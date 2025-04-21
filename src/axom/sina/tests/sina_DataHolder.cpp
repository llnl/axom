// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <stdexcept>
#include <typeinfo>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "axom/sina/core/DataHolder.hpp"
#include "axom/sina/tests/SinaMatchers.hpp"

namespace axom
{
namespace sina
{
namespace testing
{
namespace
{

using ::testing::Contains;
using ::testing::DoubleEq;
using ::testing::ElementsAre;
using ::testing::HasSubstr;
using ::testing::Key;
using ::testing::Not;

char const EXPECTED_DATA_KEY[] = "data";
char const EXPECTED_CURVE_SETS_KEY[] = "curve_sets";
char const EXPECTED_LIBRARY_DATA_KEY[] = "library_data";
char const EXPECTED_USER_DEFINED_KEY[] = "user_defined";

// Controlled from Record in real usage, but we set it here to keep
// the test cases readable, since it's not usually the focus
CurveSet::CurveOrder const defaultCurveOrder = CurveSet::CurveOrder::REGISTRATION_OLDEST_FIRST;

TEST(DataHolder, add_data_existing_key)
{
  DataHolder dh {};
  dh.add("key1", Datum {"val1"});
  EXPECT_EQ("val1", dh.getData().at("key1").getValue());
  dh.add("key1", Datum {"val2"});
  EXPECT_EQ("val2", dh.getData().at("key1").getValue());
}

TEST(DataHolder, add_curve_set_existing_key)
{
  DataHolder dh {};
  CurveSet cs1 {"cs1"};
  cs1.addDependentCurve(Curve {"original", {1, 2, 3}});
  dh.add(cs1);

  auto &csAfterFirstInsert = dh.getCurveSets();
  ASSERT_THAT(csAfterFirstInsert, Contains(Key("cs1")));
  EXPECT_THAT(csAfterFirstInsert.at("cs1").getDependentCurves(), Contains(Key("original")));

  CurveSet cs2 {"cs1"};
  cs2.addDependentCurve(Curve {"new", {1, 2, 3}});
  dh.add(cs2);

  auto &csAfterSecondInsert = dh.getCurveSets();
  ASSERT_THAT(csAfterSecondInsert, Contains(Key("cs1")));
  EXPECT_THAT(csAfterSecondInsert.at("cs1").getDependentCurves(), Not(Contains(Key("original"))));
  EXPECT_THAT(csAfterSecondInsert.at("cs1").getDependentCurves(), Contains(Key("new")));
}

TEST(DataHolder, create_fromNode_userDefined)
{
  conduit::Node originalNode;
  originalNode[EXPECTED_USER_DEFINED_KEY]["k1"] = "v1";
  originalNode[EXPECTED_USER_DEFINED_KEY]["k2"] = 123;
  std::vector<int> k3_vals {1, 2, 3};
  originalNode[EXPECTED_USER_DEFINED_KEY]["k3"] = k3_vals;

  DataHolder holder {originalNode};
  auto const &userDefined = holder.getUserDefinedContent();
  EXPECT_EQ("v1", userDefined["k1"].as_string());
  EXPECT_EQ(123, userDefined["k2"].as_int());
  auto int_array = userDefined["k3"].as_int_ptr();
  std::vector<double> udef_ints(int_array,
                                int_array + userDefined["k3"].dtype().number_of_elements());
  EXPECT_THAT(udef_ints, ElementsAre(1, 2, 3));
}

TEST(DataHolder, create_fromNode_userDefined_not_object)
{
  conduit::Node originalNode;
  originalNode[EXPECTED_USER_DEFINED_KEY] = "not an object";
  EXPECT_THROW(DataHolder {originalNode}, std::invalid_argument);
}

TEST(DataHolder, getUserDefined_initialConst)
{
  DataHolder const holder;
  conduit::Node const &userDefined = holder.getUserDefinedContent();
  EXPECT_TRUE(userDefined.dtype().is_empty());
}

TEST(DataHolder, getUserDefined_initialNonConst)
{
  DataHolder holder;
  conduit::Node &initialUserDefined = holder.getUserDefinedContent();
  EXPECT_TRUE(initialUserDefined.dtype().is_empty());
  initialUserDefined["foo"] = 123;
  EXPECT_EQ(123, holder.getUserDefinedContent()["foo"].as_int());
}

TEST(DataHolder, add_new_library)
{
  DataHolder dh {};
  auto outer = dh.addLibraryData("outer");
  auto &libDataAfterFirstInsert = dh.getLibraryData();
  ASSERT_THAT(libDataAfterFirstInsert, Contains(Key("outer")));
  dh.addLibraryData("other_outer");
  auto &libDataAfterSecondInsert = dh.getLibraryData();
  ASSERT_THAT(libDataAfterSecondInsert, Contains(Key("outer")));
  ASSERT_THAT(libDataAfterSecondInsert, Contains(Key("other_outer")));
  outer->addLibraryData("inner");
  auto &libDataAfterThirdInsert = dh.getLibraryData();
  ASSERT_THAT(libDataAfterThirdInsert.at("outer")->getLibraryData(), Contains(Key("inner")));
  ASSERT_THAT(libDataAfterThirdInsert.at("other_outer")->getLibraryData(),
              Not(Contains(Key("inner"))));
}

TEST(DataHolder, add_library_existing_key)
{
  std::string libName = "outer";
  DataHolder dh {};
  auto outer = dh.addLibraryData(libName);
  outer->add("key1", Datum {"val1"});
  ASSERT_THAT(dh.getLibraryData(libName)->getData(), Contains(Key("key1")));
  dh.addLibraryData(libName);
  ASSERT_THAT(dh.getLibraryData(libName)->getData(), Not(Contains(Key("key1"))));
}

TEST(DataHolder, create_fromNode_data)
{
  conduit::Node originalNode;
  originalNode[EXPECTED_DATA_KEY];

  std::string name1 = "datum name 1";
  std::string name2 = "datum name 2/with/slash";

  conduit::Node name1_node;
  name1_node["value"] = "value 1";
  originalNode[EXPECTED_DATA_KEY][name1] = name1_node;
  conduit::Node name2_node;
  name2_node["value"] = 2.22;
  name2_node["units"] = "g/L";
  addStringsToNode(name2_node, "tags", {"tag1", "tag2"});
  name2_node["value"] = 2.22;
  originalNode[EXPECTED_DATA_KEY].add_child(name2) = name2_node;
  DataHolder dh {originalNode};
  auto &data = dh.getData();
  ASSERT_EQ(2u, data.size());
  EXPECT_EQ("value 1", data.at(name1).getValue());
  EXPECT_THAT(2.22, DoubleEq(data.at(name2).getScalar()));
  EXPECT_EQ("g/L", data.at(name2).getUnits());
  EXPECT_EQ("tag1", data.at(name2).getTags()[0]);
  EXPECT_EQ("tag2", data.at(name2).getTags()[1]);
}

TEST(DataHolder, create_fromNode_curveSets)
{
  conduit::Node dataHolderAsNode = parseJsonValue(R"({
        "curve_sets": {
            "cs1": {
                "independent": {
                    "i1": { "value": [1, 2, 3]}
                },
                "dependent": {
                    "d1": { "value": [4, 5, 6]}
                }
            }
        }
    })");
  DataHolder dh {dataHolderAsNode};
  auto &curveSets = dh.getCurveSets();
  ASSERT_THAT(curveSets, Contains(Key("cs1")));
}

TEST(DataHolder, create_fromNode_libraryData)
{
  conduit::Node dataHolderAsNode = parseJsonValue(R"({
        "library_data": {
            "outer_lib": {
                "library_data": {
                    "inner_lib": { "data": {"i2": { "value": "good morning!"}}}
                }
            }
        }
    })");
  DataHolder dh {dataHolderAsNode};
  auto &fullLibData = dh.getLibraryData();
  ASSERT_THAT(fullLibData, Contains(Key("outer_lib")));
  auto outerLibData = fullLibData.at("outer_lib")->getLibraryData();
  ASSERT_THAT(outerLibData, Contains(Key("inner_lib")));
  auto &innerData = outerLibData.at("inner_lib")->getData();
  EXPECT_EQ("good morning!", innerData.at("i2").getValue());
}

TEST(DataHolder, toNode_default_values)
{
  DataHolder dh {};
  auto asNode = dh.toNode(defaultCurveOrder);
  EXPECT_TRUE(asNode.dtype().is_object());
  // We want to be sure that unset optional fields aren't present
  EXPECT_FALSE(asNode.has_child(EXPECTED_DATA_KEY));
  EXPECT_FALSE(asNode.has_child(EXPECTED_CURVE_SETS_KEY));
  EXPECT_FALSE(asNode.has_child(EXPECTED_LIBRARY_DATA_KEY));
}

TEST(DataHolder, toNode_data)
{
  DataHolder dh {};
  std::string name1 = "name1";
  std::string value1 = "value1";
  Datum datum1 = Datum {value1};
  datum1.setUnits("some units");
  datum1.setTags({"tag1"});
  dh.add(name1, datum1);
  std::string name2 = "name2";
  dh.add(name2, Datum {2.});
  auto asNode = dh.toNode(defaultCurveOrder);
  ASSERT_EQ(2u, asNode[EXPECTED_DATA_KEY].number_of_children());
  EXPECT_EQ("value1", asNode[EXPECTED_DATA_KEY][name1]["value"].as_string());
  EXPECT_EQ("some units", asNode[EXPECTED_DATA_KEY][name1]["units"].as_string());
  EXPECT_EQ("tag1", asNode[EXPECTED_DATA_KEY][name1]["tags"][0].as_string());

  EXPECT_THAT(asNode[EXPECTED_DATA_KEY][name2]["value"].as_double(), DoubleEq(2.));
  EXPECT_TRUE(asNode[EXPECTED_DATA_KEY][name2]["units"].dtype().is_empty());
  EXPECT_TRUE(asNode[EXPECTED_DATA_KEY][name2]["tags"].dtype().is_empty());
}

TEST(DataHolder, toNode_dataWithSlashes)
{
  DataHolder dh {};
  std::string name = "name/with/slashes";
  std::string value = "the value";
  Datum datum = Datum {value};
  dh.add(name, datum);
  auto asNode = dh.toNode(defaultCurveOrder);
  ASSERT_EQ(1u, asNode[EXPECTED_DATA_KEY].number_of_children());
  EXPECT_EQ("the value", asNode[EXPECTED_DATA_KEY].child(name)["value"].as_string());
}

TEST(DataHolder, toNode_curveSets)
{
  DataHolder dh {};
  CurveSet cs {"myCurveSet/with/slash"};
  cs.addIndependentCurve(Curve {"myCurve", {1, 2, 3}});
  dh.add(cs);
  std::string expected = R"({
        "curve_sets": {
            "myCurveSet/with/slash": {
                "independent": {
                     "myCurve": {
                         "value": [1.0, 2.0, 3.0]
                     }
                 },
                 "dependent": {}
            }
        }
    })";
  EXPECT_THAT(dh.toNode(defaultCurveOrder), MatchesJsonMatcher(expected));
}

TEST(DataHolder, toNode_libraryData)
{
  DataHolder dh {};
  auto outer = dh.addLibraryData("outer");
  outer->add("scal", Datum {"goodbye!"});
  auto inner = outer->addLibraryData("inner");
  inner->add("str", Datum {"hello!"});
  std::string expected = R"({
        "library_data": {
            "outer": {
                "library_data": {
                    "inner": {
                        "data": {"str": {"value": "hello!"}}
                    }
                },
                "data": {"scal": {"value": "goodbye!"}}
            }
        }
    })";
  EXPECT_THAT(dh.toNode(defaultCurveOrder), MatchesJsonMatcher(expected));
}

TEST(DataHolder, toNode_userDefined)
{
  DataHolder holder;
  conduit::Node userDef;
  userDef["k1"] = "v1";
  userDef["k2"] = 123;
  std::vector<int> int_vals {1, 2, 3};
  userDef["k3"] = int_vals;
  holder.setUserDefinedContent(userDef);

  auto asNode = holder.toNode(defaultCurveOrder);

  auto userDefined = asNode[EXPECTED_USER_DEFINED_KEY];
  EXPECT_EQ("v1", userDefined["k1"].as_string());
  EXPECT_EQ(123, userDefined["k2"].as_int());
  auto int_array = userDefined["k3"].as_int_ptr();
  std::vector<double> udef_ints(int_array,
                                int_array + userDefined["k3"].dtype().number_of_elements());
  EXPECT_THAT(udef_ints, ElementsAre(1, 2, 3));
}

}  // namespace
}  // namespace testing
}  // namespace sina
}  // namespace axom
