// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <stdexcept>
#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "axom/sina/core/Datum.hpp"
#include "axom/sina/core/ConduitUtil.hpp"

namespace axom
{
namespace sina
{
namespace testing
{
namespace
{

using ::testing::DoubleEq;
using ::testing::ElementsAre;
using ::testing::HasSubstr;

TEST(Datum, create)
{
  std::vector<std::string> tags = {"tag1", "tag2"};
  std::string value = "value";
  std::vector<std::string> val_list = {"val1", "val2"};
  std::vector<double> scal_list = {100, 2.0};
  Datum datum1 {value};
  datum1.setUnits("some units");
  datum1.setTags(tags);
  Datum datum2 {3.14};
  Datum datum3 {val_list};
  Datum datum4 {scal_list};

  EXPECT_EQ(ValueType::String, datum1.getType());
  EXPECT_EQ("value", datum1.getValue());
  EXPECT_EQ("some units", datum1.getUnits());
  EXPECT_EQ(tags, datum1.getTags());

  EXPECT_EQ(ValueType::Scalar, datum2.getType());
  EXPECT_THAT(datum2.getScalar(), DoubleEq(3.14));

  EXPECT_EQ(ValueType::StringArray, datum3.getType());
  EXPECT_EQ(val_list, datum3.getStringArray());

  EXPECT_EQ(ValueType::ScalarArray, datum4.getType());
  EXPECT_EQ(scal_list, datum4.getScalarArray());
}

TEST(Datum, createFromNode)
{
  conduit::Node val_with_tags;
  conduit::Node scalar_with_units;
  conduit::Node val_list_node;
  conduit::Node scal_list_node;
  conduit::Node empty_list_node;
  conduit::Node int_array_node;
  conduit::Node char_array_node;

  std::vector<std::string> tags = {"hello", "world"};
  std::vector<std::string> val_list = {"val1", "val2"};
  std::vector<double> scal_list = {100, 2.0};
  // Conduit has some special treatment of arrays
  int int_array[] = {-2, 2, 4, 8};
  std::vector<double> array_equiv = {-2, 2, 4, 8};
  char coerced_to_string[] = {'a', 'b', 'c', '\0'};
  //Empty lists are valid
  conduit::Node empty_list(conduit::DataType::list());

  val_with_tags["value"] = "the value";
  addStringsToNode(val_with_tags, "tags", tags);
  scalar_with_units["units"] = "some units";
  scalar_with_units["value"] = 3.14;
  addStringsToNode(val_list_node, "value", val_list);
  scal_list_node["value"] = scal_list;
  empty_list_node["value"] = empty_list;
  int_array_node["value"].set(int_array, 4);
  char_array_node["value"] = coerced_to_string;

  Datum datum1 {val_with_tags};
  Datum datum2 {scalar_with_units};
  Datum datum3 {val_list_node};
  Datum datum4 {scal_list_node};
  Datum datum5 {empty_list_node};
  Datum datum6 {int_array_node};
  Datum datum7 {char_array_node};

  EXPECT_EQ("the value", datum1.getValue());
  EXPECT_EQ(tags, datum1.getTags());
  EXPECT_THAT(3.14, DoubleEq(datum2.getScalar()));
  EXPECT_EQ("some units", datum2.getUnits());
  EXPECT_EQ(val_list, datum3.getStringArray());
  EXPECT_EQ(scal_list, datum4.getScalarArray());
  EXPECT_EQ(ValueType::ScalarArray, datum5.getType());
  EXPECT_EQ(array_equiv, datum6.getScalarArray());
  EXPECT_EQ("abc", datum7.getValue());
}

TEST(Datum, setUnits)
{
  std::string value = "value";
  Datum datum1 {value};
  datum1.setUnits("new units");
  EXPECT_EQ("new units", datum1.getUnits());
}

TEST(Datum, setTags)
{
  std::string value = "value";
  Datum datum1 {value};
  datum1.setTags({"new_tag"});
  EXPECT_EQ("new_tag", datum1.getTags()[0]);
}

TEST(Datum, createFromJson_missingKeys)
{
  conduit::Node object1;
  try
  {
    Datum datum1 {object1};
    FAIL() << "Should have gotten a value error";
  }
  catch(std::invalid_argument const &expected)
  {
    EXPECT_THAT(expected.what(), HasSubstr("value"));
  }
}

TEST(Datum, createFromJson_badListValue)
{
  conduit::Node object1;
  auto &mixed_scal = object1["value"].append();
  mixed_scal.set(1.0);
  auto &mixed_val = object1["value"].append();
  mixed_val.set("two");
  try
  {
    Datum datum1 {object1};
    FAIL() << "Should have gotten a value error";
  }
  catch(std::invalid_argument const &expected)
  {
    std::string warning = "it must consist of only strings or only numbers";
    EXPECT_THAT(expected.what(), HasSubstr(warning));
  }
}

TEST(Datum, toJson)
{
  std::vector<std::string> tags = {"list", "of", "tags"};
  std::string value = "Datum value";
  std::vector<double> scal_list = {-14, 22, 9};
  std::vector<std::string> val_list = {"east", "west"};
  Datum datum1 {value};
  datum1.setTags(tags);
  Datum datum2 {3.14};
  datum2.setUnits("Datum units");
  Datum datum3 {scal_list};
  Datum datum4 {val_list};
  conduit::Node datumRef1 = datum1.toNode();
  conduit::Node datumRef2 = datum2.toNode();
  conduit::Node datumRef3 = datum3.toNode();
  conduit::Node datumRef4 = datum4.toNode();
  EXPECT_EQ("Datum value", datumRef1["value"].as_string());
  std::vector<std::string> node_tags;
  auto tags_itr = datumRef1["tags"].children();
  while(tags_itr.has_next()) node_tags.emplace_back(tags_itr.next().as_string());
  EXPECT_EQ(tags, node_tags);

  EXPECT_EQ("Datum units", datumRef2["units"].as_string());
  EXPECT_THAT(3.14, DoubleEq(datumRef2["value"].value()));

  // Conduit will pack vectors of numbers into arrays, but
  // strings can only live as lists of Nodes
  auto doub_array = datumRef3["value"].as_double_ptr();
  std::vector<double> scal_child_vals(doub_array,
                                      doub_array + datumRef3["value"].dtype().number_of_elements());
  std::vector<std::string> str_child_vals;
  auto str_itr = datumRef4["value"].children();
  while(str_itr.has_next()) str_child_vals.emplace_back(str_itr.next().as_string());
  EXPECT_EQ(scal_list, scal_child_vals);
  EXPECT_EQ(val_list, str_child_vals);
}

}  // namespace
}  // namespace testing
}  // namespace sina
}  // namespace axom
