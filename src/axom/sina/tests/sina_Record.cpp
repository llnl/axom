// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <stdexcept>
#include <typeinfo>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "axom/sina/core/Record.hpp"
#include "axom/sina/core/CurveSet.hpp"

#include "axom/sina/tests/SinaMatchers.hpp"
#include "axom/sina/tests/TestRecord.hpp"

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

char const EXPECTED_TYPE_KEY[] = "type";
char const EXPECTED_LOCAL_ID_KEY[] = "local_id";
char const EXPECTED_GLOBAL_ID_KEY[] = "id";
char const EXPECTED_DATA_KEY[] = "data";
char const EXPECTED_LIBRARY_DATA_KEY[] = "library_data";
char const EXPECTED_FILES_KEY[] = "files";
char const EXPECTED_USER_DEFINED_KEY[] = "user_defined";
char const LIBRARY_DATA_ID_DATUM[] = "SINA_librarydata_id";
char const LIBRARY_DATA_TYPE_DATUM[] = "SINA_librarydata_type";

TEST(Record, create_typeMissing)
{
  conduit::Node originalNode;
  originalNode[EXPECTED_LOCAL_ID_KEY] = "the ID";
  try
  {
    Record record {originalNode};
    FAIL() << "Should have failed due to missing type";
  }
  catch(std::invalid_argument const &expected)
  {
    EXPECT_THAT(expected.what(), HasSubstr(EXPECTED_TYPE_KEY));
  }
}

TEST(Record, add_data_existing_key)
{
  Record record {ID {"the id", IDType::Local}, "test_record"};
  record.add("key1", Datum {"val1"});
  EXPECT_EQ("val1", record.getData().at("key1").getValue());
  record.add("key1", Datum {"val2"});
  EXPECT_EQ("val2", record.getData().at("key1").getValue());
}

TEST(Record, add_curve_set_existing_key)
{
  Record record {ID {"the id", IDType::Local}, "test_record"};

  CurveSet cs1 {"cs1"};
  cs1.addDependentCurve(Curve {"original", {1, 2, 3}});
  record.add(cs1);

  auto &csAfterFirstInsert = record.getCurveSets();
  ASSERT_THAT(csAfterFirstInsert, Contains(Key("cs1")));
  EXPECT_THAT(csAfterFirstInsert.at("cs1").getDependentCurves(), Contains(Key("original")));

  CurveSet cs2 {"cs1"};
  cs2.addDependentCurve(Curve {"new", {1, 2, 3}});
  record.add(cs2);

  auto &csAfterSecondInsert = record.getCurveSets();
  ASSERT_THAT(csAfterSecondInsert, Contains(Key("cs1")));
  EXPECT_THAT(csAfterSecondInsert.at("cs1").getDependentCurves(), Not(Contains(Key("original"))));
  EXPECT_THAT(csAfterSecondInsert.at("cs1").getDependentCurves(), Contains(Key("new")));
}

TEST(Record, remove_file)
{
  Record record {ID {"the id", IDType::Local}, "test_record"};
  std::string path = "the/path.txt";

  File original {path};
  original.setMimeType("txt");
  record.add(original);
  EXPECT_EQ(1u, record.getFiles().size());
  EXPECT_EQ("txt", record.getFiles().find(File {path})->getMimeType());

  record.remove(original);
  EXPECT_EQ(0u, record.getFiles().size());
}

TEST(Record, add_file_existing_key)
{
  Record record {ID {"the id", IDType::Local}, "test_record"};
  std::string path = "the/path.txt";

  File original {path};
  original.setMimeType("txt");
  record.add(original);
  EXPECT_EQ(1u, record.getFiles().size());
  EXPECT_EQ("txt", record.getFiles().find(File {path})->getMimeType());

  File replacement {path};
  replacement.setMimeType("image");
  record.add(replacement);
  EXPECT_EQ(1u, record.getFiles().size());
  EXPECT_EQ("image", record.getFiles().find(File {path})->getMimeType());
}

TEST(Record, add_child_record_as_library_data)
{
  Record parentRecord {ID {"parent id", IDType::Local}, "test_record_parent"};
  Record childRecord {ID {"child id", IDType::Local}, "test_record_child"};
  parentRecord.addRecordAsLibraryData(childRecord, "child");
  auto &parentLibData = parentRecord.getLibraryData();
  ASSERT_THAT(parentLibData, Contains(Key("child")));
  auto &childLibContents = parentLibData.at("child")->getData();
  ASSERT_THAT(childLibContents, Contains(Key(LIBRARY_DATA_ID_DATUM)));
  EXPECT_EQ("child id", childLibContents.at(LIBRARY_DATA_ID_DATUM).getValue());
  ASSERT_THAT(childLibContents, Contains(Key(LIBRARY_DATA_TYPE_DATUM)));
  EXPECT_EQ("test_record_child", childLibContents.at(LIBRARY_DATA_TYPE_DATUM).getValue());
}

TEST(Record, add_child_record_as_library_data_with_data)
{
  Record parentRecord {ID {"parent id", IDType::Local}, "test_record_parent"};
  Record childRecord {ID {"child id", IDType::Local}, "test_record_child"};
  childRecord.add("key1", Datum {"val1"});
  parentRecord.addRecordAsLibraryData(childRecord, "child");
  auto &childLibContents = parentRecord.getLibraryData().at("child")->getData();
  ASSERT_THAT(childLibContents, Contains(Key("key1")));
  EXPECT_EQ("val1", childLibContents.at("key1").getValue());
}

TEST(Record, add_child_record_as_library_data_with_files)
{
  Record parentRecord {ID {"parent id", IDType::Local}, "test_record_parent"};
  Record childRecord {ID {"child id", IDType::Local}, "test_record_child"};
  std::string path = "the/path.txt";
  File childFile {path};
  childFile.setMimeType("txt");
  childRecord.add(childFile);
  parentRecord.addRecordAsLibraryData(childRecord, "child");
  EXPECT_EQ(1u, parentRecord.getFiles().size());
  EXPECT_EQ("txt", parentRecord.getFiles().find(File {path})->getMimeType());
}

TEST(Record, add_child_record_as_library_data_with_curves)
{
  Record parentRecord {ID {"parent id", IDType::Local}, "test_record_parent"};
  Record childRecord {ID {"child id", IDType::Local}, "test_record_child"};
  CurveSet cs {"name"};
  cs.addIndependentCurve(Curve {"lime", {1.0, 2.0, 3.0}});
  cs.addIndependentCurve(Curve {"white", {1.0, 2.0, 3.0}});
  cs.addIndependentCurve(Curve {"black", {1.0, 2.0, 3.0}});
  cs.addDependentCurve(Curve {"lightgrey", {4.0, 5.0, 6.0}});
  cs.addDependentCurve(Curve {"maroon", {7.0, 8.0, 9.0}});
  cs.addDependentCurve(Curve {"brown", {1.0, 2.0, 3.0}});
  std::vector<std::string> expected_order {"lightgrey", "brown", "maroon"};
  cs.applyCustomDependentCurveOrder(expected_order);
  childRecord.add(cs);
  parentRecord.addRecordAsLibraryData(childRecord, "child");
  auto expected = R"({
      "local_id": "parent id",
      "type": "test_record_parent",
      "library_data":  {
        "child": {
          "curve_sets": {
            "name": {
              "independent": {
                  "black": { "value": [1.0, 2.0, 3.0]},
                  "lime": { "value": [1.0, 2.0, 3.0]},
                  "white": { "value": [1.0, 2.0, 3.0]}
              },
              "dependent": {
                  "lightgrey": { "value": [4.0, 5.0, 6.0]},
                  "brown": { "value": [1.0, 2.0, 3.0]},
                  "maroon": { "value": [7.0, 8.0, 9.0]}
              }
          }
      },
      "data": { "SINA_librarydata_type": { "value": "test_record_child"},
                "SINA_librarydata_id": { "value": "child id"}}
   }}})";
  // gmock's Json matcher doesn't check JSON object order, as it's not expected to be ordered; we're relying on
  // conduit to give us that, but it's still useful to make sure we didn't alter the structure of the rest of the rec
  EXPECT_THAT(parentRecord.toNode(), MatchesJsonMatcher(expected));
  EXPECT_EQ(
    parentRecord.toNode()["library_data"]["child"]["curve_sets"]["name"]["dependent"].child_names(),
    expected_order);
  std::vector<std::string> expected_reorder {"brown", "lightgrey", "maroon"};
  EXPECT_EQ(
    parentRecord
      .toNode(axom::sina::CurveSet::CurveOrder::ALPHABETIC)["library_data"]["child"]["curve_sets"]
                                                           ["name"]["dependent"]
      .child_names(),
    expected_reorder);
}

TEST(Record, create_localId_fromNode)
{
  conduit::Node originalNode;
  originalNode[EXPECTED_LOCAL_ID_KEY] = "the ID";
  originalNode[EXPECTED_TYPE_KEY] = "my type";
  Record record {originalNode};
  EXPECT_EQ("my type", record.getType());
  EXPECT_EQ("the ID", record.getId().getId());
  EXPECT_EQ(IDType::Local, record.getId().getType());
}

TEST(Record, create_globalId_fromNode)
{
  conduit::Node originalNode;
  originalNode[EXPECTED_GLOBAL_ID_KEY] = "the ID";
  originalNode[EXPECTED_TYPE_KEY] = "my type";
  Record record {originalNode};
  EXPECT_EQ("my type", record.getType());
  EXPECT_EQ("the ID", record.getId().getId());
  EXPECT_EQ(IDType::Global, record.getId().getType());
}

TEST(Record, create_globalId_withContent)
{
  conduit::Node originalNode;
  originalNode[EXPECTED_GLOBAL_ID_KEY] = "the ID";
  originalNode[EXPECTED_TYPE_KEY] = "my type";
  originalNode[EXPECTED_DATA_KEY];
  originalNode[EXPECTED_LIBRARY_DATA_KEY];

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

  std::string libName = "my_lib";
  conduit::Node libNode;
  std::string name3 = "datum name 3";
  conduit::Node name3_node;
  name3_node["value"] = "value 3";
  libNode[EXPECTED_DATA_KEY];
  libNode[EXPECTED_DATA_KEY][name3] = name3_node;
  originalNode[EXPECTED_LIBRARY_DATA_KEY][libName] = libNode;

  Record record {originalNode};
  auto &data = record.getData();
  ASSERT_EQ(2u, data.size());
  EXPECT_EQ("value 1", data.at(name1).getValue());
  EXPECT_THAT(2.22, DoubleEq(data.at(name2).getScalar()));
  EXPECT_EQ("g/L", data.at(name2).getUnits());
  EXPECT_EQ("tag1", data.at(name2).getTags()[0]);
  EXPECT_EQ("tag2", data.at(name2).getTags()[1]);

  auto &libdata = record.getLibraryData();
  EXPECT_THAT(libdata, Contains(Key(libName)));
  EXPECT_EQ("value 3", libdata.at(libName)->getData().at(name3).getValue());
}

TEST(Record, create_globalId_files)
{
  conduit::Node originalNode;
  originalNode[EXPECTED_GLOBAL_ID_KEY] = "the ID";
  originalNode[EXPECTED_TYPE_KEY] = "my type";
  originalNode[EXPECTED_FILES_KEY];

  std::string uri1 = "/some/uri.txt";
  std::string uri2 = "www.anotheruri.com";
  std::string uri3 = "yet another uri";
  originalNode[EXPECTED_FILES_KEY].add_child(uri1);
  originalNode[EXPECTED_FILES_KEY].add_child(uri2);
  originalNode[EXPECTED_FILES_KEY].add_child(uri3);
  Record record {originalNode};
  auto &files = record.getFiles();
  ASSERT_EQ(3u, files.size());
  EXPECT_EQ(1, files.count(File {uri1}));
  EXPECT_EQ(1, files.count(File {uri2}));
  EXPECT_EQ(1, files.count(File {uri3}));
}

TEST(Record, create_fromNode_curveSets)
{
  conduit::Node recordAsNode = parseJsonValue(R"({
        "id": "myId",
        "type": "myType",
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
  Record record {recordAsNode};
  auto &curveSets = record.getCurveSets();
  ASSERT_THAT(curveSets, Contains(Key("cs1")));
}

TEST(Record, create_fromNode_userDefined)
{
  conduit::Node originalNode;
  originalNode[EXPECTED_GLOBAL_ID_KEY] = "the ID";
  originalNode[EXPECTED_TYPE_KEY] = "my type";
  originalNode[EXPECTED_USER_DEFINED_KEY]["k1"] = "v1";
  originalNode[EXPECTED_USER_DEFINED_KEY]["k2"] = 123;
  std::vector<int> k3_vals {1, 2, 3};
  originalNode[EXPECTED_USER_DEFINED_KEY]["k3"] = k3_vals;

  Record record {originalNode};
  auto const &userDefined = record.getUserDefinedContent();
  EXPECT_EQ("v1", userDefined["k1"].as_string());
  EXPECT_EQ(123, userDefined["k2"].as_int());
  auto int_array = userDefined["k3"].as_int_ptr();
  std::vector<double> udef_ints(int_array,
                                int_array + userDefined["k3"].dtype().number_of_elements());
  EXPECT_THAT(udef_ints, ElementsAre(1, 2, 3));
}

TEST(Record, getUserDefined_initialConst)
{
  ID id {"the id", IDType::Local};
  Record const record {id, "my type"};
  conduit::Node const &userDefined = record.getUserDefinedContent();
  EXPECT_TRUE(userDefined.dtype().is_empty());
}

TEST(Record, getUserDefined_initialNonConst)
{
  ID id {"the id", IDType::Local};
  Record record {id, "my type"};
  conduit::Node &initialUserDefined = record.getUserDefinedContent();
  EXPECT_TRUE(initialUserDefined.dtype().is_empty());
  initialUserDefined["foo"] = 123;
  EXPECT_EQ(123, record.getUserDefinedContent()["foo"].as_int());
}

TEST(Record, toNode_localId)
{
  ID id {"the id", IDType::Global};
  Record record {id, "my type"};
  auto asNode = record.toNode();
  EXPECT_TRUE(asNode.dtype().is_object());
  EXPECT_EQ("my type", asNode[EXPECTED_TYPE_KEY].as_string());
  EXPECT_EQ("the id", asNode[EXPECTED_GLOBAL_ID_KEY].as_string());
  EXPECT_TRUE(asNode[EXPECTED_LOCAL_ID_KEY].dtype().is_empty());
}

TEST(Record, toNode_globalId)
{
  ID id {"the id", IDType::Local};
  Record record {id, "my type"};
  auto asNode = record.toNode();
  EXPECT_TRUE(asNode.dtype().is_object());
  EXPECT_EQ("my type", asNode[EXPECTED_TYPE_KEY].as_string());
  EXPECT_EQ("the id", asNode[EXPECTED_LOCAL_ID_KEY].as_string());
  EXPECT_TRUE(asNode[EXPECTED_GLOBAL_ID_KEY].dtype().is_empty());
}

TEST(Record, toNode_default_values)
{
  ID id {"the id", IDType::Global};
  Record record {id, "my type"};
  auto asNode = record.toNode();
  EXPECT_TRUE(asNode.dtype().is_object());
  // We want to be sure that unset optional fields aren't present
  EXPECT_FALSE(asNode.has_child(EXPECTED_DATA_KEY));
  EXPECT_FALSE(asNode.has_child(EXPECTED_FILES_KEY));
  EXPECT_FALSE(asNode.has_child(EXPECTED_USER_DEFINED_KEY));
}

TEST(Record, toNode_userDefined)
{
  ID id {"the id", IDType::Local};
  Record record {id, "my type"};
  conduit::Node userDef;
  userDef["k1"] = "v1";
  userDef["k2"] = 123;
  std::vector<int> int_vals {1, 2, 3};
  userDef["k3"] = int_vals;
  record.setUserDefinedContent(userDef);

  auto asNode = record.toNode();

  auto userDefined = asNode[EXPECTED_USER_DEFINED_KEY];
  EXPECT_EQ("v1", userDefined["k1"].as_string());
  EXPECT_EQ(123, userDefined["k2"].as_int());
  auto int_array = userDefined["k3"].as_int_ptr();
  std::vector<double> udef_ints(int_array,
                                int_array + userDefined["k3"].dtype().number_of_elements());
  EXPECT_THAT(udef_ints, ElementsAre(1, 2, 3));
}

TEST(Record, toNode_data)
{
  ID id {"the id", IDType::Local};
  Record record {id, "my type"};
  std::string name1 = "name1";
  std::string value1 = "value1";
  Datum datum1 = Datum {value1};
  datum1.setUnits("some units");
  datum1.setTags({"tag1"});
  record.add(name1, datum1);
  std::string name2 = "name2";
  record.add(name2, Datum {2.});
  auto asNode = record.toNode();
  ASSERT_EQ(2u, asNode[EXPECTED_DATA_KEY].number_of_children());
  EXPECT_EQ("value1", asNode[EXPECTED_DATA_KEY][name1]["value"].as_string());
  EXPECT_EQ("some units", asNode[EXPECTED_DATA_KEY][name1]["units"].as_string());
  EXPECT_EQ("tag1", asNode[EXPECTED_DATA_KEY][name1]["tags"][0].as_string());

  EXPECT_THAT(asNode[EXPECTED_DATA_KEY][name2]["value"].as_double(), DoubleEq(2.));
  EXPECT_TRUE(asNode[EXPECTED_DATA_KEY][name2]["units"].dtype().is_empty());
  EXPECT_TRUE(asNode[EXPECTED_DATA_KEY][name2]["tags"].dtype().is_empty());
}

TEST(Record, toNode_dataWithSlashes)
{
  ID id {"the id", IDType::Local};
  Record record {id, "my type"};
  std::string name = "name/with/slashes";
  std::string value = "the value";
  Datum datum = Datum {value};
  record.add(name, datum);
  auto asNode = record.toNode();
  ASSERT_EQ(1u, asNode[EXPECTED_DATA_KEY].number_of_children());
  EXPECT_EQ("the value", asNode[EXPECTED_DATA_KEY].child(name)["value"].as_string());
}

TEST(Record, toNode_files)
{
  ID id {"the id", IDType::Local};
  Record record {id, "my type"};
  std::string uri1 = "a/file/path/foo.png";
  std::string uri2 = "uri2";
  File file {uri1};
  file.setMimeType("mt1");
  record.add(file);
  record.add(File {uri2});
  // Identical uris should overwrite
  record.add(File {uri2});
  auto asNode = record.toNode();
  ASSERT_EQ(2u, asNode[EXPECTED_FILES_KEY].number_of_children());
  auto &child_with_slashes = asNode[EXPECTED_FILES_KEY].child(uri1);
  EXPECT_EQ("mt1", child_with_slashes["mimetype"].as_string());
  EXPECT_TRUE(asNode[EXPECTED_FILES_KEY][uri2]["mimetype"].dtype().is_empty());
}

TEST(Record, toNode_curveSets)
{
  ID id {"the id", IDType::Local};
  Record record {id, "my type"};
  CurveSet cs {"myCurveSet/with/slash"};
  cs.addIndependentCurve(Curve {"myCurve", {1, 2, 3}});
  cs.addIndependentCurve(Curve {"myOtherCurve", {4, 5, 6}});
  cs.addIndependentCurve(Curve {"myThirdCurve", {7, 8, 9}});
  record.add(cs);
  std::string expected = R"({
        "local_id": "the id",
        "type": "my type",
        "curve_sets": {
            "myCurveSet/with/slash": {
                "independent": {
                     "myCurve": {
                         "value": [1.0, 2.0, 3.0]
                     },
                     "myOtherCurve": {
                         "value": [4.0, 5.0, 6.0]
                     },
                     "myThirdCurve": {
                         "value": [7.0, 8.0, 9.0]
                     }
                 },
                 "dependent": {}
            }
        }
    })";
  EXPECT_THAT(record.toNode(), MatchesJsonMatcher(expected));
  std::vector<std::string> expected_order {"myCurve", "myOtherCurve", "myThirdCurve"};
  EXPECT_EQ(record.toNode()["curve_sets"].child("myCurveSet/with/slash")["independent"].child_names(),
            expected_order);
}

TEST(Record, toNode_curveSets_customOrder)
{
  ID id {"the id", IDType::Local};
  Record record {id, "my type"};
  CurveSet cs {"reordered_curves"};
  cs.addIndependentCurve(Curve {"lime", {1, 2, 3}});
  cs.addIndependentCurve(Curve {"white", {4, 5, 6}});
  cs.addIndependentCurve(Curve {"black", {7, 8, 9}});
  cs.addDependentCurve(Curve {"cyan", {1, 2, 3}});
  cs.addDependentCurve(Curve {"yellow", {1, 2, 3}});
  cs.addDependentCurve(Curve {"pink", {1, 2, 3}});
  std::vector<std::string> expected_order {"white", "black", "lime"};
  cs.applyCustomIndependentCurveOrder(expected_order);
  record.add(cs);
  auto expected = R"({
        "local_id": "the id",
        "type": "my type",
        "curve_sets": {
            "reordered_curves": {
                "independent": {
                     "white": { "value": [4.0, 5.0, 6.0] },
                     "black": { "value": [7.0, 8.0, 9.0] },
                     "lime": { "value": [1.0, 2.0, 3.0] }
                 },
                 "dependent": {
                     "pink": { "value": [1.0, 2.0, 3.0] },
                     "cyan": { "value": [1.0, 2.0, 3.0] },
                     "yellow": { "value": [1.0, 2.0, 3.0] }
                 }
            }
        }
    })";
  EXPECT_THAT(record.toNode(), MatchesJsonMatcher(expected));
  EXPECT_EQ(record.toNode()["curve_sets"]["reordered_curves"]["independent"].child_names(),
            expected_order);
}

TEST(Record, toNode_curveSets_setDefaultOrder_fromRecord)
{
  ID id {"the id", IDType::Local};
  Record record {id, "my type"};
  CurveSet cs {"reordered_curves"};
  cs.addIndependentCurve(Curve {"lime", {1, 2, 3}});
  cs.addIndependentCurve(Curve {"white", {4, 5, 6}});
  cs.addIndependentCurve(Curve {"black", {7, 8, 9}});
  cs.addDependentCurve(Curve {"cyan", {1, 2, 3}});
  cs.addDependentCurve(Curve {"yellow", {1, 2, 3}});
  cs.addDependentCurve(Curve {"pink", {1, 2, 3}});
  record.add(cs);
  record.setDefaultCurveOrder(axom::sina::CurveSet::CurveOrder::ALPHABETIC);
  auto expected = R"({
        "local_id": "the id",
        "type": "my type",
        "curve_sets": {
            "reordered_curves": {
                "independent": {
                     "black": { "value": [7.0, 8.0, 9.0] },
                     "lime": { "value": [1.0, 2.0, 3.0] },
                     "white": { "value": [4.0, 5.0, 6.0] }
                 },
                 "dependent": {
                     "cyan": { "value": [1.0, 2.0, 3.0] },
                     "pink": { "value": [1.0, 2.0, 3.0] },
                     "yellow": { "value": [1.0, 2.0, 3.0] }
                 }
            }
        }
    })";
  EXPECT_THAT(record.toNode(), MatchesJsonMatcher(expected));
  std::vector<std::string> expected_order {"black", "lime", "white"};
  EXPECT_EQ(record.toNode()["curve_sets"]["reordered_curves"]["independent"].child_names(),
            expected_order);
  std::vector<std::string> dep_expected_order {"cyan", "pink", "yellow"};
  EXPECT_EQ(record.toNode()["curve_sets"]["reordered_curves"]["dependent"].child_names(),
            dep_expected_order);
  std::vector<std::string> expected_reverse_order {"white", "lime", "black"};
  EXPECT_EQ(
    record
      .toNode(
        CurveSet::CurveOrder::REVERSE_ALPHABETICAL)["curve_sets"]["reordered_curves"]["independent"]
      .child_names(),
    expected_reverse_order);
}

TEST(Record, toNode_curveSets_setDefaultOrder)
{
  ID id {"the id", IDType::Local};
  ID id2 {"the other id", IDType::Local};
  Record record {id, "my type"};
  Record record2 {id2, "my type"};
  CurveSet cs {"reordered_curves"};
  axom::sina::setDefaultCurveOrder(axom::sina::CurveSet::CurveOrder::ALPHABETIC);
  record2.setDefaultCurveOrder(axom::sina::CurveSet::CurveOrder::REVERSE_ALPHABETICAL);

  cs.addIndependentCurve(Curve {"lime", {1, 2, 3}});
  cs.addIndependentCurve(Curve {"white", {4, 5, 6}});
  cs.addIndependentCurve(Curve {"black", {7, 8, 9}});
  cs.addDependentCurve(Curve {"cyan", {1, 2, 3}});
  cs.addDependentCurve(Curve {"yellow", {1, 2, 3}});
  cs.addDependentCurve(Curve {"pink", {1, 2, 3}});
  record.add(cs);
  record2.add(cs);
  auto expected = R"({
        "local_id": "the id",
        "type": "my type",
        "curve_sets": {
            "reordered_curves": {
                "independent": {
                     "black": { "value": [7.0, 8.0, 9.0] },
                     "lime": { "value": [1.0, 2.0, 3.0] },
                     "white": { "value": [4.0, 5.0, 6.0] }
                 },
                 "dependent": {
                     "cyan": { "value": [1.0, 2.0, 3.0] },
                     "pink": { "value": [1.0, 2.0, 3.0] },
                     "yellow": { "value": [1.0, 2.0, 3.0] }
                 }
            }
        }
    })";
  EXPECT_THAT(record.toNode(), MatchesJsonMatcher(expected));
  std::vector<std::string> expected_order {"black", "lime", "white"};
  EXPECT_EQ(record.toNode()["curve_sets"]["reordered_curves"]["independent"].child_names(),
            expected_order);
  std::vector<std::string> expected_reverse_order {"white", "lime", "black"};
  EXPECT_EQ(record2.toNode()["curve_sets"]["reordered_curves"]["independent"].child_names(),
            expected_reverse_order);
  std::vector<std::string> dep_expected_order {"cyan", "pink", "yellow"};
  EXPECT_EQ(record.toNode()["curve_sets"]["reordered_curves"]["dependent"].child_names(),
            dep_expected_order);
}

TEST(RecordLoader, load_missingLoader)
{
  RecordLoader loader;
  conduit::Node asNode;
  asNode[EXPECTED_GLOBAL_ID_KEY] = "the ID";
  asNode[EXPECTED_TYPE_KEY] = "unknownType";
  auto loaded = loader.load(asNode);

  EXPECT_NE(loaded.get(), nullptr);
  if(loaded)
  {
    auto &loadedRef = *loaded;
    EXPECT_EQ(typeid(Record), typeid(loadedRef)) << "Type was " << typeid(loadedRef).name();
  }
}

TEST(RecordLoader, load_loaderPresent)
{
  RecordLoader loader;
  EXPECT_FALSE(loader.canLoad("TestInt"));
  EXPECT_FALSE(loader.canLoad("TestString"));

  loader.addTypeLoader("TestInt", [](conduit::Node const &value) {
    return std::make_unique<TestRecord<int>>(value);
  });
  EXPECT_TRUE(loader.canLoad("TestInt"));

  loader.addTypeLoader("TestString", [](conduit::Node const &value) {
    return std::make_unique<TestRecord<std::string>>(value);
  });
  EXPECT_TRUE(loader.canLoad("TestString"));

  conduit::Node asNode;
  asNode[EXPECTED_GLOBAL_ID_KEY] = "the ID";
  asNode[EXPECTED_TYPE_KEY] = "TestString";
  asNode[TEST_RECORD_VALUE_KEY] = "The value";
  auto loaded = loader.load(asNode);
  auto testObjPointer = dynamic_cast<TestRecord<std::string> *>(loaded.get());
  ASSERT_NE(nullptr, testObjPointer);
  EXPECT_EQ("The value", testObjPointer->getValue());
  EXPECT_EQ("TestString", testObjPointer->getType());
}

TEST(RecordLoader, createRecordLoaderWithAllKnownTypes)
{
  RecordLoader loader = createRecordLoaderWithAllKnownTypes();
  EXPECT_TRUE(loader.canLoad("run"));
}

}  // namespace
}  // namespace testing
}  // namespace sina
}  // namespace axom
