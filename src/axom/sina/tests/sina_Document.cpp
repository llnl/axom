// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core/utilities/FileUtilities.hpp"

#include "axom/sina/core/Document.hpp"
#include "axom/sina/core/Run.hpp"
#include "axom/sina/tests/TestRecord.hpp"

#include "conduit.hpp"
#include "conduit_relay.hpp"
#include "conduit_relay_io.hpp"
#ifdef AXOM_USE_HDF5
  #include "conduit_relay_io_hdf5.hpp"
#endif

#include <cstdio>
#include <fstream>
#include <iostream>
#include <type_traits>
#include <utility>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "axom/sina/tests/TestRecord.hpp"
#include "axom/sina/core/CurveSet.hpp"
#include "axom/sina/core/Record.hpp"
#include "axom/sina/tests/SinaMatchers.hpp"

namespace axom
{
namespace sina
{
namespace testing
{
namespace
{

using ::testing::ElementsAre;
using ::testing::HasSubstr;

char const TEST_RECORD_TYPE[] = "test type";
char const EXPECTED_RECORDS_KEY[] = "records";
char const EXPECTED_RELATIONSHIPS_KEY[] = "relationships";

// Simple document to append into, has examples of all our data categories
std::string SIMPLE_DOCUMENT = R"(
{
  "records": [
    {
      "type": "run",
      "application": "test",
      "local_id": "bar1",
      "data": {
        "int": {
          "value": 500,
          "units": "miles"
        },
        "string": {"value": "goodbye!"},
        "str/ings": {
          "value": ["z", "o", "o"]
        }
      },
      "files": {
        "test/test.png": {}
      },
      "user_defined": {
        "foo": "bar"
      },
      "curve_sets": {
        "set_1": {
          "dependent": {
            "0": {"value": [1, 2, 3]},
            "1": {"value": [-1, -2, -3]}
          },
          "independent": {
            "0": {"value": [4, 5, 6]},
            "1": {"value": [-4, -5, -6]}
          }
        },
        "set_2": {
          "dependent": {
            "0": {"value": [10, 10]}
          },
          "independent": {
            "0": {"value": [20, 20]}
          }        
        }
      },
      "library_data": {
        "my_lib": {
          "data": {"int": {"value": 10}},
          "library_data": {
            "my_inner_lib": {
              "user_defined": {"foo/bar": "baz/qux"}
            }
          }
        }
      }
    }
  ],
  "relationships": [{"local_subject": "bar1", "predicate": "knows", "object": "something"},
                    {"local_subject": "bar1", "predicate": "this_is", "object": "unique"}]
}
)";

std::string CURVE_ORDERED_DOCUMENT = R"(
{
  "records": [
    {
      "type": "run",
      "application": "test",
      "local_id": "bar1",
      "curve_sets": {
        "set_1": {
          "dependent": {
            "i_am_first": {"value": [1.4, 2.6, 3.2]},
            "must_be_second": {"value": [-1.4, -2.0, -3.7]},
            "ordered_third": {"value": [70.0, 80.1125, 90.0]}
          },
          "independent": {
            "time": {"value": [4, 5, 6]}
          }
        }
      }
    }
  ],
  "relationships": []
}
)";

std::string CURVE_ORDERED_DOCUMENT_APPEND = R"(
{
  "records": [{
      "type": "run",  "application": "test", "local_id": "bar1",
      "curve_sets": {
        "set_1": {
          "dependent": {
            "am_i_fourth": {"value": [-100.25, -110.18, -120.4]}
    }}}}], "relationships": []}
)";

std::string MULTI_REC_DOCUMENT = R"(
{
      "records": [
          {
              "data": {"string": {"value": "hello!"},
                       "string2": {"value": "unchanged!"},
                       "int": {"value": 20}},
              "type": "run",
              "application": "test",
              "local_id": "bar1",
              "user_defined":{"hello": "and",
                              "foo": "notbar"},
              "curve_sets": {
                  "set_1": {
                      "dependent": {
                          "0": { "value": [11.0, 12.0, 13.0] },
                          "1": { "value": [10.0, 20.0, 30.0] }
                      },
                      "independent": {
                          "0": { "value": [14.0, 15.0, 16.0] },
                          "1": { "value": [7.0, 8.0, 9.0] }
                      }
                  }
              }
          },
          {
              "type": "run",
              "application": "test",
              "id": "bar2",
              "curve_sets": {
                  "set_1": {
                      "dependent": {
                          "0": { "value": [1.0, 2.0] },
                          "1": { "value": [10.0, 20.0] }
                      },
                      "independent": {
                          "0": { "value": [4.0, 5.0] },
                          "1": { "value": [7.0, 8.0] }
                      }
                  }
              }
          }
      ],
      "relationships": [{"local_subject": "bar1", "predicate": "knows", "object": "something"},
                        {"local_subject": "bar2", "predicate": "sees", "object": "bar1"}]
  })";

// Full-featured multi-record document to test appending into.
std::string long_json = R"(
{
  "records": [
    {
      "type": "foo",
      "id": "test_1",
      "user_defined": {
        "name": "bob"
      },
      "files": {
        "foo/bar.png": {
          "mimetype": "image"
        }
      },
      "data": {
        "scalar": {
          "value": 500,
          "units": "miles"
        }
      }
    },
    {
      "type": "bar",
      "id": "test_2",
      "data": {
        "scalar_list": {
          "value": [1, 2, 3]
        },
        "string_list": {
          "value": ["a", "wonderful", "world"],
          "tags": ["observation"]
        }
      }
    },
    {
      "type": "run",
      "application": "sina_test",
      "id": "test_3",
      "data": {
        "scalar": {
          "value": 12.3,
          "units": "g/s",
          "tags": ["hi"]
        },
        "scalar_list": {
          "value": [1, 2, 3.0, 4]
        }
      }
    },
    {
      "type": "bar",
      "id": "test_4",
      "data": {
        "string": {
          "value": "yarr"
        },
        "string_list": {
          "value": ["y", "a", "r"]
        }
      },
      "files": {
        "test/test.png": {}
      },
      "user_defined": {
        "hello": "there"
      }
    }
  ],
  "relationships": [
    {
      "predicate": "completes",
      "subject": "test_2",
      "object": "test_1"
    },
    {
      "subject": "test_3",
      "predicate": "overrides",
      "object": "test_4"
    }
  ]
}
)";

// Helper function to convert Conduit Node array to std::vector<double> for HDF5 assertion
std::vector<double> node_to_double_vector(const conduit::Node &node)
{
  std::vector<double> result;

  if(node.dtype().is_number())
  {
    const double *intArray = node.as_double_ptr();
    conduit::index_t numElements = node.dtype().number_of_elements();
    for(conduit::index_t i = 0; i < numElements; ++i)
    {
      result.push_back(intArray[i]);
    }
  }
  return result;
}

// Tests
TEST(Document, create_fromNode_empty)
{
  conduit::Node documentAsNode;
  RecordLoader loader;
  Document document {documentAsNode, loader};
  EXPECT_EQ(0u, document.getRecords().size());
  EXPECT_EQ(0u, document.getRelationships().size());
}

TEST(Document, create_fromNode_wrongRecordsType)
{
  conduit::Node recordsAsNodes;
  recordsAsNodes[EXPECTED_RECORDS_KEY] = 123;
  RecordLoader loader;
  try
  {
    Document document {recordsAsNodes, loader};
    FAIL() << "Should not have been able to parse records. Have " << document.getRecords().size();
  }
  catch(std::invalid_argument const &expected)
  {
    EXPECT_THAT(expected.what(), HasSubstr(EXPECTED_RECORDS_KEY));
  }
}

TEST(Document, create_fromNode_withRecords)
{
  conduit::Node recordAsNode;
  recordAsNode["type"] = "IntTestRecord";
  recordAsNode["id"] = "the ID";
  recordAsNode[TEST_RECORD_VALUE_KEY] = 123;

  conduit::Node recordsAsNodes;
  recordsAsNodes.append().set(recordAsNode);

  conduit::Node documentAsNode;
  documentAsNode[EXPECTED_RECORDS_KEY] = recordsAsNodes;

  RecordLoader loader;
  loader.addTypeLoader("IntTestRecord", [](conduit::Node const &asNode) {
    return std::make_unique<TestRecord<int>>(asNode);
  });

  Document document {documentAsNode, loader};
  auto &records = document.getRecords();
  ASSERT_EQ(1u, records.size());
  auto testRecord = dynamic_cast<TestRecord<int> const *>(records[0].get());
  ASSERT_NE(nullptr, testRecord);
  ASSERT_EQ(123, testRecord->getValue());
}

TEST(Document, create_fromNode_withRelationships)
{
  conduit::Node relationshipAsNode;
  relationshipAsNode["subject"] = "the subject";
  relationshipAsNode["object"] = "the object";
  relationshipAsNode["predicate"] = "is related to";

  conduit::Node relationshipsAsNodes;
  relationshipsAsNodes.append().set(relationshipAsNode);

  conduit::Node documentAsNode;
  documentAsNode[EXPECTED_RELATIONSHIPS_KEY] = relationshipsAsNodes;

  Document document {documentAsNode, RecordLoader {}};
  auto &relationships = document.getRelationships();
  ASSERT_EQ(1u, relationships.size());
  EXPECT_EQ("the subject", relationships[0].getSubject().getId());
  EXPECT_EQ(IDType::Global, relationships[0].getSubject().getType());
  EXPECT_EQ("the object", relationships[0].getObject().getId());
  EXPECT_EQ(IDType::Global, relationships[0].getObject().getType());
  EXPECT_EQ("is related to", relationships[0].getPredicate());
}

TEST(Document, toNode_empty)
{
  // A sina document should always have, at minimum, both records and
  // relationships as empty arrays.
  Document const document;
  conduit::Node asNode = document.toNode();
  EXPECT_TRUE(asNode[EXPECTED_RECORDS_KEY].dtype().is_list());
  EXPECT_EQ(0, asNode[EXPECTED_RECORDS_KEY].number_of_children());
  EXPECT_TRUE(asNode[EXPECTED_RELATIONSHIPS_KEY].dtype().is_list());
  EXPECT_EQ(0, asNode[EXPECTED_RELATIONSHIPS_KEY].number_of_children());
}

TEST(Document, toNode_records)
{
  Document document;
  std::string expectedIds[] = {"id 1", "id 2", "id 3"};
  std::string expectedValues[] = {"value 1", "value 2", "value 3"};

  auto numRecords = sizeof(expectedIds) / sizeof(expectedIds[0]);
  for(std::size_t i = 0; i < numRecords; ++i)
  {
    document.add(
      std::make_unique<TestRecord<std::string>>(expectedIds[i], TEST_RECORD_TYPE, expectedValues[i]));
  }

  auto asNode = document.toNode();

  auto record_nodes = asNode[EXPECTED_RECORDS_KEY];
  ASSERT_EQ(numRecords, record_nodes.number_of_children());
  for(auto i = 0; i < record_nodes.number_of_children(); ++i)
  {
    auto &actualNode = record_nodes[i];
    EXPECT_EQ(expectedIds[i], actualNode["id"].as_string());
    EXPECT_EQ(TEST_RECORD_TYPE, actualNode["type"].as_string());
    EXPECT_EQ(expectedValues[i], actualNode[TEST_RECORD_VALUE_KEY].as_string());
  }
}

TEST(Document, toNode_relationships)
{
  Document document;
  std::string expectedSubjects[] = {"subject 1", "subject 2"};
  std::string expectedObjects[] = {"object 1", "object 2"};
  std::string expectedPredicates[] = {"predicate 1", "predicate 2"};

  auto numRecords = sizeof(expectedSubjects) / sizeof(expectedSubjects[0]);
  for(unsigned long i = 0; i < numRecords; ++i)
  {
    document.add(Relationship {
      ID {expectedSubjects[i], IDType::Global},
      expectedPredicates[i],
      ID {expectedObjects[i], IDType::Global},
    });
  }

  auto asNode = document.toNode();

  auto relationship_nodes = asNode[EXPECTED_RELATIONSHIPS_KEY];
  ASSERT_EQ(numRecords, relationship_nodes.number_of_children());
  for(auto i = 0; i < relationship_nodes.number_of_children(); ++i)
  {
    auto &actualRelationship = relationship_nodes[i];
    EXPECT_EQ(expectedSubjects[i], actualRelationship["subject"].as_string());
    EXPECT_EQ(expectedObjects[i], actualRelationship["object"].as_string());
    EXPECT_EQ(expectedPredicates[i], actualRelationship["predicate"].as_string());
  }
}

TEST(Document, create_fromJson_roundtrip_json)
{
  std::string orig_json =
    "{\"records\": [{\"type\": \"test_rec\",\"id\": "
    "\"test\"}],\"relationships\": []}";
  axom::sina::Document myDocument = Document(orig_json, createRecordLoaderWithAllKnownTypes());
  EXPECT_EQ(0, myDocument.getRelationships().size());
  ASSERT_EQ(1, myDocument.getRecords().size());
  EXPECT_EQ("test_rec", myDocument.getRecords()[0]->getType());
  std::string returned_json1 = myDocument.toJson(0, 0, "", "");
  EXPECT_EQ(orig_json, returned_json1);
}

TEST(Document, create_fromJson_full_json)
{
  axom::sina::Document myDocument = Document(long_json, createRecordLoaderWithAllKnownTypes());
  EXPECT_EQ(2, myDocument.getRelationships().size());
  auto &records1 = myDocument.getRecords();
  EXPECT_EQ(4, records1.size());
}

TEST(Document, create_fromJson_value_check_json)
{
  axom::sina::Document myDocument = Document(SIMPLE_DOCUMENT, createRecordLoaderWithAllKnownTypes());
  EXPECT_EQ(2, myDocument.getRelationships().size());
  auto &records1 = myDocument.getRecords();
  EXPECT_EQ(1, records1.size());
  EXPECT_EQ(records1[0]->getType(), "run");
  auto &data1 = records1[0]->getData();
  EXPECT_EQ(data1.at("int").getScalar(), 500.0);
  std::vector<std::string> expected_string_vals = {"z", "o", "o"};
  EXPECT_EQ(data1.at("str/ings").getStringArray(), expected_string_vals);
  EXPECT_EQ(records1[0]->getFiles().count(File {"test/test.png"}), 1);
}

TEST(Document, saveDocument_json)
{
  axom::utilities::filesystem::TempFile tempFile("", ".json");

  Document document;
  document.add(std::make_unique<Record>(ID {"the id", IDType::Global}, "the type"));

  saveDocument(document, tempFile.getPath());

  conduit::Node readContents;
  readContents.parse(tempFile.getFileContents(), "json");

  ASSERT_TRUE(readContents[EXPECTED_RECORDS_KEY].dtype().is_list());
  EXPECT_EQ(1, readContents[EXPECTED_RECORDS_KEY].number_of_children());
  auto &readRecord = readContents[EXPECTED_RECORDS_KEY][0];
  EXPECT_EQ("the id", readRecord["id"].as_string());
  EXPECT_EQ("the type", readRecord["type"].as_string());
}

TEST(Document, load_specifiedRecordLoader)
{
  using RecordType = TestRecord<int>;
  auto originalRecord = std::make_unique<RecordType>("the ID", "my type", 123);
  Document originalDocument;
  originalDocument.add(std::move(originalRecord));

  axom::utilities::filesystem::TempFile tempfile("load_specifiedRecordLoader", ".json");
  tempfile.write(originalDocument.toNode().to_json());

  RecordLoader loader;
  loader.addTypeLoader("my type", [](conduit::Node const &asNode) {
    return std::make_unique<RecordType>(
      getRequiredString("id", asNode, "Test type"),
      getRequiredString("type", asNode, "Test type"),
      static_cast<int>(getRequiredField(TEST_RECORD_VALUE_KEY, asNode, "Test type").as_int64()));
  });
  Document loadedDocument = loadDocument(tempfile.getPath(), loader);
  ASSERT_EQ(1u, loadedDocument.getRecords().size());
  auto loadedRecord = dynamic_cast<RecordType const *>(loadedDocument.getRecords()[0].get());
  ASSERT_NE(nullptr, loadedRecord);
  EXPECT_EQ(123, loadedRecord->getValue());
}

TEST(Document, load_defaultRecordLoaders)
{
  auto originalRun =
    std::make_unique<axom::sina::Run>(ID {"the ID", IDType::Global}, "the app", "1.2.3", "jdoe");
  Document originalDocument;
  originalDocument.add(std::move(originalRun));

  axom::utilities::filesystem::TempFile tempfile("load_defaultRecordLoaders", ".json");
  tempfile.write(originalDocument.toNode().to_json());

  Document loadedDocument = loadDocument(tempfile.getPath());
  ASSERT_EQ(1u, loadedDocument.getRecords().size());
  auto loadedRun = dynamic_cast<axom::sina::Run const *>(loadedDocument.getRecords()[0].get());
  EXPECT_NE(nullptr, loadedRun);
}

TEST(Document, test_validate_append_typeclash)
{
  conduit::Node appendTo = parseJsonValue(SIMPLE_DOCUMENT);
  conduit::Node appendFrom = parseJsonValue(R"({"id": "rec_1", "type": "mismatch_type"})");
  conduit::Node msgNode = validateAppendDocument(appendTo, appendFrom, "", 1, 0);
  ASSERT_EQ(msgNode.number_of_children(), 1);
  // TODO: ugly string! Is there a good way to de-escape / from conduit?
  ASSERT_EQ(msgNode.child(0).to_string(),
            "\"Failed to append record 0: type mismatch \\\"run\\\"vs \\\"mismatch_type\\\"\"");
}

TEST(Document, test_validate_append_protocol_clash)
{
  conduit::Node appendTo = parseJsonValue(SIMPLE_DOCUMENT);
  // Note lack of id or type--could be a librarydata for all the method should care
  conduit::Node appendFrom = parseJsonValue(R"({"data": { "int": {"value": 20}}})");
  conduit::Node msgNode = validateAppendDocument(appendTo, appendFrom, "", 3, 0);
  ASSERT_EQ(msgNode.number_of_children(), 1);
  ASSERT_EQ(msgNode.child(0).to_string(),
            "\"Failed to append record 0 (protocol 3): conflicting data: int\"");
  appendFrom = parseJsonValue(R"({"user_defined": { "foo": "blarp"}})");
  ASSERT_EQ(validateAppendDocument(appendTo, appendFrom, "", 3, 0).child(0).to_string(),
            "\"Failed to append record 0 (protocol 3): conflicting user_defined: foo\"");
  appendFrom = parseJsonValue(R"({"files": { "test/test.png": {}}})");
  ASSERT_EQ(validateAppendDocument(appendTo, appendFrom, "", 3, 0).child(0).to_string(),
            "\"Failed to append record 0 (protocol 3): conflicting files: test/test.png\"");
}

TEST(Document, test_validate_append_missing_curve)
{
  conduit::Node appendTo = parseJsonValue(SIMPLE_DOCUMENT);
  // Only foo is being appended to, not bar. Bad!
  conduit::Node appendFrom =
    parseJsonValue(R"({"curve_sets": {"set_1": {"independent": {"0": {"value": [4, 5, 6]}}}}})");
  ASSERT_EQ(validateAppendDocument(appendTo, appendFrom, "", 1, 0).number_of_children(), 2);
  ASSERT_EQ(validateAppendDocument(appendTo, appendFrom, "", 1, 0).child(0).to_string(),
            "\"Failed to append record 0: did not append ALL or NO pre-existing curves (causing "
            "append element count mismatch)\"");
  ASSERT_EQ(validateAppendDocument(appendTo, appendFrom, "", 1, 0).child(1).to_string(),
            "\"Failed to append record 0's curve '0': count of appended elements would differ "
            "between series\"");
}

TEST(Document, test_validate_append_mismatched_curve_lengths)
{
  conduit::Node appendTo = parseJsonValue(SIMPLE_DOCUMENT);
  conduit::Node appendFrom = parseJsonValue(R"({"curve_sets": {"set_1": {
    "independent": {"0": {"value": [4, 5, 6]}, "1": {"value": [4, 5, 6]}},
    "dependent": {"0": {"value": [7, 8, 9]}, "1": {"value": [7, 8]}}}}})");
  conduit::Node msgNode = validateAppendDocument(appendTo, appendFrom, "", 1, 0);
  ASSERT_EQ(msgNode.number_of_children(), 1);
  ASSERT_EQ(msgNode.child(0).to_string(),
            "\"Failed to append record 0's curve '1': count of appended elements would differ "
            "between series\"");
}

TEST(Document, test_validate_append_mismatched_curve_lengths_newcurve)
{
  conduit::Node appendTo = parseJsonValue(SIMPLE_DOCUMENT);
  conduit::Node appendFrom =
    parseJsonValue(R"({"curve_sets": {"set_1": {"independent": {"new": {"value": [4]}}}}})");
  conduit::Node msgNode = validateAppendDocument(appendTo, appendFrom, "", 1, 0);
  ASSERT_EQ(msgNode.number_of_children(), 1);
  ASSERT_EQ(msgNode.child(0).to_string(),
            "\"Failed to append record 0's curve 'new': count of appended elements would differ "
            "between series\"");

  appendFrom = parseJsonValue(R"({"curve_sets": {"set_1": {
    "independent": {"new": {"value": [4, 5, 6]},
                    "0": {"value": [4, 5, 6]},
                    "1": {"value": [4, 5, 6]}},
    "dependent": {"0": {"value": [4, 5, 6]},
                  "1": {"value": [4, 5, 6]}}}}})");
  msgNode = validateAppendDocument(appendTo, appendFrom, "", 1, 0);
  ASSERT_EQ(msgNode.number_of_children(), 1);
  ASSERT_EQ(msgNode.child(0).to_string(),
            "\"Failed to append record 0's curve 'new': count of appended elements would differ "
            "between series\"");
}

TEST(Document, test_validate_append_issue_in_librarydata)
{
  conduit::Node appendTo = parseJsonValue(SIMPLE_DOCUMENT);
  conduit::Node appendFrom =
    parseJsonValue(R"({"library_data": {"my_lib": {"data": {"int": {"value": 50}}}}})");
  conduit::Node msgNode = validateAppendDocument(appendTo, appendFrom, "", 3, 0);
  ASSERT_EQ(msgNode.number_of_children(), 1);
  ASSERT_EQ(
    msgNode.child(0).to_string(),
    "\"Failed to append record 0/library_data/my_lib (protocol 3): conflicting data: int\"");
}

TEST(Document, test_validate_append_issue_nested_librarydata)
{
  conduit::Node appendTo = parseJsonValue(SIMPLE_DOCUMENT);
  conduit::Node appendFrom = parseJsonValue(
    R"({"library_data": {"my_lib": {"library_data": {"my_inner_lib": {"user_defined": {"foo/bar": "baz/qux"}}}}}})");
  conduit::Node msgNode = validateAppendDocument(appendTo, appendFrom, "", 3, 0);
  ASSERT_EQ(msgNode.number_of_children(), 1);
  ASSERT_EQ(msgNode.child(0).to_string(),
            "\"Failed to append record 0/library_data/my_lib/library_data/my_inner_lib (protocol "
            "3): conflicting user_defined: foo/bar\"");
}

TEST(Document, test_validate_append_valid)
{
  conduit::Node appendTo = parseJsonValue(SIMPLE_DOCUMENT);
  conduit::Node appendFrom = parseJsonValue(SIMPLE_DOCUMENT);
  ASSERT_EQ(validateAppendDocument(appendTo, appendFrom, "", 1, 0).number_of_children(), 0);
}

void doEveryErrorTest(
  const std::string &protocol,
  std::function<conduit::Node(const std::string &, const sina::Document &, int, bool)> appendDocumentFunc,
  bool skipValidation = false)
{
  std::string append_to_file = "test." + protocol;
  axom::sina::Document append_to_doc =
    Document(SIMPLE_DOCUMENT, createRecordLoaderWithAllKnownTypes());
  Protocol enum_protocol = (protocol == "hdf5") ? Protocol::HDF5 : Protocol::JSON;
  saveDocument(append_to_doc, append_to_file, enum_protocol);
  // Every error you could want and more
  conduit::Node appendFrom = parseJsonValue(R"(
    {"records": [{"id": "bar1", "type": "notarun", "data": { "int": {"value": 20}},
     "curve_sets": {"set_1": {"independent": {"0": {"value": [4, 5, 6]}}}},
     "library_data": {"my_lib": {"library_data": {"my_inner_lib": {"user_defined": {"foo/bar": "baz/qux"}}}}}}]})");
  axom::sina::Document new_doc = Document(appendFrom, createRecordLoaderWithAllKnownTypes());
  conduit::Node resultMsg = appendDocumentFunc(append_to_file, new_doc, 3, skipValidation);
  // Make sure no data changed
  conduit::Node root;
  conduit::relay::io::load(append_to_file, root);
  conduit::Node expect_root = parseJsonValue(SIMPLE_DOCUMENT);
  EXPECT_EQ(expect_root["records"].child(0)["id"].to_string(),
            root["records"].child(0)["id"].to_string());
  int expected = skipValidation ? 0 : 5;
  EXPECT_EQ(resultMsg.number_of_children(), expected);
}

TEST(Document, test_appendErrorCodepathsJSON) { doEveryErrorTest("json", appendDocumentToJson); }

TEST(Document, test_ignoreEveryError) { doEveryErrorTest("json", appendDocumentToJson, true); }

#ifdef AXOM_USE_HDF5
TEST(Document, test_appendErrorCodepathsHDF5) { doEveryErrorTest("hdf5", appendDocumentToHDF5); }
#endif

// Appending into an empty document
void doSimpleAppendTest(
  const std::string &protocol,
  std::function<conduit::Node(const std::string &, const sina::Document &, int, bool)> appendDocumentFunc)
{
  std::string empty_file = "test." + protocol;
  axom::sina::Document empty_doc =
    Document(R"({"records": [], "relationships": []})", createRecordLoaderWithAllKnownTypes());
  Protocol enum_protocol = (protocol == "hdf5") ? Protocol::HDF5 : Protocol::JSON;
  saveDocument(empty_doc, empty_file, enum_protocol);
  axom::sina::Document new_doc = Document(SIMPLE_DOCUMENT, createRecordLoaderWithAllKnownTypes());
  conduit::Node resultMsg = appendDocumentFunc(empty_file, new_doc, 3, true);  // skip validation
  EXPECT_EQ(resultMsg.number_of_children(), 0);
  conduit::Node root;
  conduit::relay::io::load(empty_file, root);
  conduit::Node expect_root = parseJsonValue(SIMPLE_DOCUMENT);
  EXPECT_EQ(expect_root["records"].child(0)["id"].to_string(),
            root["records"].child(0)["id"].to_string());
}

TEST(Document, test_simpleAppendDocumentToJson)
{
  doSimpleAppendTest("json", appendDocumentToJson);
}

#ifdef AXOM_USE_HDF5
TEST(Document, test_simpleAppendDocumentToHDF5)
{
  doSimpleAppendTest("hdf5", appendDocumentToHDF5);
}
#endif

// One unchanged, one merged
void doFullAppendTest(
  const std::string &protocol,
  std::function<conduit::Node(const std::string &, const sina::Document &, int, bool)> appendDocumentFunc)
{
  std::string filePath = "test." + protocol;
  sina::Document testDoc = Document(MULTI_REC_DOCUMENT, createRecordLoaderWithAllKnownTypes());
  Protocol enum_protocol = (protocol == "hdf5") ? Protocol::HDF5 : Protocol::JSON;
  saveDocument(testDoc, filePath, enum_protocol);

  axom::sina::Document new_doc = Document(SIMPLE_DOCUMENT, createRecordLoaderWithAllKnownTypes());
  conduit::Node resultMsg = appendDocumentFunc(filePath, new_doc, 1, false);
  EXPECT_EQ(resultMsg.number_of_children(), 0);

  conduit::Node root;
  conduit::relay::io::load(filePath, root);
  EXPECT_EQ(root["records"].number_of_children(), 2);
  conduit::Node appendTo = parseJsonValue(MULTI_REC_DOCUMENT);
  conduit::Node appendFrom = parseJsonValue(SIMPLE_DOCUMENT);
  // One record should be unchanged, but loading it into a document means we
  // don't guarantee data order except where important (curve set order). Spot check shared val.
  std::vector<double> expected = {1.0, 2.0};
  const conduit::Node &rec1 = root["records"].child(1);
  auto actual = node_to_double_vector(rec1["curve_sets"]["set_1"]["dependent"]["0"]["value"]);
  EXPECT_EQ(expected, actual);
  // The hard one, now a blend of the prior and new record
  conduit::Node rec0;
  restoreSlashes(root["records"].child(0), rec0);

  EXPECT_EQ(rec0["type"].as_string(), "run");
  EXPECT_EQ(rec0["data"]["string"]["value"].as_string(), "goodbye!");
  EXPECT_EQ(rec0["data"].child("str/ings")["value"][0].as_string(), "z");
  EXPECT_EQ(rec0["data"]["int"]["value"].as_float64(), 500);
  EXPECT_EQ(rec0["data"]["string2"]["value"].as_string(), "unchanged!");
  EXPECT_EQ(rec0["user_defined"]["hello"].as_string(), "and");
  EXPECT_EQ(rec0["user_defined"]["foo"].as_string(), "bar");
  expected = {11, 12, 13, 1, 2, 3};
  actual = node_to_double_vector(rec0["curve_sets"]["set_1"]["dependent"]["0"]["value"]);
  EXPECT_EQ(expected, actual);
  expected = {10, 20, 30, -1, -2, -3};
  actual = node_to_double_vector(rec0["curve_sets"]["set_1"]["dependent"]["1"]["value"]);
  EXPECT_EQ(expected, actual);
  expected = {14, 15, 16, 4, 5, 6};
  actual = node_to_double_vector(rec0["curve_sets"]["set_1"]["independent"]["0"]["value"]);
  EXPECT_EQ(expected, actual);
  expected = {7, 8, 9, -4, -5, -6};
  actual = node_to_double_vector(rec0["curve_sets"]["set_1"]["independent"]["1"]["value"]);
  EXPECT_EQ(expected, actual);
  EXPECT_EQ(expected, actual);
  // Relationships are easy, we just need the union (no duplicates)
  const conduit::Node &rels = root["relationships"];
  EXPECT_EQ(rels.number_of_children(), 3);
  EXPECT_EQ(rels.to_string(),
            "\n- \n  predicate: \"knows\"\n  local_subject: \"bar1\"\n  object: \"something\"\n- "
            "\n  predicate: \"sees\"\n  local_subject: \"bar2\"\n  object: \"bar1\"\n- \n  "
            "predicate: \"this_is\"\n  local_subject: \"bar1\"\n  object: \"unique\"\n");
  return;
}

TEST(Document, test_appendDocumentToJson) { doFullAppendTest("json", appendDocumentToJson); }

#ifdef AXOM_USE_HDF5
TEST(Document, test_appendDocumentToHDF5) { doFullAppendTest("hdf5", appendDocumentToHDF5); }
#endif

// Making sure we respect curve order (Records are in charge of ordering their curves, not documents)
void doAppendOrderedCurveTest(
  const std::string &protocol,
  std::function<conduit::Node(const std::string &, const sina::Document &, int, bool)> appendDocumentFunc)
{
  std::string curvedump_file = "test_curve." + protocol;
  axom::sina::Document ordered_curves =
    Document(CURVE_ORDERED_DOCUMENT, createRecordLoaderWithAllKnownTypes());
  axom::sina::Document additional_curve =
    Document(CURVE_ORDERED_DOCUMENT_APPEND, createRecordLoaderWithAllKnownTypes());
  Protocol enum_protocol = (protocol == "hdf5") ? Protocol::HDF5 : Protocol::JSON;
  saveDocument(ordered_curves, curvedump_file, enum_protocol);
  conduit::Node resultMsg = appendDocumentFunc(curvedump_file, additional_curve, 3, false);
  EXPECT_EQ(resultMsg.number_of_children(), 0);
  conduit::Node root;
  conduit::relay::io::load(curvedump_file, root);
  conduit::Node expected_dependents =
    parseJsonValue(CURVE_ORDERED_DOCUMENT)["records"].child(0)["curve_sets"]["set_1"]["dependent"];
  conduit::Node actual_dependents = root["records"].child(0)["curve_sets"]["set_1"]["dependent"];
  auto curveIter = actual_dependents.children();
  for(int i = 0; i < expected_dependents.number_of_children(); i++)
  {
    EXPECT_EQ(actual_dependents.child(i).name(), expected_dependents.child(i).name());
    EXPECT_EQ(node_to_double_vector(actual_dependents.child(i)["value"]),
              node_to_double_vector(expected_dependents.child(i)["value"]));
  }
  EXPECT_EQ(actual_dependents.child(3).name(), "am_i_fourth");
  std::vector<double> expected = {-100.25, -110.18, -120.4};
  EXPECT_EQ(node_to_double_vector(actual_dependents.child(3)["value"]), expected);
}

TEST(Document, test_appendOrderedCurvesToJson)
{
  doAppendOrderedCurveTest("json", appendDocumentToJson);
}

#ifdef AXOM_USE_HDF5
TEST(Document, test_appendOrderedCurvesToHDF5)
{
  doAppendOrderedCurveTest("hdf5", appendDocumentToHDF5);
}

TEST(Document, create_fromJson_roundtrip_hdf5)
{
  std::string orig_json =
    "{\"records\": [{\"type\": \"test_rec\",\"id\": "
    "\"test\"}],\"relationships\": []}";
  axom::sina::Document myDocument = Document(orig_json, createRecordLoaderWithAllKnownTypes());
  saveDocument(myDocument, "round_json.hdf5", Protocol::HDF5);
  Document loadedDocument = loadDocument("round_json.hdf5", Protocol::HDF5);
  EXPECT_EQ(0, loadedDocument.getRelationships().size());
  ASSERT_EQ(1, loadedDocument.getRecords().size());
  EXPECT_EQ("test_rec", loadedDocument.getRecords()[0]->getType());
  std::string returned_json2 = loadedDocument.toJson(0, 0, "", "");
  EXPECT_EQ(orig_json, returned_json2);
}

TEST(Document, create_fromJson_full_hdf5)
{
  axom::sina::Document myDocument = Document(long_json, createRecordLoaderWithAllKnownTypes());
  saveDocument(myDocument, "long_json.hdf5", Protocol::HDF5);
  Document loadedDocument = loadDocument("long_json.hdf5", Protocol::HDF5);
  EXPECT_EQ(2, loadedDocument.getRelationships().size());
  auto &records2 = loadedDocument.getRecords();
  EXPECT_EQ(4, records2.size());
}

TEST(Document, create_fromJson_value_check_hdf5)
{
  axom::sina::Document myDocument = Document(SIMPLE_DOCUMENT, createRecordLoaderWithAllKnownTypes());
  std::vector<std::string> expected_string_vals = {"z", "o", "o"};
  saveDocument(myDocument, "data_json.hdf5", Protocol::HDF5);
  Document loadedDocument = loadDocument("data_json.hdf5", Protocol::HDF5);
  EXPECT_EQ(2, loadedDocument.getRelationships().size());
  auto &records2 = loadedDocument.getRecords();
  EXPECT_EQ(1, records2.size());
  EXPECT_EQ(records2[0]->getType(), "run");
  auto &data2 = records2[0]->getData();
  EXPECT_EQ(data2.at("int").getScalar(), 500.0);
  EXPECT_EQ(data2.at("str/ings").getStringArray(), expected_string_vals);
  EXPECT_EQ(records2[0]->getFiles().count(File {"test/test.png"}), 1);
}

TEST(Document, saveDocument_hdf5)
{
  axom::utilities::filesystem::TempFile tempFile("saveDocument", ".hdf5");

  Document document;
  document.add(std::make_unique<Record>(ID {"the id", IDType::Global}, "the type"));

  saveDocument(document, tempFile.getPath(), Protocol::HDF5);

  conduit::Node readContents;
  conduit::relay::io::load(tempFile.getPath(), "hdf5", readContents);

  ASSERT_TRUE(readContents[EXPECTED_RECORDS_KEY].dtype().is_list());
  EXPECT_EQ(1, readContents[EXPECTED_RECORDS_KEY].number_of_children());
  auto &readRecord = readContents[EXPECTED_RECORDS_KEY][0];
  EXPECT_EQ("the id", readRecord["id"].as_string());
  EXPECT_EQ("the type", readRecord["type"].as_string());
}

#endif
}  // namespace
}  // namespace testing
}  // namespace sina
}  // namespace axom
