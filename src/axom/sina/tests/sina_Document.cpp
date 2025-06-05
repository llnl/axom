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


// Simple document to append into
std::string SIMPLE_DOCUMENT = R"(
{
  "records": [
    {
      "type": "run",
      "application": "test",
      "id": "rec_1",
      "data": {
        "int": {
          "value": 500,
          "units": "miles"
        },
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
          "independent": {
            "foo": {"values": [1, 2, 3]}
          },
          "dependent": {
            "bar": {"values": [4, 5, 6]}
          }
        }
      }
    }
  ]
}
)";

std::string data_json = R"(
{
  "records": [
    {
      "type": "run",
      "application": "test",
      "id": "test_1",
      "data": {
        "int": {
          "value": 500,
          "units": "miles"
        },
        "str/ings": {
          "value": ["z", "o", "o"]
        }
      },
      "files": {
        "test/test.png": {}
      }
    }
  ]
}
)";

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

// String of a new Doc to use for appending
std::string new_data =
  "{\"records\": ["
  "   {"
  "       \"type\": \"foo1\","
  "       \"id\": \"bar1\","
  "       \"data\": {\"scal1\": {\"value\": \"goodbye!\"}, \"scal2\": {\"value\": \"goodbye!\"}},"
  "       \"user_defined\":{\"hello\": \"goodbye\"},"
  "       \"curve_sets\": {"
  "           \"1\": {"
  "               \"dependent\": {"
  "                   \"0\": {\"value\": [1, 1]},"
  "                   \"1\": {\"value\": [8, 9]}"
  "               },"
  "               \"independent\": {"
  "                   \"0\": {\"value\": [4, 5]},"
  "                   \"1\": {\"value\": [7, 9]}"
  "               }"
  "           },"
  "           \"2\": {"
  "               \"dependent\": {"
  "                   \"0\": {\"value\": [1, 1]},"
  "                   \"1\": {\"value\": [10, 10]}"
  "               },"
  "               \"independent\": {"
  "                   \"0\": {\"value\": [4, 4]},"
  "                   \"1\": {\"value\": [7, 7]}"
  "               }"
  "           }"
  "       }"
  "   },"
  "   {"
  "       \"type\": \"foo2\","
  "       \"id\": \"bar3\","
  "       \"curve_sets\": {"
  "           \"1\": {"
  "               \"dependent\": {"
  "                   \"0\": {\"value\": [1, 1]},"
  "                   \"1\": {\"value\": [10, 10]}"
  "               },"
  "               \"independent\": {"
  "                   \"0\": {\"value\": [4, 4]},"
  "                   \"1\": {\"value\": [7, 7]}"
  "               }"
  "           }"
  "       }"
  "   }"
  "]}";

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

// We need to call these asserts multiple times for the append failures
void ValidateRecordsUnchanged(const conduit::Node &root, const std::string &context = "")
{
  std::vector<double> expected, actual;
  SCOPED_TRACE("Validation context: " + context);

  ASSERT_EQ(root["records"].number_of_children(), 2);

  const conduit::Node &rec0 = root["records"].child(0);
  ASSERT_EQ(rec0["id"].as_string(), "bar1");
  ASSERT_EQ(rec0["type"].as_string(), "foo1");
  ASSERT_EQ(rec0["data"]["scal1"]["value"].as_string(), "hello!");
  ASSERT_FALSE(rec0["data"].has_child("scal2"));
  ASSERT_EQ(rec0["user_defined"]["hello"].as_string(), "and");
  expected = {1.0, 2.0};
  actual = node_to_double_vector(rec0["curve_sets"]["1"]["dependent"]["0"]["value"]);
  ASSERT_EQ(actual, expected);
  expected = {10.0, 20.0};
  actual = node_to_double_vector(rec0["curve_sets"]["1"]["dependent"]["1"]["value"]);
  ASSERT_EQ(actual, expected);
  expected = {4.0, 5.0};
  actual = node_to_double_vector(rec0["curve_sets"]["1"]["independent"]["0"]["value"]);
  ASSERT_EQ(actual, expected);
  expected = {7.0, 8.0};
  actual = node_to_double_vector(rec0["curve_sets"]["1"]["independent"]["1"]["value"]);
  ASSERT_EQ(actual, expected);
  ASSERT_FALSE(rec0["curve_sets"].has_child("2"));

  const conduit::Node &rec1 = root["records"].child(1);
  ASSERT_EQ(rec1["id"].as_string(), "bar2");
  ASSERT_EQ(rec1["type"].as_string(), "foo2");
  expected = {1.0, 2.0};
  actual = node_to_double_vector(rec1["curve_sets"]["1"]["dependent"]["0"]["value"]);
  ASSERT_EQ(actual, expected);
  expected = {10.0, 20.0};
  actual = node_to_double_vector(rec1["curve_sets"]["1"]["dependent"]["1"]["value"]);
  ASSERT_EQ(actual, expected);
  expected = {4.0, 5.0};
  actual = node_to_double_vector(rec1["curve_sets"]["1"]["independent"]["0"]["value"]);
  ASSERT_EQ(actual, expected);
  expected = {7.0, 8.0};
  actual = node_to_double_vector(rec1["curve_sets"]["1"]["independent"]["1"]["value"]);
  ASSERT_EQ(actual, expected);
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

/**
 * Instances of this class acquire a temporary file name when created
 * and delete the file when destructed.
 *
 * NOTE: This class uses unsafe methods and should only be used for testing
 * purposes. DO NOT move it to the main code!!!!
 */
class NamedTempFile
{
public:
  NamedTempFile();

  // As a resource-holding class, we don't want this to be copyable
  // (or movable since there is no reason to return it from a function)
  NamedTempFile(NamedTempFile const &) = delete;

  NamedTempFile(NamedTempFile &&) = delete;

  NamedTempFile &operator=(NamedTempFile const &) = delete;

  NamedTempFile &operator=(NamedTempFile &&) = delete;

  ~NamedTempFile();

  std::string const &getName() const { return fileName; }

private:
  std::string fileName;
};

NamedTempFile::NamedTempFile()
{
  std::vector<char> tmpFileName;
  tmpFileName.resize(L_tmpnam);
  // tmpnam is not the best way to do this, but it is standard and this is
  // only a test.
  if(!std::tmpnam(tmpFileName.data()))
  {
    throw std::ios::failure {"Could not get temporary file"};
  }
  fileName = tmpFileName.data();
}

NamedTempFile::~NamedTempFile() { axom::utilities::filesystem::removeFile(fileName.data()); }

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
  axom::sina::Document myDocument = Document(data_json, createRecordLoaderWithAllKnownTypes());
  EXPECT_EQ(0, myDocument.getRelationships().size());
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
  NamedTempFile tmpFile;

  // First, write some random stuff to the temp file to make sure it is
  // overwritten.
  {
    std::ofstream fout {tmpFile.getName()};
    fout << "Initial contents";
  }

  Document document;
  document.add(std::make_unique<Record>(ID {"the id", IDType::Global}, "the type"));

  saveDocument(document, tmpFile.getName() + ".json");

  conduit::Node readContents;
  {
    std::ifstream fin {tmpFile.getName() + ".json"};
    std::stringstream f_buf;
    f_buf << fin.rdbuf();
    readContents.parse(f_buf.str(), "json");
  }

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

  NamedTempFile file;
  {
    std::ofstream fout {file.getName()};
    fout << originalDocument.toNode().to_json();
  }

  RecordLoader loader;
  loader.addTypeLoader("my type", [](conduit::Node const &asNode) {
    return std::make_unique<RecordType>(
      getRequiredString("id", asNode, "Test type"),
      getRequiredString("type", asNode, "Test type"),
      static_cast<int>(getRequiredField(TEST_RECORD_VALUE_KEY, asNode, "Test type").as_int64()));
  });
  Document loadedDocument = loadDocument(file.getName(), loader);
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

  NamedTempFile file;
  {
    std::ofstream fout {file.getName()};
    fout << originalDocument.toNode().to_json();
  }

  Document loadedDocument = loadDocument(file.getName());
  ASSERT_EQ(1u, loadedDocument.getRecords().size());
  auto loadedRun = dynamic_cast<axom::sina::Run const *>(loadedDocument.getRecords()[0].get());
  EXPECT_NE(nullptr, loadedRun);
}

// Get us a file handle for a basic document.
conduit::relay::io::IOHandle setup_basic_file(){
  std::string jsonFilePath = "simple_sina_doc.json";
  std::ofstream testFile(jsonFilePath);
  testFile << SIMPLE_DOCUMENT;
  testFile.close();
  conduit::relay::io::IOHandle simpleDoc;
  simpleDoc.open(jsonFilePath);
  return simpleDoc;
}

//Figure out where the librarydata is going
conduit::relay::io::IOHandle do_hideous_print_append(const std::string &json_value){
  std::string jsonFilePath = "simple_sina_doc.json";
  std::ofstream testFile(jsonFilePath);
  testFile << SIMPLE_DOCUMENT;
  testFile.close();
  conduit::relay::io::IOHandle simpleDoc;
  simpleDoc.open(jsonFilePath);
  return simpleDoc;
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
  conduit::Node appendFrom = parseJsonValue(R"({"id": "rec_1", "type": "run", "data": { "int": {"value": 20}}})");
  conduit::Node msgNode = validateAppendDocument(appendTo, appendFrom, "", 3, 0);
  ASSERT_EQ(msgNode.number_of_children(), 1);
  ASSERT_EQ(msgNode.child(0).to_string(),
    "\"Failed to append record 0 (protocol 3): conflicting data: int\"");
  appendFrom = parseJsonValue(R"({"user_defined": { "foo": "blarp"}})");
  ASSERT_EQ(validateAppendDocument(appendTo, appendFrom, "", 3, 0).child(0).to_string(),
           "\"Failed to append record 0 (protocol 3): conflicting user_defined: foo\""
  );
  appendFrom = parseJsonValue(R"({"files": { "test/test.png": {}}})");
  ASSERT_EQ(validateAppendDocument(appendTo, appendFrom, "", 3, 0).child(0).to_string(),
            "\"Failed to append record 0 (protocol 3): conflicting files: test/test.png\""
  );
}

TEST(Document, test_validate_append_missing_curve)
{
  conduit::Node appendTo = parseJsonValue(SIMPLE_DOCUMENT);
  // Only foo is being appended to, not bar. Bad!
  conduit::Node appendFrom = parseJsonValue(R"({"curve_sets": {"set_1": {"independent": {"foo": {"values": [4, 5, 6]}}}}})");
  ASSERT_EQ(validateAppendDocument(appendTo, appendFrom, "", 1, 0).number_of_children(), 1);
}

TEST(Document, test_validate_append_mismatched_curve_lengths)
{
  conduit::Node appendTo = parseJsonValue(SIMPLE_DOCUMENT);
  conduit::Node appendFrom = parseJsonValue(R"({"curve_sets": {"set_1": {
    "independent": {"foo": {"values": [4, 5, 6]}},
    "dependent": {"bar": {"values": [7, 8]}}}}})");
  ASSERT_EQ(validateAppendDocument(appendTo, appendFrom, "", 1, 0).number_of_children(), 1);
}

TEST(Document, test_validate_append_mismatched_curve_lengths_newcurve)
{
  conduit::Node appendTo = parseJsonValue(SIMPLE_DOCUMENT);
  conduit::Node appendFrom = parseJsonValue(R"({"curve_sets": {"set_1": {"independent": {"new": {"values": [4]}}}}})");
  ASSERT_EQ(validateAppendDocument(appendTo, appendFrom, "", 1, 0).number_of_children(), 0);
  appendFrom = parseJsonValue(R"({"curve_sets": {"set_1": {
    "independent": {"new": {"values": [4, 5, 6]}},
    "dependent": {"bar": {"values": [4, 5, 6]}}}}})");
  ASSERT_EQ(validateAppendDocument(appendTo, appendFrom, "", 1, 0).number_of_children(), 1);
}

TEST(Document, test_validate_append_valid)
{
  conduit::Node appendTo = parseJsonValue(SIMPLE_DOCUMENT);
  conduit::Node appendFrom = parseJsonValue(SIMPLE_DOCUMENT);
  ASSERT_EQ(validateAppendDocument(appendTo, appendFrom, "", 1, 0).number_of_children(), 0);
}

TEST(Document, test_basic_append)
{
  std::string jsonFilePath = "find-me.json";
  std::ofstream testFile(jsonFilePath);
  // Create a near-empty document
  testFile << R"(
  {
    "records": [{
      "type": "foo",
      "id": "i_am_empty"
    }]
  })";
  testFile.close();
  axom::sina::Document new_doc = Document(SIMPLE_DOCUMENT, createRecordLoaderWithAllKnownTypes());
  /*conduit::Node resultMsg = appendDocumentToJson(jsonFilePath, new_doc, 1);
  ASSERT_EQ(resultMsg.number_of_children(), 0);
  conduit::Node append_root;
  append_root.load(jsonFilePath, "json");
  conduit::Node expect_root;
  append_root.load(SIMPLE_DOCUMENT, "json");
  ASSERT_EQ(expect_root.child(0)["id"].to_string(), append_root.child(1)["id"].to_string());*/
}

TEST(Document, test_appendDocumentToJson)
{
  std::string jsonFilePath = "test.json";

  std::ofstream testFile(jsonFilePath);
  testFile << R"({
      "records": [
          {
              "data": {"scal1": {"value": "hello!"}},
              "type": "foo1",
              "id": "bar1",
              "user_defined":{"hello": "and"},
              "curve_sets": {
                  "1": {
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
          },
          {
              "type": "foo2",
              "id": "bar2",
              "curve_sets": {
                  "1": {
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
      ]
  })";
  testFile.close();

  axom::sina::Document new_doc = Document(new_data, createRecordLoaderWithAllKnownTypes());

  conduit::Node resultMsg = appendDocumentToJson(jsonFilePath, new_doc, 1);
  ASSERT_EQ(resultMsg.number_of_children(), 0);

  conduit::Node root;
  root.load(jsonFilePath, "json");
  std::vector<double> expected, actual;

  ASSERT_EQ(root["records"].number_of_children(), 3);

  // ----- Record 0 -----
  const conduit::Node &rec0 = root["records"].child(0);
  ASSERT_EQ(rec0["id"].as_string(), "bar1");
  ASSERT_EQ(rec0["type"].as_string(), "foo1");
  ASSERT_EQ(rec0["data"]["scal1"]["value"].as_string(), "goodbye!");
  ASSERT_EQ(rec0["data"]["scal2"]["value"].as_string(), "goodbye!");
  ASSERT_EQ(rec0["user_defined"]["hello"].as_string(), "goodbye");
  expected = {1, 2, 1.0, 1.0};
  actual = node_to_double_vector(rec0["curve_sets"]["1"]["dependent"]["0"]["value"]);
  ASSERT_EQ(actual, expected);
  expected = {10, 20, 8.0, 9.0};
  actual = node_to_double_vector(rec0["curve_sets"]["1"]["dependent"]["1"]["value"]);
  ASSERT_EQ(actual, expected);
  expected = {4, 5, 4.0, 5.0};
  actual = node_to_double_vector(rec0["curve_sets"]["1"]["independent"]["0"]["value"]);
  ASSERT_EQ(actual, expected);
  expected = {7, 8, 7.0, 9.0};
  actual = node_to_double_vector(rec0["curve_sets"]["1"]["independent"]["1"]["value"]);
  ASSERT_EQ(actual, expected);
  expected = {1.0, 1.0};
  actual = node_to_double_vector(rec0["curve_sets"]["2"]["dependent"]["0"]["value"]);
  ASSERT_EQ(actual, expected);
  expected = {10.0, 10.0};
  actual = node_to_double_vector(rec0["curve_sets"]["2"]["dependent"]["1"]["value"]);
  ASSERT_EQ(actual, expected);
  expected = {4.0, 4.0};
  actual = node_to_double_vector(rec0["curve_sets"]["2"]["independent"]["0"]["value"]);
  ASSERT_EQ(actual, expected);
  expected = {7.0, 7.0};
  actual = node_to_double_vector(rec0["curve_sets"]["2"]["independent"]["1"]["value"]);
  ASSERT_EQ(actual, expected);

  // ----- Record 1 -----
  const conduit::Node &rec1 = root["records"].child(1);
  ASSERT_EQ(rec1["id"].as_string(), "bar2");
  ASSERT_EQ(rec1["type"].as_string(), "foo2");

  expected = {1, 2};
  actual = node_to_double_vector(rec1["curve_sets"]["1"]["dependent"]["0"]["value"]);
  ASSERT_EQ(actual, expected);
  expected = {10, 20};
  actual = node_to_double_vector(rec1["curve_sets"]["1"]["dependent"]["1"]["value"]);
  ASSERT_EQ(actual, expected);
  expected = {4, 5};
  actual = node_to_double_vector(rec1["curve_sets"]["1"]["independent"]["0"]["value"]);
  ASSERT_EQ(actual, expected);
  expected = {7, 8};
  actual = node_to_double_vector(rec1["curve_sets"]["1"]["independent"]["1"]["value"]);
  ASSERT_EQ(actual, expected);

  // ----- Record 2 -----
  const conduit::Node &rec2 = root["records"].child(2);
  ASSERT_EQ(rec2["id"].as_string(), "bar3");
  ASSERT_EQ(rec2["type"].as_string(), "foo2");

  expected = {1.0, 1.0};
  actual = node_to_double_vector(rec2["curve_sets"]["1"]["dependent"]["0"]["value"]);
  ASSERT_EQ(actual, expected);
  expected = {10.0, 10.0};
  actual = node_to_double_vector(rec2["curve_sets"]["1"]["dependent"]["1"]["value"]);
  ASSERT_EQ(actual, expected);
  expected = {4.0, 4.0};
  actual = node_to_double_vector(rec2["curve_sets"]["1"]["independent"]["0"]["value"]);
  ASSERT_EQ(actual, expected);
  expected = {7.0, 7.0};
  actual = node_to_double_vector(rec2["curve_sets"]["1"]["independent"]["1"]["value"]);
  ASSERT_EQ(actual, expected);
}

TEST(Document, test_appendDocumentToJson_failures)
{
  // We have a lot of potential failures
  // 1. Fail validation due to an appended-to curve not being uniform
  // 2. Fail validation due to a new curve not being uniform
  // 3. Fail because Dependents won't be uniform with Independents
  // 4. Fail because Protocol 3 Duplicate Data
  // 5. Fail because Protocol 3 Duplicate UDC

  std::string jsonFilePath = "test_failure.json";

  std::ofstream testFile(jsonFilePath);
  testFile << R"({
    "records": [
        {
            "data": {"scal1": {"value": "hello!"}},
            "type": "foo1",
            "id": "bar1",
            "user_defined":{"hello": "and"},
            "curve_sets": {
                "1": {
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
        },
        {
            "type": "foo2",
            "id": "bar2",
            "curve_sets": {
                "1": {
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
    ]
  })";
  testFile.close();

  // ---- Failure 1 ----
  SCOPED_TRACE("Attempting Case 1");
  std::string f1_data = ([](std::string s) {
    conduit::Node tmp;
    tmp.parse(s, "json");
    conduit::Node &valNode = tmp["records"].child(0)["curve_sets"]["1"]["dependent"]["0"]["value"];
    valNode.reset();
    std::vector<double> new_val = {1, 1, 1};
    valNode.set(new_val);
    return tmp.to_json();
  })(new_data);

  axom::sina::Document new_doc = Document(f1_data, createRecordLoaderWithAllKnownTypes());
  conduit::Node resultMsg = appendDocumentToJson(jsonFilePath, new_doc, 1);
  ASSERT_EQ(resultMsg.number_of_children(), 4);
  conduit::Node j;
  j.load(jsonFilePath, "json");
  ValidateRecordsUnchanged(j, "Case 1: Checking if File Was Changed");

  // ---- Failure 2 ----
  SCOPED_TRACE("Attempting Case 2");
  std::string f2_data = ([](std::string s) {
    conduit::Node tmp;
    tmp.parse(s, "json");
    conduit::Node &curve_set1 = tmp["records"].child(0)["curve_sets"]["1"];
    if(curve_set1.has_child("dependent")) curve_set1["dependent"].reset();
    conduit::Node &newDep = curve_set1["dependent"]["2"];
    newDep.reset();
    std::vector<double> newVal = {1};
    newDep["value"].set(newVal);
    if(curve_set1.has_child("independent")) curve_set1["independent"].reset();
    return tmp.to_json();
  })(new_data);

  /*new_doc = Document(f2_data, createRecordLoaderWithAllKnownTypes());
  result = appendDocumentToJson(jsonFilePath, new_doc, 1);
  ASSERT_FALSE(result);
  j.load(jsonFilePath, "json");
  ValidateRecordsUnchanged(j, "Case 2: Checking if File Was Changed");

  // ---- Failure 3 ----
  SCOPED_TRACE("Attempting Case 3");
  std::string f3_data = ([](std::string s) {
    conduit::Node tmp;
    tmp.parse(s, "json");
    conduit::Node &valNode =
      tmp["records"].child(0)["curve_sets"]["1"]["independent"]["1"]["value"];
    valNode.reset();
    std::vector<double> new_val = {1, 1, 1};
    valNode.set(new_val);
    return tmp.to_json();
  })(new_data);

  new_doc = Document(f3_data, createRecordLoaderWithAllKnownTypes());
  result = appendDocumentToJson(jsonFilePath, new_doc, 1);
  ASSERT_FALSE(result);
  j.load(jsonFilePath, "json");
  ValidateRecordsUnchanged(j, "Case 3: Checking if File Was Changed");

  // ---- Failure 4 ----
  SCOPED_TRACE("Attempting Case 4");
  new_doc = Document(new_data, createRecordLoaderWithAllKnownTypes());
  result = appendDocumentToJson(jsonFilePath, new_doc, 3);
  ASSERT_FALSE(result);
  j.load(jsonFilePath, "json");
  ValidateRecordsUnchanged(j, "Case 4: Checking if File Was Changed");

  // ---- Failure 5 ----
  SCOPED_TRACE("Attempting Case 5");
  new_doc = Document(new_data, createRecordLoaderWithAllKnownTypes());
  result = appendDocumentToJson(jsonFilePath, new_doc, 1);
  ASSERT_FALSE(result);
  j.load(jsonFilePath, "json");
  ValidateRecordsUnchanged(j, "Case 5: Checking if File Was Changed");

  // ---- Assert Error Messages ----
  std::cerr.rdbuf(origBuffer);
  std::string expected =
  ASSERT_EQ(msgNode.child(0),
           "Error validating dependent: Record bar1, Curve Set 1, Curve 1 size mismatch "
           "(expected 5, got 4).");
  ASSERT_EQ(msgNode.child(1),
    "Error validating dependent: Record bar1, Curve Set 1, Curve 2 size mismatch "
    "(expected 2, got 1).");
  ASSERT_EQ(msgNode.child(2),
    "Error validating independent: Record bar1, Curve Set 1, Curve 1 size mismatch "
    "(expected 4, got 5).");
  ASSERT_EQ(msgNode.child(3),"Found a duplicate data entry, protocol 3 dictates "
    "append cancellation. Found duplicate UDC, protocol 3 dictates append cancellation.";)*/
}

#ifdef AXOM_USE_HDF5
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
  axom::sina::Document myDocument = Document(data_json, createRecordLoaderWithAllKnownTypes());
  std::vector<std::string> expected_string_vals = {"z", "o", "o"};
  saveDocument(myDocument, "data_json.hdf5", Protocol::HDF5);
  Document loadedDocument = loadDocument("data_json.hdf5", Protocol::HDF5);
  EXPECT_EQ(0, loadedDocument.getRelationships().size());
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
  NamedTempFile tmpFile;

  // First, write some random stuff to the temp file to make sure it is
  // overwritten.
  {
    std::ofstream fout {tmpFile.getName()};
    fout << "Initial contents";
  }

  Document document;
  document.add(std::make_unique<Record>(ID {"the id", IDType::Global}, "the type"));

  saveDocument(document, tmpFile.getName() + ".hdf5", Protocol::HDF5);

  conduit::Node readContents;
  conduit::relay::io::load(tmpFile.getName() + ".hdf5", "hdf5", readContents);

  ASSERT_TRUE(readContents[EXPECTED_RECORDS_KEY].dtype().is_list());
  EXPECT_EQ(1, readContents[EXPECTED_RECORDS_KEY].number_of_children());
  auto &readRecord = readContents[EXPECTED_RECORDS_KEY][0];
  EXPECT_EQ("the id", readRecord["id"].as_string());
  EXPECT_EQ("the type", readRecord["type"].as_string());
}

TEST(Document, test_append_to_hdf5)
{
  std::string hdf5FilePath = "test.hdf5";
  conduit::Node initialData;
  initialData["records/0/data/scal1/value"].set("hello!");
  initialData["records/0/type"].set("foo1");
  initialData["records/0/id"].set("bar1");
  initialData["records/0/user_defined/hello"].set("and");
  initialData["records/0/curve_sets/1/dependent/0/value"].set(std::vector<double> {1.0, 2.0});
  initialData["records/0/curve_sets/1/dependent/1/value"].set(std::vector<double> {10.0, 20.0});
  initialData["records/0/curve_sets/1/independent/0/value"].set(std::vector<double> {4.0, 5.0});
  initialData["records/0/curve_sets/1/independent/1/value"].set(std::vector<double> {7.0, 8.0});
  initialData["records/1/type"].set("foo2");
  initialData["records/1/id"].set("bar2");
  initialData["records/1/curve_sets/1/dependent/0/value"].set(std::vector<double> {1.0, 2.0});
  initialData["records/1/curve_sets/1/dependent/1/value"].set(std::vector<double> {10.0, 20.0});
  initialData["records/1/curve_sets/1/independent/0/value"].set(std::vector<double> {4.0, 5.0});
  initialData["records/1/curve_sets/1/independent/1/value"].set(std::vector<double> {7.0, 8.0});
  conduit::relay::io::hdf5_write(initialData, hdf5FilePath);

  axom::sina::Document new_doc = Document(new_data, createRecordLoaderWithAllKnownTypes());
  conduit::Node msgNode = append_to_hdf5(hdf5FilePath, new_doc, 1);
  ASSERT_EQ(msgNode.number_of_children(), 0);
  conduit::Node root;
  conduit::relay::io::hdf5_read(hdf5FilePath, root);

  // ----- Record 0 -----
  EXPECT_EQ(root["records/0/type"].as_string(), "foo1");
  EXPECT_EQ(root["records/0/id"].as_string(), "bar1");
  EXPECT_EQ(root["records/0/data/scal1/value"].as_string(), "goodbye!");
  EXPECT_EQ(root["records/0/data/scal2/value"].as_string(), "goodbye!");
  EXPECT_EQ(root["records/0/user_defined/hello"].as_string(), "goodbye");

  const conduit::Node &cs0 = root["records/0/curve_sets/1"];
  std::vector<double> dep0_vec = node_to_double_vector(cs0["dependent/0/value"]);
  std::vector<double> expected_dep0 = {1.0, 2.0, 1.0, 1.0};
  EXPECT_EQ(dep0_vec, expected_dep0);
  std::vector<double> dep1_vec = node_to_double_vector(cs0["dependent/1/value"]);
  std::vector<double> expected_dep1 = {10.0, 20.0, 8.0, 9.0};
  EXPECT_EQ(dep1_vec, expected_dep1);
  std::vector<double> ind0_vec = node_to_double_vector(cs0["independent/0/value"]);
  std::vector<double> expected_ind0 = {4.0, 5.0, 4.0, 5.0};
  EXPECT_EQ(ind0_vec, expected_ind0);
  std::vector<double> ind1_vec = node_to_double_vector(cs0["independent/1/value"]);
  std::vector<double> expected_ind1 = {7.0, 8.0, 7.0, 9.0};
  const conduit::Node &csnew = root["records/0/curve_sets/2"];
  std::vector<double> new_dep0 = node_to_double_vector(csnew["dependent/0/value"]);
  std::vector<double> expected_new_dep0 = {1.0, 1.0};
  EXPECT_EQ(new_dep0, expected_new_dep0);
  std::vector<double> new_dep1 = node_to_double_vector(csnew["dependent/1/value"]);
  std::vector<double> expected_new_dep1 = {10.0, 10.0};
  EXPECT_EQ(new_dep1, expected_new_dep1);
  std::vector<double> new_ind0 = node_to_double_vector(csnew["independent/0/value"]);
  std::vector<double> expected_new_ind0 = {4.0, 4.0};
  EXPECT_EQ(new_ind0, expected_new_ind0);
  std::vector<double> new_ind1 = node_to_double_vector(csnew["independent/1/value"]);
  std::vector<double> expected_new_ind1 = {7.0, 7.0};
  EXPECT_EQ(new_ind1, expected_new_ind1);

  // ----- Record 1 -----
  EXPECT_EQ(root["records/1/type"].as_string(), "foo2");
  EXPECT_EQ(root["records/1/id"].as_string(), "bar2");
  const conduit::Node &cs1 = root["records/1/curve_sets/1"];

  std::vector<double> rec1_dep0 = node_to_double_vector(cs1["dependent/0/value"]);
  std::vector<double> expected_rec1_dep0 = {1.0, 2.0};
  EXPECT_EQ(rec1_dep0, expected_rec1_dep0);
  std::vector<double> rec1_dep1 = node_to_double_vector(cs1["dependent/1/value"]);
  std::vector<double> expected_rec1_dep1 = {10.0, 20.0};
  EXPECT_EQ(rec1_dep1, expected_rec1_dep1);
  std::vector<double> rec1_ind0 = node_to_double_vector(cs1["independent/0/value"]);
  std::vector<double> expected_rec1_ind0 = {4.0, 5.0};
  EXPECT_EQ(rec1_ind0, expected_rec1_ind0);
  std::vector<double> rec1_ind1 = node_to_double_vector(cs1["independent/1/value"]);
  std::vector<double> expected_rec1_ind1 = {7.0, 8.0};
  EXPECT_EQ(rec1_ind1, expected_rec1_ind1);

  // ----- Record 2 -----
  EXPECT_EQ(root["records/2/type"].as_string(), "foo2");
  EXPECT_EQ(root["records/2/id"].as_string(), "bar3");

  const conduit::Node &cs2 = root["records/2/curve_sets/1"];
  std::vector<double> rec2_dep0 = node_to_double_vector(cs2["dependent/0/value"]);
  std::vector<double> expected_rec2_dep0 = {1.0, 1.0};
  EXPECT_EQ(rec2_dep0, expected_rec2_dep0);
  std::vector<double> rec2_dep1 = node_to_double_vector(cs2["dependent/1/value"]);
  std::vector<double> expected_rec2_dep1 = {10.0, 10.0};
  EXPECT_EQ(rec2_dep1, expected_rec2_dep1);
  std::vector<double> rec2_ind0 = node_to_double_vector(cs2["independent/0/value"]);
  std::vector<double> expected_rec2_ind0 = {4.0, 4.0};
  EXPECT_EQ(rec2_ind0, expected_rec2_ind0);
  std::vector<double> rec2_ind1 = node_to_double_vector(cs2["independent/1/value"]);
  std::vector<double> expected_rec2_ind1 = {7.0, 7.0};
  EXPECT_EQ(rec2_ind1, expected_rec2_ind1);
}

TEST(Document, test_append_to_hdf5_failures)
{
  std::streambuf *origBuffer = std::cerr.rdbuf();
  std::ostringstream capturedCerr;
  std::cerr.rdbuf(capturedCerr.rdbuf());

  std::string hdf5FilePath = "test_failure.hdf5";
  conduit::Node initialData;
  initialData["records/0/data/scal1/value"].set("hello!");
  initialData["records/0/type"].set("foo1");
  initialData["records/0/id"].set("bar1");
  initialData["records/0/user_defined/hello"].set("and");
  initialData["records/0/curve_sets/1/dependent/0/value"].set(std::vector<double> {1.0, 2.0});
  initialData["records/0/curve_sets/1/dependent/1/value"].set(std::vector<double> {10.0, 20.0});
  initialData["records/0/curve_sets/1/independent/0/value"].set(std::vector<double> {4.0, 5.0});
  initialData["records/0/curve_sets/1/independent/1/value"].set(std::vector<double> {7.0, 8.0});
  initialData["records/1/type"].set("foo2");
  initialData["records/1/id"].set("bar2");
  initialData["records/1/curve_sets/1/dependent/0/value"].set(std::vector<double> {1.0, 2.0});
  initialData["records/1/curve_sets/1/dependent/1/value"].set(std::vector<double> {10.0, 20.0});
  initialData["records/1/curve_sets/1/independent/0/value"].set(std::vector<double> {4.0, 5.0});
  initialData["records/1/curve_sets/1/independent/1/value"].set(std::vector<double> {7.0, 8.0});
  conduit::relay::io::hdf5_write(initialData, hdf5FilePath);

  // ---- Failure 1 ----
  SCOPED_TRACE("Attempting Case 1");
  std::string f1_data = ([](std::string s) {
    conduit::Node tmp;
    tmp.parse(s, "json");
    conduit::Node &valNode = tmp["records"].child(0)["curve_sets"]["1"]["dependent"]["0"]["value"];
    valNode.reset();
    std::vector<double> new_val = {1, 1, 1};
    valNode.set(new_val);
    return tmp.to_json();
  })(new_data);

  axom::sina::Document new_doc = Document(f1_data, createRecordLoaderWithAllKnownTypes());

  conduit::Node msgNode = append_to_hdf5(hdf5FilePath, new_doc, 1);
  ASSERT_FALSE(msgNode.number_of_children() == 0);
  conduit::Node root;
  conduit::relay::io::hdf5_read(hdf5FilePath, root);
  ValidateRecordsUnchanged(root, "Case 1: Checking if File Was Changed");

  // ---- Failure 2 ----
  SCOPED_TRACE("Attempting Case 2");
  std::string f2_data = ([](std::string s) {
    conduit::Node tmp;
    tmp.parse(s, "json");
    conduit::Node &curve_set1 = tmp["records"].child(0)["curve_sets"]["1"];
    if(curve_set1.has_child("dependent")) curve_set1["dependent"].reset();
    conduit::Node &newDep = curve_set1["dependent"]["2"];
    newDep.reset();
    std::vector<double> newVal = {1};
    newDep["value"].set(newVal);
    if(curve_set1.has_child("independent")) curve_set1["independent"].reset();
    return tmp.to_json();
  })(new_data);

  /*new_doc = Document(f2_data, createRecordLoaderWithAllKnownTypes());

  result = append_to_hdf5(hdf5FilePath, new_doc, 1);
  ASSERT_FALSE(result);
  conduit::relay::io::hdf5_read(hdf5FilePath, root);
  ValidateRecordsUnchanged(root, "Case 2: Checking if File Was Changed");

  // ---- Failure 3 ----
  SCOPED_TRACE("Attempting Case 3");
  std::string f3_data = ([](std::string s) {
    conduit::Node tmp;
    tmp.parse(s, "json");
    conduit::Node &valNode =
      tmp["records"].child(0)["curve_sets"]["1"]["independent"]["1"]["value"];
    valNode.reset();
    std::vector<double> new_val = {1, 1, 1};
    valNode.set(new_val);
    return tmp.to_json();
  })(new_data);

  new_doc = Document(f3_data, createRecordLoaderWithAllKnownTypes());

  result = append_to_hdf5(hdf5FilePath, new_doc, 1);
  ASSERT_FALSE(result);
  conduit::relay::io::hdf5_read(hdf5FilePath, root);
  ValidateRecordsUnchanged(root, "Case 3: Checking if File Was Changed");

  // ---- Failure 4 ----
  SCOPED_TRACE("Attempting Case 4");
  new_doc = Document(new_data, createRecordLoaderWithAllKnownTypes());

  result = append_to_hdf5(hdf5FilePath, new_doc, 3);
  ASSERT_FALSE(result);
  conduit::relay::io::hdf5_read(hdf5FilePath, root);
  ValidateRecordsUnchanged(root, "Case 4: Checking if File Was Changed");

  // ---- Failure 5 ----
  SCOPED_TRACE("Attempting Case 5");
  new_doc = Document(new_data, createRecordLoaderWithAllKnownTypes());

  result = append_to_hdf5(hdf5FilePath, new_doc, 1);
  ASSERT_FALSE(result);
  conduit::relay::io::hdf5_read(hdf5FilePath, root);
  ValidateRecordsUnchanged(root, "Case 5: Checking if File Was Changed");

  // ---- Assert Error Messages ----
  std::cerr.rdbuf(origBuffer);
  std::string expected =
    "Error validating dependent: Record bar1, Curve Set 1, Curve 1 size mismatch (expected 5, got "
    "4).\n"
    "Error validating dependent: Record bar1, Curve Set 1, Curve 2 size mismatch (expected 2, got "
    "1).\n"
    "Error validating independent: Record bar1, Curve Set 1, Curve 1 size mismatch (expected 4, "
    "got 5).\n"
    "Found a duplicate data entry, protocol 3 dictates append cancellation.\n"
    "Found duplicate UDC, protocol 3 dictates append cancellation.\n";
  ASSERT_EQ(capturedCerr.str(), expected);*/
}
#endif

}  // namespace
}  // namespace testing
}  // namespace sina
}  // namespace axom
