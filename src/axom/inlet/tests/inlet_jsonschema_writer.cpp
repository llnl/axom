// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic.hpp"

#include "axom/config.hpp"
#include "axom/inlet/Inlet.hpp"
#include "axom/inlet/JSONSchemaWriter.hpp"
#include "axom/inlet/tests/inlet_test_utils.hpp"

#include "gtest/gtest.h"

#include <cstdlib>
#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <streambuf>

using axom::inlet::Inlet;
using axom::inlet::JSONSchemaWriter;

bool hasSchemaUtility()
{
  const static bool hasUtility = std::system(JSONSCHEMA_EXECUTABLE " --version") == 0;
  return hasUtility;
}

bool validateString(Inlet& inlet, const std::string& luaString)
{
  const static auto instanceFile = "test_input.json";
  const static auto schemaFile = "input_schema.json";

  std::ofstream jsonInstance(instanceFile);
  jsonInstance << axom::inlet::detail::fromLuaTo<axom::inlet::JSONReader>(luaString);
  jsonInstance.close();

  inlet.write(JSONSchemaWriter(schemaFile));

  const std::string cmd =
    axom::fmt::format(JSONSCHEMA_EXECUTABLE " -i {0} {1}", instanceFile, schemaFile);
  return std::system(cmd.c_str()) == 0;
}

template <typename InletReader>
Inlet createBasicInlet(const std::string& luaString)
{
  std::unique_ptr<InletReader> reader(new InletReader());
  reader->parseString(axom::inlet::detail::fromLuaTo<InletReader>(luaString));
  const bool enableDocs = true;
  return Inlet(std::move(reader), enableDocs);
}

template <typename InletReader>
class inlet_jsonschema_writer : public ::testing::Test
{
protected:
  void SetUp() override
  {
    if(!hasSchemaUtility())
    {
      GTEST_SKIP() << "jsonschema tool was not found";
    }
  }
};

TYPED_TEST_SUITE(inlet_jsonschema_writer, axom::inlet::detail::ReaderTypes);

TYPED_TEST(inlet_jsonschema_writer, top_level_bools)
{
  std::string testString = "foo = true; bar = false";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  // Check for existing fields
  inlet.addBool("foo", "foo's description");
  inlet.addBool("bar", "bar's description");

  // Check one that doesn't exist and doesn't have a default value
  inlet.addBool("non/existant", "nothing");

  const bool jsonVerified = validateString(inlet, testString);

  EXPECT_TRUE(jsonVerified);
  EXPECT_EQ(inlet.verify(), jsonVerified);
}

TYPED_TEST(inlet_jsonschema_writer, top_level_bools_reqd)
{
  std::string testString = "foo = true";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  // Check for existing fields
  inlet.addBool("foo", "foo's description");
  inlet.addBool("bar", "bar's description").required();

  // Check one that doesn't exist and doesn't have a default value
  inlet.addBool("non/existant", "nothing");

  const bool jsonVerified = validateString(inlet, testString);

  EXPECT_FALSE(jsonVerified);
  EXPECT_EQ(inlet.verify(), jsonVerified);
}

TYPED_TEST(inlet_jsonschema_writer, top_level_ints)
{
  std::string testString = "foo = 12; bar = 16";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addInt("foo", "foo's description");
  inlet.addInt("bar", "bar's description");

  const bool jsonVerified = validateString(inlet, testString);

  EXPECT_TRUE(jsonVerified);
  EXPECT_EQ(inlet.verify(), jsonVerified);
}

TYPED_TEST(inlet_jsonschema_writer, top_level_ints_wrong_type)
{
  std::string testString = "foo = 'first'; bar = 'second'";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  // We specify the default as a "hint" to the JSONSchemaWriter
  // since expected type information is not stored explicitly in the datastore
  inlet.addInt("foo", "foo's description").defaultValue(11);
  inlet.addInt("bar", "bar's description").defaultValue(10);

  const bool jsonVerified = validateString(inlet, testString);

  EXPECT_FALSE(jsonVerified);
  EXPECT_EQ(inlet.verify(), jsonVerified);
}

TYPED_TEST(inlet_jsonschema_writer, top_level_ints_range_pass)
{
  std::string testString = "foo = 5; bar = 12";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addInt("foo", "foo's description").range(2, 10);
  inlet.addInt("bar", "bar's description").range(10, 15);

  const bool jsonVerified = validateString(inlet, testString);

  EXPECT_TRUE(jsonVerified);
  EXPECT_EQ(inlet.verify(), jsonVerified);
}

TYPED_TEST(inlet_jsonschema_writer, top_level_ints_range_fail)
{
  std::string testString = "foo = 12; bar = 5";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addInt("foo", "foo's description").range(2, 10);
  inlet.addInt("bar", "bar's description").range(10, 15);

  const bool jsonVerified = validateString(inlet, testString);

  EXPECT_FALSE(jsonVerified);
  EXPECT_EQ(inlet.verify(), jsonVerified);
}

TYPED_TEST(inlet_jsonschema_writer, top_level_ints_valid_set)
{
  std::string testString = "foo = 5; bar = 8";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addInt("foo", "foo's description").validValues({5, 7, 9});
  inlet.addInt("bar", "bar's description").validValues({4, 6, 8});

  const bool jsonVerified = validateString(inlet, testString);

  EXPECT_TRUE(jsonVerified);
  EXPECT_EQ(inlet.verify(), jsonVerified);
}

TYPED_TEST(inlet_jsonschema_writer, top_level_ints_valid_set_fail)
{
  std::string testString = "foo = 8; bar = 5";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addInt("foo", "foo's description").validValues({5, 7, 9});
  inlet.addInt("bar", "bar's description").validValues({4, 6, 8});

  const bool jsonVerified = validateString(inlet, testString);

  EXPECT_FALSE(jsonVerified);
  EXPECT_EQ(inlet.verify(), jsonVerified);
}

// Valid string values use different logic than other "simple" primitives
// so we test this separately
TYPED_TEST(inlet_jsonschema_writer, top_level_strings_valid_set)
{
  std::string testString = "foo = 'first'; bar = 'second'";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addString("foo", "foo's description").validValues({"first", "second"});
  inlet.addString("bar", "bar's description").validValues({"first", "second"});

  const bool jsonVerified = validateString(inlet, testString);

  EXPECT_TRUE(jsonVerified);
  EXPECT_EQ(inlet.verify(), jsonVerified);
}

TYPED_TEST(inlet_jsonschema_writer, top_level_strings_valid_set_fail)
{
  std::string testString = "foo = 'first'; bar = 'third'";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addString("foo", "foo's description").validValues({"first", "second"});
  inlet.addString("bar", "bar's description").validValues({"first", "second"});

  const bool jsonVerified = validateString(inlet, testString);

  EXPECT_FALSE(jsonVerified);
  EXPECT_EQ(inlet.verify(), jsonVerified);
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
