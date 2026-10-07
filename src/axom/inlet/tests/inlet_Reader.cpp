// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/inlet/tests/inlet_test_utils.hpp"

#include "gtest/gtest.h"

#include <algorithm>
#include <string>
#include <vector>
#include <memory>
#include <variant>

template <typename InletReader>
class inlet_Reader : public testing::Test
{ };

TYPED_TEST_SUITE(inlet_Reader, axom::inlet::detail::ReaderTypes);

using axom::inlet::ReaderResult;
using axom::inlet::detail::fromLuaTo;

TYPED_TEST(inlet_Reader, getTopLevelBools)
{
  TypeParam reader;
  reader.parseString(fromLuaTo<TypeParam>("foo = true; bar = false"));

  ReaderResult retValue;
  bool value;

  value = false;
  retValue = reader.getBool("foo", value);
  EXPECT_EQ(retValue, ReaderResult::Success);
  EXPECT_EQ(value, true);

  value = true;
  retValue = reader.getBool("bar", value);
  EXPECT_EQ(retValue, ReaderResult::Success);
  EXPECT_EQ(value, false);
}

TYPED_TEST(inlet_Reader, getTopLevelBoolsWrongType)
{
  TypeParam reader;
  reader.parseString(fromLuaTo<TypeParam>("foo = true; bar = false"));

  ReaderResult retValue;
  double value;

  retValue = reader.getDouble("foo", value);
  EXPECT_EQ(retValue, ReaderResult::WrongType);

  value = true;
  retValue = reader.getDouble("bar", value);
  EXPECT_EQ(retValue, ReaderResult::WrongType);
}

TYPED_TEST(inlet_Reader, getInsideBools)
{
  TypeParam reader;
  reader.parseString(fromLuaTo<TypeParam>("foo = { bar = false; baz = true }"));

  ReaderResult retValue;
  bool value;

  value = true;
  retValue = reader.getBool("foo/bar", value);
  EXPECT_EQ(retValue, ReaderResult::Success);
  EXPECT_EQ(value, false);

  value = false;
  retValue = reader.getBool("foo/baz", value);
  EXPECT_EQ(retValue, ReaderResult::Success);
  EXPECT_EQ(value, true);
}

TYPED_TEST(inlet_Reader, getTopLevelStrings)
{
  TypeParam reader;
  reader.parseString(
    fromLuaTo<TypeParam>("foo = \"this is a test string\"; bar = \"TesT StrInG\""));

  ReaderResult retValue;
  std::string value;

  value = "";
  retValue = reader.getString("foo", value);
  EXPECT_EQ(retValue, ReaderResult::Success);
  EXPECT_EQ(value, "this is a test string");

  value = "";
  retValue = reader.getString("bar", value);
  EXPECT_EQ(retValue, ReaderResult::Success);
  EXPECT_EQ(value, "TesT StrInG");
}

TYPED_TEST(inlet_Reader, getInsideStrings)
{
  TypeParam reader;
  reader.parseString(
    fromLuaTo<TypeParam>("foo = { bar = \"this is a test string\"; baz = \"TesT StrInG\" }"));

  ReaderResult retValue;
  std::string value;

  value = "";
  retValue = reader.getString("foo/bar", value);
  EXPECT_EQ(retValue, ReaderResult::Success);
  EXPECT_EQ(value, "this is a test string");

  value = "";
  retValue = reader.getString("foo/baz", value);
  EXPECT_EQ(retValue, ReaderResult::Success);
  EXPECT_EQ(value, "TesT StrInG");
}

TYPED_TEST(inlet_Reader, mixLevelContainers)
{
  TypeParam reader;
  reader.parseString(
    fromLuaTo<TypeParam>("t = { innerT = { foo = 1 }, anotherInnerT = {baz = 3}}"));

  ReaderResult retValue;
  int value;

  value = 0;
  retValue = reader.getInt("t/innerT/foo", value);
  EXPECT_EQ(retValue, ReaderResult::Success);
  EXPECT_EQ(value, 1);

  value = 0;
  retValue = reader.getInt("t/doesntexist", value);
  EXPECT_EQ(retValue, ReaderResult::NotFound);
  EXPECT_EQ(value, 0);

  value = 0;
  retValue = reader.getInt("t/anotherInnerT/baz", value);
  EXPECT_EQ(retValue, ReaderResult::Success);
  EXPECT_EQ(value, 3);
}

TYPED_TEST(inlet_Reader, getMap)
{
  // Keep this contiguous in order to test all supported input languages
  std::string testString =
    "luaArray = { [0] = 4, [1] = 5, [2] = 6 , [3] = true, [4] = false, [5] = "
    "2.4, [6] = 'hello', [7] = 'bye' }";
  TypeParam reader;
  reader.parseString(fromLuaTo<TypeParam>(testString));

  std::unordered_map<int, int> ints;
  ReaderResult retValue = reader.getIntMap("luaArray", ints);
  EXPECT_EQ(retValue, ReaderResult::NotHomogeneous);
  std::unordered_map<int, int> expectedInts {{0, 4}, {1, 5}, {2, 6}, {5, 2}};
  EXPECT_EQ(expectedInts, ints);

  std::unordered_map<int, double> doubles;
  retValue = reader.getDoubleMap("luaArray", doubles);
  EXPECT_EQ(retValue, ReaderResult::NotHomogeneous);
  std::unordered_map<int, double> expectedDoubles {{0, 4}, {1, 5}, {2, 6}, {5, 2.4}};
  EXPECT_EQ(expectedDoubles, doubles);

  std::unordered_map<int, bool> bools;
  retValue = reader.getBoolMap("luaArray", bools);
  EXPECT_EQ(retValue, ReaderResult::NotHomogeneous);
  std::unordered_map<int, bool> expectedBools {{3, true}, {4, false}};
  EXPECT_EQ(expectedBools, bools);

  // Conduit's YAML parser doesn't distinguish boolean literals from strings
  // so the YAML version will extract the "true" and "false" here
  std::unordered_map<int, std::string> strs;
  retValue = reader.getStringMap("luaArray", strs);
  EXPECT_EQ(retValue, ReaderResult::NotHomogeneous);
  // std::unordered_map<int, std::string> expectedStrs {{6, "hello"},
  //                                                    {7, "bye"}};
  // EXPECT_EQ(expectedStrs, strs);
}

TYPED_TEST(inlet_Reader, getVariantMap)
{
  std::string testString = "luaArray = { [0] = 42, [1] = 'hello', [2] = true, [3] = 3.14 }";
  TypeParam reader;
  reader.parseString(fromLuaTo<TypeParam>(testString));

  std::unordered_map<int, axom::inlet::VariantValue> values;
  ReaderResult retValue = reader.getVariantMap("luaArray", values);
  EXPECT_EQ(retValue, ReaderResult::Success);

  std::unordered_map<int, axom::inlet::VariantValue> expected {
    {0, axom::inlet::VariantValue {42}},
    {1, axom::inlet::VariantValue {std::string {"hello"}}},
    {2, axom::inlet::VariantValue {true}},
    {3, axom::inlet::VariantValue {3.14}}};
  EXPECT_EQ(expected, values);
}

TEST(inlet_Reader_JSON, getVariantMapBoolArray)
{
  axom::inlet::JSONReader reader;
  bool result = reader.parseString("{\"bools\": [true, false, true]}");
  EXPECT_TRUE(result);

  std::unordered_map<int, axom::inlet::VariantValue> values;
  ReaderResult retValue = reader.getVariantMap("bools", values);
  EXPECT_EQ(retValue, ReaderResult::Success);

  std::unordered_map<int, axom::inlet::VariantValue> expected {
    {0, axom::inlet::VariantValue {true}},
    {1, axom::inlet::VariantValue {false}},
    {2, axom::inlet::VariantValue {true}}};
  EXPECT_EQ(expected, values);
}

TEST(inlet_Reader_JSON, getVariantMapStringArray)
{
  axom::inlet::JSONReader reader;
  bool result = reader.parseString("{\"strings\": [\"red\", \"green\", \"blue\"]}");
  EXPECT_TRUE(result);

  std::unordered_map<int, axom::inlet::VariantValue> values;
  ReaderResult retValue = reader.getVariantMap("strings", values);
  EXPECT_EQ(retValue, ReaderResult::Success);

  std::unordered_map<int, axom::inlet::VariantValue> expected {
    {0, axom::inlet::VariantValue {std::string {"red"}}},
    {1, axom::inlet::VariantValue {std::string {"green"}}},
    {2, axom::inlet::VariantValue {std::string {"blue"}}}};
  EXPECT_EQ(expected, values);
}

TYPED_TEST(inlet_Reader, emptyCollections)
{
  TypeParam reader;
  reader.parseString(fromLuaTo<TypeParam>("arr = { }"));

  ReaderResult retValue;
  std::vector<axom::inlet::VariantKey> indices;
  std::vector<axom::inlet::VariantKey> expected_indices;

  retValue = reader.getIndices("arr", indices);
  EXPECT_EQ(retValue, ReaderResult::Success);
  EXPECT_EQ(indices, expected_indices);

  retValue = reader.getIndices("nonexistent_arr", indices);
  EXPECT_EQ(retValue, ReaderResult::NotFound);
  EXPECT_EQ(indices, expected_indices);
}

TYPED_TEST(inlet_Reader, simple_name_retrieval)
{
  TypeParam reader;
  reader.parseString(
    fromLuaTo<TypeParam>("t = { innerT = { foo = 1 }, anotherInnerT = {baz = 3}}"));

  auto found_names = reader.getAllNames();
  std::vector<std::string> expected_names {"t",
                                           "t/innerT",
                                           "t/innerT/foo",
                                           "t/anotherInnerT",
                                           "t/anotherInnerT/baz"};
  std::sort(found_names.begin(), found_names.end());
  std::sort(expected_names.begin(), expected_names.end());
  EXPECT_EQ(found_names, expected_names);
}

TYPED_TEST(inlet_Reader, simple_name_retrieval_arrays)
{
  TypeParam reader;
  reader.parseString(
    fromLuaTo<TypeParam>("t = { [0] = { foo = 1, bar = 2}, [1] = { foo = 3, bar = 4} }"));

  auto found_names = reader.getAllNames();
  std::vector<std::string> expected_names {
    "t",
    "t/0",
    "t/0/foo",
    "t/0/bar",
    "t/1",
    "t/1/foo",
    "t/1/bar",
  };
  std::sort(found_names.begin(), found_names.end());
  std::sort(expected_names.begin(), expected_names.end());
  EXPECT_EQ(found_names, expected_names);
}

TYPED_TEST(inlet_Reader, intReadsTruncateTowardZero)
{
  // Every reader narrows a non-integral number to an int the same way
  TypeParam reader;
  reader.parseString(fromLuaTo<TypeParam>("up = 2.7; down = -2.7; half = 2.5; whole = 7"));

  int value = 0;
  EXPECT_EQ(ReaderResult::Success, reader.getInt("up", value));
  EXPECT_EQ(2, value);
  EXPECT_EQ(ReaderResult::Success, reader.getInt("down", value));
  EXPECT_EQ(-2, value);
  EXPECT_EQ(ReaderResult::Success, reader.getInt("half", value));
  EXPECT_EQ(2, value);
  EXPECT_EQ(ReaderResult::Success, reader.getInt("whole", value));
  EXPECT_EQ(7, value);
}

TEST(inlet_Reader_YAML, getInsideBools)
{
  axom::inlet::YAMLReader reader;
  bool result = reader.parseString(
    "foo:\n"
    "  bar: false\n"
    "  baz: true");
  EXPECT_TRUE(result);

  ReaderResult retValue;
  bool value;

  value = true;
  retValue = reader.getBool("foo/bar", value);
  EXPECT_EQ(retValue, ReaderResult::Success);
  EXPECT_EQ(value, false);

  value = false;
  retValue = reader.getBool("foo/baz", value);
  EXPECT_EQ(retValue, ReaderResult::Success);
  EXPECT_EQ(value, true);
}

TEST(inlet_Reader_YAML, mixLevelContainers)
{
  axom::inlet::YAMLReader reader;
  bool result = reader.parseString(
    "t:\n"
    "  innerT:\n"
    "    foo: 1\n"
    "  anotherInnerT:\n"
    "    baz: 3");
  EXPECT_TRUE(result);

  ReaderResult retValue;
  int value;

  value = 0;
  retValue = reader.getInt("t/innerT/foo", value);
  EXPECT_EQ(retValue, ReaderResult::Success);
  EXPECT_EQ(value, 1);

  value = 0;
  retValue = reader.getInt("t/doesntexist", value);
  EXPECT_EQ(retValue, ReaderResult::NotFound);
  EXPECT_EQ(value, 0);

  value = 0;
  retValue = reader.getInt("t/anotherInnerT/baz", value);
  EXPECT_EQ(retValue, ReaderResult::Success);
  EXPECT_EQ(value, 3);
}

TEST(inlet_Reader_YAML, mixLevelContainers_invalid)
{
  axom::inlet::YAMLReader reader;
  bool result = reader.parseString(
    "t:\n"
    "  innerT: foo: 1\n"
    "  anotherInnerT:\n"
    "    baz: 3");

  EXPECT_FALSE(result);
}

TEST(inlet_Reader_JSON, getInsideBools)
{
  axom::inlet::JSONReader reader;
  bool result = reader.parseString(
    "{\n"
    "  foo: {\n"
    "    bar: false,\n"
    "    baz: true\n"
    "  }\n"
    "}");
  EXPECT_TRUE(result);

  ReaderResult retValue;
  bool value;

  value = true;
  retValue = reader.getBool("foo/bar", value);
  EXPECT_EQ(retValue, ReaderResult::Success);
  EXPECT_EQ(value, false);

  value = false;
  retValue = reader.getBool("foo/baz", value);
  EXPECT_EQ(retValue, ReaderResult::Success);
  EXPECT_EQ(value, true);
}

TEST(inlet_Reader_JSON, mixLevelContainers)
{
  axom::inlet::JSONReader reader;
  bool result = reader.parseString(
    "{\n"
    "  t: {\n"
    "    innerT: {\n"
    "      foo: 1\n"
    "    },\n"
    "    anotherInnerT: {\n"
    "      baz: 3\n"
    "    }\n"
    "  }\n"
    "}");
  EXPECT_TRUE(result);

  ReaderResult retValue;
  int value;

  value = 0;
  retValue = reader.getInt("t/innerT/foo", value);
  EXPECT_EQ(retValue, ReaderResult::Success);
  EXPECT_EQ(value, 1);

  value = 0;
  retValue = reader.getInt("t/doesntexist", value);
  EXPECT_EQ(retValue, ReaderResult::NotFound);
  EXPECT_EQ(value, 0);

  value = 0;
  retValue = reader.getInt("t/anotherInnerT/baz", value);
  EXPECT_EQ(retValue, ReaderResult::Success);
  EXPECT_EQ(value, 3);
}

TEST(inlet_Reader_JSON, mixLevelContainers_invalid)
{
  axom::inlet::JSONReader reader;
  bool result = reader.parseString(
    "{\n"
    "  t: {\n"
    "    innerT: {\n"
    "      foo: 1\n"
    "    }\n"
    "    anotherInnerT: {\n"
    "      baz: 3\n"
    "    }\n"
    "  }\n"
    "}");

  EXPECT_FALSE(result);
}

#ifdef AXOM_USE_SOL
// Checks that LuaReader parses array information as expected
// Discontiguous arrays are lua-specific
TEST(inlet_Reader_lua, getDiscontiguousMap)
{
  std::string testString =
    "luaArray = { [1] = 4, [2] = 5, [3] = 6 , [4] = true, [8] = false, [12] = "
    "2.4, [33] = 'hello', [200] = 'bye' }";
  axom::inlet::LuaReader reader;
  reader.parseString(testString);

  std::unordered_map<int, int> ints;
  ReaderResult retValue = reader.getIntMap("luaArray", ints);
  EXPECT_EQ(retValue, ReaderResult::NotHomogeneous);
  std::unordered_map<int, int> expectedInts {{1, 4}, {2, 5}, {3, 6}, {12, 2}};
  EXPECT_EQ(expectedInts, ints);

  std::unordered_map<int, double> doubles;
  retValue = reader.getDoubleMap("luaArray", doubles);
  EXPECT_EQ(retValue, ReaderResult::NotHomogeneous);
  std::unordered_map<int, double> expectedDoubles {{1, 4}, {2, 5}, {3, 6}, {12, 2.4}};
  EXPECT_EQ(expectedDoubles, doubles);

  std::unordered_map<int, bool> bools;
  retValue = reader.getBoolMap("luaArray", bools);
  EXPECT_EQ(retValue, ReaderResult::NotHomogeneous);
  std::unordered_map<int, bool> expectedBools {{4, true}, {8, false}};
  EXPECT_EQ(expectedBools, bools);

  std::unordered_map<int, std::string> strs;
  retValue = reader.getStringMap("luaArray", strs);
  EXPECT_EQ(retValue, ReaderResult::NotHomogeneous);
  std::unordered_map<int, std::string> expectedStrs {{33, "hello"}, {200, "bye"}};
  EXPECT_EQ(expectedStrs, strs);
}

TEST(inlet_Reader_lua, objectLookupReportsConsistentReaderResults)
{
  axom::inlet::LuaReader reader;
  reader.parseString(R"(
    callback = function() return {1, 2} end
    nested = {[7] = {values = {[2] = 42, [5] = "five"}}}
  )");

  double scalar = 0.0;
  EXPECT_EQ(ReaderResult::WrongType, reader.getDouble("callback", scalar));
  EXPECT_EQ(ReaderResult::NotFound, reader.getDouble("callback/value", scalar));

  std::unordered_map<int, double> typedValues {{99, 99.0}};
  EXPECT_EQ(ReaderResult::WrongType, reader.getDoubleMap("callback", typedValues));
  EXPECT_TRUE(typedValues.empty());

  std::unordered_map<int, axom::inlet::VariantValue> values {{99, axom::inlet::VariantValue {99}}};
  EXPECT_EQ(ReaderResult::WrongType, reader.getVariantMap("callback", values));
  EXPECT_TRUE(values.empty());
  EXPECT_EQ(ReaderResult::NotFound, reader.getVariantMap("missing", values));
  EXPECT_TRUE(values.empty());

  EXPECT_EQ(ReaderResult::Success, reader.getVariantMap("nested/7/values", values));
  const std::unordered_map<int, axom::inlet::VariantValue> expectedValues {
    {2, axom::inlet::VariantValue {42}},
    {5, axom::inlet::VariantValue {std::string {"five"}}}};
  EXPECT_EQ(expectedValues, values);

  std::vector<int> indices {99};
  EXPECT_EQ(ReaderResult::WrongType, reader.getIndices("callback", indices));
  EXPECT_TRUE(indices.empty());
  EXPECT_EQ(ReaderResult::Success, reader.getIndices("nested/7/values", indices));
  std::sort(indices.begin(), indices.end());
  EXPECT_EQ((std::vector<int> {2, 5}), indices);
}

TEST(inlet_Reader_lua, nonTableValuesAreWrongTypeForCollections)
{
  // Regression test against read that terminated the process with an uncaught sol::error
  axom::inlet::LuaReader reader;
  reader.parseString("scalar = 3.0; label = 'text'; callback = function() return 1 end");

  for(const std::string name : {"scalar", "label", "callback"})
  {
    std::unordered_map<int, double> doubles {{99, 99.0}};
    EXPECT_EQ(ReaderResult::WrongType, reader.getDoubleMap(name, doubles)) << name;
    EXPECT_TRUE(doubles.empty()) << name;

    std::unordered_map<axom::inlet::VariantKey, std::string> strings {{99, "stale"}};
    EXPECT_EQ(ReaderResult::WrongType, reader.getStringMap(name, strings)) << name;
    EXPECT_TRUE(strings.empty()) << name;

    std::unordered_map<axom::inlet::VariantKey, axom::inlet::VariantValue> variants {
      {99, axom::inlet::VariantValue {99}}};
    EXPECT_EQ(ReaderResult::WrongType, reader.getVariantMap(name, variants)) << name;
    EXPECT_TRUE(variants.empty()) << name;

    std::vector<axom::inlet::VariantKey> indices {99};
    EXPECT_EQ(ReaderResult::WrongType, reader.getIndices(name, indices)) << name;
    EXPECT_TRUE(indices.empty()) << name;
  }
}

TEST(inlet_Reader_lua, pathsThroughNonTablesAreNotFound)
{
  // Check that indexing through a non-table doesn't terminate the process with an uncaught sol::error
  axom::inlet::LuaReader reader;
  reader.parseString("scalar = 3.0; label = 'text'; callback = function() return 1 end");

  for(const std::string parent : {"scalar", "label", "callback"})
  {
    const std::string name = parent + "/child";
    double value = -1.0;
    EXPECT_EQ(ReaderResult::NotFound, reader.getDouble(name, value)) << name;
    EXPECT_DOUBLE_EQ(-1.0, value) << name;

    std::unordered_map<int, double> doubles {{99, 99.0}};
    EXPECT_EQ(ReaderResult::NotFound, reader.getDoubleMap(name, doubles)) << name;
    EXPECT_TRUE(doubles.empty()) << name;

    std::vector<int> indices {99};
    EXPECT_EQ(ReaderResult::NotFound, reader.getIndices(name, indices)) << name;
    EXPECT_TRUE(indices.empty()) << name;

    EXPECT_FALSE(reader.getFunction(name, axom::inlet::FunctionTag::Double, {})) << name;
  }
}

TEST(inlet_Reader_lua, nonFunctionValuesAreNotFunctions)
{
  // Check that a top-level non-function doesn't terminate the process with an uncaught sol::error
  axom::inlet::LuaReader reader;
  reader.parseString("scalar = 3.0; label = 'text'; group = {scalar = 3.0}");

  for(const std::string name : {"scalar", "label", "group", "group/scalar", "missing"})
  {
    EXPECT_FALSE(reader.getFunction(name, axom::inlet::FunctionTag::Double, {})) << name;
  }
}

TEST(inlet_Reader_lua, numericPathComponentsPreferIntegerKeys)
{
  // Every component of a path, including the last, is looked up as an integer key
  // when it is numeric, and as a string key otherwise or if no integer key exists
  axom::inlet::LuaReader reader;
  reader.parseString(R"(
    values = {5.0, 6.0}
    both = {[1] = 1.0, ["1"] = 2.0}
    stringKeyed = {["1"] = 7.0}
    callbacks = {function() return 9.0 end}
    nested = {[3] = {[4] = 8.0}}
  )");

  double value = -1.0;
  EXPECT_EQ(ReaderResult::Success, reader.getDouble("values/2", value));
  EXPECT_DOUBLE_EQ(6.0, value);
  EXPECT_EQ(ReaderResult::Success, reader.getDouble("both/1", value));
  EXPECT_DOUBLE_EQ(1.0, value);
  EXPECT_EQ(ReaderResult::Success, reader.getDouble("stringKeyed/1", value));
  EXPECT_DOUBLE_EQ(7.0, value);
  EXPECT_EQ(ReaderResult::Success, reader.getDouble("nested/3/4", value));
  EXPECT_DOUBLE_EQ(8.0, value);
  EXPECT_EQ(ReaderResult::NotFound, reader.getDouble("values/3", value));

  auto callback = reader.getFunction("callbacks/1", axom::inlet::FunctionTag::Double, {});
  ASSERT_TRUE(callback);
  EXPECT_DOUBLE_EQ(9.0, callback.call<double>());
}

TEST(inlet_Reader_lua, scalarReadsCheckTheLuaType)
{
  // Lua coerces numeric strings to numbers, but Inlet reads do not
  axom::inlet::LuaReader reader;
  reader.parseString("numeric = '3'; flag = true; number = 3; text = 'text'");

  int intValue = -1;
  double doubleValue = -1.0;
  bool boolValue = false;
  std::string stringValue = "unset";

  EXPECT_EQ(ReaderResult::WrongType, reader.getInt("numeric", intValue));
  EXPECT_EQ(ReaderResult::WrongType, reader.getDouble("numeric", doubleValue));
  EXPECT_EQ(ReaderResult::WrongType, reader.getInt("flag", intValue));
  EXPECT_EQ(ReaderResult::WrongType, reader.getDouble("flag", doubleValue));
  EXPECT_EQ(ReaderResult::WrongType, reader.getBool("number", boolValue));
  EXPECT_EQ(ReaderResult::WrongType, reader.getString("number", stringValue));
  EXPECT_EQ(ReaderResult::WrongType, reader.getBool("text", boolValue));
  EXPECT_EQ(-1, intValue);
  EXPECT_DOUBLE_EQ(-1.0, doubleValue);
  EXPECT_EQ("unset", stringValue);

  EXPECT_EQ(ReaderResult::Success, reader.getString("numeric", stringValue));
  EXPECT_EQ("3", stringValue);
  EXPECT_EQ(ReaderResult::Success, reader.getBool("flag", boolValue));
  EXPECT_TRUE(boolValue);
}

TEST(inlet_Reader_lua, getIndicesClearsOutputWhenNotFound)
{
  axom::inlet::LuaReader reader;
  reader.parseString("values = {1, 2}");

  std::vector<int> indices {99};
  EXPECT_EQ(ReaderResult::NotFound, reader.getIndices("missing", indices));
  EXPECT_TRUE(indices.empty());
}
#endif

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
