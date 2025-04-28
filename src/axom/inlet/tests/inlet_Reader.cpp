// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/inlet/tests/inlet_test_utils.hpp"

#include "gtest/gtest.h"

#include <string>
#include <vector>
#include <memory>

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
