// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/Path.hpp"
#include "axom/slic/core/SimpleLogger.hpp"
#include "axom/sidre.hpp"

#include "axom/inlet/Inlet.hpp"
#include "axom/inlet/tests/inlet_test_utils.hpp"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <array>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>

using axom::Path;
using axom::inlet::Inlet;
using axom::inlet::InletType;
using axom::inlet::VariantKey;
using axom::inlet::VerificationError;

using ::testing::Contains;
using ::testing::Truly;

template <typename InletReader>
Inlet createBasicInlet(const std::string& luaString, bool enableDocs = true)
{
  std::unique_ptr<InletReader> reader(new InletReader());
  reader->parseString(axom::inlet::detail::fromLuaTo<InletReader>(luaString));
  return Inlet(std::move(reader), enableDocs);
}

struct Foo
{
  bool bar;
  bool baz;

  bool operator==(const Foo& other) const { return bar == other.bar && baz == other.baz; }
};

template <>
struct FromInlet<Foo>
{
  Foo operator()(const axom::inlet::Container& base)
  {
    Foo f {base["bar"], base["baz"]};
    return f;
  }
};

template <typename InletReader>
class inlet_object : public ::testing::Test
{ };

TYPED_TEST_SUITE(inlet_object, axom::inlet::detail::ReaderTypes);

TYPED_TEST(inlet_object, simple_struct_by_value)
{
  std::string testString = "foo = { bar = true; baz = false }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  // Define schema
  inlet.addBool("foo/bar", "bar's description");
  inlet.addBool("foo/baz", "baz's description");

  Foo foo;

  foo = inlet.get<Foo>("foo");
  EXPECT_TRUE(foo.bar);
  EXPECT_FALSE(foo.baz);
}

TYPED_TEST(inlet_object, simple_struct_verify_pass)
{
  std::string testString = "";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& foo_table = inlet.addStruct("foo");
  foo_table.addBool("bar", "bar's description").required(true);
  foo_table.addBool("baz", "baz's description").required(true);

  // Should pass verification as the struct is not present
  EXPECT_TRUE(inlet.verify());
}

TYPED_TEST(inlet_object, simple_struct_verify_fail)
{
  std::string testString = "foo = { baz = true }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& foo_table = inlet.addStruct("foo");
  foo_table.addBool("bar", "bar's description").required(true);
  foo_table.addBool("baz", "baz's description").required(true);

  // It should fail because a) the struct exists, and
  // b) one of the required fields is not present
  // Note that the only reason the struct is marked as present
  // is because at least one of its expected fields is present, i.e.,
  // having a foo = { quux = true } will pass verification
  EXPECT_FALSE(inlet.verify());
}

TYPED_TEST(inlet_object, simple_array_of_struct_by_value)
{
  std::string testString =
    "foo = { [0] = { bar = true; baz = false}, "
    "        [1] = { bar = false; baz = true} }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& arr_container = inlet.addStructArray("foo");

  arr_container.addBool("bar", "bar's description");
  arr_container.addBool("baz", "baz's description");
  std::unordered_map<int, Foo> expected_foos = {{0, {true, false}}, {1, {false, true}}};
  auto foos = inlet["foo"].get<std::unordered_map<int, Foo>>();
  EXPECT_EQ(foos, expected_foos);
}

TYPED_TEST(inlet_object, simple_array_of_struct_implicit_idx)
{
  std::string testString =
    "foo = { { bar = true; baz = false}, "
    "        { bar = false; baz = true} }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& arr_container = inlet.addStructArray("foo");

  arr_container.addBool("bar", "bar's description");
  arr_container.addBool("baz", "baz's description");
  // Lua is 1-indexed
  const int base_idx = TypeParam::baseIndex;
  std::unordered_map<int, Foo> expected_foos = {{base_idx, {true, false}},
                                                {base_idx + 1, {false, true}}};
  auto foos = inlet["foo"].get<std::unordered_map<int, Foo>>();
  EXPECT_EQ(foos, expected_foos);
}

TYPED_TEST(inlet_object, simple_array_of_struct_verify_optional)
{
  std::string testString =
    "foo = { [4] = { bar = true;}, "
    "        [7] = { bar = false;} }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& arr_container = inlet.addStructArray("foo");

  arr_container.addBool("bar", "bar's description").required(true);
  arr_container.addBool("baz", "baz's description").required(false);

  EXPECT_TRUE(inlet.verify());
}

TYPED_TEST(inlet_object, simple_array_of_struct_verify_reqd)
{
  std::string testString =
    "foo = { [4] = { bar = true;}, "
    "        [7] = { bar = false;} }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& arr_container = inlet.addStructArray("foo");

  arr_container.addBool("bar", "bar's description").required(true);
  arr_container.addBool("baz", "baz's description").required(true);

  EXPECT_FALSE(inlet.verify());
}

TYPED_TEST(inlet_object, simple_array_of_struct_optional_empty_pass)
{
  std::string testString = "foo = { }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& arr_container = inlet.addStructArray("foo");
  // Even though these are required, the input file should still
  // "verify" as the array is empty
  arr_container.addBool("bar", "bar's description").required(true);
  arr_container.addBool("baz", "baz's description").required(true);

  EXPECT_TRUE(inlet.verify());
}

TYPED_TEST(inlet_object, simple_array_of_struct_required_empty_pass)
{
  std::string testString = "foo = { }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  // Verification should pass because the array exists but is empty
  auto& arr_container = inlet.addStructArray("foo").required(true);
  arr_container.addBool("bar", "bar's description").required(true);
  arr_container.addBool("baz", "baz's description").required(true);

  EXPECT_TRUE(inlet.verify());
}

TYPED_TEST(inlet_object, simple_array_of_primitive_required_empty_pass)
{
  std::string testString = "foo = { }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  // Verification should pass because the array exists but is empty
  inlet.addIntArray("foo").required(true);

  EXPECT_TRUE(inlet.verify());
}

TYPED_TEST(inlet_object, simple_array_of_struct_nonexistent_fail)
{
  std::string testString = "";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  // Verification should fail because the array does not exist
  auto& arr_container = inlet.addStructArray("foo").required(true);
  arr_container.addBool("bar", "bar's description").required(true);
  arr_container.addBool("baz", "baz's description").required(true);

  EXPECT_FALSE(inlet.verify());
}

TYPED_TEST(inlet_object, simple_array_of_primitive_nonexistent_fail)
{
  std::string testString = "";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  // Verification should fail because the array does not exist
  inlet.addIntArray("foo").required(true);

  EXPECT_FALSE(inlet.verify());
}

TYPED_TEST(inlet_object, simple_array_of_struct_verify_lambda)
{
  std::string testString =
    "foo = { [4] = { bar = true;}, "
    "        [7] = { bar = false;} }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& arr_container = inlet.addStructArray("foo");

  arr_container.addBool("bar", "bar's description").required();

  // Ensuring that foo is the element of the collection
  arr_container.registerVerifier(
    [](const axom::inlet::Container& foo) { return foo.contains("bar"); });

  EXPECT_TRUE(inlet.verify());
}

TYPED_TEST(inlet_object, simple_array_of_struct_verify_lambda_pass)
{
  std::string testString =
    "foo = { [4] = { bar = true;}, "
    "        [7] = { baz = false;} }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& arr_container = inlet.addStructArray("foo");

  arr_container.addBool("bar", "bar's description");
  arr_container.addBool("baz", "baz's description");

  // Can specify either "bar" or "baz" but not both
  arr_container.registerVerifier([](const axom::inlet::Container& foo) {
    return !(foo.contains("bar") && foo.contains("baz"));
  });

  EXPECT_TRUE(inlet.verify());
}

TYPED_TEST(inlet_object, simple_array_of_struct_verify_lambda_fail)
{
  std::string testString =
    "foo = { [4] = { bar = true;}, "
    "        [7] = { bar = false; baz = true} }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& arr_container = inlet.addStructArray("foo");

  arr_container.addBool("bar", "bar's description");
  arr_container.addBool("baz", "baz's description");

  // Can specify either "bar" or "baz" but not both
  arr_container.registerVerifier(
    [](const axom::inlet::Container& foo, std::vector<VerificationError>* errors) {
      INLET_VERIFICATION_WARNING("foo", "No bar or baz", errors);
      return !(foo.contains("bar") && foo.contains("baz"));
    });

  std::vector<VerificationError> errors;
  EXPECT_FALSE(inlet.verify(&errors));
  EXPECT_THAT(errors, Contains(Truly([](const VerificationError& err) {
                return err.path == Path("foo") && err.message == "No bar or baz";
              })));
}

struct FooWithArray
{
  std::unordered_map<int, int> arr;
  bool operator==(const FooWithArray& other) const { return arr == other.arr; }
};

template <>
struct FromInlet<FooWithArray>
{
  FooWithArray operator()(const axom::inlet::Container& base)
  {
    FooWithArray f = {base["arr"]};
    return f;
  }
};

TYPED_TEST(inlet_object, array_of_struct_containing_array)
{
  std::string testString =
    "foo = { [0] = { arr = { [0] = 3 }; }, "
    "        [1] = { arr = { [0] = 2 }; } }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& arr_container = inlet.addStructArray("foo");

  arr_container.addIntArray("arr", "arr's description");
  // Contiguous indexing for generality
  std::unordered_map<int, FooWithArray> expected_foos = {{0, {{{0, 3}}}}, {1, {{{0, 2}}}}};
  std::unordered_map<int, FooWithArray> foos_with_arr;
  foos_with_arr = inlet["foo"].get<std::unordered_map<int, FooWithArray>>();
  EXPECT_EQ(foos_with_arr, expected_foos);
}

struct MoveOnlyFoo
{
  MoveOnlyFoo() = delete;
  MoveOnlyFoo(const MoveOnlyFoo&) = delete;
  MoveOnlyFoo(MoveOnlyFoo&&) = default;
  MoveOnlyFoo(bool first, bool second) : bar(first), baz(second) { }
  bool bar;
  bool baz;
};

template <>
struct FromInlet<MoveOnlyFoo>
{
  MoveOnlyFoo operator()(const axom::inlet::Container& base)
  {
    MoveOnlyFoo f(base["bar"], base["baz"]);
    return f;
  }
};

TYPED_TEST(inlet_object, simple_moveonly_struct_by_value)
{
  std::string testString = "foo = { bar = true; baz = false }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  // Define schema
  // Check for existing fields
  inlet.addBool("foo/bar", "bar's description");

  inlet.addBool("foo/baz", "baz's description");

  MoveOnlyFoo foo = inlet.get<MoveOnlyFoo>("foo");
  EXPECT_TRUE(foo.bar);
  EXPECT_FALSE(foo.baz);
}

TYPED_TEST(inlet_object, simple_value_from_bracket)
{
  std::string testString = "foo = true";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  // Define schema
  // Check for existing fields
  inlet.addBool("foo", "foo's description");

  bool foo = inlet["foo"];
  EXPECT_TRUE(foo);
}

TYPED_TEST(inlet_object, simple_struct_from_bracket)
{
  std::string testString = "foo = { bar = true; baz = false }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  // Define schema
  // Check for existing fields
  inlet.addBool("foo/bar", "bar's description");

  inlet.addBool("foo/baz", "baz's description");

  auto foo = inlet["foo"].get<Foo>();
  EXPECT_TRUE(foo.bar);
  EXPECT_FALSE(foo.baz);
}

TYPED_TEST(inlet_object, contains_from_bracket)
{
  std::string testString = "foo = { bar = true; baz = false }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  // Define schema
  // Check for existing fields
  inlet.addBool("foo/bar", "bar's description");

  inlet.addBool("foo/baz", "baz's description");

  EXPECT_TRUE(inlet["foo"].contains("bar"));
  EXPECT_TRUE(inlet["foo"].contains("baz"));
}

TYPED_TEST(inlet_object, array_from_bracket)
{
  std::string testString =
    "luaArrays = { arr1 = { [0] = 4}, "
    "              arr2 = {[0] = true, [1] = false}, "
    "              arr3 = {[0] = 'hello', [1] = 'bye'}, "
    "              arr4 = { [0] = 2.4 } }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  std::unordered_map<int, int> intMap;
  std::unordered_map<int, bool> boolMap;
  std::unordered_map<int, std::string> strMap;
  std::unordered_map<int, double> doubleMap;
  inlet.addIntArray("luaArrays/arr1");
  inlet.addBoolArray("luaArrays/arr2");
  inlet.addStringArray("luaArrays/arr3");
  inlet.addDoubleArray("luaArrays/arr4");

  std::unordered_map<int, int> expectedInts {{0, 4}};
  std::unordered_map<int, bool> expectedBools {{0, true}, {1, false}};
  std::unordered_map<int, std::string> expectedStrs {{0, "hello"}, {1, "bye"}};
  std::unordered_map<int, double> expectedDoubles {{0, 2.4}};

  intMap = inlet["luaArrays/arr1"].get<std::unordered_map<int, int>>();
  EXPECT_EQ(intMap, expectedInts);

  boolMap = inlet["luaArrays/arr2"].get<std::unordered_map<int, bool>>();
  EXPECT_EQ(boolMap, expectedBools);

  strMap = inlet["luaArrays/arr3"].get<std::unordered_map<int, std::string>>();
  EXPECT_EQ(strMap, expectedStrs);

  doubleMap = inlet["luaArrays/arr4"].get<std::unordered_map<int, double>>();
  EXPECT_EQ(doubleMap, expectedDoubles);
}

TYPED_TEST(inlet_object, primitive_type_checks)
{
  std::string testString = " bar = true; baz = 12; quux = 2.5; corge = 'hello' ";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  // Define schema
  // Check for existing fields
  inlet.addBool("bar", "bar's description");

  inlet.addInt("baz", "baz's description");

  inlet.addDouble("quux", "quux's description");

  inlet.addString("corge", "corge's description");

  EXPECT_EQ(inlet["bar"].type(), InletType::Bool);
  bool bar = inlet["bar"];
  EXPECT_EQ(bar, true);

  EXPECT_EQ(inlet["baz"].type(), InletType::Integer);
  int baz = inlet["baz"];
  EXPECT_EQ(baz, 12);

  EXPECT_EQ(inlet["quux"].type(), InletType::Double);
  double quux = inlet["quux"];
  EXPECT_DOUBLE_EQ(quux, 2.5);

  EXPECT_EQ(inlet["corge"].type(), InletType::String);
  std::string corge = inlet["corge"];
  EXPECT_EQ(corge, "hello");
}

TYPED_TEST(inlet_object, composite_type_checks)
{
  std::string testString =
    "luaArrays = { arr1 = { [1] = 4}, "
    "              arr2 = {[4] = true, [8] = false} }; "
    "foo = { bar = true; baz = false }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  // Define schema
  // Check for existing fields
  inlet.addBool("foo/bar", "bar's description");

  inlet.addBool("foo/baz", "baz's description");

  inlet.addIntArray("luaArrays/arr1");
  inlet.addBoolArray("luaArrays/arr2");

  auto arr_container = inlet["luaArrays"];
  // The container containing the two arrays is not an array, but an object
  EXPECT_EQ(arr_container.type(), InletType::Object);

  // But the things it contains are arrays
  EXPECT_EQ(arr_container["arr1"].type(), InletType::Collection);
  EXPECT_EQ(arr_container["arr2"].type(), InletType::Collection);

  auto foo_container = inlet["foo"];
  // Similarly, the container containing the two bools is an object
  EXPECT_EQ(foo_container.type(), InletType::Object);

  // But the things it contains are booleans
  EXPECT_EQ(foo_container["bar"].type(), InletType::Bool);
  EXPECT_EQ(foo_container["baz"].type(), InletType::Bool);
}

TYPED_TEST(inlet_object, implicit_conversion_primitives)
{
  std::string testString =
    " bar = true; baz = 12; quux = 2.5; corge = 'hello'; arr = { [0] = 4, [1] "
    "= 6, [2] = 10}";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  // Define schema
  inlet.addBool("bar", "bar's description");
  inlet.addInt("baz", "baz's description");
  inlet.addDouble("quux", "quux's description");
  inlet.addString("corge", "corge's description");

  inlet.addIntArray("arr");

  // Attempt both construction and assignment
  bool bar = inlet["bar"];
  EXPECT_EQ(bar, true);
  bar = inlet["bar"];
  EXPECT_EQ(bar, true);

  int baz = inlet["baz"];
  EXPECT_EQ(baz, 12);
  baz = inlet["baz"];
  EXPECT_EQ(baz, 12);

  double quux = inlet["quux"];
  EXPECT_DOUBLE_EQ(quux, 2.5);
  quux = inlet["quux"];
  EXPECT_DOUBLE_EQ(quux, 2.5);

  std::string corge = inlet["corge"];
  EXPECT_EQ(corge, "hello");
  corge = inlet["corge"];
  EXPECT_EQ(corge, "hello");

  std::unordered_map<int, int> expected_arr {{0, 4}, {1, 6}, {2, 10}};
  std::unordered_map<int, int> arr = inlet["arr"];
  EXPECT_EQ(arr, expected_arr);
  arr = inlet["arr"];
  EXPECT_EQ(arr, expected_arr);
}

struct BarWithFooWithArray
{
  FooWithArray foo;
  bool operator==(const BarWithFooWithArray& other) const { return foo == other.foo; }
};

template <>
struct FromInlet<BarWithFooWithArray>
{
  BarWithFooWithArray operator()(const axom::inlet::Container& base)
  {
    BarWithFooWithArray b;
    b.foo = base["foo"].get<FooWithArray>();
    return b;
  }
};

TYPED_TEST(inlet_object, nested_array_of_struct_containing_array)
{
  std::string testString =
    "bars = { [0] = { foo = { arr = { [0] = 3 }; } }, "
    "         [1] = { foo = { arr = { [0] = 2 }; } } }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& bar_container = inlet.addStructArray("bars");
  auto& foo_container = bar_container.addStruct("foo");
  foo_container.addIntArray("arr", "arr's description");

  // Contiguous indexing for generality
  std::unordered_map<int, BarWithFooWithArray> expected_bars = {{0, {{{{0, 3}}}}}, {1, {{{{0, 2}}}}}};
  std::unordered_map<int, BarWithFooWithArray> bars_with_foo;
  bars_with_foo = inlet["bars"].get<std::unordered_map<int, BarWithFooWithArray>>();
  EXPECT_EQ(bars_with_foo, expected_bars);
}

struct QuuxWithFooArray
{
  std::unordered_map<int, Foo> arr;
  bool operator==(const QuuxWithFooArray& other) const { return arr == other.arr; }
};

template <>
struct FromInlet<QuuxWithFooArray>
{
  QuuxWithFooArray operator()(const axom::inlet::Container& base)
  {
    QuuxWithFooArray q;
    q.arr = base["arr"].get<std::unordered_map<int, Foo>>();
    return q;
  }
};

TYPED_TEST(inlet_object, nested_array_of_struct)
{
  std::string testString =
    "quux = { [0] = { arr = { [0] = { bar = true; baz = false}, "
    "                         [1] = { bar = false; baz = true} } }, "
    "         [1] = { arr = { [0] = { bar = false; baz = false}, "
    "                         [1] = { bar = true; baz = true} } } }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& quux_container = inlet.addStructArray("quux");
  auto& foo_container = quux_container.addStructArray("arr");

  foo_container.addBool("bar", "bar's description");
  foo_container.addBool("baz", "baz's description");

  // Contiguous indexing for generality
  std::unordered_map<int, QuuxWithFooArray> expected_quuxs = {
    {0, {{{0, {true, false}}, {1, {false, true}}}}},
    {1, {{{0, {false, false}}, {1, {true, true}}}}}};
  std::unordered_map<int, QuuxWithFooArray> quuxs_with_arr;
  quuxs_with_arr = inlet["quux"].get<std::unordered_map<int, QuuxWithFooArray>>();
  EXPECT_EQ(quuxs_with_arr, expected_quuxs);
}

TYPED_TEST(inlet_object, nested_dict_of_array_of_struct)
{
  std::string testString =
    "quux  = { ['first'] = { arr = { [0] = { bar = true; baz = false}, "
    "                                [1] = { bar = false; baz = true} } }, "
    "         ['second'] = { arr = { [0] = { bar = false; baz = false}, "
    "                                [1] = { bar = true; baz = true} } } }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& quux_container = inlet.addStructDictionary("quux");
  auto& foo_container = quux_container.addStructArray("arr");

  foo_container.addBool("bar", "bar's description");
  foo_container.addBool("baz", "baz's description");

  // Contiguous indexing for generality
  std::unordered_map<std::string, QuuxWithFooArray> expected_quuxs = {
    {"first", {{{0, {true, false}}, {1, {false, true}}}}},
    {"second", {{{0, {false, false}}, {1, {true, true}}}}}};
  std::unordered_map<std::string, QuuxWithFooArray> quuxs_with_arr;
  quuxs_with_arr = inlet["quux"].get<std::unordered_map<std::string, QuuxWithFooArray>>();
  EXPECT_EQ(quuxs_with_arr, expected_quuxs);
}

struct CorgeWithQuuxDictionary
{
  std::unordered_map<std::string, QuuxWithFooArray> dict;
  bool operator==(const CorgeWithQuuxDictionary& other) const { return dict == other.dict; }
};

template <>
struct FromInlet<CorgeWithQuuxDictionary>
{
  CorgeWithQuuxDictionary operator()(const axom::inlet::Container& base)
  {
    CorgeWithQuuxDictionary c;
    c.dict = base["dict"].get<std::unordered_map<std::string, QuuxWithFooArray>>();
    return c;
  }
};

TYPED_TEST(inlet_object, nested_array_of_dict_of_array_of_struct)
{
  // clang-format off
  std::string testString =
    "corge = { [0] = { dict = { ['first'] = { arr = { [0] = { bar = true; baz = false}, "
    "                                                 [1] = { bar = false; baz = true} } }, "
    "                          ['second'] = { arr = { [0] = { bar = false; baz = false}, "
    "                                                 [1] = { bar = true; baz = true} } } } }, "
    "          [1] = { dict = { ['third'] = { arr = { [0] = { bar = true; baz = false}, "
    "                                                 [1] = { bar = false; baz = true} } }, "
    "                          ['fourth'] = { arr = { [0] = { bar = false; baz = false}, "
    "                                                 [1] = { bar = true; baz = true} } } } } }";
  // clang-format on
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& corge_container = inlet.addStructArray("corge");
  auto& quux_container = corge_container.addStructDictionary("dict");
  auto& foo_container = quux_container.addStructArray("arr");

  foo_container.addBool("bar", "bar's description");
  foo_container.addBool("baz", "baz's description");

  // Contiguous indexing for generality
  std::unordered_map<int, CorgeWithQuuxDictionary> expected_corges = {
    {0,
     {{{"first", {{{0, {true, false}}, {1, {false, true}}}}},
       {"second", {{{0, {false, false}}, {1, {true, true}}}}}}}},
    {1,
     {{{"third", {{{0, {true, false}}, {1, {false, true}}}}},
       {"fourth", {{{0, {false, false}}, {1, {true, true}}}}}}}}};
  std::unordered_map<int, CorgeWithQuuxDictionary> corges_with_dict;
  corges_with_dict = inlet["corge"].get<std::unordered_map<int, CorgeWithQuuxDictionary>>();
  EXPECT_EQ(corges_with_dict, expected_corges);
}

struct QuuxWithFooWithArray
{
  std::unordered_map<int, FooWithArray> arr;
  bool operator==(const QuuxWithFooWithArray& other) const { return arr == other.arr; }
};

template <>
struct FromInlet<QuuxWithFooWithArray>
{
  QuuxWithFooWithArray operator()(const axom::inlet::Container& base)
  {
    QuuxWithFooWithArray q;
    q.arr = base["outer_arr"].get<std::unordered_map<int, FooWithArray>>();
    return q;
  }
};

TYPED_TEST(inlet_object, nested_dict_of_array_of_struct_with_array)
{
  std::string testString =
    "quux = { ['first'] = { outer_arr = {  [0] = { arr = { [0] = 1 } }, "
    "                                      [1] = { arr = { [0] = 2 } } } }, "
    "         ['second'] = { outer_arr = { [0] = { arr = { [0] = 3 } }, "
    "                                      [1] = { arr = { [0] = 4 } } } } }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& quux_container = inlet.addStructDictionary("quux");
  auto& foo_container = quux_container.addStructArray("outer_arr");

  foo_container.addIntArray("arr", "arr's description");

  // Contiguous indexing for generality
  std::unordered_map<std::string, QuuxWithFooWithArray> expected_quuxs = {
    {"first", {{{0, {{{0, 1}}}}, {1, {{{0, 2}}}}}}},
    {"second", {{{0, {{{0, 3}}}}, {1, {{{0, 4}}}}}}}};
  std::unordered_map<std::string, QuuxWithFooWithArray> quuxs_with_arr;
  quuxs_with_arr = inlet["quux"].get<std::unordered_map<std::string, QuuxWithFooWithArray>>();
  EXPECT_EQ(quuxs_with_arr, expected_quuxs);
}

struct QuuxWithSingleFoo
{
  Foo foo;
  bool operator==(const QuuxWithSingleFoo& other) const { return foo == other.foo; }
};

template <>
struct FromInlet<QuuxWithSingleFoo>
{
  QuuxWithSingleFoo operator()(const axom::inlet::Container& base)
  {
    QuuxWithSingleFoo q;
    q.foo = base["foo"].get<Foo>();
    return q;
  }
};

TYPED_TEST(inlet_object, nested_array_of_nested_structs)
{
  std::string testString =
    "quux = { [0] = { foo = { bar = true; baz = false } }, "
    "         [1] = { foo = { bar = false; baz = true } } }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& quux_schema = inlet.addStructArray("quux");
  auto& foo_schema = quux_schema.addStruct("foo");

  foo_schema.addBool("bar", "bar's description");
  foo_schema.addBool("baz", "baz's description");

  // Contiguous indexing for generality
  std::unordered_map<int, QuuxWithSingleFoo> expected_quuxs = {{0, {true, false}},
                                                               {1, {false, true}}};
  std::unordered_map<int, QuuxWithSingleFoo> quuxs_with_foo;
  quuxs_with_foo = inlet["quux"].get<std::unordered_map<int, QuuxWithSingleFoo>>();
  EXPECT_EQ(quuxs_with_foo, expected_quuxs);
}

TYPED_TEST(inlet_object, nested_array_of_struct_verify_pass)
{
  std::string testString = "quux = { }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& quux_schema = inlet.addStructArray("quux");

  // Simple required field
  quux_schema.addString("corge").required();

  // The struct itself is required, but its fields are not
  auto& foo_schema = quux_schema.addStruct("foo").required();
  foo_schema.addBool("bar", "bar's description");
  foo_schema.addBool("baz", "baz's description");

  // The struct itself is not required, but its fields are
  auto& grault_schema = quux_schema.addStruct("grault");
  grault_schema.addDouble("thud").required();

  EXPECT_TRUE(inlet.verify());
}

TYPED_TEST(inlet_object, primitive_arrays_as_std_vector)
{
  std::string testString = " arr = { [0] = 4, [1] = 6, [2] = 10}";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  // Define schema
  inlet.addIntArray("arr");

  // Attempt both construction and assignment
  std::vector<int> expected_arr {4, 6, 10};
  std::vector<int> arr = inlet["arr"];
  EXPECT_EQ(arr, expected_arr);
  arr = inlet["arr"];
  EXPECT_EQ(arr, expected_arr);

  // Also possible to retrieve with indices
  std::unordered_map<int, int> expected_arr_w_indices {{0, 4}, {1, 6}, {2, 10}};
  std::unordered_map<int, int> arr_w_indices = inlet["arr"];
  EXPECT_EQ(arr_w_indices, expected_arr_w_indices);
}

TYPED_TEST(inlet_object, primitive_arrays_as_std_vector_wrong_type)
{
  std::string testString = " arr = { [0] = 'a', [1] = 'b', [2] = 'c'}";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  // Define schema
  inlet.addIntArray("arr");

  // Even though the array was not required, the presence of string elements
  // should trigger a verification failure
  EXPECT_FALSE(inlet.verify());
}

TYPED_TEST(inlet_object, primitive_arrays_as_std_vector_wrong_type_reqd_fail)
{
  std::string testString = " arr = { [0] = 'a', [1] = 'b', [2] = 'c'}";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  // Define schema
  inlet.addIntArray("arr").required();

  // The array was empty (same as if it didn't exist), but marked as required
  EXPECT_FALSE(inlet.verify());
}

TYPED_TEST(inlet_object, primitive_arrays_as_std_vector_mixed_type)
{
  std::string testString = " arr = { [0] = 4, [1] = 6, [2] = 'a', [3] = 'b'}";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  // Define schema
  inlet.addIntArray("arr");

  // Even though the array was not required and some double elements are present,
  // the presence of string elements should trigger a verification failure
  EXPECT_FALSE(inlet.verify());
}

TYPED_TEST(inlet_object, struct_arrays_as_std_vector)
{
  std::string testString =
    "foo = { [0] = { bar = true; baz = false}, "
    "        [1] = { bar = false; baz = true} }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& arr_container = inlet.addStructArray("foo");

  arr_container.addBool("bar", "bar's description");
  arr_container.addBool("baz", "baz's description");
  std::vector<Foo> expected_foos = {{true, false}, {false, true}};
  auto foos = inlet["foo"].get<std::vector<Foo>>();
  EXPECT_EQ(foos, expected_foos);
}

TYPED_TEST(inlet_object, default_scalar_user_provided)
{
  std::string testString = " ";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& scalar = inlet.addInt("foo").defaultValue(2);
  // The field itself exists but was not provided by the user
  auto& field = static_cast<axom::inlet::Field&>(scalar);
  // "foo" exists because it has the default value of 2
  EXPECT_TRUE(field.exists());
  // It can also be said that the hierarchy "contains" foo - again because it has a default
  EXPECT_TRUE(inlet.contains("foo"));
  // But because the input string was empty "foo" was not provided by the user
  EXPECT_FALSE(field.isUserProvided());
  EXPECT_FALSE(inlet.isUserProvided("foo"));

  // ...but it should still be possible to retrieve the default
  const int foo = inlet["foo"];
  EXPECT_EQ(foo, 2);

  // and it should not impede verification
  EXPECT_TRUE(inlet.verify());
}

TYPED_TEST(inlet_object, default_struct_field_user_provided)
{
  std::string testString = " ";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& foo_container = inlet.addStruct("foo");
  foo_container.addBool("bar", "bar's description").defaultValue(true);
  foo_container.addBool("baz", "baz's description").defaultValue(false);

  // The container itself exists (has defaults) but was not provided by the user
  EXPECT_TRUE(foo_container.exists());
  EXPECT_TRUE(inlet.contains("foo"));
  // But no struct elements were provided by the user b/c empty input string
  EXPECT_FALSE(foo_container.isUserProvided());
  EXPECT_FALSE(inlet.isUserProvided("foo"));

  // ...but it should still be possible to retrieve the default
  const Foo expected_foo {true, false};
  const auto foo = inlet["foo"].get<Foo>();
  EXPECT_EQ(foo, expected_foo);

  // and it should not impede verification
  EXPECT_TRUE(inlet.verify());
}

TYPED_TEST(inlet_object, default_struct_field_user_provided_reqd)
{
  std::string testString = " ";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& foo_container = inlet.addStruct("foo").required();
  foo_container.addBool("bar", "bar's description").defaultValue(true);
  foo_container.addBool("baz", "baz's description").defaultValue(false);

  // The container itself exists (has defaults) but was not provided by the user
  EXPECT_TRUE(foo_container.exists());
  EXPECT_TRUE(inlet.contains("foo"));
  // But no struct elements were provided by the user b/c empty input string
  EXPECT_FALSE(foo_container.isUserProvided());
  EXPECT_FALSE(inlet.isUserProvided("foo"));

  // ...but it should still be possible to retrieve the default
  const Foo expected_foo {true, false};
  const auto foo = inlet["foo"].get<Foo>();
  EXPECT_EQ(foo, expected_foo);

  // and it should not impede verification
  // EXPECT_TRUE(inlet.verify()); // fails here...
  // FIXME: This one is a bit of a degenerate case where all fields have defaults
  // and the struct is marked as required - need to determine how this can be
  // handled in the general case
}

TYPED_TEST(inlet_object, default_struct_field_marked_true)
{
  std::string testString = "foo = { bar = true }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& foo_container = inlet.addStruct("foo");
  foo_container.addBool("bar", "bar's description").defaultValue(true);
  foo_container.addBool("baz", "baz's description").defaultValue(false);

  // The container itself exists and was provided by the user
  EXPECT_TRUE(foo_container.exists());
  EXPECT_TRUE(inlet.contains("foo"));
  EXPECT_TRUE(foo_container.isUserProvided());
  EXPECT_TRUE(inlet.isUserProvided("foo"));

  // If we examine the individual struct fields, "bar" will exist and be user-provided
  EXPECT_TRUE(inlet.contains("foo/bar"));
  EXPECT_TRUE(foo_container.contains("bar"));
  EXPECT_TRUE(foo_container.isUserProvided("bar"));
  EXPECT_TRUE(inlet.isUserProvided("foo/bar"));
  // ...but "baz" will exist and not be user-provided
  EXPECT_TRUE(inlet.contains("foo/baz"));
  EXPECT_TRUE(foo_container.contains("baz"));
  EXPECT_FALSE(foo_container.isUserProvided("baz"));
  EXPECT_FALSE(inlet.isUserProvided("foo/baz"));
  // ... while a nonexistent "quux" will be neither
  EXPECT_FALSE(inlet.contains("foo/quux"));
  EXPECT_FALSE(foo_container.contains("quux"));
  EXPECT_FALSE(foo_container.isUserProvided("quux"));
  EXPECT_FALSE(inlet.isUserProvided("foo/quux"));

  // ...and it should still be possible to retrieve the full struct
  const Foo expected_foo {true, false};
  const auto foo = inlet["foo"].get<Foo>();
  EXPECT_EQ(foo, expected_foo);

  // and it verification should succeed
  EXPECT_TRUE(inlet.verify());
}

TYPED_TEST(inlet_object, basic_unused_names)
{
  std::string testString =
    "foo = { [0] = { bar = true; baz = false}, "
    "        [1] = { bar = false; baz = true} }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& arr_container = inlet.addStructArray("foo");

  arr_container.addBool("bar", "bar's description");
  // Baz is left unused

  // Should still verify - unexpected fields do not mean invalid
  EXPECT_TRUE(inlet.verify());

  std::vector<std::string> expected_unused {"foo/0/baz", "foo/1/baz"};
  EXPECT_EQ(expected_unused, inlet.unexpectedNames());
}

TYPED_TEST(inlet_object, basic_unused_names_strict)
{
  std::string testString =
    "foo = { [0] = { bar = true; baz = false}, "
    "        [1] = { bar = false; baz = true} }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& arr_container = inlet.addStructArray("foo").strict();

  arr_container.addBool("bar", "bar's description");
  // Baz is left unused

  // Should fail verification, marked as strict
  EXPECT_FALSE(inlet.verify());
}

TYPED_TEST(inlet_object, basic_unused_names_substring)
{
  std::string testString =
    "foo = { [0] = { bar = true; barz = false}, "
    "        [1] = { bar = false; barz = true} }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& arr_container = inlet.addStructArray("foo");

  arr_container.addBool("barz", "barz's description");
  // Baz is left unused

  // Should still verify - unexpected fields do not mean invalid
  EXPECT_TRUE(inlet.verify());

  // Check to make sure that a naive substring is not used and that checks are path-aware
  std::vector<std::string> expected_unused {"foo/0/bar", "foo/1/bar"};
  EXPECT_EQ(expected_unused, inlet.unexpectedNames());
}

TYPED_TEST(inlet_object, basic_unused_names_substring_strict)
{
  std::string testString =
    "foo = { [0] = { bar = true; barz = false}, "
    "        [1] = { bar = false; barz = true} }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& arr_container = inlet.addStructArray("foo").strict();

  arr_container.addBool("barz", "barz's description");
  // Baz is left unused

  // Should fail verification, marked as strict
  EXPECT_FALSE(inlet.verify());
}

template <typename InletReader>
class inlet_object_dict : public ::testing::Test
{ };

TYPED_TEST_SUITE(inlet_object_dict, axom::inlet::detail::ReaderTypes);

TYPED_TEST(inlet_object_dict, basic_dicts)
{
  std::string testString = "foo = { ['key1'] = 4, ['key3'] = 6, ['key2'] = 10}";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addIntDictionary("foo", "foo's description");
  std::unordered_map<std::string, int> dict = inlet["foo"];
  std::unordered_map<std::string, int> correct = {{"key1", 4}, {"key3", 6}, {"key2", 10}};
  EXPECT_EQ(dict, correct);
}

TYPED_TEST(inlet_object_dict, simple_dict_of_struct_by_value)
{
  std::string testString =
    "foo = { ['key1'] = { bar = true; baz = false}, "
    "        ['key2'] = { bar = false; baz = true} }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& dict_container = inlet.addStructDictionary("foo");

  dict_container.addBool("bar", "bar's description");
  dict_container.addBool("baz", "baz's description");
  std::unordered_map<std::string, Foo> expected_foos = {{"key1", {true, false}},
                                                        {"key2", {false, true}}};
  std::unordered_map<std::string, Foo> foos;
  foos = inlet["foo"].get<std::unordered_map<std::string, Foo>>();
  EXPECT_EQ(foos, expected_foos);
}

struct FooWithDict
{
  std::unordered_map<std::string, int> arr;
  bool operator==(const FooWithDict& other) const { return arr == other.arr; }
};

template <>
struct FromInlet<FooWithDict>
{
  FooWithDict operator()(const axom::inlet::Container& base)
  {
    FooWithDict f = {base["arr"]};
    return f;
  }
};

TYPED_TEST(inlet_object_dict, dict_of_struct_containing_dict)
{
  std::string testString =
    "foo = { ['key3'] = { arr = { ['key1'] = 3 }; }, "
    "        ['key4'] = { arr = { ['key2'] = 2 }; } }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& dict_container = inlet.addStructDictionary("foo");

  dict_container.addIntDictionary("arr", "arr's description");
  std::unordered_map<std::string, FooWithDict> expected_foos = {{"key3", {{{"key1", 3}}}},
                                                                {"key4", {{{"key2", 2}}}}};
  std::unordered_map<std::string, FooWithDict> foos_with_dict;
  foos_with_dict = inlet["foo"].get<std::unordered_map<std::string, FooWithDict>>();
  EXPECT_EQ(foos_with_dict, expected_foos);
}

TYPED_TEST(inlet_object_dict, dict_of_struct_containing_array)
{
  std::string testString =
    "foo = { ['key3'] = { arr = { [0] = 3 }; }, "
    "        ['key4'] = { arr = { [0] = 2 }; } }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& dict_container = inlet.addStructDictionary("foo");

  dict_container.addIntArray("arr", "arr's description");
  std::unordered_map<std::string, FooWithArray> expected_foos = {{"key3", {{{0, 3}}}},
                                                                 {"key4", {{{0, 2}}}}};
  std::unordered_map<std::string, FooWithArray> foos_with_array;
  foos_with_array = inlet["foo"].get<std::unordered_map<std::string, FooWithArray>>();
  EXPECT_EQ(foos_with_array, expected_foos);
}

TYPED_TEST(inlet_object_dict, array_of_struct_containing_dict)
{
  std::string testString =
    "foo = { [0] = { arr = { ['key1'] = 3 }; }, "
    "        [1] = { arr = { ['key2'] = 2 }; } }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& arr_container = inlet.addStructArray("foo");

  arr_container.addIntDictionary("arr", "arr's description");
  std::unordered_map<int, FooWithDict> expected_foos = {{0, {{{"key1", 3}}}}, {1, {{{"key2", 2}}}}};
  std::unordered_map<int, FooWithDict> foos_with_dict;
  foos_with_dict = inlet["foo"].get<std::unordered_map<int, FooWithDict>>();
  EXPECT_EQ(foos_with_dict, expected_foos);
}

/*
FIXME: These are currently error conditions.  If these should be supported
or handled differently these tests can be re-enabled.

TEST(inlet_dict, mixed_keys_primitive_duplicated)
{
  std::string testString = "foo = { ['1'] = 4, [1] = 6 }";
  auto inlet = createBasicInlet(&ds, testString);

  inlet.addIntDictionary("foo", "foo's description");
  std::unordered_map<VariantKey, int> dict = inlet["foo"];
  std::unordered_map<VariantKey, int> correct_dict = {{"1", 4}, {1, 6}};
  EXPECT_EQ(dict, correct_dict);
}

TEST(inlet_dict, key_with_slash)
{
  std::string testString =
    "foo = { ['key1/subkey1'] = 4, ['key3'] = 6, ['key2'] = 10}";
  Inlet inlet = createBasicInlet(&ds, testString);

  inlet.addIntDictionary("foo", "foo's description");
  std::unordered_map<std::string, int> dict = inlet["foo"];
  std::unordered_map<std::string, int> correct = {{"key1/subkey1", 4},
                                                  {"key3", 6},
                                                  {"key2", 10}};
  EXPECT_EQ(dict, correct);
}
*/

// Noncontiguous array tests that are only valid in Lua
#ifdef AXOM_USE_SOL

TEST(inlet_object_lua, array_of_struct_containing_array)
{
  std::string testString =
    "foo = { [4] = { arr = { [1] = 3 }; }, "
    "        [7] = { arr = { [6] = 2 }; } }";
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(testString);

  auto& arr_container = inlet.addStructArray("foo");

  arr_container.addIntArray("arr", "arr's description");
  std::unordered_map<int, FooWithArray> expected_foos = {{4, {{{1, 3}}}}, {7, {{{6, 2}}}}};
  std::unordered_map<int, FooWithArray> foos_with_arr;
  foos_with_arr = inlet["foo"].get<std::unordered_map<int, FooWithArray>>();
  EXPECT_EQ(foos_with_arr, expected_foos);
}

TEST(inlet_object_lua, array_from_bracket)
{
  std::string testString =
    "luaArrays = { arr1 = { [1] = 4}, "
    "              arr2 = {[4] = true, [8] = false}, "
    "              arr3 = {[33] = 'hello', [2] = 'bye'}, "
    "              arr4 = { [12] = 2.4 } }";
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(testString);

  std::unordered_map<int, int> intMap;
  std::unordered_map<int, bool> boolMap;
  std::unordered_map<int, std::string> strMap;
  std::unordered_map<int, double> doubleMap;
  inlet.addIntArray("luaArrays/arr1");
  inlet.addBoolArray("luaArrays/arr2");
  inlet.addStringArray("luaArrays/arr3");
  inlet.addDoubleArray("luaArrays/arr4");

  std::unordered_map<int, int> expectedInts {{1, 4}};
  std::unordered_map<int, bool> expectedBools {{4, true}, {8, false}};
  std::unordered_map<int, std::string> expectedStrs {{33, "hello"}, {2, "bye"}};
  std::unordered_map<int, double> expectedDoubles {{12, 2.4}};

  intMap = inlet["luaArrays/arr1"].get<std::unordered_map<int, int>>();
  EXPECT_EQ(intMap, expectedInts);

  boolMap = inlet["luaArrays/arr2"].get<std::unordered_map<int, bool>>();
  EXPECT_EQ(boolMap, expectedBools);

  strMap = inlet["luaArrays/arr3"].get<std::unordered_map<int, std::string>>();
  EXPECT_EQ(strMap, expectedStrs);

  doubleMap = inlet["luaArrays/arr4"].get<std::unordered_map<int, double>>();
  EXPECT_EQ(doubleMap, expectedDoubles);
}

TEST(inlet_object_lua, primitive_arrays_as_std_vector_discontiguous)
{
  std::string testString = " arr = { [0] = 4, [8] = 6, [12] = 10}";
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(testString);

  // Define schema
  inlet.addIntArray("arr");

  // Attempt both construction and assignment
  std::vector<int> expected_arr {4, 6, 10};
  std::vector<int> arr = inlet["arr"];
  EXPECT_EQ(arr, expected_arr);
  arr = inlet["arr"];
  EXPECT_EQ(arr, expected_arr);
}

TEST(inlet_object_lua, struct_arrays_as_std_vector_discontiguous)
{
  std::string testString =
    "foo = { [6] = { bar = true; baz = false}, "
    "        [11] = { bar = false; baz = true} }";
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(testString);

  auto& arr_container = inlet.addStructArray("foo");

  arr_container.addBool("bar", "bar's description");
  arr_container.addBool("baz", "baz's description");
  std::vector<Foo> expected_foos = {{true, false}, {false, true}};
  auto foos = inlet["foo"].get<std::vector<Foo>>();
  EXPECT_EQ(foos, expected_foos);
}

TEST(inlet_object_lua, primitive_arrays_as_std_vector_implicit_idx)
{
  std::string testString = " arr = { 4, 6, 10}";
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(testString);

  // Define schema
  inlet.addIntArray("arr");

  // Attempt both construction and assignment
  std::vector<int> expected_arr {4, 6, 10};
  std::vector<int> arr = inlet["arr"];
  EXPECT_EQ(arr, expected_arr);
  arr = inlet["arr"];
  EXPECT_EQ(arr, expected_arr);
}

TEST(inlet_object_lua, struct_arrays_as_std_vector_implicit_idx)
{
  std::string testString =
    "foo = { { bar = true; baz = false}, "
    "        { bar = false; baz = true} }";
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(testString);

  auto& arr_container = inlet.addStructArray("foo");

  arr_container.addBool("bar", "bar's description");
  arr_container.addBool("baz", "baz's description");
  std::vector<Foo> expected_foos = {{true, false}, {false, true}};
  auto foos = inlet["foo"].get<std::vector<Foo>>();
  EXPECT_EQ(foos, expected_foos);
}

TEST(inlet_object_lua_dict, dict_of_struct_containing_array)
{
  std::string testString =
    "foo = { ['key3'] = { arr = { [1] = 3 }; }, "
    "        ['key4'] = { arr = { [6] = 2 }; } }";
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(testString);

  auto& dict_container = inlet.addStructDictionary("foo");

  dict_container.addIntArray("arr", "arr's description");
  std::unordered_map<std::string, FooWithArray> expected_foos = {{"key3", {{{1, 3}}}},
                                                                 {"key4", {{{6, 2}}}}};
  std::unordered_map<std::string, FooWithArray> foos_with_array;
  foos_with_array = inlet["foo"].get<std::unordered_map<std::string, FooWithArray>>();
  EXPECT_EQ(foos_with_array, expected_foos);
}

TEST(inlet_object_lua_dict, array_of_struct_containing_dict)
{
  std::string testString =
    "foo = { [7] = { arr = { ['key1'] = 3 }; }, "
    "        [4] = { arr = { ['key2'] = 2 }; } }";
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(testString);

  auto& arr_container = inlet.addStructArray("foo");

  arr_container.addIntDictionary("arr", "arr's description");
  std::unordered_map<int, FooWithDict> expected_foos = {{7, {{{"key1", 3}}}}, {4, {{{"key2", 2}}}}};
  std::unordered_map<int, FooWithDict> foos_with_dict;
  foos_with_dict = inlet["foo"].get<std::unordered_map<int, FooWithDict>>();
  EXPECT_EQ(foos_with_dict, expected_foos);
}

TEST(inlet_object_lua_dict, mixed_keys_primitive)
{
  std::string testString = "foo = { ['key1'] = 4, [1] = 6 }";
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(testString);

  inlet.addIntDictionary("foo", "foo's description");
  std::unordered_map<VariantKey, int> dict = inlet["foo"];
  std::unordered_map<VariantKey, int> correct_dict = {{"key1", 4}, {1, 6}};
  EXPECT_EQ(dict, correct_dict);
}

TEST(inlet_object_lua_dict, mixed_keys_primitive_ignore_string_only)
{
  std::string testString = "foo = { ['key1'] = 4, [1] = 6 }";
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(testString);

  inlet.addIntDictionary("foo", "foo's description");
  std::unordered_map<std::string, int> dict = inlet["foo"];
  std::unordered_map<std::string, int> correct_dict = {{"key1", 4}};
  EXPECT_EQ(dict, correct_dict);
}

TEST(inlet_object_lua_dict, mixed_keys_primitive_ignore_int_only)
{
  std::string testString = "foo = { ['key1'] = 4, [1] = 6 }";
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(testString);

  inlet.addIntArray("foo", "foo's description");
  std::unordered_map<int, int> array = inlet["foo"];
  std::unordered_map<int, int> correct_array = {{1, 6}};
  EXPECT_EQ(array, correct_array);
}

TEST(inlet_object_lua_dict, mixed_keys_object)
{
  std::string testString =
    "foo = { ['key1'] = { bar = true; baz = false}, "
    "        [1] = { bar = false; baz = true} }";
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(testString);

  auto& dict_container = inlet.addStructDictionary("foo");

  dict_container.addBool("bar", "bar's description");
  dict_container.addBool("baz", "baz's description");
  std::unordered_map<VariantKey, Foo> expected_foos = {{"key1", {true, false}}, {1, {false, true}}};
  std::unordered_map<VariantKey, Foo> foos;
  foos = inlet["foo"].get<std::unordered_map<VariantKey, Foo>>();
  EXPECT_EQ(foos, expected_foos);
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
