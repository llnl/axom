// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include <string>
#include <vector>
#include <unordered_map>

#include <iostream>

#include "axom/sidre.hpp"

#include "axom/inlet/Inlet.hpp"

#include "axom/inlet/tests/inlet_test_utils.hpp"

using axom::inlet::Container;
using axom::inlet::Field;
using axom::inlet::Inlet;
using axom::inlet::InletType;
using axom::inlet::Proxy;
using axom::sidre::DataStore;

template <typename InletReader>
Inlet createBasicInlet(DataStore* ds,
                       const std::string& luaString = {},
                       bool enableDocs = true,
                       bool reconstruct = true)
{
  std::unique_ptr<InletReader> reader(new InletReader());
  if(!luaString.empty())
  {
    reader->parseString(axom::inlet::detail::fromLuaTo<InletReader>(luaString));
  }

  return Inlet(std::move(reader), ds->getRoot(), enableDocs, reconstruct);
}

template <typename InletReader>
class inlet_restart : public ::testing::Test
{ };

TYPED_TEST_SUITE(inlet_restart, axom::inlet::detail::ReaderTypes);

TYPED_TEST(inlet_restart, simple_scalars)
{
  std::string testString = "foo = true; bar = false";
  DataStore ds;
  // Separate scope blocks for each Inlet to guarantee separation
  {
    Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

    // Check for existing fields
    inlet.addBool("foo", "foo's description");
    inlet.addBool("bar", "bar's description");
  }

  {
    // No input provided - the datastore should already contain all the data
    Inlet restartInlet = createBasicInlet<TypeParam>(&ds);

    ASSERT_TRUE(restartInlet.contains("foo"));
    ASSERT_TRUE(restartInlet.contains("bar"));

    bool value = false;
    // Check for existing fields
    value = restartInlet["foo"];
    EXPECT_TRUE(value);

    value = restartInlet.get<bool>("foo");
    EXPECT_TRUE(value);

    value = restartInlet["bar"];
    EXPECT_FALSE(value);

    value = restartInlet.get<bool>("bar");
    EXPECT_FALSE(value);
  }
}

TYPED_TEST(inlet_restart, simple_scalars_repeat_schema)
{
  std::string testString = "foo = true; bar = false";
  DataStore ds;

  {
    Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

    // Check for existing fields
    inlet.addBool("foo", "foo's description");
    inlet.addBool("bar", "bar's description");
  }

  {
    // No input provided - the datastore should already contain all the data
    Inlet restartInlet = createBasicInlet<TypeParam>(&ds);

    // Unlike the simple_scalars test, this tests adds foo and bar to the schema,
    // which is redundant but should not be an error
    restartInlet.addBool("foo", "foo's description");
    restartInlet.addBool("bar", "bar's description");

    bool value = false;
    // Check for existing fields
    value = restartInlet["foo"];
    EXPECT_TRUE(value);

    value = restartInlet.get<bool>("foo");
    EXPECT_TRUE(value);

    value = restartInlet["bar"];
    EXPECT_FALSE(value);

    value = restartInlet.get<bool>("bar");
    EXPECT_FALSE(value);
  }
}

struct Foo
{
  bool bar;
  bool baz;
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

TYPED_TEST(inlet_restart, simple_struct)
{
  std::string testString = "foo = { bar = true; baz = false }";
  DataStore ds;
  {
    Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

    // Check for existing fields
    auto& foo = inlet.addStruct("foo", "foo's description");
    foo.addBool("bar", "bar's description");
    foo.addBool("baz", "baz's description");
  }

  {
    // No input provided - the datastore should already contain all the data
    Inlet restartInlet = createBasicInlet<TypeParam>(&ds);

    ASSERT_TRUE(restartInlet.contains("foo"));

    Foo f = restartInlet["foo"].get<Foo>();
    EXPECT_TRUE(f.bar);
    EXPECT_FALSE(f.baz);
  }
}

TYPED_TEST(inlet_restart, simple_struct_repeat_schema)
{
  std::string testString = "foo = { bar = true; baz = false }";
  DataStore ds;
  {
    Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

    // Check for existing fields
    auto& foo = inlet.addStruct("foo", "foo's description");
    foo.addBool("bar", "bar's description");
    foo.addBool("baz", "baz's description");
  }

  {
    // No input provided - the datastore should already contain all the data
    Inlet restartInlet = createBasicInlet<TypeParam>(&ds);

    ASSERT_TRUE(restartInlet.contains("foo"));

    auto& restartFoo = restartInlet.addStruct("foo", "foo's description");
    restartFoo.addBool("bar", "bar's description");
    restartFoo.addBool("baz", "baz's description");

    Foo f = restartInlet["foo"].get<Foo>();
    EXPECT_TRUE(f.bar);
    EXPECT_FALSE(f.baz);
  }
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

TYPED_TEST(inlet_restart, nested_array_of_struct_containing_array)
{
  std::string testString =
    "bars = { [0] = { foo = { arr = { [0] = 3 }; } }, "
    "         [1] = { foo = { arr = { [0] = 2 }; } } }";
  DataStore ds;
  {
    Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

    auto& bar_container = inlet.addStructArray("bars");
    auto& foo_container = bar_container.addStruct("foo");
    foo_container.addIntArray("arr", "arr's description");
  }

  {
    // No input provided - the datastore should already contain all the data
    Inlet restartInlet = createBasicInlet<TypeParam>(&ds);
    // Contiguous indexing for generality
    std::unordered_map<int, BarWithFooWithArray> expected_bars = {{0, {{{{0, 3}}}}},
                                                                  {1, {{{{0, 2}}}}}};
    ASSERT_TRUE(restartInlet.contains("bars"));
    auto bars_with_foo = restartInlet["bars"].get<std::unordered_map<int, BarWithFooWithArray>>();
    EXPECT_EQ(bars_with_foo, expected_bars);
  }
}

TYPED_TEST(inlet_restart, nested_array_of_struct_containing_array_repeat_schema)
{
  std::string testString =
    "bars = { [0] = { foo = { arr = { [0] = 3 }; } }, "
    "         [1] = { foo = { arr = { [0] = 2 }; } } }";
  DataStore ds;
  {
    Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

    auto& bar_container = inlet.addStructArray("bars");
    auto& foo_container = bar_container.addStruct("foo");
    foo_container.addIntArray("arr", "arr's description");
  }

  {
    // No input provided - the datastore should already contain all the data
    Inlet restartInlet = createBasicInlet<TypeParam>(&ds);
    auto& bar_container = restartInlet.addStructArray("bars");
    auto& foo_container = bar_container.addStruct("foo");
    foo_container.addIntArray("arr", "arr's description");
    // Contiguous indexing for generality
    std::unordered_map<int, BarWithFooWithArray> expected_bars = {{0, {{{{0, 3}}}}},
                                                                  {1, {{{{0, 2}}}}}};
    ASSERT_TRUE(restartInlet.contains("bars"));
    auto bars_with_foo = restartInlet["bars"].get<std::unordered_map<int, BarWithFooWithArray>>();
    EXPECT_EQ(bars_with_foo, expected_bars);
  }
}
