// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/LRUCache.hpp"

#include <string>

TEST(core_lru_cache, basic_insert_find_and_clear)
{
  axom::LRUCache<std::string, int> cache(2);
  EXPECT_TRUE(cache.empty());
  EXPECT_EQ(cache.capacity(), 2u);

  cache.insert("one", 1);
  cache.insert("two", 2);

  EXPECT_EQ(cache.size(), 2u);
  EXPECT_TRUE(cache.contains("one"));
  EXPECT_TRUE(cache.contains("two"));

  auto* value = cache.find("one");
  ASSERT_NE(value, nullptr);
  EXPECT_EQ(*value, 1);
  EXPECT_EQ(cache.find("missing"), nullptr);

  cache.clear();
  EXPECT_TRUE(cache.empty());
  EXPECT_FALSE(cache.contains("one"));
}

TEST(core_lru_cache, evicts_least_recently_used_entry)
{
  axom::LRUCache<int, int> cache(3);
  cache.insert(1, 10);
  cache.insert(2, 20);
  cache.insert(3, 30);

  cache.insert(4, 40);

  EXPECT_FALSE(cache.contains(1));
  EXPECT_TRUE(cache.contains(2));
  EXPECT_TRUE(cache.contains(3));
  EXPECT_TRUE(cache.contains(4));
  EXPECT_EQ(cache.size(), 3u);
}

TEST(core_lru_cache, find_promotes_entry_to_most_recently_used)
{
  axom::LRUCache<int, int> cache(3);
  cache.insert(1, 10);
  cache.insert(2, 20);
  cache.insert(3, 30);

  auto* value = cache.find(1);
  ASSERT_NE(value, nullptr);
  EXPECT_EQ(*value, 10);

  cache.insert(4, 40);

  EXPECT_TRUE(cache.contains(1));
  EXPECT_FALSE(cache.contains(2));
  EXPECT_TRUE(cache.contains(3));
  EXPECT_TRUE(cache.contains(4));
}

TEST(core_lru_cache, insert_replaces_and_promotes_existing_entry)
{
  axom::LRUCache<int, int> cache(3);
  cache.insert(1, 10);
  cache.insert(2, 20);
  cache.insert(3, 30);

  cache.insert(1, 100);
  cache.insert(4, 40);

  auto* value = cache.find(1);
  ASSERT_NE(value, nullptr);
  EXPECT_EQ(*value, 100);
  EXPECT_FALSE(cache.contains(2));
  EXPECT_TRUE(cache.contains(3));
  EXPECT_TRUE(cache.contains(4));
  EXPECT_EQ(cache.size(), 3u);
}

TEST(core_lru_cache, peek_does_not_change_eviction_order)
{
  axom::LRUCache<int, int> cache(3);
  cache.insert(1, 10);
  cache.insert(2, 20);
  cache.insert(3, 30);

  const int* value = cache.peek(1);
  ASSERT_NE(value, nullptr);
  EXPECT_EQ(*value, 10);

  cache.insert(4, 40);

  EXPECT_FALSE(cache.contains(1));
  EXPECT_TRUE(cache.contains(2));
  EXPECT_TRUE(cache.contains(3));
  EXPECT_TRUE(cache.contains(4));
}
