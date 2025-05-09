// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/StackArray.hpp"
#include "axom/core/Map.hpp"
#include "gtest/gtest.h"
#include <string>

namespace experimental = axom::experimental;

namespace

{
// Note: axom::Map does not (currently) rehash during insertion operations.
// Most of the tests below assume that a Map with N slots can store all keys
// from 0 to N-1 without a rehash. This is possible with the identity hash function
struct IdentityHash
{
  std::size_t operator()(int i) const { return i; }
};

using TestKey = int;
using TestVal = int;
using TestMap = experimental::Map<TestKey, TestVal, IdentityHash>;

TestMap init(int N, int len)
{
  TestMap test(N, len);
  EXPECT_EQ(N * len, test.max_size());
  EXPECT_EQ(0, test.size());
  EXPECT_EQ(true, test.empty());
  EXPECT_EQ(len, test.bucket_size());
  EXPECT_EQ(N, test.bucket_count());
  return test;
}

void test_storage(TestMap& test)
{
  using Key = TestMap::key_type;
  using T = TestMap::mapped_type;

  EXPECT_TRUE(test.empty());
  for(int i = 0; i < test.max_size(); i++)
  {
    Key key = i;
    T value = key * 27;
    auto ret_test = test.insert(key, value);
    EXPECT_TRUE(ret_test.second);
  }
  EXPECT_FALSE(test.empty());
  for(int i = 0; i < test.max_size(); i++)
  {
    Key key = i;
    EXPECT_EQ(key * 27, test.find(key).value);
  }
  //This should fail, since we're at capacity.
  auto ret = test.insert(test.max_size(), 900);
  EXPECT_FALSE(ret.second);
}

void test_subscript(TestMap& test)
{
  using Key = TestMap::key_type;
  for(int i = 0; i < test.size(); i++)
  {
    Key key = i;
    EXPECT_EQ(key * 27, test[key]);
  }
}

void test_insert_assign(TestMap& test)
{
  using Key = TestMap::key_type;
  using T = TestMap::mapped_type;

  for(int i = 0; i < test.max_size(); i++)
  {
    Key key = i;
    T value = key * 27;
    auto ret_test = test.insert_or_assign(key, value);
    EXPECT_EQ(true, ret_test.second);
  }

  EXPECT_EQ(false, test.empty());
  for(int i = 0; i < test.max_size(); i++)
  {
    Key key = i;
    EXPECT_EQ(key * 27, test.find(key).value);
  }

  for(int i = 0; i < test.max_size(); i++)
  {
    Key key = i;
    T value = key * 28;
    auto ret_test = test.insert_or_assign(key, value);
    EXPECT_EQ(false, ret_test.second);
    EXPECT_EQ(ret_test.first->key, key);
  }
  EXPECT_EQ(test.size(), test.max_size());
  for(int i = 0; i < test.max_size(); i++)
  {
    Key key = i;
    EXPECT_EQ(key * 28, test.find(key).value);
  }
}

void test_remove(TestMap& test)
{
  using Key = TestMap::key_type;

  std::size_t to_erase = test.size();
  for(std::size_t i = 0; i < to_erase; i++)
  {
    Key key = (Key)i;
    bool erased = test.erase(key);
    EXPECT_EQ(erased, true);
    EXPECT_EQ(test.find(key), test.end());
  }
  EXPECT_EQ(test.size(), 0);
  test.insert(0, 900);
  auto ret = test.find(0);
  EXPECT_EQ(900, ret.value);
  EXPECT_EQ(test.size(), 1);
}

void test_rehash(TestMap& test, int num, int fact)
{
  using Key = TestMap::key_type;
  using T = TestMap::mapped_type;

  auto original_size = test.size();
  test.rehash(num, fact);

  for(int i = 0; i < original_size; i++)
  {
    Key key = i;
    EXPECT_EQ(key * 27, test.find(key).value);
  }

  for(int i = original_size; i < test.max_size(); i++)
  {
    Key key = i;
    T value = key * 27;
    auto ret_test = test.insert(key, value);
    EXPECT_EQ(true, ret_test.second);
  }

  for(int i = original_size; i < test.max_size(); i++)
  {
    Key key = i;
    EXPECT_EQ(key * 27, test.find(key).value);
  }
  auto ret = test.insert(test.max_size(), 900);
  EXPECT_EQ(false, ret.second);
}

template <typename Key, typename T>
experimental::axom_map::Node<Key, T> init_node(axom::IndexType next = -1, Key key = 5, T value = 15)
{
  axom::experimental::axom_map::Node<Key, T> n;
  n.next = next;
  n.key = key;
  n.value = value;

  return n;
}

template <typename Key, typename T>
experimental::axom_map::Pair<Key, T> init_pair(experimental::axom_map::Node<Key, T>* node = nullptr,
                                               bool status = false)
{
  axom::experimental::axom_map::Pair<Key, T> p(node, status);

  return p;
}

template <typename Key, typename T>
experimental::axom_map::Bucket<Key, T> init_bucket(int length)
{
  axom::experimental::axom_map::Bucket<Key, T> b(length);

  return b;
}

// The first parameter is the number of buckets;
// the second is a lambda to fill the buckets
template <typename Key, typename T>
experimental::axom_map::Bucket<Key, T> init_filled_bucket(int length, std::function<T(Key)>&& fn)
{
  axom::experimental::axom_map::Bucket<Key, T> b(length);

  for(int i = 0; i < length; ++i)
  {
    b.insert_update(i, fn(i));
  }

  return b;
}

template <typename MapType>
MapType init_filled_map(int numBucket,
                        int bucketLength,
                        std::function<typename MapType::mapped_type(typename MapType::key_type)>&& fn)
{
  using T = typename MapType::mapped_type;

  MapType map(numBucket, bucketLength);

  const int sz = numBucket * bucketLength;
  bool validState = true;
  do
  {
    if(!validState)
    {
      validState = map.rehash();
    }

    if(validState)
    {
      for(int i = 0; i < sz; ++i)
      {
        auto res = map.insert_or_assign(i, fn(i));
        const bool validInsert = (res.second == true || res.first->key == i);
        if(!validInsert)
        {
          validState = false;
          break;
        }
      }
    }
  } while(!validState);

  // Check that the values actually got inserted
  for(int i = 0; i < sz; ++i)
  {
    T exp_value = fn(i);
    T value = map[i];
    EXPECT_EQ(exp_value, value);
  }

  return map;
}

}  // end anonymous namespace

TEST(core_map, initialization)
{
  for(int i : {1, 2, 5, 10, 20, 100})
  {
    for(int j : {1, 2, 5, 10})
    {
      TestMap test = ::init(i, j);
    }
  }
}

TEST(core_map, insertion)
{
  for(int i : {1, 2, 5, 10, 20, 100})
  {
    for(int j : {1, 2, 5, 10})
    {
      TestMap test = ::init(i, j);
      ::test_storage(test);
    }
  }
}

TEST(core_map, insert_or_assign)
{
  for(int i : {1, 2, 5, 10, 20, 100})
  {
    for(int j : {1, 2, 5, 10})
    {
      TestMap test = ::init(i, j);
      ::test_insert_assign(test);
    }
  }
}

TEST(core_map, subscript)
{
  for(int i : {1, 2, 5, 10, 20, 100})
  {
    for(int j : {1, 2, 5, 10})
    {
      TestMap test = ::init(i, j);
      ::test_storage(test);
      ::test_subscript(test);
    }
  }
}

TEST(core_map, removal)
{
  for(int i : {1, 2, 5, 10, 20, 100})
  {
    for(int j : {1, 2, 5, 10})
    {
      TestMap test = ::init(i, j);
      ::test_storage(test);
      ::test_remove(test);
    }
  }
}

TEST(core_map, rehash)
{
  for(int i : {1, 2, 5, 10, 20, 100})
  {
    for(int j : {1, 2, 5, 10})
    {
      for(int k : {2, 4, 8})
      {
        TestMap test = ::init(i, j);
        ::test_storage(test);
        ::test_rehash(test, -1, k);
        EXPECT_EQ(test.max_size(), k * i * j);
        ::test_remove(test);
      }
    }
  }
  for(int i : {1, 2, 5, 10, 20, 100})
  {
    for(int j : {1, 2, 5, 10})
    {
      for(int k = 1; k < 4; k++)
      {
        TestMap test = ::init(i, j);
        ::test_storage(test);
        ::test_rehash(test, test.max_size() + 20 * k, -1);
        EXPECT_EQ((i * j + 20 * k) * j, test.max_size());
        ::test_remove(test);
      }
    }
  }
}

TEST(core_map, node_return_value)
{
  using Key = int;
  using T = int;
  using NodeType = experimental::axom_map::Node<Key, T>;

  axom::IndexType n = -1;
  Key k = 5;
  T v = 15;
  NodeType node = ::init_node<Key, T>(n, k, v);

  EXPECT_EQ(n, node.next);
  EXPECT_EQ(k, node.key);
  EXPECT_EQ(v, node.value);
}

TEST(core_map, pair_return_value)
{
  using Key = int;
  using T = int;
  using NodeType = experimental::axom_map::Node<Key, T>;
  using PairType = experimental::axom_map::Pair<Key, T>;

  {
    axom::IndexType n = -1;
    Key k = 5;
    T v = 15;
    NodeType node = ::init_node<Key, T>(n, k, v);
    PairType pair = ::init_pair<Key, T>(&node, true);

    EXPECT_EQ(node, *pair.first);
    EXPECT_TRUE(pair.second);
  }
  {
    PairType pair = ::init_pair<Key, T>(nullptr, false);

    EXPECT_EQ(nullptr, pair.first);
    EXPECT_FALSE(pair.second);
  }
}

TEST(core_map, bucket_return_value)
{
  using Key = int;
  using T = int;
  using BucketType = experimental::axom_map::Bucket<Key, T>;

  {
    const int length = 5;
    BucketType bucket = ::init_bucket<Key, T>(length);

    EXPECT_EQ(0, bucket.get_size());
    EXPECT_EQ(length, bucket.get_capacity());

    // Note: Bucket's not intended to be used on its own; explicitly free the memory
    axom::deallocate(bucket.m_list);
  }

  {
    const int length = 5;
    auto fn = [](Key i) { return i * i; };
    BucketType bucket = ::init_filled_bucket<Key, T>(length, fn);

    EXPECT_EQ(length, bucket.get_size());
    EXPECT_EQ(length, bucket.get_capacity());

    for(int i = 0; i < length; ++i)
    {
      auto& node = bucket.find(i);
      EXPECT_EQ(i, node.key);
      EXPECT_EQ(fn(i), node.value);
    }

    // Note: Bucket's not intended to be used on its own; explicitly free the memory
    axom::deallocate(bucket.m_list);
  }

  {
    const int length = 15;
    auto fn = [](Key i) { return i * i * i - 1; };
    BucketType bucket = ::init_filled_bucket<Key, T>(length, fn);

    EXPECT_EQ(length, bucket.get_size());
    EXPECT_EQ(length, bucket.get_capacity());

    for(int i = 0; i < length; ++i)
    {
      auto& node = bucket.find(i);
      EXPECT_EQ(i, node.key);
      EXPECT_EQ(fn(i), node.value);
    }

    // Note: Bucket's not intended to be used on its own; explicitly free the memory
    axom::deallocate(bucket.m_list);
  }
}

TEST(core_map, hashmap_return_value)
{
  {
    const int buckets_num = 100;
    const int buckets_len = 3;
    auto fn = [](TestKey i) { return i * i + i; };

    // Use the test map, with the identity hash function
    {
      auto map = ::init_filled_map<TestMap>(buckets_num, buckets_len, fn);

      // Check that all expected values are present
      const int sz = buckets_num;
      for(int i = 0; i < sz; ++i)
      {
        TestKey k = i;
        TestVal exp_value = fn(k);
        TestVal value = map[k];
        EXPECT_EQ(exp_value, value);
      }
    }
    // Use a map with the compiler's default std::hash function
    {
      using DefaultHasherMap = experimental::Map<TestKey, TestVal>;
      auto map = ::init_filled_map<DefaultHasherMap>(buckets_num, buckets_len, fn);

      // Check that all expected values are present
      const int sz = buckets_num;
      for(int i = 0; i < sz; ++i)
      {
        TestKey k = i;
        TestVal exp_value = fn(k);
        TestVal value = map[k];
        EXPECT_EQ(exp_value, value);
      }
    }
  }
}
