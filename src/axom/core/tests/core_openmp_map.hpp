// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/Map.hpp"
#include "gtest/gtest.h" /* for TEST and EXPECT_* macros */
#include <string>

namespace axom
{
namespace internal
{
template <typename Key, typename T, typename Hash, typename Policy>
experimental::Map<Key, T, Hash, Policy> init(int N, int len)
{
  experimental::Map<Key, T, Hash, Policy> test(N, len);
  EXPECT_EQ(N * len, test.max_size());
  EXPECT_EQ(0, test.size());
  EXPECT_EQ(true, test.empty());
  EXPECT_EQ(len, test.bucket_size());
  EXPECT_EQ(N, test.bucket_count());
  return test;
}

template <typename Key, typename T, typename Hash, typename Policy>
void test_storage(experimental::Map<Key, T, Hash, Policy> &test)
{
  experimental::Map<Key, T, Hash, Policy> *test2 = &test;
  axom::for_all<Policy>(0, test.max_size(), [=](IndexType idx) {
    Key key = idx;
    T value = key * 27;
    auto ret_test = test2->insert(key, value);
    EXPECT_EQ(true, ret_test.second);
  });
  EXPECT_EQ(false, test.empty());
  axom::for_all<Policy>(0, test.max_size(), [=](IndexType idx) {
    Key key = idx;
    EXPECT_EQ(key * 27, test2->find(key).value);
  });
  //This should fail, since we're at capacity.
  auto ret = test.insert(test.max_size(), 900);
  EXPECT_EQ(false, ret.second);
}

template <typename Key, typename T, typename Hash, typename Policy>
void test_subscript(experimental::Map<Key, T, Hash, Policy> &test)
{
  experimental::Map<Key, T, Hash, Policy> *test2 = &test;
  axom::for_all<Policy>(0, test.size(), [=](IndexType idx) {
    Key key = idx;
    EXPECT_EQ(key * 27, (*test2)[key]);
  });
}

template <typename Key, typename T, typename Hash, typename Policy>
void test_insert_assign(experimental::Map<Key, T, Hash, Policy> &test)
{
  experimental::Map<Key, T, Hash, Policy> *test2 = &test;
  axom::for_all<Policy>(0, test.max_size(), [=](IndexType idx) {
    Key key = idx;
    T value = key * 27;
    auto ret_test = test2->insert_or_assign(key, value);
    EXPECT_EQ(true, ret_test.second);
  });

  EXPECT_EQ(false, test.empty());
  axom::for_all<Policy>(0, test.max_size(), [=](IndexType idx) {
    Key key = idx;
    EXPECT_EQ(key * 27, test2->find(key).value);
  });

  axom::for_all<Policy>(0, test.max_size(), [=](IndexType idx) {
    Key key = idx;
    T value = key * 28;
    auto ret_test = test2->insert_or_assign(key, value);
    EXPECT_EQ(false, ret_test.second);
    EXPECT_EQ(ret_test.first->key, key);
  });

  EXPECT_EQ(test.size(), test.max_size());
  axom::for_all<Policy>(0, test.max_size(), [=](IndexType idx) {
    Key key = idx;
    EXPECT_EQ(key * 28, test2->find(key).value);
  });
}

template <typename Key, typename T, typename Hash, typename Policy>
void test_remove(experimental::Map<Key, T, Hash, Policy> &test)
{
  std::size_t to_erase = test.size();
  experimental::Map<Key, T, Hash, Policy> *test2 = &test;

  axom::for_all<Policy>(0, to_erase, [=](IndexType idx) {
    Key key = (Key)idx;
    bool erased = test2->erase(key);
    EXPECT_EQ(erased, true);
    EXPECT_EQ(test2->find(key), test2->end());
  });
  EXPECT_EQ(test.size(), 0);
  test.insert(0, 900);
  auto ret = test.find(0);
  EXPECT_EQ(900, ret.value);
  EXPECT_EQ(test.size(), 1);
}

template <typename Key, typename T, typename Hash, typename Policy>
void test_rehash(experimental::Map<Key, T, Hash, Policy> &test, int num, int fact)
{
  auto original_size = test.size();
  experimental::Map<Key, T, Hash, Policy> *test2 = &test;
  test.rehash(num, fact);

  axom::for_all<Policy>(0, original_size, [=](IndexType idx) {
    Key key = idx;
    EXPECT_EQ(key * 27, test2->find(key).value);
  });

  axom::for_all<Policy>(original_size, test.max_size(), [=](IndexType idx) {
    Key key = idx;
    T value = key * 27;
    auto ret_test = test2->insert(key, value);
    EXPECT_EQ(true, ret_test.second);
  });

  axom::for_all<Policy>(original_size, test.max_size(), [=](IndexType idx) {
    Key key = idx;
    EXPECT_EQ(key * 27, test2->find(key).value);
  });
  auto ret = test.insert(test.max_size(), 900);
  EXPECT_EQ(false, ret.second);
}

}  // namespace internal

#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
TEST(core_map, initialization)
{
  for(int i : {1, 2, 5, 10, 20, 100})
  {
    for(int j : {1, 2, 5, 10})
    {
      experimental::Map<int, int, std::hash<int>, axom::OMP_EXEC> test =
        internal::init<int, int, std::hash<int>, axom::OMP_EXEC>(i, j);
    }
  }
}

TEST(core_map, insertion)
{
  for(int i : {1, 2, 5, 10, 20, 100})
  {
    for(int j : {1, 2, 5, 10})
    {
      experimental::Map<int, int, std::hash<int>, axom::OMP_EXEC> test =
        internal::init<int, int, std::hash<int>, axom::OMP_EXEC>(i, j);
      internal::test_storage<int, int, std::hash<int>, axom::OMP_EXEC>(test);
    }
  }
}

TEST(core_map, insert_or_assign)
{
  for(int i : {1, 2, 5, 10, 20, 100})
  {
    for(int j : {1, 2, 5, 10})
    {
      experimental::Map<int, int, std::hash<int>, axom::OMP_EXEC> test =
        internal::init<int, int, std::hash<int>, axom::OMP_EXEC>(i, j);
      internal::test_insert_assign<int, int, std::hash<int>, axom::OMP_EXEC>(test);
    }
  }
}

TEST(core_map, subscript)
{
  for(int i : {1, 2, 5, 10, 20, 100})
  {
    for(int j : {1, 2, 5, 10})
    {
      experimental::Map<int, int, std::hash<int>, axom::OMP_EXEC> test =
        internal::init<int, int, std::hash<int>, axom::OMP_EXEC>(i, j);
      internal::test_storage<int, int, std::hash<int>, axom::OMP_EXEC>(test);
      internal::test_subscript<int, int, std::hash<int>, axom::OMP_EXEC>(test);
    }
  }
}

TEST(core_map, removal)
{
  for(int i : {1, 2, 5, 10, 20, 100})
  {
    for(int j : {1, 2, 5, 10})
    {
      experimental::Map<int, int, std::hash<int>, axom::OMP_EXEC> test =
        internal::init<int, int, std::hash<int>, axom::OMP_EXEC>(i, j);
      internal::test_storage<int, int, std::hash<int>, axom::OMP_EXEC>(test);
      internal::test_remove<int, int, std::hash<int>, axom::OMP_EXEC>(test);
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
        experimental::Map<int, int, std::hash<int>, axom::OMP_EXEC> test =
          internal::init<int, int, std::hash<int>, axom::OMP_EXEC>(i, j);
        internal::test_storage<int, int, std::hash<int>, axom::OMP_EXEC>(test);
        internal::test_rehash<int, int, std::hash<int>, axom::OMP_EXEC>(test, -1, k);
        EXPECT_EQ(test.max_size(), k * i * j);
        internal::test_remove<int, int, std::hash<int>, axom::OMP_EXEC>(test);
      }
    }
  }
  for(int i : {1, 2, 5, 10, 20, 100})
  {
    for(int j : {1, 2, 5, 10})
    {
      for(int k = 1; k < 4; k++)
      {
        experimental::Map<int, int, std::hash<int>, axom::OMP_EXEC> test =
          internal::init<int, int, std::hash<int>, axom::OMP_EXEC>(i, j);
        internal::test_storage<int, int, std::hash<int>, axom::OMP_EXEC>(test);
        internal::test_rehash<int, int, std::hash<int>, axom::OMP_EXEC>(test,
                                                                        test.max_size() + 20 * k,
                                                                        -1);
        EXPECT_EQ((i * j + 20 * k) * j, test.max_size());
        internal::test_remove<int, int, std::hash<int>, axom::OMP_EXEC>(test);
      }
    }
  }
}
#endif  // AXOM_USE_OPENMP && AXOM_USE_RAJA

} /* namespace axom */
