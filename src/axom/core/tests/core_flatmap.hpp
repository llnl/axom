// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/FlatMap.hpp"
#include "axom/core/execution/runtime_policy.hpp"

// gtest includes
#include "gtest/gtest.h"

// Unit test for QuadraticProbing
TEST(core_flatmap_unit, quadratic_probing)
{
  axom::detail::flat_map::QuadraticProbing instance {};
  for(int array_size = 4; array_size <= 8192; array_size *= 2)
  {
    // Generate probe indices.
    std::vector<int> probe_index(array_size);
    std::vector<int> probe_index_mod(array_size);
    int curr_offset = 0;
    // Compute probe indices for an array of size N.
    for(int i = 0; i < array_size; i++)
    {
      probe_index[i] = curr_offset;
      probe_index_mod[i] = curr_offset % array_size;
      curr_offset += instance.getNext(i);
    }

    for(int i = 0; i < array_size; i++)
    {
      // Each probe index should match the formula:
      // H(i) = H_0 + i/2 + i^2/2 mod m.
      int expected_probe_index = i * (i + 1) / 2;
      EXPECT_EQ(probe_index[i], expected_probe_index);
      EXPECT_EQ(probe_index_mod[i], expected_probe_index % array_size);
    }

    // Probe indices should be a permutation of the range [0, array_size - 1).
    std::vector<int> expected_permutation(array_size);
    std::iota(expected_permutation.begin(), expected_permutation.end(), 0);

    // Sort the modular probe indexes and compare.
    std::sort(probe_index_mod.begin(), probe_index_mod.end());
    EXPECT_EQ(probe_index_mod, expected_permutation);
  }
}

inline void flatmap_get_value(double key, std::string& out) { out = std::to_string(key); }

inline void flatmap_get_value(int key, std::string& out) { out = std::to_string(key); }

inline void flatmap_get_value(const std::string& key, double& out) { out = std::stod(key); }

inline void flatmap_get_value(const std::string& key, int& out) { out = std::stoi(key); }

template <typename T, typename U>
inline void flatmap_get_value(T key, U& out)
{
  out = key;
}

template <typename FlatMapType>
class core_flatmap : public ::testing::Test
{
public:
  using MapType = FlatMapType;
  using KeyType = typename FlatMapType::key_type;
  using ValueType = typename FlatMapType::mapped_type;

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  using DeviceExec = axom::HIP_EXEC<256>;
#elif defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  using DeviceExec = axom::CUDA_EXEC<256>;
#else
  using DeviceExec = axom::SEQ_EXEC;
#endif
  using HostExec = axom::SEQ_EXEC;

  template <typename T>
  KeyType getKey(T input)
  {
    KeyType key;
    flatmap_get_value(input, key);
    return key;
  }

  template <typename T>
  ValueType getValue(T input)
  {
    ValueType val;
    flatmap_get_value(input, val);
    return val;
  }

  ValueType getDefaultValue() { return ValueType(); }
};

using MyTypes = ::testing::Types<axom::FlatMap<int, double>,
                                 axom::FlatMap<int, std::string>,
                                 axom::FlatMap<std::string, double>,
                                 axom::FlatMap<std::string, std::string>>;

TYPED_TEST_SUITE(core_flatmap, MyTypes);

AXOM_TYPED_TEST(core_flatmap, default_init)
{
  using MapType = typename TestFixture::MapType;
  MapType test_map;
  EXPECT_EQ(0, test_map.size());
  EXPECT_EQ(true, test_map.empty());
}

AXOM_TYPED_TEST(core_flatmap, iterators)
{
  using MapType = typename TestFixture::MapType;
  using MapIterType = typename MapType::iterator;

  MapIterType empty_iterator;
  MapIterType empty_iterator_2;

  // Default-constructed iterators should all be equivalent.
  EXPECT_EQ(empty_iterator, empty_iterator_2);

  MapType test_map;

  const int NUM_ELEMS = 100;

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);

    auto it = test_map.emplace(key, value);

    // Default-constructed iterators should not be equivalent to an actual map's
    // iterators.
    EXPECT_NE(empty_iterator, it.first);
  }

  MapType test_map_2 = test_map;
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);

    auto it = test_map.find(key);
    auto it2 = test_map_2.find(key);

    // Iterators between two different maps should not be equivalent.
    EXPECT_NE(it, it2);
  }
}

AXOM_TYPED_TEST(core_flatmap, prealloc_buckets)
{
  using MapType = typename TestFixture::MapType;

  const std::vector<int> sizes = {100, 400, 1000, 4000, 10000};

  for(int size : sizes)
  {
    MapType test_map(size);
    EXPECT_EQ(0, test_map.size());
    EXPECT_LE(size, test_map.bucket_count() * test_map.max_load_factor());
  }
}

AXOM_TYPED_TEST(core_flatmap, insert_only)
{
  using MapType = typename TestFixture::MapType;
  MapType test_map;

  const int NUM_ELEMS = 100;

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);
    const auto expected_size = i + 1;
    // Initial insertion of a given key should succeed.
    auto initial_insert = test_map.insert({key, value});
    EXPECT_EQ(test_map.size(), expected_size);
    EXPECT_EQ(initial_insert.first, test_map.find(key));
    EXPECT_EQ(value, test_map.at(key));
    EXPECT_TRUE(initial_insert.second);

    int current_bucket_capacity = test_map.bucket_count();

    // Inserting a duplicate key should not change the value.
    auto value_dup = this->getValue(i * 10.0 + 5.0);
    auto duplicate_insert = test_map.insert({key, value_dup});
    EXPECT_EQ(test_map.size(), expected_size);
    EXPECT_EQ(duplicate_insert.first, test_map.find(key));
    EXPECT_EQ(value, test_map.at(key));
    EXPECT_FALSE(duplicate_insert.second);

    // Using operator[] with an already-existing key should return the
    // existing value and not add a value.
    auto value_indexed = test_map[key];
    EXPECT_EQ(value_indexed, value);
    EXPECT_EQ(test_map.size(), expected_size);

    // Check that a rehash didn't occur on the second insertion.
    EXPECT_EQ(duplicate_insert.first, initial_insert.first);
    EXPECT_EQ(current_bucket_capacity, test_map.bucket_count());
  }
}

AXOM_TYPED_TEST(core_flatmap, insert_or_assign)
{
  using MapType = typename TestFixture::MapType;
  MapType test_map;

  const int NUM_ELEMS = 100;

  // Test insert behavior of FlatMap::insert_or_assign.
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);

    auto result = test_map.insert_or_assign(key, value);
    EXPECT_EQ(i + 1, test_map.size());
    EXPECT_EQ(value, test_map.at(key));
    EXPECT_EQ(result.first, test_map.find(key));
    EXPECT_TRUE(result.second);
  }

  // Test assign behavior of FlatMap::insert_or_assign.
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 7.0);

    auto result = test_map.insert_or_assign(key, value);
    EXPECT_EQ(value, test_map.at(key));
    EXPECT_EQ(result.first, test_map.find(key));
    EXPECT_FALSE(result.second);
  }

  // Assignments should not change size of FlatMap.
  EXPECT_EQ(NUM_ELEMS, test_map.size());
}

TEST(core_flatmap_moveonly, try_emplace)
{
  axom::FlatMap<int, std::unique_ptr<double>> test_map;

  const int NUM_ELEMS = 40;

  // Test behavior when key does not exist.
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    std::unique_ptr<double> value {new double {i + 10.0}};
    auto result = test_map.try_emplace(i, std::move(value));
    EXPECT_EQ(*(test_map[i]), i + 10.0);
    EXPECT_EQ(result.first, test_map.find(i));
    EXPECT_TRUE(result.second);
    // Value should have been moved.
    EXPECT_EQ(value.get(), nullptr);
  }

  // Test behavior when key already exists.
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    std::unique_ptr<double> value {new double {i + 20.0}};
    auto result = test_map.try_emplace(i, std::move(value));
    EXPECT_EQ(*(test_map[i]), i + 10.0);
    EXPECT_EQ(result.first, test_map.find(i));
    EXPECT_FALSE(result.second);
    // Since key already exists, value should NOT be moved.
    EXPECT_NE(value.get(), nullptr);
    EXPECT_EQ(*value, i + 20.0);
  }
}

AXOM_TYPED_TEST(core_flatmap, initializer_list)
{
  using MapType = typename TestFixture::MapType;
  MapType test_map {{this->getKey(0), this->getValue(10.0)},
                    {this->getKey(1), this->getValue(20.0)},
                    {this->getKey(2), this->getValue(30.0)}};

  EXPECT_EQ(3, test_map.size());

  // Check consistency of added values.
  const double expected_str[3] {10.0, 20.0, 30.0};
  for(int i = 0; i < 3; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(expected_str[i]);

    auto iterator = test_map.find(key);
    EXPECT_NE(iterator, test_map.end());
    EXPECT_EQ(iterator->first, key);
    EXPECT_EQ(iterator->second, value);

    // Using operator[] with an already-existing key should return the
    // existing value and not add a value.
    auto indexed_value = test_map[key];
    EXPECT_EQ(indexed_value, value);
    EXPECT_EQ(test_map.size(), 3);
  }
}

AXOM_TYPED_TEST(core_flatmap, index_operator_default)
{
  using MapType = typename TestFixture::MapType;
  MapType test_map;

  const int NUM_ELEMS = 100;

  auto expected_default_value = this->getDefaultValue();
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto default_value = test_map[key];

    EXPECT_EQ(default_value, expected_default_value);

    auto new_value = this->getValue(i * 10.0 + 5.0);
    test_map[key] = new_value;
  }

  EXPECT_EQ(NUM_ELEMS, test_map.size());

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);

    auto iterator = test_map.find(key);
    EXPECT_EQ(iterator->second, value);
  }
}

AXOM_TYPED_TEST(core_flatmap, init_and_clear)
{
  using MapType = typename TestFixture::MapType;
  MapType test_map;

  // Insert enough elements to trigger a resize of the buckets.
  // This allows us to test that a clear() doesn't reset the allocated buckets.
  int NUM_ELEMS_RESIZE = 100;
  EXPECT_GT(NUM_ELEMS_RESIZE, test_map.bucket_count());

  for(int i = 0; i < NUM_ELEMS_RESIZE; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i + 10.0);

    test_map[key] = value;
  }

  EXPECT_EQ(NUM_ELEMS_RESIZE, test_map.size());

  int buckets_before_clear = test_map.bucket_count();

  test_map.clear();

  EXPECT_EQ(test_map.size(), 0);
  EXPECT_EQ(test_map.load_factor(), 0.0);
  EXPECT_EQ(test_map.bucket_count(), buckets_before_clear);
  for(int i = 0; i < NUM_ELEMS_RESIZE; i++)
  {
    auto key = this->getKey(i);
    auto iterator = test_map.find(key);
    EXPECT_EQ(iterator, test_map.end());
  }
}

AXOM_TYPED_TEST(core_flatmap, init_and_move)
{
  using MapType = typename TestFixture::MapType;
  MapType test_map;
  int NUM_ELEMS = 40;

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i + 10.0);

    test_map[key] = value;
  }

  MapType moved_to_map = std::move(test_map);

  EXPECT_EQ(test_map.size(), 0);
  EXPECT_EQ(test_map.load_factor(), 0);
  EXPECT_EQ(moved_to_map.size(), NUM_ELEMS);
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i + 10.0);
    EXPECT_EQ(moved_to_map[key], value);

    auto old_it = test_map.find(key);
    EXPECT_EQ(old_it, test_map.end());
  }
}

TEST(core_flatmap_moveonly, init_and_move_moveonly)
{
  axom::FlatMap<int, std::unique_ptr<double>> test_map;
  int NUM_ELEMS = 40;

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    test_map.emplace(i, new double {i + 10.0});
  }

  axom::FlatMap<int, std::unique_ptr<double>> int_to_dbl_move = std::move(test_map);

  EXPECT_EQ(test_map.size(), 0);
  EXPECT_EQ(test_map.load_factor(), 0);
  EXPECT_EQ(int_to_dbl_move.size(), NUM_ELEMS);
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    EXPECT_EQ(*(int_to_dbl_move[i]), i + 10.0);

    auto old_it = test_map.find(i);
    EXPECT_EQ(old_it, test_map.end());
  }
}

AXOM_TYPED_TEST(core_flatmap, init_and_copy)
{
  using MapType = typename TestFixture::MapType;
  MapType test_map;
  int NUM_ELEMS = 40;

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i + 10.0);

    test_map[key] = value;
  }

  int expected_buckets = test_map.bucket_count();

  MapType int_to_dbl_copy = test_map;

  EXPECT_EQ(test_map.size(), NUM_ELEMS);
  EXPECT_EQ(test_map.bucket_count(), expected_buckets);
  EXPECT_EQ(int_to_dbl_copy.size(), NUM_ELEMS);
  EXPECT_EQ(int_to_dbl_copy.bucket_count(), expected_buckets);
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i + 10.0);

    EXPECT_EQ(test_map[key], value);
    EXPECT_EQ(int_to_dbl_copy[key], value);
  }
}

AXOM_TYPED_TEST(core_flatmap, insert_until_rehash)
{
  using MapType = typename TestFixture::MapType;
  MapType test_map;

  const int INIT_CAPACITY = test_map.bucket_count();
  const double LOAD_FACTOR = test_map.max_load_factor();
  const int SIZE_NO_REHASH = LOAD_FACTOR * INIT_CAPACITY;

  for(int i = 0; i < SIZE_NO_REHASH; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(2. * i + 1);

    test_map.insert({key, value});
  }
  EXPECT_EQ(test_map.bucket_count(), INIT_CAPACITY);
  EXPECT_EQ(test_map.size(), SIZE_NO_REHASH);

  // Next insert should trigger a rehash.
  {
    auto key_rehash = this->getKey(SIZE_NO_REHASH);
    auto value_rehash = this->getValue(2. * SIZE_NO_REHASH + 1);

    test_map.insert({key_rehash, value_rehash});
  }
  EXPECT_GT(test_map.bucket_count(), INIT_CAPACITY);
  EXPECT_EQ(test_map.size(), SIZE_NO_REHASH + 1);

  // Check consistency of values.
  for(int i = 0; i < SIZE_NO_REHASH + 1; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(2. * i + 1);

    auto iterator = test_map.find(key);
    EXPECT_NE(iterator, test_map.end());
    EXPECT_EQ(iterator->first, key);
    EXPECT_EQ(iterator->second, value);
  }
}

AXOM_TYPED_TEST(core_flatmap, insert_then_delete)
{
  using MapType = typename TestFixture::MapType;
  MapType test_map;

  const int INIT_CAPACITY = test_map.bucket_count();
  const double LOAD_FACTOR = test_map.max_load_factor();
  const int NUM_INSERTS = LOAD_FACTOR * INIT_CAPACITY * 4;

  for(int i = 0; i < NUM_INSERTS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(2. * i + 1);

    test_map.insert({key, value});
  }
  EXPECT_EQ(test_map.size(), NUM_INSERTS);
  EXPECT_GE(test_map.bucket_count(), NUM_INSERTS);

  int num_deletions = 0;
  for(int i = 0; i < NUM_INSERTS; i += 3)
  {
    auto key = this->getKey(i);

    auto iterator_to_remove = test_map.find(key);
    auto one_after_elem = iterator_to_remove;
    one_after_elem++;
    // Delete every third entry starting from 0, inclusive.
    // (i.e. keys 0, 3, 6, ...)
    auto deleted_iterator = test_map.erase(iterator_to_remove);

    EXPECT_EQ(deleted_iterator, one_after_elem);
    num_deletions++;
  }

  EXPECT_EQ(test_map.size(), NUM_INSERTS - num_deletions);

  // Check consistency of values.
  for(int i = 0; i < NUM_INSERTS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(2. * i + 1);

    auto iterator = test_map.find(key);
    if(i % 3 == 0)
    {
      EXPECT_EQ(iterator, test_map.end());
      EXPECT_EQ(0, test_map.count(key));
      EXPECT_EQ(false, test_map.contains(key));
    }
    else
    {
      EXPECT_NE(iterator, test_map.end());
      EXPECT_EQ(iterator->first, key);
      EXPECT_EQ(iterator->second, value);
      EXPECT_EQ(1, test_map.count(key));
      EXPECT_EQ(true, test_map.contains(key));
    }
  }
}

AXOM_TYPED_TEST(core_flatmap, iterator_loop)
{
  using MapType = typename TestFixture::MapType;
  MapType test_map;

  const int NUM_ELEMS = 100;

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);

    test_map.insert({key, value});
  }

  std::vector<int> have_iterator(NUM_ELEMS, 0);

  int iter_count = 0;
  // Test constant iteration
  for(typename MapType::const_iterator it = test_map.begin(); it != test_map.end(); ++it)
  {
    auto pair = *it;
    auto iter_key = pair.first;
    auto iter_value = pair.second;

    // Get the original integer value of the key.
    int iter_key_int;
    flatmap_get_value(iter_key, iter_key_int);

    // Check that the key value is in range.
    EXPECT_GE(iter_key_int, 0);
    EXPECT_LT(iter_key_int, NUM_ELEMS);

    // Count the key that we got.
    have_iterator[iter_key_int]++;

    // Check that the value is what we expect for the given key..
    auto expected_value = this->getValue(iter_key_int * 10.0 + 5.0);
    EXPECT_EQ(iter_value, expected_value);

    // Count the number of iterations.
    iter_count++;
  }
  EXPECT_EQ(iter_count, NUM_ELEMS);

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    // We should have iterated through every index exactly once.
    EXPECT_EQ(have_iterator[i], 1);
  }
}

AXOM_TYPED_TEST(core_flatmap, iterator_loop_write)
{
  using MapType = typename TestFixture::MapType;
  MapType test_map;

  const int NUM_ELEMS = 100;

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);

    test_map.insert({key, value});
  }

  // Test mutable iteration
  for(typename MapType::iterator it = test_map.begin(); it != test_map.end(); ++it)
  {
    auto pair = *it;
    auto iter_key = pair.first;

    // Get the original integer value of the key.
    int iter_key_int;
    flatmap_get_value(iter_key, iter_key_int);

    // Modify the stored value.
    it->second = this->getValue(iter_key_int * 10.0 + 7.0);
  }

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    // All values should be set to the new value.
    const auto key = this->getKey(i);
    const auto expected_value = this->getValue(i * 10.0 + 7.0);
    EXPECT_EQ(test_map[key], expected_value);
  }
}

AXOM_TYPED_TEST(core_flatmap, range_for_loop)
{
  using MapType = typename TestFixture::MapType;
  MapType test_map;

  const int NUM_ELEMS = 100;

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);

    test_map.insert({key, value});
  }

  std::vector<int> have_key(NUM_ELEMS, 0);

  int iter_count = 0;
  // Test constant iteration
  for(const auto& pair : test_map)
  {
    auto iter_key = pair.first;
    auto iter_value = pair.second;

    // Get the original integer value of the key.
    int iter_key_int;
    flatmap_get_value(iter_key, iter_key_int);

    // Check that the key value is in range.
    EXPECT_GE(iter_key_int, 0);
    EXPECT_LT(iter_key_int, NUM_ELEMS);

    // Count the key that we got.
    have_key[iter_key_int]++;

    // Check that the value is what we expect for the given key..
    auto expected_value = this->getValue(iter_key_int * 10.0 + 5.0);
    EXPECT_EQ(iter_value, expected_value);

    // Count the number of iterations.
    iter_count++;
  }
  EXPECT_EQ(iter_count, NUM_ELEMS);

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    // We should have iterated through every index exactly once.
    EXPECT_EQ(have_key[i], 1);
  }
}

AXOM_TYPED_TEST(core_flatmap, range_for_loop_write)
{
  using MapType = typename TestFixture::MapType;
  MapType test_map;

  const int NUM_ELEMS = 100;

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);

    test_map.insert({key, value});
  }

  // Test mutable iteration
  for(auto& pair : test_map)
  {
    auto iter_key = pair.first;

    // Get the original integer value of the key.
    int iter_key_int;
    flatmap_get_value(iter_key, iter_key_int);

    // Modify the stored value.
    pair.second = this->getValue(iter_key_int * 10.0 + 7.0);
  }

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    // All values should be set to the new value.
    const auto key = this->getKey(i);
    const auto expected_value = this->getValue(i * 10.0 + 7.0);
    EXPECT_EQ(test_map[key], expected_value);
  }
}
#if defined(AXOM_USE_HIP) || defined(AXOM_USE_CUDA)

AXOM_TYPED_TEST(core_flatmap, copy_host_device)
{
  using MapType = typename TestFixture::MapType;

  using HostExec = typename TestFixture::HostExec;
  using DeviceExec = typename TestFixture::DeviceExec;

  // CUDA failure - Skip tests if key or value is of type std::string
  #if defined(AXOM_USE_CUDA)
  if constexpr(std::is_same<typename TestFixture::KeyType, std::string>::value ||
               std::is_same<typename TestFixture::ValueType, std::string>::value)
  {
    return;
  }
  #endif

  const int host_alloc_id = axom::execution_space<HostExec>::allocatorID();
  const int device_alloc_id = axom::execution_space<DeviceExec>::allocatorID();

  MapType test_map(axom::Allocator {host_alloc_id});

  const int NUM_ELEMS = 100;
  const int EMPTY_CHECKS = 100;

  // Create map on the host.
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);

    test_map.insert({key, value});
  }

  // Copy from the host to the device.
  MapType test_map_device = MapType(test_map, axom::Allocator {device_alloc_id});

  // Check that maps are equivalent with a different allocator ID.
  EXPECT_EQ(test_map_device.getAllocator().getID(), device_alloc_id);
  EXPECT_EQ(test_map_device.bucket_count(), test_map.bucket_count());
  EXPECT_EQ(test_map_device.size(), test_map.size());

  // Copy back to the host.
  MapType test_map_host = MapType(test_map_device, axom::Allocator {host_alloc_id});

  // Check that maps are equivalent with the same allocator ID.
  EXPECT_EQ(test_map_host.getAllocator().getID(), host_alloc_id);
  EXPECT_EQ(test_map_host.bucket_count(), test_map.bucket_count());
  EXPECT_EQ(test_map_host.size(), test_map.size());

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);

    EXPECT_EQ(test_map_host.contains(key), true);
    EXPECT_EQ(test_map_host.count(key), 1);
    EXPECT_EQ(test_map_host.find(key)->first, test_map.find(key)->first);
    EXPECT_EQ(test_map_host.find(key)->second, test_map.find(key)->second);
  }

  for(int i = NUM_ELEMS; i < NUM_ELEMS + EMPTY_CHECKS; i++)
  {
    auto key = this->getKey(i);

    EXPECT_EQ(test_map_host.contains(key), false);
    EXPECT_EQ(test_map_host.count(key), 0);
    EXPECT_EQ(test_map_host.find(key), test_map_host.end());
  }
}

#endif

AXOM_TYPED_TEST(core_flatmap, empty_view)
{
  using MapType = typename TestFixture::MapType;
  MapType test_map;

  const int EMPTY_CHECKS = 100;

  // Create a read-only view from the FlatMap
  auto view = test_map.view();

  EXPECT_EQ(view.empty(), test_map.empty());
  EXPECT_EQ(view.size(), test_map.size());
  EXPECT_EQ(view.bucket_count(), test_map.bucket_count());

  for(int i = 0; i < EMPTY_CHECKS; i++)
  {
    auto key = this->getKey(i);

    EXPECT_EQ(view.contains(key), false);
    EXPECT_EQ(view.count(key), 0);
    EXPECT_EQ(view.find(key), view.end());
  }
}

AXOM_TYPED_TEST(core_flatmap, view_iterators)
{
  using MapType = typename TestFixture::MapType;
  using MapViewType = typename MapType::View;
  using MapIterType = typename MapViewType::iterator;

  MapIterType empty_iterator;
  MapIterType empty_iterator_2;

  // Default-constructed iterators should all be equivalent.
  EXPECT_EQ(empty_iterator, empty_iterator_2);

  MapType test_map;

  const int NUM_ELEMS = 100;

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);

    test_map.emplace(key, value);
  }

  // Get a view over the test map.
  auto test_map_view = test_map.view();

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);

    auto it = test_map_view.find(key);

    // Default-constructed iterators should not be equivalent to an actual map's
    // iterators.
    EXPECT_NE(it, empty_iterator);
  }

  // Get a second map view over the same test map.
  auto test_map_view_2 = test_map.view();

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);

    auto it = test_map_view.find(key);
    auto it2 = test_map_view_2.find(key);

    // Iterators between two view instances should compare equal if the two views
    // point to the same underlying FlatMap.
    EXPECT_EQ(it, it2);
  }
}

AXOM_TYPED_TEST(core_flatmap, init_view)
{
  using MapType = typename TestFixture::MapType;
  MapType test_map;

  const int NUM_ELEMS = 100;
  const int EMPTY_CHECKS = 100;

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);

    test_map.insert({key, value});
  }

  // Create a read-only view from the FlatMap
  auto view = test_map.view();

  EXPECT_EQ(view.empty(), test_map.empty());
  EXPECT_EQ(view.size(), test_map.size());
  EXPECT_EQ(view.bucket_count(), test_map.bucket_count());

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);

    EXPECT_EQ(view.contains(key), true);
    EXPECT_EQ(view.count(key), 1);
    EXPECT_EQ(view.find(key)->first, test_map.find(key)->first);
    EXPECT_EQ(view.find(key)->second, test_map.find(key)->second);
  }

  for(int i = NUM_ELEMS; i < NUM_ELEMS + EMPTY_CHECKS; i++)
  {
    auto key = this->getKey(i);

    EXPECT_EQ(view.contains(key), false);
    EXPECT_EQ(view.count(key), 0);
    EXPECT_EQ(view.find(key), view.end());
  }
}
