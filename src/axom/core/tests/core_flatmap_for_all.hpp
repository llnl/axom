// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/FlatMap.hpp"
#include "axom/core/FlatMapView.hpp"

// gtest includes
#include "gtest/gtest.h"

#include <random>

template <typename FlatMapType, typename ExecSpace, axom::MemorySpace SPACE = axom::MemorySpace::Dynamic>
struct FlatMapTestParams
{
  using ViewExecSpace = ExecSpace;
  using MapType = FlatMapType;
  static constexpr axom::MemorySpace KernelSpace = SPACE;
};

template <typename ExecParams>
class core_flatmap_for_all : public ::testing::Test
{
public:
  using MapType = typename ExecParams::MapType;
  using MapViewType = typename MapType::View;
  using MapViewConstType = typename MapType::ConstView;
  using KeyType = typename MapType::key_type;
  using ValueType = typename MapType::mapped_type;
  using ExecSpace = typename ExecParams::ViewExecSpace;

  template <typename T>
  KeyType getKey(T input)
  {
    return (KeyType)input;
  }

  template <typename T>
  ValueType getValue(T input)
  {
    return (ValueType)input;
  }

  ValueType getDefaultValue() { return ValueType(); }

  static int getKernelAllocatorID()
  {
    return axom::detail::getAllocatorID<ExecParams::KernelSpace>();
  }
  static int getHostAllocatorID()
  {
#ifdef AXOM_USE_UMPIRE
    static constexpr axom::MemorySpace HostSpace = axom::MemorySpace::Host;
#else
    static constexpr axom::MemorySpace HostSpace = axom::MemorySpace::Dynamic;
#endif
    return axom::detail::getAllocatorID<HostSpace>();
  }
};

/**
 * Hash that returns the same value for all elements.
 * Even with worst-case hash collision behavior, the FlatMap should
 * nevertheless be correctly constructible and queryable.
 */
template <typename KeyType>
struct ConstantHash
{
  using argument_type = KeyType;
  using result_type = axom::IndexType;

  AXOM_HOST_DEVICE axom::IndexType operator()(KeyType) const { return 0; }
};

using ViewTypes = ::testing::Types<
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP)
  FlatMapTestParams<axom::FlatMap<int, double>, axom::OMP_EXEC>,
  FlatMapTestParams<axom::FlatMap<int, double, ConstantHash<int>>, axom::OMP_EXEC>,
#endif
// TODO - Batch insertion for CUDA pinned memory fails (less than NUM_ELEMS insertions or deadlock)
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
  FlatMapTestParams<axom::FlatMap<int, double>, axom::CUDA_EXEC<256>, axom::MemorySpace::Device>,
  FlatMapTestParams<axom::FlatMap<int, double>, axom::CUDA_EXEC<256>, axom::MemorySpace::Unified>,
  // FlatMapTestParams<axom::FlatMap<int, double>, axom::CUDA_EXEC<256>, axom::MemorySpace::Pinned>,
  FlatMapTestParams<axom::FlatMap<int, double, ConstantHash<int>>, axom::CUDA_EXEC<256>, axom::MemorySpace::Device>,
  FlatMapTestParams<axom::FlatMap<int, double, ConstantHash<int>>, axom::CUDA_EXEC<256>, axom::MemorySpace::Unified>,
// FlatMapTestParams<axom::FlatMap<int, double, ConstantHash<int>>, axom::CUDA_EXEC<256>, axom::MemorySpace::Pinned>,
#endif
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(AXOM_USE_UMPIRE)
  FlatMapTestParams<axom::FlatMap<int, double>, axom::HIP_EXEC<256>, axom::MemorySpace::Device>,
  FlatMapTestParams<axom::FlatMap<int, double>, axom::HIP_EXEC<256>, axom::MemorySpace::Unified>,
  FlatMapTestParams<axom::FlatMap<int, double>, axom::HIP_EXEC<256>, axom::MemorySpace::Pinned>,
  FlatMapTestParams<axom::FlatMap<int, double, ConstantHash<int>>, axom::HIP_EXEC<256>, axom::MemorySpace::Device>,
  FlatMapTestParams<axom::FlatMap<int, double, ConstantHash<int>>, axom::HIP_EXEC<256>, axom::MemorySpace::Unified>,
  FlatMapTestParams<axom::FlatMap<int, double, ConstantHash<int>>, axom::HIP_EXEC<256>, axom::MemorySpace::Pinned>,
#endif
#if defined(AXOM_USE_UMPIRE)
  FlatMapTestParams<axom::FlatMap<int, double>, axom::SEQ_EXEC, axom::MemorySpace::Host>,
  FlatMapTestParams<axom::FlatMap<int, double, ConstantHash<int>>, axom::SEQ_EXEC, axom::MemorySpace::Host>,
#endif
  FlatMapTestParams<axom::FlatMap<int, double>, axom::SEQ_EXEC>,
  FlatMapTestParams<axom::FlatMap<int, double, ConstantHash<int>>, axom::SEQ_EXEC>>;

TYPED_TEST_SUITE(core_flatmap_for_all, ViewTypes);

AXOM_TYPED_TEST(core_flatmap_for_all, insert_and_find)
{
  using MapType = typename TestFixture::MapType;
  using MapViewConstType = typename TestFixture::MapViewConstType;
  using ExecSpace = typename TestFixture::ExecSpace;

  MapType test_map;

  const int NUM_ELEMS = 100;
  const int EXTRA_THREADS = 100;

  // First do insertions of elements.
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);

    test_map.insert({key, value});
  }

  MapType test_map_gpu(test_map, axom::Allocator {this->getKernelAllocatorID()});
  MapViewConstType test_map_view(test_map_gpu);

  const int TOTAL_NUM_THREADS = NUM_ELEMS + EXTRA_THREADS;
  axom::Array<int> valid_vec(TOTAL_NUM_THREADS, TOTAL_NUM_THREADS, this->getKernelAllocatorID());
  axom::Array<int> keys_vec(NUM_ELEMS, NUM_ELEMS, this->getKernelAllocatorID());
  axom::Array<double> values_vec(NUM_ELEMS, NUM_ELEMS, this->getKernelAllocatorID());
  axom::Array<double> values_vec_bracket(TOTAL_NUM_THREADS,
                                         TOTAL_NUM_THREADS,
                                         this->getKernelAllocatorID());
  const auto valid_out = valid_vec.view();
  const auto keys_out = keys_vec.view();
  const auto values_out = values_vec.view();
  const auto values_out_bracket = values_vec_bracket.view();

  // Read values out in a captured lambda.
  axom::for_all<ExecSpace>(
    NUM_ELEMS + EXTRA_THREADS,
    AXOM_LAMBDA(axom::IndexType idx) {
      auto it = test_map_view.find(idx);
      if(it != test_map_view.end())
      {
        keys_out[idx] = it->first;
        values_out[idx] = it->second;
        valid_out[idx] = true;
      }
      else
      {
        valid_out[idx] = false;
      }
      values_out_bracket[idx] = test_map_view[idx];
    });

  axom::Array<int> valid_host(valid_vec, this->getHostAllocatorID());
  axom::Array<int> keys_host(keys_vec, this->getHostAllocatorID());
  axom::Array<double> values_host(values_vec, this->getHostAllocatorID());
  axom::Array<double> values_host_bracket(values_vec_bracket, this->getHostAllocatorID());

  // Check contents on the host
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    EXPECT_EQ(valid_host[i], true);
    EXPECT_EQ(keys_host[i], this->getKey(i));
    EXPECT_EQ(values_host[i], this->getValue(i * 10.0 + 5.0));
    EXPECT_EQ(values_host_bracket[i], this->getValue(i * 10.0 + 5.0));
  }
  for(int i = NUM_ELEMS; i < NUM_ELEMS + EXTRA_THREADS; i++)
  {
    EXPECT_EQ(valid_host[i], false);
    EXPECT_EQ(values_host_bracket[i], this->getValue(0));
  }
}

AXOM_TYPED_TEST(core_flatmap_for_all, insert_and_modify)
{
  using MapType = typename TestFixture::MapType;
  using MapViewType = typename TestFixture::MapViewType;
  using ExecSpace = typename TestFixture::ExecSpace;

  MapType test_map;

  const int NUM_ELEMS = 100;
  const int EXTRA_THREADS = 100;

  // First do insertions of elements.
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);

    test_map.insert({key, value});
  }

  MapType test_map_gpu(test_map, axom::Allocator {this->getKernelAllocatorID()});
  MapViewType test_map_view(test_map_gpu);

  // Write new values into the flat map, where existing keys are.
  // This should work from a map view because we are not inserting
  // existing keys, which would potentially trigger rehashes.
  axom::for_all<ExecSpace>(
    NUM_ELEMS + EXTRA_THREADS,
    AXOM_LAMBDA(axom::IndexType idx) {
      auto it = test_map_view.find(idx);
      if(it != test_map_view.end())
      {
        it->second = idx * 11.0 + 7.0;
      }
    });

  test_map = MapType(test_map_gpu, axom::Allocator {this->getHostAllocatorID()});

  // Check contents of the original map on the host
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    EXPECT_EQ(test_map.count(i), true);
    EXPECT_EQ(test_map.find(i)->first, this->getKey(i));
    // Ensure that this k-v pair has an updated value
    EXPECT_EQ(test_map.find(i)->second, this->getValue(i * 11.0 + 7.0));
    EXPECT_NE(test_map.find(i)->second, this->getValue(i * 10.0 + 5.0));
  }
  for(int i = NUM_ELEMS; i < NUM_ELEMS + EXTRA_THREADS; i++)
  {
    EXPECT_EQ(test_map.count(i), false);
  }
}

AXOM_TYPED_TEST(core_flatmap_for_all, insert_batched)
{
  using MapType = typename TestFixture::MapType;
  using ExecSpace = typename TestFixture::ExecSpace;

  const int NUM_ELEMS = 100;

  axom::Array<int> keys_vec(NUM_ELEMS);
  axom::Array<double> values_vec(NUM_ELEMS);
  // Create batch of array elements
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);

    keys_vec[i] = key;
    values_vec[i] = value;
  }

  // Copy keys and values to GPU space.
  axom::Array<int> keys_gpu(keys_vec, this->getKernelAllocatorID());
  axom::Array<double> values_gpu(values_vec, this->getKernelAllocatorID());

  // Construct a flat map with the key-value pairs.
  MapType test_map_gpu =
    MapType::template create<ExecSpace>(keys_gpu,
                                        values_gpu,
                                        axom::Allocator {this->getKernelAllocatorID()});

  // Copy back flat map to host for testing.
  MapType test_map(test_map_gpu, axom::Allocator {this->getHostAllocatorID()});

  // Check contents on the host
  EXPECT_EQ(NUM_ELEMS, test_map.size());

  // Check that every element we inserted is in the map
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto expected_key = this->getKey(i);
    auto expected_val = this->getValue(i * 10.0 + 5.0);
    EXPECT_EQ(1, test_map.count(expected_key));
    EXPECT_EQ(expected_val, test_map.at(expected_key));
  }
}

AXOM_TYPED_TEST(core_flatmap_for_all, insert_batched_with_existing)
{
  using MapType = typename TestFixture::MapType;
  using ExecSpace = typename TestFixture::ExecSpace;

  const int NUM_ELEMS_INIT = 100;
  const int NUM_ELEMS_INSERT = 100;
  const int NUM_ELEMS = NUM_ELEMS_INIT + NUM_ELEMS_INSERT;

  axom::Array<int> keys_vec(NUM_ELEMS_INIT);
  axom::Array<double> values_vec(NUM_ELEMS_INIT);
  // Create batch of array elements
  for(int i = 0; i < NUM_ELEMS_INIT; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);

    keys_vec[i] = key;
    values_vec[i] = value;
  }

  // Copy keys and values to GPU space.
  axom::Array<int> keys_gpu(keys_vec, this->getKernelAllocatorID());
  axom::Array<double> values_gpu(values_vec, this->getKernelAllocatorID());

  // Construct a flat map with the key-value pairs.
  MapType test_map_gpu =
    MapType::template create<ExecSpace>(keys_gpu,
                                        values_gpu,
                                        axom::Allocator {this->getKernelAllocatorID()});

  // Create batch of pairs.
  axom::Array<std::pair<int, double>> kv_insert_vec(NUM_ELEMS_INSERT);
  for(int i = 0; i < NUM_ELEMS_INSERT; i++)
  {
    int offset_i = i + NUM_ELEMS_INIT;
    auto key = this->getKey(offset_i);
    auto value = this->getValue(offset_i * 10.0 + 5.0);

    kv_insert_vec[i] = {key, value};
  }

  // Copy pairs to GPU space.
  axom::Array<std::pair<int, double>> kv_insert_gpu(kv_insert_vec, this->getKernelAllocatorID());

  // Insert pairs into existing flatmap.
  test_map_gpu.template insert<ExecSpace>(kv_insert_gpu.data(),
                                          kv_insert_gpu.data() + NUM_ELEMS_INSERT);

  // Copy back flat map to host for testing.
  MapType test_map(test_map_gpu, axom::Allocator {this->getHostAllocatorID()});

  // Check contents on the host
  EXPECT_EQ(NUM_ELEMS, test_map.size());

  // Check that every element we inserted is in the map
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto expected_key = this->getKey(i);
    auto expected_val = this->getValue(i * 10.0 + 5.0);
    EXPECT_EQ(1, test_map.count(expected_key));
    EXPECT_EQ(expected_val, test_map.at(expected_key));
  }
}

AXOM_TYPED_TEST(core_flatmap_for_all, insert_batched_with_dups)
{
  using MapType = typename TestFixture::MapType;
  using ExecSpace = typename TestFixture::ExecSpace;

  const int NUM_ELEMS = 100;

  axom::Array<int> keys_vec(NUM_ELEMS * 2);
  axom::Array<double> values_vec(NUM_ELEMS * 2);
  // Create batch of array elements
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);

    keys_vec[i] = key;
    values_vec[i] = value;
  }

  // Add some duplicate key values
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 7.0);

    keys_vec[i + NUM_ELEMS] = key;
    values_vec[i + NUM_ELEMS] = value;
  }

  // Copy keys and values to GPU space.
  axom::Array<int> keys_gpu(keys_vec, this->getKernelAllocatorID());
  axom::Array<double> values_gpu(values_vec, this->getKernelAllocatorID());

  // Construct a flat map with the key-value pairs.
  MapType test_map_gpu =
    MapType::template create<ExecSpace>(keys_gpu,
                                        values_gpu,
                                        axom::Allocator {this->getKernelAllocatorID()});

  // Copy back flat map to host for testing.
  MapType test_map(test_map_gpu, axom::Allocator {this->getHostAllocatorID()});

  // Check contents on the host. Only one of the duplicate keys should remain.
  EXPECT_EQ(NUM_ELEMS, test_map.size());

  // Check that every element we inserted is in the map
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto expected_key = this->getKey(i);
    auto expected_val1 = this->getValue(i * 10.0 + 5.0);
    auto expected_val2 = this->getValue(i * 10.0 + 7.0);
    EXPECT_EQ(1, test_map.count(expected_key));
    // Second key-value pair in batch-order should overwrite first pair with
    // same key.
    EXPECT_EQ(expected_val2, test_map.at(expected_key));
    EXPECT_NE(expected_val1, test_map.at(expected_key));
  }

  // Check that we only have one instance of every key in the map
  axom::Array<int> dedup_keys(NUM_ELEMS);
  for(auto &pair : test_map)
  {
    // Check that we haven't seen another K-V pair with the same key.
    EXPECT_EQ(dedup_keys[pair.first], 0);
    dedup_keys[pair.first]++;

    // Check that we got the second KV pair, not the first.
    auto expected_val1 = this->getValue(pair.first * 10.0 + 5.0);
    auto expected_val2 = this->getValue(pair.first * 10.0 + 7.0);
    EXPECT_EQ(expected_val2, pair.second);
    EXPECT_NE(expected_val1, pair.second);
  }
}

AXOM_TYPED_TEST(core_flatmap_for_all, insert_multiple_batch_with_dups)
{
  using MapType = typename TestFixture::MapType;
  using ExecSpace = typename TestFixture::ExecSpace;

  const int NUM_ELEMS = 100;

  axom::Array<int> keys_vec(NUM_ELEMS);
  axom::Array<double> values_vec(NUM_ELEMS);
  // Create batch of array elements
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);

    keys_vec[i] = key;
    values_vec[i] = value;
  }

  // Copy keys and values to GPU space.
  axom::Array<int> keys_gpu(keys_vec, this->getKernelAllocatorID());
  axom::Array<double> values_gpu(values_vec, this->getKernelAllocatorID());

  // Construct a flat map with the key-value pairs.
  MapType test_map_gpu =
    MapType::template create<ExecSpace>(keys_gpu,
                                        values_gpu,
                                        axom::Allocator {this->getKernelAllocatorID()});

  axom::Array<std::pair<int, double>> second_batch_pairs(NUM_ELEMS);
  // Add some duplicate key values through the batched interface.
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 7.0);

    second_batch_pairs[i] = {key, value};
  }
  // Copy pairs to GPU space.
  axom::Array<std::pair<int, double>> second_batch_gpu(second_batch_pairs,
                                                       this->getKernelAllocatorID());

  test_map_gpu.template insert<ExecSpace>(second_batch_gpu.data(),
                                          second_batch_gpu.data() + NUM_ELEMS);

  // Copy back flat map to host for testing.
  MapType test_map(test_map_gpu, axom::Allocator {this->getHostAllocatorID()});

  // Check contents on the host. Only one of the duplicate keys should remain.
  EXPECT_EQ(NUM_ELEMS, test_map.size());

  // Check that every element we inserted is in the map
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto expected_key = this->getKey(i);
    auto expected_val1 = this->getValue(i * 10.0 + 5.0);
    auto expected_val2 = this->getValue(i * 10.0 + 7.0);
    EXPECT_EQ(1, test_map.count(expected_key));
    // Second key-value pair in batch-order should overwrite first pair with
    // same key.
    EXPECT_EQ(expected_val2, test_map.at(expected_key));
    EXPECT_NE(expected_val1, test_map.at(expected_key));
  }

  // Check that we only have one instance of every key in the map
  axom::Array<int> dedup_keys(NUM_ELEMS);
  for(auto &pair : test_map)
  {
    // Check that we haven't seen another K-V pair with the same key.
    EXPECT_EQ(dedup_keys[pair.first], 0);
    dedup_keys[pair.first]++;

    // Check that we got the second KV pair, not the first.
    auto expected_val1 = this->getValue(pair.first * 10.0 + 5.0);
    auto expected_val2 = this->getValue(pair.first * 10.0 + 7.0);
    EXPECT_EQ(expected_val2, pair.second);
    EXPECT_NE(expected_val1, pair.second);
  }
}

/**
 * Rigorous test of FlatMap probing insert behavior.
 * High level steps:
 *  - Insert elements, along common probing path
 *  - Delete a few of the elements on the same probing path
 *  - Perform a batched insert, with duplicates.
 */
AXOM_TYPED_TEST(core_flatmap_for_all, insert_batch_with_gaps_and_dups)
{
  using MapType = typename TestFixture::MapType;
  using ExecSpace = typename TestFixture::ExecSpace;

  const int NUM_ELEMS = 200;

  // Seed triplet observed to trigger an insertion regression in an OpenMP configuration
  const std::mt19937::result_type seed_insert = 1087231065u;
  const std::mt19937::result_type seed_batch_1 = 1693880942u;
  const std::mt19937::result_type seed_batch_2 = 3511204532u;

  // Repeat the test as necessary to catch scheduling dependent regression
  int num_trials = 20;
  if(const char *env_iters = std::getenv("AXOM_FLATMAP_TEST_STRESS_ITERS"))
  {
    num_trials = std::max(1, std::atoi(env_iters));
  }

  for(int trial = 0; trial < num_trials; ++trial)
  {
    // Allocate enough space to ensure rehashes don't eliminate probing sequence gaps.
    MapType test_map(NUM_ELEMS * 4);

    // Shuffle inserted elements.
    std::vector<int> shuffled_indexes(NUM_ELEMS);
    std::iota(shuffled_indexes.begin(), shuffled_indexes.end(), 0);
    std::shuffle(shuffled_indexes.begin(), shuffled_indexes.end(), std::mt19937 {seed_insert});

    // Insert initial elements, in shuffled order.
    for(int i = 0; i < NUM_ELEMS; i++)
    {
      int shuf_i = shuffled_indexes[i];
      auto key = this->getKey(shuf_i);
      auto value = this->getValue(shuf_i * 10.0 + 5.0);

      test_map.emplace(key, value);
    }

    // Delete every third element.
    // This is intended to creates gaps in probing sequences.
    int num_erases = 0;
    for(int i = 0; i < NUM_ELEMS / 3; i++)
    {
      const int index = i * 3 + 2;
      const auto key = this->getKey(index);
      const auto it = test_map.find(key);
      EXPECT_TRUE(it != test_map.end());

      test_map.erase(it);
      num_erases++;
    }

    EXPECT_EQ(test_map.size(), NUM_ELEMS - num_erases);

    MapType test_map_gpu(test_map, axom::Allocator {this->getKernelAllocatorID()});

    axom::Array<std::pair<int, double>> second_batch_pairs(NUM_ELEMS * 2);

    // Add some duplicate key values through the batched interface.
    std::shuffle(shuffled_indexes.begin(), shuffled_indexes.end(), std::mt19937 {seed_batch_1});
    for(int i = 0; i < NUM_ELEMS; i++)
    {
      const int shuf_i = shuffled_indexes[i];
      const auto key = this->getKey(shuf_i);
      const auto value = this->getValue(shuf_i * 10.0 + 6.0);

      second_batch_pairs[i] = {key, value};
    }

    // Add a second set of duplicates.
    // Since these are at the end of the sequence, they should be the ones written into the map.
    std::shuffle(shuffled_indexes.begin(), shuffled_indexes.end(), std::mt19937 {seed_batch_2});
    for(int i = 0; i < NUM_ELEMS; i++)
    {
      const int shuf_i = shuffled_indexes[i];
      const auto key = this->getKey(shuf_i);
      const auto value = this->getValue(shuf_i * 10.0 + 7.0);

      second_batch_pairs[i + NUM_ELEMS] = {key, value};
    }
    // Copy pairs to GPU space.
    axom::Array<std::pair<int, double>> second_batch_gpu(second_batch_pairs,
                                                         this->getKernelAllocatorID());

    // Perform batched insert.
    test_map_gpu.template insert<ExecSpace>(second_batch_gpu.data(),
                                            second_batch_gpu.data() + NUM_ELEMS * 2);

    // Copy back flat map to host for testing.
    test_map = MapType(test_map_gpu, axom::Allocator {this->getHostAllocatorID()});

    // Check contents on the host. Only one of the duplicate keys should remain.
    EXPECT_EQ(NUM_ELEMS, test_map.size());

    // Check that every element we inserted is in the map
    for(int i = 0; i < NUM_ELEMS; i++)
    {
      auto expected_key = this->getKey(i);
      auto expected_val1 = this->getValue(i * 10.0 + 5.0);
      auto expected_val2 = this->getValue(i * 10.0 + 7.0);
      EXPECT_EQ(1, test_map.count(expected_key));
      // Second key-value pair in batch-order should overwrite first pair with same key
      EXPECT_EQ(expected_val2, test_map.at(expected_key));
      EXPECT_NE(expected_val1, test_map.at(expected_key));
    }

    // Check that we only have one instance of every key in the map
    axom::Array<int> dedup_keys(NUM_ELEMS);
    for(auto &pair : test_map)
    {
      // Check that we haven't seen another K-V pair with the same key.
      EXPECT_EQ(dedup_keys[pair.first], 0);
      dedup_keys[pair.first]++;

      // Check that we got the second KV pair, not the first.
      auto expected_val1 = this->getValue(pair.first * 10.0 + 5.0);
      auto expected_val2 = this->getValue(pair.first * 10.0 + 7.0);
      EXPECT_EQ(expected_val2, pair.second);
      EXPECT_NE(expected_val1, pair.second);
    }
  }
}
