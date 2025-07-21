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

template <typename FlatMapType, typename ExecSpace, axom::MemorySpace SPACE = axom::MemorySpace::Dynamic>
struct FlatMapTestParams
{
  using ViewExecSpace = ExecSpace;
  using MapType = FlatMapType;
  static constexpr axom::MemorySpace KernelSpace = SPACE;
};

template <typename ExecParams>
class core_flatmap_forall : public ::testing::Test
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

using ViewTypes = ::testing::Types<
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP)
  FlatMapTestParams<axom::FlatMap<int, double>, axom::OMP_EXEC>,
#endif
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
  FlatMapTestParams<axom::FlatMap<int, double>, axom::CUDA_EXEC<256>, axom::MemorySpace::Device>,
  FlatMapTestParams<axom::FlatMap<int, double>, axom::CUDA_EXEC<256>, axom::MemorySpace::Unified>,
  FlatMapTestParams<axom::FlatMap<int, double>, axom::CUDA_EXEC<256>, axom::MemorySpace::Pinned>,
#endif
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(AXOM_USE_UMPIRE)
  FlatMapTestParams<axom::FlatMap<int, double>, axom::HIP_EXEC<256>, axom::MemorySpace::Device>,
  FlatMapTestParams<axom::FlatMap<int, double>, axom::HIP_EXEC<256>, axom::MemorySpace::Unified>,
  FlatMapTestParams<axom::FlatMap<int, double>, axom::HIP_EXEC<256>, axom::MemorySpace::Pinned>,
#endif
#if defined(AXOM_USE_UMPIRE)
  FlatMapTestParams<axom::FlatMap<int, double>, axom::SEQ_EXEC, axom::MemorySpace::Host>,
#endif
  FlatMapTestParams<axom::FlatMap<int, double>, axom::SEQ_EXEC>>;

TYPED_TEST_SUITE(core_flatmap_forall, ViewTypes);

AXOM_TYPED_TEST(core_flatmap_forall, insert_and_find)
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

AXOM_TYPED_TEST(core_flatmap_forall, insert_and_modify)
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

AXOM_TYPED_TEST(core_flatmap_forall, insert_batched)
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

AXOM_TYPED_TEST(core_flatmap_forall, insert_batched_with_dups)
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
  axom::Array<std::pair<int, double>> kv_out(NUM_ELEMS);
  int index = 0;
  for(auto &pair : test_map)
  {
    EXPECT_LT(index, NUM_ELEMS);
    kv_out[index++] = {pair.first, pair.second};
  }

  std::sort(kv_out.begin(),
            kv_out.end(),
            [](const std::pair<int, double> &first, const std::pair<int, double> &second) -> bool {
              return first.first < second.first;
            });

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto expected_key = this->getKey(i);
    auto expected_val1 = this->getValue(i * 10.0 + 5.0);
    auto expected_val2 = this->getValue(i * 10.0 + 7.0);
    EXPECT_EQ(kv_out[i].first, expected_key);
    EXPECT_EQ(expected_val2, test_map.at(expected_key));
    EXPECT_NE(expected_val1, test_map.at(expected_key));
  }
}

template <typename KeyType>
struct ConstantHash
{
  using argument_type = KeyType;
  using result_type = axom::IndexType;

  AXOM_HOST_DEVICE axom::IndexType operator()(KeyType) const { return 0; }
};

AXOM_TYPED_TEST(core_flatmap_forall, insert_batched_constant_hash)
{
  using ExecSpace = typename TestFixture::ExecSpace;
  using KeyType = typename TestFixture::KeyType;
  using ValueType = typename TestFixture::ValueType;

  using MapType = axom::FlatMap<KeyType, ValueType, ConstantHash<KeyType>>;

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
