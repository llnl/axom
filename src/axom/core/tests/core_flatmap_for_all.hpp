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

template <typename FlatMapType>
class core_flatmap_forall : public ::testing::Test
{
public:
  using MapType = FlatMapType;
  using MapViewType = typename FlatMapType::View;
  using KeyType = typename FlatMapType::key_type;
  using ValueType = typename FlatMapType::mapped_type;
  using ExecSpace = axom::SEQ_EXEC;

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
};

using ViewTypes = ::testing::Types<axom::FlatMap<int, double>>;

TYPED_TEST_SUITE(core_flatmap_forall, ViewTypes);

AXOM_TYPED_TEST(core_flatmap_forall, insert_and_find)
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

  MapViewType test_map_view(test_map);

  axom::Array<int> valid_vec(NUM_ELEMS + EXTRA_THREADS, 0);
  axom::Array<int> keys_vec(NUM_ELEMS);
  axom::Array<double> values_vec(NUM_ELEMS);
  const auto valid_out = valid_vec.view();
  const auto keys_out = keys_vec.view();
  const auto values_out = values_vec.view();

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
    });

  // Check contents on the host
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    EXPECT_EQ(valid_out[i], true);
    EXPECT_EQ(keys_out[i], this->getKey(i));
    EXPECT_EQ(values_out[i], this->getValue(i * 10.0 + 5.0));
  }
  for(int i = NUM_ELEMS; i < NUM_ELEMS + EXTRA_THREADS; i++)
  {
    EXPECT_EQ(valid_out[i], false);
  }
}
