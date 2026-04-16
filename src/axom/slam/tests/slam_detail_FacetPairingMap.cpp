// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/slam/mesh_struct/detail/FacetPairingMap.hpp"

#include <unordered_set>

using namespace axom;
using namespace axom::slam;

//------------------------------------------------------------------------------
// FacetKey Tests
//------------------------------------------------------------------------------

TEST(slam_detail_FacetKey, construction_2d)
{
  using KeyType = axom::slam::detail::FacetKey<2, int>;

  // Default construction
  KeyType key1;
  EXPECT_EQ(key1.key0, KeyType::INVALID);
  EXPECT_EQ(key1.key1, KeyType::INVALID);

  // Construction with single key (2D)
  KeyType key2(42);
  EXPECT_EQ(key2.key0, 42);
  EXPECT_EQ(key2.key1, KeyType::INVALID);

  // key1 parameter should be ignored in 2D
  KeyType key3(42, 100);
  EXPECT_EQ(key3.key0, 42);
  // key1 is present but not used for matching in 2D
}

TEST(slam_detail_FacetKey, construction_3d)
{
  using KeyType = axom::slam::detail::FacetKey<3, int>;

  // Default construction
  KeyType key1;
  EXPECT_EQ(key1.key0, KeyType::INVALID);
  EXPECT_EQ(key1.key1, KeyType::INVALID);

  // Construction with sorted keys
  KeyType key2(10, 20);
  EXPECT_EQ(key2.key0, 10);
  EXPECT_EQ(key2.key1, 20);

  // Construction with unsorted keys (should auto-sort)
  KeyType key3(20, 10);
  EXPECT_EQ(key3.key0, 10);  // Sorted
  EXPECT_EQ(key3.key1, 20);  // Sorted
}

TEST(slam_detail_FacetKey, equality_2d)
{
  using KeyType = axom::slam::detail::FacetKey<2, int>;

  KeyType key1(42);
  KeyType key2(42);
  KeyType key3(43);

  EXPECT_TRUE(key1 == key2);
  EXPECT_FALSE(key1 == key3);
  EXPECT_TRUE(key1 != key3);
}

TEST(slam_detail_FacetKey, equality_3d_with_sorting)
{
  using KeyType = axom::slam::detail::FacetKey<3, int>;

  // Same keys in different order should match
  KeyType key1(10, 20);
  KeyType key2(20, 10);  // Will be sorted to (10, 20)
  KeyType key3(10, 21);

  EXPECT_TRUE(key1 == key2);  // Should match after sorting
  EXPECT_FALSE(key1 == key3);
}

//------------------------------------------------------------------------------
// FacetData Tests
//------------------------------------------------------------------------------

TEST(slam_detail_FacetData, construction)
{
  using DataType = axom::slam::detail::FacetData<int>;

  // Default construction
  DataType data1;
  EXPECT_EQ(data1.element_idx, DataType::INVALID);
  EXPECT_EQ(data1.face_idx, DataType::INVALID);

  // Construction with values
  DataType data2(100, 2);
  EXPECT_EQ(data2.element_idx, 100);
  EXPECT_EQ(data2.face_idx, 2);
}

//------------------------------------------------------------------------------
// FacetPairingMap 2D Tests
//------------------------------------------------------------------------------

TEST(slam_detail_FacetPairingMap, basic_2d_insert_and_match)
{
  using MapType = axom::slam::detail::FacetPairingMap<2, int>;
  using KeyType = typename MapType::KeyType;
  using DataType = typename MapType::DataType;

  MapType map;
  map.prepareForInsertions(10);

  // Insert a facet
  KeyType key1(42);
  DataType data1(100, 1);
  map.insert(key1, data1);

  EXPECT_EQ(map.pendingCount(), 1);
  EXPECT_FALSE(map.allFacetsPaired());

  // Match with same key
  KeyType key2(42);
  auto match = map.findAndRemove(key2);

  ASSERT_TRUE(match.has_value());
  EXPECT_EQ(match->element_idx, 100);
  EXPECT_EQ(match->face_idx, 1);
  EXPECT_EQ(map.pendingCount(), 0);
  EXPECT_TRUE(map.allFacetsPaired());
}

TEST(slam_detail_FacetPairingMap, multiple_2d_facets)
{
  using MapType = axom::slam::detail::FacetPairingMap<2, int>;
  using KeyType = typename MapType::KeyType;
  using DataType = typename MapType::DataType;

  MapType map;
  map.prepareForInsertions(20);

  // Insert several facets
  for(int i = 0; i < 10; ++i)
  {
    KeyType key(i * 10);
    DataType data(i, 0);
    map.insert(key, data);
  }

  EXPECT_EQ(map.pendingCount(), 10);

  // Match them all
  for(int i = 0; i < 10; ++i)
  {
    KeyType key(i * 10);
    auto match = map.findAndRemove(key);

    ASSERT_TRUE(match.has_value());
    EXPECT_EQ(match->element_idx, i);
    EXPECT_EQ(match->face_idx, 0);
  }

  EXPECT_EQ(map.pendingCount(), 0);
  EXPECT_TRUE(map.allFacetsPaired());
}

TEST(slam_detail_FacetPairingMap, not_found_2d)
{
  using MapType = axom::slam::detail::FacetPairingMap<2, int>;
  using KeyType = typename MapType::KeyType;
  using DataType = typename MapType::DataType;

  MapType map;
  map.prepareForInsertions(10);

  // Insert one facet
  KeyType key1(42);
  DataType data1(100, 1);
  map.insert(key1, data1);

  // Try to find a different key
  KeyType key2(43);
  auto match = map.findAndRemove(key2);

  EXPECT_FALSE(match.has_value());
  EXPECT_EQ(map.pendingCount(), 1);  // Original still there
}

//------------------------------------------------------------------------------
// FacetPairingMap 3D Tests
//------------------------------------------------------------------------------

TEST(slam_detail_FacetPairingMap, basic_3d_insert_and_match)
{
  using MapType = axom::slam::detail::FacetPairingMap<3, int>;
  using KeyType = typename MapType::KeyType;
  using DataType = typename MapType::DataType;

  MapType map;
  map.prepareForInsertions(10);

  // Insert a facet
  KeyType key1(10, 20);
  DataType data1(100, 1);
  map.insert(key1, data1);

  EXPECT_EQ(map.pendingCount(), 1);

  // Match with equivalent key (different order)
  KeyType key2(20, 10);  // Should match (10, 20) after sorting
  auto match = map.findAndRemove(key2);

  ASSERT_TRUE(match.has_value());
  EXPECT_EQ(match->element_idx, 100);
  EXPECT_EQ(match->face_idx, 1);
  EXPECT_EQ(map.pendingCount(), 0);
}

TEST(slam_detail_FacetPairingMap, multiple_3d_facets)
{
  using MapType = axom::slam::detail::FacetPairingMap<3, int>;
  using KeyType = typename MapType::KeyType;
  using DataType = typename MapType::DataType;

  MapType map;
  map.prepareForInsertions(50);

  // Insert many facets to stress-test collision handling
  for(int i = 0; i < 50; ++i)
  {
    KeyType key(i, i + 1000);
    DataType data(i, 2);
    map.insert(key, data);
  }

  EXPECT_EQ(map.pendingCount(), 50);

  // Match them all (including with reversed key order)
  for(int i = 0; i < 50; ++i)
  {
    // Use reversed order to test sorting
    KeyType key(i + 1000, i);
    auto match = map.findAndRemove(key);

    ASSERT_TRUE(match.has_value()) << "Failed to find key (" << i << ", " << (i + 1000) << ")";
    EXPECT_EQ(match->element_idx, i);
    EXPECT_EQ(match->face_idx, 2);
  }

  EXPECT_EQ(map.pendingCount(), 0);
  EXPECT_TRUE(map.allFacetsPaired());
}

TEST(slam_detail_FacetPairingMap, 3d_key_ordering_invariant)
{
  using MapType = axom::slam::detail::FacetPairingMap<3, int>;
  using KeyType = typename MapType::KeyType;
  using DataType = typename MapType::DataType;

  MapType map;
  map.prepareForInsertions(10);

  // Insert with keys in sorted order
  KeyType key1(5, 10);
  DataType data1(100, 1);
  map.insert(key1, data1);

  // Should match with reversed order
  KeyType key2(10, 5);
  auto match = map.findAndRemove(key2);

  ASSERT_TRUE(match.has_value());
  EXPECT_EQ(match->element_idx, 100);

  // Insert again with reversed order
  KeyType key3(15, 10);  // Will be sorted to (10, 15)
  DataType data2(200, 2);
  map.insert(key3, data2);

  // Match with normal order
  KeyType key4(10, 15);
  auto match2 = map.findAndRemove(key4);

  ASSERT_TRUE(match2.has_value());
  EXPECT_EQ(match2->element_idx, 200);
}

//------------------------------------------------------------------------------
// Collision and Tombstone Tests
//------------------------------------------------------------------------------

TEST(slam_detail_FacetPairingMap, collision_handling)
{
  using MapType = axom::slam::detail::FacetPairingMap<3, int>;
  using KeyType = typename MapType::KeyType;
  using DataType = typename MapType::DataType;

  MapType map;
  map.prepareForInsertions(100);

  // Insert many facets to force collisions
  std::vector<KeyType> keys;
  for(int i = 0; i < 100; ++i)
  {
    KeyType key(i, i + 5000);
    keys.push_back(key);
    DataType data(i, 0);
    map.insert(key, data);
  }

  EXPECT_EQ(map.pendingCount(), 100);

  // Match all keys in random order
  std::vector<int> indices(100);
  for(int i = 0; i < 100; ++i) indices[i] = i;
  // Simple shuffle
  for(int i = 0; i < 100; ++i)
  {
    int j = (i * 7 + 13) % 100;
    std::swap(indices[i], indices[j]);
  }

  for(int idx : indices)
  {
    auto match = map.findAndRemove(keys[idx]);
    ASSERT_TRUE(match.has_value()) << "Failed to find key at index " << idx;
    EXPECT_EQ(match->element_idx, idx);
  }

  EXPECT_EQ(map.pendingCount(), 0);
}

TEST(slam_detail_FacetPairingMap, tombstone_reuse)
{
  using MapType = axom::slam::detail::FacetPairingMap<2, int>;
  using KeyType = typename MapType::KeyType;
  using DataType = typename MapType::DataType;

  MapType map;
  map.prepareForInsertions(10);

  // Insert and remove to create tombstones
  for(int i = 0; i < 5; ++i)
  {
    KeyType key1(i * 10);
    DataType data1(i, 0);
    map.insert(key1, data1);

    KeyType key2(i * 10);
    auto match = map.findAndRemove(key2);
    ASSERT_TRUE(match.has_value());
  }

  EXPECT_EQ(map.pendingCount(), 0);

  // Now insert new keys - should reuse tombstone slots
  for(int i = 0; i < 5; ++i)
  {
    KeyType key(i * 10 + 5);  // Different keys
    DataType data(i + 100, 1);
    map.insert(key, data);
  }

  EXPECT_EQ(map.pendingCount(), 5);

  // Verify new keys are findable
  for(int i = 0; i < 5; ++i)
  {
    KeyType key(i * 10 + 5);
    auto match = map.findAndRemove(key);
    ASSERT_TRUE(match.has_value());
    EXPECT_EQ(match->element_idx, i + 100);
  }
}

//------------------------------------------------------------------------------
// Generation and Reuse Tests
//------------------------------------------------------------------------------

TEST(slam_detail_FacetPairingMap, multiple_prepare_cycles)
{
  using MapType = axom::slam::detail::FacetPairingMap<2, int>;
  using KeyType = typename MapType::KeyType;
  using DataType = typename MapType::DataType;

  MapType map;

  // Run multiple prepare/insert/match cycles
  for(int cycle = 0; cycle < 10; ++cycle)
  {
    map.prepareForInsertions(10);

    // Insert facets
    for(int i = 0; i < 5; ++i)
    {
      KeyType key(cycle * 100 + i);
      DataType data(i, cycle);
      map.insert(key, data);
    }

    EXPECT_EQ(map.pendingCount(), 5);

    // Match them
    for(int i = 0; i < 5; ++i)
    {
      KeyType key(cycle * 100 + i);
      auto match = map.findAndRemove(key);
      ASSERT_TRUE(match.has_value());
      EXPECT_EQ(match->element_idx, i);
      EXPECT_EQ(match->face_idx, cycle);
    }

    EXPECT_EQ(map.pendingCount(), 0);
  }
}

TEST(slam_detail_FacetPairingMap, generation_isolation)
{
  using MapType = axom::slam::detail::FacetPairingMap<2, int>;
  using KeyType = typename MapType::KeyType;
  using DataType = typename MapType::DataType;

  MapType map;

  // First generation
  map.prepareForInsertions(10);
  KeyType key1(42);
  DataType data1(100, 1);
  map.insert(key1, data1);
  EXPECT_EQ(map.pendingCount(), 1);

  // Second generation (should not see first generation's data)
  map.prepareForInsertions(10);
  EXPECT_EQ(map.pendingCount(), 0);  // Reset

  // Try to find old key - should not be found
  KeyType key2(42);
  auto match = map.findAndRemove(key2);
  EXPECT_FALSE(match.has_value());
}

//------------------------------------------------------------------------------
// Table Sizing Tests
//------------------------------------------------------------------------------

TEST(slam_detail_FacetPairingMap, table_sizing_2d)
{
  using MapType = axom::slam::detail::FacetPairingMap<2, int>;

  MapType map;

  // Small count
  map.prepareForInsertions(10);
  EXPECT_GE(map.tableSize(), 8);  // Minimum size

  // Larger count - should scale appropriately
  map.prepareForInsertions(100);
  EXPECT_GE(map.tableSize(), 400);  // 4x in 2D

  // Table size should be power of 2
  EXPECT_EQ(map.tableSize() & (map.tableSize() - 1), 0);
}

TEST(slam_detail_FacetPairingMap, table_sizing_3d)
{
  using MapType = axom::slam::detail::FacetPairingMap<3, int>;

  MapType map;

  // Small count
  map.prepareForInsertions(10);
  EXPECT_GE(map.tableSize(), 8);

  // Larger count - should scale appropriately
  map.prepareForInsertions(100);
  EXPECT_GE(map.tableSize(), 800);  // 8x in 3D

  // Table size should be power of 2
  EXPECT_EQ(map.tableSize() & (map.tableSize() - 1), 0);
}

//------------------------------------------------------------------------------
// Edge Cases
//------------------------------------------------------------------------------

TEST(slam_detail_FacetPairingMap, empty_map_operations)
{
  using MapType = axom::slam::detail::FacetPairingMap<2, int>;
  using KeyType = typename MapType::KeyType;

  MapType map;
  map.prepareForInsertions(10);

  // Try to find in empty map
  KeyType key(42);
  auto match = map.findAndRemove(key);

  EXPECT_FALSE(match.has_value());
  EXPECT_EQ(map.pendingCount(), 0);
  EXPECT_TRUE(map.allFacetsPaired());
}

TEST(slam_detail_FacetPairingMap, single_facet)
{
  using MapType = axom::slam::detail::FacetPairingMap<2, int>;
  using KeyType = typename MapType::KeyType;
  using DataType = typename MapType::DataType;

  MapType map;
  map.prepareForInsertions(1);

  KeyType key(42);
  DataType data(100, 1);
  map.insert(key, data);

  EXPECT_EQ(map.pendingCount(), 1);

  auto match = map.findAndRemove(key);
  ASSERT_TRUE(match.has_value());
  EXPECT_EQ(map.pendingCount(), 0);
}

TEST(slam_detail_FacetPairingMap, large_vertex_ids)
{
  using MapType = axom::slam::detail::FacetPairingMap<3, int>;
  using KeyType = typename MapType::KeyType;
  using DataType = typename MapType::DataType;

  MapType map;
  map.prepareForInsertions(10);

  // Use large vertex IDs to test hash collision behavior
  const int LARGE_ID = 1000000;

  for(int i = 0; i < 10; ++i)
  {
    KeyType key(LARGE_ID + i, LARGE_ID + i + 1);
    DataType data(i, 0);
    map.insert(key, data);
  }

  EXPECT_EQ(map.pendingCount(), 10);

  for(int i = 0; i < 10; ++i)
  {
    KeyType key(LARGE_ID + i, LARGE_ID + i + 1);
    auto match = map.findAndRemove(key);
    ASSERT_TRUE(match.has_value());
  }

  EXPECT_TRUE(map.allFacetsPaired());
}

//------------------------------------------------------------------------------
// main
//------------------------------------------------------------------------------

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  int result = RUN_ALL_TESTS();

  return result;
}
