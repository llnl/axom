// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/DeviceHash.hpp"

// gtest includes
#include "gtest/gtest.h"

TEST(core_device_hash, hash_int)
{
  axom::DeviceHash<int> device_hasher;

  constexpr int NUM_HASHES = 4;

  int things_to_hash[NUM_HASHES] {0, 1, 37, 1100};

  axom::IndexType computed_hashes[NUM_HASHES];

  // Compute hashes.
  for(int i = 0; i < NUM_HASHES; i++)
  {
    computed_hashes[i] = device_hasher(things_to_hash[i]);
  }

  for(int i = 0; i < NUM_HASHES; i++)
  {
    // Invocations of the hash function should be idempotent.
    EXPECT_EQ(computed_hashes[i], device_hasher(things_to_hash[i]));

    // Check that we don't have hash collisions with other values.
    for(int j = i + 1; j < NUM_HASHES; j++)
    {
      EXPECT_NE(computed_hashes[i], computed_hashes[j]);
    }
  }
}

TEST(core_device_hash, hash_float)
{
  axom::DeviceHash<float> device_hasher;

  constexpr int NUM_HASHES = 4;

  float things_to_hash[NUM_HASHES] {0.f, 1.f, 37.f, 1100.f};

  axom::IndexType computed_hashes[NUM_HASHES];

  // Compute hashes.
  for(int i = 0; i < NUM_HASHES; i++)
  {
    computed_hashes[i] = device_hasher(things_to_hash[i]);
  }

  for(int i = 0; i < NUM_HASHES; i++)
  {
    // Invocations of the hash function should be idempotent.
    EXPECT_EQ(computed_hashes[i], device_hasher(things_to_hash[i]));

    // Check that we don't have hash collisions with other values.
    for(int j = i + 1; j < NUM_HASHES; j++)
    {
      EXPECT_NE(computed_hashes[i], computed_hashes[j]);
    }
  }

  // Since 0.0f == -0.0f, they should hash to the same value.
  EXPECT_EQ(device_hasher(0.0f), device_hasher(-0.0f));
}

TEST(core_device_hash, hash_string)
{
  axom::DeviceHash<std::string> device_hasher;

  constexpr int NUM_HASHES = 4;

  std::string things_to_hash[NUM_HASHES] {"0", "1", "37", "1100"};

  axom::IndexType computed_hashes[NUM_HASHES];

  // Compute hashes.
  for(int i = 0; i < NUM_HASHES; i++)
  {
    computed_hashes[i] = device_hasher(things_to_hash[i]);
  }

  for(int i = 0; i < NUM_HASHES; i++)
  {
    // Invocations of the hash function should be idempotent.
    EXPECT_EQ(computed_hashes[i], device_hasher(things_to_hash[i]));

    // Check that we don't have hash collisions with other values.
    for(int j = i + 1; j < NUM_HASHES; j++)
    {
      EXPECT_NE(computed_hashes[i], computed_hashes[j]);
    }
  }
}

enum class TestEnumHash
{
  Zero,
  One,
  Two,
  Three
};

TEST(core_device_hash, hash_enum)
{
  axom::DeviceHash<TestEnumHash> device_hasher;

  constexpr int NUM_HASHES = 4;

  TestEnumHash things_to_hash[NUM_HASHES] {TestEnumHash::Zero,
                                           TestEnumHash::One,
                                           TestEnumHash::Two,
                                           TestEnumHash::Three};

  axom::IndexType computed_hashes[NUM_HASHES];

  // Compute hashes.
  for(int i = 0; i < NUM_HASHES; i++)
  {
    computed_hashes[i] = device_hasher(things_to_hash[i]);
  }

  for(int i = 0; i < NUM_HASHES; i++)
  {
    // Invocations of the hash function should be idempotent.
    EXPECT_EQ(computed_hashes[i], device_hasher(things_to_hash[i]));

    // Check that we don't have hash collisions with other values.
    for(int j = i + 1; j < NUM_HASHES; j++)
    {
      EXPECT_NE(computed_hashes[i], computed_hashes[j]);
    }
  }
}
