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

template <typename TheExecSpace>
class core_device_hash : public ::testing::Test
{
public:
  using ExecSpace = TheExecSpace;
};

using HashTestTypes = ::testing::Types<
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
  axom::CUDA_EXEC<256>,
#endif
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(AXOM_USE_UMPIRE)
  axom::HIP_EXEC<256>,
#endif
  axom::SEQ_EXEC>;

TYPED_TEST_SUITE(core_device_hash, HashTestTypes);

AXOM_TYPED_TEST(core_device_hash, hash_int)
{
  using ExecSpace = typename TestFixture::ExecSpace;

  axom::DeviceHash<int> device_hasher;

  constexpr int NUM_HASHES = 4;

  int things_to_hash[NUM_HASHES] {0, 1, 37, 1100};

  // Allocate space for hash results.
  int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
  axom::IndexType *computed_hashes = axom::allocate<axom::IndexType>(NUM_HASHES, allocatorID);

  // Compute hashes.
  axom::for_all<ExecSpace>(
    NUM_HASHES,
    AXOM_LAMBDA(int i) { computed_hashes[i] = device_hasher(things_to_hash[i]); });

  // Copy back to host.
  axom::IndexType computed_hashes_host[NUM_HASHES];
  axom::copy(computed_hashes_host, computed_hashes, sizeof(axom::IndexType) * NUM_HASHES);
  axom::deallocate(computed_hashes);

  for(int i = 0; i < NUM_HASHES; i++)
  {
    // Invocations of the hash function should be idempotent.
    EXPECT_EQ(computed_hashes_host[i], device_hasher(things_to_hash[i]));

    // Check that we don't have hash collisions with other values.
    for(int j = i + 1; j < NUM_HASHES; j++)
    {
      EXPECT_NE(computed_hashes_host[i], computed_hashes_host[j]);
    }
  }
}

AXOM_TYPED_TEST(core_device_hash, hash_float)
{
  using ExecSpace = typename TestFixture::ExecSpace;

  axom::DeviceHash<float> device_hasher;

  constexpr int NUM_HASHES = 4;

  float things_to_hash[NUM_HASHES] {0.f, 1.f, 37.f, 1100.f};

  // Allocate space for hash results.
  int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
  axom::IndexType *computed_hashes = axom::allocate<axom::IndexType>(NUM_HASHES, allocatorID);

  // Compute hashes.
  axom::for_all<ExecSpace>(
    NUM_HASHES,
    AXOM_LAMBDA(int i) { computed_hashes[i] = device_hasher(things_to_hash[i]); });

  // Copy back to host.
  axom::IndexType computed_hashes_host[NUM_HASHES];
  axom::copy(computed_hashes_host, computed_hashes, sizeof(axom::IndexType) * NUM_HASHES);
  axom::deallocate(computed_hashes);

  for(int i = 0; i < NUM_HASHES; i++)
  {
    // Invocations of the hash function should be idempotent.
    EXPECT_EQ(computed_hashes_host[i], device_hasher(things_to_hash[i]));

    // Check that we don't have hash collisions with other values.
    for(int j = i + 1; j < NUM_HASHES; j++)
    {
      EXPECT_NE(computed_hashes_host[i], computed_hashes_host[j]);
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

AXOM_TYPED_TEST(core_device_hash, hash_enum)
{
  using ExecSpace = typename TestFixture::ExecSpace;

  axom::DeviceHash<TestEnumHash> device_hasher;

  constexpr int NUM_HASHES = 4;

  TestEnumHash things_to_hash[NUM_HASHES] {TestEnumHash::Zero,
                                           TestEnumHash::One,
                                           TestEnumHash::Two,
                                           TestEnumHash::Three};

  // Allocate space for hash results.
  int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
  axom::IndexType *computed_hashes = axom::allocate<axom::IndexType>(NUM_HASHES, allocatorID);

  // Compute hashes.
  axom::for_all<ExecSpace>(
    NUM_HASHES,
    AXOM_LAMBDA(int i) { computed_hashes[i] = device_hasher(things_to_hash[i]); });

  // Copy back to host.
  axom::IndexType computed_hashes_host[NUM_HASHES];
  axom::copy(computed_hashes_host, computed_hashes, sizeof(axom::IndexType) * NUM_HASHES);
  axom::deallocate(computed_hashes);

  for(int i = 0; i < NUM_HASHES; i++)
  {
    // Invocations of the hash function should be idempotent.
    EXPECT_EQ(computed_hashes_host[i], device_hasher(things_to_hash[i]));

    // Check that we don't have hash collisions with other values.
    for(int j = i + 1; j < NUM_HASHES; j++)
    {
      EXPECT_NE(computed_hashes_host[i], computed_hashes_host[j]);
    }
  }
}

namespace axom_testing
{
template <typename T>
struct UserVector
{
  T x, y, z;
};
}  // namespace axom_testing

// Test that we can correctly specialize a device hash for a user-defined type.
namespace axom
{
template <typename T>
struct DeviceHash<axom_testing::UserVector<T>>
{
  using argument_type = axom_testing::UserVector<T>;
  using result_type = axom::IndexType;

  AXOM_HOST_DEVICE axom::IndexType operator()(axom_testing::UserVector<T> value) const
  {
    // Copy byte representation over
    constexpr int NWORDS = sizeof(axom_testing::UserVector<T>) / sizeof(int);
    alignas(axom_testing::UserVector<T>) int bytes[NWORDS];
    // NOTE: Separating these statements fixes a warning about strict-aliasing.
    auto ptr = reinterpret_cast<axom_testing::UserVector<T> *>(bytes);
    *ptr = value;

    axom::IndexType hash_result {};
    for(int i = 0; i < NWORDS; i++)
    {
      hash_result ^= (bytes[i] + 0x853c49e6);
    }
    return hash_result;
  }
};
}  // namespace axom

AXOM_TYPED_TEST(core_device_hash, hash_user_defined)
{
  using ExecSpace = typename TestFixture::ExecSpace;

  axom::DeviceHash<axom_testing::UserVector<float>> device_hasher;

  constexpr int NUM_HASHES = 4;

  axom_testing::UserVector<float> things_to_hash[NUM_HASHES] = {{0.0, 0.0, 0.0},
                                                                {1.0, 3.0, 5.0},
                                                                {2.0, 5.0, 8.0},
                                                                {10.0, 20.0, 30.0}};

  // Allocate space for hash results.
  int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
  axom::IndexType *computed_hashes = axom::allocate<axom::IndexType>(NUM_HASHES, allocatorID);

  // Compute hashes.
  axom::for_all<ExecSpace>(
    NUM_HASHES,
    AXOM_LAMBDA(int i) { computed_hashes[i] = device_hasher(things_to_hash[i]); });

  // Copy back to host.
  axom::IndexType computed_hashes_host[NUM_HASHES];
  axom::copy(computed_hashes_host, computed_hashes, sizeof(axom::IndexType) * NUM_HASHES);
  axom::deallocate(computed_hashes);

  for(int i = 0; i < NUM_HASHES; i++)
  {
    // Invocations of the hash function should be idempotent.
    EXPECT_EQ(computed_hashes_host[i], device_hasher(things_to_hash[i]));

    // Check that we don't have hash collisions with other values.
    for(int j = i + 1; j < NUM_HASHES; j++)
    {
      EXPECT_NE(computed_hashes_host[i], computed_hashes_host[j]);
    }
  }
}
