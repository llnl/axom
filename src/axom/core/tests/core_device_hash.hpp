// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

// Axom includes
#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/DeviceHash.hpp"

// gtest includes
#include "gtest/gtest.h"

// C++ includes
#include <cstdint>
#include <set>
#include <type_traits>

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
  using HashResult = typename decltype(device_hasher)::result_type;

  constexpr int NUM_HASHES = 4;

  int things_to_hash[NUM_HASHES] {0, 1, 37, 1100};

  // Allocate space for hash results.
  int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
  HashResult* computed_hashes = axom::allocate<HashResult>(NUM_HASHES, allocatorID);

  // Compute hashes.
  axom::for_all<ExecSpace>(
    NUM_HASHES,
    AXOM_LAMBDA(int i) { computed_hashes[i] = device_hasher(things_to_hash[i]); });

  // Copy back to host.
  HashResult computed_hashes_host[NUM_HASHES];
  axom::copy(computed_hashes_host, computed_hashes, sizeof(HashResult) * NUM_HASHES);
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
  using HashResult = typename decltype(device_hasher)::result_type;

  constexpr int NUM_HASHES = 4;

  float things_to_hash[NUM_HASHES] {0.f, 1.f, 37.f, 1100.f};

  // Allocate space for hash results.
  int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
  HashResult* computed_hashes = axom::allocate<HashResult>(NUM_HASHES, allocatorID);

  // Compute hashes.
  axom::for_all<ExecSpace>(
    NUM_HASHES,
    AXOM_LAMBDA(int i) { computed_hashes[i] = device_hasher(things_to_hash[i]); });

  // Copy back to host.
  HashResult computed_hashes_host[NUM_HASHES];
  axom::copy(computed_hashes_host, computed_hashes, sizeof(HashResult) * NUM_HASHES);
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
  using HashResult = typename decltype(device_hasher)::result_type;

  constexpr int NUM_HASHES = 4;

  std::string things_to_hash[NUM_HASHES] {"0", "1", "37", "1100"};

  HashResult computed_hashes[NUM_HASHES];

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
  using HashResult = typename decltype(device_hasher)::result_type;

  constexpr int NUM_HASHES = 4;

  TestEnumHash things_to_hash[NUM_HASHES] {TestEnumHash::Zero,
                                           TestEnumHash::One,
                                           TestEnumHash::Two,
                                           TestEnumHash::Three};

  // Allocate space for hash results.
  int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
  HashResult* computed_hashes = axom::allocate<HashResult>(NUM_HASHES, allocatorID);

  // Compute hashes.
  axom::for_all<ExecSpace>(
    NUM_HASHES,
    AXOM_LAMBDA(int i) { computed_hashes[i] = device_hasher(things_to_hash[i]); });

  // Copy back to host.
  HashResult computed_hashes_host[NUM_HASHES];
  axom::copy(computed_hashes_host, computed_hashes, sizeof(HashResult) * NUM_HASHES);
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
    auto ptr = reinterpret_cast<axom_testing::UserVector<T>*>(bytes);
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
  using HashResult = typename decltype(device_hasher)::result_type;

  constexpr int NUM_HASHES = 4;

  axom_testing::UserVector<float> things_to_hash[NUM_HASHES] = {{0.0, 0.0, 0.0},
                                                                {1.0, 3.0, 5.0},
                                                                {2.0, 5.0, 8.0},
                                                                {10.0, 20.0, 30.0}};

  // Allocate space for hash results.
  int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
  HashResult* computed_hashes = axom::allocate<HashResult>(NUM_HASHES, allocatorID);

  // Compute hashes.
  axom::for_all<ExecSpace>(
    NUM_HASHES,
    AXOM_LAMBDA(int i) { computed_hashes[i] = device_hasher(things_to_hash[i]); });

  // Copy back to host.
  HashResult computed_hashes_host[NUM_HASHES];
  axom::copy(computed_hashes_host, computed_hashes, sizeof(HashResult) * NUM_HASHES);
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

AXOM_TYPED_TEST(core_device_hash, hash_uint64_distinguishes_high_bits)
{
  using ExecSpace = typename TestFixture::ExecSpace;

  axom::DeviceHash<std::uint64_t> device_hasher;
  using HashResult = typename decltype(device_hasher)::result_type;

  constexpr int NUM_HASHES = 3;
  std::uint64_t things_to_hash[NUM_HASHES] = {std::uint64_t {1},
                                              std::uint64_t {1} + (std::uint64_t {1} << 32),
                                              std::uint64_t {1} + (std::uint64_t {1} << 33)};

  int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
  HashResult* computed_hashes = axom::allocate<HashResult>(NUM_HASHES, allocatorID);

  axom::for_all<ExecSpace>(
    NUM_HASHES,
    AXOM_LAMBDA(int i) { computed_hashes[i] = device_hasher(things_to_hash[i]); });

  HashResult computed_hashes_host[NUM_HASHES];
  axom::copy(computed_hashes_host, computed_hashes, sizeof(HashResult) * NUM_HASHES);
  axom::deallocate(computed_hashes);

  EXPECT_NE(computed_hashes_host[0], computed_hashes_host[1]);
  EXPECT_NE(computed_hashes_host[0], computed_hashes_host[2]);
  EXPECT_NE(computed_hashes_host[1], computed_hashes_host[2]);
}

AXOM_TYPED_TEST(core_device_hash, hash_fractional_float_device)
{
  using ExecSpace = typename TestFixture::ExecSpace;

  axom::DeviceHash<float> device_hasher;
  using HashResult = typename decltype(device_hasher)::result_type;

  constexpr int NUM_HASHES = 8;
  float things_to_hash[NUM_HASHES] = {0.25f, 0.75f, -0.5f, 0.5f, 0.125f, 0.625f, 0.875f, 1.25f};

  int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
  HashResult* computed_hashes = axom::allocate<HashResult>(NUM_HASHES, allocatorID);

  axom::for_all<ExecSpace>(
    NUM_HASHES,
    AXOM_LAMBDA(int i) { computed_hashes[i] = device_hasher(things_to_hash[i]); });

  HashResult computed_hashes_host[NUM_HASHES];
  axom::copy(computed_hashes_host, computed_hashes, sizeof(HashResult) * NUM_HASHES);
  axom::deallocate(computed_hashes);

  // Idempotence and pairwise distinctness for these chosen values.
  for(int i = 0; i < NUM_HASHES; i++)
  {
    EXPECT_EQ(computed_hashes_host[i], device_hasher(things_to_hash[i]));
    for(int j = i + 1; j < NUM_HASHES; j++)
    {
      EXPECT_NE(computed_hashes_host[i], computed_hashes_host[j]);
    }
  }

  EXPECT_EQ(device_hasher(0.0f), device_hasher(-0.0f));
}

AXOM_TYPED_TEST(core_device_hash, hash_fractional_double_device)
{
  using ExecSpace = typename TestFixture::ExecSpace;

  axom::DeviceHash<double> device_hasher;
  using HashResult = typename decltype(device_hasher)::result_type;

  constexpr int NUM_HASHES = 8;
  double things_to_hash[NUM_HASHES] = {0.25, 0.75, -0.5, 0.5, 0.125, 0.625, 0.875, 1.25};

  int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
  HashResult* computed_hashes = axom::allocate<HashResult>(NUM_HASHES, allocatorID);

  axom::for_all<ExecSpace>(
    NUM_HASHES,
    AXOM_LAMBDA(int i) { computed_hashes[i] = device_hasher(things_to_hash[i]); });

  HashResult computed_hashes_host[NUM_HASHES];
  axom::copy(computed_hashes_host, computed_hashes, sizeof(HashResult) * NUM_HASHES);
  axom::deallocate(computed_hashes);

  for(int i = 0; i < NUM_HASHES; i++)
  {
    EXPECT_EQ(computed_hashes_host[i], device_hasher(things_to_hash[i]));
    for(int j = i + 1; j < NUM_HASHES; j++)
    {
      EXPECT_NE(computed_hashes_host[i], computed_hashes_host[j]);
    }
  }

  EXPECT_EQ(device_hasher(0.0), device_hasher(-0.0));
}

TEST(core_device_hash, hash_width_decoupled_from_indextype)
{
  // The hash result must be 64 bits wide regardless of the configured
  // axom::IndexType. When the result type was IndexType, builds with
  // AXOM_USE_64BIT_INDEXTYPE=OFF truncated integer keys to 32 bits before
  // the FlatMap bit mixer ran, so keys equal mod 2^32 (e.g. deep Morton codes)
  // produced identical hashes. The type assertions catch the coupling in every
  // build configuration; the value checks fail in the truncating configuration itself.
  static_assert(std::is_same<axom::DeviceHash<std::uint64_t>::result_type, std::uint64_t>::value,
                "integral hash result must be std::uint64_t");
  static_assert(std::is_same<axom::DeviceHash<int>::result_type, std::uint64_t>::value,
                "integral hash result must be std::uint64_t");
  static_assert(std::is_same<axom::DeviceHash<double>::result_type, std::uint64_t>::value,
                "floating-point hash result must be std::uint64_t");
  static_assert(std::is_same<axom::DeviceHash<int*>::result_type, std::uint64_t>::value,
                "pointer hash result must be std::uint64_t");
  static_assert(std::is_same<axom::DeviceHash<std::string>::result_type, std::uint64_t>::value,
                "catch-all (std::hash) result must be std::uint64_t");
  static_assert(
    std::is_same<decltype(axom::DeviceHash<std::uint64_t> {}(std::uint64_t {})), std::uint64_t>::value,
    "integral hash operator() must return std::uint64_t");

  axom::DeviceHash<std::uint64_t> device_hasher;
  const std::uint64_t base = 1;
  const std::uint64_t plus_2_32 = base + (std::uint64_t {1} << 32);
  const std::uint64_t plus_2_33 = base + (std::uint64_t {1} << 33);
  EXPECT_NE(device_hasher(base), device_hasher(plus_2_32));
  EXPECT_NE(device_hasher(base), device_hasher(plus_2_33));
  EXPECT_NE(device_hasher(plus_2_32), device_hasher(plus_2_33));
}

TEST(core_device_hash, hash_float_bit_pattern)
{
  // Floating-point keys must be hashed by bit pattern, not by integer value conversion.
  // This is a regression test for a previous implementation where the conversion collapsed
  // every key with the same integer value, e.g. all numbers between -1 and 1 converted to integer 0
  // so a FlatMap keyed on fractional floats degenerated into a single probe chain.
  axom::DeviceHash<float> float_hasher;
  axom::DeviceHash<double> double_hasher;

  EXPECT_NE(float_hasher(0.25f), float_hasher(0.75f));
  EXPECT_NE(float_hasher(0.25f), std::uint64_t {0});
  EXPECT_NE(double_hasher(0.25), double_hasher(0.75));

  // A spread of fractional keys must be collision-free at this scale
  std::set<std::uint64_t> float_hashes, double_hashes;
  const int NUM_KEYS = 1000;
  for(int i = 1; i <= NUM_KEYS; i++)
  {
    float_hashes.insert(
      float_hasher(static_cast<float>(i) / static_cast<float>(NUM_KEYS + 1)));
    double_hashes.insert(double_hasher(i / static_cast<double>(NUM_KEYS + 1)));
  }
  EXPECT_EQ(float_hashes.size(), NUM_KEYS);
  EXPECT_EQ(double_hashes.size(), NUM_KEYS);

  // Signed zeros compare equal and must hash equal
  EXPECT_EQ(float_hasher(0.0f), float_hasher(-0.0f));
  EXPECT_EQ(double_hasher(0.0), double_hasher(-0.0));

  // Magnitudes beyond any integer type's range are now well-defined and distinct
  EXPECT_NE(double_hasher(1e300), double_hasher(2e300));
  EXPECT_NE(float_hasher(-0.5f), float_hasher(0.5f));
}

TEST(core_device_hash, hash_long_double_has_stable_equal_value_hash)
{
  axom::DeviceHash<long double> long_double_hasher;

  static_assert(std::is_same<axom::DeviceHash<long double>::result_type, std::uint64_t>::value,
                "long double hash result must be std::uint64_t");

  EXPECT_EQ(long_double_hasher(0.0L), long_double_hasher(-0.0L));
  EXPECT_EQ(long_double_hasher(0.25L), long_double_hasher(static_cast<long double>(0.25)));
  EXPECT_NE(long_double_hasher(0.25L), long_double_hasher(0.75L));
}
