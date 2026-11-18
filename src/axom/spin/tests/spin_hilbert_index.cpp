// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file spin_hilbert_index.cpp
 * \brief Tests the Hilbert curve convention and reversible integer conversions.
 */

#include "axom/spin/HilbertIndex.hpp"
#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/slic.hpp"
#include "gtest/gtest.h"

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <random>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{
// Fixed seed for the sampled points and keys, so a failure reproduces.
// The value is arbitrary.
constexpr std::uint64_t TEST_SEED = 42;

template <int DIM>
void exhaustive(int bits)
{
  using H = axom::spin::Hilbertizer<std::int64_t, std::uint64_t, DIM>;
  const std::uint64_t side = std::uint64_t {1} << bits;
  const std::uint64_t count = std::uint64_t {1} << (bits * DIM);
  std::vector<bool> seen(count, false);

  // Enumerate the Cartesian grid independently of the Hilbert decoder.
  for(std::uint64_t ordinal = 0; ordinal < count; ++ordinal)
  {
    SCOPED_TRACE(::testing::Message() << "dimension " << DIM << ", point " << ordinal);
    typename H::PointType point;
    auto remaining = ordinal;
    for(int axis = 0; axis < DIM; ++axis)
    {
      point[axis] = static_cast<std::int64_t>(remaining % side);
      remaining /= side;
    }
    const auto key = H::toIndex(point, bits);
    ASSERT_LT(key, count);
    ASSERT_FALSE(seen[key]);
    seen[key] = true;
    EXPECT_EQ(H::toPoint(key, bits), point);
  }

  typename H::PointType last;
  for(std::uint64_t key = 0; key < count; ++key)
  {
    SCOPED_TRACE(::testing::Message() << "dimension " << DIM << ", key " << key);
    EXPECT_TRUE(seen[key]);
    const auto point = H::toPoint(key, bits);
    ASSERT_TRUE(H::isValidPoint(point, bits));
    EXPECT_EQ(H::toIndex(point, bits), key);
    if(key > 0)
    {
      std::int64_t distance = 0;
      for(int axis = 0; axis < DIM; ++axis)
      {
        distance += std::abs(point[axis] - last[axis]);
      }
      EXPECT_EQ(distance, 1);
    }
    if(bits > 1)
    {
      auto parent = point;
      for(int axis = 0; axis < DIM; ++axis)
      {
        parent[axis] >>= 1;
      }
      EXPECT_EQ(H::toIndex(parent, bits - 1), key >> DIM);
    }
    last = point;
  }
}

template <typename Coord, typename Index, int DIM>
void allDepths()
{
  using H = axom::spin::Hilbertizer<Coord, Index, DIM>;
  // Keys depend only on the point and the depth, so each type pair must agree
  // with the widest pair at every depth that both support.
  using Reference = axom::spin::Hilbertizer<std::uint64_t, std::uint64_t, DIM>;
  const auto widen = [](const typename H::PointType& point) {
    typename Reference::PointType wide;
    for(int axis = 0; axis < DIM; ++axis)
    {
      wide[axis] = static_cast<std::uint64_t>(point[axis]);
    }
    return wide;
  };
  constexpr int coord_bits = std::numeric_limits<Coord>::digits;
  constexpr int index_bits = std::numeric_limits<Index>::digits;
  constexpr int capacity = coord_bits < index_bits / DIM ? coord_bits : index_bits / DIM;
  EXPECT_EQ(H::COORD_BITS, coord_bits);
  EXPECT_EQ(H::INDEX_BITS, index_bits);
  EXPECT_EQ(H::MAX_UNIQUE_BITS, capacity);
  EXPECT_EQ(H::maxBitsPerCoord(), capacity);
  std::mt19937_64 random(TEST_SEED);
  for(int bits = 1; bits <= capacity; ++bits)
  {
    SCOPED_TRACE(::testing::Message()
                 << "dimension " << DIM << ", depth " << bits << ", coordinate digits "
                 << coord_bits << ", index digits " << index_bits);
    const std::uint64_t max_coord = (std::uint64_t {1} << bits) - 1;
    const auto max_key = std::numeric_limits<std::uint64_t>::max() >> (64 - bits * DIM);
    EXPECT_TRUE(H::isValidDepth(bits));
    EXPECT_TRUE(H::isValidIndex(static_cast<Index>(max_key), bits));
    for(int iteration = 0; iteration < 100; ++iteration)
    {
      typename H::PointType point;
      for(int axis = 0; axis < DIM; ++axis)
      {
        // Check every corner, then sample points with a fixed seed.
        point[axis] =
          static_cast<Coord>(iteration < (1 << DIM) ? (((iteration >> axis) & 1) ? max_coord : 0)
                                                    : (random() & max_coord));
      }
      ASSERT_TRUE(H::isValidPoint(point, bits));
      const auto encoded = H::toIndex(point, bits);
      ASSERT_TRUE(H::isValidIndex(encoded, bits));
      EXPECT_EQ(H::toPoint(encoded, bits), point);
      EXPECT_EQ(axom::spin::convertPointToHilbert<Index>(point, bits), encoded);
      EXPECT_EQ((axom::spin::convertHilbertToPoint<Coord, DIM>(encoded, bits)), point);
      EXPECT_EQ(static_cast<std::uint64_t>(encoded), Reference::toIndex(widen(point), bits));

      // Independent keys: zero, maximum, every single bit, then random keys.
      const auto key =
        static_cast<Index>(iteration == 0                  ? 0
                             : iteration == 1              ? max_key
                             : iteration <= bits * DIM + 1 ? (std::uint64_t {1} << (iteration - 2))
                                                           : random() & max_key);
      const auto decoded = H::toPoint(key, bits);
      ASSERT_TRUE(H::isValidPoint(decoded, bits));
      EXPECT_EQ(H::toIndex(decoded, bits), key);
      EXPECT_EQ(widen(decoded), Reference::toPoint(static_cast<std::uint64_t>(key), bits));
    }

    // Validate the first out-of-range value whenever the input type can hold it.
    if(bits < coord_bits)
    {
      for(int axis = 0; axis < DIM; ++axis)
      {
        typename H::PointType outside;
        outside[axis] = static_cast<Coord>(max_coord + 1);
        EXPECT_FALSE(H::isValidPoint(outside, bits));
      }
    }
    if(bits * DIM < index_bits)
    {
      EXPECT_FALSE(H::isValidIndex(static_cast<Index>(max_key + 1), bits));
    }
  }
  for(int bits :
      {std::numeric_limits<int>::min(), -1, 0, capacity + 1, std::numeric_limits<int>::max()})
  {
    EXPECT_FALSE(H::isValidDepth(bits));
    EXPECT_FALSE(H::isValidPoint(typename H::PointType(), bits));
    EXPECT_FALSE(H::isValidIndex(0, bits));
  }
  if constexpr(std::is_signed<Coord>::value)
  {
    for(int axis = 0; axis < DIM; ++axis)
    {
      typename H::PointType negative;
      negative[axis] = -1;
      EXPECT_FALSE(H::isValidPoint(negative, capacity));
    }
  }
  if constexpr(std::is_signed<Index>::value)
  {
    EXPECT_FALSE(H::isValidIndex(-1, capacity));
  }
}

using IntegerTypes =
  std::tuple<std::int8_t, std::uint8_t, std::int16_t, std::uint16_t, std::int32_t, std::uint32_t, std::int64_t, std::uint64_t>;
constexpr std::size_t NUM_INTEGER_TYPES = std::tuple_size_v<IntegerTypes>;

// Register every coordinate/index pair as its own GoogleTest type parameter.
template <std::size_t... I>
auto integerPairs(std::index_sequence<I...>)
  -> ::testing::Types<std::pair<std::tuple_element_t<I / NUM_INTEGER_TYPES, IntegerTypes>,
                                std::tuple_element_t<I % NUM_INTEGER_TYPES, IntegerTypes>>...>;

using IntegerPairs =
  decltype(integerPairs(std::make_index_sequence<NUM_INTEGER_TYPES * NUM_INTEGER_TYPES> {}));

struct IntegerPairNames
{
  template <typename Pair>
  static std::string GetName(int)
  {
    using Coord = typename Pair::first_type;
    using Index = typename Pair::second_type;
    return std::string("Coord") + (std::is_signed_v<Coord> ? "Int" : "UInt") +
      std::to_string(sizeof(Coord) * 8) + "_Index" + (std::is_signed_v<Index> ? "Int" : "UInt") +
      std::to_string(sizeof(Index) * 8);
  }
};

template <int DIM>
void checkWideCoordinates()
{
  using H = axom::spin::Hilbertizer<std::uint64_t, std::uint8_t, DIM>;
  const std::uint64_t max_coord = (std::uint64_t {1} << H::MAX_UNIQUE_BITS) - 1;
  EXPECT_TRUE(H::isValidPoint(typename H::PointType(max_coord), H::MAX_UNIQUE_BITS));
  for(int axis = 0; axis < DIM; ++axis)
  {
    typename H::PointType point;
    point[axis] = std::uint64_t {1} << 40;
    EXPECT_FALSE(H::isValidPoint(point, H::MAX_UNIQUE_BITS));
    point[axis] |= 1;
    EXPECT_FALSE(H::isValidPoint(point, H::MAX_UNIQUE_BITS));
  }
}

#if AXOM_SLIC_DEBUG_MACROS_ENABLED
template <int DIM>
void checkInvalidPreconditions()
{
  using H = axom::spin::Hilbertizer<int, std::int32_t, DIM>;
  axom::slic::ScopedAbortToThrow abort_guard;
  typename H::PointType negative;
  negative[0] = -1;
  typename H::PointType outside;
  outside[DIM - 1] = 4;
  EXPECT_THROW(H::toIndex(negative, 2), axom::slic::SlicAbortException);
  EXPECT_THROW(H::toIndex(outside, 2), axom::slic::SlicAbortException);
  EXPECT_THROW(H::toPoint(-1, 2), axom::slic::SlicAbortException);
  EXPECT_THROW(H::toPoint(1 << (2 * DIM), 2), axom::slic::SlicAbortException);
  EXPECT_THROW(H::toIndex(typename H::PointType(), 0), axom::slic::SlicAbortException);
  EXPECT_THROW(H::toPoint(0, H::MAX_UNIQUE_BITS + 1), axom::slic::SlicAbortException);
}
#endif

template <typename ExecSpace, int DIM>
void checkExecutionSpace()
{
  using H = axom::spin::Hilbertizer<std::uint32_t, std::uint64_t, DIM>;
  constexpr int count = 128;
  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int allocator = axom::execution_space<ExecSpace>::allocatorID();

  // Sample a depth, a key and a point on the host, cycling through every depth.
  axom::Array<int> depths(count, count, host_allocator);
  axom::Array<std::uint64_t> keys(count, count, host_allocator);
  axom::Array<typename H::PointType> points(count, count, host_allocator);
  std::mt19937_64 random(TEST_SEED);
  for(int i = 0; i < count; ++i)
  {
    depths[i] = 1 + i % H::MAX_UNIQUE_BITS;
    keys[i] = random() & (std::numeric_limits<std::uint64_t>::max() >> (64 - depths[i] * DIM));
    for(int axis = 0; axis < DIM; ++axis)
    {
      points[i][axis] = static_cast<std::uint32_t>(random() & ((std::uint64_t {1} << depths[i]) - 1));
    }
  }

  const axom::Array<int> device_depths(depths, allocator);
  const axom::Array<std::uint64_t> device_keys(keys, allocator);
  const axom::Array<typename H::PointType> device_points(points, allocator);
  axom::Array<std::uint64_t> results(count * 3, count * 3, allocator);
  axom::Array<typename H::PointType> round_trips(count, count, allocator);
  const auto depths_view = device_depths.view();
  const auto keys_view = device_keys.view();
  const auto points_view = device_points.view();
  const auto results_view = results.view();
  const auto round_trips_view = round_trips.view();

  axom::for_all<ExecSpace>(count, [=] AXOM_HOST_DEVICE(axom::IndexType i) {
    const int bits = depths_view[i];
    const auto key = keys_view[i];
    const auto point = H::toPoint(key, bits);
    results_view[3 * i] = H::toIndex(point, bits);
    results_view[3 * i + 1] = axom::spin::convertPointToHilbert<std::uint64_t>(
      axom::spin::convertHilbertToPoint<std::uint32_t, DIM>(key, bits),
      bits);
    results_view[3 * i + 2] = H::isValidPoint(point, bits) && H::isValidIndex(key, bits);
    round_trips_view[i] = H::toPoint(H::toIndex(points_view[i], bits), bits);
  });

  const axom::Array<std::uint64_t> host(results, host_allocator);
  const axom::Array<typename H::PointType> host_round_trips(round_trips, host_allocator);
  for(int i = 0; i < count; ++i)
  {
    EXPECT_EQ(host[3 * i], keys[i]);
    EXPECT_EQ(host[3 * i + 1], keys[i]);
    EXPECT_EQ(host[3 * i + 2], 1);
    EXPECT_EQ(host_round_trips[i], points[i]);
  }
}

template <int DIM>
void checkExecutionSpaces()
{
  checkExecutionSpace<axom::SEQ_EXEC, DIM>();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA)
  checkExecutionSpace<axom::CUDA_EXEC<256>, DIM>();
  #endif

  #if defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP)
  checkExecutionSpace<axom::HIP_EXEC<256>, DIM>();
  #endif
#endif
}
}  // namespace

template <typename Pair>
class Hilbertizer2DIntegerTest : public ::testing::Test
{ };

template <typename Pair>
class Hilbertizer3DIntegerTest : public ::testing::Test
{ };

template <typename Pair>
class Hilbertizer4DIntegerTest : public ::testing::Test
{ };

TYPED_TEST_SUITE(Hilbertizer2DIntegerTest, IntegerPairs, IntegerPairNames);
TYPED_TEST_SUITE(Hilbertizer3DIntegerTest, IntegerPairs, IntegerPairNames);
TYPED_TEST_SUITE(Hilbertizer4DIntegerTest, IntegerPairs, IntegerPairNames);

TYPED_TEST(Hilbertizer2DIntegerTest, every_depth)
{
  allDepths<typename TypeParam::first_type, typename TypeParam::second_type, 2>();
}

TYPED_TEST(Hilbertizer3DIntegerTest, every_depth)
{
  allDepths<typename TypeParam::first_type, typename TypeParam::second_type, 3>();
}

TYPED_TEST(Hilbertizer4DIntegerTest, every_depth)
{
  allDepths<typename TypeParam::first_type, typename TypeParam::second_type, 4>();
}

TEST(Hilbertizer2DTest, known_orientation)
{
  using H2 = axom::spin::Hilbertizer<std::int64_t, std::uint64_t, 2>;
  const H2::PointType square[] = {{0, 0}, {0, 1}, {1, 1}, {1, 0}};
  for(std::uint64_t key = 0; key < 4; ++key)
  {
    EXPECT_EQ(H2::toIndex(square[key], 1), key);
    EXPECT_EQ(H2::toPoint(key, 1), square[key]);
  }

  const H2::PointType curve[] = {{0, 0},
                                 {1, 0},
                                 {1, 1},
                                 {0, 1},
                                 {0, 2},
                                 {0, 3},
                                 {1, 3},
                                 {1, 2},
                                 {2, 2},
                                 {2, 3},
                                 {3, 3},
                                 {3, 2},
                                 {3, 1},
                                 {2, 1},
                                 {2, 0},
                                 {3, 0}};
  for(std::uint64_t key = 0; key < 16; ++key)
  {
    EXPECT_EQ(H2::toIndex(curve[key], 2), key);
    EXPECT_EQ(H2::toPoint(key, 2), curve[key]);
  }
}

TEST(Hilbertizer3DTest, known_orientation)
{
  using H3 = axom::spin::Hilbertizer<std::int64_t, std::int64_t, 3>;
  EXPECT_EQ(H3::toIndex({5, 10, 20}, 5), 7865);
  EXPECT_EQ(H3::toPoint(7865, 5), (H3::PointType {5, 10, 20}));
  const H3::PointType cube[] =
    {{0, 0, 0}, {0, 0, 1}, {0, 1, 1}, {0, 1, 0}, {1, 1, 0}, {1, 1, 1}, {1, 0, 1}, {1, 0, 0}};
  for(std::int64_t key = 0; key < 8; ++key)
  {
    EXPECT_EQ(H3::toIndex(cube[key], 1), key);
    EXPECT_EQ(H3::toPoint(key, 1), cube[key]);
  }
}

TEST(Hilbertizer4DTest, known_orientation)
{
  using H = axom::spin::Hilbertizer<std::uint32_t, std::uint64_t, 4>;
  // At depth one, use binary-reflected Gray codes with axis zero first.
  for(std::uint64_t key = 0; key < 16; ++key)
  {
    const auto gray = key ^ (key >> 1);
    H::PointType point;
    for(int axis = 0; axis < 4; ++axis)
    {
      point[axis] = static_cast<std::uint32_t>((gray >> (3 - axis)) & 1);
    }
    EXPECT_EQ(H::toIndex(point, 1), key);
    EXPECT_EQ(H::toPoint(key, 1), point);
  }
}

TEST(HilbertizerTest, max_depth_reference_keys)
{
  // Test roundtrip on unsigned types: point -> Hilbert index -> point
  {
    using H2 = axom::spin::Hilbertizer<std::uint64_t, std::uint64_t, 2>;
    using H3 = axom::spin::Hilbertizer<std::uint64_t, std::uint64_t, 3>;
    using H4 = axom::spin::Hilbertizer<std::uint64_t, std::uint64_t, 4>;
    const H2::PointType p2 {3141592653ULL, 2718281828ULL};
    const H3::PointType p3 {1234567, 765432, 2097151};
    const H4::PointType p4 {12345, 54321, 65535, 1};

    EXPECT_EQ(H2::toIndex(p2, 32), 9893194776509318881ULL);
    EXPECT_EQ(H2::toPoint(9893194776509318881ULL, 32), p2);
    EXPECT_EQ(H3::toIndex(p3, 21), 7299960257221567926ULL);
    EXPECT_EQ(H3::toPoint(7299960257221567926ULL, 21), p3);
    EXPECT_EQ(H4::toIndex(p4, 16), 4954868739613205852ULL);
    EXPECT_EQ(H4::toPoint(4954868739613205852ULL, 16), p4);
  }

  // Test roundtrip on signed types: point -> Hilbert index -> point
  {
    using S2 = axom::spin::Hilbertizer<std::int64_t, std::int64_t, 2>;
    using S4 = axom::spin::Hilbertizer<std::int64_t, std::int64_t, 4>;
    const S2::PointType s2 {1234567890, 987654321};
    const S4::PointType s4 {12345, 4321, 32767, 1};

    EXPECT_EQ(S2::toIndex(s2, 31), 3831921454093212845LL);
    EXPECT_EQ(S2::toPoint(3831921454093212845LL, 31), s2);
    EXPECT_EQ(S4::toIndex(s4, 15), 264373023102641500LL);
    EXPECT_EQ(S4::toPoint(264373023102641500LL, 15), s4);
  }
}

TEST(Hilbertizer2DTest, exhaustive_bijection_adjacency_and_prefixes)
{
  for(int bits = 1; bits <= 5; ++bits)
  {
    SCOPED_TRACE(bits);
    ASSERT_NO_FATAL_FAILURE(exhaustive<2>(bits));
  }
}

TEST(Hilbertizer3DTest, exhaustive_bijection_adjacency_and_prefixes)
{
  for(int bits = 1; bits <= 5; ++bits)
  {
    SCOPED_TRACE(bits);
    ASSERT_NO_FATAL_FAILURE(exhaustive<3>(bits));
  }
}

TEST(Hilbertizer4DTest, exhaustive_bijection_adjacency_and_prefixes)
{
  for(int bits = 1; bits <= 4; ++bits)
  {
    SCOPED_TRACE(bits);
    ASSERT_NO_FATAL_FAILURE(exhaustive<4>(bits));
  }
}

TEST(Hilbertizer2DTest, wide_coordinates_are_validated_before_narrowing)
{
  checkWideCoordinates<2>();
}

TEST(Hilbertizer3DTest, wide_coordinates_are_validated_before_narrowing)
{
  checkWideCoordinates<3>();
}

TEST(Hilbertizer4DTest, wide_coordinates_are_validated_before_narrowing)
{
  checkWideCoordinates<4>();
}

#if AXOM_SLIC_DEBUG_MACROS_ENABLED
TEST(Hilbertizer2DTest, invalid_preconditions) { checkInvalidPreconditions<2>(); }

TEST(Hilbertizer3DTest, invalid_preconditions) { checkInvalidPreconditions<3>(); }

TEST(Hilbertizer4DTest, invalid_preconditions) { checkInvalidPreconditions<4>(); }
#endif

AXOM_CUDA_TEST(Hilbertizer2DTest, execution_spaces) { checkExecutionSpaces<2>(); }

AXOM_CUDA_TEST(Hilbertizer3DTest, execution_spaces) { checkExecutionSpaces<3>(); }

AXOM_CUDA_TEST(Hilbertizer4DTest, execution_spaces) { checkExecutionSpaces<4>(); }

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;
  return RUN_ALL_TESTS();
}
