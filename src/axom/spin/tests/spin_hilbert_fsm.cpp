// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file spin_hilbert_fsm.cpp
 * \brief Tests Hilbert FSM conversion, orientation, adjacency, and prefixes.
 */

#include "axom/spin/internal/HilbertFSM.hpp"
#include "axom/spin/internal/HilbertIndexImpl.hpp"
#include "axom/spin/HilbertIndex.hpp"
#include "axom/slic.hpp"

#include "gtest/gtest.h"

#include <cstdint>
#include <cstdlib>
#include <limits>
#include <random>
#include <utility>
#include <vector>

namespace
{
namespace hilbert = axom::spin::internal::hilbert;
constexpr std::uint64_t TEST_SEED = 42;
constexpr int SAMPLES_PER_DEPTH = 100;

/*!
 * \brief Compares finite state machine (FSM) encoding and decoding against
 * the scalar transform at every valid depth.
 *
 * Covers corners, individual coordinate bits, and deterministic random points,
 * along with independently sampled keys and round trips in both directions.
 *
 * \tparam Coord Integer type for point coordinates.
 * \tparam Index Integer type for Hilbert indices.
 * \tparam DIM Number of coordinates in each point.
 * \tparam LEVELS Number of bits per coordinate processed by each table lookup.
 */
template <typename Coord, typename Index, int DIM, int LEVELS>
void allDepths()
{
  using FSM = hilbert::HilbertFSMTransform<Coord, Index, DIM, LEVELS>;
  using Scalar = hilbert::HilbertTransform<Coord, Index, DIM>;
  using Public = axom::spin::Hilbertizer<Coord, Index, DIM>;
  using Work = typename FSM::WorkType;
  using Morton = typename FSM::MortonizerType;
  std::mt19937_64 random(TEST_SEED + DIM);
  for(int bits = 1; bits <= Public::MAX_UNIQUE_BITS; ++bits)
  {
    SCOPED_TRACE(::testing::Message()
                 << "dimension " << DIM << ", group " << LEVELS << ", depth " << bits);
    const std::uint64_t coord_mask = (std::uint64_t {1} << bits) - 1;
    const auto key_mask = std::numeric_limits<std::uint64_t>::max() >> (64 - DIM * bits);
    const auto check_point = [&](const typename FSM::PointType& point) {
      const auto expected = Scalar::toIndex(point, bits);
      const auto morton = Morton::mortonize(point);
      EXPECT_EQ(FSM::toIndex(point, bits), expected);
      EXPECT_EQ(FSM::toPoint(expected, bits), point);
      EXPECT_EQ((hilbert::convert<DIM, LEVELS, false>(morton, bits)), static_cast<Work>(expected));
      EXPECT_EQ((hilbert::convert<DIM, LEVELS, true>(static_cast<Work>(expected), bits)), morton);
    };

    // Single coordinate bits exercise every input position, including the high bit.
    for(int axis = 0; axis < DIM; ++axis)
    {
      for(int bit = 0; bit < bits; ++bit)
      {
        typename FSM::PointType point;
        point[axis] = static_cast<Coord>(std::uint64_t {1} << bit);
        check_point(point);
      }
    }
    for(int sample = 0; sample < SAMPLES_PER_DEPTH; ++sample)
    {
      typename FSM::PointType point;
      for(int axis = 0; axis < DIM; ++axis)
      {
        point[axis] = static_cast<Coord>(
          sample < (1 << DIM) ? (((sample >> axis) & 1) ? coord_mask : 0) : random() & coord_mask);
      }
      check_point(point);

      // Decode independent keys: zero, maximum, single bits, and random samples.
      const auto key =
        static_cast<Index>(sample == 0                  ? 0
                             : sample == 1              ? key_mask
                             : sample <= DIM * bits + 1 ? std::uint64_t {1} << (sample - 2)
                                                        : random() & key_mask);
      const auto decoded = FSM::toPoint(key, bits);
      EXPECT_EQ(decoded, Scalar::toPoint(key, bits));
      EXPECT_TRUE(Public::isValidPoint(decoded, bits));
      EXPECT_EQ(FSM::toIndex(decoded, bits), key);
    }
  }
}

/*!
 * \brief Exhaustively checks conversion, adjacency, and prefixes on small grids.
 *
 * Enumerates every point and key to check scalar agreement, unique keys, and round trips.
 * Consecutive keys must map to neighboring grid points.
 * Removing the lowest bit from each coordinate must remove the lowest DIM bits from its key.
 *
 * \tparam DIM Number of coordinates in each point.
 * \tparam LEVELS Number of bits per coordinate processed by each table lookup.
 */
template <int DIM, int LEVELS>
void exhaustive()
{
  using FSM = hilbert::HilbertFSMTransform<std::int64_t, std::uint64_t, DIM, LEVELS>;
  using Scalar = hilbert::HilbertTransform<std::int64_t, std::uint64_t, DIM>;
  constexpr int max_depth = DIM == 4 ? 4 : 5;
  for(int bits = 1; bits <= max_depth; ++bits)
  {
    SCOPED_TRACE(::testing::Message()
                 << "dimension " << DIM << ", group " << LEVELS << ", depth " << bits);
    const std::uint64_t side = std::uint64_t {1} << bits;
    const std::uint64_t count = std::uint64_t {1} << (DIM * bits);
    std::vector<bool> seen(count, false);
    for(std::uint64_t ordinal = 0; ordinal < count; ++ordinal)
    {
      typename FSM::PointType point;
      auto remaining = ordinal;
      for(int axis = 0; axis < DIM; ++axis)
      {
        point[axis] = remaining % side;
        remaining /= side;
      }
      const auto key = FSM::toIndex(point, bits);
      ASSERT_LT(key, count);
      ASSERT_FALSE(seen[key]);
      seen[key] = true;
      EXPECT_EQ(key, Scalar::toIndex(point, bits));
      EXPECT_EQ(FSM::toPoint(key, bits), point);
    }
    typename FSM::PointType previous;
    for(std::uint64_t key = 0; key < count; ++key)
    {
      const auto point = FSM::toPoint(key, bits);
      EXPECT_EQ(point, Scalar::toPoint(key, bits));
      EXPECT_EQ(FSM::toIndex(point, bits), key);
      if(key > 0)
      {
        std::int64_t distance = 0;
        for(int axis = 0; axis < DIM; ++axis)
        {
          distance += std::abs(point[axis] - previous[axis]);
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
        EXPECT_EQ(FSM::toIndex(parent, bits - 1), key >> DIM);
      }
      previous = point;
    }
  }
}

/*!
 * \brief Checks the Hilbert curve's orientation against known point/key pairs.
 *
 * Checks both conversion directions using depth-one Gray-code fixtures,
 * a fixed 2D traversal at depth two, and selected pairs at greater depths.
 *
 * \tparam DIM Number of coordinates in each point.
 * \tparam LEVELS Number of bits per coordinate processed by each table lookup.
 */
template <int DIM, int LEVELS>
void knownOrientation()
{
  using FSM = hilbert::HilbertFSMTransform<std::uint64_t, std::uint64_t, DIM, LEVELS>;
  const auto check = [](const typename FSM::PointType& point, std::uint64_t key, int bits) {
    EXPECT_EQ(FSM::toIndex(point, bits), key);
    EXPECT_EQ(FSM::toPoint(key, bits), point);
  };
  // Depth-one fixtures follow the binary-reflected Gray code, with axis zero first.
  for(std::uint64_t key = 0; key < (1u << DIM); ++key)
  {
    const auto gray = key ^ (key >> 1);
    typename FSM::PointType point;
    for(int axis = 0; axis < DIM; ++axis)
    {
      point[axis] = (gray >> (DIM - 1 - axis)) & 1;
    }
    check(point, key, 1);
  }
  if constexpr(DIM == 2)
  {
    const typename FSM::PointType curve[] = {{0, 0},
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
      check(curve[key], key, 2);
    }
    check({3141592653ULL, 2718281828ULL}, 9893194776509318881ULL, 32);
  }
  else if constexpr(DIM == 3)
  {
    check({5, 10, 20}, 7865, 5);
    check({1234567, 765432, 2097151}, 7299960257221567926ULL, 21);
  }
  else
  {
    check({12345, 54321, 65535, 1}, 4954868739613205852ULL, 16);
  }
}

using IntegerPairs = ::testing::Types<std::pair<std::int8_t, std::int8_t>,
                                      std::pair<std::uint8_t, std::uint8_t>,
                                      std::pair<std::int16_t, std::int16_t>,
                                      std::pair<std::uint16_t, std::uint16_t>,
                                      std::pair<std::int32_t, std::int32_t>,
                                      std::pair<std::uint32_t, std::uint32_t>,
                                      std::pair<std::int64_t, std::int64_t>,
                                      std::pair<std::uint64_t, std::uint64_t>,
                                      std::pair<std::int8_t, std::uint64_t>,
                                      std::pair<std::uint32_t, std::uint64_t>,
                                      std::pair<std::int64_t, std::uint32_t>,
                                      std::pair<std::uint32_t, std::int64_t>>;
}  // namespace

template <typename Pair>
class HilbertFSMIntegerTest : public ::testing::Test
{ };

TYPED_TEST_SUITE(HilbertFSMIntegerTest, IntegerPairs);

TYPED_TEST(HilbertFSMIntegerTest, every_depth_and_group_size)
{
  using Coord = typename TypeParam::first_type;
  using Index = typename TypeParam::second_type;
  allDepths<Coord, Index, 2, 1>();
  allDepths<Coord, Index, 2, 2>();
  allDepths<Coord, Index, 2, 3>();
  allDepths<Coord, Index, 3, 1>();
  allDepths<Coord, Index, 3, 2>();
  allDepths<Coord, Index, 4, 1>();
}

TEST(HilbertFSMTest, exhaustive_bijection_adjacency_and_prefixes)
{
  ASSERT_NO_FATAL_FAILURE((exhaustive<2, 1>()));
  ASSERT_NO_FATAL_FAILURE((exhaustive<2, 2>()));
  ASSERT_NO_FATAL_FAILURE((exhaustive<2, 3>()));
  ASSERT_NO_FATAL_FAILURE((exhaustive<3, 1>()));
  ASSERT_NO_FATAL_FAILURE((exhaustive<3, 2>()));
  ASSERT_NO_FATAL_FAILURE((exhaustive<4, 1>()));
}

TEST(HilbertFSMTest, known_orientation)
{
  knownOrientation<2, 1>();
  knownOrientation<2, 2>();
  knownOrientation<2, 3>();
  knownOrientation<3, 1>();
  knownOrientation<3, 2>();
  knownOrientation<4, 1>();
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;
  return RUN_ALL_TESTS();
}
