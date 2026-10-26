// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file spin_hilbert_fsm.cpp
 * \brief Tests Hilbert FSM conversion against the scalar transform.
 */

#include "axom/spin/internal/HilbertFSM.hpp"
#include "axom/spin/internal/HilbertIndexImpl.hpp"
#include "axom/spin/HilbertIndex.hpp"
#include "axom/slic.hpp"

#include "gtest/gtest.h"

#include <cstdint>
#include <limits>
#include <random>

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
 * along with independently sampled keys.
 *
 * \tparam Index Integer type for Hilbert indices.
 * \tparam DIM Number of coordinates in each point.
 * \tparam LEVELS Number of bits per coordinate processed by each table lookup.
 */
template <typename Index, int DIM, int LEVELS>
void allDepths()
{
  using FSM = hilbert::HilbertFSMTransform<std::uint32_t, Index, DIM, LEVELS>;
  using Scalar = hilbert::HilbertTransform<std::uint32_t, Index, DIM>;
  using Public = axom::spin::Hilbertizer<std::uint32_t, Index, DIM>;
  std::mt19937_64 random(TEST_SEED + DIM);
  for(int bits = 1; bits <= Public::MAX_UNIQUE_BITS; ++bits)
  {
    SCOPED_TRACE(::testing::Message()
                 << "dimension " << DIM << ", group " << LEVELS << ", depth " << bits);
    const std::uint64_t coord_mask = (std::uint64_t {1} << bits) - 1;
    const auto key_mask = std::numeric_limits<std::uint64_t>::max() >> (64 - DIM * bits);
    const auto check_point = [&](const typename FSM::PointType& point) {
      const auto expected = Scalar::toIndex(point, bits);
      EXPECT_EQ(FSM::toIndex(point, bits), expected);
      EXPECT_EQ(FSM::toPoint(expected, bits), point);
    };

    // Single coordinate bits exercise every input position, including the high bit.
    for(int axis = 0; axis < DIM; ++axis)
    {
      for(int bit = 0; bit < bits; ++bit)
      {
        typename FSM::PointType point;
        point[axis] = static_cast<std::uint32_t>(std::uint64_t {1} << bit);
        check_point(point);
      }
    }
    for(int sample = 0; sample < SAMPLES_PER_DEPTH; ++sample)
    {
      typename FSM::PointType point;
      for(int axis = 0; axis < DIM; ++axis)
      {
        point[axis] = static_cast<std::uint32_t>(
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
    }
  }
}

/*!
 * \brief Compares FSM and scalar encoding and decoding for every point and key on small grids.
 *
 * \tparam DIM Number of coordinates in each point.
 * \tparam LEVELS Number of bits per coordinate processed by each table lookup.
 */
template <int DIM, int LEVELS>
void exhaustive()
{
  using FSM = hilbert::HilbertFSMTransform<std::uint32_t, std::uint64_t, DIM, LEVELS>;
  using Scalar = hilbert::HilbertTransform<std::uint32_t, std::uint64_t, DIM>;
  constexpr int max_depth = DIM == 4 ? 4 : 5;
  for(int bits = 1; bits <= max_depth; ++bits)
  {
    SCOPED_TRACE(::testing::Message()
                 << "dimension " << DIM << ", group " << LEVELS << ", depth " << bits);
    const std::uint64_t side = std::uint64_t {1} << bits;
    const std::uint64_t count = std::uint64_t {1} << (DIM * bits);
    for(std::uint64_t key = 0; key < count; ++key)
    {
      // Enumerate coordinates independently of either Hilbert decoder.
      typename FSM::PointType point;
      auto remaining = key;
      for(int axis = 0; axis < DIM; ++axis)
      {
        point[axis] = static_cast<std::uint32_t>(remaining % side);
        remaining /= side;
      }
      EXPECT_EQ(FSM::toIndex(point, bits), Scalar::toIndex(point, bits));
      EXPECT_EQ(FSM::toPoint(key, bits), Scalar::toPoint(key, bits));
    }
  }
}

// Public API tests cover signed and narrow types, orientation, adjacency, and prefixes.
// Here, exercise every lookup group with both working widths.
using IndexTypes = ::testing::Types<std::uint32_t, std::uint64_t>;
}  // namespace

template <typename Index>
class HilbertFSMIndexTest : public ::testing::Test
{ };

TYPED_TEST_SUITE(HilbertFSMIndexTest, IndexTypes);

TYPED_TEST(HilbertFSMIndexTest, every_depth_and_group_size)
{
  using Index = TypeParam;
  allDepths<Index, 2, 1>();
  allDepths<Index, 2, 2>();
  allDepths<Index, 2, 3>();
  allDepths<Index, 3, 1>();
  allDepths<Index, 3, 2>();
  allDepths<Index, 4, 1>();
}

TEST(HilbertFSMTest, exhaustive_scalar_agreement)
{
  ASSERT_NO_FATAL_FAILURE((exhaustive<2, 1>()));
  ASSERT_NO_FATAL_FAILURE((exhaustive<2, 2>()));
  ASSERT_NO_FATAL_FAILURE((exhaustive<2, 3>()));
  ASSERT_NO_FATAL_FAILURE((exhaustive<3, 1>()));
  ASSERT_NO_FATAL_FAILURE((exhaustive<3, 2>()));
  ASSERT_NO_FATAL_FAILURE((exhaustive<4, 1>()));
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;
  return RUN_ALL_TESTS();
}
