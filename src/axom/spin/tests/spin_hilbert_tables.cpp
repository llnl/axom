// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file spin_hilbert_tables.cpp
 * \brief Checks generated Hilbert table outputs, inverses, and state transitions.
 */

#include "axom/spin/internal/HilbertTables.hpp"
#include "axom/spin/internal/HilbertIndexImpl.hpp"
#include "axom/slic.hpp"

#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>

namespace
{
namespace tables = axom::spin::internal::hilbert;

// Unpack Morton digits and apply the incoming permutation and reflections before scalar encoding.
template <int DIM>
std::uint64_t referenceIndex(std::uint64_t digits,
                             int bits,
                             const std::array<int, DIM>& permutation,
                             unsigned flips)
{
  using Scalar = tables::HilbertTransform<std::uint32_t, std::uint64_t, DIM>;
  typename Scalar::PointType point;
  for(int axis = 0; axis < DIM; ++axis)
  {
    for(int bit = 0; bit < bits; ++bit)
    {
      const auto value =
        ((digits >> (bit * DIM + permutation[axis])) ^ (flips >> permutation[axis])) & 1;
      point[axis] |= static_cast<std::uint32_t>(value) << bit;
    }
  }
  return Scalar::toIndex(point, bits);
}

template <int LEVELS, typename Entry, std::size_t N>
void check2DTable(const Entry (&forward)[N], const Entry (&inverse)[N])
{
  constexpr unsigned dim = 2;
  constexpr unsigned states = 4;
  constexpr unsigned plane_mask = (1u << dim) - 1;
  constexpr unsigned width = dim * LEVELS;
  constexpr unsigned size = 1u << width;
  constexpr unsigned mask = size - 1;
  static_assert(N == (LEVELS == 3 ? size : states * size),
                "2D tables store four states for each input group");
  const auto lookup = [](const Entry(&table)[N], unsigned state, unsigned digit) {
    if constexpr(LEVELS == 3)
    {
      // Three-plane tables pack the four state bytes into one word per digit.
      return (table[digit] >> (state * 8)) & 0xffu;
    }
    else
    {
      return static_cast<unsigned>(table[state * size + digit]);
    }
  };

  // Find a short input prefix reaching each state. Compare each prefix and group
  // with the scalar transform to check every row's orientation.
  std::array<std::uint64_t, states> prefixes {};
  std::array<int, states> depths {};
  std::array<bool, states> reached {};
  std::array<unsigned, states> pending {};
  reached[0] = true;
  unsigned count = 1;
  for(unsigned head = 0; head < count; ++head)
  {
    const auto state = pending[head];
    for(unsigned digit = 0; digit <= plane_mask; ++digit)
    {
      const auto next = tables::table2D1[(state << dim) | digit] >> dim;
      ASSERT_LT(next, states);
      if(!reached[next])
      {
        reached[next] = true;
        prefixes[next] = (prefixes[state] << dim) | digit;
        depths[next] = depths[state] + 1;
        pending[count++] = next;
      }
    }
  }
  ASSERT_EQ(count, states);

  for(unsigned state = 0; state < states; ++state)
  {
    for(unsigned digits = 0; digits < size; ++digits)
    {
      SCOPED_TRACE(::testing::Message()
                   << "planes " << LEVELS << ", state " << state << ", digits " << digits);
      const unsigned entry = lookup(forward, state, digits);
      const auto output = entry & mask;
      const auto next = entry >> width;
      ASSERT_LT(next, states);
      const auto decoded = lookup(inverse, state, output);
      EXPECT_EQ(decoded & mask, digits);
      EXPECT_EQ(decoded >> width, next);

      const auto source = (prefixes[state] << width) | digits;
      const int bits = depths[state] + LEVELS;
      EXPECT_EQ(output, referenceIndex<2>(source, bits, {0, 1}, 0) & mask);
      // Every possible following digit tests the orientation left by this group.
      for(unsigned tail = 0; tail <= plane_mask; ++tail)
      {
        const auto suffix = tables::table2D1[(next << dim) | tail] & plane_mask;
        const auto expected = referenceIndex<2>((source << dim) | tail, bits + 1, {0, 1}, 0);
        EXPECT_EQ((output << dim) | suffix, expected & ((size << dim) - 1));
      }
    }
  }
}

template <int DIM, int LEVELS, std::size_t N, std::size_t SINGLE_N>
void checkFactoredTable(const std::uint16_t (&forward)[N],
                        const std::uint16_t (&inverse)[N],
                        const std::uint16_t (&single)[SINGLE_N])
{
  // DIM! lexicographically ordered axis permutations select the stored rows.
  constexpr unsigned permutations = DIM == 3 ? 6 : 24;
  constexpr unsigned width = DIM * LEVELS;
  constexpr unsigned size = 1u << width;
  constexpr unsigned mask = size - 1;
  constexpr unsigned plane_mask = (1u << DIM) - 1;
  constexpr unsigned repeat = mask / plane_mask;
  static_assert(N == permutations * size, "Factored tables have DIM! rows");
  static_assert(SINGLE_N == permutations * (1u << DIM), "One-plane table size");

  std::array<int, DIM> permutation;
  for(int axis = 0; axis < DIM; ++axis)
  {
    permutation[axis] = axis;
  }
  unsigned row = 0;
  do
  {
    for(unsigned flips = 0; flips <= plane_mask; ++flips)
    {
      const auto state = (row << width) | flips * repeat;
      for(unsigned digits = 0; digits < size; ++digits)
      {
        SCOPED_TRACE(::testing::Message() << "dimension " << DIM << ", planes " << LEVELS << ", row "
                                          << row << ", flips " << flips << ", digits " << digits);
        const unsigned entry = forward[state ^ digits];
        const auto output = entry & mask;
        const auto next = state ^ (entry >> width);
        const auto next_permutation = next >> width;
        const auto next_flips = next & plane_mask;
        ASSERT_LT(next_permutation, permutations);
        EXPECT_EQ(next & mask, next_flips * repeat);

        const unsigned decoded = inverse[(row << width) | output];
        EXPECT_EQ((decoded ^ state) & mask, digits);
        EXPECT_EQ(state ^ (decoded >> width), next);

        EXPECT_EQ(output, referenceIndex<DIM>(digits, LEVELS, permutation, flips));

        // Append every possible digit and compare with the scalar transform
        // to check the next orientation.
        for(unsigned tail = 0; tail <= plane_mask; ++tail)
        {
          const auto suffix = single[(next_permutation << DIM) | (tail ^ next_flips)] & plane_mask;
          const auto expected =
            referenceIndex<DIM>((std::uint64_t {digits} << DIM) | tail, LEVELS + 1, permutation, flips);
          EXPECT_EQ((output << DIM) | suffix, expected);
        }
      }
    }
    ++row;
  } while(std::next_permutation(permutation.begin(), permutation.end()));
  EXPECT_EQ(row, permutations);
}
}  // namespace

TEST(HilbertTablesTest, two_dimensional_rows_and_packed_states)
{
  ASSERT_NO_FATAL_FAILURE(check2DTable<1>(tables::table2D1, tables::inverseTable2D1));
  ASSERT_NO_FATAL_FAILURE(check2DTable<2>(tables::table2D2, tables::inverseTable2D2));
  ASSERT_NO_FATAL_FAILURE(check2DTable<3>(tables::table2D3, tables::inverseTable2D3));
}

TEST(HilbertTablesTest, three_dimensional_rows_and_reflections)
{
  ASSERT_NO_FATAL_FAILURE(
    (checkFactoredTable<3, 1>(tables::table3D1, tables::inverseTable3D1, tables::table3D1)));
  ASSERT_NO_FATAL_FAILURE(
    (checkFactoredTable<3, 2>(tables::table3D2, tables::inverseTable3D2, tables::table3D1)));
}

TEST(HilbertTablesTest, four_dimensional_rows_and_reflections)
{
  ASSERT_NO_FATAL_FAILURE(
    (checkFactoredTable<4, 1>(tables::table4D1, tables::inverseTable4D1, tables::table4D1)));
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;
  return RUN_ALL_TESTS();
}
