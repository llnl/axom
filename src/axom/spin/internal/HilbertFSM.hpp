// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/*!
 * \file HilbertFSM.hpp
 * \brief Host Hilbert encoding and decoding with grouped finite state machine (FSM) lookups.
 */

#include "axom/spin/MortonIndex.hpp"
#include "axom/spin/internal/HilbertTables.hpp"

#include <cstdint>
#include <type_traits>

namespace axom::spin::internal::hilbert
{
/// Returns converted bits and a state update for one group.
/// The row is a 2D orientation state or a 3D/4D axis permutation.
template <int DIM, int LEVELS, bool INVERSE>
inline unsigned int lookup(unsigned int state, unsigned int digits)
{
  static_assert(DIM >= 2 && DIM <= 4, "Supported dimensions are 2, 3, and 4");
  static_assert(LEVELS >= 1 && LEVELS <= 5 - DIM,
                "Supported groups are 1..3 in 2D, 1..2 in 3D, and 1 in 4D");
  const unsigned int offset = (state << (DIM * LEVELS)) | digits;
  if constexpr(DIM == 2)
  {
    if constexpr(LEVELS == 1)
    {
      return INVERSE ? inverseTable2D1[offset] : table2D1[offset];
    }
    else if constexpr(LEVELS == 2)
    {
      return INVERSE ? inverseTable2D2[offset] : table2D2[offset];
    }
    else
    {
      // Load all four state bytes using only the input digits, then select the current state's byte.
      // This keeps the load address independent of the preceding lookup's state update.
      const auto entry = INVERSE ? inverseTable2D3[digits] : table2D3[digits];
      return (entry >> (state * 8)) & 0xffu;
    }
  }
  else if constexpr(DIM == 3)
  {
    if constexpr(LEVELS == 1)
    {
      return INVERSE ? inverseTable3D1[offset] : table3D1[offset];
    }
    else
    {
      return INVERSE ? inverseTable3D2[offset] : table3D2[offset];
    }
  }
  else
  {
    return INVERSE ? inverseTable4D1[offset] : table4D1[offset];
  }
}

/*!
 * \brief Converts Morton to Hilbert keys, or Hilbert to Morton when INVERSE is true.
 *
 * Each full lookup processes LEVELS bits per coordinate, starting with the most significant bits.
 * A partial leading group uses a smaller table. Padding it with zeros changes the curve's orientation.
 * In 2D, an entry holds the output bits in its low bits and the next state in its high bits.
 * In 3D/4D, the state holds the axis permutation and a reflection mask repeated for each input plane.
 * An entry's low bits hold the output; its high bits hold an XOR mask for updating the state.
 *
 * \pre bits is positive and DIM * bits does not exceed the bit width of U.
 * \pre source uses only the lowest \c DIM*bits bits.
 */
template <int DIM, int LEVELS, bool INVERSE, typename U>
inline U convert(U source, int bits)
{
  static_assert(std::is_unsigned<U>::value, "Key conversion requires unsigned arithmetic");
  constexpr int group_bits = DIM * LEVELS;
  constexpr unsigned int mask = (1u << group_bits) - 1;
  unsigned int state = 0;
  U result = 0;
  if constexpr(DIM >= 3)
  {
    if constexpr(LEVELS == 2)
    {
      if(bits & 1)
      {
        --bits;
        const auto entry = lookup<DIM, 1, INVERSE>(0, static_cast<unsigned>(source >> (bits * DIM)));
        result = entry & ((1u << DIM) - 1);
        const auto next = entry >> DIM;
        // Repeat the reflection mask for both planes and move the permutation bits above it.
        state =
          ((next >> DIM) << group_bits) | ((next & ((1u << DIM) - 1)) * (mask / ((1u << DIM) - 1)));
      }
    }
    while(bits > 0)
    {
      bits -= LEVELS;
      const auto digit = static_cast<unsigned>((source >> (bits * DIM)) & mask);
      const auto offset = INVERSE ? ((state & ~mask) | digit) : (state ^ digit);
      const auto entry = lookup<DIM, LEVELS, INVERSE>(offset >> group_bits, offset & mask);
      result = (result << group_bits) | ((INVERSE ? (entry ^ state) : entry) & mask);
      state ^= entry >> group_bits;
    }
    return result;
  }
  else
  {
    if constexpr(LEVELS > 1)
    {
      if(bits % LEVELS == 1)
      {
        bits -= 1;
        const auto entry =
          lookup<DIM, 1, INVERSE>(0, static_cast<unsigned int>(source >> (bits * DIM)));
        result = static_cast<U>(entry & ((1u << DIM) - 1));
        state = entry >> DIM;
      }
    }
    if constexpr(LEVELS > 2)
    {
      if(bits % LEVELS == 2)
      {
        bits -= 2;
        const auto entry =
          lookup<DIM, 2, INVERSE>(0, static_cast<unsigned int>(source >> (bits * DIM)));
        result = static_cast<U>(entry & ((1u << (2 * DIM)) - 1));
        state = entry >> (2 * DIM);
      }
    }
    while(bits > 0)
    {
      bits -= LEVELS;
      const auto entry =
        lookup<DIM, LEVELS, INVERSE>(state,
                                     static_cast<unsigned int>((source >> (bits * DIM)) & mask));
      result = static_cast<U>((result << group_bits) | (entry & mask));
      state = entry >> group_bits;
    }
    return result;
  }
}

/*!
 * \brief Converts integer points and Hilbert indices through Morton packing and FSM tables.
 *
 * \tparam LEVELS Bits per coordinate processed by each full lookup group.
 * \note Host only. Callers must satisfy Hilbertizer's point, index, and depth preconditions.
 *       These routines do not validate inputs.
 */
template <typename CoordType, typename IndexType, int DIM, int LEVELS>
struct HilbertFSMTransform
{
  using PointType = primal::Point<CoordType, DIM>;
  // Unsigned words of at least 32 bits avoid signed shifts and narrow integer promotion.
  using WorkType = std::conditional_t<(sizeof(IndexType) <= 4), std::uint32_t, std::uint64_t>;
  using MortonizerType = Mortonizer<CoordType, WorkType, DIM>;

  static IndexType toIndex(const PointType& point, int bits)
  {
    return static_cast<IndexType>(convert<DIM, LEVELS, false>(MortonizerType::mortonize(point), bits));
  }

  static PointType toPoint(IndexType index, int bits)
  {
    return MortonizerType::demortonize(convert<DIM, LEVELS, true>(static_cast<WorkType>(index), bits));
  }
};
}  // namespace axom::spin::internal::hilbert
