// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/*! 
 * \file hilbert_reference.hpp
 * \brief Provides a scalar Hilbert reference for tests and benchmarks.
 */

#include "axom/core/Macros.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/spin/MortonIndex.hpp"

#include <cstdint>
#include <type_traits>

namespace axom::spin::testing
{
/*!
 * \brief Converts between integer grid points and Hilbert indices.
 *
 * A Hilbert index gives a point's position on a path through every grid point.
 * Successive indices visit points one grid step apart. Skilling's algorithm
 * computes this order with XORs and with flips and exchanges of lower coordinate bits.
 * These operations account for the path's orientation within each region.
 * Encoding applies the orientation operations before the inverse Gray-code transform.
 * Decoding applies the Gray-code transform first, then restores the coordinates.
 * A Gray code orders bit patterns so that successive patterns differ in one bit.
 *
 * \note Skilling's AxestoTranspose and TransposetoAxes store the index bits in
 * DIM unsigned integers, one per axis. This layout is called the transpose.
 * For example, three two-bit words a1 a0, b1 b0, and c1 c0 represent the single
 * index a1 b1 c1 a0 b0 c0, written from most to least significant bit.
 * This class also packs and unpacks that single integer using Morton interleaving.
 */
template <typename CoordType, typename IndexType, int DIM>
struct HilbertReference
{
  using PointType = primal::Point<CoordType, DIM>;
  // Use unsigned words of at least 32 bits to avoid integer promotion in shifts.
  using WorkType =
    typename std::conditional_t<(sizeof(IndexType) <= 4), std::uint32_t, std::uint64_t>;
  using WorkPoint = primal::Point<WorkType, DIM>;
  using MortonizerType = Mortonizer<WorkType, WorkType, DIM>;

  /*!
   * \brief Encodes an integer lattice point as a Hilbert index.
   *
   * \param[in] point Nonnegative integer coordinates.
   * \param[in] bits Bits per coordinate. Use the same depth to decode.
   * \return Hilbert index for \p point at depth \p bits.
   * \pre Hilbertizer::isValidPoint(point, bits).
   * \note Does not validate inputs.
   */
  AXOM_HOST_DEVICE static IndexType toIndex(const PointType& point, int bits)
  {
    WorkPoint axes = convertPoint<WorkType>(point);

    // Undo lower-bit flips and axis exchanges, from the largest regions downwards.
    for(int bit = bits - 1; bit > 0; --bit)
    {
      const WorkType q = WorkType {1} << bit;
      const WorkType p = q - 1;
      for(int axis = 0; axis < DIM; ++axis)
      {
        if(axes[axis] & q)
        {
          axes[0] ^= p;
        }
        else
        {
          const WorkType exchange = (axes[0] ^ axes[axis]) & p;
          axes[0] ^= exchange;
          axes[axis] ^= exchange;
        }
      }
    }

    // Apply the inverse Gray-code transform across all index bits.
    for(int axis = 1; axis < DIM; ++axis)
    {
      axes[axis] ^= axes[axis - 1];
    }
    WorkType correction = 0;
    for(int bit = bits - 1; bit > 0; --bit)
    {
      const WorkType q = WorkType {1} << bit;
      if(axes[DIM - 1] & q)
      {
        correction ^= q - 1;
      }
    }
    for(int axis = 0; axis < DIM; ++axis)
    {
      axes[axis] ^= correction;
    }

    // Put axis zero's index bits first in each group that Mortonizer packs.
    reverseAxes(axes);
    return static_cast<IndexType>(MortonizerType::mortonize(axes));
  }

  /*!
   * \brief Decodes a Hilbert index to an integer lattice point.
   *
   * \param[in] index Hilbert index to decode.
   * \param[in] bits Bits per coordinate. Must match the depth used to encode.
   * \return Integer lattice point for \p index at depth \p bits.
   * \pre Hilbertizer::isValidIndex(index, bits).
   * \note Does not validate inputs.
   */
  AXOM_HOST_DEVICE static PointType toPoint(IndexType index, int bits)
  {
    WorkPoint axes = MortonizerType::demortonize(static_cast<WorkType>(index));
    // Restore axis zero's word from the most significant bit of each group.
    reverseAxes(axes);

    // Apply h ^ (h >> 1) to the index bits spread across these words.
    // Axis zero uses the last word shifted right, rather than a previous word.
    const WorkType correction = axes[DIM - 1] >> 1;
    for(int axis = DIM - 1; axis > 0; --axis)
    {
      axes[axis] ^= axes[axis - 1];
    }
    axes[0] ^= correction;

    // Correct lower-bit flips and axis exchanges, from the smallest regions upwards.
    WorkType axis0 = axes[0];
    for(int bit = 1; bit < bits; ++bit)
    {
      const WorkType q = WorkType {1} << bit;
      const WorkType p = q - 1;
      for(int axis = DIM - 1; axis > 0; --axis)
      {
        if(axes[axis] & q)
        {
          axis0 ^= p;
        }
        else
        {
          const WorkType exchange = (axis0 ^ axes[axis]) & p;
          axis0 ^= exchange;
          axes[axis] ^= exchange;
        }
      }
      // Axis zero's exchange is a no-op. This formulation helps the compiler to use a conditional move.
      axis0 = (axis0 & q) ? (axis0 ^ p) : axis0;
    }
    axes[0] = axis0;

    return convertPoint<CoordType>(axes);
  }

private:
  /// Casts point coordinates without a conversion loop or temporary array.
  template <typename To, typename From>
  AXOM_HOST_DEVICE static primal::Point<To, DIM> convertPoint(const primal::Point<From, DIM>& point)
  {
    // Explicit assignments help avoid compiler checks for vector instruction support.
    primal::Point<To, DIM> converted;
    converted[0] = static_cast<To>(point[0]);
    converted[1] = static_cast<To>(point[1]);
    if constexpr(DIM > 2)
    {
      converted[2] = static_cast<To>(point[2]);
    }
    if constexpr(DIM > 3)
    {
      converted[3] = static_cast<To>(point[3]);
    }
    return converted;
  }

  /// Reverses axis order in place for the supported dimensions, 2, 3, and 4.
  AXOM_HOST_DEVICE static void reverseAxes(WorkPoint& axes)
  {
    if constexpr(DIM > 1)
    {
      axom::utilities::swap(axes[0], axes[DIM - 1]);
    }
    if constexpr(DIM > 3)
    {
      axom::utilities::swap(axes[1], axes[DIM - 2]);
    }
  }
};
}  // namespace axom::spin::testing
