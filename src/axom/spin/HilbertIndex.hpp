// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/*! 
 * \file HilbertIndex.hpp
 * \brief Converts bounded integer grid points to and from Hilbert indices.
 */

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/NumericLimits.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/slic/interface/slic.hpp"
#include "axom/spin/internal/HilbertCodec.hpp"

#include <type_traits>

namespace axom::spin
{
/*!
 * \brief Maps 2D, 3D, or 4D integer lattice points to and from a Hilbert index.
 *
 * The index is a point's position on a path that visits every point in the grid once.
 * Consecutive indices visit points one grid step apart along one axis.
 * At depth 1, indices 0 through 3 visit the 2D points (0,0), (0,1), (1,1), (1,0).
 *
 * At a valid bit depth, each point in [0, 2^bits - 1]^DIM maps to exactly one
 * index in [0, 2^(DIM * bits) - 1], and vice versa. Use the same depth to encode
 * and decode, and store the dimension and depth with the key.
 * Every conversion and check takes the depth explicitly, since a point's index
 * changes with the depth.
 *
 * Quantize floating-point positions to integer grid coordinates before encoding.
 * Check that the resulting coordinates are nonnegative and fit the selected depth.
 * Conversion functions require valid inputs and don't clamp them.
 *
 * Signed types use only nonnegative values. With 64-bit coordinates, uint64_t
 * keys allow 32, 21, and 16 bits per axis in 2D, 3D, and 4D.
 * The corresponding limits for int64_t keys are 31, 21, and 15 bits per axis.
 * Narrower coordinate types can reduce these limits.
 *
 * Encoding and decoding use a table-based codec. Host tables are constexpr,
 * while GPU tables use read-only global device memory.
 * The tables follow J. Skilling's algorithm in "Programming the Hilbert curve",
 * AIP Conference Proceedings 707, 381 (2004), https://doi.org/10.1063/1.1751381.
 * The scalar Skilling codec is a test and benchmark reference, not a library option.
 * \see internal::hilbert::HilbertFSMTransform for the table implementation.
 *
 * \tparam CoordType Integral coordinate type other than bool, at most 64 bits.
 * \tparam IndexType Integral index type other than bool, at most 64 bits.
 * \tparam DIM Supported spatial dimension, 2, 3, or 4.
 */
template <typename CoordType, typename IndexType, int DIM>
struct Hilbertizer
{
  static_assert(DIM >= 2 && DIM <= 4, "Hilbert indexing supports dimensions 2, 3, and 4");
  static_assert(std::is_integral<CoordType>::value && !std::is_same<CoordType, bool>::value,
                "Hilbert coordinates must be non-bool integral values");
  static_assert(std::is_integral<IndexType>::value && !std::is_same<IndexType, bool>::value,
                "Hilbert indices must be non-bool integral values");
  static_assert(sizeof(CoordType) <= 8 && sizeof(IndexType) <= 8,
                "Hilbert indexing supports integral types of at most 64 bits");

  using PointType = primal::Point<CoordType, DIM>;
  using UnsignedIndex = typename std::make_unsigned<IndexType>::type;
  static constexpr int NDIM = DIM;
  static constexpr int INDEX_BITS = axom::numeric_limits<IndexType>::digits;
  static constexpr int COORD_BITS = axom::numeric_limits<CoordType>::digits;
  /// Maximum bits per coordinate without loss.
  static constexpr int MAX_UNIQUE_BITS =
    COORD_BITS < INDEX_BITS / DIM ? COORD_BITS : INDEX_BITS / DIM;

  /// Whether a depth is positive and fits the coordinate/index type combination.
  AXOM_HOST_DEVICE static bool isValidDepth(int bits)
  {
    return bits > 0 && bits <= MAX_UNIQUE_BITS;
  }

  /// Whether all coordinates are nonnegative and fit in the requested depth.
  AXOM_HOST_DEVICE static bool isValidPoint(const PointType& point, int bits)
  {
    using UnsignedCoord = typename std::make_unsigned<CoordType>::type;

    if(!isValidDepth(bits))
    {
      return false;
    }
    // DIM >= 2 keeps the shift strictly below the index type's width.
    const UnsignedIndex limit = UnsignedIndex {1} << bits;
    for(int axis = 0; axis < DIM; ++axis)
    {
      if constexpr(std::is_signed<CoordType>::value)
      {
        if(point[axis] < 0)
        {
          return false;
        }
      }
      // CoordType may be wider than IndexType, compare before narrowing.
      if(static_cast<UnsignedCoord>(point[axis]) >= limit)
      {
        return false;
      }
    }
    return true;
  }

  /// Whether an index is nonnegative and fits in DIM * bits bits.
  AXOM_HOST_DEVICE static bool isValidIndex(IndexType index, int bits)
  {
    if(!isValidDepth(bits))
    {
      return false;
    }
    if constexpr(std::is_signed<IndexType>::value)
    {
      if(index < 0)
      {
        return false;
      }
    }
    constexpr int unsigned_bits = axom::numeric_limits<UnsignedIndex>::digits;
    const auto mask = axom::numeric_limits<UnsignedIndex>::max() >> (unsigned_bits - DIM * bits);
    return static_cast<UnsignedIndex>(index) <= mask;
  }

  /*!
   * \brief Converts an integer lattice point to a Hilbert index.
   * \param[in] point Nonnegative integer coordinates.
   * \param[in] bits Bits per coordinate, in [1, maxBitsPerCoord()].
   * \pre isValidPoint(point, bits).
   */
  AXOM_HOST_DEVICE static IndexType toIndex(const PointType& point, int bits)
  {
    SLIC_ASSERT(isValidPoint(point, bits));
    return HilbertFSMTransform::toIndex(point, bits);
  }

  /*!
   * \brief Converts a Hilbert index back to an integer lattice point.
   * \param[in] index Hilbert index to decode.
   * \param[in] bits Depth used to encode the index, in [1, maxBitsPerCoord()].
   * \pre isValidIndex(index, bits).
   */
  AXOM_HOST_DEVICE static PointType toPoint(IndexType index, int bits)
  {
    SLIC_ASSERT(isValidIndex(index, bits));
    return HilbertFSMTransform::toPoint(index, bits);
  }

  /// Maximum supported depth.
  AXOM_HOST_DEVICE static constexpr int maxBitsPerCoord() { return MAX_UNIQUE_BITS; }

private:
  // Bit planes per table lookup, shared by host and GPU calls.
  static constexpr int FSM_LEVELS = DIM == 2 ? 3 : DIM == 3 ? 2 : 1;
  using HilbertFSMTransform =
    internal::hilbert::HilbertFSMTransform<CoordType, IndexType, DIM, FSM_LEVELS>;
};

/// Encodes a valid integer lattice point; see Hilbertizer::toIndex().
template <typename IndexType, typename CoordType, int DIM>
AXOM_HOST_DEVICE inline IndexType convertPointToHilbert(const primal::Point<CoordType, DIM>& point,
                                                        int bits)
{
  return Hilbertizer<CoordType, IndexType, DIM>::toIndex(point, bits);
}

/// Decodes a valid Hilbert index at the same depth; see Hilbertizer::toPoint().
template <typename CoordType, int DIM, typename IndexType>
AXOM_HOST_DEVICE inline primal::Point<CoordType, DIM> convertHilbertToPoint(IndexType index, int bits)
{
  return Hilbertizer<CoordType, IndexType, DIM>::toPoint(index, bits);
}
}  // namespace axom::spin
