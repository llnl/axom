// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file MortonIndex
 *
 * \brief Classes and functions to convert between points on an integer grid and
 *  their unidimensional MortonIndex.
 *
 * Also has some utility functions for 'mortonizing' and 'demortonizing' points and a PointHash
 * functor class that can be used as a std::hash for unordered_maps and axom::FlatMap
 */

#ifndef AXOM_SPIN_MORTON_INDEX_HPP_
#define AXOM_SPIN_MORTON_INDEX_HPP_

#include "axom/config.hpp"
#include "axom/core/Types.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/NumericLimits.hpp"
#include "axom/primal/geometry/Point.hpp"

#include <type_traits>

namespace
{
/*!
 * \brief A helper type trait class for the Mortonizer class
 *
 * This type trait helps determine the number of iterations required
 * to convert between Points and a 1D Morton index.
 */
template <typename IntegerType>
struct NumReps
{
  static constexpr int value = 5;
};

template <>
struct NumReps<std::int64_t>
{
  static constexpr int value = 5;
};
template <>
struct NumReps<std::uint64_t>
{
  static constexpr int value = 5;
};

template <>
struct NumReps<std::int32_t>
{
  static constexpr int value = 4;
};
template <>
struct NumReps<std::uint32_t>
{
  static constexpr int value = 4;
};

template <>
struct NumReps<std::int16_t>
{
  static constexpr int value = 3;
};
template <>
struct NumReps<std::uint16_t>
{
  static constexpr int value = 3;
};

template <>
struct NumReps<std::int8_t>
{
  static constexpr int value = 2;
};
template <>
struct NumReps<std::uint8_t>
{
  static constexpr int value = 2;
};
}  // namespace

namespace axom
{
namespace spin
{
/*!
 * \class
 * \brief Base class for Dimension independent Morton indexing
 *
 * \note Uses CRTP to access dimension-dependent data from the derived class
 * \note This class only works for integral CoordTypes
 */
template <typename CoordType, typename MortonIndexType, typename Derived>
struct MortonBase
{
  // static assert to ensure we only instantiate on integral types
  static_assert(std::is_integral<CoordType>::value,
                "Coordtype must be integral for Morton indexing");
  static_assert(std::is_integral<MortonIndexType>::value,
                "MortonIndexType must be integral for Morton indexing");

  // we get undefined behavior w/ our bit shifts if using signed types
  static_assert(std::is_unsigned<MortonIndexType>::value,
                "MortonIndexType must be an unsigned type");

private:
  // Magic numbers for efficient base-2 log-like function -- maxSetBit()
  static const CoordType MaxBit_B[];
  static const int MaxBit_S[];

protected:
  /*!
   * \brief Expands bits in bitwise representation of an integral type
   *        and zero-fills the holes
   *
   * \param [in] x The integer type that we are expanding
   * \return A zero-filled expanded MortonIndex
   * In dimension D, it adds (D-1) zeros between each bit,
   * so, e.g. in 2D, 6 == 0b0110 becomes 0b*0*1*1*0 == 0b00010100 == 20
   */
  AXOM_HOST_DEVICE
  static MortonIndexType expandBits(MortonIndexType x)
  {
    for(int i = Derived::EXPAND_MAX_ITER; i >= 0; --i)
    {
      x = (x | (x << Derived::GetS(i))) & Derived::GetB(i);
    }

    return x;
  }

  /*!
   * \brief Contracts bits in bitwise representation of x
   *
   * \param [in] x The Morton index that we are contracting
   * \return A contracted MortonIndex
   *
   * In dimension D, it retains every (D-1)\f$^th\f$ bit,
   * so, e.g. in 2D, 20 = 0b00010100 == 0b*0*1*1*0 becomes  0b0110 = 6
   */
  AXOM_HOST_DEVICE
  static MortonIndexType contractBits(MortonIndexType x)
  {
    for(int i = 0; i < Derived::CONTRACT_MAX_ITER; ++i)
    {
      x = (x | (x >> Derived::GetS(i))) & Derived::GetB(i + 1);
    }

    return x;
  }

public:
  /*! \brief Finds the index of the maximum set bit (MSB) in an integral type */
  static int maxSetBit(CoordType x)
  {
    CoordType res = 0;
    for(int i = Derived::CONTRACT_MAX_ITER; i >= 0; --i)
    {
      if(x & MaxBit_B[i])
      {
        x >>= MaxBit_S[i];
        res |= MaxBit_S[i];
      }
    }

    return static_cast<int>(res);
  }
};

/*!
 * \class Mortonizer
 * \brief Helper class for MortonIndexing of a point's coordinate
 *
 * The Morton index of a point interleaves the bits of its coordinates
 * (with the least significant bit coming from the x-coordinate
 * E.g. if we have a point in 2D (6,3) == (0b0110, 0b0011) in its binary
 * representation.If we label the individual set bits we get (0b0ab0, 0b00yz).
 * When we expand the bits by inserting a 0 (denoted with *),
 * we get (0b*0*a*b*0,0b*0*0*y*z).
 * Finally, by shifting the y-coordinate and interleaving, we get
 * the Morton index of this point: 0b000yazb0 == 0b00011110 ==  30
 */
template <typename CoordType, typename MortonIndexType, int DIM>
struct Mortonizer;

/*!
 * \class
 * \brief A 2D specialization of Mortonizer
 *
 * This class converts between 2D integer Grid points and MortonIndexes
 * where a grid point is a point whose coordinates are integral types
 * Expand bits will add a 0 between every bit and contract bits
 * will remove every other bit
 * \see Mortonizer
 */
template <typename CoordType, typename MortonIndexType>
struct Mortonizer<CoordType, MortonIndexType, 2>
  : public MortonBase<CoordType, MortonIndexType, Mortonizer<CoordType, MortonIndexType, 2>>
{
  using self = Mortonizer<CoordType, MortonIndexType, 2>;
  using Base = MortonBase<CoordType, MortonIndexType, self>;

  /*! The dimension of the Mortonizer */
  static constexpr int NDIM = 2;

  /*! The number of bits in a CoordType  */
  static constexpr int COORD_BITS = axom::numeric_limits<CoordType>::digits;

  /*! The number of bits in a MortonIndex  */
  static constexpr int MORTON_BITS = axom::numeric_limits<MortonIndexType>::digits;

  /*! The number of representable Morton bits per dimension */
  static constexpr int MB_PER_DIM = MORTON_BITS / NDIM;

  /*!
     * The maximum number of unique bits from each coordinate of type CoordType
     *  that can be represented in a MortonIndex.
     *
     * \note If we are use Mortonizer as a (one-way) hash function, it is ok to
     *  use more bits. But, if we would like to be able to reverse the
     *  MortonIndex, then we cannot safely use more than MAX_UNIQUE_BITS per
     *  coordinate.
     *
     */
  static constexpr int MAX_UNIQUE_BITS = (MB_PER_DIM < COORD_BITS) ? MB_PER_DIM : COORD_BITS;

  /*!
     * The number of iterations required for converting from MortonIndexes to
     * CoordType using the bit interleaving algorithm in MortonBase.
     */
  static constexpr int CONTRACT_MAX_ITER = NumReps<MortonIndexType>::value;

  /*!
     * The number of iterations required for converting between CoordTypes
     * and MortonIndexes using the bit interleaving algorithm in MortonBase.
     */
  static constexpr int EXPAND_MAX_ITER = NumReps<MortonIndexType>::value;

  // Magic numbers in 2D
  AXOM_HOST_DEVICE static MortonIndexType GetB(int i)
  {
    constexpr MortonIndexType B[] = {static_cast<MortonIndexType>(0x5555555555555555),  // 0101'0101
                                     static_cast<MortonIndexType>(0x3333333333333333),  // 0011'0011
                                     static_cast<MortonIndexType>(0x0F0F0F0F0F0F0F0F),  // 0000'1111
                                     static_cast<MortonIndexType>(0x00FF00FF00FF00FF),  // 0x8  1x8
                                     static_cast<MortonIndexType>(0x0000FFFF0000FFFF),  // 0x16 1x16
                                     static_cast<MortonIndexType>(0x00000000FFFFFFFF)};  // 0x32 1x32;
    return B[i];
  }

  AXOM_HOST_DEVICE static int GetS(int i)
  {
    constexpr int S[] = {1, 2, 4, 8, 16, 32};
    return S[i];
  }

  /*!
   * \brief A function to convert a 2D point to a Morton index
   *
   * Morton indexing interleaves the bits of the point's coordinates
   * \param [in] x The x-coordinate of the point
   * \param [in] y The y-coordinate of the point
   * \pre CoordType must be an integral type
   * \pre x and y must be positive with value less than \f$ 2^{16} \f$
   * \note These preconditions can easily be relaxed, if desired
   *       (e.g. to 32 bits for 64 bit Morton indices)
   * \return The MortonIndex of the 2D point
   */
  AXOM_HOST_DEVICE
  static inline MortonIndexType mortonize(CoordType x, CoordType y)
  {
    return (Base::expandBits(x) | (Base::expandBits(y) << 1));
  }

  /*!
   * \brief A function to convert a 2D point to a Morton index
   *
   * \see mortonize(CoordType, CoordType)
   */
  AXOM_HOST_DEVICE
  static inline MortonIndexType mortonize(const primal::Point<CoordType, NDIM>& pt)
  {
    return (Base::expandBits(pt[0]) | (Base::expandBits(pt[1]) << 1));
  }

  /*!
   * \brief A function to convert a Morton index back to a 2D point
   *
   * \param [in] morton The MortonIndex of the desired point
   * \param [out] x The x-coordinate of the point
   * \param [out] y The y-coordinate of the point
   *  Morton indexing interleaves the bits of the point's coordinates
   * \note The point's coordinates are returned in the x and y parameters
   */
  AXOM_HOST_DEVICE
  static inline void demortonize(MortonIndexType morton, CoordType& x, CoordType& y)
  {
    const MortonIndexType b0 = GetB(0);

    x = static_cast<CoordType>(Base::contractBits(morton & b0));
    y = static_cast<CoordType>(Base::contractBits((morton >> 1) & b0));
  }

  /*!
   * \brief A function to convert a Morton index back to a 2D point
   *
   * \see demortonize(MortonIndex,CoordType,CoordType)
   */
  AXOM_HOST_DEVICE
  static inline primal::Point<CoordType, NDIM> demortonize(MortonIndexType morton)
  {
    primal::Point<CoordType, NDIM> pt;
    demortonize(morton, pt[0], pt[1]);
    return pt;
  }

  /*!
   * \brief returns the maximum number of bits per coordinate
   */
  static int maxBitsPerCoord() { return MAX_UNIQUE_BITS; }
};

/*!
 * \class
 * \brief A 3D specialization of Mortonizer
 *
 * This class converts between 3D integer grid points and MortonIndexes
 * where a grid point is a point whose coordinates are integral types.
 * Expand bits will add a 0 between every bit
 * and contract bits will remove every other bit
 * \see Mortonizer
 */
template <typename CoordType, typename MortonIndexType>
struct Mortonizer<CoordType, MortonIndexType, 3>
  : public MortonBase<CoordType, MortonIndexType, Mortonizer<CoordType, MortonIndexType, 3>>
{
  using self = Mortonizer<CoordType, MortonIndexType, 3>;
  using Base = MortonBase<CoordType, MortonIndexType, self>;

  /*! The dimension of the Mortonizer */
  constexpr static int NDIM = 3;

  /*! The number of bits in a CoordType  */
  constexpr static int COORD_BITS = axom::numeric_limits<CoordType>::digits;

  /*! The number of bits in a MortonIndex  */
  constexpr static int MORTON_BITS = axom::numeric_limits<MortonIndexType>::digits;

  /*! The number of representable morton bits per dimension */
  constexpr static int MB_PER_DIM = MORTON_BITS / NDIM;

  /*!
    * The maximum number of unique bits from each coordinate of type CoordType
    *  that can be represented in a MortonIndex.
    *  \note If we are use Mortonizer as a (one-way) hash function,
    *        it is ok to use more bits. But, if we would like to be
    *        able to reverse the MortonIndex, then we cannot safely use
    *        more than MAX_UNIQUE_BITS per coordinate.
    */
  constexpr static int MAX_UNIQUE_BITS = (MB_PER_DIM < COORD_BITS) ? MB_PER_DIM : COORD_BITS;

  /*!
    * The number of iterations required for converting from MortonIndexes
    * to CoordType using the bit interleaving algorithm in MortonBase.
    */
  constexpr static int CONTRACT_MAX_ITER = NumReps<MortonIndexType>::value;

  /*!
    * The number of iterations required for converting between CoordTypes
    * and MortonIndexes using the bit interleaving algorithm in MortonBase.
    */
  constexpr static int EXPAND_MAX_ITER = NumReps<MortonIndexType>::value - 1;

  // Magic numbers in 3D from C. Ericson's Real Time Collision Detection book
  AXOM_HOST_DEVICE static MortonIndexType GetB(int i)
  {
    constexpr MortonIndexType B[] = {
      static_cast<MortonIndexType>(0x9249249249249249),   // 0010'0100'1001'0010'0100'1001
      static_cast<MortonIndexType>(0x30C30C30C30C30C3),   // 0000'1100'0011'0000'1100'0011
      static_cast<MortonIndexType>(0xF00F00F00F00F00F),   // 0000'0000'1111'0000'0000'1111
      static_cast<MortonIndexType>(0x00FF0000FF0000FF),   // 0000'0000'0000'0000'1111'1111
      static_cast<MortonIndexType>(0xFFFF00000000FFFF),   // x16
      static_cast<MortonIndexType>(0x00000000FFFFFFFF)};  // x32
    return B[i];
  }

  AXOM_HOST_DEVICE static int GetS(int i)
  {
    constexpr int S[] = {2, 4, 8, 16, 32, 0};
    return S[i];
  }

  /*!
   * \brief A function to convert a 3D point to a Morton index
   *
   *  Morton indexing interleaves the bits of the point's coordinates
   * \param [in] x The x-coordinate of the point
   * \param [in] y The y-coordinate of the point
   * \param [in] z The z-coordinate of the point
   * \pre CoordType must be an integral type
   * \pre x and y must be positive with value less than \f$ 2^{10} \f$
   * \note These preconditions can easily be relaxed, if desired
   * (e.g. to 21 bits for 64 bit Morton indices)
   * \return The MortonIndex of the 2D point
   */
  AXOM_HOST_DEVICE
  static inline MortonIndexType mortonize(CoordType x, CoordType y, CoordType z)
  {
    const MortonIndexType b5 = GetB(5);

    return (Base::expandBits(x & b5) | (Base::expandBits(y & b5) << 1) |
            (Base::expandBits(z & b5) << 2));
  }

  /*!
   * \brief A function to convert a 3D point to a Morton index
   *
   * \see mortonize(CoordType, CoordType, CoordType)
   */
  AXOM_HOST_DEVICE
  static inline MortonIndexType mortonize(const primal::Point<CoordType, NDIM>& pt)
  {
    const MortonIndexType b5 = GetB(5);

    return (Base::expandBits(pt[0] & b5) | (Base::expandBits(pt[1] & b5) << 1) |
            (Base::expandBits(pt[2] & b5) << 2));
  }

  /*!
   * \brief A function to convert a Morton index back to a 3D point
   *
   * \param [in] morton The MortonIndex of the desired point
   * \param [out] x The x-coordinate of the point
   * \param [out] y The y-coordinate of the point
   * \param [out] z The z-coordinate of the point
   *  Morton indexing interleaves the bits of the point's coordinates
   * \note The point's coordinates are returned in the x, y and z parameters
   */
  AXOM_HOST_DEVICE
  static inline void demortonize(MortonIndexType morton, CoordType& x, CoordType& y, CoordType& z)
  {
    const MortonIndexType b0 = GetB(0);

    x = static_cast<CoordType>(Base::contractBits(morton & b0));
    y = static_cast<CoordType>(Base::contractBits((morton >> 1) & b0));
    z = static_cast<CoordType>(Base::contractBits((morton >> 2) & b0));
  }

  /*!
   * \brief A function to convert a Morton index back to a 3D point
   *
   * \see demortonize(MortonIndex,CoordType,CoordType,CoordType)
   */
  AXOM_HOST_DEVICE
  static inline primal::Point<CoordType, NDIM> demortonize(MortonIndexType morton)
  {
    primal::Point<CoordType, NDIM> pt;
    demortonize(morton, pt[0], pt[1], pt[2]);
    return pt;
  }

  /*!
   * \brief returns the maximum number of bits per coordinate
   */
  static int maxBitsPerCoord() { return MAX_UNIQUE_BITS; }
};

template <typename CoordType, typename MortonIndexType, typename Derived>
const CoordType MortonBase<CoordType, MortonIndexType, Derived>::MaxBit_B[] = {
  static_cast<CoordType>(0x2),
  static_cast<CoordType>(0xC),
  static_cast<CoordType>(0xF0),
  static_cast<CoordType>(0xFF00),
  static_cast<CoordType>(0xFFFF0000),
  static_cast<CoordType>(0xFFFFFFFF00000000)};

template <typename CoordType, typename MortonIndexType, typename Derived>
const int MortonBase<CoordType, MortonIndexType, Derived>::MaxBit_S[] = {1, 2, 4, 8, 16, 32};

/*!
 * \brief A helper function to convert a point directly to a MortonIndex
 *
 * \return The Morton index of the point
 */
template <typename MortonIndexType, typename CoordType, int DIM>
AXOM_HOST_DEVICE inline MortonIndexType convertPointToMorton(const primal::Point<CoordType, DIM>& pt)
{
  return Mortonizer<CoordType, MortonIndexType, DIM>::mortonize(pt);
}

/*!
 * \brief A helper function to convert a MortonIndex back to a point
 *
 * \return The demortonized Point
 */
template <typename CoordType, int DIM, typename MortonIndexType>
inline primal::Point<CoordType, DIM> convertMortonToPoint(MortonIndexType idx)
{
  return Mortonizer<CoordType, MortonIndexType, DIM>::demortonize(idx);
}

/*!
 * \class
 * \brief A functor class for Mortonizing points.
 *
 * This can be used as a hashing function for Points in dimensions 1..4
 */
template <typename CoordType>
struct PointHash
{
  using MortonIndex = std::size_t;
  using result_type = std::size_t;

  /*!
   * \brief Mortonizes a coordinate (viewed as a 1D point)
   *
   * \note This is a no-op and is provided for genericity in point dimension
   * \param [in] coord The coordinate of the
   * \returns The morton index of the 1D point
   */
  std::size_t operator()(CoordType& coord) const { return static_cast<std::size_t>(coord); }

  /*!
   * \brief Mortonizes a 1D point
   *
   * \note This is a no-op and is provided for genericity in point dimension
   * \param [in] pt The 1D point
   * \returns The morton index of the point
   */
  std::size_t operator()(primal::Point<CoordType, 1> const& pt) const
  {
    return static_cast<std::size_t>(pt[0]);
  }

  /*!
   * \brief Mortonizes a 2D point
   *
   * \param [in] pt The 2D point
   * \returns The morton index of the point
   */
  std::size_t operator()(primal::Point<CoordType, 2> const& pt) const
  {
    return Mortonizer<CoordType, MortonIndex, 2>::mortonize(pt);
  }

  /*!
   * \brief Mortonizes a 3D point
   *
   * \param [in] pt The 3D point
   * \returns The morton index of the point
   */
  std::size_t operator()(primal::Point<CoordType, 3> const& pt) const
  {
    return Mortonizer<CoordType, MortonIndex, 3>::mortonize(pt);
  }

  /*!
   * \brief Mortonizes a 4D point
   *
   * \param [in] pt The 4D point
   * \returns A morton index of the point
   */
  std::size_t operator()(primal::Point<CoordType, 4> const& pt) const
  {
    // Mortonize two morton indices for a 4D morton index

    using Pt2M = primal::Point<MortonIndex, 2>;

    Pt2M pMorton = Pt2M::make_point(Mortonizer<CoordType, MortonIndex, 2>::mortonize(pt[0], pt[2]),
                                    Mortonizer<CoordType, MortonIndex, 2>::mortonize(pt[1], pt[3]));

    return Mortonizer<MortonIndex, MortonIndex, 2>::mortonize(pMorton);
  }
};

}  // end namespace spin
}  // end namespace axom

#endif  // AXOM_SPIN_MORTON_INDEX_HPP_
