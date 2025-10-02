// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_BROOD__HPP_
#define AXOM_SPIN_BROOD__HPP_

#include "axom/config.hpp"
#include "axom/spin/MortonIndex.hpp"

namespace axom
{
namespace spin
{
/**
 * \class
 * \brief Helper class to handle subindexing of block data within octree siblings
 *
 * \note A brood is a collection of siblings that are generated simultaneously.
 * \note This class converts a grid point at the given level into a brood index of the point.
 *
 * The base brood is the MortonIndex of the grid point's octree parent and its offset index
 *  is obtained by interleaving the least significant bit of its coordinates.
 */
template <typename GridPt, typename MortonIndexType>
struct Brood
{
  enum
  {
    DIM = GridPt::DIMENSION,
    BROOD_BITMASK = (1 << DIM) - 1
  };

  using MortonizerType = spin::Mortonizer<typename GridPt::CoordType, MortonIndexType, DIM>;

  /**
   * \brief Constructor for a brood offset relative to the given grid point pt
   *
   * \param [in] pt The grid point within the octree level
   */
  Brood(const GridPt& pt)
  {
    m_broodIdx = MortonizerType::mortonize(pt);
    m_offset = m_broodIdx & BROOD_BITMASK;
    m_broodIdx >>= DIM;
  }

  /**
   * \brief Accessor for the base point of the entire brood
   *
   * \return MortonIndex of the base point
   */
  const MortonIndexType& base() const { return m_broodIdx; }

  /** \brief Accessor for the offset of the point within the brood */
  const int& offset() const { return m_offset; }

  /// \brief Reconstruct a grid point from a brood's Morton index and an offset
  static GridPt reconstructGridPt(MortonIndexType morton, int offset)
  {
    return static_cast<GridPt>(
      MortonizerType::demortonize(static_cast<MortonIndexType>((morton << DIM) + offset)));
  }

private:
  MortonIndexType m_broodIdx;  /// MortonIndex of the base point of all blocks within the brood
  int m_offset;                /// Index of the block within the brood. Value is in [0, 2^DIM)
};

/**
 * \class
 * \brief Template specialization of Brood which does not use MortonIndexing
 *
 * \note A brood is a collection of siblings that are generated simultaneously.
 * \note This class converts a grid point at the given level into a brood index
 *  of the point.The base brood point has the coordinates of the grid point's
 *  octree parent and its offset index is obtained by interleaving the least
 *  significant bit of its coordinates in each dimension.
 *
 *  \see Brood
 */
template <typename GridPt>
struct Brood<GridPt, GridPt>
{
  /**
   * \brief Constructor for a brood offset relative to the given grid point pt
   *
   * \param [in] pt The grid point within the octree level
   */
  Brood(const GridPt& pt) : m_broodPt(pt.array() / 2), m_offset(0)
  {
    for(int i = 0; i < GridPt::DIMENSION; ++i)
    {
      m_offset |= (pt[i] & 1) << i;  // interleave the least significant bits
    }
  }

  /// \brief Accessor for the base point of the entire brood
  const GridPt& base() const { return m_broodPt; }

  /// \brief Accessor for the index of the point within the brood
  const int& offset() const { return m_offset; }

  /// \brief Reconstruct a grid point from a brood's base point and an offset
  static GridPt reconstructGridPt(const GridPt& pt, int offset)
  {
    // shift and add offset to each coordinate
    GridPt retPt;
    for(int i = 0; i < GridPt::DIMENSION; ++i)
    {
      retPt[i] = (pt[i] << 1) + (offset & (1 << i) ? 1 : 0);
    }

    return retPt;
  }

private:
  GridPt m_broodPt;  /// Base point of all blocks within the brood
  int m_offset;      /// Index of the block within the brood. Value is in [0, 2^DIM)
};

}  // end namespace spin
}  // end namespace axom

#endif  // AXOM_SPIN_BROOD__HPP_
