// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_SPATIAL_OCTREE__HPP_
#define AXOM_SPIN_SPATIAL_OCTREE__HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"

#include "axom/spin/OctreeBase.hpp"

namespace axom
{
namespace spin
{
/**
 * \class SpatialOctree
 * \brief Adds spatial extents to an OctreeBase, allowing point location
 */
template <int DIM, typename BlockDataType>
class SpatialOctree : public OctreeBase<DIM, BlockDataType>
{
public:
  using GeometricBoundingBox = primal::BoundingBox<double, DIM>;
  using SpacePt = primal::Point<double, DIM>;
  using SpaceVector = primal::Vector<double, DIM>;

  using BaseOctree = OctreeBase<DIM, BlockDataType>;

  using GridPt = typename BaseOctree::GridPt;
  using CoordType = typename GridPt::CoordType;

  using BlockIndex = typename BaseOctree::BlockIndex;

  using SpaceVectorLevelMap = slam::Map<SpaceVector>;

public:
  /**
   * \brief Construct a spatial octree from a spatial bounding box
   *
   * \param [in] bb The spatial extent to be indexed by the octree
   */
  SpatialOctree(const GeometricBoundingBox& bb)
    : BaseOctree()
    , m_deltaLevelMap(&this->m_levels)
    , m_invDeltaLevelMap(&this->m_levels)
    , m_boundingBox(bb)
  {
    // Cache the extents of a grid cell at each level of resolution
    const SpaceVector bbRange = m_boundingBox.range();
    for(int lev = 0; lev < this->m_levels.size(); ++lev)
    {
      m_deltaLevelMap[lev] = bbRange / static_cast<double>(CoordType(1) << lev);

      for(int dim = 0; dim < DIM; ++dim)
      {
        m_invDeltaLevelMap[lev][dim] = 1. / m_deltaLevelMap[lev][dim];
      }
    }
  }

  /**
   * \brief Returns a reference to the octree's bounding box
   *  (i.e. the bounding box of the root block)
   */
  const GeometricBoundingBox& boundingBox() const { return m_boundingBox; }

  /**
   * \brief Return the spatial bounding box of a grid cell
   * at the given level or resolution
   */
  GeometricBoundingBox blockBoundingBox(const BlockIndex& block) const
  {
    return blockBoundingBox(block.pt(), block.level());
  }

  /**
   * \brief Return the spatial bounding box of a grid cell
   * at the given level or resolution
   */
  GeometricBoundingBox blockBoundingBox(const GridPt& gridPt, int level) const
  {
    const SpaceVector& deltaVec = m_deltaLevelMap[level];

    SpacePt lower(m_boundingBox.getMin());
    SpacePt upper(m_boundingBox.getMin());
    for(int i = 0; i < DIM; ++i)
    {
      lower[i] += gridPt[i] * deltaVec[i];
      upper[i] += (gridPt[i] + 1) * deltaVec[i];
    }

    return GeometricBoundingBox(lower, upper);
  }

  /**
   * Returns the width of an octree block at level of resolution level
   * \param level The level of resolution that we are checking
   */
  const SpaceVector& spacingAtLevel(int level) const { return m_deltaLevelMap[level]; }

  /**
   * \brief Finds the index of the leaf block covering the query point pt
   *
   * \param [in] pt The query point in space
   * \param [in] startingLevel (Optional) starting level for the query
   * \pre pt must be in the bounding box of the octree (i.e. boundingBox.contains(pt) == true )
   * \note The collection of leaves covers the bounding box, and the interiors
   * of the leaves do not intersect, so every point in the bounding box 
   * should be located in a unique leaf block.
   * \note We are assuming a half-open interval on the bounding boxes.
   * \return The block index (i.e. grid point and level) of the leaf block containing the query point
   */
  BlockIndex findLeafBlock(const SpacePt& pt, int startingLevel = -1) const
  {
    SLIC_ASSERT_MSG(
      m_boundingBox.contains(pt),
      "SpatialOctree::findLeafNode -- Did not find " << pt << " in bounding box " << m_boundingBox);

    // Perform binary search on levels to find the leaf block containing the point
    CoordType minLev = 0;
    CoordType maxLev = this->maxLeafLevel();
    CoordType lev = (startingLevel == -1) ? maxLev >> 1 : startingLevel;

    while(minLev <= maxLev)
    {
      GridPt gridPt = findGridCellAtLevel(pt, static_cast<int>(lev));
      switch(this->blockStatus(gridPt, lev))
      {
      case BlockNotInTree:
        // Block must be in coarser levels -- update upper bound
        maxLev = lev - 1;
        lev = (maxLev + minLev) >> 1;
        break;
      case InternalBlock:
        // Block must be in deeper levels -- update lower bound
        minLev = lev + 1;
        lev = (maxLev + minLev) >> 1;
        break;
      case LeafBlock:
        return BlockIndex(gridPt, static_cast<int>(lev));
      }
    }

    SLIC_ASSERT_MSG(false, "Point " << pt << " not found in a leaf block of the octree");

    return BlockIndex::invalid_index();
  }

  /**
   * \brief Utility function to find the quantized grid cell at level lev for query point pt
   *
   * \param [in] pt The point at which we are querying.
   * \param [in] lev The level or resolution.
   * \pre \f$ 0 \le lev < octree.maxLeafLevel() \f$
   * \post Each coordinate of the returned gridPt is in range  \f$ [0, 2^{lev}) \f$
   * \return The grid point of the block covering this point at this level
   * \internal
   * \todo KW: Should this function be protected? Is it generally useful?
   */
  GridPt findGridCellAtLevel(const SpacePt& pt, int lev) const
  {
    GridPt quantizedPt;

    const SpacePt& bbMin = m_boundingBox.getMin();
    const SpaceVector& invDelta = m_invDeltaLevelMap[lev];
    const CoordType highestCell = this->maxCoordAtLevel(lev);

    for(int i = 0; i < DIM; ++i)
    {
      // Note: quantCell is always positive, and within range of CoordType
      //       so truncating is equivalent to the floor function
      // Note: we need to clamp to avoid setting coordinates past the upper
      // boundaries
      const CoordType quantCell = static_cast<CoordType>((pt[i] - bbMin[i]) * invDelta[i]);
      quantizedPt[i] = std::min(quantCell, highestCell);
    }

    return quantizedPt;
  }

private:
  DISABLE_COPY_AND_ASSIGNMENT(SpatialOctree);
  DISABLE_MOVE_AND_ASSIGNMENT(SpatialOctree);

protected:
  SpaceVectorLevelMap m_deltaLevelMap;     // The width of a cell at each level or resolution
  SpaceVectorLevelMap m_invDeltaLevelMap;  // Its inverse is useful for quantizing
  GeometricBoundingBox m_boundingBox;
};

}  // end namespace spin
}  // end namespace axom

#endif  // AXOM_SPIN_SPATIAL_OCTREE__HPP_
