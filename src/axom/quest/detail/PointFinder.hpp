// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_POINT_IN_CELL_POINT_FINDER_HPP_
#define AXOM_QUEST_POINT_IN_CELL_POINT_FINDER_HPP_

#include "axom/spin/ImplicitGrid.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"

namespace axom
{
namespace quest
{
// Predeclare mesh traits class
template <typename mesh_tag>
struct PointInCellTraits;

namespace detail
{
// Predeclare mesh wrapper class
template <typename mesh_tag>
class PointInCellMeshWrapper;

/*!
 * \class PointFinder
 *
 * \brief A class to encapsulate locating points within the cells
 * of a computational mesh.
 *
 * \tparam NDIMS The dimension of the mesh
 * \tparam mesh_tag A tag struct used to identify the mesh
 *
 * \note This class implements part of the functionality of \a PointInCell
 * \note This class assumes the existence of specialized implementations of
 * the following two classes for the provided \a mesh_tag:
 *   \arg axom::quest::PointInCellTraits
 *   \arg axom::quest::detail::PointInCellMeshWrapper
 */
template <int NDIMS, typename mesh_tag, typename ExecSpace>
class PointFinder
{
public:
  using MeshWrapperType = PointInCellMeshWrapper<mesh_tag>;
  using IndexType = typename MeshWrapperType::IndexType;

  using GridType = spin::ImplicitGrid<NDIMS, ExecSpace, IndexType>;

  using SpacePoint = typename GridType::SpacePoint;
  using SpatialBoundingBox = typename GridType::SpatialBoundingBox;

private:
  constexpr static bool DeviceExec = axom::execution_space<ExecSpace>::onDevice();

public:
  /*!
   * Constructor for PointFinder
   *
   * \param meshWrapper A non-null MeshWrapperType
   * \param res The grid resolution for the spatial acceleration structure
   * \param bboxScaleFactor A number slightly larger than 1 by which to expand
   * cell bounding boxes
   *
   * \sa constructors in PointInCell class for more details about parameters
   */
  PointFinder(const MeshWrapperType* meshWrapper, const int* res, double bboxScaleFactor, int allocatorID)
    : m_meshWrapper(meshWrapper)
    , m_allocatorID(allocatorID)
  {
    SLIC_ASSERT(m_meshWrapper != nullptr);
    SLIC_ASSERT(bboxScaleFactor >= 1.);

    const axom::IndexType numCells = m_meshWrapper->numElements();

    // setup bounding boxes -- Slightly scaled for robustness

    SpatialBoundingBox meshBBox;
#ifdef AXOM_USE_UMPIRE
    axom::Array<SpatialBoundingBox, 1, MemorySpace::Host> cellBBoxesHost(numCells);
#else
    axom::Array<SpatialBoundingBox> cellBBoxesHost(numCells);
#endif
    m_meshWrapper->template computeBoundingBoxes<NDIMS>(bboxScaleFactor,
                                                        cellBBoxesHost.data(),
                                                        meshBBox);
    if(DeviceExec)
    {
      // Copy the host-side bounding boxes to GPU memory.
      m_cellBBoxes = axom::Array<SpatialBoundingBox>(cellBBoxesHost, m_allocatorID);
    }
    else
    {
      m_cellBBoxes = std::move(cellBBoxesHost);
    }

    // initialize implicit grid, handle case where resolution is a NULL pointer
    if(res != nullptr)
    {
      using GridResolution = axom::primal::Point<int, NDIMS>;
      GridResolution gridRes(res);
      m_grid.initialize(meshBBox, &gridRes, numCells, allocatorID);
    }
    else
    {
      m_grid.initialize(meshBBox, nullptr, numCells, allocatorID);
    }

    // add mesh elements to grid
    m_grid.insert(numCells, m_cellBBoxes.data());
  }

  /*!
   * Query to find the mesh cell containing query point with coordinates \a pos
   *
   * \sa PointInCell::locatePoint() for more details about parameters
   */
  IndexType locatePoint(const double* pos, double* isoparametric) const
  {
    IndexType containingCell = PointInCellTraits<mesh_tag>::NO_CELL;

    SLIC_ASSERT(pos != nullptr);
    SpacePoint pt(pos);
    SpacePoint isopar;

    if(DeviceExec)
    {
      axom::Array<SpacePoint> dev_ptr(axom::ArrayView<const SpacePoint>(&pt, 1), m_allocatorID);
      locatePoints(dev_ptr, &containingCell, &isopar);
    }
    else
    {
      locatePoints(axom::ArrayView<const SpacePoint>(&pt, 1), &containingCell, &isopar);
    }

    // Copy data back to input parameter isoparametric, if necessary
    if(isoparametric != nullptr)
    {
      isopar.array().to_array(isoparametric);
    }

    return containingCell;
  }

  void locatePoints(axom::ArrayView<const SpacePoint> pts,
                    IndexType* outCellIds,
                    SpacePoint* outIsoparametricCoords) const
  {
    using IndexArray = axom::Array<IndexType>;

#ifdef AXOM_USE_RAJA
    using IndexView = axom::ArrayView<IndexType>;
  #ifdef AXOM_USE_UMPIRE
    using HostIndexArray = axom::Array<IndexType, 1, axom::MemorySpace::Host>;
    using HostPointArray = axom::Array<SpacePoint, 1, axom::MemorySpace::Host>;

    using HostIndexView = axom::ArrayView<IndexType, 1, axom::MemorySpace::Host>;
    using HostPointView = axom::ArrayView<SpacePoint, 1, axom::MemorySpace::Host>;
    using ConstHostPointView = axom::ArrayView<const SpacePoint, 1, axom::MemorySpace::Host>;
  #else
    using HostIndexArray = IndexArray;
    using HostPointArray = axom::Array<SpacePoint>;

    using HostIndexView = IndexView;
    using HostPointView = axom::ArrayView<SpacePoint>;
    using ConstHostPointView = axom::ArrayView<const SpacePoint>;
  #endif  // AXOM_USE_UMPIRE
#endif    // AXOM_USE_RAJA

    auto gridQuery = m_grid.getQueryObject();

    axom::IndexType npts = pts.size();

    IndexArray offsets(npts, npts, m_allocatorID);
    IndexArray counts(npts, npts, m_allocatorID);

#ifdef AXOM_USE_RAJA
    IndexView countsPtr = counts;

    axom::ReduceSum<ExecSpace, IndexType> totalCountReduce(0);
    // Step 1: count number of candidate intersections for each point
    for_all<ExecSpace>(
      npts,
      AXOM_LAMBDA(IndexType i) {
        countsPtr[i] = gridQuery.countCandidates(pts[i]);
        totalCountReduce += countsPtr[i];
      });

    // Step 2: exclusive scan for offsets in candidate array
    axom::exclusive_scan<ExecSpace>(counts, offsets);

    axom::IndexType totalCount = totalCountReduce.get();

    // Step 3: allocate memory for all candidates
    IndexArray candidates(totalCount, totalCount, m_allocatorID);
    IndexView candidatesPtr = candidates;
    IndexView offsetsPtr = offsets;
    const SpatialBoundingBox* cellBBoxes = m_cellBBoxes.data();

    // Step 4: fill candidate array for each query box
    for_all<ExecSpace>(
      npts,
      AXOM_LAMBDA(IndexType i) {
        int startIdx = offsetsPtr[i];
        int currCount = 0;
        const SpacePoint& pt = pts[i];
        auto onCandidate = [&](int candidateIdx) -> bool {
          // Check that point is in bounding box of candidate element
          if(cellBBoxes[candidateIdx].contains(pt))
          {
            candidatesPtr[startIdx] = candidateIdx;
            currCount++;
            startIdx++;
          }
          return currCount >= countsPtr[i];
        };
        gridQuery.visitCandidates(pts[i], onCandidate);
        countsPtr[i] = currCount;
      });

    // Temporary host arrays we copy device-side data into when the candidate
    // search is conducted on the GPU
    HostPointArray ptsHost, outIsoparHost;
    HostIndexArray outCellIdsHost;
    HostIndexArray candidatesHost, offsetsHost, countsHost;

    // For sequential/OpenMP execution, just use the argument pointers
    // directly.
    HostIndexView outCellIdsPtr(outCellIds, pts.size());
    HostPointView outIsoparPtr(outIsoparametricCoords, pts.size());

    // If the candidate search takes place on the GPU, we need to copy the
    // device-side data first, then set these array views to point to the
    // intermediate arrays. Otherwise, we can set these to point to the result
    // arrays directly.
    ConstHostPointView ptsHostPtr;
    HostIndexView candidatesHostPtr, offsetsHostPtr, countsHostPtr;

    if(DeviceExec)
    {
      // Copy points and candidate intersections to host memory.
      ptsHost = pts;
      candidatesHost = candidates;
      offsetsHost = offsets;
      countsHost = counts;
      // Set up views from intermediate host arrays
      ptsHostPtr = ptsHost;
      candidatesHostPtr = candidatesHost;
      offsetsHostPtr = offsetsHost;
      countsHostPtr = countsHost;
      // Allocate intermediate output buffers on the host side.
      outCellIdsHost.resize(pts.size());
      if(outIsoparametricCoords)
      {
        outIsoparHost.resize(pts.size());
      }
      outCellIdsPtr = outCellIdsHost;
      outIsoparPtr = outIsoparHost;
    }
    else
    {
      ptsHostPtr = pts;
      candidatesHostPtr = candidates;
      offsetsHostPtr = offsets;
      countsHostPtr = counts;
    }

    // Step 5: Check each candidate
    // TODO: This only supports sequential execution right now, because we
    // don't build MFEM in a thread-safe manner.
    const MeshWrapperType* meshWrapperPtr = m_meshWrapper;
    for_all<SEQ_EXEC>(
      npts,
      AXOM_HOST_LAMBDA(IndexType i) {
        outCellIdsPtr[i] = PointInCellTraits<mesh_tag>::NO_CELL;
        const SpacePoint& pt = ptsHostPtr[i];
        SpacePoint isopar;
        for(int icell = 0; icell < countsHostPtr[i]; icell++)
        {
          const int cellIdx = candidatesHostPtr[icell + offsetsHostPtr[i]];
          // if isopar is in the proper range
          if(meshWrapperPtr->locatePointInCell(cellIdx, pt.data(), isopar.data()))
          {
            // then we have found the cellID
            outCellIdsPtr[i] = cellIdx;
            break;
          }
        }
        if(outIsoparametricCoords != nullptr)
        {
          outIsoparPtr[i] = isopar;
        }
      });

    if(DeviceExec)
    {
      // Copy back to GPU memory.
      axom::copy(outCellIds, outCellIdsHost.data(), outCellIdsHost.size() * sizeof(IndexType));
      axom::copy(outIsoparametricCoords,
                 outIsoparHost.data(),
                 outIsoparHost.size() * sizeof(SpacePoint));
    }
#else   // AXOM_USE_RAJA
    for(int i = 0; i < npts; i++)
    {
      const SpacePoint& pt = pts[i];
      SpacePoint isopar;
      outCellIds[i] = PointInCellTraits<mesh_tag>::NO_CELL;
      gridQuery.visitCandidates(pt, [&](int candidateIdx) -> bool {
        if(m_cellBBoxes[candidateIdx].contains(pt))
        {
          if(m_meshWrapper->locatePointInCell(candidateIdx, pt.data(), isopar.data()))
          {
            outCellIds[i] = candidateIdx;
            return true;
          }
        }
        return false;
      });
      if(outIsoparametricCoords != nullptr)
      {
        outIsoparametricCoords[i] = isopar;
      }
    }
#endif  // AXOM_USE_RAJA
  }

  /*! Returns a const reference to the given cells's bounding box */
  const SpatialBoundingBox& cellBoundingBox(IndexType cellIdx) const
  {
    return m_cellBBoxes[cellIdx];
  }

  std::vector<IndexType> getCandidates(SpacePoint const& pt) const
  {
    return m_grid.getCandidatesAsArray(pt);
  }

private:
  GridType m_grid;
  const MeshWrapperType* m_meshWrapper;
  axom::Array<SpatialBoundingBox> m_cellBBoxes;
  int m_allocatorID;
};

}  // end namespace detail
}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_POINT_IN_CELL_POINT_FINDER_HPP_
