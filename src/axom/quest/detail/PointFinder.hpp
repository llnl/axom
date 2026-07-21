// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/core/memory_management.hpp"
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
  PointFinder(const MeshWrapperType* meshWrapper,
              const int* res,
              double bboxScaleFactor,
              int allocatorID,
              HostAllocator hostAllocator)
    : m_meshWrapper(meshWrapper)
    , m_allocatorID(allocatorID)
    , m_hostAllocator(hostAllocator)
  {
    SLIC_ASSERT(m_meshWrapper != nullptr);
    SLIC_ASSERT(bboxScaleFactor >= 1.);

    const axom::IndexType numCells = m_meshWrapper->numElements();

    // setup bounding boxes -- Slightly scaled for robustness

    SpatialBoundingBox meshBBox;
    axom::Array<SpatialBoundingBox> cellBBoxesHost(numCells,
                                                   numCells,
                                                   m_hostAllocator.getID(),
                                                   m_hostAllocator);
    m_meshWrapper->template computeBoundingBoxes<NDIMS>(bboxScaleFactor,
                                                        cellBBoxesHost.data(),
                                                        meshBBox);
    if(DeviceExec)
    {
      // Copy the host-side bounding boxes to GPU memory.
      m_cellBBoxes = axom::Array<SpatialBoundingBox>(cellBBoxesHost, m_allocatorID, m_hostAllocator);
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
      axom::Array<SpacePoint> dev_ptr(axom::ArrayView<const SpacePoint>(&pt, 1),
                                      m_allocatorID,
                                      m_hostAllocator);
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
    using HostIndexArray = axom::Array<IndexType>;
    using HostPointArray = axom::Array<SpacePoint>;

    using HostIndexView = axom::ArrayView<IndexType>;
    using HostPointView = axom::ArrayView<SpacePoint>;
    using ConstHostPointView = axom::ArrayView<const SpacePoint>;
#endif  // AXOM_USE_RAJA

    auto gridQuery = m_grid.getQueryObject();

    axom::IndexType npts = pts.size();

    IndexArray offsets(npts, npts, m_allocatorID, m_hostAllocator);
    IndexArray counts(npts, npts, m_allocatorID, m_hostAllocator);

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
    IndexArray candidates(totalCount, totalCount, m_allocatorID, m_hostAllocator);
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

    auto locateCandidatesHost = [&](ConstHostPointView ptsHostPtr,
                                    HostIndexView candidatesHostPtr,
                                    HostIndexView offsetsHostPtr,
                                    HostIndexView countsHostPtr,
                                    HostIndexView outCellIdsPtr,
                                    HostPointView outIsoparPtr) {
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
            if(meshWrapperPtr->locatePointInCell(cellIdx, pt.data(), isopar.data()))
            {
              outCellIdsPtr[i] = cellIdx;
              break;
            }
          }
          if(outIsoparametricCoords != nullptr)
          {
            outIsoparPtr[i] = isopar;
          }
        });
    };

    if(DeviceExec)
    {
      HostPointArray ptsHost(pts, m_hostAllocator.getID(), m_hostAllocator);
      HostIndexArray candidatesHost(candidates, m_hostAllocator.getID(), m_hostAllocator);
      HostIndexArray offsetsHost(offsets, m_hostAllocator.getID(), m_hostAllocator);
      HostIndexArray countsHost(counts, m_hostAllocator.getID(), m_hostAllocator);
      HostIndexArray outCellIdsHost(npts, npts, m_hostAllocator.getID(), m_hostAllocator);
      HostPointArray outIsoparHost(0, 0, m_hostAllocator.getID(), m_hostAllocator);
      if(outIsoparametricCoords != nullptr)
      {
        outIsoparHost.resize(npts);
      }

      locateCandidatesHost(ptsHost, candidatesHost, offsetsHost, countsHost, outCellIdsHost, outIsoparHost);

      // Copy back to GPU memory.
      axom::copy(outCellIds, outCellIdsHost.data(), outCellIdsHost.size() * sizeof(IndexType));
      if(outIsoparametricCoords != nullptr)
      {
        axom::copy(outIsoparametricCoords,
                   outIsoparHost.data(),
                   outIsoparHost.size() * sizeof(SpacePoint));
      }
    }
    else
    {
      locateCandidatesHost(pts,
                           candidates,
                           offsets,
                           counts,
                           HostIndexView(outCellIds, pts.size()),
                           HostPointView(outIsoparametricCoords, pts.size()));
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
  HostAllocator m_hostAllocator;
};

}  // end namespace detail
}  // end namespace quest
}  // end namespace axom
