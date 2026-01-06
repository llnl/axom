// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/quest/MeshClipper.hpp"
#include "axom/quest/detail/clipping/MeshClipperImpl.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/runtime_policy.hpp"
#include "axom/slic/interface/slic_macros.hpp"
#include "axom/fmt.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{

MeshClipper::MeshClipper(quest::experimental::ShapeMesh& shapeMesh,
                         const std::shared_ptr<quest::experimental::MeshClipperStrategy>& strategy)
  : m_shapeMesh(shapeMesh)
  , m_strategy(strategy)
  , m_impl(newImpl())
  , m_verbose(false)
{ }

void MeshClipper::clip(axom::Array<double>& ovlap)
{
  const int allocId = m_shapeMesh.getAllocatorID();
  const axom::IndexType cellCount = m_shapeMesh.getCellCount();

  // Resize output array and use appropriate allocator.
  if(ovlap.size() < cellCount || ovlap.getAllocatorID() != allocId)
  {
    AXOM_ANNOTATE_SCOPE("MeshClipper::clip_alloc");
    ovlap = axom::Array<double>(ArrayOptions::Uninitialized(), cellCount, cellCount, allocId);
  }
  clip(ovlap.view());
}

/**
 * @brief Orchestrates the geometry clipping by using the capabilities of the
 * MeshClipperStrategy implementation.
 *
 * If the strategy can label cells as inside/on/outside geometry
 * boundary, use that to reduce reliance on expensive clipping methods.
 *
 * Regardless of labeling, try to use specialized clipping first.
 * If specialized methods aren't implemented, resort to discretizing
 * geometry into tets or octs for clipping against mesh cells.
 */
void MeshClipper::clip(axom::ArrayView<double> ovlap)
{
  const int allocId = m_shapeMesh.getAllocatorID();
  const axom::IndexType cellCount = m_shapeMesh.getCellCount();
  SLIC_ASSERT(ovlap.size() == cellCount);

  // Try to label cells as inside, outside or on shape boundary
  axom::Array<LabelType> cellLabels;
  bool withCellInOut = m_strategy->labelCellsInOut(m_shapeMesh, cellLabels);

  bool done = false;

  if(withCellInOut)
  {
    SLIC_ERROR_IF(
      cellLabels.size() != m_shapeMesh.getCellCount(),
      axom::fmt::format("MeshClipperStrategy '{}' did not return the correct array size of {}",
                        m_strategy->name(),
                        m_shapeMesh.getCellCount()));
    SLIC_ERROR_IF(cellLabels.getAllocatorID() != allocId,
                  axom::fmt::format("MeshClipperStrategy '{}' failed to provide cellLabels data "
                                    "with the required allocator id {}",
                                    m_strategy->name(),
                                    allocId));

    if(m_verbose)
    {
      logLabelStats(cellLabels, "cells");
    }

    AXOM_ANNOTATE_BEGIN("MeshClipper::processInOut");

    m_impl->initVolumeOverlaps(cellLabels.view(), ovlap);

    axom::Array<axom::IndexType> cellsOnBdry;
    m_impl->collectOnIndices(cellLabels.view(), cellsOnBdry);

    axom::Array<LabelType> tetLabels;
    bool withTetInOut = m_strategy->labelTetsInOut(m_shapeMesh, cellsOnBdry.view(), tetLabels);

    axom::Array<axom::IndexType> tetsOnBdry;

    if(withTetInOut)
    {
      if(m_verbose)
      {
        logLabelStats(tetLabels, "tets");
      }
      m_impl->collectOnIndices(tetLabels.view(), tetsOnBdry);
      m_impl->remapTetIndices(cellsOnBdry, tetsOnBdry);

      SLIC_ASSERT(tetsOnBdry.getAllocatorID() == m_shapeMesh.getAllocatorID());
      SLIC_ASSERT(tetsOnBdry.size() <= cellsOnBdry.size() * NUM_TETS_PER_HEX);

      m_impl->addVolumesOfInteriorTets(cellsOnBdry.view(), tetLabels.view(), ovlap);
    }

    AXOM_ANNOTATE_END("MeshClipper::processInOut");

    //
    // If implementation has a specialized clip, use it.
    //
    if(withTetInOut)
    {
      done = m_strategy->specializedClipTets(m_shapeMesh, ovlap, tetsOnBdry);
    }
    else
    {
      done = m_strategy->specializedClipCells(m_shapeMesh, ovlap, cellsOnBdry);
    }

    if(!done)
    {
      AXOM_ANNOTATE_SCOPE("MeshClipper::clip3D_limited");
      if(withTetInOut)
      {
        m_impl->computeClipVolumes3DTets(tetsOnBdry.view(), ovlap);
      }
      else
      {
        m_impl->computeClipVolumes3D(cellsOnBdry.view(), ovlap);
      }
    }

    m_localCellInCount = cellsOnBdry.size();
  }
  else  // !withCellInOut
  {
    std::string msg =
      axom::fmt::format("MeshClipper strategy '{}' did not provide in/out cell labels.\n",
                        m_strategy->name());
    SLIC_INFO(msg);
    m_impl->zeroVolumeOverlaps(ovlap);
    done = m_strategy->specializedClipCells(m_shapeMesh, ovlap);

    if(!done)
    {
      AXOM_ANNOTATE_SCOPE("MeshClipper::clip3D");
      m_impl->computeClipVolumes3D(ovlap);
    }

    m_localCellInCount = cellCount;
  }
}

/*!
 * @brief Allocate an Impl for the execution-space computations
 * of this clipper.
 */
std::unique_ptr<MeshClipper::Impl> MeshClipper::newImpl()
{
  using RuntimePolicy = axom::runtime_policy::Policy;

  auto runtimePolicy = m_shapeMesh.getRuntimePolicy();

  std::unique_ptr<Impl> impl;
  if(runtimePolicy == RuntimePolicy::seq)
  {
    impl.reset(new detail::MeshClipperImpl<axom::SEQ_EXEC>(*this));
  }
#ifdef AXOM_RUNTIME_POLICY_USE_OPENMP
  else if(runtimePolicy == RuntimePolicy::omp)
  {
    impl.reset(new detail::MeshClipperImpl<axom::OMP_EXEC>(*this));
  }
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_CUDA
  else if(runtimePolicy == RuntimePolicy::cuda)
  {
    impl.reset(new detail::MeshClipperImpl<axom::CUDA_EXEC<256>>(*this));
  }
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_HIP
  else if(runtimePolicy == RuntimePolicy::hip)
  {
    impl.reset(new detail::MeshClipperImpl<axom::HIP_EXEC<256>>(*this));
  }
#endif
  else
  {
    SLIC_ERROR(axom::fmt::format("MeshClipper has no impl for runtime policy {}", runtimePolicy));
  }
  return impl;
}

void MeshClipper::logLabelStats(axom::ArrayView<const LabelType> labels, const std::string& labelType)
{
  axom::IndexType countsa[4];
  axom::IndexType countsb[4];
  getLabelCounts(labels, countsa[0], countsa[1], countsa[2]);
  countsa[3] = m_shapeMesh.getCellCount();
#ifdef AXOM_USE_MPI
  MPI_Reduce(countsa, countsb, 4, axom::mpi_traits<axom::IndexType>::type, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
  std::string msg = axom::fmt::format(
    "MeshClipper strategy '{}' globally labeled {} {} inside, {} on and {} outside, for mesh with "
    "{} cells ({} tets)\n",
    m_strategy->name(),
    labelType,
    countsb[0],
    countsb[1],
    countsb[2],
    countsb[3],
    countsb[3] * NUM_TETS_PER_HEX);
  SLIC_INFO(msg);
}

void MeshClipper::getClippingStats(axom::IndexType& localCellInCount,
                                   axom::IndexType& globalCellInCount,
                                   axom::IndexType& maxLocalCellInCount) const
{
  localCellInCount = m_localCellInCount;
#ifdef AXOM_USE_MPI
  MPI_Reduce(&localCellInCount,
             &globalCellInCount,
             1,
             axom::mpi_traits<axom::IndexType>::type,
             MPI_SUM,
             0,
             MPI_COMM_WORLD);
  MPI_Reduce(&localCellInCount,
             &maxLocalCellInCount,
             1,
             axom::mpi_traits<axom::IndexType>::type,
             MPI_MAX,
             0,
             MPI_COMM_WORLD);
#else
  maxLocalCellInCount = localCellInCount;
  globalCellInCount = localCellInCount;
#endif
}

}  // namespace experimental
}  // end namespace quest
}  // end namespace axom
