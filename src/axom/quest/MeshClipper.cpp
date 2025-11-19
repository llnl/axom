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
      axom::fmt::format("MeshClipperStrategy '{}' did not return the correct cell label array size of {}",
                        m_strategy->name(),
                        m_shapeMesh.getCellCount()));
    SLIC_ERROR_IF(cellLabels.getAllocatorID() != allocId,
                  axom::fmt::format("MeshClipperStrategy '{}' failed to provide cellLabels data "
                                    "with the required allocator id {}",
                                    m_strategy->name(),
                                    allocId));

    // Counting labels is non-essential and presumed to be relatively very fast.
    getLabelCounts(cellLabels, m_cellsInCount, m_cellsOnCount, m_cellsOutCount);
    if(m_verbose) { logClippingStats(); }

    AXOM_ANNOTATE_BEGIN("MeshClipper::processInOut");

    m_impl->initVolumeOverlaps(cellLabels.view(), ovlap);

    axom::Array<axom::IndexType> cellsOnBdry;
    m_impl->collectOnIndices(cellLabels.view(), cellsOnBdry);
    SLIC_ASSERT(cellsOnBdry.size() == m_cellsOnCount);

    axom::Array<LabelType> tetLabels;
    bool withTetInOut = m_strategy->labelTetsInOut(m_shapeMesh, cellsOnBdry.view(), tetLabels);

    axom::Array<axom::IndexType> tetsOnBdry;

    if(withTetInOut)
    {
      SLIC_ERROR_IF(
        tetLabels.size() != NUM_TETS_PER_HEX*cellsOnBdry.size(),
        axom::fmt::format("MeshClipperStrategy '{}' did not return the correct tet label array size of {}",
                          m_strategy->name(),
                          NUM_TETS_PER_HEX*cellsOnBdry.size()));
      SLIC_ERROR_IF(tetLabels.getAllocatorID() != allocId,
                    axom::fmt::format("MeshClipperStrategy '{}' failed to provide tetLabels data "
                                      "with the required allocator id {}",
                                      m_strategy->name(),
                                      allocId));

      // Counting labels is non-essential and presumed to be very fast.
      getLabelCounts(tetLabels, m_tetsInCount, m_tetsOnCount, m_tetsOutCount);
      if(m_verbose) { logClippingStats(); }

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
      if(done)
      {
        m_tetsClipped = tetsOnBdry.size();
      }
    }
    else
    {
      done = m_strategy->specializedClipCells(m_shapeMesh, ovlap, cellsOnBdry);
      if(done)
      {
        m_hexesClipped = cellsOnBdry.size();
      }
    }

    if(!done)
    {
      AXOM_ANNOTATE_SCOPE("MeshClipper::clip3D_limited");
      if(withTetInOut)
      {
        m_impl->computeClipVolumes3DTets(tetsOnBdry.view(), ovlap);
        m_tetsClipped = tetsOnBdry.size();
      }
      else
      {
        m_impl->computeClipVolumes3D(cellsOnBdry.view(), ovlap);
        m_hexesClipped = cellsOnBdry.size();
      }
    }
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

    m_cellsOnCount = cellCount;
    m_hexesClipped = cellCount;
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

template<typename T>
void globalReduce(axom::Array<T>& values, int reduceOp)
{
#ifdef AXOM_USE_MPI
  axom::Array<T> localValues(values);
  MPI_Allreduce(localValues.data(),
             values.data(),
             values.size(),
             axom::mpi_traits<T>::type,
             reduceOp,
             MPI_COMM_WORLD);
#endif
}

conduit::Node MeshClipper::getClippingStats() const
{
  axom::Array<std::int64_t> sums {
    m_cellsInCount, m_cellsOnCount, m_cellsOutCount, m_hexesClipped,
    m_tetsInCount, m_tetsOnCount, m_tetsOutCount, m_tetsClipped};
  axom::Array<std::int64_t> maxs {
    m_cellsInCount, m_cellsOnCount, m_cellsOutCount, m_hexesClipped,
    m_tetsInCount, m_tetsOnCount, m_tetsOutCount, m_tetsClipped};

  conduit::Node stats;

  auto& locNode = stats["loc"];
  locNode["cellsIn"].set(sums[0]);
  locNode["cellsOn"].set(sums[1]);
  locNode["cellsOut"].set(sums[2]);
  locNode["hexesClipped"].set(sums[3]);
  locNode["tetsIn"].set(sums[4]);
  locNode["tetsOn"].set(sums[5]);
  locNode["tetsOut"].set(sums[6]);
  locNode["tetsClipped"].set(sums[7]);

  auto& sumNode = stats["sum"];
  globalReduce(sums, MPI_SUM);
  sumNode["cellsIn"].set(sums[0]);
  sumNode["cellsOn"].set(sums[1]);
  sumNode["cellsOut"].set(sums[2]);
  sumNode["hexesClipped"].set(sums[3]);
  sumNode["tetsIn"].set(sums[4]);
  sumNode["tetsOn"].set(sums[5]);
  sumNode["tetsOut"].set(sums[6]);
  sumNode["tetsClipped"].set(sums[7]);

  auto& maxNode = stats["max"];
  globalReduce(maxs, MPI_MAX);
  maxNode["cellsIn"].set(maxs[0]);
  maxNode["cellsOn"].set(maxs[1]);
  maxNode["cellsOut"].set(maxs[2]);
  maxNode["hexesClipped"].set(maxs[3]);
  maxNode["tetsIn"].set(maxs[4]);
  maxNode["tetsOn"].set(maxs[5]);
  maxNode["tetsOut"].set(maxs[6]);
  maxNode["tetsClipped"].set(maxs[7]);

  return stats;
}

void MeshClipper::logClippingStats() const
{
  conduit::Node stats = getClippingStats();
  SLIC_INFO(std::string("MeshClipper loc-stats: ") + stats["loc"].to_string("yaml", 2, 0, "", " "));
  SLIC_INFO(std::string("MeshClipper sum-stats: ") + stats["sum"].to_string("yaml", 2, 0, "", " "));
  SLIC_INFO(std::string("MeshClipper max-stats: ") + stats["max"].to_string("yaml", 2, 0, "", " "));
}

}  // namespace experimental
}  // end namespace quest
}  // end namespace axom
