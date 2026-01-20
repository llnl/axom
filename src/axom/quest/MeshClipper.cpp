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
  , m_screenLevel(3)
{
  // Initialize statistics used by this class.
  m_counterStats["cellsIn"].set_int64(0);
  m_counterStats["cellsOn"].set_int64(0);
  m_counterStats["cellsOut"].set_int64(0);
  m_counterStats["tetsIn"].set_int64(0);
  m_counterStats["tetsOn"].set_int64(0);
  m_counterStats["tetsOut"].set_int64(0);
}

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
  [[maybe_unused]] const axom::IndexType cellCount = m_shapeMesh.getCellCount();
  SLIC_ASSERT(ovlap.size() == cellCount);

  auto& cellsInCount = *m_counterStats["cellsIn"].as_int64_ptr();
  auto& cellsOnCount = *m_counterStats["cellsOn"].as_int64_ptr();
  auto& cellsOutCount = *m_counterStats["cellsOut"].as_int64_ptr();
  auto& tetsInCount = *m_counterStats["tetsIn"].as_int64_ptr();
  auto& tetsOnCount = *m_counterStats["tetsOn"].as_int64_ptr();
  auto& tetsOutCount = *m_counterStats["tetsOut"].as_int64_ptr();

  // Try to label cells as inside, outside or on shape boundary
  axom::Array<LabelType> cellLabels;
  bool withCellInOut = false;
  if(m_screenLevel >= 1)
  {
    AXOM_ANNOTATE_BEGIN("MeshClipper:label_cells");
    withCellInOut = m_strategy->labelCellsInOut(m_shapeMesh, cellLabels);
    AXOM_ANNOTATE_END("MeshClipper:label_cells");
  }

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
      getLabelCounts(cellLabels, cellsInCount, cellsOnCount, cellsOutCount);
      logClippingStats();
    }

    AXOM_ANNOTATE_BEGIN("MeshClipper::processInOut");

    m_impl->initVolumeOverlaps(cellLabels.view(), ovlap);

    axom::Array<axom::IndexType> cellsOnBdry;
    m_impl->collectOnIndices(cellLabels.view(), cellsOnBdry);

    axom::Array<LabelType> tetLabels;
    bool withTetInOut = false;
    if(m_screenLevel >= 2)
    {
      AXOM_ANNOTATE_BEGIN("MeshClipper:label_tets");
      withTetInOut = m_strategy->labelTetsInOut(m_shapeMesh, cellsOnBdry.view(), tetLabels);
      AXOM_ANNOTATE_END("MeshClipper:label_tets");
    }

    axom::Array<axom::IndexType> tetsOnBdry;

    if(withTetInOut)
    {
      if(m_verbose)
      {
        getLabelCounts(tetLabels, tetsInCount, tetsOnCount, tetsOutCount);
        logClippingStats();
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
      done = m_strategy->specializedClipTets(m_shapeMesh, ovlap, tetsOnBdry, m_counterStats);
    }
    else
    {
      done = m_strategy->specializedClipCells(m_shapeMesh, ovlap, cellsOnBdry, m_counterStats);
    }

    if(!done)
    {
      AXOM_ANNOTATE_SCOPE("MeshClipper::clip3D_limited");
      if(withTetInOut)
      {
        m_impl->computeClipVolumes3DTets(tetsOnBdry.view(), ovlap, m_counterStats);
      }
      else
      {
        m_impl->computeClipVolumes3D(cellsOnBdry.view(), ovlap, m_counterStats);
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
    AXOM_ANNOTATE_BEGIN("MeshClipper:specialized_clip");
    done = m_strategy->specializedClipCells(m_shapeMesh, ovlap, m_counterStats);
    AXOM_ANNOTATE_END("MeshClipper:specialized_clip");

    if(!done)
    {
      AXOM_ANNOTATE_SCOPE("MeshClipper:clip_fcn");
      m_impl->computeClipVolumes3D(ovlap, m_counterStats);
    }
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

#if defined(AXOM_USE_MPI)
template <typename T>
void globalReduce(axom::Array<T>& values, int reduceOp)
{
  axom::Array<T> localValues(values);
  MPI_Allreduce(localValues.data(),
                values.data(),
                values.size(),
                axom::mpi_traits<T>::type,
                reduceOp,
                MPI_COMM_WORLD);
}
#endif

void MeshClipper::accumulateClippingStats(conduit::Node& curStats, const conduit::Node& newStats)
{
  for(int i = 0; i < newStats.number_of_children(); ++i)
  {
    const auto& newStat = newStats.child(i);
    SLIC_ERROR_IF(!newStat.dtype().is_integer(),
                  "MeshClipper statistic must be integer"
                  " (at least until a need for floats arises).");
    auto& currentStat = curStats[newStat.name()];
    if(currentStat.dtype().is_empty())
    {
      currentStat.set_int64(newStat.as_int64());
    }
    else
    {
      *currentStat.as_int64_ptr() += newStat.as_int64();
    }
  }
}

conduit::Node MeshClipper::getGlobalClippingStats() const
{
  conduit::Node stats;
  auto& locNode = stats["loc"];
  auto& maxNode = stats["max"];
  auto& sumNode = stats["sum"];

  locNode.set(m_counterStats);
  sumNode.set(m_counterStats);
  maxNode.set(m_counterStats);

#if defined(AXOM_USE_MPI)
  // Do sum and max reductions.
  axom::Array<std::int64_t> sums(0, sumNode.number_of_children());
  for(int i = 0; i < sumNode.number_of_children(); ++i)
  {
    sums.push_back(locNode.child(i).as_int64());
  }
  axom::Array<std::int64_t> maxs(sums);
  globalReduce(maxs, MPI_MAX);
  globalReduce(sums, MPI_SUM);

  for(int i = 0; i < sumNode.number_of_children(); ++i)
  {
    *maxNode.child(i).as_int64_ptr() = maxs[i];
    *sumNode.child(i).as_int64_ptr() = sums[i];
  }
#endif

  return stats;
}

void MeshClipper::logClippingStats(bool local, bool sum, bool max) const
{
  conduit::Node stats = getGlobalClippingStats();
  if(local)
  {
    SLIC_INFO(std::string("MeshClipper loc-stats: ") +
              stats["loc"].to_string("yaml", 2, 0, "", " "));
  }
  if(sum)
  {
    SLIC_INFO(std::string("MeshClipper sum-stats: ") +
              stats["sum"].to_string("yaml", 2, 0, "", " "));
  }
  if(max)
  {
    SLIC_INFO(std::string("MeshClipper max-stats: ") +
              stats["max"].to_string("yaml", 2, 0, "", " "));
  }
}

}  // namespace experimental
}  // end namespace quest
}  // end namespace axom
