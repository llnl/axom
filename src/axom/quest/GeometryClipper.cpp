// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/quest/GeometryClipper.hpp"
#include "axom/quest/detail/GeometryClipperDelegateExec.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/runtime_policy.hpp"
#include "axom/slic/interface/slic_macros.hpp"
#include "axom/fmt.hpp"

namespace axom
{
namespace quest
{

GeometryClipper::GeometryClipper(quest::ShapeeMesh& shapeeMesh,
                                 const std::shared_ptr<quest::GeometryClipperStrategy>& strategy)
  : m_shapeeMesh(shapeeMesh)
  , m_strategy(strategy)
  , m_delegate(newDelegate())
  , m_verbose(false)
{ }

/*
  If the strategy can label cells as inside/on/outside geometry
  boundary, use that to reduce the use of expensive clipping methods.

  Regardless of labling, try to use specialized clipping first.
  If those methods aren't implemented, resort to discretizing
  geomety into tets or octs for brute-force clipping.
*/
void GeometryClipper::clip(axom::Array<double>& ovlap)
{
  const int allocId = m_shapeeMesh.getAllocatorID();
  const axom::IndexType cellCount = m_shapeeMesh.getCellCount();

  // Resize output array and use appropriate allocator.
  if(ovlap.size() < cellCount || ovlap.getAllocatorID() != allocId)
  {
    AXOM_ANNOTATE_SCOPE("GeometryClipper::clip_alloc");
    ovlap = axom::Array<double>(ArrayOptions::Uninitialized(), cellCount, cellCount, allocId);
  }
  clip(ovlap.view());
}

void GeometryClipper::clip(axom::ArrayView<double> ovlap)
{
  const int allocId = m_shapeeMesh.getAllocatorID();
  const axom::IndexType cellCount = m_shapeeMesh.getCellCount();
  SLIC_ASSERT(ovlap.size() == cellCount);

  // Try to label cells as inside, outside or on shape boundary
  axom::Array<char> labels;
  bool withInOut = m_strategy->labelInOut(m_shapeeMesh, labels);

  bool done = false;

  if(withInOut)
  {
    SLIC_ERROR_IF(
      labels.size() != m_shapeeMesh.getCellCount(),
      axom::fmt::format("GeometryClipperStrategy '{}' did not return the correct array size of {}",
                        m_strategy->name(),
                        m_shapeeMesh.getCellCount()));
    SLIC_ERROR_IF(labels.getAllocatorID() != allocId,
                  axom::fmt::format("GeometryClipperStrategy '{}' failed to provide labels data "
                                    "with the required allocator id {}",
                                    m_strategy->name(),
                                    allocId));

    if(m_verbose)
    {
      axom::IndexType inCount;
      axom::IndexType onCount;
      axom::IndexType outCount;
      getLabelCounts(labels.view(), inCount, onCount, outCount);
      std::string msg = axom::fmt::format(
        "GeometryClipper with strategy '{}' labeled cells {} inside, {} on and {} outside, out of {} "
        "cells",
        m_strategy->name(),
        inCount,
        onCount,
        outCount,
        m_shapeeMesh.getCellCount());
      SLIC_INFO(msg);
    }

    AXOM_ANNOTATE_BEGIN("GeometryClipper::processInOut");
    m_delegate->initVolumeOverlaps(labels.view(), ovlap);

    axom::Array<axom::IndexType> unlabeledCells;
    m_delegate->collectUnlabeledCellIndices(labels.view(), unlabeledCells);
    AXOM_ANNOTATE_END("GeometryClipper::processInOut");

    done = m_strategy->specializedClip(m_shapeeMesh, ovlap, unlabeledCells);

    if(!done)
    {
      AXOM_ANNOTATE_SCOPE("GeometryClipper::clip3D_limited");
      m_delegate->computeClipVolumes3D(unlabeledCells.view(), ovlap);
    }
  }
  else  // !withInOut
  {
    done = m_strategy->specializedClip(m_shapeeMesh, ovlap);

    if(!done)
    {
      AXOM_ANNOTATE_SCOPE("GeometryClipper::clip3D");
      m_delegate->computeClipVolumes3D(ovlap);
    }
  }
}

/*!
  @brief Allocate a Delegate for the runtime policy of m_shapeeMesh.
*/
std::unique_ptr<GeometryClipper::Delegate> GeometryClipper::newDelegate()
{
  using RuntimePolicy = axom::runtime_policy::Policy;

  auto runtimePolicy = m_shapeeMesh.getRuntimePolicy();

  std::unique_ptr<Delegate> delegate;
  if(runtimePolicy == RuntimePolicy::seq)
  {
    delegate.reset(new detail::GeometryClipperDelegateExec<axom::SEQ_EXEC>(*this));
  }
#ifdef AXOM_RUNTIME_POLICY_USE_OPENMP
  else if(runtimePolicy == RuntimePolicy::omp)
  {
    delegate.reset(new detail::GeometryClipperDelegateExec<axom::OMP_EXEC>(*this));
  }
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_CUDA
  else if(runtimePolicy == RuntimePolicy::cuda)
  {
    delegate.reset(new detail::GeometryClipperDelegateExec<axom::CUDA_EXEC<256>>(*this));
  }
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_HIP
  else if(runtimePolicy == RuntimePolicy::hip)
  {
    delegate.reset(new detail::GeometryClipperDelegateExec<axom::HIP_EXEC<256>>(*this));
  }
#endif
  else
  {
    SLIC_ERROR(
      axom::fmt::format("GeometryClipper has no delegate for runtime policy {}", runtimePolicy));
  }
  return delegate;
}

}  // end namespace quest
}  // end namespace axom
