// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#ifndef AXOM_USE_RAJA
  #error "quest::GeometryClipper requires RAJA."
#endif

#include "axom/quest/GeometryClipper.hpp"
#include "axom/quest/detail/GeometryClipperDelegateExec.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/runtime_policy.hpp"
#include "axom/fmt.hpp"

namespace axom
{
namespace quest
{

GeometryClipper::GeometryClipper(
  quest::ShapeeMesh& shapeeMesh,
  const std::shared_ptr<quest::GeometryClipperStrategy>& strategy)
  : m_shapeeMesh(shapeeMesh)
  , m_strategy(strategy)
  , m_delegate(newDelegate())
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
  const int allocId = m_shapeeMesh.getAllocatorId();
  const axom::IndexType cellCount = m_shapeeMesh.getCellCount();

  axom::Array<char> labels;
  bool withInOut = m_strategy->labelInOut(m_shapeeMesh, labels);

  if(ovlap.size() < cellCount || ovlap.getAllocatorID() != allocId)
  {
    ovlap = axom::Array<double>(ArrayOptions::Uninitialized(),
                                cellCount,
                                cellCount,
                                allocId);
  }

  bool done = false;

  if(withInOut)
  {
    SLIC_ERROR_IF(labels.getAllocatorID() != allocId,
                  "GeometryClipperStrategy '" + m_strategy->name() + "' failed to provide 'labels' data with the required allocator id " + std::to_string(allocId));

    m_delegate->setCleanVolumeOverlaps(labels.view(), ovlap);

    axom::Array<axom::IndexType> unlabeledCells;
    m_delegate->collectUnlabeledCellIndices(labels.view(), unlabeledCells);

    done =
      m_strategy->specializedClip(m_shapeeMesh, ovlap.view(), unlabeledCells);

    if(!done)
    {
      m_delegate->computeClipVolumes3D(unlabeledCells.view(), ovlap.view());
    }
  }
  else
  {
    done = m_strategy->specializedClip(m_shapeeMesh, ovlap.view());

    if(!done)
    {
      m_delegate->computeClipVolumes3D(ovlap);
    }
  }
  if(done)
  {
    SLIC_ERROR_IF(labels.getAllocatorID() != allocId,
                  "GeometryClipperStrategy '" + m_strategy->name() + "' failed to provide 'ovlap' data with the required allocator id " + std::to_string(allocId));
  }

  if(!done)
  {
    if(m_shapeeMesh.dimension() == 3)
    {
      axom::Array<GeometryClipperStrategy::TetrahedronType> geomAsTets;
      axom::Array<GeometryClipperStrategy::OctahedronType> geomAsOcts;
      if(m_strategy->getShapeAsTets(m_shapeeMesh, geomAsTets))
      {
        SLIC_ASSERT("Code for clipping tets is incomplete.");
      }
      else if(m_strategy->getShapeAsOcts(m_shapeeMesh, geomAsOcts))
      {
        SLIC_ASSERT("Code for clipping octs is incomplete.");
      }
      else
      {
        SLIC_ERROR(
          "GeometryClipperStrategy subclass is missing required "
          "implementations.  See the documentation for "
          "GeometryClipperStrategy.");
      }
    }
    if(m_shapeeMesh.dimension() == 2)
    {
      SLIC_ERROR("GeometryClipper for 2D mesh is not implemented yet.");
    }
  }
}

/*!
  @brief Allocate a Delegate, for the runtime policy.
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
    delegate.reset(
      new detail::GeometryClipperDelegateExec<axom::CUDA_EXEC<256>>(*this));
  }
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_HIP
  else if(runtimePolicy == RuntimePolicy::hip)
  {
    delegate.reset(
      new detail::GeometryClipperDelegateExec<axom::HIP_EXEC<256>>(*this));
  }
#endif
  else
  {
    SLIC_ERROR(
      axom::fmt::format("GeometryClipper has no delegate for runtime policy {}",
                        runtimePolicy));
  }
  return delegate;
}

}  // end namespace quest
}  // end namespace axom
