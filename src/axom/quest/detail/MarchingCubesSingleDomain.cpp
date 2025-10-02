// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifndef AXOM_USE_CONDUIT
  #error "MarchingCubes.cpp requires conduit"
#endif
#include "conduit_blueprint.hpp"

#include "axom/core/execution/execution_space.hpp"
#include "axom/quest/detail/MarchingCubesSingleDomain.hpp"
#include "axom/quest/detail/MarchingCubesImpl.hpp"
#include "axom/fmt.hpp"

namespace axom
{
namespace quest
{
namespace detail
{
namespace marching_cubes
{
MarchingCubesSingleDomain::MarchingCubesSingleDomain(MarchingCubes& mc)
  : m_mc(mc)
  , m_runtimePolicy(mc.m_runtimePolicy)
  , m_allocatorID(mc.m_allocatorID)
  , m_dataParallelism(mc.m_dataParallelism)
  , m_dom(nullptr)
  , m_ndim(0)
  , m_topologyName()
  , m_fcnFieldName()
  , m_fcnPath()
  , m_maskFieldName()
  , m_maskPath()
{
  return;
}

void MarchingCubesSingleDomain::setDomain(const conduit::Node& dom,
                                          const std::string& topologyName,
                                          const std::string& maskField)
{
  m_topologyName = topologyName;

  SLIC_ASSERT_MSG(!conduit::blueprint::mesh::is_multi_domain(dom),
                  "Internal error.  Attempt to set a multi-domain mesh in "
                  "MarchingCubesSingleDomain.");
  SLIC_ASSERT(dom.fetch_existing("topologies/" + m_topologyName + "/type").as_string() ==
              "structured");

  const std::string coordsetPath =
    "coordsets/" + dom.fetch_existing("topologies/" + m_topologyName + "/coordset").as_string();
  SLIC_ASSERT(dom.has_path(coordsetPath));

  if(!m_maskPath.empty())
  {
    m_maskPath = maskField.empty() ? std::string() : "fields/" + maskField;
    SLIC_ASSERT(dom.has_path(m_maskPath + "/values"));
  }
  else
  {
    m_maskPath.clear();
  }

  m_dom = &dom;

  m_ndim = conduit::blueprint::mesh::topology::dims(
    dom.fetch_existing(axom::fmt::format("topologies/{}", m_topologyName)));
  SLIC_ASSERT(m_ndim >= 2 && m_ndim <= 3);

  SLIC_ASSERT_MSG(
    !conduit::blueprint::mcarray::is_interleaved(dom.fetch_existing(coordsetPath + "/values")),
    "MarchingCubes currently requires contiguous coordinates layout.");

  m_impl = newMarchingCubesImpl();

  m_impl->setDomain(dom, topologyName, maskField);
  m_impl->setDataParallelism(m_dataParallelism);
}

/*!
  @brief Allocate a MarchingCubesImpl object, template-specialized
  for caller-specified runtime policy and physical dimension.
*/
std::unique_ptr<MarchingCubesSingleDomain::ImplBase> MarchingCubesSingleDomain::newMarchingCubesImpl()
{
  SLIC_ASSERT(m_ndim >= 2 && m_ndim <= 3);
  std::unique_ptr<ImplBase> impl;
  if(m_runtimePolicy == MarchingCubes::RuntimePolicy::seq)
  {
    impl = m_ndim == 2
      ? std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<2, axom::SEQ_EXEC, axom::SEQ_EXEC>(m_mc.m_allocatorID,
                                                                   m_mc.m_caseIdsFlat,
                                                                   m_mc.m_crossingFlags,
                                                                   m_mc.m_scannedFlags,
                                                                   m_mc.m_facetIncrs))
      : std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<3, axom::SEQ_EXEC, axom::SEQ_EXEC>(m_mc.m_allocatorID,
                                                                   m_mc.m_caseIdsFlat,
                                                                   m_mc.m_crossingFlags,
                                                                   m_mc.m_scannedFlags,
                                                                   m_mc.m_facetIncrs));
  }
#ifdef AXOM_RUNTIME_POLICY_USE_OPENMP
  else if(m_runtimePolicy == MarchingCubes::RuntimePolicy::omp)
  {
    impl = m_ndim == 2
      ? std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<2, axom::OMP_EXEC, axom::SEQ_EXEC>(m_mc.m_allocatorID,
                                                                   m_mc.m_caseIdsFlat,
                                                                   m_mc.m_crossingFlags,
                                                                   m_mc.m_scannedFlags,
                                                                   m_mc.m_facetIncrs))
      : std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<3, axom::OMP_EXEC, axom::SEQ_EXEC>(m_mc.m_allocatorID,
                                                                   m_mc.m_caseIdsFlat,
                                                                   m_mc.m_crossingFlags,
                                                                   m_mc.m_scannedFlags,
                                                                   m_mc.m_facetIncrs));
  }
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_CUDA
  else if(m_runtimePolicy == MarchingCubes::RuntimePolicy::cuda)
  {
    impl = m_ndim == 2
      ? std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<2, axom::CUDA_EXEC<256>, axom::CUDA_EXEC<1>>(m_mc.m_allocatorID,
                                                                             m_mc.m_caseIdsFlat,
                                                                             m_mc.m_crossingFlags,
                                                                             m_mc.m_scannedFlags,
                                                                             m_mc.m_facetIncrs))
      : std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<3, axom::CUDA_EXEC<256>, axom::CUDA_EXEC<1>>(m_mc.m_allocatorID,
                                                                             m_mc.m_caseIdsFlat,
                                                                             m_mc.m_crossingFlags,
                                                                             m_mc.m_scannedFlags,
                                                                             m_mc.m_facetIncrs));
  }
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_HIP
  else if(m_runtimePolicy == MarchingCubes::RuntimePolicy::hip)
  {
    impl = m_ndim == 2
      ? std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<2, axom::HIP_EXEC<256>, axom::HIP_EXEC<1>>(m_mc.m_allocatorID,
                                                                           m_mc.m_caseIdsFlat,
                                                                           m_mc.m_crossingFlags,
                                                                           m_mc.m_scannedFlags,
                                                                           m_mc.m_facetIncrs))
      : std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<3, axom::HIP_EXEC<256>, axom::HIP_EXEC<1>>(m_mc.m_allocatorID,
                                                                           m_mc.m_caseIdsFlat,
                                                                           m_mc.m_crossingFlags,
                                                                           m_mc.m_scannedFlags,
                                                                           m_mc.m_facetIncrs));
  }
#endif
  else
  {
    SLIC_ERROR(
      axom::fmt::format("MarchingCubesSingleDomain has no implementation for runtime policy {}",
                        m_runtimePolicy));
  }
  return impl;
}

int32_t MarchingCubesSingleDomain::getDomainId(int32_t defaultId) const
{
  int rval = defaultId;
  if(m_dom->has_path("state/domain_id"))
  {
    rval = m_dom->fetch_existing("state/domain_id").to_int32();
  }
  return rval;
}

}  // namespace marching_cubes
}  // namespace detail
}  // end namespace quest
}  // end namespace axom
