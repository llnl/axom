// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifndef AXOM_USE_CONDUIT
  #error "MarchingCubesSingleDomain.cpp requires conduit"
#endif
#include "conduit_blueprint.hpp"

#include "axom/quest/detail/MarchingCubesSingleDomain.hpp"
#include "axom/fmt.hpp"

namespace axom::quest::detail::marching_cubes
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
  // The Bump implementation validates its supported topology types.
  if(!m_mc.m_useBumpBackend)
  {
    SLIC_ASSERT(dom.fetch_existing("topologies/" + m_topologyName + "/type").as_string() ==
                "structured");
  }

  const std::string coordsetPath =
    "coordsets/" + dom.fetch_existing("topologies/" + m_topologyName + "/coordset").as_string();
  SLIC_ASSERT(dom.has_path(coordsetPath));

  m_maskFieldName = maskField;
  if(!m_maskFieldName.empty())
  {
    m_maskPath = "fields/" + m_maskFieldName;
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

  // The Bump coordset dispatcher validates its supported layouts.
  if(!m_mc.m_useBumpBackend)
  {
    SLIC_ASSERT_MSG(
      !conduit::blueprint::mcarray::is_interleaved(dom.fetch_existing(coordsetPath + "/values")),
      "The legacy MarchingCubes backend requires a contiguous coordinate layout.");
  }

  m_impl = newMarchingCubesImpl();

  m_impl->setDomain(dom, topologyName, maskField);
  m_impl->setDataParallelism(m_dataParallelism);
}

std::unique_ptr<MarchingCubesSingleDomain::ImplBase> MarchingCubesSingleDomain::newMarchingCubesImpl()
{
  SLIC_ASSERT(m_ndim >= 2 && m_ndim <= 3);
  auto makeImplementation = [this](auto dimension) -> std::unique_ptr<ImplBase> {
    if(m_runtimePolicy == MarchingCubes::RuntimePolicy::seq)
    {
      return newMarchingCubesSeqImpl(dimension);
    }
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
    else if(m_runtimePolicy == MarchingCubes::RuntimePolicy::omp)
    {
      return newMarchingCubesOpenMPImpl(dimension);
    }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
    else if(m_runtimePolicy == MarchingCubes::RuntimePolicy::cuda)
    {
      return newMarchingCubesCudaImpl(dimension);
    }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
    else if(m_runtimePolicy == MarchingCubes::RuntimePolicy::hip)
    {
      return newMarchingCubesHipImpl(dimension);
    }
#endif
    SLIC_ERROR(
      axom::fmt::format("MarchingCubesSingleDomain has no implementation for runtime policy {}",
                        m_runtimePolicy));
    return nullptr;
  };

  return m_ndim == 2 ? makeImplementation(std::integral_constant<int, 2> {})
                     : makeImplementation(std::integral_constant<int, 3> {});
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

}  // end namespace axom::quest::detail::marching_cubes
