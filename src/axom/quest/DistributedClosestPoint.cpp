// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/runtime_policy.hpp"
#include "axom/core/NumericLimits.hpp"
#include "axom/quest/DistributedClosestPoint.hpp"
#include "axom/quest/detail/DistributedClosestPointImpl.hpp"

#include "axom/fmt.hpp"

#include "conduit_blueprint.hpp"
#include "conduit_blueprint_mpi.hpp"
#include "conduit_relay_mpi.hpp"

#include <cstdlib>

#ifndef AXOM_USE_MPI
  #error This file requires Axom to be configured with MPI
#endif
#include "mpi.h"

namespace axom
{
namespace quest
{
DistributedClosestPoint::DistributedClosestPoint()
  : m_mpiComm(MPI_COMM_WORLD)
  , m_mpiCommIsPrivate(false)
  , m_allocatorID(axom::INVALID_ALLOCATOR_ID)
  , m_sqDistanceThreshold(axom::numeric_limits<double>::max())
{
  setDefaultAllocatorID();
  setMpiCommunicator(MPI_COMM_WORLD);
}

DistributedClosestPoint::~DistributedClosestPoint()
{
  if(m_mpiCommIsPrivate)
  {
    int mpiIsFinalized = 0;
    MPI_Finalized(&mpiIsFinalized);
    if(!mpiIsFinalized)
    {
      MPI_Comm_free(&m_mpiComm);
    }
  }
}

void DistributedClosestPoint::setDefaultAllocatorID()
{
  int defaultAllocatorID = axom::INVALID_ALLOCATOR_ID;
  switch(m_runtimePolicy)
  {
  case RuntimePolicy::seq:
    defaultAllocatorID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
    break;

#ifdef AXOM_RUNTIME_POLICY_USE_OPENMP
  case RuntimePolicy::omp:
    defaultAllocatorID = axom::execution_space<axom::OMP_EXEC>::allocatorID();
    break;
#endif

#ifdef __CUDACC__
  #ifdef AXOM_RUNTIME_POLICY_USE_CUDA
  case RuntimePolicy::cuda:
    defaultAllocatorID = axom::execution_space<axom::CUDA_EXEC<256>>::allocatorID();
    break;
  #endif
#endif

#ifdef AXOM_RUNTIME_POLICY_USE_HIP
  case RuntimePolicy::hip:
    defaultAllocatorID = axom::execution_space<axom::HIP_EXEC<256>>::allocatorID();
    break;
#endif
  }
  if(defaultAllocatorID == axom::INVALID_ALLOCATOR_ID)
  {
    SLIC_ERROR(
      axom::fmt::format("There is no default allocator for runtime policy {}", m_runtimePolicy));
  }
  setAllocatorID(defaultAllocatorID);
}

void DistributedClosestPoint::setAllocatorID(int allocatorID)
{
  SLIC_ASSERT_MSG(allocatorID != axom::INVALID_ALLOCATOR_ID, "Invalid allocator id.");
  m_allocatorID = allocatorID;

  if(m_impl)
  {
    m_impl->setAllocatorID(m_allocatorID);
  }
}

void DistributedClosestPoint::setMpiCommunicator(MPI_Comm mpiComm, bool duplicate)
{
  if(m_mpiCommIsPrivate)
  {
    MPI_Comm_free(&m_mpiComm);
  }

  if(duplicate)
  {
    MPI_Comm_dup(mpiComm, &m_mpiComm);
  }
  else
  {
    m_mpiComm = mpiComm;
  }
  m_mpiCommIsPrivate = duplicate;
}

void DistributedClosestPoint::setDimension(int dim)
{
  SLIC_ERROR_IF(dim < 2 || dim > 3, "DistributedClosestPoint query only supports 2D or 3D queries");
  m_dimension = dim;
}

void DistributedClosestPoint::setDistanceThreshold(double threshold)
{
  SLIC_ERROR_IF(threshold < 0.0, "Distance threshold must be non-negative.");
  m_sqDistanceThreshold = threshold * threshold;
}

void DistributedClosestPoint::setOutput(const std::string& field, bool on)
{
  // clang-format off
  bool* f
    = field == "cp_rank" ? &m_outputRank
    : field == "cp_index" ? &m_outputIndex
    : field == "cp_distance" ? &m_outputDistance
    : field == "cp_coords" ? &m_outputCoords
    : field == "cp_domain_index" ? &m_outputDomainIndex
    : nullptr;
  // clang-format on
  SLIC_ERROR_IF(f == nullptr,
                axom::fmt::format("Invalid field '{}' should be one of these: "
                                  "cp_rank, cp_index, cp_distance, cp_coords, cp_domain_index",
                                  field));
  *f = on;
}

void DistributedClosestPoint::setObjectMesh(const conduit::Node& meshNode,
                                            const std::string& topologyName)
{
  if(m_impl)
  {
    SLIC_ERROR("Object mesh already created.");
  }

  const bool isMultidomain = conduit::blueprint::mesh::is_multi_domain(meshNode);

  // If meshNode isn't multidomain, create a temporary multidomain representation.
  std::shared_ptr<conduit::Node> tmpNode;
  if(!isMultidomain)
  {
    tmpNode = std::make_shared<conduit::Node>();
    conduit::blueprint::mesh::to_multi_domain(meshNode, *tmpNode);
  }
  const conduit::Node& mdMeshNode(isMultidomain ? meshNode : *tmpNode);
  verifyTopologyName(mdMeshNode, topologyName);

  auto domainCount = conduit::blueprint::mesh::number_of_domains(mdMeshNode);

  // Extract the dimension from the coordinate values group
  // use allreduce since some ranks might be empty
  {
    int localDim = -1;
    if(domainCount > 0)
    {
      const conduit::Node& domain0 = mdMeshNode.child(0);
      const conduit::Node& topology = domain0.fetch_existing("topologies/" + topologyName);
      const std::string coordsetName = topology.fetch_existing("coordset").as_string();
      const conduit::Node& coordset = domain0.fetch_existing("coordsets/" + coordsetName);
      const conduit::Node& coordsetValues = coordset.fetch_existing("values");
      localDim = internal::extractDimension(coordsetValues);
    }
    int dim = -1;
    MPI_Allreduce(&localDim, &dim, 1, MPI_INT, MPI_MAX, m_mpiComm);
    setDimension(dim);
  }

  allocateQueryInstance();

  m_impl->importObjectPoints(mdMeshNode, topologyName);

  return;
}

bool DistributedClosestPoint::generateBVHTree()
{
  SLIC_ASSERT_MSG(m_impl, "Must call 'setObjectMesh' before calling generateBVHTree");

  bool success = m_impl->generateBVHTree();
  return success;
}

void DistributedClosestPoint::computeClosestPoints(conduit::Node& query_node,
                                                   const std::string& topology)
{
  SLIC_ASSERT_MSG(m_impl, "Must call 'setObjectMesh' before calling generateBVHTree");

  SLIC_ASSERT(this->isValidBlueprint(query_node));

  m_impl->setSquaredDistanceThreshold(m_sqDistanceThreshold);
  m_impl->setMpiCommunicator(m_mpiComm);
  m_impl->setOutputSwitches(m_outputRank,
                            m_outputIndex,
                            m_outputDistance,
                            m_outputCoords,
                            m_outputDomainIndex);
  m_impl->computeClosestPoints(query_node, topology);
}

void DistributedClosestPoint::allocateQueryInstance()
{
  switch(m_runtimePolicy)
  {
  case RuntimePolicy::seq:
    m_dimension == 2 ? allocateQueryInstance<2, axom::SEQ_EXEC>()
                     : allocateQueryInstance<3, axom::SEQ_EXEC>();
    break;

#ifdef AXOM_RUNTIME_POLICY_USE_OPENMP
  case RuntimePolicy::omp:
    m_dimension == 2 ? allocateQueryInstance<2, axom::OMP_EXEC>()
                     : allocateQueryInstance<3, axom::OMP_EXEC>();
    break;
#endif

#ifdef AXOM_RUNTIME_POLICY_USE_CUDA
  case RuntimePolicy::cuda:
    m_dimension == 2 ? allocateQueryInstance<2, axom::CUDA_EXEC<256>>()
                     : allocateQueryInstance<3, axom::CUDA_EXEC<256>>();
    break;
#endif

#ifdef AXOM_RUNTIME_POLICY_USE_HIP
  case RuntimePolicy::hip:
    m_dimension == 2 ? allocateQueryInstance<2, axom::HIP_EXEC<256>>()
                     : allocateQueryInstance<3, axom::HIP_EXEC<256>>();
    break;
#endif

  default:
    SLIC_ERROR(
      axom::fmt::format("DistriburedClosestPoint: axom was not built for runtime policy {}."
                        "  Please select another policy.",
                        axom::runtime_policy::s_policyToName.at(m_runtimePolicy)));
  }
}

template <int DIM, typename ExecSpace>
void DistributedClosestPoint::allocateQueryInstance()
{
  m_impl = std::make_unique<internal::DistributedClosestPointExec<DIM, ExecSpace>>(m_allocatorID,
                                                                                   m_isVerbose);
}

bool DistributedClosestPoint::isValidBlueprint(const conduit::Node& mesh_node) const
{
  bool success = true;
  conduit::Node info;
  if(!conduit::blueprint::mpi::verify("mesh", mesh_node, info, m_mpiComm))
  {
    SLIC_INFO("Invalid blueprint for particle mesh: \n" << info.to_yaml());
    success = false;
  }

  return success;
}

void DistributedClosestPoint::verifyTopologyName(const conduit::Node& meshNode,
                                                 const std::string& topologyName)
{
  std::string coordsetPath;
  const std::string topologyPath = axom::fmt::format("topologies/{}", topologyName);
  for(axom::IndexType d = 0; d < meshNode.number_of_children(); ++d)
  {
    const auto& domain = meshNode.child(d);
    if(!domain.has_path(topologyPath))
    {
      auto errMsg = fmt::format("No such topology '{}' found.", topologyName);
      if(domain.has_path("coordsets/" + topologyName))
      {
        errMsg += fmt::format(
          "  You may have mistakenly specified a coordset name."
          "  The interface has changed to use topology name"
          " instead of coordset.");
      }
      SLIC_ERROR(errMsg);
    }
  }
}

}  // end namespace quest
}  // end namespace axom
