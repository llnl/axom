// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_DISTRIBUTED_CLOSEST_POINT_H_
#define QUEST_DISTRIBUTED_CLOSEST_POINT_H_

#include "axom/config.hpp"
#include "axom/core/execution/runtime_policy.hpp"
#include "axom/slic.hpp"

#include "conduit_node.hpp"

#include <memory>
#include <cstdlib>

#ifndef AXOM_USE_MPI
  #error This file requires Axom to be configured with MPI
#endif
#include "mpi.h"

namespace axom
{
namespace quest
{
namespace internal
{
class DistributedClosestPointImpl;
}

/**
 * \brief Encapsulated the Distributed closest point query for a collection of query points
 * over an "object mesh"
 *
 * The object mesh and the query mesh are provided as conduit nodes
 * using the mesh blueprint schema.  Each of these are distributed
 * over the same MPI rank space.  Ranks are allowed to have any number
 * of domains, including zero.  This class orchestrates passing the
 * query points to all ranks whose object meshes might contain a
 * closest point.
 *
 * \note The class currently supports object meshes that are comprised
 * of a collection of points.  In the future, we'd like to consider
 * more general object meshes, e.g. triangle meshes.
 *
 * To use this class, first set some parameters, such as the runtime
 * execution policy, then pass in the object mesh and build a spatial
 * index over this mesh.  Finally, compute the closest points in the
 * object mesh to each point in a query mesh using the \a
 * computeClosestPoints() function.
 *
 * \note To prevent mixing unrelated MPI communications, you can set a
 * custom MPI Communicator using setMpiCommunicator().
 */
class DistributedClosestPoint
{
public:
  using RuntimePolicy = axom::runtime_policy::Policy;

public:
  DistributedClosestPoint();

  ~DistributedClosestPoint();

  /*!
    @brief Set runtime execution policy for local queries

    See axom::runtime_policy.
  */
  void setRuntimePolicy(RuntimePolicy policy)
  {
    SLIC_ASSERT_MSG(!m_impl, "Runtime policy may not change after setObjectMesh()");
    m_runtimePolicy = policy;
  }

  /*!  @brief Sets the allocator ID to the default associated with the
    execution policy
  */
  void setDefaultAllocatorID();

  /*!  @brief Sets the allocator ID.

    If not explitly set, the allocator ID is the default is the id
    associated with the runtime policy.
  */
  void setAllocatorID(int allocatorID);

  /**
   * \brief Set the MPI communicator.
   *
   * By default, the communicator is MPI_COMM_WORLD.
   *
   * \param [i] mpiComm The MPI communicator to use.
   * \param [i] duplicate Whether to duplicate mpiComm for exclusive use
   */
  void setMpiCommunicator(MPI_Comm mpiComm, bool duplicate = false);

  /**
   * \brief Sets the threshold for the query
   *
   * \param [in] threshold Ignore distances greater than this value.
   */
  void setDistanceThreshold(double threshold);

  /*!
    @brief Set what fields to output.

    @param [i] field Must be one of these:
      "cp_rank", "cp_index", "cp_distance", "cp_coords", "cp_domain_index".
    @param [i] on Whether to enable outputing @c field.

    By default, all are on.
  */
  void setOutput(const std::string& field, bool on);

  /// Sets the logging verbosity of the query. By default the query is not verbose
  void setVerbosity(bool isVerbose) { m_isVerbose = isVerbose; }

  /**
   * \brief Sets the object mesh for the query
   *
   * \param [in] meshNode Conduit node for the object mesh
   * \param [in] topologyName The name of the topology for the object mesh
   *
   * \pre \a meshNode must follow the mesh blueprint convention.
   * \pre Dimension of the mesh must be 2D or 3D
   *
   * \post The physical dimension and runtime policy are locked in,
   * and attempts to change them are disallowed.
   *
   * \internal This method also allocates the templated query instance,
   * which locks in the physical dimension and the runtime policy.
   *
   */
  void setObjectMesh(const conduit::Node& meshNode, const std::string& topologyName);

  /**
   * \brief Generates a BVH tree over the object mesh using the runtime execution policy
   *
   * \pre Users must set the object mesh before generating the BVH tree
   * \sa setObjectMesh()
   */
  bool generateBVHTree();

  /**
   * \brief Computes the closest point on the object mesh for each point
   * on the provided query mesh
   *
   * \param [in] query_node The root node of a mesh blueprint for the query points,
   * which can be empty if there are no query points for the calling rank
   * \param [in] topology The name of the topology within query_node
   *
   * @c queryMesh should have data on the host, regardless of the runtime
   * policy setting.  Data will be copied to device as needed.
   *
   * On completion, the query mesh contains the following fields:
   *   - cp_rank: will hold the rank of the object point containing the closest point
   *   - cp_domain_index: will hold the index of the object domain containing
   *     the closest points.
   *   - cp_index: Will hold the index of the closest object points.
   *     For multiple object mesh domains on a rank, cp_index is relative to
   *     each domain.
   *   - cp_coords: Will hold the coordinates of the closest points
   *     interleaved in a 1D array.
   *   - cp_distance: Will hold the distances to the closest points.
   * See setOutput() to toggle these outputs.
   *
   * \note The current implementation assumes that the mesh coordinates
   * are interleaved or contiguous.  The output cp_coords will be contiguous.
   */
  void computeClosestPoints(conduit::Node& query_node, const std::string& topology);

private:
  /*!
    @brief Allocate the DistributedClosestPointImpl object, which actually does the work.

    \post The implementation object is for a specific dimension and runtime policy,
    so changes to those are disallowed.
  */
  void allocateQueryInstance();

  //!@brief Allocate the DistributedClosestPointImpl for compile-time dimension and execution space.
  template <int DIM, typename ExecSpace>
  void allocateQueryInstance();

  /// Check validity of blueprint group
  bool isValidBlueprint(const conduit::Node& mesh_node) const;

  void verifyTopologyName(const conduit::Node& meshNode, const std::string& topologyName);

  void setDimension(int dim);

  RuntimePolicy m_runtimePolicy {RuntimePolicy::seq};
  MPI_Comm m_mpiComm;
  bool m_mpiCommIsPrivate;
  int m_allocatorID;
  int m_dimension {-1};
  bool m_isVerbose {false};
  double m_sqDistanceThreshold;

  bool m_outputRank = true;
  bool m_outputIndex = true;
  bool m_outputDistance = true;
  bool m_outputCoords = true;
  bool m_outputDomainIndex = true;

  //!@brief Instantiated implementation, created by setting object mesh.
  std::unique_ptr<internal::DistributedClosestPointImpl> m_impl;
};

}  // end namespace quest
}  // end namespace axom

#endif  //  QUEST_DISTRIBUTED_CLOSEST_POINT_H_
