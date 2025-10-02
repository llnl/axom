// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_DISTRIBUTED_CLOSEST_POINT_IMPL_H_
#define QUEST_DISTRIBUTED_CLOSEST_POINT_IMPL_H_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/core/NumericLimits.hpp"
#include "axom/core/execution/runtime_policy.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/spin.hpp"

#include "axom/fmt.hpp"

#include "conduit_blueprint.hpp"
#include "conduit_blueprint_mcarray.hpp"
#include "conduit_blueprint_mpi.hpp"
#include "conduit_relay_mpi.hpp"
#include "conduit_relay_io.hpp"

#include <memory>
#include <cstdlib>
#include <cmath>
#include <list>
#include <vector>

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
// Utility function to dump a conduit node on each rank, e.g. for debugging
inline void dump_node(const conduit::Node& n,
                      const std::string&& fname,
                      const std::string& protocol = "json")
{
  conduit::relay::io::save(n, fname, protocol);
}

/**
 * \brief Utility function to get a typed pointer to the beginning of an array
 * stored by a conduit::Node
 */
template <typename T>
T* getPointer(conduit::Node& node)
{
  T* ptr = node.value();
  return ptr;
}

/**
 * \brief Utility function to create an axom::ArrayView over the array
 * of native types stored by a conduit::Node
 */
template <typename T>
axom::ArrayView<T> ArrayView_from_Node(conduit::Node& node, int sz)
{
  T* ptr = node.value();
  return axom::ArrayView<T>(ptr, sz);
}

/**
 * \brief Template specialization of ArrayView_from_Node for Point<double,2>
 *
 * \warning Assumes the underlying data is an MCArray with stride 2 access
 */
template <>
inline axom::ArrayView<primal::Point<double, 2>> ArrayView_from_Node(conduit::Node& node, int sz)
{
  using PointType = primal::Point<double, 2>;

  PointType* ptr = static_cast<PointType*>(node.data_ptr());
  return axom::ArrayView<PointType>(ptr, sz);
}

/**
 * \brief Template specialization of ArrayView_from_Node for Point<double,3>
 *
 * \warning Assumes the underlying data is an MCArray with stride 3 access
 */
template <>
inline axom::ArrayView<primal::Point<double, 3>> ArrayView_from_Node(conduit::Node& node, int sz)
{
  using PointType = primal::Point<double, 3>;

  PointType* ptr = static_cast<PointType*>(node.data_ptr());
  return axom::ArrayView<PointType>(ptr, sz);
}

/**
 * \brief Put BoundingBox into a Conduit Node.
 */
template <int NDIMS>
void put_bounding_box_to_conduit_node(const primal::BoundingBox<double, NDIMS>& bb,
                                      conduit::Node& node)
{
  node["dim"].set(bb.dimension());
  if(bb.isValid())
  {
    node["lo"].set(bb.getMin().data(), bb.dimension());
    node["hi"].set(bb.getMax().data(), bb.dimension());
  }
}

/**
 * \brief Get BoundingBox from a Conduit Node.
 */
template <int NDIMS>
void get_bounding_box_from_conduit_node(primal::BoundingBox<double, NDIMS>& bb,
                                        const conduit::Node& node)
{
  using PointType = primal::Point<double, NDIMS>;

  SLIC_ASSERT(NDIMS == node.fetch_existing("dim").as_int());

  bb.clear();

  if(node.has_child("lo"))
  {
    bb.addPoint(PointType(node.fetch_existing("lo").as_double_ptr(), NDIMS));
    bb.addPoint(PointType(node.fetch_existing("hi").as_double_ptr(), NDIMS));
  }
}

/// Helper function to extract the dimension from the coordinate values group
/// of a mesh blueprint coordset
inline int extractDimension(const conduit::Node& values_node)
{
  SLIC_ASSERT(values_node.has_child("x"));
  return values_node.has_child("z") ? 3 : (values_node.has_child("y") ? 2 : 1);
}

/// Helper function to extract the number of points from the coordinate values group
/// of a mesh blueprint coordset
inline int extractSize(const conduit::Node& values_node)
{
  SLIC_ASSERT(values_node.has_child("x"));
  return values_node["x"].dtype().number_of_elements();
}

namespace relay
{
namespace mpi
{
/**
 * \brief Sends a conduit node along with its schema using MPI_Isend
 *
 * \param [in] node node to send
 * \param [in] dest ID of MPI rank to send to
 * \param [in] tag tag for MPI message
 * \param [in] comm MPI communicator to use
 * \param [in] request object holding state for the sent data
 * \note Adapted from conduit's relay::mpi's \a send_using_schema and \a isend
 * to use non-blocking \a MPI_Isend instead of blocking \a MPI_Send
 */
inline int isend_using_schema(conduit::Node& node,
                              int dest,
                              int tag,
                              MPI_Comm comm,
                              conduit::relay::mpi::Request* request)
{
  conduit::Schema s_data_compact;

  // schema will only be valid if compact and contig
  if(node.is_compact() && node.is_contiguous())
  {
    s_data_compact = node.schema();
  }
  else
  {
    node.schema().compact_to(s_data_compact);
  }
  const std::string snd_schema_json = s_data_compact.to_json();

  conduit::Schema s_msg;
  s_msg["schema_len"].set(conduit::DataType::int64());
  s_msg["schema"].set(conduit::DataType::char8_str(snd_schema_json.size() + 1));
  s_msg["data"].set(s_data_compact);

  // create a compact schema to use
  conduit::Schema s_msg_compact;
  s_msg.compact_to(s_msg_compact);
  request->m_buffer.reset();
  request->m_buffer.set_schema(s_msg_compact);

  // set up the message's node using this schema
  request->m_buffer["schema_len"].set((std::int64_t)snd_schema_json.length());
  request->m_buffer["schema"].set(snd_schema_json);
  request->m_buffer["data"].update(node);

  // for wait_all,  this must always be NULL except for
  // the irecv cases where copy out is necessary
  // isend case must always be NULL
  request->m_rcv_ptr = nullptr;

  auto msg_data_size = request->m_buffer.total_bytes_compact();
  int mpi_error = MPI_Isend(const_cast<void*>(request->m_buffer.data_ptr()),
                            static_cast<int>(msg_data_size),
                            MPI_BYTE,
                            dest,
                            tag,
                            comm,
                            &(request->m_request));

  // Error checking -- Note: expansion of CONDUIT_CHECK_MPI_ERROR
  if(static_cast<int>(mpi_error) != MPI_SUCCESS)
  {
    char check_mpi_err_str_buff[MPI_MAX_ERROR_STRING];
    int check_mpi_err_str_len = 0;
    MPI_Error_string(mpi_error, check_mpi_err_str_buff, &check_mpi_err_str_len);

    SLIC_ERROR(fmt::format("MPI call failed: error code = {} error message = {}",
                           mpi_error,
                           check_mpi_err_str_buff));
  }

  return mpi_error;
}

}  // namespace mpi
}  // namespace relay

/*!
  @brief Non-templated base class for the distributed closest point
  implementation.

  This class provides an abstract base class handle for
  DistributedClosestPointExec, which generically implements the code
  for templated dimensions and execution spaces.
  This class implements the non-templated parts of the implementation.
  The two are highly coupled.
*/
class DistributedClosestPointImpl
{
public:
  DistributedClosestPointImpl(int allocatorID, bool isVerbose)
    : m_allocatorID(allocatorID)
    , m_isVerbose(isVerbose)
    , m_mpiComm(MPI_COMM_NULL)
    , m_rank(-1)
    , m_nranks(-1)
    , m_sqDistanceThreshold(axom::numeric_limits<double>::max())
  { }

  virtual ~DistributedClosestPointImpl() { }

  virtual int getDimension() const = 0;

  /*!  @brief Sets the allocator ID to the default associated with the
    execution policy
  */
  void setAllocatorID(int allocatorID)
  {
    SLIC_ASSERT(allocatorID != axom::INVALID_ALLOCATOR_ID);
    // TODO: If appropriate, how to check for compatibility with runtime policy?
    m_allocatorID = allocatorID;
  }

  /*!
   @brief Import object mesh points from the object blueprint mesh into internal memory.

   @param [in] mdMeshNode The blueprint mesh containing the object points.
   @param [in] topologyName Name of the blueprint topology in \a mdMeshNode.
   @note This function currently supports mesh blueprints with the "point" topology
  */
  virtual void importObjectPoints(const conduit::Node& mdMeshNode,
                                  const std::string& topologyName) = 0;

  //! @brief Generates the BVH tree for the classes execution space
  virtual bool generateBVHTree() = 0;

  /*!
   @brief Set the MPI communicator.
  */
  void setMpiCommunicator(MPI_Comm mpiComm)
  {
    m_mpiComm = mpiComm;
    MPI_Comm_rank(m_mpiComm, &m_rank);
    MPI_Comm_size(m_mpiComm, &m_nranks);
  }

  /*!
   @brief Sets the threshold for the query

   @param [in] threshold Ignore distances greater than this value.
  */
  void setSquaredDistanceThreshold(double sqThreshold)
  {
    SLIC_ERROR_IF(sqThreshold < 0.0, "Squared distance-threshold must be non-negative.");
    m_sqDistanceThreshold = sqThreshold;
  }

  /*!
    @brief Set which output data fields to generate.
  */
  void setOutputSwitches(bool outputRank,
                         bool outputIndex,
                         bool outputDistance,
                         bool outputCoords,
                         bool outputDomainIndex)
  {
    m_outputRank = outputRank;
    m_outputIndex = outputIndex;
    m_outputDistance = outputDistance;
    m_outputCoords = outputCoords;
    m_outputDomainIndex = outputDomainIndex;
  }

  /*!
   * Copy parts of query mesh partition to a conduit::Node for
   * computation and communication.
   * queryNode must be a blueprint multidomain mesh.
   */
  void node_copy_query_to_xfer(conduit::Node& queryNode,
                               conduit::Node& xferNode,
                               const std::string& topologyName) const
  {
    xferNode["homeRank"] = m_rank;
    xferNode["is_first"] = 1;

    const bool isMultidomain = conduit::blueprint::mesh::is_multi_domain(queryNode);
    const auto domainCount = conduit::blueprint::mesh::number_of_domains(queryNode);
    conduit::Node& xferDoms = xferNode["xferDoms"];
    for(conduit::index_t domainNum = 0; domainNum < domainCount; ++domainNum)
    {
      auto& queryDom = isMultidomain ? queryNode.child(domainNum) : queryNode;

      const std::string coordsetName =
        queryDom.fetch_existing(axom::fmt::format("topologies/{}/coordset", topologyName)).as_string();
      const std::string& domName = queryDom.name();
      conduit::Node& xferDom = xferDoms[domName];
      conduit::Node& queryCoords = queryDom.fetch_existing(fmt::format("coordsets/{}", coordsetName));
      conduit::Node& queryCoordsValues = queryCoords.fetch_existing("values");

      const int dim = internal::extractDimension(queryCoordsValues);
      const int qPtCount = internal::extractSize(queryCoordsValues);
      xferDom["qPtCount"] = qPtCount;
      xferDom["dim"] = dim;

      copy_components_to_interleaved(queryCoordsValues, xferDom["coords"]);

      constexpr bool isInt32 = std::is_same<axom::IndexType, std::int32_t>::value;
      auto dtype = isInt32 ? conduit::DataType::int32() : conduit::DataType::int64();
      dtype.set_number_of_elements(qPtCount);
      xferDom["cp_index"].set_dtype(dtype);
      xferDom["cp_rank"].set_dtype(dtype);
      xferDom["cp_domain_index"].set_dtype(dtype);
      xferDom["debug/cp_distance"].set_dtype(conduit::DataType::float64(qPtCount));
      xferDom["cp_coords"].set_dtype(conduit::DataType::float64(dim * qPtCount));
    }
  }

  /// Copy xferNode back to query mesh partition.
  void node_copy_xfer_to_query(conduit::Node& xferNode,
                               conduit::Node& queryNode,
                               const std::string& topologyName) const
  {
    const bool isMultidomain = conduit::blueprint::mesh::is_multi_domain(queryNode);
    const auto domainCount = conduit::blueprint::mesh::number_of_domains(queryNode);
    conduit::Node& xferDoms = xferNode.fetch_existing("xferDoms");
    SLIC_ASSERT(xferDoms.number_of_children() == domainCount);
    for(conduit::index_t domainNum = 0; domainNum < domainCount; ++domainNum)
    {
      auto& queryDom = isMultidomain ? queryNode.child(domainNum) : queryNode;
      conduit::Node& xferDom = xferDoms.child(domainNum);
      conduit::Node& fields = queryDom.fetch_existing("fields");

      conduit::Node genericHeaders;
      genericHeaders["association"] = "vertex";
      genericHeaders["topology"] = topologyName;

      if(m_outputRank)
      {
        auto& src = xferDom.fetch_existing("cp_rank");
        auto& dst = fields["cp_rank"];
        dst.set_node(genericHeaders);
        dst["values"].move(src);
      }

      if(m_outputIndex)
      {
        auto& src = xferDom.fetch_existing("cp_index");
        auto& dst = fields["cp_index"];
        dst.set_node(genericHeaders);
        dst["values"].move(src);
      }

      if(m_outputDomainIndex)
      {
        auto& src = xferDom.fetch_existing("cp_domain_index");
        auto& dst = fields["cp_domain_index"];
        dst.set_node(genericHeaders);
        dst["values"].move(src);
      }

      if(m_outputDistance)
      {
        auto& src = xferDom.fetch_existing("debug/cp_distance");
        auto& dst = fields["cp_distance"];
        dst.set_node(genericHeaders);
        dst["values"].move(src);
      }

      if(m_outputCoords)
      {
        auto& dst = fields["cp_coords"];
        dst.set_node(genericHeaders);
        auto& dstValues = dst["values"];
        copy_interleaved_to_components(xferDom.fetch_existing("cp_coords"), dstValues);
      }
    }
  }

  /*
    Special copy from coordinates (in a format that's not
    necessarily interleaved) to a 1D array of interleaved values).
    If coordinates are already interleaved, copy pointer.
  */
  void copy_components_to_interleaved(conduit::Node& components, conduit::Node& interleaved) const
  {
    const int dim = getDimension();
    const int qPtCount = internal::extractSize(components);
    bool interleavedSrc = conduit::blueprint::mcarray::is_interleaved(components);
    if(interleavedSrc)
    {
      interleaved.set_external(internal::getPointer<double>(components.child(0)), dim * qPtCount);
    }
    else
    {
      // Copy from component-wise src to 1D-interleaved dst.
      interleaved.reset();
      interleaved.set_dtype(conduit::DataType::float64(dim * qPtCount));
      for(int d = 0; d < dim; ++d)
      {
        auto src = components.child(d).as_float64_array();
        double* dst = interleaved.as_float64_ptr() + d;
        for(int i = 0; i < qPtCount; ++i)
        {
          dst[i * dim] = src[i];
        }
      }
    }
  }

  /*
    Special copy from 1D interleaved coordinate values back to
    component-wise storage.
    This is a nop if they point to the same data.
  */
  void copy_interleaved_to_components(const conduit::Node& interleaved, conduit::Node& components) const
  {
    const int dim = getDimension();
    const int qPtCount = interleaved.dtype().number_of_elements() / dim;
    components.reset();
    // Copy from 1D-interleaved src to component-wise dst.
    for(int d = 0; d < dim; ++d)
    {
      const double* src = interleaved.as_float64_ptr() + d;
      auto& dstNode = components.append();
      dstNode.set_dtype(conduit::DataType(interleaved.dtype().id(), qPtCount));
      double* dst = dstNode.as_float64_ptr();
      for(int i = 0; i < qPtCount; ++i)
      {
        dst[i] = src[i * dim];
      }
    }
  }

  /// Wait for some non-blocking sends (if any) to finish.
  void check_send_requests(std::list<conduit::relay::mpi::Request>& isendRequests, bool atLeastOne) const
  {
    std::vector<MPI_Request> reqs;
    for(auto& isr : isendRequests)
    {
      reqs.push_back(isr.m_request);
    }

    int inCount = static_cast<int>(reqs.size());
    int outCount = 0;
    std::vector<int> indices(reqs.size(), -1);
    if(atLeastOne)
    {
      MPI_Waitsome(inCount, reqs.data(), &outCount, indices.data(), MPI_STATUSES_IGNORE);
    }
    else
    {
      MPI_Testsome(inCount, reqs.data(), &outCount, indices.data(), MPI_STATUSES_IGNORE);
    }
    indices.resize(outCount);

    auto reqIter = isendRequests.begin();
    int prevIdx = 0;
    for(const int idx : indices)
    {
      for(; prevIdx < idx; ++prevIdx)
      {
        ++reqIter;
      }
      reqIter = isendRequests.erase(reqIter);
      ++prevIdx;
    }
  }

  virtual void computeClosestPoints(conduit::Node& queryMesh,
                                    const std::string& topologyName) const = 0;

protected:
  int m_allocatorID;
  bool m_isVerbose;

  MPI_Comm m_mpiComm;
  int m_rank;
  int m_nranks;

  double m_sqDistanceThreshold;

  bool m_outputRank = true;
  bool m_outputIndex = true;
  bool m_outputDistance = true;
  bool m_outputCoords = true;
  bool m_outputDomainIndex = true;

  struct MinCandidate
  {
    /// Squared distance to query point
    double sqDist {numerics::floating_point_limits<double>::max()};
    /// Index of domain of closest element
    int domainIdx {-1};
    /// Index within domain of closest element
    int pointIdx {-1};
    /// MPI rank of closest element
    int rank {-1};
  };
};

/*!
  \brief Implements the DistributedClosestPoint query for
  compile-time dimension and execution space.

  This class implements closest point search parts that depend
  on dimension and execution space.

  \tparam NDIMS The dimension of the object mesh and query points
  \tparam ExecSpace The general execution space, such as axom::SEQ_EXEC and
  axom::CUDA_EXEC<256>.
*/
template <int NDIMS, typename ExecSpace>
class DistributedClosestPointExec : public DistributedClosestPointImpl
{
public:
  static constexpr int DIM = NDIMS;
  using LoopPolicy = typename execution_space<ExecSpace>::loop_policy;
  using ReducePolicy = typename execution_space<ExecSpace>::reduce_policy;
  using PointType = primal::Point<double, DIM>;
  using BoxType = primal::BoundingBox<double, DIM>;
  using PointArray = axom::Array<PointType>;
  using BoxArray = axom::Array<BoxType>;
  using BVHTreeType = spin::BVH<DIM, ExecSpace>;

  /*!
    @brief Constructor

    @param [i] allocatorID Allocator ID, which must be compatible with
      @c ExecSpace.  See axom::allocate and axom::reallocate.
      Also see setAllocatorID().
    @param [i[ isVerbose
  */
  DistributedClosestPointExec(int allocatorID, bool isVerbose)
    : DistributedClosestPointImpl(allocatorID, isVerbose)
    , m_objectPtCoords(0, 0, allocatorID)
    , m_objectPtDomainIds(0, 0, allocatorID)
  {
    SLIC_ASSERT(allocatorID != axom::INVALID_ALLOCATOR_ID);

    setMpiCommunicator(MPI_COMM_WORLD);
  }

  int getDimension() const override { return DIM; }

  void importObjectPoints(const conduit::Node& mdMeshNode, const std::string& topologyName) override
  {
    // TODO: See if some of the copies in this method can be optimized out.

    SLIC_ASSERT(sizeof(double) * DIM == sizeof(PointType));

    // Count points in the mesh.
    int ptCount = 0;
    for(const conduit::Node& domain : mdMeshNode.children())
    {
      const std::string coordsetName =
        domain.fetch_existing(axom::fmt::format("topologies/{}/coordset", topologyName)).as_string();
      const std::string valuesPath = axom::fmt::format("coordsets/{}/values", coordsetName);
      auto& values = domain.fetch_existing(valuesPath);
      const int N = internal::extractSize(values);
      ptCount += N;
    }

    // Copy points to internal memory
    PointArray coords(ptCount, ptCount);
    axom::Array<axom::IndexType> domIds(ptCount, ptCount);
    std::size_t copiedCount = 0;
    conduit::Node tmpValues;
    for(axom::IndexType d = 0; d < mdMeshNode.number_of_children(); ++d)
    {
      const conduit::Node& domain = mdMeshNode.child(d);

      axom::IndexType domainId = d;
      if(domain.has_path("state/domain_id"))
      {
        domainId = domain.fetch_existing("state/domain_id").to_int32();
      }

      const std::string coordsetName =
        domain.fetch_existing(axom::fmt::format("topologies/{}/coordset", topologyName)).as_string();
      const std::string valuesPath = axom::fmt::format("coordsets/{}/values", coordsetName);

      auto& values = domain.fetch_existing(valuesPath);

      bool isInterleaved = conduit::blueprint::mcarray::is_interleaved(values);
      if(!isInterleaved)
      {
        conduit::blueprint::mcarray::to_interleaved(values, tmpValues);
      }
      const conduit::Node& copySrc = isInterleaved ? values : tmpValues;

      const int N = internal::extractSize(copySrc);
      const std::size_t nBytes = sizeof(double) * DIM * N;

      axom::copy(coords.data() + copiedCount, copySrc.fetch_existing("x").data_ptr(), nBytes);
      tmpValues.reset();

      domIds.fill(domainId, N, copiedCount);

      copiedCount += N;
    }
    // copy computed data to ExecSpace
    m_objectPtCoords = PointArray(coords, m_allocatorID);
    m_objectPtDomainIds = axom::Array<axom::IndexType>(domIds, m_allocatorID);
  }

  bool generateBVHTree() override
  {
    // Delegates to generateBVHTreeImpl<> which uses
    // the execution space templated bvh tree

    SLIC_ASSERT_MSG(!m_bvh, "BVH tree already initialized");

    // In case user changed the allocator after setObjectMesh,
    // move the object point data to avoid repetitive page faults.
    if(m_objectPtCoords.getAllocatorID() != m_allocatorID)
    {
      PointArray tmpPoints(m_objectPtCoords, m_allocatorID);
      m_objectPtCoords.swap(tmpPoints);
    }

    m_bvh = std::make_unique<BVHTreeType>();
    return generateBVHTreeImpl(m_bvh.get());
  }

  /// Get local copy of all ranks BVH root bounding boxes.
  void gatherBVHRoots()
  {
    SLIC_ASSERT_MSG(m_bvh, "BVH tree must be initialized before calling 'gatherBVHRoots");

    BoxType local_bb = m_bvh->getBounds();
    gatherBoundingBoxes(local_bb, m_objectPartitionBbs);
  }

  /// Allgather one bounding box from each rank.
  void gatherBoundingBoxes(const BoxType& aabb, BoxArray& all_aabbs) const
  {
    axom::Array<double> sendbuf(2 * DIM);
    aabb.getMin().to_array(&sendbuf[0]);
    aabb.getMax().to_array(&sendbuf[DIM]);
    axom::Array<double> recvbuf(m_nranks * sendbuf.size());
    // Note: Using axom::Array<double,2> may reduce clutter a tad.
    int errf = MPI_Allgather(sendbuf.data(),
                             2 * DIM,
                             mpi_traits<double>::type,
                             recvbuf.data(),
                             2 * DIM,
                             mpi_traits<double>::type,
                             m_mpiComm);
    SLIC_ASSERT(errf == MPI_SUCCESS);
    AXOM_UNUSED_VAR(errf);

    all_aabbs.clear();
    all_aabbs.reserve(m_nranks);
    for(int i = 0; i < m_nranks; ++i)
    {
      PointType lower(&recvbuf[i * 2 * DIM]);
      PointType upper(&recvbuf[i * 2 * DIM + DIM]);
      all_aabbs.emplace_back(BoxType(lower, upper, false));
    }
  }

  /// Compute bounding box for local part of a mesh.
  BoxType computeMeshBoundingBox(conduit::Node& xferNode) const
  {
    BoxType rval;

    conduit::Node& xferDoms = xferNode.fetch_existing("xferDoms");
    for(conduit::Node& xferDom : xferDoms.children())
    {
      const int qPtCount = xferDom.fetch_existing("qPtCount").value();

      /// Extract fields from the input node as ArrayViews
      auto queryPts = ArrayView_from_Node<PointType>(xferDom.fetch_existing("coords"), qPtCount);
      for(const auto& p : queryPts)
      {
        rval.addPoint(p);
      }
    }

    return rval;
  }

  /**
   * \brief Implementation of the user-facing
   * DistributedClosestPoint::computeClosestPoints() method.
   *
   * We use non-blocking sends for performance and deadlock avoidance.
   * The worst case could incur nranks^2 sends.  To avoid excessive
   * buffer usage, we occasionally check the sends for completion,
   * using check_send_requests().
   */
  void computeClosestPoints(conduit::Node& queryMesh, const std::string& topologyName) const override
  {
    SLIC_ASSERT_MSG(m_bvh, "BVH tree must be initialized before calling 'computeClosestPoints");

    std::map<int, std::shared_ptr<conduit::Node>> xferNodes;

    // create conduit Node containing data that has to xfer between ranks.
    // The node will be mostly empty if there are no domains on this rank
    {
      xferNodes[m_rank] = std::make_shared<conduit::Node>();
      conduit::Node& xferNode = *xferNodes[m_rank];
      node_copy_query_to_xfer(queryMesh, xferNode, topologyName);
      xferNode["homeRank"] = m_rank;
    }

    BoxType myQueryBb = computeMeshBoundingBox(*xferNodes[m_rank]);
    put_bounding_box_to_conduit_node(myQueryBb, xferNodes[m_rank]->fetch("aabb"));
    BoxArray allQueryBbs;
    gatherBoundingBoxes(myQueryBb, allQueryBbs);

    {
      conduit::Node& xferNode = *xferNodes[m_rank];
      computeLocalClosestPoints(xferNode);
    }

    const auto& myObjectBb = m_objectPartitionBbs[m_rank];
    int remainingRecvs = 0;
    for(int r = 0; r < m_nranks; ++r)
    {
      if(r != m_rank)
      {
        const auto& otherQueryBb = allQueryBbs[r];
        double sqDistance = axom::primal::squared_distance(otherQueryBb, myObjectBb);
        if(sqDistance <= m_sqDistanceThreshold)
        {
          ++remainingRecvs;
        }
      }
    }

    // arbitrary tags for send/recv xferNode.
    const int tag = 987342;

    std::list<conduit::relay::mpi::Request> isendRequests;

    {
      /*
        Send local query mesh to next rank with close-enough object
        partition, if any.  Increase remainingRecvs, because this data
        will come back.
      */
      int firstRecipForMyQuery = next_recipient(*xferNodes[m_rank]);
      if(m_nranks == 1)
      {
        SLIC_ASSERT(firstRecipForMyQuery == -1);
      }

      if(firstRecipForMyQuery == -1)
      {
        // No need to send anywhere.  Put computed data back into queryMesh.
        node_copy_xfer_to_query(*xferNodes[m_rank], queryMesh, topologyName);
        xferNodes.erase(m_rank);
      }
      else
      {
        isendRequests.emplace_back(conduit::relay::mpi::Request());
        auto& req = isendRequests.back();
        relay::mpi::isend_using_schema(*xferNodes[m_rank], firstRecipForMyQuery, tag, m_mpiComm, &req);
        ++remainingRecvs;
      }
    }

    while(remainingRecvs > 0)
    {
      SLIC_INFO_IF(m_isVerbose,
                   fmt::format("=======  {} receives remaining =======", remainingRecvs));

      // Receive the next xferNode
      std::shared_ptr<conduit::Node> recvXferNodePtr = std::make_shared<conduit::Node>();
      conduit::relay::mpi::recv_using_schema(*recvXferNodePtr, MPI_ANY_SOURCE, tag, m_mpiComm);

      const int homeRank = recvXferNodePtr->fetch_existing("homeRank").as_int();
      --remainingRecvs;
      xferNodes[homeRank] = recvXferNodePtr;
      conduit::Node& xferNode = *xferNodes[homeRank];

      if(homeRank == m_rank)
      {
        node_copy_xfer_to_query(xferNode, queryMesh, topologyName);
      }
      else
      {
        computeLocalClosestPoints(xferNode);

        isendRequests.emplace_back(conduit::relay::mpi::Request());
        auto& isendRequest = isendRequests.back();
        int nextRecipient = next_recipient(xferNode);
        SLIC_ASSERT(nextRecipient != -1);
        relay::mpi::isend_using_schema(xferNode, nextRecipient, tag, m_mpiComm, &isendRequest);

        // Check non-blocking sends to free memory.
        check_send_requests(isendRequests, false);
      }

    }  // remainingRecvs loop

    // Complete remaining non-blocking sends.
    while(!isendRequests.empty())
    {
      check_send_requests(isendRequests, true);
    }

    MPI_Barrier(m_mpiComm);
    slic::flushStreams();
  }

private:
  /**
    Determine the next rank (in ring order) with an object partition
    close to the query points in xferNode.  The intent is to send
    xferNode there next.
  */
  int next_recipient(const conduit::Node& xferNode) const
  {
    int homeRank = xferNode.fetch_existing("homeRank").value();
    BoxType bb;
    get_bounding_box_from_conduit_node(bb, xferNode.fetch_existing("aabb"));
    for(int i = 1; i < m_nranks; ++i)
    {
      int maybeNextRecip = (m_rank + i) % m_nranks;
      if(maybeNextRecip == homeRank)
      {
        return maybeNextRecip;
      }
      double sqDistance = primal::squared_distance(bb, m_objectPartitionBbs[maybeNextRecip]);
      if(sqDistance <= m_sqDistanceThreshold)
      {
        return maybeNextRecip;
      }
    }
    return -1;
  }

  // Note: following should be private, but nvcc complains about lambdas in private scope
public:
  /// Templated implementation of generateBVHTree function
  bool generateBVHTreeImpl(BVHTreeType* bvh)
  {
    SLIC_ASSERT(bvh != nullptr);

    const int npts = m_objectPtCoords.size();
    axom::Array<BoxType> boxesArray(npts, npts, m_allocatorID);
    auto boxesView = boxesArray.view();
    auto pointsView = m_objectPtCoords.view();

    axom::for_all<ExecSpace>(
      npts,
      AXOM_LAMBDA(axom::IndexType i) { boxesView[i] = BoxType {pointsView[i]}; });

    // Build bounding volume hierarchy
    bvh->setAllocatorID(m_allocatorID);
    int result = bvh->initialize(boxesView, npts);

    gatherBVHRoots();

    return (result == spin::BVH_BUILD_OK);
  }

  void computeLocalClosestPoints(conduit::Node& xferNode) const
  {
    using axom::primal::squared_distance;

    // Note: There is some additional computation the first time this function
    // is called for a query node, even if the local object mesh is empty
    const bool hasObjectPoints = m_objectPtCoords.size() > 0;
    const bool is_first = xferNode.has_path("is_first");
    if(!hasObjectPoints && !is_first)
    {
      return;
    }
    conduit::Node& xferDoms = xferNode["xferDoms"];
    for(conduit::Node& xferDom : xferDoms.children())
    {
      // --- Set up arrays and views in the execution space
      // Arrays are initialized in that execution space the first time
      // they are processed and are copied in during subsequent
      // processing

      // Check dimension and extract the number of points
      SLIC_ASSERT(xferDom.fetch_existing("dim").as_int() == DIM);
      const int qPtCount = xferDom.fetch_existing("qPtCount").value();

      /// Extract fields from the input node as ArrayViews
      auto queryPts = ArrayView_from_Node<PointType>(xferDom.fetch_existing("coords"), qPtCount);
      auto cpIndexes =
        ArrayView_from_Node<axom::IndexType>(xferDom.fetch_existing("cp_index"), qPtCount);
      auto cpDomainIndexes =
        ArrayView_from_Node<axom::IndexType>(xferDom.fetch_existing("cp_domain_index"), qPtCount);
      auto cpRanks =
        ArrayView_from_Node<axom::IndexType>(xferDom.fetch_existing("cp_rank"), qPtCount);
      auto cpCoords = ArrayView_from_Node<PointType>(xferDom.fetch_existing("cp_coords"), qPtCount);

      /// Create ArrayViews in ExecSpace that are compatible with fields
      // This deep-copies host memory in xferDom to device memory.
      // TODO: Avoid copying arrays (here and at the end) if both are on the host
      auto cp_idx = is_first ? axom::Array<axom::IndexType>(qPtCount, qPtCount, m_allocatorID)
                             : axom::Array<axom::IndexType>(cpIndexes, m_allocatorID);
      auto cp_domidx = is_first ? axom::Array<axom::IndexType>(qPtCount, qPtCount, m_allocatorID)
                                : axom::Array<axom::IndexType>(cpDomainIndexes, m_allocatorID);
      auto cp_rank = is_first ? axom::Array<axom::IndexType>(qPtCount, qPtCount, m_allocatorID)
                              : axom::Array<axom::IndexType>(cpRanks, m_allocatorID);

      /// PROBLEM: The striding does not appear to be retained by conduit relay
      ///          We might need to transform it? or to use a single array w/ pointers into it?
      auto cp_pos = is_first ? axom::Array<PointType>(qPtCount, qPtCount, m_allocatorID)
                             : axom::Array<PointType>(cpCoords, m_allocatorID);

      // DEBUG
      const bool has_cp_distance = xferDom.has_path("debug/cp_distance");
      auto minDist = has_cp_distance
        ? ArrayView_from_Node<double>(xferDom.fetch_existing("debug/cp_distance"), qPtCount)
        : ArrayView<double>();

      auto cp_dist = has_cp_distance
        ? (is_first ? axom::Array<double>(qPtCount, qPtCount, m_allocatorID)
                    : axom::Array<double>(minDist, m_allocatorID))
        : axom::Array<double>(0, 0, m_allocatorID);
      // END DEBUG

      if(is_first)
      {
        cp_rank.fill(-1);
        cp_idx.fill(-1);
        cp_domidx.fill(-1);
        const PointType nowhere(axom::numeric_limits<double>::signaling_NaN());
        cp_pos.fill(nowhere);
        cp_dist.fill(axom::numeric_limits<double>::signaling_NaN());
      }
      auto query_inds = cp_idx.view();
      auto query_doms = cp_domidx.view();
      auto query_ranks = cp_rank.view();
      auto query_pos = cp_pos.view();
      auto query_min_dist = cp_dist.view();

      /// Create an ArrayView in ExecSpace that is compatible with queryPts
      PointArray execPoints(queryPts, m_allocatorID);
      auto query_pts = execPoints.view();

      if(hasObjectPoints)
      {
        // Get a device-useable iterator
        auto it = m_bvh->getTraverser();
        const int rank = m_rank;

        axom::Array<double> sqDistThresh_host(1,
                                              1,
                                              axom::execution_space<axom::SEQ_EXEC>::allocatorID());
        sqDistThresh_host[0] = m_sqDistanceThreshold;
        axom::Array<double> sqDistThresh_device =
          axom::Array<double>(sqDistThresh_host, m_allocatorID);
        auto sqDistThresh_device_view = sqDistThresh_device.view();

        auto ptCoordsView = m_objectPtCoords.view();
        auto ptDomainIdsView = m_objectPtDomainIds.view();

        {
          AXOM_ANNOTATE_SCOPE("ComputeClosestPoints");
          axom::for_all<ExecSpace>(
            qPtCount,
            AXOM_LAMBDA(std::int32_t idx) mutable {
              PointType qpt = query_pts[idx];

              MinCandidate curr_min {};
              // Preset cur_min to the closest point found so far.
              if(query_ranks[idx] >= 0)
              {
                curr_min.sqDist = squared_distance(qpt, query_pos[idx]);
                curr_min.pointIdx = query_inds[idx];
                curr_min.domainIdx = query_doms[idx];
                curr_min.rank = query_ranks[idx];
              }

              auto checkMinDist = [&](std::int32_t current_node, const std::int32_t* leaf_nodes) {
                const int candidate_point_idx = leaf_nodes[current_node];
                const int candidate_domain_idx = ptDomainIdsView[candidate_point_idx];
                const PointType candidate_pt = ptCoordsView[candidate_point_idx];
                const double sq_dist = squared_distance(qpt, candidate_pt);

                if(sq_dist < curr_min.sqDist)
                {
                  curr_min.sqDist = sq_dist;
                  curr_min.pointIdx = candidate_point_idx;
                  curr_min.domainIdx = candidate_domain_idx;
                  curr_min.rank = rank;
                }
              };

              auto traversePredicate = [&](const PointType& p, const BoxType& bb) -> bool {
                auto sqDist = squared_distance(p, bb);
                return sqDist <= curr_min.sqDist && sqDist <= sqDistThresh_device_view[0];
              };

              // Traverse the tree, searching for the point with minimum distance.
              it.traverse_tree(qpt, checkMinDist, traversePredicate);

              // If modified, update the fields that changed
              if(curr_min.rank == rank)
              {
                query_inds[idx] = curr_min.pointIdx;
                query_doms[idx] = curr_min.domainIdx;
                query_ranks[idx] = curr_min.rank;
                query_pos[idx] = ptCoordsView[curr_min.pointIdx];

                //DEBUG
                if(has_cp_distance)
                {
                  query_min_dist[idx] = sqrt(curr_min.sqDist);
                }
              }
            });
        }
      }

      axom::copy(cpIndexes.data(), query_inds.data(), cpIndexes.size() * sizeof(axom::IndexType));
      axom::copy(cpDomainIndexes.data(),
                 query_doms.data(),
                 cpDomainIndexes.size() * sizeof(axom::IndexType));
      axom::copy(cpRanks.data(), query_ranks.data(), cpRanks.size() * sizeof(axom::IndexType));
      axom::copy(cpCoords.data(), query_pos.data(), cpCoords.size() * sizeof(PointType));

      // DEBUG
      if(has_cp_distance)
      {
        axom::copy(minDist.data(), query_min_dist.data(), minDist.size() * sizeof(double));
      }
    }

    // Data has now been initialized
    if(is_first)
    {
      xferNode.remove_child("is_first");
    }
  }

private:
  /*!
    @brief Object point coordindates array.

    Points from all local object mesh domains are flattened here.
  */
  PointArray m_objectPtCoords;

  axom::Array<axom::IndexType> m_objectPtDomainIds;

  /*!  @brief Object partition bounding boxes, one per rank.
    All are in physical space, not index space.
  */
  BoxArray m_objectPartitionBbs;

  std::unique_ptr<BVHTreeType> m_bvh;
};  // DistributedClosestPointExec

}  // namespace internal

}  // end namespace quest
}  // end namespace axom

#endif  //  QUEST_DISTRIBUTED_CLOSEST_POINT_IMPL_H_
