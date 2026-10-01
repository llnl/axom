// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/config.hpp"
#include "axom/core.hpp"
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
#include <algorithm>
#include <functional>
#include <optional>

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
    , m_mpiAllocatorID(MALLOC_ALLOCATOR_ID)
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

  void setMpiAllocatorID(int mpiAllocatorID)
  {
    SLIC_ASSERT(mpiAllocatorID != axom::INVALID_ALLOCATOR_ID);
    // TODO: If appropriate, how to check for compatibility with runtime policy?
    m_mpiAllocatorID = mpiAllocatorID;
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
   @brief Enables dynamic filtering of object ranks using current closest distances.
  */
  void setDynamicDistanceFiltering(bool on) { m_dynamicDistanceFiltering = on; }

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

  virtual void computeClosestPoints(conduit::Node& queryMesh,
                                    const std::string& topologyName) const = 0;

protected:
  int m_allocatorID;
  int m_mpiAllocatorID;
  bool m_isVerbose;

  MPI_Comm m_mpiComm;
  int m_rank;
  int m_nranks;

  double m_sqDistanceThreshold;

  bool m_dynamicDistanceFiltering = true;

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
 * \class DCPTransferNode
 *
 * \brief Holds data for a query mesh in a contiguous buffer.
 */
template <int NDIMS, typename ExecSpace>
struct DCPTransferNode
{
private:
  using PointType = primal::Point<double, NDIMS>;
  using BoxType = primal::BoundingBox<double, NDIMS>;

  IndexType computeSize(IndexType numPoints) const
  {
    constexpr IndexType PerNodeSize = sizeof(Metadata);
    constexpr IndexType PerPointSize = sizeof(PointType) * 2 + 3 * sizeof(IndexType);

    return PerNodeSize + PerPointSize * numPoints;
  }

  void Allocate(IndexType numPoints, int allocatorID)
  {
    IndexType total_size = computeSize(numPoints);
    buffer = axom::Array<std::uint8_t>(total_size, total_size, allocatorID);
    UpdateView();
  }

  void UpdateView()
  {
    int numPoints = metadata.numPoints;

    auto* data = buffer.data() + sizeof(Metadata);
    points = axom::ArrayView<PointType>(reinterpret_cast<PointType*>(data), numPoints);
    data += sizeof(PointType) * numPoints;

    cp_coords = axom::ArrayView<PointType>(reinterpret_cast<PointType*>(data), numPoints);
    data += sizeof(PointType) * numPoints;

    cp_index = axom::ArrayView<IndexType>(reinterpret_cast<IndexType*>(data), numPoints);
    data += sizeof(IndexType) * numPoints;

    cp_rank = axom::ArrayView<IndexType>(reinterpret_cast<IndexType*>(data), numPoints);
    data += sizeof(IndexType) * numPoints;

    cp_domain_index = axom::ArrayView<IndexType>(reinterpret_cast<IndexType*>(data), numPoints);
  }

public:
  struct Metadata
  {
    int homeRank {-1};
    int dims {0};
    int numPoints {0};
    bool isFirst {true};
    BoxType aabb;
  } metadata;

  axom::ArrayView<PointType> points;
  axom::ArrayView<PointType> cp_coords;
  axom::ArrayView<IndexType> cp_index;
  axom::ArrayView<IndexType> cp_rank;
  axom::ArrayView<IndexType> cp_domain_index;

  axom::Array<std::uint8_t> buffer;

  /*!
   * \brief Copy constructor for transfer node.
   *
   *  This constructor copies the data from a transfer node in one allocator
   *  into memory for another allocator.
   *
   * \param [in] from the node to copy from
   * \param [in] allocatorID the allocator to copy to
   */
  DCPTransferNode(const DCPTransferNode& from, int allocatorID) : metadata(from.metadata)
  {
    IndexType validSize = computeSize(metadata.numPoints);
    auto validBuffer = from.buffer.view().subspan(0, validSize);
    buffer = axom::Array<std::uint8_t>(validBuffer, allocatorID);
    UpdateView();
  }

  /*!
   * \brief Allocate an empty transfer node with backing memroy.
   *
   *  This constructor is used to allocate empty memory for incoming transfer node
   *  receives.
   *
   * \param [in] numPoints the max number of points to fit
   * \param [in] allocatorID the allocator to allocate in
   */
  DCPTransferNode(IndexType numPoints, int allocatorID) { Allocate(numPoints, allocatorID); }

  DCPTransferNode(DCPTransferNode&&) noexcept = default;
  DCPTransferNode& operator=(DCPTransferNode&&) noexcept = default;

  void Isend(int dst, int tag, MPI_Comm comm, MPI_Request& request)
  {
    IndexType total_size = computeSize(metadata.numPoints);
    axom::copy(buffer.data(), reinterpret_cast<std::uint8_t*>(&metadata), sizeof(Metadata));

    const int mpi_err =
      MPI_Isend(buffer.data(), static_cast<int>(total_size), MPI_BYTE, dst, tag, comm, &request);
    SLIC_ASSERT(mpi_err == MPI_SUCCESS);
    AXOM_UNUSED_VAR(mpi_err);
  }

  void Irecv(int src, int tag, MPI_Comm comm, MPI_Request& request)
  {
    const int mpi_err =
      MPI_Irecv(buffer.data(), static_cast<int>(buffer.size()), MPI_BYTE, src, tag, comm, &request);
    SLIC_ASSERT(mpi_err == MPI_SUCCESS);
    AXOM_UNUSED_VAR(mpi_err);
  }

  void UpdateMetadataFromBuffer()
  {
    axom::copy(reinterpret_cast<std::uint8_t*>(&metadata), buffer.data(), sizeof(Metadata));
    UpdateView();
  }

  DCPTransferNode() = default;
  DCPTransferNode(const DCPTransferNode&) = delete;
  DCPTransferNode& operator=(const DCPTransferNode&) = delete;

  /*!
   * \brief Copy a query mesh in Conduit Blueprint format into transfer buffer.
   *
   * \param [in] queryNode the blueprint mesh to copy from
   * \param [in] topologyName name of the coordinate set to use in the blueprint mesh
   * \param [in] rank current rank of this processor
   * \param [in] allocatorID where to allocate the storage for the transfer node
   */
  void copyFromConduitNode(conduit::Node& queryNode,
                           const std::string& topologyName,
                           int rank,
                           int allocatorID)
  {
    const bool isMultidomain = conduit::blueprint::mesh::is_multi_domain(queryNode);
    const auto domainCount = conduit::blueprint::mesh::number_of_domains(queryNode);
    metadata.homeRank = rank;
    metadata.isFirst = true;
    metadata.dims = NDIMS;
    metadata.numPoints = 0;
    for(conduit::index_t domainNum = 0; domainNum < domainCount; ++domainNum)
    {
      auto& queryDom = isMultidomain ? queryNode.child(domainNum) : queryNode;
      const std::string coordsetName =
        queryDom.fetch_existing(axom::fmt::format("topologies/{}/coordset", topologyName)).as_string();
      conduit::Node& queryCoords = queryDom.fetch_existing(fmt::format("coordsets/{}", coordsetName));
      conduit::Node& queryCoordsValues = queryCoords.fetch_existing("values");
      const int dim = internal::extractDimension(queryCoordsValues);
      const int qPtCount = internal::extractSize(queryCoordsValues);
      SLIC_ASSERT(dim == metadata.dims);
      metadata.numPoints += qPtCount;
    }

    Allocate(metadata.numPoints, allocatorID);

    axom::IndexType pointOffset = 0;
    for(conduit::index_t domainNum = 0; domainNum < domainCount; ++domainNum)
    {
      auto& queryDom = isMultidomain ? queryNode.child(domainNum) : queryNode;
      const std::string coordsetName =
        queryDom.fetch_existing(axom::fmt::format("topologies/{}/coordset", topologyName)).as_string();
      conduit::Node& values =
        queryDom.fetch_existing(fmt::format("coordsets/{}/values", coordsetName));
      const int qPtCount = internal::extractSize(values);
      copyFromConduitPoints(values, pointOffset);
      pointOffset += qPtCount;
    }
  }

  /*!
   * \brief Write out query results in Conduit Blueprint format
   *
   * \param [out] queryNode the blueprint mesh to write to
   * \param [in] topologyName name of the coordinate set used in the blueprint mesh
   * \param [in] allocatorID where to allocate the storage for the transfer node
   */
  void copyToConduitNode(conduit::Node& queryNode,
                         const std::string& topologyName,
                         int allocatorID,
                         bool outputRank,
                         bool outputIndex,
                         bool outputDomainIndex,
                         bool outputDistance,
                         bool outputCoords) const
  {
    const bool isMultidomain = conduit::blueprint::mesh::is_multi_domain(queryNode);
    const auto domainCount = conduit::blueprint::mesh::number_of_domains(queryNode);
    axom::IndexType pointOffset = 0;
    for(conduit::index_t domainNum = 0; domainNum < domainCount; ++domainNum)
    {
      auto& queryDom = isMultidomain ? queryNode.child(domainNum) : queryNode;
      conduit::Node& fields = queryDom.fetch_existing("fields");
      const std::string coordsetName =
        queryDom.fetch_existing(axom::fmt::format("topologies/{}/coordset", topologyName)).as_string();
      const conduit::Node& values =
        queryDom.fetch_existing(fmt::format("coordsets/{}/values", coordsetName));
      const axom::IndexType qPtCount = internal::extractSize(values);

      conduit::Node genericHeaders;
      genericHeaders["association"] = "vertex";
      genericHeaders["topology"] = topologyName;

      if(outputRank)
      {
        auto& dst = fields["cp_rank"];
        dst.set_node(genericHeaders);
        dst["values"].set(cp_rank.data() + pointOffset, qPtCount);
      }

      if(outputIndex)
      {
        auto& dst = fields["cp_index"];
        dst.set_node(genericHeaders);
        dst["values"].set(cp_index.data() + pointOffset, qPtCount);
      }

      if(outputDomainIndex)
      {
        auto& dst = fields["cp_domain_index"];
        dst.set_node(genericHeaders);
        dst["values"].set(cp_domain_index.data() + pointOffset, qPtCount);
      }

      if(outputDistance)
      {
        // Distance to closest point is a derived quantity
        auto points_v = points;
        auto cp_coords_v = cp_coords;
        axom::Array<double> cp_distance(qPtCount, qPtCount, allocatorID);
        auto cp_distance_v = cp_distance.view();

        axom::for_all<ExecSpace>(qPtCount, [=] AXOM_HOST_DEVICE(axom::IndexType idx) {
          double squared_dist = axom::primal::squared_distance(points_v[idx + pointOffset],
                                                               cp_coords_v[idx + pointOffset]);
          cp_distance_v[idx] = sqrt(squared_dist);
        });

        auto& dst = fields["cp_distance"];
        dst.set_node(genericHeaders);
        dst["values"].set(cp_distance.data(), qPtCount);
      }

      if(outputCoords)
      {
        auto& dst = fields["cp_coords"];
        dst.set_node(genericHeaders);
        auto& dstValues = dst["values"];

        copyToConduitPoints(pointOffset, qPtCount, dstValues);
      }
      pointOffset += qPtCount;
    }
  }

  /*!
   * \brief Helper to read points from a Conduit coordinate set.
   *
   *  Coordinates may be interleaved in memory or stored as separate arrays for each
   *  dimension.
   *
   * \param [in] components conduit node the coords are stored in
   * \param [in] pointOffset offset in the transfer node for the current domain
   */
  void copyFromConduitPoints(conduit::Node& components, axom::IndexType pointOffset) const
  {
    const int dim = NDIMS;
    const int qPtCount = internal::extractSize(components);
    SLIC_ASSERT(dim == metadata.dims);

    auto* dst = reinterpret_cast<double*>(points.data() + pointOffset);
    const bool interleavedSrc = conduit::blueprint::mcarray::is_interleaved(components);
    if(interleavedSrc)
    {
      axom::copy(dst,
                 internal::getPointer<double>(components.child(0)),
                 dim * qPtCount * sizeof(double));
    }
    else
    {
      // Copy from component-wise source to the interleaved destination.
      for(int d = 0; d < dim; ++d)
      {
        auto src = components.child(d).as_float64_array();
        for(int i = 0; i < qPtCount; ++i)
        {
          dst[i * dim + d] = src[i];
        }
      }
    }
  }

  /*!
   * \brief Helper to write out result closest points into Conduit Blueprint.
   *
   *  Points are always stored as component-wise.
   *
   * \param [in] pointOffset offset in the transfer node for the current domain
   * \param [in] qptCount number of points from the transfer node to write out
   * \param [out] components conduit node to write out to
   */
  void copyToConduitPoints(IndexType pointOffset, IndexType qPtCount, conduit::Node& components) const
  {
    const int dim = NDIMS;
    components.reset();

    const auto* interleaved = reinterpret_cast<const double*>(cp_coords.data() + pointOffset);
    // Copy from 1D-interleaved src to component-wise dst.
    for(int d = 0; d < dim; ++d)
    {
      const double* src = interleaved + d;
      auto& dstNode = components.append();
      dstNode.set_dtype(conduit::DataType::float64(qPtCount));
      double* dst = dstNode.as_float64_ptr();
      for(int i = 0; i < qPtCount; ++i)
      {
        dst[i] = src[i * dim];
      }
    }
  }

  /// Wait for a receive or one or more non-blocking sends to finish.
  static void WaitMPIRequests(std::deque<std::pair<DCPTransferNode, MPI_Request>>& isendRequests,
                              MPI_Request& irecvRequest)
  {
    AXOM_ANNOTATE_SCOPE("WaitMPIRequests");
    std::vector<MPI_Request> reqs;
    reqs.reserve(isendRequests.size() + 1);
    reqs.push_back(irecvRequest);
    for(auto const& isr : isendRequests)
    {
      reqs.push_back(isr.second);
    }

    int inCount = static_cast<int>(reqs.size());
    if(irecvRequest == MPI_REQUEST_NULL)
    {
      // No outstanding receive requests, just wait on all completed sends.
      const int mpi_err = MPI_Waitall(inCount, reqs.data(), MPI_STATUSES_IGNORE);
      SLIC_ASSERT(mpi_err == MPI_SUCCESS);
      AXOM_UNUSED_VAR(mpi_err);

      // Free all allocated send buffers.
      isendRequests.clear();
      return;
    }

    std::vector<int> finished_requests;
    // Keep waiting until a query is received.
    while(reqs[0] != MPI_REQUEST_NULL)
    {
      std::vector<int> indices(reqs.size(), -1);
      int outCount = 0;
      const int mpi_err =
        MPI_Waitsome(inCount, reqs.data(), &outCount, indices.data(), MPI_STATUSES_IGNORE);
      SLIC_ASSERT(mpi_err == MPI_SUCCESS);
      AXOM_UNUSED_VAR(mpi_err);

      indices.resize(outCount);
      finished_requests.insert(finished_requests.end(), indices.begin(), indices.end());
    }

    // MPI does not guarantee indices are in order.
    // Removing in descending index order ensures that moving the last element
    // cannot invalidate an index that remains to be processed.
    std::sort(finished_requests.begin(), finished_requests.end(), std::greater<int>());

    for(int request_index : finished_requests)
    {
      if(request_index == 0)
      {
        // Irecv special case: just set the irecv request to MPI_REQUEST_NULL
        irecvRequest = MPI_REQUEST_NULL;
        continue;
      }
      // Remove completed send requests to free memory
      int completed_send_index = request_index - 1;
      if(completed_send_index + 1 != isendRequests.size())
      {
        isendRequests[completed_send_index] = std::move(isendRequests.back());
      }
      isendRequests.pop_back();
    }
  }
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

  using TransferNode = DCPTransferNode<NDIMS, ExecSpace>;

public:
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
  BoxType computeMeshBoundingBox(const TransferNode& xferNode) const
  {
    BoxType rval;
    SLIC_ASSERT(xferNode.metadata.dims == DIM);
    for(const auto& p : xferNode.points)
    {
      rval.addPoint(p);
    }

    return rval;
  }

  /**
   * \brief Implementation of the user-facing
   * DistributedClosestPoint::computeClosestPoints() method.
   *
   * We use non-blocking sends for performance and deadlock avoidance.
   * The worst case could incur nranks^2 sends.  To avoid excessive
   * buffer usage, we wait for receives and sends together, freeing send
   * buffers as their requests complete.
   */
  void computeClosestPoints(conduit::Node& queryMesh, const std::string& topologyName) const override
  {
    SLIC_ASSERT_MSG(m_bvh, "BVH tree must be initialized before calling 'computeClosestPoints");

    // MPI guarantees MPI_TAG_UB is at least 32767, so use a tag below that
    constexpr int tag = 9873;

    int maxParticlesToRecv = 0;
    int remainingRecvs = 0;
    int fullXferRecvs = 0;

    std::deque<std::pair<TransferNode, MPI_Request>> isendRequests;

    TransferNode localXferNode;

    {
      TransferNode& xferNode = localXferNode;
      // create conduit Node containing data that has to xfer between ranks.
      // The node will be mostly empty if there are no domains on this rank
      xferNode.copyFromConduitNode(queryMesh, topologyName, m_rank, m_allocatorID);

      BoxType myQueryBb = computeMeshBoundingBox(xferNode);
      xferNode.metadata.aabb = myQueryBb;

      // Get maximum number of particles from any query rank.
      // This allows us to pre-alloocate the required receive buffer size,
      // which may allow corresponding MPI_Sends to be performed "eagerly."
      maxParticlesToRecv = xferNode.metadata.numPoints;
      {
        int mpi_err =
          MPI_Allreduce(MPI_IN_PLACE, &maxParticlesToRecv, 1, MPI_INT, MPI_MAX, m_mpiComm);
        SLIC_ASSERT(mpi_err == MPI_SUCCESS);
        AXOM_UNUSED_VAR(mpi_err);
      }

      double currentMaxSqDistance = computeLocalClosestPoints(xferNode);

      const auto myObjectBb = m_objectPartitionBbs[m_rank];
      axom::Array<int> send_counts(m_nranks);
      send_counts.fill(0);
      for(int r = 0; r < m_nranks; ++r)
      {
        if(r != m_rank)
        {
          const auto& otherQueryBb = m_objectPartitionBbs[r];
          if(is_statically_eligible(otherQueryBb, myQueryBb))
          {
            send_counts[r]++;
          }
        }
      }

      // We compute the receive counts by summing the send counts for a given
      // object rank from all query processors.
      {
        int mpi_err =
          MPI_Reduce_scatter_block(send_counts.data(), &remainingRecvs, 1, MPI_INT, MPI_SUM, m_mpiComm);
        SLIC_ASSERT(mpi_err == MPI_SUCCESS);
        AXOM_UNUSED_VAR(mpi_err);
      }

      /*
        Send local query mesh to next rank with close-enough object
        partition, if any.  Increase remainingRecvs, because this data
        will come back.
      */
      int firstRecipForMyQuery = next_recipient(xferNode, currentMaxSqDistance, tag, isendRequests);
      if(m_nranks == 1)
      {
        SLIC_ASSERT(firstRecipForMyQuery == -1);
      }

      if(firstRecipForMyQuery != -1)
      {
        if(xferNode.buffer.getAllocatorID() == m_mpiAllocatorID)
        {
          // Transfer node is allocated in MPI communication space, don't allocate
          // a staging buffer.
          isendRequests.emplace_back(std::move(xferNode), MPI_Request {});
        }
        else
        {
          // Copy node to selected MPI memory pool, then do the MPI communication.
          TransferNode sendNode(xferNode, m_mpiAllocatorID);
          isendRequests.emplace_back(std::move(sendNode), MPI_Request {});
        }
        auto& req = isendRequests.back();
        req.first.Isend(firstRecipForMyQuery, tag, m_mpiComm, req.second);
        ++remainingRecvs;
      }
    }

#if defined(AXOM_USE_UMPIRE)
    bool gpuAwareMpi =
      axom::isDeviceAllocator(m_mpiAllocatorID) && axom::isDeviceAllocator(m_allocatorID);
#else
    bool gpuAwareMpi = false;
#endif

    // Allocate our receive node persistently. This may improve performance by
    // avoiding repeated memory registrations in the MPI implementation.
    // TODO: test this on the CPU
    TransferNode recvXferNode(maxParticlesToRecv, m_mpiAllocatorID);

    const int totalExpectedRecvs = remainingRecvs;
    while(remainingRecvs > 0)
    {
      if(m_isVerbose)
      {
        if(m_dynamicDistanceFiltering)
        {
          const int skipTokenRecvs = totalExpectedRecvs - remainingRecvs - fullXferRecvs;
          SLIC_INFO(fmt::format("=======  receives: remaining={}, full={}, skip={} =======",
                                remainingRecvs,
                                fullXferRecvs,
                                skipTokenRecvs));
        }
        else
        {
          SLIC_INFO(fmt::format("=======  {} receives remaining =======", remainingRecvs));
        }
      }

      // Receive the next xferNode
      MPI_Request recv_req = MPI_REQUEST_NULL;
      recvXferNode.Irecv(MPI_ANY_SOURCE, tag, m_mpiComm, recv_req);

      // Wait for receive to complete
      TransferNode::WaitMPIRequests(isendRequests, recv_req);
      recvXferNode.UpdateMetadataFromBuffer();

      --remainingRecvs;
      const int homeRank = recvXferNode.metadata.homeRank;
      if(homeRank < 0)
      {
        continue;
      }

      TransferNode xferNode;
      {
        AXOM_ANNOTATE_SCOPE("CopyQueryFromMPI");
        // We need to copy the received node from the persistent MPI allocation
        // into a temporary node.
        xferNode = TransferNode(recvXferNode, m_allocatorID);
        if(gpuAwareMpi)
        {
          // Device-to-device copies are treated as asynchronous with respect to
          // the host. Needed on HIP since we access metadata from the host.
          axom::synchronize<ExecSpace>();
        }
      }
      ++fullXferRecvs;
      if(homeRank == m_rank)
      {
        localXferNode = std::move(xferNode);
      }
      else
      {
        double currentMaxSqDistance = computeLocalClosestPoints(xferNode);

        int nextRecipient = next_recipient(xferNode, currentMaxSqDistance, tag, isendRequests);
        SLIC_ASSERT(nextRecipient != -1);
        if(gpuAwareMpi)
        {
          // Transfer node is allocated in MPI communication space, don't allocate
          // a staging buffer.
          isendRequests.emplace_back(std::move(xferNode), MPI_Request {});
        }
        else
        {
          AXOM_ANNOTATE_SCOPE("CopyQueryToMPI");
          // Copy node to selected MPI memory pool, then do the MPI communication.
          TransferNode sendNode(xferNode, m_mpiAllocatorID);
          isendRequests.emplace_back(std::move(sendNode), MPI_Request {});
        }
        auto& req = isendRequests.back();
        req.first.Isend(nextRecipient, tag, m_mpiComm, req.second);
      }

    }  // remainingRecvs loop

    // Complete remaining non-blocking sends.
    MPI_Request recv_req = MPI_REQUEST_NULL;
    TransferNode::WaitMPIRequests(isendRequests, recv_req);
    SLIC_ASSERT(isendRequests.empty());

    {
      localXferNode.copyToConduitNode(queryMesh,
                                      topologyName,
                                      m_allocatorID,
                                      m_outputRank,
                                      m_outputIndex,
                                      m_outputDomainIndex,
                                      m_outputDistance,
                                      m_outputCoords);
    }

    MPI_Barrier(m_mpiComm);
    slic::flushStreams();
  }

  /*!
   * \brief Check whether static rank bounding boxes are close enough to search.
   *
   * \param [out] sqDistance Optional reference to receive the squared distance
   * between the boxes.
   */
  bool is_statically_eligible(const BoxType& queryBb,
                              const BoxType& objectBb,
                              std::optional<std::reference_wrapper<double>> sqDistance = std::nullopt) const
  {
    if(!queryBb.isValid() || !objectBb.isValid())
    {
      return false;
    }

    const double sqDist = primal::squared_distance(queryBb, objectBb);
    if(sqDistance)
    {
      sqDistance->get() = sqDist;
    }
    return sqDist <= m_sqDistanceThreshold;
  }

  /*!
   * \brief Send a minimal message indicating this rank has no query data for the destination rank.
   */
  void send_skip_token(int dest,
                       int tag,
                       std::deque<std::pair<TransferNode, MPI_Request>>& isendRequests) const
  {
    TransferNode skipToken(0, m_mpiAllocatorID);
    skipToken.metadata.homeRank = -1;
    skipToken.metadata.dims = DIM;

    isendRequests.emplace_back(std::move(skipToken), MPI_Request {});
    auto& req = isendRequests.back();
    req.first.Isend(dest, tag, m_mpiComm, req.second);
  }

  /*!
   * \brief Determine the next object rank in ring order to receive the transfer node.
   *
   * Sends skip tokens to statically eligible ranks that dynamic distance
   * filtering can prove will not improve the current closest point results.
   */
  int next_recipient(const TransferNode& xferNode,
                     double currentMaxSqDistance,
                     int tag,
                     std::deque<std::pair<TransferNode, MPI_Request>>& isendRequests) const
  {
    int homeRank = xferNode.metadata.homeRank;
    const BoxType& bb = xferNode.metadata.aabb;

    for(int i = 1; i < m_nranks; ++i)
    {
      int maybeNextRecip = (m_rank + i) % m_nranks;
      if(maybeNextRecip == homeRank)
      {
        return maybeNextRecip;
      }
      if(double sqDistance = 0.0;
         is_statically_eligible(bb, m_objectPartitionBbs[maybeNextRecip], sqDistance))
      {
        /// Use the next recipient if that rank may be able to update one of the
        /// points in the query, otherwise skip that rank. An update is possible
        /// if the minimum distance between the bounding boxes is less than
        /// the current distance of any of the query points.
        if(sqDistance <= currentMaxSqDistance)
        {
          return maybeNextRecip;
        }
        send_skip_token(maybeNextRecip, tag, isendRequests);
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

    axom::for_all<ExecSpace>(npts, [=] AXOM_HOST_DEVICE(axom::IndexType i) {
      boxesView[i] = BoxType {pointsView[i]};
    });

    // Build bounding volume hierarchy
    bvh->setAllocatorID(m_allocatorID);
    int result = bvh->initialize(boxesView, npts);

    gatherBVHRoots();

    return (result == spin::BVH_BUILD_OK);
  }

  double computeLocalClosestPoints(TransferNode& xferNode) const
  {
    using axom::primal::squared_distance;

    // Note: There is some additional computation the first time this function
    // is called for a query node, even if the local object mesh is empty
    const bool hasObjectPoints = m_objectPtCoords.size() > 0;
    const bool is_first = xferNode.metadata.isFirst;
    if(!hasObjectPoints && !is_first)
    {
      return m_sqDistanceThreshold;
    }

    double currentMaxSqDistance = -1.0;

    {
      // Check dimension and extract the number of points
      SLIC_ASSERT(xferNode.metadata.dims == DIM);
      const int qPtCount = xferNode.metadata.numPoints;

      // Extract fields from the input node as ArrayViews
      // These are allocated with the user-specified allocator ID, which is expected
      // to be accessible from the requested execution space.
      SLIC_ASSERT(xferNode.buffer.getAllocatorID() == m_allocatorID);
      auto query_pts = xferNode.points;
      auto query_inds = xferNode.cp_index;
      auto query_doms = xferNode.cp_domain_index;
      auto query_ranks = xferNode.cp_rank;
      auto query_pos = xferNode.cp_coords;

      if(is_first)
      {
        query_ranks.fill(-1);
        query_inds.fill(-1);
        query_doms.fill(-1);
        const PointType nowhere(axom::numeric_limits<double>::signaling_NaN());
        query_pos.fill(nowhere);
      }

      if(hasObjectPoints)
      {
        auto query_order = mortonSortQueryPoints(query_pts);
        auto query_order_view = query_order.view();
        const double sqDistThreshold = m_sqDistanceThreshold;
        auto it = m_bvh->getTraverser();
        const int rank = m_rank;
        auto ptCoordsView = m_objectPtCoords.view();
        auto ptDomainIdsView = m_objectPtDomainIds.view();

        /// Update the query by finding the closest points in the local mesh
        if(m_dynamicDistanceFiltering)
        {
          /// Dynamic filtering does the same distance calculation as static but
          /// additionally calculates the furthest distance of any query point
          /// after the local update
          AXOM_ANNOTATE_SCOPE("ComputeClosestPointsDynamic");
          axom::ReduceMax<ExecSpace, double> maxSqDistance(currentMaxSqDistance);
          axom::for_all<ExecSpace>(qPtCount, [=] AXOM_HOST_DEVICE(std::int32_t sorted_idx) {
            const auto idx = query_order_view[sorted_idx];
            PointType qpt = query_pts[idx];

            MinCandidate curr_min {};
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
              return sqDist <= curr_min.sqDist && sqDist <= sqDistThreshold;
            };

            it.template traverseTreeShared<ExecSpace>(qpt, checkMinDist, traversePredicate);

            if(curr_min.rank == rank)
            {
              query_inds[idx] = curr_min.pointIdx;
              query_doms[idx] = curr_min.domainIdx;
              query_ranks[idx] = curr_min.rank;
              query_pos[idx] = ptCoordsView[curr_min.pointIdx];
            }

            maxSqDistance.max(curr_min.rank >= 0 ? curr_min.sqDist : sqDistThreshold);
          });

          currentMaxSqDistance = maxSqDistance.get();
        }
        else
        {
          AXOM_ANNOTATE_SCOPE("ComputeClosestPoints");
          axom::for_all<ExecSpace>(qPtCount, [=] AXOM_HOST_DEVICE(std::int32_t sorted_idx) {
            const auto idx = query_order_view[sorted_idx];
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
              return sqDist <= curr_min.sqDist && sqDist <= sqDistThreshold;
            };

            // Traverse the tree, searching for the point with minimum distance.
            it.template traverseTreeShared<ExecSpace>(qpt, checkMinDist, traversePredicate);

            // If modified, update the fields that changed
            if(curr_min.rank == rank)
            {
              query_inds[idx] = curr_min.pointIdx;
              query_doms[idx] = curr_min.domainIdx;
              query_ranks[idx] = curr_min.rank;
              query_pos[idx] = ptCoordsView[curr_min.pointIdx];
            }
          });
        }
      }
    }

    // Data has now been initialized
    if(is_first)
    {
      xferNode.metadata.isFirst = false;
    }

    return m_dynamicDistanceFiltering && currentMaxSqDistance >= 0.0 ? currentMaxSqDistance
                                                                     : m_sqDistanceThreshold;
  }

  /*! \brief Returns query point indices ordered by their Morton codes. */
  axom::Array<axom::IndexType> mortonSortQueryPoints(const axom::ArrayView<PointType>& queryPoints) const
  {
    IndexType queryPointCount = queryPoints.size();
    axom::Array<axom::IndexType> queryOrder(queryPointCount, queryPointCount, m_allocatorID);
    if(queryPointCount == 0)
    {
      return queryOrder;
    }

    PointType minPoint;
    PointType inverseExtent;
    for(int dim = 0; dim < DIM; ++dim)
    {
      axom::ReduceMin<ExecSpace, double> minCoord(axom::numeric_limits<double>::max());
      axom::ReduceMax<ExecSpace, double> maxCoord(axom::numeric_limits<double>::lowest());
      axom::for_all<ExecSpace>(queryPointCount, [=] AXOM_HOST_DEVICE(axom::IndexType idx) {
        minCoord.min(queryPoints[idx][dim]);
        maxCoord.max(queryPoints[idx][dim]);
      });

      minPoint[dim] = minCoord.get();
      const double extent = maxCoord.get() - minPoint[dim];
      inverseExtent[dim] = extent > 0.0 ? 1.0 / extent : 0.0;
    }

    axom::Array<std::uint32_t> mortonCodes(queryPointCount, queryPointCount, m_allocatorID);
    auto morton_codes = mortonCodes.view();
    auto query_order = queryOrder.view();
    axom::for_all<ExecSpace>(queryPointCount, [=] AXOM_HOST_DEVICE(axom::IndexType idx) {
      constexpr int bits_per_dimension = 32 / DIM;
      constexpr double coordinate_scale = 1 << bits_per_dimension;
      constexpr double coordinate_max = coordinate_scale - 1.0;

      primal::Point<std::int32_t, DIM> gridPoint;
      for(int dim = 0; dim < DIM; ++dim)
      {
        const double coordinate = (queryPoints[idx][dim] - minPoint[dim]) * inverseExtent[dim];
        gridPoint[dim] = static_cast<std::int32_t>(
          axom::utilities::clampVal(coordinate * coordinate_scale, 0.0, coordinate_max));
      }

      morton_codes[idx] = spin::convertPointToMorton<std::uint32_t>(gridPoint);
      query_order[idx] = idx;
    });

    axom::stable_sort_pairs<ExecSpace>(morton_codes, query_order);
    return queryOrder;
  }

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
