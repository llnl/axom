// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/primal.hpp"
#include "axom/quest/ShapeeMesh.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/fmt.hpp"

#include "conduit_blueprint.hpp"

namespace axom
{
namespace quest
{
#if defined(AXOM_USE_64BIT_INDEXTYPE) && !defined(AXOM_NO_INT64_T)
  #if defined(AXOM_USE_CONDUIT)
static constexpr conduit::DataType::TypeID conduitDataIdOfAxomIndexType =
  conduit::DataType::INT64_ID;
  #endif
#else
  #if defined(AXOM_USE_CONDUIT)
static constexpr conduit::DataType::TypeID conduitDataIdOfAxomIndexType =
  conduit::DataType::INT32_ID;
  #endif
#endif

constexpr int NUM_TETS_PER_HEX = 24;

ShapeeMesh::ShapeeMesh(RuntimePolicy runtimePolicy,
                       int allocatorId,
                       conduit::Node& bpMesh,
                       const std::string& topo)
  : m_runtimePolicy(runtimePolicy)
  , m_allocId(allocatorId != axom::INVALID_ALLOCATOR_ID
                ? allocatorId
                : axom::policyToDefaultAllocatorID(runtimePolicy))
  , m_bpTopo(topo.empty() ? bpMesh.fetch_existing("topologies").child(0).name()
                          : topo)
  , m_bpNodeExt(&bpMesh)
  , m_bpNodeInt()
{
  // We want unstructured topo but can accomodate structured.
  const std::string topoType = m_bpNodeExt->fetch_existing("topologies")
                                 .fetch_existing(m_bpTopo)
                                 .fetch_existing("type")
                                 .as_string();
  SLIC_ASSERT_MSG(
    topoType == "unstructured",
    "ShapeeMesh currently only works with unstructured mesh, not " + topoType +
      ".");

  const conduit::Node& topoNode =
    m_bpNodeExt->fetch_existing("topologies").fetch_existing(m_bpTopo);

  const std::string coordsetName =
    topoNode.fetch_existing("coordset").as_string();
  const conduit::Node& coordsetNode =
    m_bpNodeExt->fetch_existing("coordsets").fetch_existing(coordsetName);

  // Check for unstructured and hexahedral
  SLIC_ERROR_IF(topoNode["type"].as_string() != "unstructured",
                "topology type must be 'unstructured'");
  SLIC_ERROR_IF(topoNode["elements/shape"].as_string() != "hex",
                "element shape must be 'hex'");

  m_dim = conduit::blueprint::mesh::topology::dims(topoNode);
  SLIC_ASSERT(m_dim == 3); // Will allow 2D when we support it.

  m_cellCount = conduit::blueprint::mesh::topology::length(topoNode);

  m_vertexCount = conduit::blueprint::mesh::coordset::length(coordsetNode);

  const conduit::Node& coordsValues = coordsetNode.fetch_existing("values");
  const bool isInterleaved =
    conduit::blueprint::mcarray::is_interleaved(coordsValues);
  const int stride = isInterleaved ? m_dim : 1;
  const char* dirNames[] = {"x", "y", "z"};
  for(int d = 0; d < m_dim; ++d)
  {
    m_vertCoordsViews3D[d] =
      axom::ArrayView<const double>(coordsValues[dirNames[d]].as_double_ptr(),
                                    {m_vertexCount},
                                    stride);
  }
}

axom::ArrayView<const ShapeeMesh::TetrahedronType> ShapeeMesh::getCellsAsTets()
{
  if(m_cellsAsTets.size() != m_cellCount * NUM_TETS_PER_HEX)
  {
    computeCellsAsTets();
  }
  return m_cellsAsTets;
}

axom::ArrayView<const ShapeeMesh::HexahedronType> ShapeeMesh::getCellsAsHexes()
{
  if(m_cellsAsHexes.size() != m_cellCount)
  {
    computeCellsAsHexes();
  }
  return m_cellsAsHexes;
}

axom::ArrayView<const double> ShapeeMesh::getCellVolumes()
{
  if(m_hexVolumes.size() != m_cellCount)
  {
    computeHexVolumes();
  }
  return m_hexVolumes.view();
}

axom::ArrayView<const ShapeeMesh::BoundingBox3dType> ShapeeMesh::getCellBoundingBoxes()
{
  if(m_hexBbs.size() != m_cellCount)
  {
    computeHexBbs();
  }
  return m_hexBbs.view();
}

axom::ArrayView<const axom::IndexType, 2> ShapeeMesh::getConnectivity()
{
  if(m_connectivity.size() != m_cellCount)
  {
    computeConnectivity();
  }
  return m_connectivity;
}

void ShapeeMesh::computeCellsAsHexes()
{
  AXOM_ANNOTATE_SCOPE("ShapeeMesh::computeCellsAsHexes");
  switch(m_runtimePolicy)
  {
  case RuntimePolicy::seq:
    computeCellsAsHexesImpl<axom::SEQ_EXEC>();
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case RuntimePolicy::omp:
    computeCellsAsHexesImpl<axom::OMP_EXEC>();
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case RuntimePolicy::cuda:
    computeCellsAsHexesImpl<axom::CUDA_EXEC<256>>();
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case RuntimePolicy::hip:
    computeCellsAsHexesImpl<axom::HIP_EXEC<256>>();
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
}

void ShapeeMesh::computeCellsAsTets()
{
  AXOM_ANNOTATE_SCOPE("ShapeeMesh::computeCellsAsTets");
  switch(m_runtimePolicy)
  {
  case RuntimePolicy::seq:
    computeCellsAsTetsImpl<axom::SEQ_EXEC>();
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case RuntimePolicy::omp:
    computeCellsAsTetsImpl<axom::OMP_EXEC>();
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case RuntimePolicy::cuda:
    computeCellsAsTetsImpl<axom::CUDA_EXEC<256>>();
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case RuntimePolicy::hip:
    computeCellsAsTetsImpl<axom::HIP_EXEC<256>>();
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
}

void ShapeeMesh::computeHexVolumes()
{
  AXOM_ANNOTATE_SCOPE("ShapeeMesh::computeHexVolumes");
  switch(m_runtimePolicy)
  {
  case RuntimePolicy::seq:
    computeHexVolumesImpl<axom::SEQ_EXEC>();
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case RuntimePolicy::omp:
    computeHexVolumesImpl<axom::OMP_EXEC>();
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case RuntimePolicy::cuda:
    computeHexVolumesImpl<axom::CUDA_EXEC<256>>();
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case RuntimePolicy::hip:
    computeHexVolumesImpl<axom::HIP_EXEC<256>>();
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
}

void ShapeeMesh::computeHexBbs()
{
  AXOM_ANNOTATE_SCOPE("ShapeeMesh::computeHexBoundingBoxes");
  switch(m_runtimePolicy)
  {
  case RuntimePolicy::seq:
    computeHexBbsImpl<axom::SEQ_EXEC>();
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case RuntimePolicy::omp:
    computeHexBbsImpl<axom::OMP_EXEC>();
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case RuntimePolicy::cuda:
    computeHexBbsImpl<axom::CUDA_EXEC<256>>();
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case RuntimePolicy::hip:
    computeHexBbsImpl<axom::HIP_EXEC<256>>();
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
}

void ShapeeMesh::computeConnectivity()
{
  SLIC_ASSERT(m_dim == 3);  // 2D support not done yet.

  constexpr int NUM_VERTS_PER_HEX = 8;

  const conduit::Node& topoNode =
    m_bpNodeExt->fetch_existing("topologies").fetch_existing(m_bpTopo);
  const auto& connNode = topoNode.fetch_existing("elements/connectivity");
  SLIC_ERROR_IF(connNode.dtype().id() != conduitDataIdOfAxomIndexType,
                "IntersectionShaper error: connectivity data type must be "
                "axom::IndexType.");
  const auto* connPtr = static_cast<const axom::IndexType*>(connNode.data_ptr());
  m_connectivity = axom::ArrayView<const axom::IndexType, 2> {connPtr,
                                                              m_cellCount,
                                                              NUM_VERTS_PER_HEX};
}

template <typename ExecSpace>
void ShapeeMesh::computeCellsAsHexesImpl()
{
  constexpr int NUM_VERTS_PER_HEX = 8;
  constexpr int NDIM = 3;

  SLIC_ASSERT(m_dim == NDIM);  // or we shouldn't be here.

  auto vertexCoords = getVertexCoords3D();
  const auto& vX = vertexCoords[0];
  const auto& vY = vertexCoords[1];
  const auto& vZ = vertexCoords[2];

  axom::ArrayView<const IndexType, 2> connView = getConnectivity();

  m_cellsAsHexes = axom::Array<HexahedronType>(ArrayOptions::Uninitialized(),
                                               m_cellCount,
                                               m_cellCount,
                                               m_allocId);
  axom::ArrayView<HexahedronType> cellsAsHexesView = m_cellsAsHexes.view();

  constexpr double ZERO_THRESHOLD = 1.e-10;

  axom::for_all<ExecSpace>(
    m_cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
      // Set each hexahedral element vertices
      auto cellVertIds = connView[cellId];
      auto& hex = cellsAsHexesView[cellId];

      for(int vi = 0; vi < NUM_VERTS_PER_HEX; ++vi)
      {
        int vertIndex = cellVertIds[vi];
        primal::Point3D vCoords {vX[vertIndex], vY[vertIndex], vZ[vertIndex]};

        // Snap coordinates to zero.
        for(int d = 0; d < NDIM; ++d)
        {
          if(axom::utilities::isNearlyEqual(vCoords[d], 0.0, ZERO_THRESHOLD))
          {
            vCoords[d] = 0;
          }
        }

        hex[vi] = vCoords;
      }
    });  // end of loop to initialize hexahedral elements and bounding boxes
}

template <typename ExecSpace>
void ShapeeMesh::computeCellsAsTetsImpl()
{
  constexpr int NUM_TETS_PER_HEX = primal::Hexahedron<double, 3>::NUM_TRIANGULATE;

  SLIC_ASSERT(m_dim == 3);  // or we shouldn't be here.

  using TetsInHex = axom::StackArray<TetrahedronType, NUM_TETS_PER_HEX>;

  m_cellsAsTets = axom::Array<TetrahedronType>(ArrayOptions::Uninitialized(),
                                               m_cellCount,
                                               m_cellCount,
                                               m_allocId);
  auto cellsAsTetsView = m_cellsAsTets.view();

  auto cellsAsHexesView = getCellsAsHexes();

  axom::for_all<ExecSpace>(
    m_cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
      const auto& hex = cellsAsHexesView[cellId];

      TetsInHex tetsInHex;
      hex.triangulate(tetsInHex);

      for(int vi = 0; vi < NUM_TETS_PER_HEX; ++vi)
      {
        axom::IndexType tetId = vi + cellId * NUM_TETS_PER_HEX;
        cellsAsTetsView[tetId] = tetsInHex[vi];
      }
    });
}

template <typename ExecSpace>
void ShapeeMesh::computeHexVolumesImpl()
{
  m_hexVolumes = axom::Array<double>(ArrayOptions::Uninitialized(),
                                     m_cellCount,
                                     m_cellCount,
                                     m_allocId);

  auto cellsAsHexes = getCellsAsHexes();

  auto hexVolumesView = m_hexVolumes.view();
  axom::for_all<ExecSpace>(
    m_cellCount,
    AXOM_LAMBDA(axom::IndexType i) {
      hexVolumesView[i] = cellsAsHexes[i].volume();
    });
}

template <typename ExecSpace>
void ShapeeMesh::computeHexBbsImpl()
{
  m_hexBbs = axom::Array<BoundingBox3dType>(ArrayOptions::Uninitialized(),
                                            m_cellCount,
                                            m_cellCount,
                                            m_allocId);

  auto cellsAsHexes = getCellsAsHexes();

  auto hexBbsView = m_hexBbs.view();
  axom::for_all<ExecSpace>(
    m_cellCount,
    AXOM_LAMBDA(axom::IndexType i) {
      hexBbsView[i] = primal::compute_bounding_box<double, 3>(cellsAsHexes[i]);
    });
}

}  // end namespace quest
}  // end namespace axom
