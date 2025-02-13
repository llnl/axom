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
static constexpr conduit::DataType::TypeID conduitDataIdOfAxomIndexType =
  conduit::DataType::INT64_ID;
#else
static constexpr conduit::DataType::TypeID conduitDataIdOfAxomIndexType =
  conduit::DataType::INT32_ID;
#endif

constexpr int NUM_TETS_PER_HEX = 24;

ShapeeMesh::ShapeeMesh(RuntimePolicy runtimePolicy,
                       int allocatorId,
                       conduit::Node& bpMesh,
                       const std::string& topoName,
                       const std::string& matsetName)
  : m_runtimePolicy(runtimePolicy)
  , m_allocId(allocatorId != axom::INVALID_ALLOCATOR_ID
                ? allocatorId
                : axom::policyToDefaultAllocatorID(runtimePolicy))
  , m_bpTopo(topoName.empty() && bpMesh["topologies"].number_of_children() > 0
               ? bpMesh["topologies"].child(0).name()
               : topoName)
  , m_bpMatset(matsetName.empty() && bpMesh["matsets"].number_of_children() > 0
                 ? bpMesh.fetch("matsets").child(0).name()
                 : matsetName)
  , m_bpNodeExt(&bpMesh)
  , m_bpNodeInt()
{
  SLIC_ERROR_IF(
    m_bpTopo.empty(),
    "Topology name was not provided, and no default topology was found.");
  SLIC_ERROR_IF(
    m_bpMatset.empty(),
    "Matset name was not provided, and no default matset was found.");

  // We want unstructured topo but can accomodate structured.
  const std::string topoType = m_bpNodeExt->fetch_existing("topologies")
                                 .fetch_existing(m_bpTopo)
                                 .fetch_existing("type")
                                 .as_string();
  SLIC_ERROR_IF(topoType != "unstructured",
                "ShapeeMesh currently only works with unstructured mesh, not " +
                  topoType + ".");

  const conduit::Node& topoNode =
    m_bpNodeExt->fetch_existing("topologies").fetch_existing(m_bpTopo);
  const std::string coordsetName =
    topoNode.fetch_existing("coordset").as_string();
  const conduit::Node& coordsetNode =
    m_bpNodeExt->fetch_existing("coordsets").fetch_existing(coordsetName);
  conduit::Node& matsetNode = m_bpNodeExt->fetch("matsets").fetch(m_bpMatset);

  if(!matsetNode.has_child("topology"))
  {
    matsetNode.fetch("topology").set_string(m_bpTopo);
  }

  // Input checks.
  SLIC_ERROR_IF(topoNode["type"].as_string() != "unstructured",
                "topology type must be 'unstructured'");
  SLIC_ERROR_IF(topoNode["elements/shape"].as_string() != "hex",
                "element shape must be 'hex'");
  SLIC_ERROR_IF(matsetNode["topology"].as_string() != m_bpTopo,
                "matset's topology doesn't match the specified topology");

  m_dim = conduit::blueprint::mesh::topology::dims(topoNode);
  SLIC_ASSERT(m_dim == 3);  // Will allow 2D when we support it.

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

bool ShapeeMesh::isValidForShaping(std::string& whyNot) const
{
  bool rval = true;

  /*
    Check Blueprint-validity.
    Conduit's verify should work even if m_bpNodeExt has array data on
    devices.  The verification doesn't dereference array data.
    If this changes in the future, this code must adapt.
  */
  conduit::Node info;
  rval = conduit::blueprint::mesh::verify(*m_bpNodeExt, info);

  if(rval)
  {
    std::string topoType =
      m_bpNodeExt->fetch("topologies")[m_bpTopo]["type"].as_string();
    rval = topoType == "unstructured";
    info[0].set_string("Topology is not unstructured.");
  }

  if(rval)
  {
    std::string elemShape =
      m_bpNodeExt->fetch("topologies")[m_bpTopo]["elements/shape"].as_string();
    rval = elemShape == "hex";
    info[0].set_string("Topology elements are not hex.");
  }

  whyNot = info.to_summary_string();

  return rval;
}

void ShapeeMesh::setMatsetFromVolume(const std::string& materialName,
                                     const axom::ArrayView<double>& volumes,
                                     bool isFraction)
{
  conduit::Node& matsetNode = m_bpNodeExt->fetch("matsets")[m_bpMatset];
  conduit::Node& vfValues = matsetNode["volume_fractions"][materialName];
  vfValues.set(conduit::DataType::float64(m_cellCount));
  double* vfPtr = vfValues.as_double_ptr();
  axom::copy(vfPtr, volumes.data(), m_cellCount * sizeof(double));
  if(!isFraction)
  {
    const double* cellVols = getCellVolumes().data();
    switch(m_runtimePolicy)
    {
    case RuntimePolicy::seq:
      elementwiseDivideImpl<axom::SEQ_EXEC>(vfPtr, cellVols, vfPtr, m_cellCount);
      break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
    case RuntimePolicy::omp:
      elementwiseDivideImpl<axom::OMP_EXEC>(vfPtr, cellVols, vfPtr, m_cellCount);
      break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
    case RuntimePolicy::cuda:
      elementwiseDivideImpl<axom::CUDA_EXEC<256>>(vfPtr,
                                                  cellVols,
                                                  vfPtr,
                                                  m_cellCount);
      break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
    case RuntimePolicy::hip:
      elementwiseDivideImpl<axom::HIP_EXEC<256>>(vfPtr,
                                                 cellVols,
                                                 vfPtr,
                                                 m_cellCount);
      break;
#endif
    default:
      SLIC_ERROR("ShapeeMesh internal error: Unhandled execution policy.");
    }
  }
}

void ShapeeMesh::setFreeVolumeFractions(const std::string& freeName)
{
  conduit::Node& vfsNode =
    m_bpNodeExt->fetch("matsets")[m_bpMatset]["volume_fractions"];

  SLIC_ERROR_IF(vfsNode.has_child(freeName),
                "Matset '" + m_bpMatset + "' already has a material named '" +
                  freeName + "'");

  conduit::Node& newVfNode = vfsNode[freeName];
  newVfNode.set(conduit::DataType::float64(m_cellCount));
  axom::ArrayView<double> newVfView(newVfNode.as_double_ptr(), {m_cellCount});

  fillNImpl(newVfView, 0.0);

  for(auto& vfNode : vfsNode.children())
  {
    if(vfNode.name() == newVfNode.name()) continue;
    axom::ArrayView<double> vfView(vfNode.as_double_ptr(), {m_cellCount});
    elementwiseAddImpl(newVfView, vfView, newVfView);
  }

  elementwiseComplementImpl(newVfView, 1.0, newVfView);
}

template <typename T>
void ShapeeMesh::fillNImpl(axom::ArrayView<T> a, const T& val) const
{
  auto kern = AXOM_LAMBDA(axom::IndexType i) { a[i] = val; };

  // Zero the new data for use as VF accumulation space.
  switch(m_runtimePolicy)
  {
  case RuntimePolicy::seq:
    axom::for_all<axom::SEQ_EXEC>(a.size(), kern);
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case RuntimePolicy::omp:
    axom::for_all<axom::OMP_EXEC>(a.size(), kern);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case RuntimePolicy::cuda:
    axom::for_all<axom::CUDA_EXEC<256>>(a.size(), kern);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case RuntimePolicy::hip:
    axom::for_all<axom::HIP_EXEC<256>>(a.size(), kern);
    break;
#endif
  default:
    SLIC_ERROR("ShapeeMesh internal error: Unhandled execution policy.");
  }
}

template <typename T>
void ShapeeMesh::elementwiseAddImpl(const axom::ArrayView<T> a,
                                    const axom::ArrayView<T> b,
                                    axom::ArrayView<T> result) const
{
  auto kern = AXOM_LAMBDA(axom::IndexType i) { result[i] = a[i] + b[i]; };

  switch(m_runtimePolicy)
  {
  case RuntimePolicy::seq:
    axom::for_all<axom::SEQ_EXEC>(result.size(), kern);
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case RuntimePolicy::omp:
    axom::for_all<axom::OMP_EXEC>(result.size(), kern);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case RuntimePolicy::cuda:
    axom::for_all<axom::CUDA_EXEC<256>>(result.size(), kern);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case RuntimePolicy::hip:
    axom::for_all<axom::HIP_EXEC<256>>(result.size(), kern);
    break;
#endif
  default:
    SLIC_ERROR("ShapeeMesh internal error: Unhandled execution policy.");
  }
}

template <typename T>
void ShapeeMesh::elementwiseComplementImpl(const axom::ArrayView<T> a,
                                           const T& val,
                                           axom::ArrayView<T> results) const
{
  auto kern = AXOM_LAMBDA(axom::IndexType i)
  {
    results[i] = val >= a[i] ? val - a[i] : 0.0;
  };

  switch(m_runtimePolicy)
  {
  case RuntimePolicy::seq:
    axom::for_all<axom::SEQ_EXEC>(a.size(), kern);
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case RuntimePolicy::omp:
    axom::for_all<axom::OMP_EXEC>(a.size(), kern);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case RuntimePolicy::cuda:
    axom::for_all<axom::CUDA_EXEC<256>>(a.size(), kern);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case RuntimePolicy::hip:
    axom::for_all<axom::HIP_EXEC<256>>(a.size(), kern);
    break;
#endif
  default:
    SLIC_ERROR("ShapeeMesh internal error: Unhandled execution policy.");
  }
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

template <typename ExecSpace, typename T>
void ShapeeMesh::elementwiseDivideImpl(const T* numerator,
                                       const T* denominator,
                                       T* quotient,
                                       axom::IndexType n)
{
  axom::for_all<ExecSpace>(
    n,
    AXOM_LAMBDA(axom::IndexType i) {
      quotient[i] = numerator[i] / denominator[i];
    });
}

}  // end namespace quest
}  // end namespace axom
