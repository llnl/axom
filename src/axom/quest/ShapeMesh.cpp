// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/primal.hpp"
#include "axom/quest/ShapeMesh.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/fmt.hpp"

#include "conduit_blueprint.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{
#if defined(AXOM_USE_64BIT_INDEXTYPE) && !defined(AXOM_NO_INT64_T)
static constexpr conduit::DataType::TypeID conduitDataIdOfAxomIndexType = conduit::DataType::INT64_ID;
#else
static constexpr conduit::DataType::TypeID conduitDataIdOfAxomIndexType = conduit::DataType::INT32_ID;
#endif

ShapeMesh::ShapeMesh(RuntimePolicy runtimePolicy,
                     int allocatorId,
                     conduit::Node& bpMesh,
                     const std::string& topoName,
                     const std::string& matsetName)
  : m_runtimePolicy(runtimePolicy)
  , m_allocId(allocatorId != axom::INVALID_ALLOCATOR_ID
                ? allocatorId
                : axom::policyToDefaultAllocatorID(runtimePolicy))
  , m_topoName(topoName.empty() && bpMesh["topologies"].number_of_children() > 0
                 ? bpMesh["topologies"].child(0).name()
                 : topoName)
  , m_matsetName(matsetName.empty() && bpMesh["matsets"].number_of_children() > 0
                   ? bpMesh.fetch("matsets").child(0).name()
                   : matsetName)
  , m_bpGrpExt(nullptr)
  , m_bpNodeExt(&bpMesh)
  , m_zeroThreshold(1e-10)
{
  SLIC_ERROR_IF(m_topoName.empty(),
                "Topology name was not provided, and no default topology was found.");

  const int hostAllocId = axom::execution_space<axom::SEQ_EXEC>::allocatorID();

  // We currently support only unstructured topo.
  const auto& typeNode =
    m_bpNodeExt->fetch_existing("topologies").fetch_existing(m_topoName).fetch_existing("type");
  const std::string topoType = typeNode.as_string();
  SLIC_ERROR_IF(topoType != "unstructured",
                "ShapeMesh currently only works with unstructured mesh, not " + topoType + ".");

  const conduit::Node& topoNode =
    m_bpNodeExt->fetch_existing("topologies").fetch_existing(m_topoName);
  const std::string coordsetName = topoNode.fetch_existing("coordset").as_string();
  const conduit::Node& coordsetNode =
    m_bpNodeExt->fetch_existing("coordsets").fetch_existing(coordsetName);

  if(!m_matsetName.empty())
  {
    conduit::Node& matsetNode = m_bpNodeExt->fetch("matsets").fetch(m_matsetName);

    // If matsetName was given, but topology data isn't set up yet, set it up.
    if(!matsetNode.has_child("topology"))
    {
      matsetNode.set_allocator(sidre::ConduitMemory::axomAllocIdToConduit(hostAllocId));
      matsetNode.fetch("topology").set_string(m_topoName);
    }

    SLIC_ERROR_IF(matsetNode["topology"].as_string() != m_topoName,
                  "matset's topology doesn't match the specified topology");
  }

  // Input checks.
  SLIC_ERROR_IF(topoNode["type"].as_string() != "unstructured",
                "topology type must be 'unstructured'");
  SLIC_ERROR_IF(topoNode["elements/shape"].as_string() != "hex", "element shape must be 'hex'");

  m_dim = conduit::blueprint::mesh::topology::dims(topoNode);
  SLIC_ASSERT(m_dim == 3);  // Will allow 2D when we support it.

  m_cellCount = conduit::blueprint::mesh::topology::length(topoNode);

  m_vertexCount = conduit::blueprint::mesh::coordset::length(coordsetNode);

  const conduit::Node& coordsValues = coordsetNode.fetch_existing("values");
  const bool isInterleaved = conduit::blueprint::mcarray::is_interleaved(coordsValues);
  const int stride = isInterleaved ? m_dim : 1;
  const char* dirNames[] = {"x", "y", "z"};
  for(int d = 0; d < m_dim; ++d)
  {
    m_vertCoordsViews3D[d] = axom::ArrayView<const double>(coordsValues[dirNames[d]].as_double_ptr(),
                                                           {m_vertexCount},
                                                           stride);
  }
}

#ifdef AXOM_USE_SIDRE
ShapeMesh::ShapeMesh(RuntimePolicy runtimePolicy,
                     int allocatorId,
                     sidre::Group* bpMesh,
                     const std::string& topoName,
                     const std::string& matsetName)
  : m_runtimePolicy(runtimePolicy)
  , m_allocId(allocatorId != axom::INVALID_ALLOCATOR_ID
                ? allocatorId
                : axom::policyToDefaultAllocatorID(runtimePolicy))
  , m_topoName(topoName.empty() && bpMesh->hasGroup("topologies") &&
                   bpMesh->getGroup("topologies")->getNumGroups() > 0
                 ? bpMesh->getGroup("topologies")->getGroup(0)->getName()
                 : topoName)
  , m_matsetName(matsetName.empty() && bpMesh->hasGroup("matsets") &&
                     bpMesh->getGroup("matsets")->getNumGroups() > 0
                   ? bpMesh->getGroup("matsets")->getGroup(0)->getName()
                   : matsetName)
  , m_bpGrpExt(bpMesh)
  , m_bpNodeInt()
  , m_zeroThreshold(1e-10)
{
  SLIC_ASSERT(m_topoName != sidre::InvalidName);
  SLIC_ERROR_IF(m_topoName.empty(),
                "Topology name was not provided, and no default topology was found.");

  m_bpGrpExt->createNativeLayout(m_bpNodeInt);

  // We want unstructured topo but can accomodate structured.
  const std::string topoType =
    m_bpNodeInt.fetch_existing("topologies").fetch_existing(m_topoName).fetch_existing("type").as_string();
  SLIC_ERROR_IF(topoType != "unstructured",
                "ShapeMesh currently only works with unstructured mesh, not " + topoType + ".");

  const conduit::Node& topoNode = m_bpNodeInt.fetch_existing("topologies").fetch_existing(m_topoName);
  const std::string coordsetName = topoNode.fetch_existing("coordset").as_string();
  const conduit::Node& coordsetNode =
    m_bpNodeInt.fetch_existing("coordsets").fetch_existing(coordsetName);

  if(!m_matsetName.empty())
  {
    conduit::Node& matsetNode = m_bpNodeInt.fetch("matsets").fetch(m_matsetName);

    // If matsetName was given, but not data isn't set up yet, set it up.
    if(!matsetNode.has_child("topology"))
    {
      matsetNode.fetch("topology").set_string(m_topoName);
    }

    SLIC_ERROR_IF(matsetNode["topology"].as_string() != m_topoName,
                  "matset's topology doesn't match the specified topology");
  }

  // Input checks.
  SLIC_ERROR_IF(topoNode["type"].as_string() != "unstructured",
                "topology type must be 'unstructured'");
  SLIC_ERROR_IF(topoNode["elements/shape"].as_string() != "hex", "element shape must be 'hex'");

  m_dim = conduit::blueprint::mesh::topology::dims(topoNode);
  SLIC_ASSERT(m_dim == 3);  // Will allow 2D when we support it.

  m_cellCount = conduit::blueprint::mesh::topology::length(topoNode);

  m_vertexCount = conduit::blueprint::mesh::coordset::length(coordsetNode);

  const conduit::Node& coordsValues = coordsetNode.fetch_existing("values");
  const bool isInterleaved = conduit::blueprint::mcarray::is_interleaved(coordsValues);
  const int stride = isInterleaved ? m_dim : 1;
  const char* dirNames[] = {"x", "y", "z"};
  for(int d = 0; d < m_dim; ++d)
  {
    m_vertCoordsViews3D[d] = axom::ArrayView<const double>(coordsValues[dirNames[d]].as_double_ptr(),
                                                           {m_vertexCount},
                                                           stride);
  }
}
#endif

void ShapeMesh::precomputeMeshData()
{
  AXOM_ANNOTATE_SCOPE("ShapeMesh::precomputeMeshData");
  getCellsAsHexes();
  getCellsAsTets();
  getCellVolumes();
  getTetVolumes();
  getCellBoundingBoxes();
  getCellLengths();
  getCellNodeConnectivity();
  getVertexPoints();
}

axom::ArrayView<const ShapeMesh::TetrahedronType> ShapeMesh::getCellsAsTets()
{
  if(m_cellsAsTets.size() != m_cellCount * NUM_TETS_PER_HEX)
  {
    computeCellsAsTets();
  }
  return m_cellsAsTets;
}

axom::ArrayView<const ShapeMesh::HexahedronType> ShapeMesh::getCellsAsHexes()
{
  if(m_cellsAsHexes.size() != m_cellCount)
  {
    computeCellsAsHexes();
  }
  return m_cellsAsHexes;
}

axom::ArrayView<const double> ShapeMesh::getCellVolumes()
{
  if(m_hexVolumes.size() != m_cellCount)
  {
    computeHexVolumes();
  }
  return m_hexVolumes.view();
}

axom::ArrayView<const double> ShapeMesh::getTetVolumes()
{
  if(m_tetVolumes.size() != m_cellCount * NUM_TETS_PER_HEX)
  {
    computeTetVolumes();
  }
  return m_tetVolumes.view();
}

axom::ArrayView<const ShapeMesh::BoundingBox3DType> ShapeMesh::getCellBoundingBoxes()
{
  if(m_hexBbs.size() != m_cellCount)
  {
    computeHexBbs();
  }
  return m_hexBbs.view();
}

axom::ArrayView<const double> ShapeMesh::getCellLengths()
{
  if(m_cellLengths.size() != m_cellCount)
  {
    computeCellLengths();
  }
  return m_cellLengths.view();
}

axom::ArrayView<const ShapeMesh::Point3DType> ShapeMesh::getVertexPoints()
{
  if(m_vertPoints3D.size() != m_vertexCount)
  {
    computeVertPoints();
  }
  return m_vertPoints3D.view();
}

axom::ArrayView<const axom::IndexType, 2> ShapeMesh::getCellNodeConnectivity()
{
  if(m_connectivity.size() != m_cellCount)
  {
    computeConnectivity();
  }
  return m_connectivity;
}

bool ShapeMesh::isValidForShaping(std::string& whyNot) const
{
  bool rval = true;

  // We run the check on the Conduit form of the mesh.
  const conduit::Node& bpMesh = m_bpNodeExt ? *m_bpNodeExt : m_bpNodeInt;

  /*
   * Check Blueprint-validity.
   * Conduit's verify should work even if m_bpNodeExt has array data on
   * devices.  The verification doesn't dereference array data.
   * If this changes in the future, this code must adapt.
   */
  conduit::Node info;
  rval = conduit::blueprint::mesh::verify(bpMesh, info);

  if(rval)
  {
    std::string topoType = bpMesh.fetch("topologies")[m_topoName]["type"].as_string();
    rval = topoType == "unstructured";
    if(!rval)
    {
      info[0].set_string("Topology is not unstructured.");
    }
  }

  if(rval)
  {
    std::string elemShape = bpMesh.fetch("topologies")[m_topoName]["elements/shape"].as_string();
    rval = elemShape == "hex";
    if(!rval)
    {
      info[0].set_string("Topology elements are not hex.");
    }
  }

  whyNot = info.to_summary_string();

  return rval;
}

void ShapeMesh::setMatsetFromVolume(const std::string& materialName,
                                    const axom::ArrayView<double>& volumes,
                                    bool isFraction)
{
  SLIC_ERROR_IF(m_matsetName.empty(),
                "Cannot use material set in ShapeMesh: Matset name was not provided, and no "
                "default matset was found.");
  SLIC_ERROR_IF(materialName.empty(), "Cannot have an empty materialName.");

  double* vfPtr = nullptr;

  auto dataType = conduit::DataType::float64(m_cellCount);

  if(m_bpNodeExt != nullptr)
  {
    conduit::Node& matsetNode = m_bpNodeExt->fetch("matsets")[m_matsetName];
    if(matsetNode.has_child("topology"))
    {
      SLIC_ASSERT(matsetNode.fetch_existing("topology").as_string() == m_topoName);
    }
    else
    {
      matsetNode["topology"].set_string(m_topoName);
    }
    const std::string vfValuesPath = "matsets/" + m_matsetName + "/volume_fractions/" + materialName;
    conduit::Node& vfValues = getMeshConduitPath(*m_bpNodeExt, vfValuesPath, dataType);
    vfPtr = vfValues.as_double_ptr();
  }
#if defined(AXOM_USE_SIDRE)
  if(m_bpGrpExt != nullptr)
  {
    std::string viewPath = "matsets/" + m_matsetName + "/volume_fractions/" + materialName;
    sidre::View* vfValues = m_bpGrpExt->hasView(viewPath)
      ? m_bpGrpExt->getView(viewPath)
      : m_bpGrpExt->createView("matsets/" + m_matsetName + "/volume_fractions/" + materialName);
    if(!vfValues->isAllocated())
    {
      vfValues->allocate(dataType, m_allocId);
    }
    else
    {
      SLIC_ASSERT(vfValues->getSchema().dtype().id() == dataType.id());
      SLIC_ASSERT(vfValues->getSchema().dtype().number_of_elements() == dataType.number_of_elements());
    }
    vfPtr = (double*)(vfValues->getVoidPtr());
  }
#endif

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
      elementwiseDivideImpl<axom::CUDA_EXEC<256>>(vfPtr, cellVols, vfPtr, m_cellCount);
      break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
    case RuntimePolicy::hip:
      elementwiseDivideImpl<axom::HIP_EXEC<256>>(vfPtr, cellVols, vfPtr, m_cellCount);
      break;
#endif
    default:
      SLIC_ERROR("ShapeMesh internal error: Unhandled execution policy.");
    }
  }
}

void ShapeMesh::setFreeVolumeFractions(const std::string& freeName)
{
  SLIC_ERROR_IF(m_matsetName.empty(),
                "Cannot use material set in ShapeMesh: Matset name was not provided, and no "
                "default matset was found.");
  SLIC_ERROR_IF(freeName.empty(), "Cannot have an empty material name.");

  auto dataType = conduit::DataType::float64(m_cellCount);

  double* newVfPtr = nullptr;
  axom::ArrayView<double> newVfView(newVfPtr, {m_cellCount});

  if(m_bpNodeExt != nullptr)
  {
    conduit::Node& matsetNode = m_bpNodeExt->fetch("matsets")[m_matsetName];
    if(matsetNode.has_child("topology"))
    {
      SLIC_ASSERT(matsetNode.fetch_existing("topology").as_string() == m_topoName);
    }
    else
    {
      matsetNode["topology"].set_string(m_topoName);
    }

    conduit::Node& vfsNode = matsetNode["volume_fractions"];

    conduit::Node& freeVfNode = getMeshConduitPath(vfsNode, freeName, dataType);
    newVfPtr = freeVfNode.as_double_ptr();
    newVfView = axom::ArrayView<double>(newVfPtr, {m_cellCount});
    fillNImpl(newVfView, 0.0);

    for(auto& vfNode : vfsNode.children())
    {
      if(vfNode.name() == freeVfNode.name()) continue;
      axom::ArrayView<double> vfView(vfNode.as_double_ptr(), {m_cellCount});
      elementwiseAddImpl(newVfView, vfView, newVfView);
    }
  }
#if defined(AXOM_USE_SIDRE)
  if(m_bpGrpExt != nullptr)
  {
    sidre::Group* matsetGrp = m_bpGrpExt->createGroup("matsets/" + m_matsetName, false, true);
    if(matsetGrp->hasView("topology"))
    {
      SLIC_ERROR_IF(
        matsetGrp->getView("topology")->getString() != m_topoName,
        "Material set '" + m_matsetName + "' doesn't have expected topology '" + m_topoName + "'");
    }
    else
    {
      matsetGrp->createView("topology")->setString(m_topoName);
    }

    sidre::Group* vfsGrp = matsetGrp->createGroup("volume_fractions", false, true);

    sidre::View* newVfVu = vfsGrp->createView(freeName);
    newVfVu->allocate(dataType, m_allocId);
    newVfPtr = (double*)(newVfVu->getVoidPtr());
    newVfView = axom::ArrayView<double>(newVfPtr, {m_cellCount});
    fillNImpl(newVfView, 0.0);

    for(auto& vfVu : vfsGrp->views())
    {
      if(vfVu.getName() == freeName) continue;
      axom::ArrayView<double> vfView((double*)(vfVu.getVoidPtr()), {m_cellCount});
      elementwiseAddImpl(newVfView, vfView, newVfView);
    }
  }
#endif

  elementwiseComplementImpl(newVfView, 1.0, newVfView);
}

template <typename T>
void ShapeMesh::fillNImpl(axom::ArrayView<T> a, const T& val) const
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
    SLIC_ERROR("ShapeMesh internal error: Unhandled execution policy.");
  }
}

template <typename T>
void ShapeMesh::elementwiseAddImpl(const axom::ArrayView<T> a,
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
    SLIC_ERROR("ShapeMesh internal error: Unhandled execution policy.");
  }
}

template <typename T>
void ShapeMesh::elementwiseComplementImpl(const axom::ArrayView<T> a,
                                          const T& val,
                                          axom::ArrayView<T> results) const
{
  auto kern = AXOM_LAMBDA(axom::IndexType i) { results[i] = val >= a[i] ? val - a[i] : 0.0; };

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
    SLIC_ERROR("ShapeMesh internal error: Unhandled execution policy.");
  }
}

void ShapeMesh::computeCellsAsHexes()
{
  AXOM_ANNOTATE_SCOPE("ShapeMesh::computeCellsAsHexes");
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

void ShapeMesh::computeCellsAsTets()
{
  AXOM_ANNOTATE_SCOPE("ShapeMesh::computeCellsAsTets");
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

void ShapeMesh::computeHexVolumes()
{
  AXOM_ANNOTATE_SCOPE("ShapeMesh::computeHexVolumes");
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

void ShapeMesh::computeTetVolumes()
{
  AXOM_ANNOTATE_SCOPE("ShapeMesh::computeTetVolumes");
  switch(m_runtimePolicy)
  {
  case RuntimePolicy::seq:
    computeTetVolumesImpl<axom::SEQ_EXEC>();
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case RuntimePolicy::omp:
    computeTetVolumesImpl<axom::OMP_EXEC>();
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case RuntimePolicy::cuda:
    computeTetVolumesImpl<axom::CUDA_EXEC<256>>();
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case RuntimePolicy::hip:
    computeTetVolumesImpl<axom::HIP_EXEC<256>>();
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
}

void ShapeMesh::computeHexBbs()
{
  AXOM_ANNOTATE_SCOPE("ShapeMesh::computeHexBoundingBoxes");
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

void ShapeMesh::computeVertPoints()
{
  AXOM_ANNOTATE_SCOPE("ShapeMesh::computeVertPoints(");
  switch(m_runtimePolicy)
  {
  case RuntimePolicy::seq:
    computeVertPointsImpl<axom::SEQ_EXEC>();
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case RuntimePolicy::omp:
    computeVertPointsImpl<axom::OMP_EXEC>();
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case RuntimePolicy::cuda:
    computeVertPointsImpl<axom::CUDA_EXEC<256>>();
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case RuntimePolicy::hip:
    computeVertPointsImpl<axom::HIP_EXEC<256>>();
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
}

void ShapeMesh::computeCellLengths()
{
  AXOM_ANNOTATE_SCOPE("ShapeMesh::computeCellLengths");
  switch(m_runtimePolicy)
  {
  case axom::runtime_policy::Policy::seq:
    computeCellLengthsImpl<axom::SEQ_EXEC>();
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    computeCellLengthsImpl<axom::OMP_EXEC>();
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case axom::runtime_policy::Policy::cuda:
    computeCellLengthsImpl<axom::CUDA_EXEC<256>>();
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    computeCellLengthsImpl<axom::HIP_EXEC<256>>();
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
}

void ShapeMesh::computeConnectivity()
{
  SLIC_ASSERT(m_dim == 3);  // 2D support not done yet.

  constexpr int NUM_VERTS_PER_HEX = 8;

  conduit::Node& bpMesh = m_bpNodeExt ? *m_bpNodeExt : m_bpNodeInt;

  conduit::Node& topoNode = bpMesh.fetch_existing("topologies").fetch_existing(m_topoName);
  auto& connNode = topoNode.fetch_existing("elements/connectivity");
  SLIC_ERROR_IF(connNode.dtype().id() != conduitDataIdOfAxomIndexType,
                "IntersectionShaper error: connectivity data type must be "
                "axom::IndexType.");
  const auto* connPtr = static_cast<const axom::IndexType*>(connNode.data_ptr());
  m_connectivity =
    axom::ArrayView<const axom::IndexType, 2> {connPtr, m_cellCount, NUM_VERTS_PER_HEX};
}

template <typename ExecSpace>
void ShapeMesh::computeCellsAsHexesImpl()
{
  constexpr static int NUM_VERTS_PER_HEX = 8;
  constexpr static int NDIM = 3;

  SLIC_ASSERT(m_dim == NDIM);  // or we shouldn't be here.

  const auto& vertexCoords = getVertexCoords3D();
  const auto& vX = vertexCoords[0];
  const auto& vY = vertexCoords[1];
  const auto& vZ = vertexCoords[2];

  axom::ArrayView<const IndexType, 2> connView = getCellNodeConnectivity();

  m_cellsAsHexes =
    axom::Array<HexahedronType>(ArrayOptions::Uninitialized(), m_cellCount, m_cellCount, m_allocId);
  axom::ArrayView<HexahedronType> cellsAsHexesView = m_cellsAsHexes.view();
  SLIC_ASSERT(cellsAsHexesView.data() == m_cellsAsHexes.data());

  const auto zeroThreshold = m_zeroThreshold;
  axom::for_all<ExecSpace>(
    m_cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
      auto& hex = cellsAsHexesView[cellId];

      for(int vi = 0; vi < NUM_VERTS_PER_HEX; ++vi)
      {
        axom::IndexType vertIndex = connView(cellId, vi);
        primal::Point3D vCoords(
          axom::NumericArray<double, NDIM> {vX[vertIndex], vY[vertIndex], vZ[vertIndex]});

        // Snap coordinates to zero.
        for(int d = 0; d < NDIM; ++d)
        {
          if(axom::utilities::isNearlyEqual(vCoords[d], 0.0, zeroThreshold))
          {
            vCoords[d] = 0;
          }
        }

        hex[vi] = vCoords;
      }
    });  // end of loop to initialize hexahedral elements and bounding boxes
}

template <typename ExecSpace>
void ShapeMesh::computeCellsAsTetsImpl()
{
  SLIC_ASSERT(m_dim == 3);  // or we shouldn't be here.

  m_cellsAsTets = axom::Array<TetrahedronType>(ArrayOptions::Uninitialized(),
                                               NUM_TETS_PER_HEX * m_cellCount,
                                               NUM_TETS_PER_HEX * m_cellCount,
                                               m_allocId);
  auto cellsAsTetsView = m_cellsAsTets.view();

  auto cellsAsHexesView = getCellsAsHexes();

  axom::for_all<ExecSpace>(
    m_cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
      const auto& hex = cellsAsHexesView[cellId];
      auto* firstTetPtr = &cellsAsTetsView[cellId * NUM_TETS_PER_HEX];
      hexToTets(hex, firstTetPtr);
    });
}

template <typename ExecSpace>
void ShapeMesh::computeHexVolumesImpl()
{
  m_hexVolumes =
    axom::Array<double>(ArrayOptions::Uninitialized(), m_cellCount, m_cellCount, m_allocId);

  auto cellsAsHexes = getCellsAsHexes();

  auto hexVolumesView = m_hexVolumes.view();
  axom::for_all<ExecSpace>(
    m_cellCount,
    AXOM_LAMBDA(axom::IndexType i) { hexVolumesView[i] = cellsAsHexes[i].volume(); });
}

template <typename ExecSpace>
void ShapeMesh::computeTetVolumesImpl()
{
  axom::IndexType tetCount = m_cellCount * NUM_TETS_PER_HEX;
  m_tetVolumes = axom::Array<double>(ArrayOptions::Uninitialized(), tetCount, tetCount, m_allocId);

  auto cellsAsTets = getCellsAsTets();

  auto tetVolumesView = m_tetVolumes.view();
  axom::for_all<ExecSpace>(
    tetCount,
    AXOM_LAMBDA(axom::IndexType i) { tetVolumesView[i] = cellsAsTets[i].volume(); });
}

template <typename ExecSpace>
void ShapeMesh::computeHexBbsImpl()
{
  m_hexBbs =
    axom::Array<BoundingBox3DType>(ArrayOptions::Uninitialized(), m_cellCount, m_cellCount, m_allocId);

  auto cellsAsHexes = getCellsAsHexes();

  auto hexBbsView = m_hexBbs.view();
  axom::for_all<ExecSpace>(
    m_cellCount,
    AXOM_LAMBDA(axom::IndexType i) {
      hexBbsView[i] = primal::compute_bounding_box<double, 3>(cellsAsHexes[i]);
    });
}

template <typename ExecSpace>
void ShapeMesh::computeCellLengthsImpl()
{
  m_cellLengths =
    axom::Array<double>(ArrayOptions::Uninitialized(), m_cellCount, m_cellCount, m_allocId);

  auto cellBbs = getCellBoundingBoxes();

  auto lengthsView = m_cellLengths.view();
  axom::for_all<ExecSpace>(
    m_cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) { lengthsView[cellId] = cellBbs[cellId].range().norm(); });
}

template <typename ExecSpace>
void ShapeMesh::computeVertPointsImpl()
{
  m_vertPoints3D = axom::Array<Point3DType>(m_vertexCount, m_vertexCount, m_allocId);

  auto& vertCoords = getVertexCoords3D();
  const auto& vX = vertCoords[0];
  const auto& vY = vertCoords[1];
  const auto& vZ = vertCoords[2];

  auto vertPointsView = m_vertPoints3D.view();
  axom::for_all<ExecSpace>(
    m_vertexCount,
    AXOM_LAMBDA(axom::IndexType vi) { vertPointsView[vi] = Point3DType {vX[vi], vY[vi], vZ[vi]}; });
}

template <typename ExecSpace, typename T>
void ShapeMesh::elementwiseDivideImpl(const T* numerator,
                                      const T* denominator,
                                      T* quotient,
                                      axom::IndexType n)
{
  axom::for_all<ExecSpace>(
    n,
    AXOM_LAMBDA(axom::IndexType i) { quotient[i] = numerator[i] / denominator[i]; });
}

conduit::Node& ShapeMesh::getMeshConduitPath(conduit::Node& node,
                                             const std::string& path,
                                             const conduit::DataType& dtype)
{
  conduit::Node* rval = nullptr;

  if(node.has_path(path))
  {
    rval = &node.fetch_existing(path);
    SLIC_ERROR_IF(
      rval->dtype().id() != dtype.id() ||
        rval->dtype().number_of_elements() != dtype.number_of_elements() ||
        rval->dtype().offset() != dtype.offset() || rval->dtype().stride() != dtype.stride() ||
        rval->dtype().offset() != dtype.offset() || rval->dtype().endianness() != dtype.endianness(),
      "Blueprint mesh doesn't have correct type for path '" + path + "'");
  }
  else
  {
    node.set_allocator(sidre::ConduitMemory::axomAllocIdToConduit(m_allocId));
    rval = &node.fetch(path);
    rval->set_allocator(sidre::ConduitMemory::axomAllocIdToConduit(m_allocId));
    rval->set(dtype);
  }

  // Verify that data is in memory space usable by m_runtimePolicy.
  void* dataPtr = rval->data_ptr();
  int allocId = axom::getAllocatorIDFromPointer(dataPtr);

  bool memoryOkForPolicy =
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
    m_runtimePolicy == axom::runtime_policy::Policy::omp
    ? axom::execution_space<axom::OMP_EXEC>::usesAllocId(allocId)
    :
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
    m_runtimePolicy == axom::runtime_policy::Policy::cuda
    ? axom::execution_space<axom::CUDA_EXEC<256>>::usesAllocId(allocId)
    :
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
    m_runtimePolicy == axom::runtime_policy::Policy::hip
    ? axom::execution_space<axom::HIP_EXEC<256>>::usesAllocId(allocId)
    :
#endif
    axom::execution_space<axom::SEQ_EXEC>::usesAllocId(allocId);

  SLIC_WARNING_IF(!memoryOkForPolicy,
                  "Blueprint mesh data at " + axom::fmt::format("{}", dataPtr) + " from path '" +
                    path + "' is not accessible by execution policy " +
                    axom::runtime_policy::policyToName(m_runtimePolicy));

  return *rval;
}

}  // end namespace experimental
}  // end namespace quest
}  // end namespace axom
