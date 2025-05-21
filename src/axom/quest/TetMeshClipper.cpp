// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifdef AXOM_USE_CONDUIT
  #include "conduit_blueprint.hpp"
#endif

#include "axom/mint/mesh/Mesh.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/quest/Discretize.hpp"
#include "axom/quest/TetMeshClipper.hpp"
#include "axom/fmt.hpp"

namespace axom
{
namespace quest
{

TetMeshClipper::TetMeshClipper(const klee::Geometry& kGeom, const std::string& name)
  : GeometryClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("TetMesh") : name)
  , m_topoName(kGeom.getBlueprintTopology())
  , m_cellCount(0)
  , m_transformer(m_transMat)
{
  SLIC_ASSERT(!m_topoName.empty());

  extractClipperInfo();

  transformCoordset();
}

bool TetMeshClipper::labelInOut(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  AXOM_UNUSED_VAR(shapeeMesh);
  AXOM_UNUSED_VAR(labels);
  AXOM_ANNOTATE_BEGIN("TetMeshClipper::labelInOut");
  switch(shapeeMesh.getRuntimePolicy())
  {
  case axom::runtime_policy::Policy::seq:
    labelInOutImpl<axom::SEQ_EXEC>(shapeeMesh, labels);
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    labelInOutImpl<axom::OMP_EXEC>(shapeeMesh, labels);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case axom::runtime_policy::Policy::cuda:
    labelInOutImpl<axom::CUDA_EXEC<256>>(shapeeMesh, labels);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    labelInOutImpl<axom::HIP_EXEC<256>>(shapeeMesh, labels);
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
  AXOM_ANNOTATE_END("TetMeshClipper::labelInOut");
  return true;
}

template <typename ExecSpace>
void TetMeshClipper::labelInOutImpl(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  SLIC_ERROR_IF(shapeeMesh.dimension() != 3, "TetMeshClipper requires a 3D mesh.");

  constexpr int NUM_VERTS_PER_CELL = 8;

  int allocId = shapeeMesh.getAllocatorId();
  auto cellCount = shapeeMesh.getCellCount();
  auto vertCount = shapeeMesh.getVertexCount();

  const auto& vertCoords = shapeeMesh.getVertexCoords3D();
  const auto& vX = vertCoords[0];
  const auto& vY = vertCoords[1];
  const auto& vZ = vertCoords[2];

  /*
    Compute whether mesh vertices are inside the mesh's bounding box.
    This is more overly conservative compared to checking against the
    mesh itself, but much faster.
  */
  axom::Array<bool> vertIsInside {ArrayOptions::Uninitialized(), vertCount, vertCount, allocId};
  auto vertIsInsideView = vertIsInside.view();
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vX.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vY.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vZ.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vertIsInsideView.getAllocatorID()));

  auto bb = m_bb;
  axom::for_all<ExecSpace>(
    vertCount,
    AXOM_LAMBDA(axom::IndexType vertId) {
      primal::Point3D vert {vX[vertId], vY[vertId], vZ[vertId]};
      vertIsInsideView[vertId] = bb.contains(vert);
    });

  if(labels.size() < cellCount || labels.getAllocatorID() != shapeeMesh.getAllocatorId())
  {
    labels = axom::Array<LabelType>(ArrayOptions::Uninitialized(), cellCount, cellCount, allocId);
  }

  /*
    Label cell outside if no vertex is in the bounding box.
    Otherwise, label it on boundary, because we don't know.
  */
  axom::ArrayView<const axom::IndexType, 2> connView = shapeeMesh.getConnectivity();
  SLIC_ASSERT(connView.shape() ==
              (axom::StackArray<axom::IndexType, 2> {cellCount, NUM_VERTS_PER_CELL}));

  auto labelsView = labels.view();

  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
      auto cellVertIds = connView[cellId];
      bool hasIn = vertIsInsideView[cellVertIds[0]];
      for(int vi = 0; vi < NUM_VERTS_PER_CELL; ++vi)
      {
        int vertId = cellVertIds[vi];
        bool isIn = vertIsInsideView[vertId];
        hasIn |= isIn;
      }
      labelsView[cellId] = hasIn ? LABEL_ON : LABEL_OUT;
    });

  return;
}

bool TetMeshClipper::getGeometryAsTets(quest::ShapeeMesh& shapeeMesh, axom::Array<TetrahedronType>& tets)
{
  AXOM_ANNOTATE_BEGIN("TetMeshClipper::getGeometryAsTets");
  const int hostAllocId = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int allocId = shapeeMesh.getAllocatorId();

  if(tets.getAllocatorID() != allocId || tets.size() != m_cellCount)
  {
    tets = axom::Array<TetrahedronType>(m_cellCount, m_cellCount, allocId);
  }

  /*
    1. Initialize a mint Mesh intermediary from the blueprint mesh.
       mint::getMesh() utility for this saves some coding.
    2. Populate a tet array on the host (because mint data works only for host).
       Host array could be tets or a temporary array.
    3. If needed, copy temporary array to tets.
  */

  axom::Array<TetrahedronType> tmpTets(0, 0, hostAllocId);

  // Host loop writes directly to tets if possible, else a temporary array.
  axom::Array<TetrahedronType>& tetsOnHost =
    axom::execution_space<axom::SEQ_EXEC>::usesAllocId(tets.getAllocatorID()) ? tets : tmpTets;

  if(&tetsOnHost == &tmpTets)
  {
    tmpTets.resize(m_cellCount, TetrahedronType());
  }

  axom::sidre::DataStore ds;
  auto* bpMeshGrp = ds.getRoot()->createGroup("blueprintMesh");
  bpMeshGrp->importConduitTree(*m_bpMesh);

  const bool addExtraDataForMint = true;
  if(addExtraDataForMint)
  {
    /*
      Constructing a mint mesh from meshGrp fails unless we add some
      extra data.  Blueprint doesn't require this extra data.  (The mesh
      passes conduit's Blueprint verification.)  This should be fixed,
      or we should write better blueprint support utilities.
    */
    auto* topoGrp = bpMeshGrp->getGroup("topologies")->getGroup(m_topoName);
    auto* coordValuesGrp =
      bpMeshGrp->getGroup("coordsets")->getGroup(m_coordsetName)->getGroup("values");
    /*
      Make the coordinate arrays 2D to use mint::Mesh.
      For some reason, mint::Mesh requires the arrays to be
      2D, even though the second dimension is always 1.
    */
    axom::IndexType curShape[2];
    int curDim;
    curDim = coordValuesGrp->getView("x")->getShape(2, curShape);
    assert(curDim == 1);
    const axom::IndexType vertsShape[2] = {curShape[0], 1};
    coordValuesGrp->getView("x")->reshapeArray(2, vertsShape);
    coordValuesGrp->getView("y")->reshapeArray(2, vertsShape);
    coordValuesGrp->getView("z")->reshapeArray(2, vertsShape);

    // Make connectivity array 2D for the same reason.
    auto* elementsGrp = topoGrp->getGroup("elements");
    auto* connView = elementsGrp->getView("connectivity");
    curDim = connView->getShape(2, curShape);
    constexpr axom::IndexType NUM_VERTS_PER_TET = 4;
    SLIC_ASSERT(curDim == 1 || curDim == 2);
    if(curDim == 1)
    {
      SLIC_ASSERT(curShape[0] % NUM_VERTS_PER_TET == 0);
      axom::IndexType connShape[2] = {curShape[0] / NUM_VERTS_PER_TET, NUM_VERTS_PER_TET};
      connView->reshapeArray(2, connShape);
    }

    // mint::Mesh requires connectivity strides, even though Blueprint doesn't.
    if(!elementsGrp->hasView("stride"))
    {
      elementsGrp->createViewScalar("stride", NUM_VERTS_PER_TET, hostAllocId);
    }

    // mint::Mesh requires field group, even though Blueprint doesn't.
    if(!bpMeshGrp->hasGroup("fields"))
    {
      bpMeshGrp->createGroup("fields");
    }
  }
  std::shared_ptr<mint::UnstructuredMesh<axom::mint::Topology::SINGLE_SHAPE>> mintMesh {
    (mint::UnstructuredMesh<axom::mint::Topology::SINGLE_SHAPE>*)axom::mint::getMesh(bpMeshGrp,
                                                                                     m_topoName)};

  // Initialize tetrahedra
  IndexType nodeIds[4];
  Point3D pts[4];
  for(int i = 0; i < m_cellCount; i++)
  {
    mintMesh->getCellNodeIDs(i, nodeIds);

    mintMesh->getNode(nodeIds[0], pts[0].data());
    mintMesh->getNode(nodeIds[1], pts[1].data());
    mintMesh->getNode(nodeIds[2], pts[2].data());
    mintMesh->getNode(nodeIds[3], pts[3].data());

    tetsOnHost[i] = TetrahedronType({pts[0], pts[1], pts[2], pts[3]});
  }

  if(&tets != &tetsOnHost)
  {
    axom::copy(tets.data(), tetsOnHost.data(), sizeof(TetrahedronType) * tets.size());
  }

  AXOM_ANNOTATE_END("TetMeshClipper::getGeometryAsTets");
  return true;
}

void TetMeshClipper::extractClipperInfo()
{
  m_topoName = m_info.fetch_existing("topologyName").as_string();

  m_bpMesh = &m_info.fetch_existing("klee::Geometry:tetMesh");

  SLIC_ASSERT(
    m_bpMesh->fetch_existing("topologies").fetch_existing(m_topoName).fetch_existing("type").as_string() ==
    "unstructured");

  conduit::Node& topoNode = m_bpMesh->fetch_existing("topologies").fetch_existing(m_topoName);

  bool isMultiDomain = conduit::blueprint::mesh::is_multi_domain(*m_bpMesh);
  SLIC_ERROR_IF(isMultiDomain, "TetMeshClipper does not support multi-domain tet meshes yet.");

  SLIC_ASSERT(conduit::blueprint::mesh::topology::dims(topoNode) == 3);

  m_cellCount = conduit::blueprint::mesh::topology::length(topoNode);

  m_coordsetName = topoNode.fetch_existing("coordset").as_string();
}

void TetMeshClipper::transformCoordset()
{
  // Apply transformations
  auto& oldCoordset = m_bpMesh->fetch_existing("coordsets").fetch_existing(m_coordsetName);
  const std::string newCoordsetName = m_coordsetName + ".trans";
  conduit::Node& coordset = m_bpMesh->fetch("coordsets")[newCoordsetName];
  coordset.set_node(oldCoordset);
  auto transformer = m_transformer;
  conduit::index_t count = conduit::blueprint::mesh::coordset::length(coordset);
  axom::ArrayView<double> xV(coordset.fetch_existing("values/x").as_double_ptr(), count);
  axom::ArrayView<double> yV(coordset.fetch_existing("values/y").as_double_ptr(), count);
  axom::ArrayView<double> zV(coordset.fetch_existing("values/z").as_double_ptr(), count);
  axom::for_all<axom::SEQ_EXEC>(count, AXOM_LAMBDA(axom::IndexType i) {
      transformer.transform(xV[i], yV[i], zV[i]);
      m_bb.addPoint(Point3DType{xV[i], yV[i], zV[i]});
    });
  m_bpMesh->fetch_existing("topologies").fetch_existing(m_topoName).fetch_existing("coordset").set_string(newCoordsetName);
  m_coordsetName = newCoordsetName;
}

}  // end namespace quest
}  // end namespace axom
