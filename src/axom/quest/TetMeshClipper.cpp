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
#include "axom/spin/BVH.hpp"
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
  , m_tetCount(0)
  , m_transformer(m_transMat)
{
  SLIC_ASSERT(!m_topoName.empty());

  extractClipperInfo();

  transformCoordset();

  computeTets();
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

/*
  1. Compute whether vertices are in or out of the tet mesh.
  2. Determine whether cells are in, out or on the tet mesh
     boundary.
  Unlike the TetClipper, this doesn't check edge-tet intersections,
  so it has errors.  These errors should shrink with mesh resolution.
  If needed, we can implement edge-tet detection for TetMeshClipper.
*/
template <typename ExecSpace>
void TetMeshClipper::labelInOutImpl(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  SLIC_ERROR_IF(shapeeMesh.dimension() != 3, "TetMeshClipper requires a 3D mesh.");

  int allocId = shapeeMesh.getAllocatorID();
  auto vertCount = shapeeMesh.getVertexCount();

  /*
    Tets (as bounding boxes) in BVH.
  */
  axom::Array<BoundingBox3DType> tetsAsBbs(m_tets.size(), m_tets.size(), allocId);
  auto tetsAsBbsView = tetsAsBbs.view();
  auto tetsView = m_tets.view();
  axom::for_all<ExecSpace>(tetsView.size(),
                           AXOM_LAMBDA(axom::IndexType vi) {
                             const auto&tet = tetsView[vi];
                             tetsAsBbsView[vi] = axom::primal::compute_bounding_box(tet);
                           });

  spin::BVH<3, ExecSpace, double> bvh;
  bvh.initialize(tetsAsBbsView, tetsAsBbsView.size());

  /*
    Compute whether vertices are inside tet mesh.
    (Use BVH to narrow the search to nearby candidate tets.)
  */
  axom::Array<bool> vertIsInside {ArrayOptions::Uninitialized(), vertCount, vertCount, allocId};
  vertIsInside.fill(false);
  auto vertIsInsideView = vertIsInside.view();

  axom::ArrayView<const Point3DType> vertPointsView = shapeeMesh.getVertexPoints();

  axom::Array<IndexType> offsets(vertCount, vertCount, allocId);
  axom::Array<IndexType> counts(vertCount, vertCount, allocId);
  axom::Array<IndexType> candidates;
  bvh.findPoints(offsets, counts, candidates, vertCount, vertPointsView);

  axom::for_all<ExecSpace>(
    vertCount,
    AXOM_LAMBDA(axom::IndexType vertId) {
      auto candidateCount = counts[vertId];
      bool& isInside = vertIsInsideView[vertId];
      if (!isInside && candidateCount > 0)
      {
        auto candidateIds = &candidates[offsets[vertId]];
        auto& vertex = vertPointsView[vertId];
        for(int ci = 0; ci < candidateCount && !isInside; ++ci)
        {
          axom::IndexType tetId = candidateIds[ci];
          const auto& tet = tetsView[tetId];
          isInside |= tet.contains(vertex);
        }
      }
    });

  vertexInsideToCellLabel<ExecSpace>(shapeeMesh, vertIsInsideView, labels);
}

/*
  Label cell outside if no vertex is in the bounding box.
  Otherwise, label it on boundary, because we don't know.
*/
template <typename ExecSpace>
void TetMeshClipper::vertexInsideToCellLabel(
  quest::ShapeeMesh& shapeeMesh,
  axom::ArrayView<bool>& vertIsInside,
  axom::Array<LabelType>& labels)
{
  axom::ArrayView<const axom::IndexType, 2> connView = shapeeMesh.getConnectivity();
  SLIC_ASSERT(connView.shape() ==
              (axom::StackArray<axom::IndexType, 2> {shapeeMesh.getCellCount(), HexahedronType::NUM_HEX_VERTS}));

  if(labels.size() < shapeeMesh.getCellCount() || labels.getAllocatorID() != shapeeMesh.getAllocatorID())
  {
    labels = axom::Array<LabelType>(ArrayOptions::Uninitialized(), shapeeMesh.getCellCount(), shapeeMesh.getCellCount(), shapeeMesh.getAllocatorID());
  }
  auto labelsView = labels.view();

  axom::for_all<ExecSpace>(
    shapeeMesh.getCellCount(),
    AXOM_LAMBDA(axom::IndexType cellId) {
      auto cellVertIds = connView[cellId];
      bool hasIn = vertIsInside[cellVertIds[0]];
      bool hasOut = !hasIn;
      for(int vi = 0; vi < HexahedronType::NUM_HEX_VERTS; ++vi)
      {
        int vertId = cellVertIds[vi];
        bool isIn = vertIsInside[vertId];
        hasIn |= isIn;
        hasOut |= !isIn;
      }
      labelsView[cellId] = !hasOut ? LABEL_IN : !hasIn ? LABEL_OUT : LABEL_ON;
    });

  return;
}

bool TetMeshClipper::getGeometryAsTets(quest::ShapeeMesh& shapeeMesh, axom::Array<TetrahedronType>& tets)
{
  AXOM_ANNOTATE_SCOPE("TetMeshClipper::getGeometryAsTets");
  tets = axom::Array<TetrahedronType>(m_tets, shapeeMesh.getAllocatorID());
  return true;
}

// Compute m_tets.  Keep data on host.  We don't know what allocator the ShapeeMesh uses.
void TetMeshClipper::computeTets()
{
  AXOM_ANNOTATE_SCOPE("TetMeshClipper::computeTets");
  const int hostAllocId = axom::execution_space<axom::SEQ_EXEC>::allocatorID();

  m_tets = axom::Array<TetrahedronType>(m_tetCount, m_tetCount, hostAllocId);

  /*
    1. Initialize a mint Mesh intermediary from the blueprint mesh.
       mint::getMesh() utility for this saves some coding.
    2. Populate a tet array on the host (because mint data works only for host).
  */

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
  for(int i = 0; i < m_tetCount; i++)
  {
    mintMesh->getCellNodeIDs(i, nodeIds);

    mintMesh->getNode(nodeIds[0], pts[0].data());
    mintMesh->getNode(nodeIds[1], pts[1].data());
    mintMesh->getNode(nodeIds[2], pts[2].data());
    mintMesh->getNode(nodeIds[3], pts[3].data());

    m_tets[i] = TetrahedronType({pts[0], pts[1], pts[2], pts[3]});
  }
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

  m_tetCount = conduit::blueprint::mesh::topology::length(topoNode);

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
