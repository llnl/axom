// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/mint/mesh/Mesh.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/spin/BVH.hpp"
#include "axom/quest/TetMeshClipper.hpp"
#include "axom/bump.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/execution/scans.hpp"
#include "axom/core/WhereMacro.hpp"

namespace axom
{
namespace quest
{

TetMeshClipper::TetMeshClipper(const klee::Geometry& kGeom, const std::string& name)
  : GeometryClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("TetMesh") : name)
  , m_topoName(kGeom.getBlueprintTopology())
  , m_tetCount(0)
  , m_transformer(m_extTrans)
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

  SLIC_ERROR_IF(shapeeMesh.dimension() != 3, "TetMeshClipper requires a 3D mesh.");

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

#if 1
/*
  Alternative way:
  - Put surface triangles in BVH.
  - Create a bounding box and a ray for every hex.  The ray
    originates from the bounding box center and points away from
    the center of m_tetMeshBb.
  - Use BVH::findBoundingBoxes and BVH::findRay to get surface
    triangle near the bounding boxes and rays.  We won't need
    both for most of the hexes, but it may be faster than building
    index lists of where we need them.
  - Loop through the hexes.
    - If hex bb intersects any surface triangle bb, hex is ON.
    - Else, the hex is either IN or OUT.  It can't possibly by ON.
      Count number of surface triangles the hex's ray intersects.
      use bool intersect(const Triangle<T, 3>& tri, const Ray<T, 3>& ray)
      If count is odd, hex is IN, if even, OUT.
*/
template <typename ExecSpace>
void TetMeshClipper::labelInOutImpl(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  int allocId = shapeeMesh.getAllocatorID();
  auto cellCount = shapeeMesh.getCellCount();

  if(labels.size() < cellCount || labels.getAllocatorID() != allocId)
  {
    labels = axom::Array<LabelType>(ArrayOptions::Uninitialized(), cellCount, 0, allocId);
  }
  auto labelsView = labels.view();

  /*
    Compute surface triangles of the tet mesh.
  */
  axom::Array<Triangle3DType> surfTris = computeGeometrySurface(shapeeMesh.getRuntimePolicy(), allocId);
  writeTrianglesToVTK(surfTris, "surfaceTris.vtk");
  auto surfTrisView = surfTris.view();

  /*
    Surface triangles (as bounding boxes) in BVH.
  */
  axom::Array<BoundingBox3DType> surfTrisAsBbs(surfTris.size(), 0, allocId);
  auto surfTrisAsBbsView = surfTrisAsBbs.view();

  axom::for_all<ExecSpace>(surfTris.size(),
                           AXOM_LAMBDA(axom::IndexType fi) {
                             const auto&surfTri = surfTrisView[fi];
                             surfTrisAsBbsView[fi] = axom::primal::compute_bounding_box(surfTri);
                           });

  spin::BVH<3, ExecSpace, double> bvh;
  bvh.initialize(surfTrisAsBbsView, surfTrisAsBbsView.size());

  /*
    Compute rays.  Each ray originates from its hex center and point
    away from the center of the tet mesh.
  */
  Point3DType geomCenter = m_tetMeshBb.getCentroid(); // Estimate of tet mesh center.
  axom::ArrayView<const BoundingBox3DType> hexBbs = shapeeMesh.getCellBoundingBoxes();
  axom::ArrayView<const HexahedronType> hexes = shapeeMesh.getCellsAsHexes();
  axom::Array<Ray3DType> hexRays(hexes.size(), 0, allocId);
  auto hexRaysView = hexRays.view();
  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType cellIdx)
    {
      Point3DType hexCenter = hexBbs[cellIdx].getCentroid();
      Vector3DType direction(geomCenter, hexCenter);
      hexRaysView[cellIdx] = Ray3DType(hexCenter, direction);
    });

  /*
    Find candidate surface triangles near the cells' bounding boxes and rays.
  */
  axom::Array<IndexType> bbOffsets(cellCount, 0, allocId);
  axom::Array<IndexType> bbCounts(cellCount, 0, allocId);
  axom::Array<IndexType> bbCandidates;
  bvh.findBoundingBoxes(bbOffsets, bbCounts, bbCandidates, hexBbs.size(), hexBbs);

  axom::Array<IndexType> rayOffsets(cellCount, 0, allocId);
  axom::Array<IndexType> rayCounts(cellCount, 0, allocId);
  axom::Array<IndexType> rayCandidates;
  bvh.findRays(rayOffsets, rayCounts, rayCandidates, hexRaysView.size(), hexRaysView);

  auto bbCountsView = bbCounts.view();
  auto bbOffsetsView = bbOffsets.view();
  auto bbCandidatesView = bbCandidates.view();

  auto rayCountsView = rayCounts.view();
  auto rayOffsetsView = rayOffsets.view();
  auto rayCandidatesView = rayCandidates.view();

  const double eps = 1e-12;

  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
      LabelType& label = labelsView[cellId];

      {
        // Label cell as ON boundary if it's near the boundary.
        const auto& hexBb = hexBbs[cellId];
        auto candidateCount = bbCountsView[cellId];
        auto candidateOffset = bbOffsetsView[cellId];
        auto* candidateIds = bbCandidatesView.data() + candidateOffset;
        for(int ci = 0; ci < candidateCount; ++ci)
        {
          axom::IndexType candidateId = candidateIds[ci];
          const Triangle3DType& candidate = surfTrisView[candidateId];
          bool intersects = axom::primal::intersect(candidate, hexBb);
          if (intersects)
          {
            label = LABEL_ON;
            return;
          }
        }
      }

      /*
        At this point, the cell must be IN or OUT.  No need to account
        for the possibility that it's ON.
      */

      {
        const auto& hexRay = hexRaysView[cellId];
        auto candidateCount = rayCountsView[cellId];
        auto candidateOffset = rayOffsetsView[cellId];
        auto* candidateIds = rayCandidatesView.data() + candidateOffset;
        axom::IndexType surfaceCrossingCount = 0;
        for(int ci = 0; ci < candidateCount; ++ci)
        {
          axom::IndexType candidateId = candidateIds[ci];
          Triangle3DType candidate = surfTrisView[candidateId];
          double contactT;
          Point3DType contactPt; // contact point in unnormalized barycentric coordinates.
          bool touches = axom::primal::intersect(candidate, hexRay, contactT, contactPt);
          if (touches)
          {
            /*
              Grazing contact requires more logic and computation to
              determine if the ray crosses the boundary.  Label the
              hex as ON the boundary and to the clipping computation
              handle it this edge case.
            */
            contactPt.array() /= contactPt[0] + contactPt[1] + contactPt[2]; // Normalize
            bool grazing = axom::utilities::isNearlyEqual(contactPt[0], eps) ||
              axom::utilities::isNearlyEqual(contactPt[1], eps) ||
              axom::utilities::isNearlyEqual(contactPt[2], eps);
            if(grazing)
            {
              label = LABEL_ON;
              return;
            }
            ++surfaceCrossingCount;
          }
        }
        label = surfaceCrossingCount%2 == 0 ? LABEL_OUT : LABEL_IN;
      }
    });
}
#else
/*
  This is an old version, kept around until performance evaluation
  determines which is better.

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

  int allocId = shapeeMesh.getAllocatorID();
  auto vertCount = shapeeMesh.getVertexCount();

  /*
    Tets (as bounding boxes) in BVH.
  */
  axom::Array<BoundingBox3DType> tetsAsBbs(m_tets.size(), m_tets.size(), allocId);
  auto tetsAsBbsView = tetsAsBbs.view();
  auto tetsView = m_tets.view();

  axom::Array<TetrahedronType> tmpTets;
  if(allocId != m_tets.getAllocatorID())
  {
    tmpTets = axom::Array<TetrahedronType>(m_tets, allocId);
    tetsView = tmpTets.view();
  }
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

  auto countsView = counts.view();
  auto offsetsView = offsets.view();
  auto candidatesView = candidates.view();
  axom::for_all<ExecSpace>(
    vertCount,
    AXOM_LAMBDA(axom::IndexType vertId) {
      auto candidateCount = countsView[vertId];
      bool& isInside = vertIsInsideView[vertId];
      if (!isInside && candidateCount > 0)
      {
        auto candidateIds = &candidatesView[offsetsView[vertId]];
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
#endif

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
  axom::ArrayView<const axom::IndexType, 2> hexConnView = shapeeMesh.getCellNodeConnectivity();
  SLIC_ASSERT(hexConnView.shape() ==
              (axom::StackArray<axom::IndexType, 2> {shapeeMesh.getCellCount(), HexahedronType::NUM_HEX_VERTS}));

  if(labels.size() < shapeeMesh.getCellCount() || labels.getAllocatorID() != shapeeMesh.getAllocatorID())
  {
    labels = axom::Array<LabelType>(ArrayOptions::Uninitialized(), shapeeMesh.getCellCount(), shapeeMesh.getCellCount(), shapeeMesh.getAllocatorID());
  }
  auto labelsView = labels.view();

  axom::for_all<ExecSpace>(
    shapeeMesh.getCellCount(),
    AXOM_LAMBDA(axom::IndexType cellId) {
      auto cellVertIds = hexConnView[cellId];
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
  bpMeshGrp->importConduitTree(m_bpMesh);

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
    auto* tetVertsView = elementsGrp->getView("connectivity");
    curDim = tetVertsView->getShape(2, curShape);
    constexpr axom::IndexType NUM_VERTS_PER_TET = 4;
    SLIC_ASSERT(curDim == 1 || curDim == 2);
    if(curDim == 1)
    {
      SLIC_ASSERT(curShape[0] % NUM_VERTS_PER_TET == 0);
      axom::IndexType connShape[2] = {curShape[0] / NUM_VERTS_PER_TET, NUM_VERTS_PER_TET};
      tetVertsView->reshapeArray(2, connShape);
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
  Point3DType pts[4];
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

  m_bpMesh = m_info.fetch_existing("klee::Geometry:tetMesh");

  SLIC_ASSERT(
    m_bpMesh.fetch_existing("topologies").fetch_existing(m_topoName).fetch_existing("type").as_string() ==
    "unstructured");

  conduit::Node& topoNode = m_bpMesh.fetch_existing("topologies").fetch_existing(m_topoName);

  bool isMultiDomain = conduit::blueprint::mesh::is_multi_domain(m_bpMesh);
  SLIC_ERROR_IF(isMultiDomain, "TetMeshClipper does not support multi-domain tet meshes yet.");

  SLIC_ASSERT(conduit::blueprint::mesh::topology::dims(topoNode) == 3);

  m_tetCount = conduit::blueprint::mesh::topology::length(topoNode);

  m_coordsetName = topoNode.fetch_existing("coordset").as_string();
}

void TetMeshClipper::transformCoordset()
{
  // Apply transformations
  auto& oldCoordset = m_bpMesh.fetch_existing("coordsets").fetch_existing(m_coordsetName);
  const std::string newCoordsetName = m_coordsetName + ".trans";
  conduit::Node& coordset = m_bpMesh.fetch("coordsets")[newCoordsetName];
  coordset.set_node(oldCoordset);
  auto transformer = m_transformer;
  conduit::index_t count = conduit::blueprint::mesh::coordset::length(coordset);
  axom::ArrayView<double> xV(coordset.fetch_existing("values/x").as_double_ptr(), count);
  axom::ArrayView<double> yV(coordset.fetch_existing("values/y").as_double_ptr(), count);
  axom::ArrayView<double> zV(coordset.fetch_existing("values/z").as_double_ptr(), count);
  axom::for_all<axom::SEQ_EXEC>(count, AXOM_LAMBDA(axom::IndexType i) {
      transformer.transform(xV[i], yV[i], zV[i]);
      m_tetMeshBb.addPoint(Point3DType{xV[i], yV[i], zV[i]});
    });
  m_bpMesh.fetch_existing("topologies").fetch_existing(m_topoName).fetch_existing("coordset").set_string(newCoordsetName);
  m_coordsetName = newCoordsetName;
}

axom::Array<TetMeshClipper::Triangle3DType> TetMeshClipper::computeGeometrySurface(axom::runtime_policy::Policy policy, int allocId)
{
  AXOM_ANNOTATE_SCOPE("TetMeshClipper::computeGeometrySurface");
  switch(policy)
  {
  case axom::runtime_policy::Policy::seq:
    return computeGeometrySurface<axom::SEQ_EXEC>(allocId);
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    return computeGeometrySurface<axom::OMP_EXEC>(allocId);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  return case axom::runtime_policy::Policy::cuda:
    computeGeometrySurface<axom::CUDA_EXEC<256>>(allocId);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    return computeGeometrySurface<axom::HIP_EXEC<256>>(allocId);
    break;
#endif
  }
  SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  return {};
}

template <typename ExecSpace>
axom::Array<TetMeshClipper::Triangle3DType> TetMeshClipper::computeGeometrySurface(int allocId)
{
  /*
    Make a view of the tet mesh topology.

    Note: Using axom::IndexType for integers.  Users should configure
    Axom with the type that they plan to use.  Mixing different
    IndexTypes in the same Axom build is not supported.
  */
  namespace bumpViews = axom::bump::views;
  namespace bumpUtils = axom::bump::utilities;
  using bumpShape = bumpViews::TetShape<axom::IndexType>;

  constexpr axom::IndexType FACES_PER_TET = 4;
  constexpr axom::IndexType VERTS_PER_TET = 4;
  constexpr axom::IndexType VERTS_PER_FACE = 3;

#if 0
  conduit::Node* bpMesh = m_bpMesh;
  if(!axom::execution_space<ExecSpace>::usesAllocId(axom::ConduitMemory::conduitAllocIdToAxom(bpMesh->allocator())))
  {
    SLIC_WARNING("Cannot guarantee Conduit memory is compatible with ExecSpace."
                 "  Possible fix: (1) shallow-copy to Group.  (2) selectively"
                 " reallocate some Views using reallocateTo, with a lambda"
                 " like hostAllocForScalarAndStringViews in quest_shape_in_memory.cpp");
    // This probably won't work because string values will also be on device.
    bpMesh = new conduit::Node();
    bpMesh->set_allocator(axom::ConduitMemory::axomAllocIdToConduit(allocId));
    bpMesh->set_node(*m_bpMesh);
  }
#endif

  conduit::Node& tetMesh = m_bpMesh;
  conduit::Node& tetMeshTopo = tetMesh["topologies"][m_topoName];
  conduit::Node& tetMeshConn = tetMeshTopo.fetch_existing("elements").fetch_existing("connectivity");
  const auto tetVertsView = bumpUtils::make_array_view<axom::IndexType>(tetMeshConn);
  bumpViews::UnstructuredTopologySingleShapeView<bumpShape> topoView(tetVertsView);

  conduit::Node& polyTopo = tetMesh["topologies"]["polyhedral"];
  make_polyhedral_topology<ExecSpace>(tetMesh, polyTopo);
  const auto& faceVertConn = polyTopo.fetch_existing("subelements/connectivity");
  const axom::IndexType faceCount = faceVertConn.dtype().number_of_elements()/VERTS_PER_FACE;
  const auto faceVertsView = bumpUtils::make_array_view<axom::IndexType>(faceVertConn);
  SLIC_ASSERT(faceVertsView.size() == faceCount*VERTS_PER_FACE);
// std::cout<<__WHERE<<"tetMesh:"<<std::endl; tetMesh.to_string_stream(std::cout, "yaml", 2, 0);

  conduit::Node& polyElemFaceConn = polyTopo["elements/connectivity"];
  const auto polyElemFaceView = bumpUtils::make_array_view<axom::IndexType>(polyElemFaceConn);
  SLIC_ASSERT(polyElemFaceView.size() == m_tetCount * FACES_PER_TET);

  const std::string coordsName = polyTopo["coordset"].as_string();
  const conduit::Node& coordsNode = m_bpMesh.fetch_existing("coordsets").fetch_existing(coordsName);
  auto xs = bumpUtils::make_array_view<double>(coordsNode["values/x"]);
  auto ys = bumpUtils::make_array_view<double>(coordsNode["values/y"]);
  auto zs = bumpUtils::make_array_view<double>(coordsNode["values/z"]);

  /*
    Compute tet faces as triangles.
    Compute rays from triangle centroid, in normal direction.
  */
  axom::Array<Triangle3DType> faceTris(faceCount, faceCount, allocId);
  axom::Array<Ray3DType> faceRays(faceCount, faceCount, allocId);
  auto faceTrisView = faceTris.view();
  auto faceRaysView = faceRays.view();
  axom::for_all<ExecSpace>(faceCount,
                           AXOM_LAMBDA(axom::IndexType faceIdx)
                           {
                             axom::IndexType* vertIds = faceVertsView.data() + faceIdx*VERTS_PER_FACE;
                             auto vId0 = vertIds[0];
                             auto vId1 = vertIds[1];
                             auto vId2 = vertIds[2];
                             Point3DType v0{xs[vId0], ys[vId0], zs[vId0]};
                             Point3DType v1{xs[vId1], ys[vId1], zs[vId1]};
                             Point3DType v2{xs[vId2], ys[vId2], zs[vId2]};
                             auto& faceTri = faceTrisView[faceIdx] = Triangle3DType(v0, v1, v2);
                             faceRaysView[faceIdx] = Ray3DType(faceTri.centroid(),
                                                               faceTri.normal());
                           });
writeTrianglesToVTK(faceTris, "allCellFaces.vtk");

  const auto tetFaceConn = polyTopo["elements/connectivity"];
  const auto tetFacesView = bumpUtils::make_array_view<axom::IndexType>(tetFaceConn);
  SLIC_ASSERT(tetFacesView.size() == 4*m_tetCount);

  /*
    Compute whether faces have tets on each side.
  */
  axom::Array<bool> hasCellOnFrontSide(faceCount, 0, allocId);
  axom::Array<bool> hasCellOnBackSide(faceCount, 0, allocId);
  hasCellOnFrontSide.fill(false);
  hasCellOnBackSide.fill(false);
  auto hasCellOnFrontSideView = hasCellOnFrontSide.view();
  auto hasCellOnBackSideView = hasCellOnBackSide.view();
  axom::for_all<ExecSpace>(
    m_tetCount, AXOM_LAMBDA(IndexType tetId) {
      Point3DType cellCentroid({0.0, 0.0, 0.0});
      axom::IndexType* vertIdxs = tetVertsView.data() + tetId*VERTS_PER_TET;
      TetrahedronType tet;
      for(int vi = 0; vi < VERTS_PER_TET; ++vi)
      {
        axom::IndexType vIdx = vertIdxs[vi];
        Point3DType vertCoords{xs[vIdx], ys[vIdx], zs[vIdx]};
        cellCentroid.array() += vertCoords.array();
        tet[vi] = vertCoords;
      }
      cellCentroid.array() /= VERTS_PER_TET;

      for(int fi = 0; fi < FACES_PER_TET; ++fi)
      {
        axom::IndexType faceIdx = polyElemFaceView[tetId*FACES_PER_TET + fi];
        const auto& faceTri = faceTrisView[faceIdx];
        Ray3DType faceRay(faceTri.centroid(), faceTri.normal());
        const Vector3DType& faceDir = faceRay.direction();
        auto cellToFace = cellCentroid.array() - faceRay.origin().array();
        double distParam = faceDir[0]*cellToFace[0] + faceDir[1]*cellToFace[1] + faceDir[2]*cellToFace[2];
        bool &hasCell = distParam >= 0 ? hasCellOnFrontSideView[faceIdx] : hasCellOnBackSideView[faceIdx];
        hasCell = true;
      }
    });

  /*
    Mark faces touching only 1 cell.  Reuse hasCellOnFrontSide space.
  */
  axom::Array<axom::IndexType> hasCellOnOneSide(ArrayOptions::Uninitialized(), faceCount, 0, allocId);
  auto hasCellOnOneSideView = hasCellOnOneSide.view();
  axom::for_all<ExecSpace>(
    faceCount,
    AXOM_LAMBDA(IndexType faceId) {
      hasCellOnOneSideView[faceId] =
        (hasCellOnFrontSideView[faceId] + hasCellOnBackSideView[faceId]) == 1;
    });

  /*
    Get running total of surface triangle count using prefix-sum scan.
    Then use the results to populate array of those faces.
  */
  axom::Array<axom::IndexType> prefixSum(faceCount + 1, 0, allocId);
  prefixSum.fill(0);
  auto prefixSumView = prefixSum.view();
  axom::inclusive_scan<ExecSpace>(
    hasCellOnOneSide,
    axom::ArrayView<axom::IndexType>(prefixSum.data() + 1, faceCount));

  axom::IndexType surfFaceCount = -1;
  axom::copy(&surfFaceCount,
             prefixSumView.data() + prefixSumView.size() - 1,
             sizeof(surfFaceCount));

  axom::Array<axom::IndexType> surfFaceIds(surfFaceCount, 0, allocId);
  axom::Array<Triangle3DType> surfTris(surfFaceCount, 0, allocId);
  auto surfFaceIdsView = surfFaceIds.view();
  auto surfTrisView = surfTris.view();
  axom::for_all<ExecSpace>(
    faceCount,
    AXOM_LAMBDA(axom::IndexType faceIdx)
    {
      if(prefixSumView[faceIdx] != prefixSumView[faceIdx + 1])
      {
        auto runningTotal = prefixSumView[faceIdx];
        surfFaceIdsView[runningTotal] = faceIdx;
        surfTrisView[runningTotal] = faceTrisView[faceIdx];
      }
    });

#if 0
  if(bpMesh != m_bpMesh) { delete bpMesh; } // Delete temporary.
#endif

  return surfTris;
}

void TetMeshClipper::writeTrianglesToVTK(
  const axom::Array<Triangle3DType>& triangles,
  const std::string& filename)
{
  std::ofstream ofs(filename);
  if(!ofs) {
    std::cerr << "Cannot open file for writing: " << filename << std::endl;
    return;
  }

  // Header
  ofs << "# vtk DataFile Version 3.0\n";
  ofs << "Triangle mesh\n";
  ofs << "ASCII\n";
  ofs << "DATASET POLYDATA\n";

  // Write points
  ofs << "POINTS " << triangles.size()*3 << " double\n";
  for(const auto& tri : triangles)
  {
    for(int i = 0; i < 3; ++i)
    {
      const auto& pt = tri[i];
      ofs << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
    }
  }

  // Write polygons (triangles)
  ofs << "POLYGONS " << triangles.size() << " " << triangles.size() * 4 << "\n";
  for(int i = 0; i < triangles.size(); ++i)
  {
    ofs << "3 " << i*3 << " " << i*3 + 1 << " " << i*3+2 << "\n";
  }

  ofs.close();
}

/*
  @param tetMesh Input unstructured tet mesh, single domain.
  @param polyTopo Output unstructured polyhedral topology.
*/

template<typename ExecSpace>
void TetMeshClipper::make_polyhedral_topology(
  const conduit::Node& tetMesh,
  conduit::Node& polyTopo)
{
  namespace bumpViews = axom::bump::views;
  namespace bumpUtils = axom::bump::utilities;
  using bumpShape = bumpViews::TetShape<axom::IndexType>;

  const conduit::Node& tetTopo = tetMesh["topologies"][m_topoName];
  SLIC_ASSERT( tetTopo["type"].as_string() == std::string("unstructured") );
  auto tetTopoView = bumpViews::make_unstructured_single_shape_topology<bumpShape>::view(tetTopo);
  using TopologyView = decltype(tetTopoView);
  using ConnectivityType = typename TopologyView::ConnectivityType;

  bump::MakePolyhedralTopology<ExecSpace, TopologyView> polyTopoMaker(tetTopoView);
  polyTopoMaker.execute(tetTopo, polyTopo);
std::cout<<__WHERE<<"polyTopo before face merge:"<<std::endl;
polyTopo.to_string_stream(std::cout, "yaml", 2, 0);
  bump::MergePolyhedralFaces<ExecSpace, ConnectivityType>::execute(polyTopo);
}

// Run cellKernel through a cell loop.
// Run faceKernel through a face loop.

}  // end namespace quest
}  // end namespace axom
