// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/core.hpp"
#include "axom/spin/BVH.hpp"
#include "axom/quest/detail/clipping/TetMeshClipper.hpp"
#include "axom/bump.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/execution/scans.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{

TetMeshClipper::TetMeshClipper(const klee::Geometry& kGeom, const std::string& name)
  : MeshClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("TetMesh") : name)
  , m_topoName(kGeom.getBlueprintTopology())
  , m_tetCount(0)
  , m_transformer(m_extTrans)
{
  SLIC_ASSERT(!m_topoName.empty());

  extractClipperInfo();

  transformCoordset();
}

bool TetMeshClipper::labelCellsInOut(quest::experimental::ShapeMesh& shapeMesh,
                                     axom::Array<LabelType>& labels)
{
  SLIC_ERROR_IF(shapeMesh.dimension() != 3, "TetMeshClipper requires a 3D mesh.");

  int allocId = shapeMesh.getAllocatorID();
  auto cellCount = shapeMesh.getCellCount();
  if(labels.size() < cellCount || labels.getAllocatorID() != allocId)
  {
    labels = axom::Array<LabelType>(ArrayOptions::Uninitialized(), cellCount, 0, allocId);
  }

  switch(shapeMesh.getRuntimePolicy())
  {
  case axom::runtime_policy::Policy::seq:
    labelCellsInOutImpl<axom::SEQ_EXEC>(shapeMesh, labels.view());
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    labelCellsInOutImpl<axom::OMP_EXEC>(shapeMesh, labels.view());
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case axom::runtime_policy::Policy::cuda:
    labelCellsInOutImpl<axom::CUDA_EXEC<256>>(shapeMesh, labels.view());
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    labelCellsInOutImpl<axom::HIP_EXEC<256>>(shapeMesh, labels.view());
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
  return true;
}

bool TetMeshClipper::labelTetsInOut(quest::experimental::ShapeMesh& shapeMesh,
                                    axom::ArrayView<const axom::IndexType> cellIds,
                                    axom::Array<LabelType>& tetLabels)
{
  SLIC_ERROR_IF(shapeMesh.dimension() != 3, "TetMeshClipper requires a 3D mesh.");

  int allocId = shapeMesh.getAllocatorID();
  const axom::IndexType tetCount = cellIds.size() * NUM_TETS_PER_HEX;
  if(tetLabels.size() < tetCount || tetLabels.getAllocatorID() != allocId)
  {
    tetLabels = axom::Array<LabelType>(ArrayOptions::Uninitialized(), tetCount, 0, allocId);
  }

  switch(shapeMesh.getRuntimePolicy())
  {
  case axom::runtime_policy::Policy::seq:
    labelTetsInOutImpl<axom::SEQ_EXEC>(shapeMesh, cellIds, tetLabels.view());
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    labelTetsInOutImpl<axom::OMP_EXEC>(shapeMesh, cellIds, tetLabels.view());
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case axom::runtime_policy::Policy::cuda:
    labelTetsInOutImpl<axom::CUDA_EXEC<256>>(shapeMesh, cellIds, tetLabels.view());
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    labelTetsInOutImpl<axom::HIP_EXEC<256>>(shapeMesh, cellIds, tetLabels.view());
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
  return true;
}

/*
 * To determine whether a cell is completely inside the tet mesh or
 * completely outside, create a ray from the cell interior in some
 * direction.  If the ray crosses the tet mesh boundary an odd number
 * of times, it's inside; if even, it's outside.  If can't efficiently
 * determined, call it on the boundary so the slower clipping function
 * can compute the intersection volume.
 *
 * We use a BVH to find candidate boundary facets that a ray can potentially
 * intersect.
 */
template <typename ExecSpace>
void TetMeshClipper::labelCellsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                                         axom::ArrayView<LabelType> labels)
{
  int allocId = shapeMesh.getAllocatorID();
  auto cellCount = shapeMesh.getCellCount();

  axom::Array<Triangle3DType> surfTris;
  spin::BVH<3, ExecSpace, double> bvh;
  computeSurfaceTrianglesAndBVH<ExecSpace>(allocId, surfTris, bvh);
  auto surfTrisView = surfTris.view();

  axom::Array<Ray3DType> hexRays;
  computeHexRays<ExecSpace>(shapeMesh, hexRays);
  auto hexRaysView = hexRays.view();

  /*
   * Find candidate surface triangles near the cells' bounding boxes and rays.
   */
  axom::ArrayView<const BoundingBox3DType> hexBbs = shapeMesh.getCellBoundingBoxes();
  AXOM_ANNOTATE_BEGIN("TetMeshClipper::get_surf_near_bbs");
  axom::Array<IndexType> bbOffsets(cellCount, 0, allocId);
  axom::Array<IndexType> bbCounts(cellCount, 0, allocId);
  axom::Array<IndexType> bbCandidates;
  bvh.findBoundingBoxes(bbOffsets, bbCounts, bbCandidates, hexBbs.size(), hexBbs);
  AXOM_ANNOTATE_END("TetMeshClipper::get_surf_near_bbs");

  AXOM_ANNOTATE_BEGIN("TetMeshClipper::get_surf_near_rays");
  axom::Array<IndexType> rayOffsets(cellCount, 0, allocId);
  axom::Array<IndexType> rayCounts(cellCount, 0, allocId);
  axom::Array<IndexType> rayCandidates;
  bvh.findRays(rayOffsets, rayCounts, rayCandidates, hexRaysView.size(), hexRaysView);
  AXOM_ANNOTATE_END("TetMeshClipper::get_surf_near_rays");

  auto bbCountsView = bbCounts.view();
  auto bbOffsetsView = bbOffsets.view();
  auto bbCandidatesView = bbCandidates.view();

  auto rayCountsView = rayCounts.view();
  auto rayOffsetsView = rayOffsets.view();
  auto rayCandidatesView = rayCandidates.view();

  const double EPS = 1e-12;

  AXOM_ANNOTATE_BEGIN("TetMeshClipper::compute_labels");
  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
      LabelType& label = labels[cellId];

      {
        // Label cell as ON boundary if it's bounding box intersects the boundary.
        const auto& hexBb = hexBbs[cellId];
        auto candidateCount = bbCountsView[cellId];
        auto candidateOffset = bbOffsetsView[cellId];
        auto* candidateIds = bbCandidatesView.data() + candidateOffset;
        for(int ci = 0; ci < candidateCount; ++ci)
        {
          axom::IndexType candidateId = candidateIds[ci];
          const Triangle3DType& candidate = surfTrisView[candidateId];
          bool intersects = axom::primal::intersect(candidate, hexBb);
          if(intersects)
          {
            label = LabelType::LABEL_ON;
            return;
          }
        }
      }

      /*
       * At this point, the cell must be IN or OUT.  No need to account
       * for the possibility that it's ON.
       */

      {
        const Ray3DType& hexRay = hexRaysView[cellId];
        auto candidateCount = rayCountsView[cellId];
        auto candidateOffset = rayOffsetsView[cellId];
        auto* candidateIds = rayCandidatesView.data() + candidateOffset;
        axom::IndexType surfaceCrossingCount = 0;
        for(int ci = 0; ci < candidateCount; ++ci)
        {
          axom::IndexType candidateId = candidateIds[ci];
          Triangle3DType candidate = surfTrisView[candidateId];
          double contactT;
          Point3DType contactPt;  // contact point in unnormalized barycentric coordinates.
          bool touches = axom::primal::intersect(candidate, hexRay, contactT, contactPt);
          if(touches)
          {
            /*
             * Grazing contact requires more logic and computation to
             * determine if the ray crosses the boundary.  Label the
             * hex as ON the boundary to have the clipping computation
             * handle this edge case.
             */
            contactPt.array() /= contactPt[0] + contactPt[1] + contactPt[2];  // Normalize
            bool grazing = axom::utilities::isNearlyEqual(contactPt[0], EPS) ||
              axom::utilities::isNearlyEqual(contactPt[1], EPS) ||
              axom::utilities::isNearlyEqual(contactPt[2], EPS);
            if(grazing)
            {
              label = LabelType::LABEL_ON;
              return;
            }
            ++surfaceCrossingCount;
          }
        }
        label = surfaceCrossingCount % 2 == 0 ? LabelType::LABEL_OUT : LabelType::LABEL_IN;
      }
    });
  AXOM_ANNOTATE_END("TetMeshClipper::compute_labels");
}

/*
 * To determine whether a tet is completely inside the tet mesh or
 * completely outside, create a ray from the tet center in some
 * direction.  If the ray crosses the tet mesh boundary an odd number
 * of times, it's inside; if even, it's outside.  If can't efficiently
 * determined, call it on the boundary so the slower clipping function
 * can compute the intersection volume.
 *
 * We use a BVH to find candidate boundary facets that a ray can potentially
 * intersect.
 */
template <typename ExecSpace>
void TetMeshClipper::labelTetsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                                        axom::ArrayView<const axom::IndexType> cellIds,
                                        axom::ArrayView<LabelType> tetLabels)
{
  int allocId = shapeMesh.getAllocatorID();
  auto cellCount = cellIds.size();
  auto tetCount = cellCount * NUM_TETS_PER_HEX;

  auto tetVolumes = shapeMesh.getTetVolumes();

  // Copy m_tetMesh array data to allocId if it's not done yet.
  copy_topo_and_coords_to(allocId);

  axom::Array<Triangle3DType> surfTris;
  spin::BVH<3, ExecSpace, double> bvh;
  computeSurfaceTrianglesAndBVH<ExecSpace>(allocId, surfTris, bvh);
  auto surfTrisView = surfTris.view();

  axom::Array<BoundingBox3DType> tetBbs;
  axom::Array<Ray3DType> tetRays;
  computeTetRays<ExecSpace>(shapeMesh, cellIds, tetRays, tetBbs);
  auto tetRaysView = tetRays.view();
  auto tetBbsView = tetBbs.view();

  /*
   * Find candidate surface triangles near the tets' bounding boxes and rays.
   */
  AXOM_ANNOTATE_BEGIN("TetMeshClipper::get_surf_near_bbs");
  axom::Array<IndexType> bbOffsets(tetCount, 0, allocId);
  axom::Array<IndexType> bbCounts(tetCount, 0, allocId);
  axom::Array<IndexType> bbCandidates;
  bvh.findBoundingBoxes(bbOffsets, bbCounts, bbCandidates, tetBbs.size(), tetBbsView);
  AXOM_ANNOTATE_END("TetMeshClipper::get_surf_near_bbs");

  AXOM_ANNOTATE_BEGIN("TetMeshClipper::get_surf_near_rays");
  axom::Array<IndexType> rayOffsets(tetCount, 0, allocId);
  axom::Array<IndexType> rayCounts(tetCount, 0, allocId);
  axom::Array<IndexType> rayCandidates;
  bvh.findRays(rayOffsets, rayCounts, rayCandidates, tetRaysView.size(), tetRaysView);
  AXOM_ANNOTATE_END("TetMeshClipper::get_surf_near_rays");

  auto bbCountsView = bbCounts.view();
  auto bbOffsetsView = bbOffsets.view();
  auto bbCandidatesView = bbCandidates.view();

  auto rayCountsView = rayCounts.view();
  auto rayOffsetsView = rayOffsets.view();
  auto rayCandidatesView = rayCandidates.view();

  const double EPS = 1e-12;

  AXOM_ANNOTATE_BEGIN("TetMeshClipper::compute_labels");
  axom::for_all<ExecSpace>(
    tetCount,
    AXOM_LAMBDA(axom::IndexType ti) {
      LabelType& label = tetLabels[ti];

      axom::IndexType ci = ti / NUM_TETS_PER_HEX;
      axom::IndexType tii = ti % NUM_TETS_PER_HEX;
      axom::IndexType tetId = ci * NUM_TETS_PER_HEX + tii;
      if(axom::utilities::isNearlyEqual(tetVolumes[tetId], 0.0, EPS))
      {
        label = LabelType::LABEL_OUT;
        return;
      }

      {
        // Label tet as ON boundary if it's near the boundary.
        const auto& tetBb = tetBbsView[ti];
        auto candidateCount = bbCountsView[ti];
        auto candidateOffset = bbOffsetsView[ti];
        auto* candidateIds = bbCandidatesView.data() + candidateOffset;
        for(int ci = 0; ci < candidateCount; ++ci)
        {
          axom::IndexType candidateId = candidateIds[ci];
          const Triangle3DType& candidate = surfTrisView[candidateId];
          bool intersects = axom::primal::intersect(candidate, tetBb);
          if(intersects)
          {
            label = LabelType::LABEL_ON;
            return;
          }
        }
      }

      /*
       * At this point, the cell must be IN or OUT.  No need to account
       * for the possibility that it's ON.
       */

      {
        const Ray3DType& tetRay = tetRaysView[ti];
        auto candidateCount = rayCountsView[ti];
        auto candidateOffset = rayOffsetsView[ti];
        auto* candidateIds = rayCandidatesView.data() + candidateOffset;
        axom::IndexType surfaceCrossingCount = 0;
        for(int ci = 0; ci < candidateCount; ++ci)
        {
          axom::IndexType candidateId = candidateIds[ci];
          Triangle3DType candidate = surfTrisView[candidateId];
          double contactT;
          Point3DType contactPt;  // contact point in unnormalized barycentric coordinates.
          bool touches = axom::primal::intersect(candidate, tetRay, contactT, contactPt);
          if(touches)
          {
            /*
             * Grazing contact requires more logic and computation to
             * determine if the ray crosses the boundary.  Label the
             * tet as ON the boundary so the clipping computation
             * will handle this edge case.
             */
            contactPt.array() /= contactPt[0] + contactPt[1] + contactPt[2];  // Normalize
            bool grazing = axom::utilities::isNearlyEqual(contactPt[0], EPS) ||
              axom::utilities::isNearlyEqual(contactPt[1], EPS) ||
              axom::utilities::isNearlyEqual(contactPt[2], EPS);
            if(grazing)
            {
              label = LabelType::LABEL_ON;
              return;
            }
            ++surfaceCrossingCount;
          }
        }
        label = surfaceCrossingCount % 2 == 0 ? LabelType::LABEL_OUT : LabelType::LABEL_IN;
      }
    });
  AXOM_ANNOTATE_END("TetMeshClipper::compute_labels");
}

template <typename ExecSpace>
void TetMeshClipper::computeHexRays(quest::experimental::ShapeMesh& shapeMesh,
                                    axom::Array<Ray3DType>& hexRays)
{
  AXOM_ANNOTATE_SCOPE("TetMeshClipper::computeHexRays");
  /*
   * Each ray originates from a hex interior point and oriented
   * away from the center of the tet mesh.  To get an interior
   * point for a potentially concave hex, find the biggest of the
   * tets the hex decomposes into and use that tet's centroid.
   */
  Point3DType geomCenter = m_tetMeshBb.getCentroid();  // Estimate of tet mesh center.
  auto meshHexes = shapeMesh.getCellsAsHexes();
  auto cellCount = shapeMesh.getCellCount();
  hexRays = axom::Array<Ray3DType>(cellCount, 0, shapeMesh.getAllocatorID());
  auto hexRaysView = hexRays.view();
  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType cellIdx) {
      TetrahedronType cellTets[ShapeMesh::NUM_TETS_PER_HEX];
      ShapeMesh::hexToTets(meshHexes[cellIdx], cellTets);
      IndexType iMaxTetVol = 0;
      double maxTetVol = cellTets[iMaxTetVol].signedVolume();
      for(IndexType i = 1; i < ShapeMesh::NUM_TETS_PER_HEX; ++i)
      {
        double tmpTetVol = cellTets[i].signedVolume();
        if(maxTetVol < tmpTetVol)
        {
          maxTetVol = tmpTetVol;
          iMaxTetVol = i;
        }
      }
      const auto& bigTet = cellTets[iMaxTetVol];
      Point3DType hexInterior(
        (bigTet[0].array() + bigTet[1].array() + bigTet[2].array() + bigTet[3].array()) / 4.0);
      Vector3DType direction(geomCenter, hexInterior);
      if(geomCenter == hexInterior)
      {
        direction = Vector3DType({1.0, 0.0, 0.0});
      }
      hexRaysView[cellIdx] = Ray3DType(hexInterior, direction);
    });
}

template <typename ExecSpace>
void TetMeshClipper::computeTetRays(quest::experimental::ShapeMesh& shapeMesh,
                                    axom::ArrayView<const axom::IndexType> cellIds,
                                    axom::Array<Ray3DType>& tetRays,
                                    axom::Array<BoundingBox3DType>& tetBbs)
{
  AXOM_ANNOTATE_SCOPE("TetMeshClipper::computeTetRays");
  int allocId = shapeMesh.getAllocatorID();
  auto cellCount = cellIds.size();
  auto tetCount = cellCount * NUM_TETS_PER_HEX;

  /*
   * Each ray originates from a tet center and point
   * away from the center of the tet mesh.
   */
  Point3DType geomCenter = m_tetMeshBb.getCentroid();  // Estimate of tet mesh center.
  tetBbs = axom::Array<BoundingBox3DType>(axom::ArrayOptions::Uninitialized(), tetCount, 0, allocId);
  auto tetBbsView = tetBbs.view();
  tetRays = axom::Array<Ray3DType>(axom::ArrayOptions::Uninitialized(), tetCount, 0, allocId);
  auto tetRaysView = tetRays.view();
  const auto meshTets = shapeMesh.getCellsAsTets();
  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType ci) {
      auto cellId = cellIds[ci];
      auto* tetsForCell = &meshTets[cellId * NUM_TETS_PER_HEX];
      for(int ti = 0; ti < NUM_TETS_PER_HEX; ++ti)
      {
        const auto& tet = tetsForCell[ti];
        Point3DType tetCenter((tet[0].array() + tet[1].array() + tet[2].array() + tet[3].array()) /
                              4.0);
        Vector3DType direction(geomCenter, tetCenter);
        if(geomCenter == tetCenter)
        {
          direction = Vector3DType({1.0, 0.0, 0.0});
        }
        tetRaysView[ci * NUM_TETS_PER_HEX + ti] = Ray3DType(tetCenter, direction);
        tetBbsView[ci * NUM_TETS_PER_HEX + ti] = BoundingBox3DType {tet[0], tet[1], tet[2]};
      }
    });
}

template <typename ExecSpace>
void TetMeshClipper::computeSurfaceTrianglesAndBVH(int allocId,
                                                   axom::Array<Triangle3DType>& surfTris,
                                                   spin::BVH<3, ExecSpace, double>& bvh)
{
  /*
    Compute surface triangles of the tet mesh.
  */
  AXOM_ANNOTATE_BEGIN("TetMeshClipper:compute_surface");
  surfTris = computeGeometrySurface<ExecSpace>(allocId);
  AXOM_ANNOTATE_END("TetMeshClipper:compute_surface");
  auto surfTrisView = surfTris.view();

  /*
    Surface triangles (as bounding boxes) in BVH.
  */
  AXOM_ANNOTATE_BEGIN("TetMeshClipper::make_surf_bvh");
  axom::Array<BoundingBox3DType> surfTrisAsBbs(surfTris.size(), 0, allocId);
  auto surfTrisAsBbsView = surfTrisAsBbs.view();

  axom::for_all<ExecSpace>(
    surfTris.size(),
    AXOM_LAMBDA(axom::IndexType fi) {
      const auto& surfTri = surfTrisView[fi];
      surfTrisAsBbsView[fi] = axom::primal::compute_bounding_box(surfTri);
    });

  bvh.initialize(surfTrisAsBbsView, surfTrisAsBbsView.size());
  AXOM_ANNOTATE_END("TetMeshClipper::make_surf_bvh");
}

/*
 * Label cell outside if no vertex is in the bounding box.
 * Otherwise, label it on boundary, because we don't know.
 */
template <typename ExecSpace>
void TetMeshClipper::vertexInsideToCellLabel(quest::experimental::ShapeMesh& shapeMesh,
                                             axom::ArrayView<bool>& vertIsInside,
                                             axom::Array<LabelType>& labels)
{
  axom::ArrayView<const axom::IndexType, 2> hexConnView = shapeMesh.getCellNodeConnectivity();
  SLIC_ASSERT(hexConnView.shape() ==
              (axom::StackArray<axom::IndexType, 2> {shapeMesh.getCellCount(),
                                                     HexahedronType::NUM_HEX_VERTS}));

  if(labels.size() < shapeMesh.getCellCount() ||
     labels.getAllocatorID() != shapeMesh.getAllocatorID())
  {
    labels = axom::Array<LabelType>(ArrayOptions::Uninitialized(),
                                    shapeMesh.getCellCount(),
                                    shapeMesh.getCellCount(),
                                    shapeMesh.getAllocatorID());
  }
  auto labelsView = labels.view();

  axom::for_all<ExecSpace>(
    shapeMesh.getCellCount(),
    AXOM_LAMBDA(axom::IndexType cellId) {
      auto cellVertIds = hexConnView[cellId];

      bool hasIn = vertIsInside[cellVertIds[0]];
      bool hasOut = !hasIn;
      for(int vi = 1; vi < HexahedronType::NUM_HEX_VERTS; ++vi)
      {
        int vertId = cellVertIds[vi];
        bool isIn = vertIsInside[vertId];
        hasIn |= isIn;
        hasOut |= !isIn;
      }
      labelsView[cellId] = !hasOut ? LabelType::LABEL_IN
        : !hasIn                   ? LabelType::LABEL_OUT
                                   : LabelType::LABEL_ON;
    });

  return;
}

bool TetMeshClipper::getGeometryAsTets(quest::experimental::ShapeMesh& shapeMesh,
                                       axom::Array<TetrahedronType>& tets)
{
  int allocId = shapeMesh.getAllocatorID();
  if(tets.size() < m_tetCount || tets.getAllocatorID() != allocId)
  {
    tets = axom::Array<TetrahedronType>(ArrayOptions::Uninitialized(), m_tetCount, 0, allocId);
  }

  switch(shapeMesh.getRuntimePolicy())
  {
  case axom::runtime_policy::Policy::seq:
    computeTets<axom::SEQ_EXEC>(tets.view());
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    computeTets<axom::OMP_EXEC>(tets.view());
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case axom::runtime_policy::Policy::cuda:
    computeTets<axom::CUDA_EXEC<256>>(tets.view());
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    computeTets<axom::HIP_EXEC<256>>(tets.view());
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }

  return true;
}

/*
 * Copy the tets from m_tetMesh tetsView.
 */
template <typename ExecSpace>
void TetMeshClipper::computeTets(axom::ArrayView<TetrahedronType> tetsView)
{
  AXOM_ANNOTATE_SCOPE("TetMeshClipper::computeTets");

  const int allocId = tetsView.getAllocatorID();
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(allocId));
  SLIC_ASSERT(tetsView.size() == m_tetCount);

  copy_topo_and_coords_to(allocId);

  const std::string topoName = axom::fmt::format("{}.{}", m_topoName, allocId);
  conduit::Node& tetTopo = m_tetMesh["topologies"][topoName];
  conduit::Node& tetVertConn = tetTopo.fetch_existing("elements").fetch_existing("connectivity");
  auto expectedConnectivityType = axom::sidre::detail::SidreTT<axom::IndexType>::id;
  SLIC_ERROR_IF(
    tetVertConn.dtype().id() != expectedConnectivityType,
    fmt::format(
      "Tet mesh connectivity type {} doesn't match Axom IndexType {}.  Please use the correct "
      "connectivity type for tet mesh or configure Axom with a compatible IndexType.",
      tetVertConn.dtype().id(),
      expectedConnectivityType));
  constexpr IndexType VERTS_PER_TET = 4;
  axom::ArrayView<const IndexType, 2> tetVertConnView(
    static_cast<const IndexType*>(tetVertConn.data_ptr()),
    m_tetCount,
    VERTS_PER_TET);

  const std::string coordsName = tetTopo["coordset"].as_string();
  const conduit::Node& coordsNode = m_tetMesh.fetch_existing("coordsets").fetch_existing(coordsName);
  namespace bumpUtils = axom::bump::utilities;
  auto xs = bumpUtils::make_array_view<double>(coordsNode["values/x"]);
  auto ys = bumpUtils::make_array_view<double>(coordsNode["values/y"]);
  auto zs = bumpUtils::make_array_view<double>(coordsNode["values/z"]);

  axom::for_all<ExecSpace>(
    m_tetCount,
    AXOM_LAMBDA(IndexType ti) {
      auto vertices = tetVertConnView[ti];
      Point3DType a({xs[vertices[0]], ys[vertices[0]], zs[vertices[0]]});
      Point3DType b({xs[vertices[1]], ys[vertices[1]], zs[vertices[1]]});
      Point3DType c({xs[vertices[2]], ys[vertices[2]], zs[vertices[2]]});
      Point3DType d({xs[vertices[3]], ys[vertices[3]], zs[vertices[3]]});
      tetsView[ti] = TetrahedronType(a, b, c, d);
    });
}

void TetMeshClipper::extractClipperInfo()
{
  m_topoName = m_info.fetch_existing("topologyName").as_string();

  m_tetMesh = m_info.fetch_existing("klee::Geometry:tetMesh");

  SLIC_ASSERT(
    m_tetMesh.fetch_existing("topologies").fetch_existing(m_topoName).fetch_existing("type").as_string() ==
    "unstructured");

  if(m_tetMesh.has_child("fields") && m_tetMesh["fields"].number_of_children() == 0)
  {
    // Remove empty "fields" node:  Blueprint check will fail if "fields" is empty,
    m_tetMesh.remove_child("fields");
  }

  {
    std::string whyBad;
    bool good = isValidTetMesh(m_tetMesh, whyBad);
    if(!good)
    {
      SLIC_ERROR(axom::fmt::format("TetMeshClipper given bad tet mesh: {}", whyBad));
    }
  }

  conduit::Node& topoNode = m_tetMesh.fetch_existing("topologies").fetch_existing(m_topoName);

  bool isMultiDomain = conduit::blueprint::mesh::is_multi_domain(m_tetMesh);
  SLIC_ERROR_IF(isMultiDomain, "TetMeshClipper does not support multi-domain tet meshes yet.");

  SLIC_ASSERT(conduit::blueprint::mesh::topology::dims(topoNode) == 3);
  SLIC_ASSERT(topoNode.fetch_existing("elements/shape").as_string() == "tet");

  m_tetCount = conduit::blueprint::mesh::topology::length(topoNode);

  m_coordsetName = topoNode.fetch_existing("coordset").as_string();

  bool fixOrientation = false;
  if(m_info.has_child("fixOrientation"))
  {
    fixOrientation = bool(m_info.fetch_existing("fixOrientation").as_int());
  }

#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  using HOST_EXEC = axom::OMP_EXEC;
#else
  using HOST_EXEC = axom::SEQ_EXEC;
#endif
  auto& connNode = topoNode.fetch_existing("elements/connectivity");
  if(connNode.dtype().id() != bump::utilities::cpp2conduit<IndexType>::id)
  {
    // Change connectivity index type to IndexType.
    const auto srcIndexId = connNode.dtype().id();
    const IndexType connLen = connNode.dtype().number_of_elements();
    conduit::DataType dstConnType(axom::bump::utilities::cpp2conduit<IndexType>::id, connLen);
    conduit::Node tmpConnNode(dstConnType);
    IndexType* newPtr = static_cast<IndexType*>(tmpConnNode.data_ptr());

    if(connNode.dtype().id() == bump::utilities::cpp2conduit<std::int32_t>::id)
    {
      auto curPtr = connNode.as_int32_ptr();
      axom::for_all<HOST_EXEC>(connLen, AXOM_LAMBDA(IndexType i) { newPtr[i] = curPtr[i]; });
    }
    else if(connNode.dtype().id() == bump::utilities::cpp2conduit<std::int64_t>::id)
    {
      auto curPtr = connNode.as_int64_ptr();
      axom::for_all<HOST_EXEC>(connLen, AXOM_LAMBDA(IndexType i) { newPtr[i] = curPtr[i]; });
    }
    else
    {
      SLIC_ERROR(axom::fmt::format("Unrecognized connectivity index type with conduit type id {}",
                                   srcIndexId));
    }
    connNode.swap(tmpConnNode);
  }

  auto& coordSet = m_tetMesh.fetch_existing("coordsets").fetch_existing(m_coordsetName);
  checkTetOrientations<HOST_EXEC>(connNode, coordSet, fixOrientation);
}

template <typename ExecSpace>
void TetMeshClipper::checkTetOrientations(conduit::Node& connNode,
                                          const conduit::Node& coordSet,
                                          bool fixOrientation)
{
  const IndexType tetCount = connNode.dtype().number_of_elements() / 4;
  axom::ArrayView<IndexType, 2> connArray(static_cast<IndexType*>(connNode.data_ptr()),
                                          axom::StackArray<IndexType, 2> {tetCount, 4});

  const IndexType pointCount = coordSet.fetch_existing("values/x").dtype().number_of_elements();
  axom::ArrayView<const double> x(coordSet.fetch_existing("values/x").as_double_ptr(), pointCount);
  axom::ArrayView<const double> y(coordSet.fetch_existing("values/y").as_double_ptr(), pointCount);
  axom::ArrayView<const double> z(coordSet.fetch_existing("values/z").as_double_ptr(), pointCount);

  constexpr double eps = 1e-12;
  axom::for_all<ExecSpace>(
    tetCount,
    AXOM_LAMBDA(IndexType iTet) {
      Point3DType verts[4];
      for(int vi = 0; vi < 4; ++vi)
      {
        int iPoint = connArray(iTet, vi);
        verts[vi] = Point3DType {axom::NumericArray<double, 3> {x[iPoint], y[iPoint], z[iPoint]}};
      }

      TetrahedronType tet(verts[0], verts[1], verts[2], verts[3]);
      const auto signedVol = tet.signedVolume();
      if(signedVol < -eps)
      {
        if(fixOrientation)
        {
          axom::utilities::swap(connArray(iTet, 2), connArray(iTet, 3));
        }
        else
        {
#if defined(AXOM_DEVICE_CODE)
          SLIC_ASSERT(!(signedVol < -eps));
#else
          SLIC_ERROR(axom::fmt::format(
            "TetMeshClipper: input tet[{}] is inverted, with volume {}.  Please fix or add "
            "'fixOrientation: 1' to the input klee::Geometry's hierarchical data.",
            iTet,
            signedVol));
#endif
        }
      }
    });
}

bool TetMeshClipper::isValidTetMesh(const conduit::Node& tetMesh, std::string& whyBad) const
{
  bool rval = true;

  conduit::Node info;
  rval = conduit::blueprint::mesh::verify(tetMesh, info);

  if(rval)
  {
    std::string topoType = tetMesh.fetch("topologies")[m_topoName]["type"].as_string();
    rval = topoType == "unstructured";
    info[0].set_string("Topology is not unstructured.");
  }

  if(rval)
  {
    std::string elemShape = tetMesh.fetch("topologies")[m_topoName]["elements/shape"].as_string();
    rval = elemShape == "tet";
    info[0].set_string("Topology elements are not tet.");
  }

  whyBad = info.to_summary_string();

  return rval;
}

void TetMeshClipper::transformCoordset()
{
  // Apply transformations
  auto& oldCoordset = m_tetMesh.fetch_existing("coordsets").fetch_existing(m_coordsetName);
  const std::string newCoordsetName = m_coordsetName + ".trans";
  conduit::Node& coordset = m_tetMesh.fetch("coordsets")[newCoordsetName];
  coordset.set_node(oldCoordset);
  auto transformer = m_transformer;
  conduit::index_t count = conduit::blueprint::mesh::coordset::length(coordset);
  axom::ArrayView<double> xV(coordset.fetch_existing("values/x").as_double_ptr(), count);
  axom::ArrayView<double> yV(coordset.fetch_existing("values/y").as_double_ptr(), count);
  axom::ArrayView<double> zV(coordset.fetch_existing("values/z").as_double_ptr(), count);
  const axom::IndexType numPts = static_cast<axom::IndexType>(count);
  for(axom::IndexType i = 0; i < numPts; ++i)
  {
    transformer.transform(xV[i], yV[i], zV[i]);
    m_tetMeshBb.addPoint(Point3DType {xV[i], yV[i], zV[i]});
  }
  m_tetMesh.fetch_existing("topologies")
    .fetch_existing(m_topoName)
    .fetch_existing("coordset")
    .set_string(newCoordsetName);
  m_coordsetName = newCoordsetName;
}

/*
 * Compute the surface of the tet mesh, using bump utilities.
 */
template <typename ExecSpace>
axom::Array<TetMeshClipper::Triangle3DType> TetMeshClipper::computeGeometrySurface(int allocId)
{
  // Copy some m_tetMesh data to allocId for accessing in ExecSpace.
  copy_topo_and_coords_to(allocId);

  /*
   * Make a view of the tet mesh topology.
   *
   * Note: Using axom::IndexType for integers.  Users should configure
   * Axom with the type that they plan to use.  Mixing different
   * IndexTypes in the same Axom build is not supported.
   */
  namespace bumpViews = axom::bump::views;
  namespace bumpUtils = axom::bump::utilities;

  constexpr axom::IndexType FACES_PER_TET = 4;
  constexpr axom::IndexType VERTS_PER_TET = 4;
  constexpr axom::IndexType VERTS_PER_FACE = 3;

  const std::string topoName = axom::fmt::format("{}.{}", m_topoName, allocId);
  conduit::Node& tetTopo = m_tetMesh["topologies"][topoName];
  conduit::Node& tetVertConn = tetTopo.fetch_existing("elements").fetch_existing("connectivity");
  const auto tetVertConnView = bumpUtils::make_array_view<axom::IndexType>(tetVertConn);

  conduit::Node& polyTopo = m_tetMesh["topologies"]["polyhedral"];
  make_polyhedral_topology<ExecSpace>(tetTopo, polyTopo);
  const auto& faceVertConn = polyTopo.fetch_existing("subelements/connectivity");
  const axom::IndexType faceCount = faceVertConn.dtype().number_of_elements() / VERTS_PER_FACE;
  const auto faceVertsView = bumpUtils::make_array_view<axom::IndexType>(faceVertConn);
  SLIC_ASSERT(faceVertsView.size() == faceCount * VERTS_PER_FACE);

  conduit::Node& polyElemFaceConn = polyTopo["elements/connectivity"];
  const auto polyElemFaceView = bumpUtils::make_array_view<axom::IndexType>(polyElemFaceConn);
  SLIC_ASSERT(polyElemFaceView.size() == m_tetCount * FACES_PER_TET);

  const std::string coordsName = polyTopo["coordset"].as_string();
  const conduit::Node& coordsNode = m_tetMesh.fetch_existing("coordsets").fetch_existing(coordsName);
  auto xs = bumpUtils::make_array_view<double>(coordsNode["values/x"]);
  auto ys = bumpUtils::make_array_view<double>(coordsNode["values/y"]);
  auto zs = bumpUtils::make_array_view<double>(coordsNode["values/z"]);

  /*
   * Compute tet faces as triangles.
   * Compute rays from triangle centroid, in normal direction.
   */
  axom::Array<Triangle3DType> faceTris(faceCount, faceCount, allocId);
  axom::Array<Ray3DType> faceRays(faceCount, faceCount, allocId);
  auto faceTrisView = faceTris.view();
  auto faceRaysView = faceRays.view();
  axom::for_all<ExecSpace>(
    faceCount,
    AXOM_LAMBDA(axom::IndexType faceIdx) {
      axom::IndexType* vertIds = faceVertsView.data() + faceIdx * VERTS_PER_FACE;
      auto vId0 = vertIds[0];
      auto vId1 = vertIds[1];
      auto vId2 = vertIds[2];
      Point3DType v0 {xs[vId0], ys[vId0], zs[vId0]};
      Point3DType v1 {xs[vId1], ys[vId1], zs[vId1]};
      Point3DType v2 {xs[vId2], ys[vId2], zs[vId2]};
      auto& faceTri = faceTrisView[faceIdx] = Triangle3DType(v0, v1, v2);
      faceRaysView[faceIdx] = Ray3DType(faceTri.centroid(), faceTri.normal());
    });

  /*
   * Compute whether faces have tets on each side.
   */
  axom::Array<bool> hasCellOnFrontSide(faceCount, 0, allocId);
  axom::Array<bool> hasCellOnBackSide(faceCount, 0, allocId);
  hasCellOnFrontSide.fill(false);
  hasCellOnBackSide.fill(false);
  auto hasCellOnFrontSideView = hasCellOnFrontSide.view();
  auto hasCellOnBackSideView = hasCellOnBackSide.view();
  axom::for_all<ExecSpace>(
    m_tetCount,
    AXOM_LAMBDA(IndexType tetId) {
      Point3DType cellCentroid({0.0, 0.0, 0.0});
      axom::IndexType* vertIdxs = tetVertConnView.data() + tetId * VERTS_PER_TET;
      TetrahedronType tet;
      for(int vi = 0; vi < VERTS_PER_TET; ++vi)
      {
        axom::IndexType vIdx = vertIdxs[vi];
        Point3DType vertCoords {xs[vIdx], ys[vIdx], zs[vIdx]};
        cellCentroid.array() += vertCoords.array();
        tet[vi] = vertCoords;
      }
      cellCentroid.array() /= VERTS_PER_TET;

      for(int fi = 0; fi < FACES_PER_TET; ++fi)
      {
        axom::IndexType faceIdx = polyElemFaceView[tetId * FACES_PER_TET + fi];
        const auto& faceTri = faceTrisView[faceIdx];
        Ray3DType faceRay(faceTri.centroid(), faceTri.normal());
        const Vector3DType& faceDir = faceRay.direction();
        auto cellToFace = cellCentroid.array() - faceRay.origin().array();
        double distParam =
          faceDir[0] * cellToFace[0] + faceDir[1] * cellToFace[1] + faceDir[2] * cellToFace[2];
        bool& hasCell =
          distParam >= 0 ? hasCellOnFrontSideView[faceIdx] : hasCellOnBackSideView[faceIdx];
        hasCell = true;
      }
    });

  /*
   * Mark faces touching only 1 cell.
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
   * Get running total of surface triangle count using prefix-sum scan.
   * Then use the results to populate array of those faces.
   */
  axom::Array<axom::IndexType> prefixSum(faceCount + 1, 0, allocId);
  prefixSum.fill(0);
  auto prefixSumView = prefixSum.view();
  axom::inclusive_scan<ExecSpace>(hasCellOnOneSide,
                                  axom::ArrayView<axom::IndexType>(prefixSum.data() + 1, faceCount));

  axom::IndexType surfFaceCount = -1;
  axom::copy(&surfFaceCount, prefixSumView.data() + prefixSumView.size() - 1, sizeof(surfFaceCount));

  axom::Array<axom::IndexType> surfFaceIds(surfFaceCount, 0, allocId);
  axom::Array<Triangle3DType> surfTris(surfFaceCount, 0, allocId);
  auto surfFaceIdsView = surfFaceIds.view();
  auto surfTrisView = surfTris.view();
  axom::for_all<ExecSpace>(
    faceCount,
    AXOM_LAMBDA(axom::IndexType faceIdx) {
      if(prefixSumView[faceIdx] != prefixSumView[faceIdx + 1])
      {
        auto runningTotal = prefixSumView[faceIdx];
        surfFaceIdsView[runningTotal] = faceIdx;
        surfTrisView[runningTotal] = faceTrisView[faceIdx];
      }
    });

  return surfTris;
}

void TetMeshClipper::writeTrianglesToVTK(const axom::Array<Triangle3DType>& triangles,
                                         const std::string& filename)
{
  std::ofstream ofs(filename);
  if(!ofs)
  {
    std::cerr << "Cannot open file for writing: " << filename << std::endl;
    return;
  }

  axom::Array<Triangle3DType> hostTriangles(triangles, axom::MALLOC_ALLOCATOR_ID);

  // Header
  ofs << "# vtk DataFile Version 3.0\n";
  ofs << "Triangle mesh\n";
  ofs << "ASCII\n";
  ofs << "DATASET POLYDATA\n";

  // Write points
  ofs << "POINTS " << hostTriangles.size() * 3 << " double\n";
  for(const auto& tri : hostTriangles)
  {
    for(int i = 0; i < 3; ++i)
    {
      const auto& pt = tri[i];
      ofs << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
    }
  }

  // Write polygons (triangles)
  ofs << "POLYGONS " << hostTriangles.size() << " " << hostTriangles.size() * 4 << "\n";
  for(int i = 0; i < hostTriangles.size(); ++i)
  {
    ofs << "3 " << i * 3 << " " << i * 3 + 1 << " " << i * 3 + 2 << "\n";
  }

  ofs.close();
}

template <typename ExecSpace>
void TetMeshClipper::make_polyhedral_topology(conduit::Node& tetTopo, conduit::Node& polyTopo)
{
  namespace bumpViews = axom::bump::views;
  namespace bumpUtils = axom::bump::utilities;
  using BumpShape = bumpViews::TetShape<axom::IndexType>;

  SLIC_ASSERT(tetTopo["type"].as_string() == std::string("unstructured"));
  auto tetTopoView = bumpViews::make_unstructured_single_shape_topology<BumpShape>::view(tetTopo);
  using TopologyView = decltype(tetTopoView);
  using ConnectivityType = typename TopologyView::ConnectivityType;

  bump::MakePolyhedralTopology<ExecSpace, TopologyView> polyTopoMaker(tetTopoView);
  polyTopoMaker.execute(tetTopo, polyTopo);
  bump::MergePolyhedralFaces<ExecSpace, ConnectivityType>::execute(polyTopo);
}

/*
 * Copy a conduit Node, with special allocator ids for new arrays.
 */
void TetMeshClipper::copy_node_with_array_allocator(conduit::Node& src,
                                                    conduit::Node& dst,
                                                    conduit::index_t conduitArrayAllocId)
{
  if(src.number_of_children() == 0)  // Leaf node
  {
    auto& srcDtype = src.dtype();
    if(srcDtype.number_of_elements() > 1 && !srcDtype.is_string())
    {
      dst.set_allocator(conduitArrayAllocId);
    }
    else
    {
      dst.set_allocator(src.allocator());
    }
    dst.set(src);
  }
  else  // Branch: recurse into children
  {
    for(const auto& name : src.child_names())
    {
      conduit::Node& childSrc = src[name];
      conduit::Node& childDst = dst[name];
      copy_node_with_array_allocator(childSrc, childDst, conduitArrayAllocId);
    }
  }
  return;
}

void TetMeshClipper::copy_hierarchy_with_array_allocator(conduit::Node& hierarchy,
                                                         const std::string& srcPath,
                                                         const std::string& dstPath,
                                                         int allocId)
{
  if(hierarchy.has_path(dstPath))
  {
    return;
  }  // Already done.

  conduit::Node& src = hierarchy.fetch_existing(srcPath);
  conduit::Node& dst = hierarchy[dstPath];

  const conduit::index_t conduitArrayAllocId = sidre::ConduitMemory::axomAllocIdToConduit(allocId);
  copy_node_with_array_allocator(src, dst, conduitArrayAllocId);
  return;
}

/*
 * Copy m_tetMesh topology and coordset into allocId if it's not there yet.
 * (Copy only the arrays this object uses: coordset and connectivity.)
 * We give new array data keys "<oldKey>.<allocId>".
 * If the new array exists, no need to redo it.
 */
void TetMeshClipper::copy_topo_and_coords_to(int allocId)
{
  const std::string oldTopoPath = "topologies/" + m_topoName;
  const std::string oldCoordsetPath = "coordsets/" + m_coordsetName;
  const std::string newTopoPath = axom::fmt::format("{}.{}", oldTopoPath, allocId);
  const std::string newCoordsetPath = axom::fmt::format("{}.{}", oldCoordsetPath, allocId);

  copy_hierarchy_with_array_allocator(m_tetMesh, oldTopoPath, newTopoPath, allocId);
  copy_hierarchy_with_array_allocator(m_tetMesh, oldCoordsetPath, newCoordsetPath, allocId);

  const std::string newCoordsetName = axom::fmt::format("{}.{}", m_coordsetName, allocId);
  m_tetMesh[newTopoPath + "/coordset"].set(newCoordsetName);
}

}  // namespace experimental
}  // end namespace quest
}  // end namespace axom
