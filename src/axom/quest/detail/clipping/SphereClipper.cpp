// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/quest/Discretize.hpp"
#include "axom/quest/detail/clipping/SphereClipper.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{

SphereClipper::SphereClipper(const klee::Geometry& kGeom, const std::string& name)
  : MeshClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("Sphere") : name)
  , m_transformer(m_extTrans)
{
  extractClipperInfo();

  transformSphere();
}

bool SphereClipper::labelCellsInOut(quest::experimental::ShapeMesh& shapeMesh, axom::Array<LabelType>& labels)
{
  AXOM_ANNOTATE_SCOPE("SphereClipper::labelCellsInOut");
  switch(shapeMesh.getRuntimePolicy())
  {
  case axom::runtime_policy::Policy::seq:
    labelCellsInOutImpl<axom::SEQ_EXEC>(shapeMesh, labels);
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    labelCellsInOutImpl<axom::OMP_EXEC>(shapeMesh, labels);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case axom::runtime_policy::Policy::cuda:
    labelCellsInOutImpl<axom::CUDA_EXEC<256>>(shapeMesh, labels);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    labelCellsInOutImpl<axom::HIP_EXEC<256>>(shapeMesh, labels);
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
  return true;
}

template <typename ExecSpace>
void SphereClipper::labelCellsInOutImplOld(quest::experimental::ShapeMesh& shapeMesh,
                                           axom::Array<LabelType>& labels)
{
  SLIC_ERROR_IF(shapeMesh.dimension() != 3, "SphereClipper requires a 3D mesh.");

  constexpr int NUM_VERTS_PER_CELL = 8;

  int allocId = shapeMesh.getAllocatorID();
  auto cellCount = shapeMesh.getCellCount();
  auto vertCount = shapeMesh.getVertexCount();

  const auto& vertCoords = shapeMesh.getVertexCoords3D();
  const auto& vX = vertCoords[0];
  const auto& vY = vertCoords[1];
  const auto& vZ = vertCoords[2];

  /*
    Compute whether vertices are inside shape.
  */
  axom::Array<double> vertDist {ArrayOptions::Uninitialized(), vertCount, vertCount, allocId};
  auto vertDistView = vertDist.view();
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vX.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vY.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vZ.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vertDistView.getAllocatorID()));

  auto sphere = m_sphere;
  axom::for_all<ExecSpace>(
    vertCount,
    AXOM_LAMBDA(axom::IndexType vertId) {
      primal::Point3D vert {vX[vertId], vY[vertId], vZ[vertId]};
      double signedDist = sphere.computeSignedDistance(vert);
      vertDistView[vertId] = signedDist;
    });

  if(labels.size() < cellCount || labels.getAllocatorID() != shapeMesh.getAllocatorID())
  {
    labels = axom::Array<LabelType>(ArrayOptions::Uninitialized(), cellCount, cellCount, allocId);
  }

  /*
   * Label cell by whether it has vertices inside, outside or both.
   * Sphere may intersect cell between its vertices, so we classify
   * the cells conservatively.
   * - If cell has vertex less than a small distance outside the
   *   sphere, the cell is classified as having a vertex inside.
   * - We don't need to do the same for vertices a small distance
   *   inside the sphere, because the sphere is convex.  There's no
   *   cell with all vertices inside the sphere and intersects the
   *   sphere.
  */
  auto cellLengths = shapeMesh.getCellLengths();
  const double lenFactor = 0.5;

  axom::ArrayView<const axom::IndexType, 2> connView = shapeMesh.getCellNodeConnectivity();
  SLIC_ASSERT(connView.shape() ==
              (axom::StackArray<axom::IndexType, 2> {cellCount, NUM_VERTS_PER_CELL}));

  auto labelsView = labels.view();

  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
      LabelType& cellLabel = labelsView[cellId];
      auto cellVertIds = connView[cellId];
      const double proximityThreshold = cellLengths[cellId] * lenFactor;
      bool hasIn = vertDistView[cellVertIds[0]] < proximityThreshold;
      bool hasOut = vertDistView[cellVertIds[0]] > 0;
      for(int vi = 1; vi < NUM_VERTS_PER_CELL; ++vi)
      {
        int vertId = cellVertIds[vi];
        bool isIn = vertDistView[vertId] < proximityThreshold;
        bool isOut = vertDistView[vertId] > 0;
        hasIn |= isIn;
        hasOut |= isOut;
      }
      cellLabel = !hasOut ? LabelType::LABEL_IN : !hasIn ? LabelType::LABEL_OUT : LabelType::LABEL_ON;
    });

  bool checkEdges = false;
  if(checkEdges)
  {
    /*
     * The vertices are all outside, Check if any edge crosses the
     * geometry.  This check rarely makes a difference.  Any errors
     * corrected is probably O(h^3) and much smaller than the error
     * from discretizing the sphere.  It's probably not worth the cost.
    */
    axom::ArrayView<const TetrahedronType> cellsAsTets = shapeMesh.getCellsAsTets();
    axom::for_all<ExecSpace>(
      cellCount,
      AXOM_LAMBDA(axom::IndexType cellId) {
        LabelType& cellLabel = labelsView[cellId];
        if(cellLabel == LabelType::LABEL_OUT)
        {
          constexpr int NUM_TETS_PER_HEX = primal::Hexahedron<double, 3>::NUM_TRIANGULATE;
          const double sqRadius = sphere.getRadius() * sphere.getRadius();

          const axom::IndexType tetIdxStart = cellId * NUM_TETS_PER_HEX;
          const axom::IndexType tetIdxEnd = (1 + cellId) * NUM_TETS_PER_HEX;
          for(axom::IndexType ti = tetIdxStart; ti < tetIdxEnd; ++ti)
          {
            const TetrahedronType& tet = cellsAsTets[ti];
            for(int vA = 0; vA < 4 && cellLabel == LabelType::LABEL_OUT; ++vA)
            {
              for(int vB = vA + 1; vB < 4 && cellLabel == LabelType::LABEL_OUT; ++vB)
              {
                const Segment3DType seg(tet[vA], tet[vB]);
                const Vector3DType vec(tet[vA], tet[vB]);
                const Plane3DType plane(vec, sphere.getCenter());
                double t;
                bool intersects = axom::primal::intersect(plane, seg, t);
                if(intersects)
                {
                  Point3DType intersectionPt = seg.at(t);
                  Vector3DType centerToIntersection(sphere.getCenter(), intersectionPt);
                  double sqNorm = centerToIntersection.squared_norm();
                  if(sqNorm < sqRadius)
                  {
                    cellLabel = LabelType::LABEL_ON;
                  }
                }
              }
            }
          }
        }
      });
  }

  return;
}

template <typename ExecSpace>
void SphereClipper::labelCellsInOutImpl(quest::experimental::ShapeMesh& shapeMesh, axom::Array<LabelType>& labels)
{
  SLIC_ERROR_IF(shapeMesh.dimension() != 3, "SphereClipper requires a 3D mesh.");

  constexpr int NUM_VERTS_PER_CELL = 8;

  int allocId = shapeMesh.getAllocatorID();
  auto cellCount = shapeMesh.getCellCount();
  auto vertCount = shapeMesh.getVertexCount();

  auto sphere = m_sphere;

  const auto& vertCoords = shapeMesh.getVertexCoords3D();
  const auto& vX = vertCoords[0];
  const auto& vY = vertCoords[1];
  const auto& vZ = vertCoords[2];

  /*
   * Compute whether vertices are inside shape.
  */

  axom::Array<bool> vertIsOutside {ArrayOptions::Uninitialized(), vertCount, vertCount, allocId};
  auto vertIsOutsideView = vertIsOutside.view();
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vX.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vY.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vZ.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vertIsOutsideView.getAllocatorID()));

  axom::for_all<ExecSpace>(
    vertCount,
    AXOM_LAMBDA(axom::IndexType vertId) {
      primal::Point3D vert {vX[vertId], vY[vertId], vZ[vertId]};
      int orientation = sphere.getOrientation(vert, 0.0);
      vertIsOutsideView[vertId] = orientation == primal::ON_POSITIVE_SIDE;
    });

  /*
   * Compute cell labels.
  */

  axom::ArrayView<const axom::IndexType, 2> connView = shapeMesh.getCellNodeConnectivity();
  SLIC_ASSERT(connView.shape() ==
              (axom::StackArray<axom::IndexType, 2> {cellCount, NUM_VERTS_PER_CELL}));

  if(labels.size() < cellCount || labels.getAllocatorID() != shapeMesh.getAllocatorID())
  {
    labels = axom::Array<LabelType>(ArrayOptions::Uninitialized(), cellCount, cellCount, allocId);
  }

  auto labelsView = labels.view();

  /*
   * Label cell:
   * - Compute cell's bounding sphere.  Use bounding box's
   *   bounding sphere as a fast conservative approximation.
   * - If bounding sphere doesn't intersect the geometry, cell is LabelType::LABEL_OUT.
   * - If all cell vertices are inside geometry, cell is LabelType::LABEL_IN.
   *   This is true because geometry is convex.
   * - If spheres intersect and not all vertices are inside,
   *   some parts of the cell may intersect boundary, so it's LabelType::LABEL_ON.
  */

  auto cellBbs = shapeMesh.getCellBoundingBoxes();

  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
      LabelType& cellLabel = labelsView[cellId];
      const auto& bb = cellBbs[cellId];
      const SphereType boundingSphere(bb.getCentroid(), bb.range().norm() / 2);
      if(sphere.intersectsWith(boundingSphere, 0.0))
      {
        auto cellVertIds = connView[cellId];
        bool hasOut = vertIsOutsideView[cellVertIds[0]];
        for(int vi = 1; vi < NUM_VERTS_PER_CELL; ++vi)
        {
          int vertId = cellVertIds[vi];
          hasOut |= vertIsOutsideView[vertId];
        }
        cellLabel = hasOut ? LabelType::LABEL_ON : LabelType::LABEL_IN;
      }
      else
      {
        cellLabel = LabelType::LABEL_OUT;
      }
    });

  return;
}

bool SphereClipper::labelTetsInOut(quest::experimental::ShapeMesh& shapeMesh,
                                   axom::ArrayView<const axom::IndexType> cellIds,
                                   axom::Array<LabelType>& tetLabels)
{
  AXOM_ANNOTATE_SCOPE("SphereClipper::labelTetsInOut");
  switch(shapeMesh.getRuntimePolicy())
  {
  case axom::runtime_policy::Policy::seq:
    labelTetsInOutImpl<axom::SEQ_EXEC>(shapeMesh, cellIds, tetLabels);
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    labelTetsInOutImpl<axom::OMP_EXEC>(shapeMesh, cellIds, tetLabels);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case axom::runtime_policy::Policy::cuda:
    labelTetsInOutImpl<axom::CUDA_EXEC<256>>(shapeMesh, cellIds, tetLabels);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    labelTetsInOutImpl<axom::HIP_EXEC<256>>(shapeMesh, cellIds, tetLabels);
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
  return true;
}

template <typename ExecSpace>
void SphereClipper::labelTetsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                                       axom::ArrayView<const axom::IndexType> cellIds,
                                       axom::Array<LabelType>& tetLabels)
{
  SLIC_ERROR_IF(shapeMesh.dimension() != 3, "SphereClipper requires a 3D mesh.");

  constexpr int NUM_VERTS_PER_TET = 4;

  const axom::IndexType cellCount = cellIds.size();

  int allocId = shapeMesh.getAllocatorID();

  if(tetLabels.size() < cellCount * NUM_TETS_PER_HEX ||
     tetLabels.getAllocatorID() != shapeMesh.getAllocatorID())
  {
    tetLabels = axom::Array<LabelType>(ArrayOptions::Uninitialized(),
                                       cellCount * NUM_TETS_PER_HEX,
                                       cellCount * NUM_TETS_PER_HEX,
                                       allocId);
  }
  auto tetLabelsView = tetLabels.view();

  /*
   * Label tets:
   * - Compute tets's bounding sphere.  Use bounding box's
   *   bounding sphere as a fast conservative approximation.
   * - If bounding sphere doesn't intersect sphere geometry,
   *   cell is LabelType::LABEL_OUT.
   * - If all cell vertices are inside sphere geometry,
   *   cell is LabelType::LABEL_IN.  This is true because geometry is convex.
   * - If spheres intersect and not all vertices are inside,
   *   some parts of the cell may intersect boundary, so it's LabelType::LABEL_ON.
   *
   * TODO: It shouldn't be hard to compute the closest point
   * of a BoundingBox to another point.  It would be less
   * conservative than using the BoundingBox circumsphere.
   * Worth doing because too many tets are unnecessarily
   * labeled ON.
  */

  auto meshHexes = shapeMesh.getCellsAsHexes();

  auto sphere = m_sphere;
  const double squaredRad = sphere.getRadius() * sphere.getRadius();

  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType ci) {
      axom::IndexType cellId = cellIds[ci];
      const HexahedronType& hex = meshHexes[cellId];

      TetrahedronType cellTets[NUM_TETS_PER_HEX];
      hex.triangulate(cellTets);

      for(IndexType ti = 0; ti < NUM_TETS_PER_HEX; ++ti)
      {
        const TetrahedronType& tet = cellTets[ti];
        LabelType& tetLabel = tetLabelsView[ci * NUM_TETS_PER_HEX + ti];

        BoundingBox3DType bb(tet[0]);
        bb.addPoint(tet[1]);
        bb.addPoint(tet[2]);
        bb.addPoint(tet[3]);
        const SphereType boundingSphere(bb.getCentroid(), bb.range().norm() / 2);

        if(sphere.intersectsWith(boundingSphere, 0.0))
        {
          bool hasOut = false;
          for(int vi = 0; vi < NUM_VERTS_PER_TET; ++vi)
          {
            double sqDistance = axom::primal::squared_distance(sphere.getCenter(), tet[vi]);
            hasOut |= sqDistance > squaredRad;
          }
          tetLabel = hasOut ? LabelType::LABEL_ON : LabelType::LABEL_IN;
        }
        else
        {
          tetLabel = LabelType::LABEL_OUT;
        }
      }
    });

  return;
}

/*
  TODO: If possible: Port to GPU.  Will need to rewrite quest/Discretize.[ch]pp.
*/
bool SphereClipper::getGeometryAsOcts(quest::experimental::ShapeMesh& shapeMesh,
                                      axom::Array<axom::primal::Octahedron<double, 3>>& octs)
{
  AXOM_ANNOTATE_SCOPE("SphereClipper::getGeometryAsOcts");
  int octCount = 0;
  axom::quest::discretize(m_sphereBeforeTrans, m_levelOfRefinement, octs, octCount);

  auto octsView = octs.view();
  auto transformer = m_transformer;
  int allocId = shapeMesh.getAllocatorID();
  axom::for_all<axom::SEQ_EXEC>(
    octCount,
    AXOM_LAMBDA(axom::IndexType iOct) {
      OctahedronType& oct = octsView[iOct];
      for(int iVert = 0; iVert < OctType::NUM_VERTS; ++iVert)
      {
        Point3DType& ptCoords = oct[iVert];
        transformer.transform(ptCoords.array());
      }
    });

  // The disretize method uses host data.  Place into proper space if needed.
  if(octs.getAllocatorID() != allocId)
  {
    octs = axom::Array<axom::primal::Octahedron<double, 3>>(octs, allocId);
  }
  return true;
}

void SphereClipper::extractClipperInfo()
{
  const auto c = m_info.fetch_existing("center").as_double_array();
  const double radius = m_info.fetch_existing("radius").as_double();
  Point3DType center;
  for(int d = 0; d < 3; ++d)
  {
    center[d] = c[d];
  }
  m_sphereBeforeTrans = SphereType(center, radius);
  m_levelOfRefinement = m_info.fetch_existing("levelOfRefinement").to_int32();
}

// Include external transformations in m_sphere.
void SphereClipper::transformSphere()
{
  const auto& centerBeforeTrans = m_sphereBeforeTrans.getCenter();
  const double radiusBeforeTrans = m_sphereBeforeTrans.getRadius();
  Point3DType surfacePtBeforeTrans {centerBeforeTrans.array() +
                                    Point3DType::NumericArray {radiusBeforeTrans, 0, 0}};

  auto center = m_transformer.getTransformed(centerBeforeTrans);
  Point3DType surfacePoint = m_transformer.getTransformed(surfacePtBeforeTrans);
  const double radius = Vector3DType(center, surfacePoint).norm();
  m_sphere = SphereType(center, radius);
}

}  // namespace experimental
}  // end namespace quest
}  // end namespace axom
