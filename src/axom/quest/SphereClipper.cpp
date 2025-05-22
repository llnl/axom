// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifdef AXOM_USE_CONDUIT
  #include "conduit_blueprint.hpp"
#endif

#include "axom/quest/Discretize.hpp"
#include "axom/quest/SphereClipper.hpp"
#include "axom/fmt.hpp"

namespace axom
{
namespace quest
{

SphereClipper::SphereClipper(const klee::Geometry& kGeom, const std::string& name)
  : GeometryClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("Sphere") : name)
  , m_transformer(m_transMat)
{
  extractClipperInfo();

  transformSphere();
}

bool SphereClipper::labelInOut(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  AXOM_ANNOTATE_BEGIN("SphereClipper::labelInOut");
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
  AXOM_ANNOTATE_END("SphereClipper::labelInOut");
  return true;
}

template <typename ExecSpace>
void SphereClipper::labelInOutImpl(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  SLIC_ERROR_IF(shapeeMesh.dimension() != 3, "SphereClipper requires a 3D mesh.");

  constexpr int NUM_VERTS_PER_CELL = 8;

  int allocId = shapeeMesh.getAllocatorID();
  auto cellCount = shapeeMesh.getCellCount();
  auto vertCount = shapeeMesh.getVertexCount();

  const auto& vertCoords = shapeeMesh.getVertexCoords3D();
  const auto& vX = vertCoords[0];
  const auto& vY = vertCoords[1];
  const auto& vZ = vertCoords[2];

  /*
    Compute whether vertices are inside shape.
  */
  axom::Array<bool> vertIsInside {ArrayOptions::Uninitialized(), vertCount, vertCount, allocId};
  auto vertIsInsideView = vertIsInside.view();
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vX.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vY.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vZ.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vertIsInsideView.getAllocatorID()));

  auto sphere = m_sphere;
  axom::for_all<ExecSpace>(
    vertCount,
    AXOM_LAMBDA(axom::IndexType vertId) {
      primal::Point3D vert {vX[vertId], vY[vertId], vZ[vertId]};
      double signedDist = sphere.computeSignedDistance(vert);
      vertIsInsideView[vertId] = signedDist < 0;
    });

  if(labels.size() < cellCount || labels.getAllocatorID() != shapeeMesh.getAllocatorID())
  {
    labels = axom::Array<LabelType>(ArrayOptions::Uninitialized(), cellCount, cellCount, allocId);
  }

  /*
    Label cell by whether it has vertices inside, outside or both.
  */
  axom::ArrayView<const axom::IndexType, 2> connView = shapeeMesh.getConnectivity();
  SLIC_ASSERT(connView.shape() ==
              (axom::StackArray<axom::IndexType, 2> {cellCount, NUM_VERTS_PER_CELL}));

  auto labelsView = labels.view();

  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
      LabelType& cellLabel = labelsView[cellId];
      auto cellVertIds = connView[cellId];
      bool hasIn = vertIsInsideView[cellVertIds[0]];
      bool hasOut = !hasIn;
      for(int vi = 0; vi < NUM_VERTS_PER_CELL; ++vi)
      {
        int vertId = cellVertIds[vi];
        bool isIn = vertIsInsideView[vertId];
        hasIn |= isIn;
        hasOut |= !isIn;
      }
      cellLabel = !hasOut ? LABEL_IN : !hasIn ? LABEL_OUT : LABEL_ON;
    });

  bool checkEdges = false;
  if(checkEdges)
  {
    /*
      The vertices are all outside, Check if any edge crosses the
      geometry.  This check rarely makes a difference.  Any errors
      corrected is probably O(h^3) and much smaller than the error
      from discretizing the sphere.  It's probably not worth the cost.
    */
    axom::ArrayView<const TetrahedronType> cellsAsTets = shapeeMesh.getCellsAsTets();
    axom::for_all<ExecSpace>(
      cellCount,
      AXOM_LAMBDA(axom::IndexType cellId) {
        LabelType& cellLabel = labelsView[cellId];
        if (cellLabel == LABEL_OUT)
        {
          constexpr int NUM_TETS_PER_HEX = primal::Hexahedron<double, 3>::NUM_TRIANGULATE;
          const double sqRadius = sphere.getRadius() * sphere.getRadius();

          const axom::IndexType tetIdxStart = cellId * NUM_TETS_PER_HEX;
          const axom::IndexType tetIdxEnd = (1 + cellId) * NUM_TETS_PER_HEX;
          for(axom::IndexType ti = tetIdxStart; ti < tetIdxEnd; ++ti)
          {
            const TetrahedronType& tet = cellsAsTets[ti];
            for(int vA = 0; vA < 4 && cellLabel == LABEL_OUT; ++vA)
            {
              for(int vB=vA + 1; vB < 4 && cellLabel == LABEL_OUT; ++vB)
              {
                const Segment3DType seg(tet[vA], tet[vB]);
                const Vector3DType vec(tet[vA], tet[vB]);
                const Plane3DType plane(vec, sphere.getCenter());
                double t;
                bool intersects = axom::primal::intersect(plane, seg, t);
                if (intersects)
                {
                  Point3DType intersectionPt = seg.at(t);
                  Vector3DType centerToIntersection(sphere.getCenter(), intersectionPt);
                  double sqNorm = centerToIntersection.squared_norm();
                  if (sqNorm < sqRadius)
                  {
                    cellLabel = LABEL_ON;
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

bool SphereClipper::getGeometryAsOcts(quest::ShapeeMesh& shapeeMesh,
                                   axom::Array<axom::primal::Octahedron<double, 3>>& octs)
{
  AXOM_ANNOTATE_BEGIN("SphereClipper::getGeometryAsOcts");
  int octCount = 0;
  axom::quest::discretize(m_sphereBeforeTrans, m_levelOfRefinement, octs, octCount);

  auto octsView = octs.view();
  auto transformer = m_transformer;
  int allocId = shapeeMesh.getAllocatorID();
  axom::for_all<axom::SEQ_EXEC>(
    octCount,
    AXOM_LAMBDA(axom::IndexType iOct) {
      OctahedronType& oct = octsView[iOct];
      for(int iVert = 0; iVert < OctType::NUM_VERTS; ++iVert)
      {
        Point3DType& ptCoords = oct[iVert];
        transformer.transform(ptCoords);
      }
    });

  // The disretize method uses host data.  Place into proper space if needed.
  if(octs.getAllocatorID() != allocId)
  {
    octs = axom::Array<axom::primal::Octahedron<double, 3>>(octs, allocId);
  }
  AXOM_ANNOTATE_END("SphereClipper::getGeometryAsOcts");
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

void SphereClipper::transformSphere()
{
  const auto& centerBeforeTrans = m_sphereBeforeTrans.getCenter();
  const double radiusBeforeTrans = m_sphereBeforeTrans.getRadius();
  Point3DType surfacePtBeforeTrans { centerBeforeTrans.array() +
                                     Point3DType::NumericArray{radiusBeforeTrans, 0, 0} };

  auto center = m_transformer.getTransform(centerBeforeTrans);
  Point3DType surfacePoint = m_transformer.getTransform(surfacePtBeforeTrans);
  const double radius = Vector3DType(center, surfacePoint).norm();
  m_sphere = SphereType(center, radius);
}

}  // end namespace quest
}  // end namespace axom
