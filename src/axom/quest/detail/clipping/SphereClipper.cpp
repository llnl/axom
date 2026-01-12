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

bool SphereClipper::labelCellsInOut(quest::experimental::ShapeMesh& shapeMesh,
                                    axom::Array<LabelType>& labels)
{
  SLIC_ERROR_IF(shapeMesh.dimension() != 3, "SphereClipper requires a 3D mesh.");

  int allocId = shapeMesh.getAllocatorID();
  auto cellCount = shapeMesh.getCellCount();
  if(labels.size() < cellCount || labels.getAllocatorID() != shapeMesh.getAllocatorID())
  {
    labels = axom::Array<LabelType>(ArrayOptions::Uninitialized(), cellCount, cellCount, allocId);
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

template <typename ExecSpace>
void SphereClipper::labelCellsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                                        axom::ArrayView<LabelType> labels)
{
  auto cellCount = shapeMesh.getCellCount();
  auto cellsAsHexes = shapeMesh.getCellsAsHexes();
  auto cellVolumes = shapeMesh.getCellVolumes();
  constexpr double EPS = 1e-10;
  auto sphere = m_sphere;
  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
      LabelType& cellLabel = labels[cellId];
      if(axom::utilities::isNearlyEqual(cellVolumes[cellId], 0.0, EPS))
      {
        cellLabel = LabelType::LABEL_OUT;
        return;
      }
      const auto& hex = cellsAsHexes[cellId];
      cellLabel = polyhedronToLabel(hex, sphere);
      // Note: cellLabel may be set to LABEL_ON if polyhedronToLabel
      // cannot efficiently determine whether the hex is IN or OUT.
      // See MeshClipperStrategy::labelCellsInOut().
    });
  return;
}

bool SphereClipper::labelTetsInOut(quest::experimental::ShapeMesh& shapeMesh,
                                   axom::ArrayView<const axom::IndexType> cellIds,
                                   axom::Array<LabelType>& tetLabels)
{
  SLIC_ERROR_IF(shapeMesh.dimension() != 3, "SphereClipper requires a 3D mesh.");

  const axom::IndexType cellCount = cellIds.size();
  const int allocId = shapeMesh.getAllocatorID();

  if(tetLabels.size() < cellCount * NUM_TETS_PER_HEX ||
     tetLabels.getAllocatorID() != shapeMesh.getAllocatorID())
  {
    tetLabels = axom::Array<LabelType>(ArrayOptions::Uninitialized(),
                                       cellCount * NUM_TETS_PER_HEX,
                                       cellCount * NUM_TETS_PER_HEX,
                                       allocId);
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

template <typename ExecSpace>
void SphereClipper::labelTetsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                                       axom::ArrayView<const axom::IndexType> cellIds,
                                       axom::ArrayView<LabelType> tetLabels)
{
  const axom::IndexType cellCount = cellIds.size();
  auto meshHexes = shapeMesh.getCellsAsHexes();
  auto tetVolumes = shapeMesh.getTetVolumes();
  constexpr double EPS = 1e-10;
  auto sphere = m_sphere;

  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType ci) {
      axom::IndexType cellId = cellIds[ci];
      const HexahedronType& hex = meshHexes[cellId];

      TetrahedronType cellTets[NUM_TETS_PER_HEX];
      ShapeMesh::hexToTets(hex, cellTets);

      for(IndexType ti = 0; ti < NUM_TETS_PER_HEX; ++ti)
      {
        LabelType& tetLabel = tetLabels[ci * NUM_TETS_PER_HEX + ti];
        const axom::IndexType tetId = cellId * NUM_TETS_PER_HEX + ti;
        if(axom::utilities::isNearlyEqual(tetVolumes[tetId], 0.0, EPS))
        {
          tetLabel = LabelType::LABEL_OUT;
          continue;
        }
        const TetrahedronType& tet = cellTets[ti];
        tetLabel = polyhedronToLabel(tet, sphere);
        // Note: cellLabel may be set to LABEL_ON if polyhedronToLabel
        // cannot efficiently determine whether the tet is IN or OUT.
        // See MeshClipperStrategy::labelTetsInOut().
      }
    });
  return;
}

template <typename Polyhedron>
AXOM_HOST_DEVICE inline MeshClipperStrategy::LabelType SphereClipper::polyhedronToLabel(
  const Polyhedron& verts,
  const SphereType& sphere) const
{
  /*
    If bounding box of polyhedron is more than the radius distance
    from center, it is LABEL_OUT.  (Comparing vertices for this check
    can miss intersections by edges and facets, so we compare bounding
    box.)

    Otherwise, polyhedron is labeled either LABEL_ON or LABEL_IN.
    Sphere is convex, so polyhedron is IN only if all vertices are inside.

    Some polyhedra may be LABEL_ON even though they are actually LABEL_OUT,
    but this is a conservative error.  The clip function will compute the
    correct overlap volume.  The purpose of labeling is bypass the
    clip function where we can do it efficiently.
  */
  BoundingBox3DType bb(verts[0]);
  auto vertCount = Polyhedron::numVertices();
  for(int i = 1; i < vertCount; ++i)
  {
    bb.addPoint(verts[i]);
  }

  const double sqRad = sphere.getRadius() * sphere.getRadius();

  double sqDistToBb = primal::squared_distance(sphere.getCenter(), bb);

  if(sqDistToBb >= sqRad)
  {
    return LabelType::LABEL_OUT;
  }

  for(int i = 0; i < vertCount; ++i)
  {
    const auto& vert = verts[i];
    double sqDistToVert = axom::primal::squared_distance(sphere.getCenter(), vert);
    if(sqDistToVert > sqRad)
    {
      return LabelType::LABEL_ON;
    }
  }
  return LabelType::LABEL_IN;
}

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

  SLIC_INFO(axom::fmt::format("SphereClipper '{}' {}-level refined got {} geometry octs.",
                              name(),
                              m_levelOfRefinement,
                              octs.size()));
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
