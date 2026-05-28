// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
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

namespace
{

constexpr double SPHERE_CLIP_EPS = 1e-10;
constexpr double SPHERE_CLIP_GRAD_EPS = 1e-20;
constexpr int SPHERE_MAX_SUBDIVISION_LEVELS = 4;
constexpr double SPHERE_LINEARIZATION_EDGE_FRACTION = 0.02;

/*!
 * @brief Construct a tetrahedron and force a positive orientation.
 *
 * Primal's clipping code expects consistently oriented tets. The midpoint
 * subdivision below creates new tets algebraically, so we normalize the
 * orientation once as each child is formed.
 */
AXOM_HOST_DEVICE
MeshClipperStrategy::TetrahedronType makePositiveTet(const MeshClipperStrategy::Point3DType& v0,
                                                     const MeshClipperStrategy::Point3DType& v1,
                                                     const MeshClipperStrategy::Point3DType& v2,
                                                     const MeshClipperStrategy::Point3DType& v3)
{
  using TetrahedronType = MeshClipperStrategy::TetrahedronType;

  TetrahedronType tet(v0, v1, v2, v3);
  tet.checkAndFixOrientation();
  return tet;
}

/*!
 * @brief Split a tetrahedron into eight children using edge midpoints.
 *
 * This reduces the local tet diameter quickly, which gives the final
 * linearized sphere clip a much better approximation than a centroid-based
 * split for the same depth.
 */
AXOM_HOST_DEVICE
void subdivideTetByMidpoints(const MeshClipperStrategy::TetrahedronType& tet,
                             MeshClipperStrategy::TetrahedronType children[8])
{
  const auto m01 = MeshClipperStrategy::Point3DType::midpoint(tet[0], tet[1]);
  const auto m02 = MeshClipperStrategy::Point3DType::midpoint(tet[0], tet[2]);
  const auto m03 = MeshClipperStrategy::Point3DType::midpoint(tet[0], tet[3]);
  const auto m12 = MeshClipperStrategy::Point3DType::midpoint(tet[1], tet[2]);
  const auto m13 = MeshClipperStrategy::Point3DType::midpoint(tet[1], tet[3]);
  const auto m23 = MeshClipperStrategy::Point3DType::midpoint(tet[2], tet[3]);

  children[0] = makePositiveTet(tet[0], m01, m02, m03);
  children[1] = makePositiveTet(m01, tet[1], m12, m13);
  children[2] = makePositiveTet(m02, m12, tet[2], m23);
  children[3] = makePositiveTet(m03, m13, m23, tet[3]);

  // Split the central octahedron along the m01-m23 diagonal.
  children[4] = makePositiveTet(m01, m02, m03, m23);
  children[5] = makePositiveTet(m01, m02, m12, m23);
  children[6] = makePositiveTet(m01, m12, m13, m23);
  children[7] = makePositiveTet(m01, m03, m13, m23);
}

/*!
 * @brief Return the squared length of the longest edge of a tetrahedron.
 *
 * The midpoint subdivision used below halves all edge lengths at each level,
 * so the longest-edge metric lets us choose how many levels are needed before
 * the final linearized clip is sufficiently local relative to sphere radius.
 */
AXOM_HOST_DEVICE
double tetMaxEdgeSquared(const MeshClipperStrategy::TetrahedronType& tet)
{
  double maxEdgeSq = 0.0;
  for(int i = 0; i < MeshClipperStrategy::TetrahedronType::NUM_VERTS; ++i)
  {
    for(int j = i + 1; j < MeshClipperStrategy::TetrahedronType::NUM_VERTS; ++j)
    {
      const double edgeSq = (tet[j] - tet[i]).squared_norm();
      maxEdgeSq = axom::utilities::max(maxEdgeSq, edgeSq);
    }
  }
  return maxEdgeSq;
}

/*!
 * @brief Choose the subdivision depth needed before local plane clipping.
 *
 * Coarse meshes need more refinement because the sphere curvature is visible
 * across a larger tet. Finer meshes can stop earlier and avoid unnecessary
 * subdivision work while still using the same GPU-friendly explicit stack.
 */
AXOM_HOST_DEVICE
int chooseSubdivisionLevels(const MeshClipperStrategy::TetrahedronType& tet,
                            const MeshClipperStrategy::SphereType& sphere)
{
  const double targetEdge = sphere.getRadius() * SPHERE_LINEARIZATION_EDGE_FRACTION;
  const double targetEdgeSq = targetEdge * targetEdge;

  double edgeSq = tetMaxEdgeSquared(tet);
  int levels = 0;
  while(levels < SPHERE_MAX_SUBDIVISION_LEVELS && edgeSq > targetEdgeSq)
  {
    edgeSq *= 0.25;
    ++levels;
  }

  return levels;
}

/*!
 * @brief Approximate the sphere/tet overlap by linearizing the signed-distance
 * field over a small tetrahedron and clipping against the resulting plane.
 *
 * The signed-distance samples at the tet vertices define an affine field whose
 * zero isosurface is the best local planar approximation available from those
 * samples. The plane is oriented so that its positive side matches the sphere
 * interior, which matches primal::clip(tet, plane).
 */
AXOM_HOST_DEVICE
double clipTetAgainstLinearizedSphere(const MeshClipperStrategy::TetrahedronType& tet,
                                      const MeshClipperStrategy::SphereType& sphere)
{
  using Plane3DType = MeshClipperStrategy::Plane3DType;
  using Point3DType = MeshClipperStrategy::Point3DType;
  using Vector3DType = MeshClipperStrategy::Vector3DType;

  const Point3DType& v0 = tet[0];
  const Point3DType& v1 = tet[1];
  const Point3DType& v2 = tet[2];
  const Point3DType& v3 = tet[3];

  const double phi0 = sphere.computeSignedDistance(v0);
  const double phi1 = sphere.computeSignedDistance(v1);
  const double phi2 = sphere.computeSignedDistance(v2);
  const double phi3 = sphere.computeSignedDistance(v3);

  const Vector3DType e1 = v1 - v0;
  const Vector3DType e2 = v2 - v0;
  const Vector3DType e3 = v3 - v0;

  const double denom = Vector3DType::scalar_triple_product(e1, e2, e3);
  if(axom::utilities::isNearlyEqual(denom, 0.0, SPHERE_CLIP_EPS))
  {
    return 0.0;
  }

  const Vector3DType grad = ((phi1 - phi0) * Vector3DType::cross_product(e2, e3) +
                             (phi2 - phi0) * Vector3DType::cross_product(e3, e1) +
                             (phi3 - phi0) * Vector3DType::cross_product(e1, e2)) /
    denom;
  const double gradSqNorm = grad.squared_norm();
  if(gradSqNorm <= SPHERE_CLIP_GRAD_EPS)
  {
    return 0.0;
  }

  const Point3DType planePoint = v0 - (phi0 / gradSqNorm) * grad;
  const Plane3DType plane(-grad, planePoint);
  const primal::Polyhedron<double, 3> overlap = primal::clip(tet, plane, SPHERE_CLIP_EPS);
  return overlap.volume();
}

/*!
 * @brief Conservatively classify a tetrahedron against the sphere.
 *
 * This mirrors SphereClipper::polyhedronToLabel() but is specialized to tets
 * so the file-local adaptive clipper can use the same inexpensive screening.
 */
AXOM_HOST_DEVICE
MeshClipperStrategy::LabelType tetToSphereLabel(const MeshClipperStrategy::TetrahedronType& tet,
                                                const MeshClipperStrategy::SphereType& sphere)
{
  using LabelType = MeshClipperStrategy::LabelType;

  MeshClipperStrategy::BoundingBox3DType bb(tet[0]);
  for(int i = 1; i < MeshClipperStrategy::TetrahedronType::NUM_VERTS; ++i)
  {
    bb.addPoint(tet[i]);
  }

  const double sqRad = sphere.getRadius() * sphere.getRadius();
  if(primal::squared_distance(sphere.getCenter(), bb) >= sqRad)
  {
    return LabelType::LABEL_OUT;
  }

  for(int i = 0; i < MeshClipperStrategy::TetrahedronType::NUM_VERTS; ++i)
  {
    if(axom::primal::squared_distance(sphere.getCenter(), tet[i]) > sqRad)
    {
      return LabelType::LABEL_ON;
    }
  }

  return LabelType::LABEL_IN;
}

/*!
 * @brief Compute a sphere/tet overlap using fixed-depth adaptive subdivision.
 *
 * Most boundary tets quickly break into fully interior or exterior subtets.
 * Only the smallest unresolved subtets pay for the local plane clip, which is
 * why this path is much cheaper than sending every boundary tet through the
 * generic sphere discretization and BVH clip pipeline.
 */
AXOM_HOST_DEVICE
double clipTetAgainstSphere(const MeshClipperStrategy::TetrahedronType& tet,
                            const MeshClipperStrategy::SphereType& sphere)
{
  using LabelType = MeshClipperStrategy::LabelType;
  using TetrahedronType = MeshClipperStrategy::TetrahedronType;
  constexpr int CHILD_COUNT = 8;
  constexpr int MAX_STACK_SIZE = 1 + (CHILD_COUNT - 1) * SPHERE_MAX_SUBDIVISION_LEVELS;

  struct WorkItem
  {
    TetrahedronType tet;
    int level;
  };

  const int maxSubdivisionLevels = chooseSubdivisionLevels(tet, sphere);
  double volume = 0.0;
  WorkItem stack[MAX_STACK_SIZE];
  int stackSize = 1;
  stack[0] = {tet, 0};

  // Use an explicit stack instead of recursion so the same code remains viable
  // in device execution spaces.
  while(stackSize > 0)
  {
    const WorkItem item = stack[--stackSize];
    const LabelType label = tetToSphereLabel(item.tet, sphere);
    if(label == LabelType::LABEL_IN)
    {
      volume += item.tet.volume();
      continue;
    }
    if(label == LabelType::LABEL_OUT)
    {
      continue;
    }

    if(item.level >= maxSubdivisionLevels)
    {
      volume += clipTetAgainstLinearizedSphere(item.tet, sphere);
      continue;
    }

    TetrahedronType children[CHILD_COUNT];
    subdivideTetByMidpoints(item.tet, children);
    for(int i = 0; i < CHILD_COUNT; ++i)
    {
      // The stack bound is fixed by the subdivision depth and fanout above.
      stack[stackSize++] = {children[i], item.level + 1};
    }
  }

  return volume;
}

}  // namespace

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
  const SphereType& sphere)
{
  /*
    If bounding box of polyhedron is more than the radius distance
    from center, it is LABEL_OUT.  (Comparing vertices for this check
    can miss intersections by edges and facets, so we compare bounding
    box.)

    Otherwise, polyhedron is labeled either LABEL_ON or LABEL_IN.
    Sphere is convex, so polyhedron is IN only if all vertices are inside.

    Some polyhedra may be labeled ON even though they are actually outside.
    This is a conservative error that is fixed ty the clip function later.
    The purpose of labeling is bypass the clip function where we can do
    it efficiently.
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
  const int allocId = shapeMesh.getAllocatorID();
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

  // The discretize method uses host data. Place into the required allocator if needed.
  if(octs.getAllocatorID() != allocId)
  {
    octs = axom::Array<axom::primal::Octahedron<double, 3>>(octs, allocId);
  }

  SLIC_DEBUG(axom::fmt::format("SphereClipper '{}' {}-level refined got {} geometry octs.",
                               name(),
                               m_levelOfRefinement,
                               octs.size()));
  return true;
}

bool SphereClipper::specializedClipTets(quest::experimental::ShapeMesh& shapeMesh,
                                        axom::ArrayView<double> ovlap,
                                        const axom::ArrayView<IndexType>& tetIds,
                                        conduit::Node& statistics)
{
  switch(shapeMesh.getRuntimePolicy())
  {
  case axom::runtime_policy::Policy::seq:
    specializedClipTetsImpl<axom::SEQ_EXEC>(shapeMesh, ovlap, tetIds, statistics);
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    specializedClipTetsImpl<axom::OMP_EXEC>(shapeMesh, ovlap, tetIds, statistics);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case axom::runtime_policy::Policy::cuda:
    specializedClipTetsImpl<axom::CUDA_EXEC<256>>(shapeMesh, ovlap, tetIds, statistics);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    specializedClipTetsImpl<axom::HIP_EXEC<256>>(shapeMesh, ovlap, tetIds, statistics);
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
  return true;
}

template <typename ExecSpace>
void SphereClipper::specializedClipTetsImpl(quest::experimental::ShapeMesh& shapeMesh,
                                            axom::ArrayView<double> ovlap,
                                            const axom::ArrayView<IndexType>& tetIds,
                                            conduit::Node& statistics)
{
  auto meshTets = shapeMesh.getCellsAsTets();
  auto sphere = m_sphere;
  const IndexType tetCount = tetIds.size();

  axom::for_all<ExecSpace>(
    tetCount,
    AXOM_LAMBDA(axom::IndexType ti) {
      const axom::IndexType tetId = tetIds[ti];
      const axom::IndexType cellId = tetId / NUM_TETS_PER_HEX;
      const auto& tet = meshTets[tetId];
      const double vol = clipTetAgainstSphere(tet, sphere);
      axom::atomicAdd<ExecSpace>(ovlap.data() + cellId, vol);
    });

  statistics["onSum"].set_int64(tetCount);
  statistics["clipsSum"].set_int64(tetCount);
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
