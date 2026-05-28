// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/core/numerics/matvecops.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/primal/operators/squared_distance.hpp"
#include "axom/quest/Discretize.hpp"
#include "axom/quest/detail/clipping/MonotonicZSORClipper.hpp"
#include "axom/quest/detail/clipping/TetrahedronClipUtils.hpp"
#include "axom/fmt.hpp"

#include <limits>

namespace axom
{
namespace quest
{
namespace experimental
{
namespace
{

struct RzBounds
{
  double zMin;
  double zMax;
  double rMinSquared;
  double rMaxSquared;
};

struct LinearSorData
{
  double z0;
  double r0;
  double zMin;
  double zMax;
  double drDz;
};

constexpr int LINEAR_SOR_MAX_SUBDIVISION_LEVELS = 4;
constexpr double LINEAR_SOR_EDGE_FRACTION = 0.08;
constexpr double LINEAR_CYLINDER_EDGE_FRACTION = 0.10;
constexpr double LINEAR_SOR_Z_EPS = 1e-14;

template <typename ExecSpace>
constexpr bool isCpuExecutionSpace()
{
  return std::is_same<ExecSpace, axom::SEQ_EXEC>::value
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
    || std::is_same<ExecSpace, axom::OMP_EXEC>::value
#endif
    ;
}

AXOM_HOST_DEVICE inline MeshClipperStrategy::LabelType labelRzBoundsAgainstLinearSor(
  const RzBounds& bounds,
  const LinearSorData& sor)
{
  const double zLo = bounds.zMin > sor.zMin ? bounds.zMin : sor.zMin;
  const double zHi = bounds.zMax < sor.zMax ? bounds.zMax : sor.zMax;
  if(zLo > zHi)
  {
    return MeshClipperStrategy::LabelType::LABEL_OUT;
  }

  const double rAtLo = sor.r0 + sor.drDz * (zLo - sor.z0);
  const double rAtHi = sor.r0 + sor.drDz * (zHi - sor.z0);
  const double rAllowedMax = rAtLo > rAtHi ? rAtLo : rAtHi;
  const double rAllowedMin = rAtLo < rAtHi ? rAtLo : rAtHi;
  const double rAllowedMaxSquared = rAllowedMax * rAllowedMax;
  const double rAllowedMinSquared = rAllowedMin * rAllowedMin;

  if(bounds.rMinSquared > rAllowedMaxSquared)
  {
    return MeshClipperStrategy::LabelType::LABEL_OUT;
  }

  if(bounds.zMin >= sor.zMin && bounds.zMax <= sor.zMax && bounds.rMaxSquared <= rAllowedMinSquared)
  {
    return MeshClipperStrategy::LabelType::LABEL_IN;
  }

  return MeshClipperStrategy::LabelType::LABEL_ON;
}

template <typename PolyhedronType>
AXOM_HOST_DEVICE inline RzBounds estimateRzBounds(const PolyhedronType& vertices)
{
  double zMin = vertices[0][0];
  double zMax = vertices[0][0];
  double yMin = vertices[0][1];
  double yMax = vertices[0][1];
  double xMin = vertices[0][2];
  double xMax = vertices[0][2];
  double rMaxSquared = 0.0;
  for(int i = 0; i < PolyhedronType::numVertices(); ++i)
  {
    const auto& vert = vertices[i];
    zMin = axom::utilities::min(zMin, vert[0]);
    zMax = axom::utilities::max(zMax, vert[0]);
    yMin = axom::utilities::min(yMin, vert[1]);
    yMax = axom::utilities::max(yMax, vert[1]);
    xMin = axom::utilities::min(xMin, vert[2]);
    xMax = axom::utilities::max(xMax, vert[2]);
    rMaxSquared = axom::utilities::max(rMaxSquared, vert[1] * vert[1] + vert[2] * vert[2]);
  }

  const double yClosest = yMin > 0.0 ? yMin : (yMax < 0.0 ? yMax : 0.0);
  const double xClosest = xMin > 0.0 ? xMin : (xMax < 0.0 ? xMax : 0.0);
  return {zMin, zMax, yClosest * yClosest + xClosest * xClosest, rMaxSquared};
}

AXOM_HOST_DEVICE inline void subdivideTetByMidpoints(const MeshClipperStrategy::TetrahedronType& tet,
                                                     MeshClipperStrategy::TetrahedronType children[8])
{
  const auto m01 = MeshClipperStrategy::Point3DType::midpoint(tet[0], tet[1]);
  const auto m02 = MeshClipperStrategy::Point3DType::midpoint(tet[0], tet[2]);
  const auto m03 = MeshClipperStrategy::Point3DType::midpoint(tet[0], tet[3]);
  const auto m12 = MeshClipperStrategy::Point3DType::midpoint(tet[1], tet[2]);
  const auto m13 = MeshClipperStrategy::Point3DType::midpoint(tet[1], tet[3]);
  const auto m23 = MeshClipperStrategy::Point3DType::midpoint(tet[2], tet[3]);

  children[0] = MeshClipperStrategy::TetrahedronType(tet[0], m01, m02, m03);
  children[1] = MeshClipperStrategy::TetrahedronType(m01, tet[1], m12, m13);
  children[2] = MeshClipperStrategy::TetrahedronType(m02, m12, tet[2], m23);
  children[3] = MeshClipperStrategy::TetrahedronType(m03, m13, m23, tet[3]);
  children[4] = MeshClipperStrategy::TetrahedronType(m01, m02, m03, m23);
  children[5] = MeshClipperStrategy::TetrahedronType(m01, m12, m02, m23);
  children[6] = MeshClipperStrategy::TetrahedronType(m01, m13, m12, m23);
  children[7] = MeshClipperStrategy::TetrahedronType(m01, m03, m13, m23);
}

AXOM_HOST_DEVICE inline double tetMaxEdgeSquared(const MeshClipperStrategy::TetrahedronType& tet)
{
  double maxEdgeSquared = 0.0;
  for(int i = 0; i < MeshClipperStrategy::TetrahedronType::NUM_VERTS; ++i)
  {
    for(int j = i + 1; j < MeshClipperStrategy::TetrahedronType::NUM_VERTS; ++j)
    {
      maxEdgeSquared = axom::utilities::max(maxEdgeSquared, (tet[j] - tet[i]).squared_norm());
    }
  }
  return maxEdgeSquared;
}

AXOM_HOST_DEVICE inline int chooseSubdivisionLevels(const MeshClipperStrategy::TetrahedronType& tet,
                                                    double targetEdgeSquared)
{
  double edgeSquared = tetMaxEdgeSquared(tet);
  int levels = 0;
  while(levels < LINEAR_SOR_MAX_SUBDIVISION_LEVELS && edgeSquared > targetEdgeSquared)
  {
    edgeSquared *= 0.25;
    ++levels;
  }
  return levels;
}

AXOM_HOST_DEVICE inline MeshClipperStrategy::LabelType classifyTetAgainstLinearSor(
  const MeshClipperStrategy::TetrahedronType& bodyTet,
  const LinearSorData& sor,
  double radialSquared[4])
{
  constexpr double eps = 1e-12;
  bool allInside = true;
  double zMin = bodyTet[0][0];
  double zMax = bodyTet[0][0];
  double yMin = bodyTet[0][1];
  double yMax = bodyTet[0][1];
  double xMin = bodyTet[0][2];
  double xMax = bodyTet[0][2];
  double rMaxSquared = 0.0;

  for(int i = 0; i < 4; ++i)
  {
    const auto& vert = bodyTet[i];
    zMin = axom::utilities::min(zMin, vert[0]);
    zMax = axom::utilities::max(zMax, vert[0]);
    yMin = axom::utilities::min(yMin, vert[1]);
    yMax = axom::utilities::max(yMax, vert[1]);
    xMin = axom::utilities::min(xMin, vert[2]);
    xMax = axom::utilities::max(xMax, vert[2]);

    radialSquared[i] = vert[1] * vert[1] + vert[2] * vert[2];
    rMaxSquared = axom::utilities::max(rMaxSquared, radialSquared[i]);

    const double radius = sor.r0 + sor.drDz * (vert[0] - sor.z0);
    const double radiusSquared = radius * radius;
    allInside = allInside && vert[0] >= sor.zMin - eps && vert[0] <= sor.zMax + eps &&
      radius >= 0.0 && radialSquared[i] <= radiusSquared + eps * (1.0 + radiusSquared);
  }
  if(allInside)
  {
    return MeshClipperStrategy::LabelType::LABEL_IN;
  }

  const double yClosest = yMin > 0.0 ? yMin : (yMax < 0.0 ? yMax : 0.0);
  const double xClosest = xMin > 0.0 ? xMin : (xMax < 0.0 ? xMax : 0.0);
  const RzBounds boundsInRz {zMin, zMax, yClosest * yClosest + xClosest * xClosest, rMaxSquared};
  const auto bbLabel = labelRzBoundsAgainstLinearSor(boundsInRz, sor);
  return bbLabel == MeshClipperStrategy::LabelType::LABEL_OUT
    ? MeshClipperStrategy::LabelType::LABEL_OUT
    : MeshClipperStrategy::LabelType::LABEL_ON;
}

AXOM_HOST_DEVICE inline double linearSorSignedDistance(const MeshClipperStrategy::Point3DType& pt,
                                                       double radialSquared,
                                                       const LinearSorData& sor)
{
  const double radius = sor.r0 + sor.drDz * (pt[0] - sor.z0);
  const double radialDistance = std::sqrt(radialSquared);
  return radius - radialDistance;
}

AXOM_HOST_DEVICE inline double clipBodyTetAgainstLinearSor(
  const MeshClipperStrategy::TetrahedronType& bodyTet,
  const LinearSorData& sor,
  const double radialSquared[4],
  double bodyTetVolume)
{
  using Plane3DType = MeshClipperStrategy::Plane3DType;
  using Point3DType = MeshClipperStrategy::Point3DType;
  using Vector3DType = MeshClipperStrategy::Vector3DType;

  constexpr double eps = 1e-10;

  double tetZMin = bodyTet[0][0];
  double tetZMax = bodyTet[0][0];
  for(int vi = 1; vi < MeshClipperStrategy::TetrahedronType::NUM_VERTS; ++vi)
  {
    tetZMin = axom::utilities::min(tetZMin, bodyTet[vi][0]);
    tetZMax = axom::utilities::max(tetZMax, bodyTet[vi][0]);
  }

  const Point3DType& v0 = bodyTet[0];
  const Point3DType& v1 = bodyTet[1];
  const Point3DType& v2 = bodyTet[2];
  const Point3DType& v3 = bodyTet[3];
  const double phi0 = linearSorSignedDistance(v0, radialSquared[0], sor);
  const double phi1 = linearSorSignedDistance(v1, radialSquared[1], sor);
  const double phi2 = linearSorSignedDistance(v2, radialSquared[2], sor);
  const double phi3 = linearSorSignedDistance(v3, radialSquared[3], sor);
  const double phi[4] = {phi0, phi1, phi2, phi3};

  if(tetZMin >= sor.zMin && tetZMax <= sor.zMax)
  {
    return detail::clipTetByVertexValues(phi, bodyTetVolume);
  }

  primal::Polyhedron<double, 3> overlap = primal::Polyhedron<double, 3>::from_primitive(bodyTet);
  if(tetZMin < sor.zMin)
  {
    overlap =
      primal::clip(overlap,
                   Plane3DType(Vector3DType {1.0, 0.0, 0.0}, Point3DType {sor.zMin, 0.0, 0.0}),
                   eps);
    if(overlap.numVertices() < 4)
    {
      return 0.0;
    }
  }

  if(tetZMax > sor.zMax)
  {
    overlap =
      primal::clip(overlap,
                   Plane3DType(Vector3DType {-1.0, 0.0, 0.0}, Point3DType {sor.zMax, 0.0, 0.0}),
                   eps);
    if(overlap.numVertices() < 4)
    {
      return 0.0;
    }
  }

  const Vector3DType e1 = v1 - v0;
  const Vector3DType e2 = v2 - v0;
  const Vector3DType e3 = v3 - v0;

  const double denom = Vector3DType::scalar_triple_product(e1, e2, e3);
  if(axom::utilities::isNearlyEqual(denom, 0.0, eps))
  {
    return 0.0;
  }

  const Vector3DType grad = ((phi1 - phi0) * Vector3DType::cross_product(e2, e3) +
                             (phi2 - phi0) * Vector3DType::cross_product(e3, e1) +
                             (phi3 - phi0) * Vector3DType::cross_product(e1, e2)) /
    denom;
  const double gradSqNorm = grad.squared_norm();
  if(gradSqNorm <= 1e-20)
  {
    return 0.0;
  }

  const Point3DType planePoint = v0 - (phi0 / gradSqNorm) * grad;
  overlap = primal::clip(overlap, Plane3DType(grad, planePoint), eps);
  return overlap.volume();
}

}  // namespace

MonotonicZSORClipper::MonotonicZSORClipper(const klee::Geometry& kGeom, const std::string& name)
  : MeshClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("FSor") : name)
  , m_maxRadius(0.0)
  , m_minRadius(numerics::floating_point_limits<double>::max())
  , m_transformer()
  , m_bodyVertexCache(std::make_shared<BodyVertexCache>())
{
  extractClipperInfo();

  combineRadialSegments(m_sorCurve);
  axom::Array<axom::IndexType> turnIndices = findZSwitchbacks(m_sorCurve.view());
  if(turnIndices.size() > 2)
  {
    // The 2 "turns" allowed are the first and last points.  Anything else is a switchback.
    SLIC_ERROR(
      "MonotonicZSORClipper does not work when a curve doubles back"
      " in the axial direction.  Use SORClipper instead.");
  }

  for(auto& pt : m_sorCurve)
  {
    m_maxRadius = fmax(m_maxRadius, pt[1]);
    m_minRadius = fmin(m_minRadius, pt[1]);
  }
  SLIC_ERROR_IF(m_minRadius < 0.0,
                axom::fmt::format("MonotonicZSORClipper '{}' has a negative radius", m_name));

  // Combine internal and external rotations into m_transformer.
  m_transformer.applyRotation(Vector3DType({1, 0, 0}), m_sorDirection);
  m_transformer.applyTranslation(m_sorOrigin.array());
  m_transformer.applyMatrix(m_extTrans);
  m_invTransformer = m_transformer.getInverse();

  for(const auto& pt : m_sorCurve)
  {
    m_curveBb.addPoint(pt);
  }
}

MonotonicZSORClipper::MonotonicZSORClipper(const klee::Geometry& kGeom,
                                           const std::string& name,
                                           axom::ArrayView<const Point2DType> discreteFunction,
                                           const Point3DType& sorOrigin,
                                           const Vector3DType& sorDirection,
                                           axom::IndexType levelOfRefinement,
                                           std::shared_ptr<BodyVertexCache> bodyVertexCache)
  : MeshClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("FSor") : name)
  , m_sorCurve(discreteFunction, axom::execution_space<axom::SEQ_EXEC>::allocatorID())
  , m_maxRadius(0.0)
  , m_minRadius(numerics::floating_point_limits<double>::max())
  , m_sorOrigin(sorOrigin)
  , m_sorDirection(sorDirection)
  , m_levelOfRefinement(levelOfRefinement)
  , m_transformer()
  , m_bodyVertexCache(bodyVertexCache ? bodyVertexCache : std::make_shared<BodyVertexCache>())
{
  combineRadialSegments(m_sorCurve);
  axom::Array<axom::IndexType> turnIndices = findZSwitchbacks(m_sorCurve.view());
  if(turnIndices.size() > 2)
  {
    // The 2 "turns" allowed are the first and last points.  Anything else is a switchback.
    SLIC_ERROR(
      "MonotonicZSORClipper does not work when a curve doubles back"
      " in the axial direction.  Use SORClipper instead.");
  }

  for(auto& pt : m_sorCurve)
  {
    m_maxRadius = fmax(m_maxRadius, pt[1]);
    m_minRadius = fmin(m_minRadius, pt[1]);
  }
  SLIC_ERROR_IF(m_minRadius < 0.0,
                axom::fmt::format("MonotonicZSORClipper '{}' has a negative radius", m_name));

  // Combine internal and external rotations into m_transformer.
  m_transformer.applyRotation(Vector3DType({1, 0, 0}), m_sorDirection);
  m_transformer.applyTranslation(m_sorOrigin.array());
  m_transformer.applyMatrix(m_extTrans);
  m_invTransformer = m_transformer.getInverse();

  for(const auto& pt : m_sorCurve)
  {
    m_curveBb.addPoint(pt);
  }
}

bool MonotonicZSORClipper::labelCellsInOut(quest::experimental::ShapeMesh& shapeMesh,
                                           axom::Array<LabelType>& labels)
{
  SLIC_ERROR_IF(shapeMesh.dimension() != 3, "MonotonicZSORClipper requires a 3D mesh.");

  const int allocId = shapeMesh.getAllocatorID();
  const auto cellCount = shapeMesh.getCellCount();
  if(labels.size() < cellCount || labels.getAllocatorID() != allocId)
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

bool MonotonicZSORClipper::labelTetsInOut(quest::experimental::ShapeMesh& shapeMesh,
                                          axom::ArrayView<const axom::IndexType> cellIds,
                                          axom::Array<LabelType>& tetLabels)
{
  SLIC_ERROR_IF(shapeMesh.dimension() != 3, "MonotonicZSORClipper requires a 3D mesh.");

  const int allocId = shapeMesh.getAllocatorID();
  const auto cellCount = cellIds.size();
  const auto tetCount = cellCount * NUM_TETS_PER_HEX;
  if(tetLabels.size() < tetCount || tetLabels.getAllocatorID() != allocId)
  {
    tetLabels = axom::Array<LabelType>(ArrayOptions::Uninitialized(), tetCount, tetCount, allocId);
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

bool MonotonicZSORClipper::specializedClipTets(quest::experimental::ShapeMesh& shapeMesh,
                                               axom::ArrayView<double> ovlap,
                                               const axom::ArrayView<IndexType>& tetIds,
                                               conduit::Node& statistics)
{
  if(m_sorCurve.size() != 2 || m_levelOfRefinement > 5)
  {
    return false;
  }

  if(axom::utilities::isNearlyEqual(m_sorCurve[1][0] - m_sorCurve[0][0], 0.0, LINEAR_SOR_Z_EPS))
  {
    return false;
  }

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
    return false;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    return false;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
  return true;
}

template <typename ExecSpace>
void MonotonicZSORClipper::specializedClipTetsImpl(quest::experimental::ShapeMesh& shapeMesh,
                                                   axom::ArrayView<double> ovlap,
                                                   const axom::ArrayView<IndexType>& tetIds,
                                                   conduit::Node& statistics)
{
  AXOM_ANNOTATE_SCOPE("MonotonicZSORClipper::adaptive_clip");
  struct WorkTet
  {
    TetrahedronType bodyTet;
    int depth;
  };

  const double z0 = m_sorCurve[0][0];
  const double r0 = m_sorCurve[0][1];
  const double z1 = m_sorCurve[1][0];
  const double r1 = m_sorCurve[1][1];
  const LinearSorData sorData {z0, r0, z0 < z1 ? z0 : z1, z0 < z1 ? z1 : z0, (r1 - r0) / (z1 - z0)};
  const auto invTransformer = m_invTransformer;
  auto meshTets = shapeMesh.getCellsAsTets();
  const IndexType tetCount = tetIds.size();

  axom::ReduceSum<ExecSpace, std::int64_t> inSum {0};
  axom::ReduceSum<ExecSpace, std::int64_t> onSum {0};
  axom::ReduceSum<ExecSpace, std::int64_t> outSum {0};

  constexpr int childCount = 8;
  constexpr int maxStackSize = 1 + (childCount - 1) * LINEAR_SOR_MAX_SUBDIVISION_LEVELS;
  const double maxRadius = axom::utilities::max(r0, r1);
  const double edgeFraction = axom::utilities::isNearlyEqual(sorData.drDz, 0.0)
    ? LINEAR_CYLINDER_EDGE_FRACTION
    : LINEAR_SOR_EDGE_FRACTION;
  const double targetEdge = maxRadius * edgeFraction;
  const double targetEdgeSquared = targetEdge * targetEdge;

  AXOM_ANNOTATE_BEGIN("MonotonicZSORClipper::transform_subdivide_clip");
  axom::for_all<ExecSpace>(tetCount, [=] AXOM_HOST_DEVICE(axom::IndexType ti) {
    const axom::IndexType tetId = tetIds[ti];
    const axom::IndexType cellId = tetId / NUM_TETS_PER_HEX;
    SLIC_ASSERT(cellId >= 0 && cellId < ovlap.size());
    const TetrahedronType worldTet = meshTets[tetId];

    TetrahedronType bodyTet = worldTet;
    for(int vi = 0; vi < 4; ++vi)
    {
      invTransformer.transform(bodyTet[vi].array());
    }
    bodyTet.checkAndFixOrientation();
    const double bodyRootVolume = bodyTet.volume();
    if(bodyRootVolume <= 0.0)
    {
      return;
    }
    const double worldRootVolume = worldTet.volume();
    const double volumeScale = worldRootVolume / bodyRootVolume;
    const int maxDepth = chooseSubdivisionLevels(bodyTet, targetEdgeSquared);
    double levelVolumes[LINEAR_SOR_MAX_SUBDIVISION_LEVELS + 1];
    double bodyLevelVolumes[LINEAR_SOR_MAX_SUBDIVISION_LEVELS + 1];
    levelVolumes[0] = worldRootVolume;
    bodyLevelVolumes[0] = bodyRootVolume;
    for(int depth = 1; depth <= maxDepth; ++depth)
    {
      levelVolumes[depth] = levelVolumes[depth - 1] * 0.125;
      bodyLevelVolumes[depth] = bodyLevelVolumes[depth - 1] * 0.125;
    }

    WorkTet stack[maxStackSize];
    int stackSize = 1;
    stack[0] = {bodyTet, 0};

    double overlap = 0.0;
    while(stackSize > 0)
    {
      const WorkTet current = stack[--stackSize];
      double radialSquared[4];
      const LabelType label = classifyTetAgainstLinearSor(current.bodyTet, sorData, radialSquared);

      if(label == LabelType::LABEL_IN)
      {
        overlap += levelVolumes[current.depth];
        inSum += 1;
        continue;
      }
      if(label == LabelType::LABEL_OUT)
      {
        outSum += 1;
        continue;
      }

      if(current.depth == maxDepth)
      {
        const double bodyOverlap = clipBodyTetAgainstLinearSor(current.bodyTet,
                                                               sorData,
                                                               radialSquared,
                                                               bodyLevelVolumes[current.depth]);
        overlap += bodyOverlap * volumeScale;
        onSum += 1;
        continue;
      }

      TetrahedronType bodyChildren[childCount];
      subdivideTetByMidpoints(current.bodyTet, bodyChildren);
      for(int child = 0; child < childCount; ++child)
      {
        stack[stackSize++] = {bodyChildren[child], current.depth + 1};
      }
    }

    detail::addToOverlapVolume<ExecSpace>(ovlap.data() + cellId, overlap);
  });
  AXOM_ANNOTATE_END("MonotonicZSORClipper::transform_subdivide_clip");

  AXOM_ANNOTATE_BEGIN("MonotonicZSORClipper::record_statistics");
  const std::int64_t clipsInCount = inSum.get();
  const std::int64_t clipsOnCount = onSum.get();
  const std::int64_t clipsOutCount = outSum.get();
  statistics["clipsIn"].set_int64(clipsInCount);
  statistics["clipsOn"].set_int64(clipsOnCount);
  statistics["clipsOut"].set_int64(clipsOutCount);
  statistics["clipsMiss"].set_int64(0);
  statistics["clipsSum"].set_int64(clipsInCount + clipsOnCount + clipsOutCount);
  statistics["clipsCandidates"].set(static_cast<IndexType>(tetCount));
  AXOM_ANNOTATE_END("MonotonicZSORClipper::record_statistics");
}

/*
 * Implementation: (reverse) transform the mesh vertices to the r-z
 * frame where the curve is defined as a r(z) function.  It's easier to
 * determine whether the point is in the sor that way.
 */
template <typename ExecSpace>
void MonotonicZSORClipper::labelCellsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                                               axom::ArrayView<LabelType> labels)
{
  if constexpr(isCpuExecutionSpace<ExecSpace>())
  {
    if(m_sorCurve.size() == 2 &&
       !axom::utilities::isNearlyEqual(m_sorCurve[1][0] - m_sorCurve[0][0], 0.0, LINEAR_SOR_Z_EPS))
    {
      const double z0 = m_sorCurve[0][0];
      const double r0 = m_sorCurve[0][1];
      const double z1 = m_sorCurve[1][0];
      const double r1 = m_sorCurve[1][1];
      const LinearSorData sorData {z0, r0, z0 < z1 ? z0 : z1, z0 < z1 ? z1 : z0, (r1 - r0) / (z1 - z0)};

      const auto cellCount = shapeMesh.getCellCount();
      const auto vertCount = shapeMesh.getVertexCount();
      const auto& vertCoords = shapeMesh.getVertexCoords3D();
      const auto& vX = vertCoords[0];
      const auto& vY = vertCoords[1];
      const auto& vZ = vertCoords[2];
      const auto connView = shapeMesh.getCellNodeConnectivity();
      auto meshCellVolumes = shapeMesh.getCellVolumes();
      auto invTransformer = m_invTransformer;
      constexpr double EPS = 1e-10;

      auto& bodyVertices = m_bodyVertexCache->vertices;
      const bool bodyVerticesCached = m_bodyVertexCache->shapeMesh == &shapeMesh &&
        m_bodyVertexCache->sourceCoordinates[0] == vX.data() &&
        m_bodyVertexCache->sourceCoordinates[1] == vY.data() &&
        m_bodyVertexCache->sourceCoordinates[2] == vZ.data() && bodyVertices.size() == vertCount &&
        bodyVertices.getAllocatorID() == shapeMesh.getAllocatorID();
      AXOM_ANNOTATE_BEGIN("MonotonicZSORClipper::transform_vertices");
      if(!bodyVerticesCached)
      {
        bodyVertices = axom::Array<Point3DType>(ArrayOptions::Uninitialized(),
                                                vertCount,
                                                vertCount,
                                                shapeMesh.getAllocatorID());
        auto bodyVerticesView = bodyVertices.view();
        axom::for_all<ExecSpace>(vertCount, [=] AXOM_HOST_DEVICE(axom::IndexType vertId) {
          Point3DType bodyVertex {vX[vertId], vY[vertId], vZ[vertId]};
          invTransformer.transform(bodyVertex.array());
          bodyVerticesView[vertId] = bodyVertex;
        });
        m_bodyVertexCache->shapeMesh = &shapeMesh;
        m_bodyVertexCache->sourceCoordinates[0] = vX.data();
        m_bodyVertexCache->sourceCoordinates[1] = vY.data();
        m_bodyVertexCache->sourceCoordinates[2] = vZ.data();
      }
      AXOM_ANNOTATE_END("MonotonicZSORClipper::transform_vertices");
      auto bodyVerticesView = bodyVertices.view();

      AXOM_ANNOTATE_BEGIN("MonotonicZSORClipper::classify_cells_linear");
      axom::for_all<ExecSpace>(cellCount, [=] AXOM_HOST_DEVICE(axom::IndexType cellId) {
        if(axom::utilities::isNearlyEqual(meshCellVolumes[cellId], 0.0, EPS))
        {
          labels[cellId] = LabelType::LABEL_OUT;
          return;
        }

        HexahedronType cellHex;
        const auto cellVertIds = connView[cellId];
        for(int vi = 0; vi < HexahedronType::NUM_HEX_VERTS; ++vi)
        {
          cellHex[vi] = bodyVerticesView[cellVertIds[vi]];
        }

        labels[cellId] = labelRzBoundsAgainstLinearSor(estimateRzBounds(cellHex), sorData);
      });
      AXOM_ANNOTATE_END("MonotonicZSORClipper::classify_cells_linear");
      return;
    }
  }

  axom::Array<BoundingBox2DType> bbOn;
  axom::Array<BoundingBox2DType> bbUnder;
  AXOM_ANNOTATE_BEGIN("MonotonicZSORClipper::compute_curve_boxes");
  computeCurveBoxes<ExecSpace>(shapeMesh, bbOn, bbUnder);
  AXOM_ANNOTATE_END("MonotonicZSORClipper::compute_curve_boxes");
  const axom::ArrayView<const BoundingBox2DType> bbOnView = bbOn.view();
  const axom::ArrayView<const BoundingBox2DType> bbUnderView = bbUnder.view();

  const auto cellCount = shapeMesh.getCellCount();
  auto meshHexes = shapeMesh.getCellsAsHexes();
  auto meshCellVolumes = shapeMesh.getCellVolumes();
  auto invTransformer = m_invTransformer;
  constexpr double EPS = 1e-10;

  AXOM_ANNOTATE_BEGIN("MonotonicZSORClipper::classify_cells_generic");
  axom::for_all<ExecSpace>(cellCount, [=] AXOM_HOST_DEVICE(axom::IndexType cellId) {
    if(axom::utilities::isNearlyEqual(meshCellVolumes[cellId], 0.0, EPS))
    {
      labels[cellId] = LabelType::LABEL_OUT;
      return;
    }
    auto cellHex = meshHexes[cellId];
    for(int vi = 0; vi < HexahedronType::NUM_HEX_VERTS; ++vi)
    {
      invTransformer.transform(cellHex[vi].array());
    }
    BoundingBox2DType cellBbInRz = estimateBoundingBoxInRz(cellHex);
    labels[cellId] = rzBbToLabel(cellBbInRz, bbOnView, bbUnderView);
  });
  AXOM_ANNOTATE_END("MonotonicZSORClipper::classify_cells_generic");
}

template <typename ExecSpace>
void MonotonicZSORClipper::labelTetsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                                              axom::ArrayView<const axom::IndexType> cellIds,
                                              axom::ArrayView<LabelType> labels)
{
  if constexpr(isCpuExecutionSpace<ExecSpace>())
  {
    if(m_sorCurve.size() == 2 &&
       !axom::utilities::isNearlyEqual(m_sorCurve[1][0] - m_sorCurve[0][0], 0.0, LINEAR_SOR_Z_EPS))
    {
      const double z0 = m_sorCurve[0][0];
      const double r0 = m_sorCurve[0][1];
      const double z1 = m_sorCurve[1][0];
      const double r1 = m_sorCurve[1][1];
      const LinearSorData sorData {z0, r0, z0 < z1 ? z0 : z1, z0 < z1 ? z1 : z0, (r1 - r0) / (z1 - z0)};

      const auto cellCount = cellIds.size();
      auto meshHexes = shapeMesh.getCellsAsHexes();
      const auto connView = shapeMesh.getCellNodeConnectivity();
      auto tetVolumes = shapeMesh.getTetVolumes();
      auto invTransformer = m_invTransformer;
      const auto& bodyVertices = m_bodyVertexCache->vertices;
      const bool bodyVerticesCached = m_bodyVertexCache->shapeMesh == &shapeMesh &&
        bodyVertices.size() == shapeMesh.getVertexCount() &&
        bodyVertices.getAllocatorID() == shapeMesh.getAllocatorID();
      const auto bodyVerticesView = bodyVertices.view();
      constexpr double EPS = 1e-10;

      AXOM_ANNOTATE_BEGIN("MonotonicZSORClipper::classify_tets_linear");
      axom::for_all<ExecSpace>(cellCount, [=] AXOM_HOST_DEVICE(axom::IndexType ci) {
        axom::IndexType cellId = cellIds[ci];

        HexahedronType hex;
        if(bodyVerticesCached)
        {
          const auto cellVertIds = connView[cellId];
          for(int vi = 0; vi < HexahedronType::NUM_HEX_VERTS; ++vi)
          {
            hex[vi] = bodyVerticesView[cellVertIds[vi]];
          }
        }
        else
        {
          hex = meshHexes[cellId];
          for(int vi = 0; vi < HexahedronType::NUM_HEX_VERTS; ++vi)
          {
            invTransformer.transform(hex[vi].array());
          }
        }

        TetrahedronType cellTets[NUM_TETS_PER_HEX];
        ShapeMesh::hexToTets(hex, cellTets);

        for(IndexType ti = 0; ti < NUM_TETS_PER_HEX; ++ti)
        {
          axom::IndexType tetId = cellId * NUM_TETS_PER_HEX + ti;
          LabelType& tetLabel = labels[ci * NUM_TETS_PER_HEX + ti];
          if(axom::utilities::isNearlyEqual(tetVolumes[tetId], 0.0, EPS))
          {
            tetLabel = LabelType::LABEL_OUT;
            continue;
          }

          const TetrahedronType& tet = cellTets[ti];
          tetLabel = labelRzBoundsAgainstLinearSor(estimateRzBounds(tet), sorData);
        }
      });
      AXOM_ANNOTATE_END("MonotonicZSORClipper::classify_tets_linear");
      return;
    }
  }

  axom::Array<BoundingBox2DType> bbOn;
  axom::Array<BoundingBox2DType> bbUnder;
  AXOM_ANNOTATE_BEGIN("MonotonicZSORClipper::compute_curve_boxes");
  computeCurveBoxes<ExecSpace>(shapeMesh, bbOn, bbUnder);
  AXOM_ANNOTATE_END("MonotonicZSORClipper::compute_curve_boxes");
  const axom::ArrayView<const BoundingBox2DType> bbOnView = bbOn.view();
  const axom::ArrayView<const BoundingBox2DType> bbUnderView = bbUnder.view();

  const auto cellCount = cellIds.size();
  auto meshHexes = shapeMesh.getCellsAsHexes();
  auto tetVolumes = shapeMesh.getTetVolumes();
  auto invTransformer = m_invTransformer;
  constexpr double EPS = 1e-10;

  AXOM_ANNOTATE_BEGIN("MonotonicZSORClipper::classify_tets_generic");
  axom::for_all<ExecSpace>(cellCount, [=] AXOM_HOST_DEVICE(axom::IndexType ci) {
    axom::IndexType cellId = cellIds[ci];

    HexahedronType hex = meshHexes[cellId];
    for(int vi = 0; vi < HexahedronType::NUM_HEX_VERTS; ++vi)
    {
      invTransformer.transform(hex[vi].array());
    }

    TetrahedronType cellTets[NUM_TETS_PER_HEX];
    ShapeMesh::hexToTets(hex, cellTets);

    for(IndexType ti = 0; ti < NUM_TETS_PER_HEX; ++ti)
    {
      axom::IndexType tetId = cellId * NUM_TETS_PER_HEX + ti;
      LabelType& tetLabel = labels[ci * NUM_TETS_PER_HEX + ti];
      if(axom::utilities::isNearlyEqual(tetVolumes[tetId], 0.0, EPS))
      {
        tetLabel = LabelType::LABEL_OUT;
        continue;
      }
      const TetrahedronType& tet = cellTets[ti];
      BoundingBox2DType bbInRz = estimateBoundingBoxInRz(tet);
      tetLabel = rzBbToLabel(bbInRz, bbOnView, bbUnderView);
    }
  });
  AXOM_ANNOTATE_END("MonotonicZSORClipper::classify_tets_generic");
}

/*
  Compute bounding box in rz space for a tet or hex geometry in
  body frame (the 3D frame with the rotation along +x).

  1. Rotate the tet or hex vertices into the rz plane.
  2. Compute bounding box for vertices.
  3. Expand 2D bounding box to contain edge that may
     intersect SOR between vertices.
*/
template <typename PolyhedronType>
AXOM_HOST_DEVICE MonotonicZSORClipper::BoundingBox2DType MonotonicZSORClipper::estimateBoundingBoxInRz(
  const PolyhedronType& vertices)
{
  MonotonicZSORClipper::BoundingBox2DType bbInRz;

  // Range of vertex angles in cylindrical coordinates.
  primal::BoundingBox<double, 1> angleRange;

  for(IndexType vi = 0; vi < vertices.numVertices(); ++vi)
  {
    auto& vert = vertices[vi];
    Point2DType vertOnXPlane {vert[1], vert[2]};
    Point2DType vertOnRz {
      vert[0],
      std::sqrt(numerics::dot_product(vertOnXPlane.data(), vertOnXPlane.data(), 2))};
    bbInRz.addPoint(vertOnRz);

    double angle = atan2(vertOnXPlane[1], vertOnXPlane[0]);
    angleRange.addPoint(primal::Point<double, 1> {angle});
  }
  /*
    The geometry can be closer to the rotation axis than its
    individual vertices are, depending on the angle (about the axis)
    between the vertices.  Given the angle, extend the bottom of bbInRz
    for the worst case.
  */
  auto angleDiff = angleRange.range()[0];
  double factor = angleDiff > M_PI ? 0.0 : cos(angleDiff / 2);
  auto newMin = bbInRz.getMin();
  newMin[1] *= factor;
  bbInRz.addPoint(newMin);
  return bbInRz;
}

/*
  Compute label based on a bounding box in rz space.

  - If bbInRz is close to any bbOn, label it ON.
  - Else if bbInRz touches any bbUnder, label it IN.
    It cannot possibly be partially outside, because it
    doesn't cross the boundary or even touch any bbOn.
  - Else, label bbInRz OUT.

  We expect bbOn and bbUnder to be small arrays, so we use
  linear searches.  If that's too slow, we can use a BVH.
*/
AXOM_HOST_DEVICE inline MeshClipperStrategy::LabelType MonotonicZSORClipper::rzBbToLabel(
  const BoundingBox2DType& bbInRz,
  const axom::ArrayView<const BoundingBox2DType>& bbOn,
  const axom::ArrayView<const BoundingBox2DType>& bbUnder)
{
  LabelType label = LabelType::LABEL_OUT;

  for(const auto& bbOn : bbOn)
  {
    double sqDist = axom::primal::squared_distance(bbInRz, bbOn);
    if(sqDist <= 0.0)
    {
      label = LabelType::LABEL_ON;
    }
  }

  if(label == LabelType::LABEL_OUT)
  {
    for(const auto& bbUnder : bbUnder)
    {
      if(bbInRz.intersectsWith(bbUnder))
      {
        label = LabelType::LABEL_IN;
      }
    }
  }

  return label;
}

/*
*/
template <typename ExecSpace>
void MonotonicZSORClipper::computeCurveBoxes(quest::experimental::ShapeMesh& shapeMesh,
                                             axom::Array<BoundingBox2DType>& bbOn,
                                             axom::Array<BoundingBox2DType>& bbUnder)
{
  /*
   * Compute bounding boxes bbOn, which cover the curve segments, and
   * bbUnder, which cover the space between bbOn and the z axis.  bbOn
   * includes end caps, the segments that join the curve to the
   * z-axis.
   */
  const int allocId = shapeMesh.getAllocatorID();
  const IndexType cellCount = shapeMesh.getCellCount();

  axom::ArrayView<const double> cellLengths = shapeMesh.getCellLengths();

  axom::ReduceSum<ExecSpace, double> sumCharLength(0.0);
  axom::for_all<ExecSpace>(cellCount, [=] AXOM_HOST_DEVICE(axom::IndexType cellId) {
    sumCharLength += cellLengths[cellId];
  });
  double avgCharLength = sumCharLength.get() / cellCount;

  /*
    Subdivide the SOR curve and place it with the correct allocator.
    Create temporary sorCurve that is equivalent to m_sorCurve but
    - with long segments subdivided into subsegments based on
      characteristic length of mesh cells.
    - with memory from allocId.
  */
  axom::Array<Point2DType> sorCurve = subdivideCurve(m_sorCurve,
                                                     3 * avgCharLength /* maxMean */,
                                                     -1 /* maxDz, negative disables */,
                                                     -1 /* minDz, negative disables */);
  sorCurve = axom::Array<Point2DType>(sorCurve, allocId);
  auto sorCurveView = sorCurve.view();

  /*
    Compute 2 sets of boxes.
    - bbOn have boxes over each segment.
    - bbUnder have boxes from the z axis to the bottom of bbOn.
    Add add to bbOn boxes representing the vertical endcaps of the curve.
  */
  auto segCount = sorCurve.size() - 1;
  bbOn = axom::Array<BoundingBox2DType>(segCount + 2, segCount + 2, allocId);
  bbUnder = axom::Array<BoundingBox2DType>(segCount, segCount, allocId);
  auto bbOnView = bbOn.view();
  auto bbUnderView = bbUnder.view();

  axom::for_all<ExecSpace>(segCount, [=] AXOM_HOST_DEVICE(axom::IndexType i) {
    BoundingBox2DType& on = bbOnView[i];
    BoundingBox2DType& under = bbUnderView[i];
    on.addPoint(sorCurveView[i]);
    on.addPoint(sorCurveView[i + 1]);
    Point2DType underMin {on.getMin()[0], 0.0};
    Point2DType underMax {on.getMax()[0], on.getMin()[1]};
    under = BoundingBox2DType(underMin, underMax);
  });

  axom::Array<BoundingBox2DType> endCaps(2, 2);
  endCaps[0].addPoint(m_sorCurve.front());
  endCaps[0].addPoint(Point2DType {m_sorCurve.front()[0], 0.0});
  endCaps[1].addPoint(m_sorCurve.back());
  endCaps[1].addPoint(Point2DType {m_sorCurve.back()[0], 0.0});
  axom::copy(&bbOn[segCount], endCaps.data(), endCaps.size() * sizeof(BoundingBox2DType));
}

/*
 * Replace SOR curve segments that have bounding boxes that overlap
 * too much beyond what the segments actually overlap.
 *
 * Goal: Split up segments with excessively large bounding boxes,
 * which reach too far beyond the SOR curve.  These are long diagonal
 * segments.  But don't split up segments aligned close to z or r
 * directions, because they don't have excessively large bounding
 * boxes for their size.  We do this by limiting the harmonic mean of
 * the r and z sides of the bounding boxes.
 */
Array<MonotonicZSORClipper::Point2DType> MonotonicZSORClipper::subdivideCurve(
  const Array<Point2DType>& sorCurveIn,
  double maxMean,
  double maxDz,
  double minDz)
{
  Array<Point2DType> sorCurveOut;

  if(sorCurveIn.empty())
  {
    return sorCurveOut;
  }

  // Reserve guessed total number of points needed
  sorCurveOut.reserve(sorCurveIn.size() * 1.2 + 10);
  sorCurveOut.push_back(sorCurveIn[0]);

  for(IndexType i = 1; i < sorCurveIn.size(); ++i)
  {
    const Point2DType& segStart = sorCurveIn[i - 1];
    const Point2DType& segEnd = sorCurveIn[i];

    const auto delta = segEnd.array() - segStart.array();
    const auto absDelta = axom::abs(delta);
    const double segDz = absDelta[0];
    const double segDr = absDelta[1];
    const double segMean = 2 * segDz * segDr / (segDz + segDr);

    int numSplitsByMean =
      maxMean <= 0 && segMean > maxMean ? 0 : static_cast<int>(std::ceil(segMean / maxMean)) - 1;
    int numSplitsByDz =
      maxDz <= 0 && segDz > maxDz ? 0 : static_cast<int>(std::ceil(segDz / maxDz)) - 1;

    // Prevent dz from falling below minDz
    int numSplitsByMinDz = minDz <= 0 && segDz > minDz ? 0 : static_cast<int>(segDz / minDz) - 1;

    int numSplits = std::min(std::max(numSplitsByMean, numSplitsByDz), numSplitsByMinDz);

    for(int j = 1; j < numSplits; ++j)
    {
      double t = static_cast<double>(j) / numSplits;
      Point2DType newPt(segStart.array() + t * delta);
      sorCurveOut.push_back(newPt);
    }
    sorCurveOut.push_back(segEnd);
  }

  return sorCurveOut;
}

bool MonotonicZSORClipper::getGeometryAsOcts(quest::experimental::ShapeMesh& shapeMesh,
                                             axom::Array<OctahedronType>& octs)
{
  AXOM_ANNOTATE_SCOPE("MonotonicZSORClipper::getGeometryAsOcts");
  switch(shapeMesh.getRuntimePolicy())
  {
  case axom::runtime_policy::Policy::seq:
    getGeometryAsOctsImpl<axom::SEQ_EXEC>(shapeMesh, octs);
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    getGeometryAsOctsImpl<axom::OMP_EXEC>(shapeMesh, octs);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case axom::runtime_policy::Policy::cuda:
    getGeometryAsOctsImpl<axom::CUDA_EXEC<256>>(shapeMesh, octs);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    getGeometryAsOctsImpl<axom::HIP_EXEC<256>>(shapeMesh, octs);
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
  return true;
}

/*
  Compute octahedral geometry representation, with an execution policy.

  Side effect: m_sorCurve data is reallocated to the shapeMesh allocator,
  if it's not there yet.
*/
template <typename ExecSpace>
bool MonotonicZSORClipper::getGeometryAsOctsImpl(quest::experimental::ShapeMesh& shapeMesh,
                                                 axom::Array<OctahedronType>& octs)
{
  const int allocId = shapeMesh.getAllocatorID();
  octs = axom::Array<OctahedronType>(0, 0, allocId);

  const auto cellCount = shapeMesh.getCellCount();

  // Compute an average characteristic length for the mesh cells.
  axom::ArrayView<const double> cellVolumes = shapeMesh.getCellVolumes();
  axom::ReduceSum<ExecSpace, double> sumVolume(0.0);
  axom::for_all<ExecSpace>(cellCount, [=] AXOM_HOST_DEVICE(axom::IndexType cellId) {
    sumVolume += cellVolumes[cellId];
  });
  double avgVolume = sumVolume.get() / cellCount;
  double avgCharLength = pow(avgVolume, 1. / 3);

  axom::Array<Point2DType> sorCurve = subdivideCurve(m_sorCurve,
                                                     3 * avgCharLength /* maxMean */,
                                                     3 * avgCharLength /* maxDz */,
                                                     2 * avgCharLength /* minDz */);

  // Generate the Octahedra
  int octCount = 0;
  const bool good = axom::quest::discretize<ExecSpace>(sorCurve.view(),
                                                       int(sorCurve.size()),
                                                       m_levelOfRefinement,
                                                       octs,
                                                       octCount);

  AXOM_UNUSED_VAR(good);
  SLIC_ASSERT(good);
  SLIC_ASSERT(octCount == octs.size());

  auto transformer = m_transformer;
  auto octsView = octs.view();
  axom::for_all<ExecSpace>(octCount, [=] AXOM_HOST_DEVICE(axom::IndexType iOct) {
    OctahedronType& oct = octsView[iOct];
    for(int iVert = 0; iVert < OctahedronType::NUM_VERTS; ++iVert)
    {
      transformer.transform(oct[iVert].array());
    }
  });

  SLIC_DEBUG(axom::fmt::format(
    "MonotonicZSORClipper '{}' {}-level refinement got {} geometry octs from {} curve points.",
    name(),
    m_levelOfRefinement,
    octs.size(),
    sorCurve.size()));

  return true;
}

/*
  Combine consecutive radial segments in SOR curve.  Change in place.
*/
void MonotonicZSORClipper::combineRadialSegments(axom::Array<Point2DType>& sorCurve)
{
  int ptCount = sorCurve.size();
  if(ptCount < 3)
  {
    return;
  }

  constexpr double eps = 1e-14;

  // Set sorCurve[j] to sorCurve[i] where j <= i, skipping points
  // joining consecutive radial segments.

  int j = 1;
  bool prevIsRadial = axom::utilities::isNearlyEqual(sorCurve[j][0] - sorCurve[j - 1][0], eps);
  bool curIsRadial = false;
  for(int i = 2; i < ptCount; ++i)
  {
    curIsRadial = axom::utilities::isNearlyEqual(sorCurve[i][0] - sorCurve[i - 1][0], eps);
    /*
      Current and previous segments share point j.  If both are
      consecutive radial segments, discard point j by overwriting it
      with point i.  Else, copy point i to a new point j.
    */
    if(!(curIsRadial && prevIsRadial))
    {
      ++j;
    }
    sorCurve[j] = sorCurve[i];
    prevIsRadial = curIsRadial;
  }
  sorCurve.resize(j + 1);
}

/*
  Find points along the r(z) curve where the z-coordinate changes direction.

  Cases 1 and 2 below show direction changes at point o.  Case 3
  shows a potential change at the radial segment, but not a real
  change.  (Radial segments have constant z and align with the radial
  direction.)  To decide between cases 2 and 3, defer until the
  segment after the radial segment.  (The next segment is not radial
  because adjacent radials have been combined by combineRadialSegments.)
  For case 2, prefer to split at the point closer to the axis of
  rotation.

     r   ^
  (or y) |    (1)         (2)         (3)
         |  Single      Radial      Radial
         |  point       segment     segment w/o
         |  change      change      change
         |
         |    \            \          \
         |     \            \          \
         |      o            |          |
         |     /             o           \
         |    /             /             \
         +-------------------------------------> z (or x)
*/
axom::Array<axom::IndexType> MonotonicZSORClipper::findZSwitchbacks(
  axom::ArrayView<const Point2DType> pts)
{
  const axom::IndexType segCount = pts.size() - 1;
  SLIC_ASSERT(segCount > 0);

  // boundaryIdx is where curve's axial direction changes, plus end points.
  axom::Array<axom::IndexType> boundaryIdx(0, 2);
  boundaryIdx.push_back(0);

  constexpr double eps = 1e-14;

  if(segCount > 1)
  {
    // Direction is whether z increases or decreases along the curve.
    // curDir is the current direction, ignoring radial segments,
    // which don't change z.
    int curDir = axom::utilities::sign_of(pts[1][0] - pts[0][0], eps);
    if(curDir == 0)
    {
      curDir = axom::utilities::sign_of(pts[2][0] - pts[1][0], eps);
    }

    // Detect where z changes direction, and note those indices.
    for(axom::IndexType i = 1; i < segCount; ++i)
    {
      int segDir = axom::utilities::sign_of(pts[i + 1][0] - pts[i][0], eps);
      if(segDir == 0)
      {
        // Radial segment may or may not indicate change. Decide with next segment.
        continue;
      }
      if(segDir != curDir)
      {
        // Direction change
        int prevSegDir = axom::utilities::sign_of(pts[i][0] - pts[i - 1][0], eps);
        if(prevSegDir != 0)
        {
          // Case 1, a clear turn not involving a radial segment.
          boundaryIdx.push_back(i);
        }
        else
        {
          // Case 2, involving a radial segment.
          // Use the radially-closer point of the segment.
          int splitI = pts[i][1] < pts[i - 1][1] ? i : i - 1;
          boundaryIdx.push_back(splitI);
        }
        curDir = segDir;
        SLIC_ASSERT(curDir != 0);  // curDir ignores radial segments.
      }
    }
  }
  boundaryIdx.push_back(pts.size() - 1);
  return boundaryIdx;
}

void MonotonicZSORClipper::extractClipperInfo()
{
  auto sorOriginArray = m_info.fetch_existing("sorOrigin").as_double_array();
  auto sorDirectionArray = m_info.fetch_existing("sorDirection").as_double_array();
  for(int d = 0; d < 3; ++d)
  {
    m_sorOrigin[d] = sorOriginArray[d];
    m_sorDirection[d] = sorDirectionArray[d];
  }

  auto discreteFunctionArray = m_info.fetch_existing("discreteFunction").as_double_array();
  auto n = discreteFunctionArray.number_of_elements();

  SLIC_ERROR_IF(
    n % 2 != 0,
    axom::fmt::format(
      "***MonotonicZSORClipper: Discrete function must have an even number of values.  It has {}.",
      n));

  m_sorCurve.resize(axom::ArrayOptions::Uninitialized(), n / 2);
  for(int i = 0; i < n / 2; ++i)
  {
    m_sorCurve[i] = Point2DType {discreteFunctionArray[i * 2], discreteFunctionArray[i * 2 + 1]};
  }

  m_levelOfRefinement = m_info.fetch_existing("levelOfRefinement").to_double();
}

}  // namespace experimental
}  // end namespace quest
}  // end namespace axom
