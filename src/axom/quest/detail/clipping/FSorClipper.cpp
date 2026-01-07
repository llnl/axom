// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/core/numerics/matvecops.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/primal/operators/squared_distance.hpp"
#include "axom/quest/Discretize.hpp"
#include "axom/quest/detail/clipping/FSorClipper.hpp"
#include "axom/fmt.hpp"

#include <limits>

namespace axom
{
namespace quest
{
namespace experimental
{

FSorClipper::FSorClipper(const klee::Geometry& kGeom, const std::string& name)
  : MeshClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("FSor") : name)
  , m_maxRadius(0.0)
  , m_minRadius(numerics::floating_point_limits<double>::max())
  , m_transformer()
{
  extractClipperInfo();

  combineRadialSegments(m_sorCurve);
  axom::Array<axom::IndexType> turnIndices = findZSwitchbacks(m_sorCurve.view());
  if(turnIndices.size() > 2)
  {
    // The 2 "turns" allowed are the first and last points.  Anything else is a switchback.
    SLIC_ERROR(
      "FSorClipper does not work when a curve doubles back"
      " in the axial direction.  Use SorClipper instead.");
  }

  for(auto& pt : m_sorCurve)
  {
    m_maxRadius = fmax(m_maxRadius, pt[1]);
    m_minRadius = fmin(m_minRadius, pt[1]);
  }
  SLIC_ERROR_IF(m_minRadius < 0.0,
                axom::fmt::format("FSorClipper '{}' has a negative radius", m_name));

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

FSorClipper::FSorClipper(const klee::Geometry& kGeom,
                         const std::string& name,
                         axom::ArrayView<const Point2DType> sorCurve,
                         const Point3DType& sorOrigin,
                         const Vector3DType& sorDirection,
                         axom::IndexType levelOfRefinement)
  : MeshClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("FSor") : name)
  , m_sorCurve(sorCurve, axom::execution_space<axom::SEQ_EXEC>::allocatorID())
  , m_maxRadius(0.0)
  , m_minRadius(numerics::floating_point_limits<double>::max())
  , m_sorOrigin(sorOrigin)
  , m_sorDirection(sorDirection)
  , m_levelOfRefinement(levelOfRefinement)
  , m_transformer()
{
  combineRadialSegments(m_sorCurve);
  axom::Array<axom::IndexType> turnIndices = findZSwitchbacks(m_sorCurve.view());
  if(turnIndices.size() > 2)
  {
    // The 2 "turns" allowed are the first and last points.  Anything else is a switchback.
    SLIC_ERROR(
      "FSorClipper does not work when a curve doubles back"
      " in the axial direction.  Use SorClipper instead.");
  }

  for(auto& pt : m_sorCurve)
  {
    m_maxRadius = fmax(m_maxRadius, pt[1]);
    m_minRadius = fmin(m_minRadius, pt[1]);
  }
  SLIC_ERROR_IF(m_minRadius < 0.0,
                axom::fmt::format("FSorClipper '{}' has a negative radius", m_name));

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

bool FSorClipper::labelCellsInOut(quest::experimental::ShapeMesh& shapeMesh,
                                  axom::Array<LabelType>& labels)
{
  SLIC_ERROR_IF(shapeMesh.dimension() != 3, "FSorClipper requires a 3D mesh.");

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

bool FSorClipper::labelTetsInOut(quest::experimental::ShapeMesh& shapeMesh,
                                 axom::ArrayView<const axom::IndexType> cellIds,
                                 axom::Array<LabelType>& tetLabels)
{
  SLIC_ERROR_IF(shapeMesh.dimension() != 3, "FSorClipper requires a 3D mesh.");

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

/*
 * Implementation: (reverse) transform the mesh vertices to the r-z
 * frame where the curve is defined as a r(z) function.  It's easier to
 * determine whether the point is in the sor that way.
 */
template <typename ExecSpace>
void FSorClipper::labelCellsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                                      axom::ArrayView<LabelType> labels)
{
  axom::Array<BoundingBox2DType> bbOn;
  axom::Array<BoundingBox2DType> bbUnder;
  computeCurveBoxes<ExecSpace>(shapeMesh, bbOn, bbUnder);
  const axom::ArrayView<const BoundingBox2DType> bbOnView = bbOn.view();
  const axom::ArrayView<const BoundingBox2DType> bbUnderView = bbUnder.view();

  const auto cellCount = shapeMesh.getCellCount();
  auto meshHexes = shapeMesh.getCellsAsHexes();
  auto meshCellVolumes = shapeMesh.getCellVolumes();
  auto invTransformer = m_invTransformer;
  constexpr double EPS = 1e-10;

  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
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
}

template <typename ExecSpace>
void FSorClipper::labelTetsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                                     axom::ArrayView<const axom::IndexType> cellIds,
                                     axom::ArrayView<LabelType> labels)
{
  axom::Array<BoundingBox2DType> bbOn;
  axom::Array<BoundingBox2DType> bbUnder;
  computeCurveBoxes<ExecSpace>(shapeMesh, bbOn, bbUnder);
  const axom::ArrayView<const BoundingBox2DType> bbOnView = bbOn.view();
  const axom::ArrayView<const BoundingBox2DType> bbUnderView = bbUnder.view();

  const auto cellCount = cellIds.size();
  auto meshHexes = shapeMesh.getCellsAsHexes();
  auto tetVolumes = shapeMesh.getTetVolumes();
  auto invTransformer = m_invTransformer;
  constexpr double EPS = 1e-10;

  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType ci) {
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
AXOM_HOST_DEVICE FSorClipper::BoundingBox2DType FSorClipper::estimateBoundingBoxInRz(
  const PolyhedronType& vertices)
{
  FSorClipper::BoundingBox2DType bbInRz;

  // Range of vertex angles in cylindrical coordinates.
  double minAngle = numerics::floating_point_limits<double>::max();
  double maxAngle = -numerics::floating_point_limits<double>::max();

  for(IndexType vi = 0; vi < vertices.numVertices(); ++vi)
  {
    auto& vert = vertices[vi];
    Point2DType vertOnXPlane {vert[1], vert[2]};
    Point2DType vertOnRz {
      vert[0],
      std::sqrt(numerics::dot_product(vertOnXPlane.data(), vertOnXPlane.data(), 2))};
    bbInRz.addPoint(vertOnRz);

    double angle = atan2(vertOnXPlane[1], vertOnXPlane[0]);
    minAngle = std::min(minAngle, angle);
    maxAngle = std::max(maxAngle, angle);
  }
#if 1
  /*
    The geometry can be closer to the rotation axis than its
    individual vertices are, depending on the angle (about the axis)
    between the vertices.  Given the angle, scale the bottom of bbInRz
    for the worst case.
  */
  double angleRange = maxAngle - minAngle;
  double factor = angleRange > M_PI ? 0.0 : cos(angleRange / 2);
  auto newMin = bbInRz.getMin();
  newMin[1] *= factor;
  bbInRz.addPoint(newMin);
#endif
#if 0
  /*
    Jeff's method to account for the angle.  Faster, but less
    discriminating, I think.
  */
  auto newMin = bbInRz.getMin();
  newMin[1] -= 0.5 * cellLength;
  if (newMin[1] < 0.0) newMin[1] = 0.0;
  bbInRz.addPoint(newMin);
#endif
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
AXOM_HOST_DEVICE inline MeshClipperStrategy::LabelType FSorClipper::rzBbToLabel(
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
void FSorClipper::computeCurveBoxes(quest::experimental::ShapeMesh& shapeMesh,
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

  using ReducePolicy = typename axom::execution_space<ExecSpace>::reduce_policy;
  using LoopPolicy = typename execution_space<ExecSpace>::loop_policy;
  RAJA::ReduceSum<ReducePolicy, double> sumCharLength(0.0);
  RAJA::forall<LoopPolicy>(
    RAJA::RangeSegment(0, cellCount),
    AXOM_LAMBDA(axom::IndexType cellId) { sumCharLength += cellLengths[cellId]; });
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

  axom::for_all<ExecSpace>(
    segCount,
    AXOM_LAMBDA(axom::IndexType i) {
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
Array<FSorClipper::Point2DType> FSorClipper::subdivideCurve(const Array<Point2DType>& sorCurveIn,
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

bool FSorClipper::getGeometryAsOcts(quest::experimental::ShapeMesh& shapeMesh,
                                    axom::Array<OctahedronType>& octs)
{
  AXOM_ANNOTATE_SCOPE("FSorClipper::getGeometryAsOcts");
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
bool FSorClipper::getGeometryAsOctsImpl(quest::experimental::ShapeMesh& shapeMesh,
                                        axom::Array<OctahedronType>& octs)
{
  const int allocId = shapeMesh.getAllocatorID();
  octs = axom::Array<OctahedronType>(0, 0, allocId);

  const auto cellCount = shapeMesh.getCellCount();

  // Compute an average characteristic length for the mesh cells.
  using ReducePolicy = typename axom::execution_space<ExecSpace>::reduce_policy;
  using LoopPolicy = typename execution_space<ExecSpace>::loop_policy;
#if 1
  axom::ArrayView<const double> cellVolumes = shapeMesh.getCellVolumes();
  RAJA::ReduceSum<ReducePolicy, double> sumVolume(0.0);
  RAJA::forall<LoopPolicy>(
    RAJA::RangeSegment(0, cellCount),
    AXOM_LAMBDA(axom::IndexType cellId) { sumVolume += cellVolumes[cellId]; });
  double avgVolume = sumVolume.get() / cellCount;
  double avgCharLength = pow(avgVolume, 1. / 3);
#else
  axom::ArrayView<const double> cellLengths = shapeMesh.getCellLengths();
  RAJA::ReduceSum<ReducePolicy, double> sumCharLength(0.0);
  RAJA::forall<LoopPolicy>(
    RAJA::RangeSegment(0, cellCount),
    AXOM_LAMBDA(axom::IndexType cellId) { sumCharLength += cellLengths[cellId]; });
  double avgCharLength = sumCharLength.get() / cellCount;
#endif

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
  axom::for_all<ExecSpace>(
    octCount,
    AXOM_LAMBDA(axom::IndexType iOct) {
      OctahedronType& oct = octsView[iOct];
      for(int iVert = 0; iVert < OctahedronType::NUM_VERTS; ++iVert)
      {
        transformer.transform(oct[iVert].array());
      }
    });

  SLIC_INFO(axom::fmt::format(
    "FSorClipper '{}' {}-level refinement got {} geometry octs from {} curve points.",
    name(),
    m_levelOfRefinement,
    octs.size(),
    sorCurve.size()));

  return true;
}

/*
  Combine consecutive radial segments in SOR curve.  Change in place.
*/
void FSorClipper::combineRadialSegments(axom::Array<Point2DType>& sorCurve)
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
axom::Array<axom::IndexType> FSorClipper::findZSwitchbacks(axom::ArrayView<const Point2DType> pts)
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

void FSorClipper::extractClipperInfo()
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
      "***FSorClipper: Discrete function must have an even number of values.  It has {}.",
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
