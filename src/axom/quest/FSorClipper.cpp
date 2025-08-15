// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/core/utilities/Utilities.hpp"
#include "axom/primal/operators/squared_distance.hpp"
#include "axom/quest/Discretize.hpp"
#include "axom/quest/FSorClipper.hpp"
#include "axom/fmt.hpp"

#include <limits>

namespace axom
{
namespace quest
{

FSorClipper::FSorClipper(const klee::Geometry& kGeom, const std::string& name)
  : GeometryClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("FSor") : name)
  , m_maxRadius(0.0)
  , m_minRadius(std::numeric_limits<double>::max())
  , m_transformer()
{
  extractClipperInfo();

  combineRadialSegments(m_sorCurve);
  axom::Array<axom::IndexType> turnIndices =
    findZSwitchbacks(m_sorCurve.view());
  if(turnIndices.size() > 2)
  {
    SLIC_ERROR("FSorClipper does not work when a curve doubles back"
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
  m_transformer.addRotation(Vector3DType({1,0,0}), m_sorDirection);
  m_transformer.addTranslation(m_sorOrigin.array());
  m_transformer.addMatrix(m_extTrans);
  m_inverseTransformer = m_transformer.getInverse();

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
  : GeometryClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("FSor") : name)
  , m_sorCurve(sorCurve, axom::execution_space<axom::SEQ_EXEC>::allocatorID())
  , m_maxRadius(0.0)
  , m_minRadius(std::numeric_limits<double>::max())
  , m_sorOrigin(sorOrigin)
  , m_sorDirection(sorDirection)
  , m_levelOfRefinement(levelOfRefinement)
  , m_transformer()
{
  combineRadialSegments(m_sorCurve);
  axom::Array<axom::IndexType> turnIndices =
    findZSwitchbacks(m_sorCurve.view());
  if(turnIndices.size() > 2)
  {
    SLIC_ERROR("FSorClipper does not work when a curve doubles back"
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
  m_transformer.addRotation(Vector3DType({1,0,0}), m_sorDirection);
  m_transformer.addTranslation(m_sorOrigin.array());
  m_transformer.addMatrix(m_extTrans);
  m_inverseTransformer = m_transformer.getInverse();

  for(const auto& pt : m_sorCurve)
  {
    m_curveBb.addPoint(pt);
  }
}

bool FSorClipper::labelInOut(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  AXOM_ANNOTATE_SCOPE("FSorClipper::labelInOut");
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
  return true;
}

/*
  Implementation: (reverse) transform the mesh vertices to the r-z
  frame where the curve is defined as a r(z) function.  It's easier to
  determine whether the point is in the sor that way.
*/
template <typename ExecSpace>
void FSorClipper::labelInOutImpl(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  SLIC_ERROR_IF(shapeeMesh.dimension() != 3, "FSorClipper requires a 3D mesh.");

  const int allocId = shapeeMesh.getAllocatorID();
  const auto cellCount = shapeeMesh.getCellCount();

  auto inverseTransformer = m_inverseTransformer;
  auto transformer = m_transformer;

  axom::ArrayView<const double> cellLengths = shapeeMesh.getCellLengths();

  using ReducePolicy = typename axom::execution_space<ExecSpace>::reduce_policy;
  using LoopPolicy = typename execution_space<ExecSpace>::loop_policy;
  RAJA::ReduceSum<ReducePolicy, double> sumCharLength(0.0);
  RAJA::forall<LoopPolicy>(
    RAJA::RangeSegment(0, cellCount),
      AXOM_LAMBDA(axom::IndexType cellId) {
      sumCharLength += cellLengths[cellId];
    });
  double avgCharLength = sumCharLength.get()/cellCount;

  // Subdivide the SOR curve and place in the right allocator.
  axom::Array<Point2DType> sorCurve = subdivideCurve(m_sorCurve, avgCharLength);
  sorCurve = axom::Array<Point2DType>(sorCurve, shapeeMesh.getAllocatorID());
  auto sorCurveView = sorCurve.view();

  /*
   * Compute bounding boxes bbOn, which cover the curve segments, and
   * bbUnder, which cover the space between bbOn and the z axis.  bbOn
   * includes end caps, the segments that join the curve to the
   * z-axis.
   */
  auto segCount = sorCurve.size() - 1;
  axom::Array<BoundingBox2DType> bbOn(segCount + 2,
                                      segCount + 2,
                                      allocId);
  axom::Array<BoundingBox2DType> bbUnder(segCount,
                                         segCount,
                                         allocId);
  auto bbOnView = bbOn.view();
  auto bbUnderView = bbUnder.view();
  axom::for_all<ExecSpace>(
    segCount,
    AXOM_LAMBDA(axom::IndexType i) {
      BoundingBox2DType& on = bbOnView[i];
      BoundingBox2DType& under = bbUnderView[i];
      // on = BoundingBox2DType{&sorCurve[i], 2};
      on.addPoint(sorCurveView[i]);
      on.addPoint(sorCurveView[i+1]);
      Point2DType underMin{on.getMin()[0], 0.0};
      Point2DType underMax{on.getMax()[0], on.getMin()[1]};
      under = BoundingBox2DType(underMin, underMax);
    });

  axom::Array<BoundingBox2DType> endCaps(2, 2);
  endCaps[0].addPoint(m_sorCurve.front());
  endCaps[0].addPoint(Point2DType{m_sorCurve.front()[0], 0.0});
  endCaps[1].addPoint(m_sorCurve.back());
  endCaps[1].addPoint(Point2DType{m_sorCurve.back()[0], 0.0});
  axom::copy(&bbOn[segCount], endCaps.data(), endCaps.size()*sizeof(BoundingBox2DType));

  const double lenFactor = 0.5;

  auto cellBbs = shapeeMesh.getCellBoundingBoxes();
  constexpr int NUM_BB_VERTS = 8;

  if(labels.size() < cellCount || labels.getAllocatorID() != shapeeMesh.getAllocatorID())
  {
    labels = axom::Array<LabelType>(ArrayOptions::Uninitialized(), cellCount, cellCount, allocId);
  }

  auto labelsView = labels.view();

  /*
    Compute cell bounding boxes in rz plane, to be checked against
    bbOn and bbUnder.
  */
  axom::Array<BoundingBox2DType> cellBbsInRz(cellBbs.size(),
                                             cellBbs.size(),
                                             shapeeMesh.getAllocatorID());
  auto cellBbsInRzView = cellBbsInRz.view();
  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
      const auto& cellBb = cellBbs[cellId];
      const auto& boxMin = cellBb.getMin();
      const auto& boxMax = cellBb.getMax();
      BoundingBox2DType& cellBbInRz = cellBbsInRzView[cellId];
      for(int vi = 0; vi < NUM_BB_VERTS; ++vi)
      {
        Point3DType vert3D({(vi & 1) ? boxMax[0] : boxMin[0],
                            (vi & 2) ? boxMax[1] : boxMin[1],
                            (vi & 4) ? boxMax[2] : boxMin[2]});
        Point3DType vert3Dt = inverseTransformer.getTransformed(vert3D);
        Point2DType vertRz({vert3Dt[0], sqrt(vert3Dt[1]*vert3Dt[1] + vert3Dt[2]*vert3Dt[2])});
        cellBbInRz.addPoint(vertRz);
      }
    });

  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
      const auto& cellBb = cellBbsInRzView[cellId];

      /*
        - If cellBb is close to any bbOn, label it ON.
        - Else if cellBb touches any bbUnder, label it IN.
          It cannot possibly be partially outside, because it
          doesn't cross the boundary or even touch any bbOn.
        - Else, label cellBb OUT.

        We expect bbOn and bbUnder to be small arrays, so we use
        linear searches.  If that's too slow, we can use a BVH.
      */

      auto& cellLabel = labelsView[cellId];
      double sqDistThreshold = lenFactor * cellLengths[cellId];
      sqDistThreshold = sqDistThreshold*sqDistThreshold;
      for(const auto& bbOn : bbOnView)
      {
        double sqDist = axom::primal::squared_distance(cellBb, bbOn);
        if (sqDist <= sqDistThreshold)
        {
          cellLabel = LABEL_ON;
          return;
        }
      }

      for(const auto& bbUnder : bbUnderView)
      {
        if(cellBb.intersectsWith(bbUnder))
        {
          cellLabel = LABEL_IN;
          return;
        }
      }

      cellLabel = LABEL_OUT;
    });
}

/*
  Reduce large segments, which have bounding boxes that
  overlap much more than the segments actually overlap.

  Use harmonic mean because it works better than dz to prevent
  dividing segments that are cylindrical or nearly cylindrical.
*/
Array<FSorClipper::Point2DType> FSorClipper::subdivideCurve(
  const Array<Point2DType>& sorCurveIn,
  double cellCharacteristicLength)
{
  Array<Point2DType> sorCurveOut;

  if(sorCurveIn.empty())
  {
    return sorCurveOut;
  }

  const double maxMean = 3 * cellCharacteristicLength;
  const double minDz = 2 * cellCharacteristicLength;

  // Reserve estimated total number of points needed
  sorCurveOut.reserve(sorCurveIn.size() * 1.2 + 10);
  sorCurveOut.push_back(sorCurveIn[0]);

  for(IndexType i = 1; i < sorCurveIn.size(); ++i)
  {
    const Point2DType& segStart = sorCurveIn[i-1];
    const Point2DType& segEnd = sorCurveIn[i];

    const auto delta = segEnd.array() - segStart.array();
    const auto absDelta = axom::abs(delta);
    const double segDz = absDelta[0];
    const double segDr = absDelta[1];

    // To drive harmonic mean below maxMean.
    double segMean = 2 * segDz * segDr / (segDz + segDr);
    int numSplitsByMean = static_cast<int>(std::ceil(segMean / maxMean));

    // To prevent dz from falling below minDz
    int numSplitsByMinDz = static_cast<int>(std::ceil(segDz / minDz));

    int numSplits = std::min(numSplitsByMean, numSplitsByMinDz);

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

bool FSorClipper::getGeometryAsOcts(quest::ShapeeMesh& shapeeMesh, axom::Array<OctahedronType>& octs)
{
  AXOM_ANNOTATE_SCOPE("FSorClipper::getGeometryAsOcts");
  switch(shapeeMesh.getRuntimePolicy())
  {
  case axom::runtime_policy::Policy::seq:
    getGeometryAsOctsImpl<axom::SEQ_EXEC>(shapeeMesh, octs);
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    getGeometryAsOctsImpl<axom::OMP_EXEC>(shapeeMesh, octs);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case axom::runtime_policy::Policy::cuda:
    getGeometryAsOctsImpl<axom::CUDA_EXEC<256>>(shapeeMesh, octs);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    getGeometryAsOctsImpl<axom::HIP_EXEC<256>>(shapeeMesh, octs);
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
  return true;
}

/*
  Compute octahedral geometry representation, with an execution policy.

  Side effect: m_sorCurve data is reallocated to the shapeeMesh allocator,
  if it's not there yet.
*/
template<typename ExecSpace>
bool FSorClipper::getGeometryAsOctsImpl(quest::ShapeeMesh& shapeeMesh, axom::Array<OctahedronType>& octs)
{
  const int allocId = shapeeMesh.getAllocatorID();
  octs = axom::Array<OctahedronType>(0, 0, allocId);

  // Generate the Octahedra
  int octCount = 0;
  const bool good = axom::quest::discretize<ExecSpace>(m_sorCurve.view(),
                                                       int(m_sorCurve.size()),
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
      for(int iVert = 0; iVert < OctType::NUM_VERTS; ++iVert)
      {
        transformer.transform(oct[iVert].array());
      }
    });

  return true;
}

/*
  Combine consecutive radial segments in SOR curve.  Change in place.
*/
void FSorClipper::combineRadialSegments(axom::Array<Point2DType>& sorCurve)
{
  int ptCount = sorCurve.size();
  if(ptCount < 3) { return; }

  constexpr double eps = 1e-14;

  // Set sorCurve[j] to sorCurve[i] where j <= i, skipping points
  // joining consecutive radial segments.

  int j = 1;
  bool prevIsRadial = axom::utilities::isNearlyEqual(sorCurve[j][0] - sorCurve[j-1][0], eps);
  bool curIsRadial = false;
  for (int i = 2; i < ptCount; ++i)
  {
    curIsRadial = axom::utilities::isNearlyEqual(sorCurve[i][0] - sorCurve[i-1][0], eps);
    /*
      Current and previous segments share point j.  If both are
      consecutive radial segments, discard point j by overwriting it
      with point i.  Else, copy point i to a new point j.
    */
    if (!(curIsRadial && prevIsRadial)) { ++j; }
    sorCurve[j] = sorCurve[i];
    prevIsRadial = curIsRadial;
  }
  sorCurve.resize(j + 1);
}

/*
  Find points along the r-z curve where the z-coordinate changes direction.

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
axom::Array<axom::IndexType> FSorClipper::findZSwitchbacks(
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
    if (curDir == 0) { curDir = axom::utilities::sign_of(pts[2][0] - pts[1][0], eps); }

    // Detect where z changes direction, and note those indices.
    for(axom::IndexType i = 1; i < segCount; ++i)
    {
      int segDir = axom::utilities::sign_of(pts[i+1][0] - pts[i][0], eps);
      if (segDir == 0) {
        // Radial segment may or may not indicate change. Decide with next segment.
        continue;
      }
      if (segDir != curDir)
      {
        // Direction change
        int prevSegDir = axom::utilities::sign_of(pts[i][0] - pts[i-1][0], eps);
        if (prevSegDir != 0)
        {
          // Case 1, a clear turn not involving a radial segment.
          boundaryIdx.push_back(i);
        }
        else
        {
          // Case 2, involving a radial segment.
          // Use the radially-closer point of the segment.
          int splitI = pts[i][1] < pts[i-1][1] ? i : i - 1;
          boundaryIdx.push_back(splitI);
        }
        curDir = segDir;
        SLIC_ASSERT( curDir != 0 ); // curDir ignores radial segments.
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
  for(int i = 0; i < n/2; ++i)
  {
    m_sorCurve[i] = Point2DType{discreteFunctionArray[i*2], discreteFunctionArray[i*2 + 1]};
  }

  m_levelOfRefinement = m_info.fetch_existing("levelOfRefinement").to_double();
}

}  // end namespace quest
}  // end namespace axom
