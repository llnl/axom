// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifdef AXOM_USE_CONDUIT
  #include "conduit_blueprint.hpp"
#endif

#include "axom/core/utilities/Utilities.hpp"
#include "axom/primal/operators/squared_distance.hpp"
#include "axom/quest/Discretize.hpp"
#include "axom/quest/FSorClipper.hpp"
#include "axom/spin/BVH.hpp"
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

  for(auto& pt : m_sorCurve)
  {
    m_maxRadius = fmax(m_maxRadius, pt[1]);
    m_minRadius = fmin(m_minRadius, pt[1]);
  }
  SLIC_ERROR_IF(m_minRadius < 0.0,
                axom::fmt::format("FSorClipper '{}' has a negative radius", m_name));

  if(!pointsAreAxiallyMonotonic(m_sorCurve.view()))
  {
    SLIC_ERROR("FSorClipper does not work when a curve doubles back in the axial direction.  Use SorClipper instead.");
  }

  // Combine internal and external rotations into m_transformer.
  m_transformer.addRotation(Vector3DType({1,0,0}), m_sorDirection);
  m_transformer.addTranslation(m_sorOrigin.array());
  m_transformer.addMatrix(m_extTrans);
  m_inverseTransformer = m_transformer.getInverse();

  for(const auto& pt : m_sorCurve)
  {
    m_curveBb.addPoint(pt);
  }

  clusterSorFunction();
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
  for(auto& pt : m_sorCurve)
  {
    m_maxRadius = fmax(m_maxRadius, pt[1]);
    m_minRadius = fmin(m_minRadius, pt[1]);
  }
  SLIC_ERROR_IF(m_minRadius < 0.0,
                axom::fmt::format("FSorClipper '{}' has a negative radius", m_name));

  if(!pointsAreAxiallyMonotonic(m_sorCurve.view()))
  {
    SLIC_ERROR("FSorClipper does not work when a curve doubles back in the axial direction.  Use SorClipper instead.");
  }

  // Combine internal and external rotations into m_transformer.
  m_transformer.addRotation(Vector3DType({1,0,0}), m_sorDirection);
  m_transformer.addTranslation(m_sorOrigin.array());
  m_transformer.addMatrix(m_extTrans);
  m_inverseTransformer = m_transformer.getInverse();

  for(const auto& pt : m_sorCurve)
  {
    m_curveBb.addPoint(pt);
  }

  clusterSorFunction();
}

bool FSorClipper::labelInOut(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  /*
    Planned implementation: (reverse) transform the mesh vertices to the
    r-z frame where the curve is defined as a 1D function.  It's easier
    to determine whether the point is in the sor that way.  I can check
    the SOR one conical section at a time.
  */
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

template <typename ExecSpace>
void FSorClipper::labelInOutImpl(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  SLIC_ERROR_IF(shapeeMesh.dimension() != 3, "FSorClipper requires a 3D mesh.");

  int allocId = shapeeMesh.getAllocatorID();
  auto cellCount = shapeeMesh.getCellCount();

  auto inverseTransformer = m_inverseTransformer;
  auto transformer = m_transformer;

  auto cellBbs = shapeeMesh.getCellBoundingBoxes();
  constexpr int NUM_BB_VERTS = 8;

  if(labels.size() < cellCount || labels.getAllocatorID() != shapeeMesh.getAllocatorID())
  {
    labels = axom::Array<LabelType>(ArrayOptions::Uninitialized(), cellCount, cellCount, allocId);
  }

  auto labelsView = labels.view();

  /*
    Compute cell bounding boxes in rz plane, to be checked against the
    2D curve.
  */

  axom::Array<BoundingBox2DType> cellBbsInRz(cellCount,
                                             cellCount,
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
        Point3DType pt3D({(vi & 1) ? boxMax[0] : boxMin[0],
                          (vi & 2) ? boxMax[1] : boxMin[1],
                          (vi & 4) ? boxMax[2] : boxMin[2]});
        Point3DType ptIt = inverseTransformer.getTransformed(pt3D);
        Point2DType ptRz({ptIt[0], sqrt(ptIt[1]*ptIt[1] + ptIt[2]*ptIt[2])});
        cellBbInRz.addPoint(ptRz);
      }
    });

  /*
    Compute function's segments and their bounding boxes.
    Construct BVH for those boxes.
  */
  axom::Array<Point2DType> sorCurve;
  axom::Array<Point2DType> sorCurveTmp;
  axom::ArrayView<Point2DType> sorCurveView =
    getCurveWithAxisPoints(sorCurveTmp) ? sorCurveTmp.view() : m_sorCurve.view();
  if(!axom::execution_space<ExecSpace>::usesAllocId(sorCurveView.getAllocatorID()))
  {
    sorCurve = axom::Array<Point2DType>(sorCurveView, allocId);
    sorCurveView = sorCurve.view();
  }

  const axom::IndexType segCount = sorCurveView.size() - 1;
  axom::Array<Segment2DType> segs(segCount, segCount, allocId);
  axom::Array<BoundingBox2DType> segBbs(segCount, segCount, allocId);
  auto segsView = segs.view();
  auto segBbsView = segBbs.view();
  axom::for_all<ExecSpace>(
    segCount,
    AXOM_LAMBDA(axom::IndexType i)
    {
      auto& seg = segsView[i];
      auto& segBb = segBbsView[i];
      seg = Segment2DType(sorCurveView[i], sorCurveView[i+1]);
      segBb.addPoint(seg[0]);
      segBb.addPoint(seg[1]);
    });

  spin::BVH<2, ExecSpace, double> bvh(segBbsView, segBbsView.size(), allocId);

  BoundingBox2DType curveBb = bvh.getBounds();

  /*
    Find cells intersecting the segment bounding boxes.
    These cells are ON the SOR boundary.
  */

  axom::Array<IndexType> offsets(cellCount, cellCount, allocId);
  axom::Array<IndexType> counts(cellCount, cellCount, allocId);
  axom::Array<IndexType> candidates;
  AXOM_ANNOTATE_BEGIN("bvh.findBoundingBoxes");
  bvh.findBoundingBoxes(offsets, counts, candidates, cellCount, cellBbsInRzView);
  AXOM_ANNOTATE_END("bvh.findBoundingBoxes");

  auto offsetsView = offsets.view();
  auto countsView = counts.view();
  auto candidatesView = candidates.view();

  axom::for_all<ExecSpace>(
    cellBbsInRzView.size(),
    AXOM_LAMBDA(axom::IndexType i)
    {
      auto& cellBbInRz = cellBbsInRzView[i];
      if (!cellBbInRz.intersectsWith(curveBb))
      {
        labelsView[i] = LABEL_OUT;
        return;
      }

      auto candidatesBegin = offsetsView[i];
      auto candidatesEnd = offsetsView[i] + countsView[i];
      for(axom::IndexType j=candidatesBegin; j<candidatesEnd; ++j)
      {
        axom::IndexType candidateIdx = candidatesView[j];
        const auto& candidateSeg = segsView[candidateIdx];
        if(axom::primal::intersect(candidateSeg, cellBbInRz))
        {
          labelsView[i] = LABEL_ON;
          return;
        }
      }
      /*
        Cell is definitely NOT ON boundary.  Count segments above
        the cell.  If it's even, cell is OUT, else it's IN.

        Use an exhaustive search through segs.
        Alternative if that's slow: Create rays in the r-direction,
        from centroids, and count number of intersecting segments.
        Use bvh.findRays to limit the search.  (This has to be done
        in a separate loop.)
      */
      auto cellCentroid = cellBbInRz.getCentroid();
      const Ray2DType radial(cellCentroid, {0.0, 1.0});
      int intersectionCount = 0;
      for(const auto& seg : segsView) {
        bool intersects = axom::primal::intersect(radial, seg);
        intersectionCount += intersects;
      }
      labelsView[i] = intersectionCount % 2 == 1 ? LABEL_IN : LABEL_OUT;
    });
}

/*
  Compute m_bbOn and m_bbUnder, the boxes that block of areas
  on and under the rz curve.  The blocking is rough but conservative.

  Currently, we have just 1 box under the curve, but we code for
  future arrays of boxes to be more discriminating.
  - m_bbOn has 3 boxes.  One for the base plane, one for the top
    plane and one for the m_sorCurve curve.
  - m_bbUnder has only the boxes under m_sorCurve.
*/
void FSorClipper::clusterSorFunction()
{
  const axom::IndexType topI = m_sorCurve.size() - 1;
  const Point2DType basePt = m_sorCurve[0];
  const Point2DType topPt = m_sorCurve[topI];

  m_bbOn.resize(3);
  m_bbOn[0].addPoint(Point2DType{basePt[0], m_minRadius});
  m_bbOn[0].addPoint(Point2DType{topPt[0], m_maxRadius});
  m_bbOn[1].addPoint(Point2DType{basePt[0], 0.0});
  m_bbOn[1].addPoint(basePt);
  m_bbOn[2].addPoint(Point2DType{topPt[0], 0.0});
  m_bbOn[2].addPoint(topPt);

  m_bbUnder.resize(m_bbOn.size() - 2);
  axom::for_all<axom::SEQ_EXEC>(m_bbUnder.size(),
                                AXOM_LAMBDA(axom::IndexType i) {
                                  // Set bbUnder to the box from the z-axis to the bottom of bbOn.
                                  auto minPt = m_bbOn[i].getMin();
                                  auto maxPt = m_bbOn[i].getMax();
                                  maxPt[1] = minPt[1];
                                  minPt[1] = 0.0;
                                  m_bbUnder[i] = BoundingBox2DType(minPt, maxPt, false);
                                });
}

// TODO: Factor out execution-space-specific computations into a template method,
// instead of computing on host and copying to GPU.
bool FSorClipper::getGeometryAsOcts(quest::ShapeeMesh& shapeeMesh, axom::Array<OctahedronType>& octs)
{
  AXOM_ANNOTATE_BEGIN("FSorClipper::getGeometryAsOcts");
  const int hostAllocId = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int allocId = shapeeMesh.getAllocatorID();

  if(octs.getAllocatorID() != allocId || octs.size() != 0)
  {
    octs = axom::Array<OctahedronType>(0, 0, allocId);
  }

  axom::Array<OctahedronType> tmpOcts(0, 0, hostAllocId);

  // Host loop writes directly to octs if possible, else a temporary array.
  axom::Array<OctahedronType>& octsOnHost =
    axom::execution_space<axom::SEQ_EXEC>::usesAllocId(octs.getAllocatorID()) ? octs : tmpOcts;

  if(&octsOnHost == &tmpOcts)
  {
    tmpOcts.resize(0, OctahedronType());
  }

  // Generate the Octahedra
  int octCount = 0;
  axom::ArrayView<Point2D> polyline = m_sorCurve.view();
  const bool good = axom::quest::discretize<axom::SEQ_EXEC>(polyline,
                                                            int(polyline.size()),
                                                            m_levelOfRefinement,
                                                            octsOnHost,
                                                            octCount);
  AXOM_UNUSED_VAR(good);
  SLIC_ASSERT(good);
  SLIC_ASSERT(octCount == octsOnHost.size());

  auto transformer = m_transformer;
  auto octsView = octsOnHost.view();
  axom::for_all<axom::SEQ_EXEC>(
    octCount,
    AXOM_LAMBDA(axom::IndexType iOct) {
      OctahedronType& oct = octsView[iOct];
      for(int iVert = 0; iVert < OctType::NUM_VERTS; ++iVert)
      {
        transformer.transform(oct[iVert].array());
      }
    });

  // The disretize method uses host data.  Place into proper space if needed.
  if(&octs != &octsOnHost)
  {
    octs.resize(octsOnHost.size());
    axom::copy(octs.data(), octsOnHost.data(), sizeof(OctahedronType) * octs.size());
  }

  AXOM_ANNOTATE_END("FSorClipper::getGeometryAsOcts");
  return true;
}

bool FSorClipper::pointsAreAxiallyMonotonic(const axom::ArrayView<Point2DType>& sorCurve)
{
  constexpr double eps = 1e-14;
  for(int i = 1; i < sorCurve.size()-1; ++i)
  {
    double prevDiff = (sorCurve[i][0] - sorCurve[i-1][0]);
    double nextDiff = (sorCurve[i+1][0] - sorCurve[i][0]);
    int prevSign = axom::utilities::sign_of(prevDiff, eps);
    int nextSign = axom::utilities::sign_of(nextDiff, eps);
    if (prevSign != nextSign && prevSign != 0 && nextSign != 0)
    {
      return false;
    }
  }
  return true;
}

bool FSorClipper::getCurveWithAxisPoints(axom::Array<Point2DType>& curveWithAxisPoints)
{
  /*
    The function is considered a loop if the first and last points are
    close to each other.  If not a loop, add a point to each end to
    bring the curve down to the axis of rotation.
  */
  const double eps = 1e-12;
  auto ptCount = m_sorCurve.size();
  Point2DType firstPt =  m_sorCurve.front();
  Point2DType lastPt = m_sorCurve.back();
  double sep = Segment2DType(firstPt, lastPt).length();
  bool isLoop = sep < eps;
  if (!isLoop) {
    bool addFirst = firstPt[1] > eps;
    bool addLast = lastPt[1] > eps;
    auto newPtCount = ptCount + addFirst + addLast;
    curveWithAxisPoints = axom::Array<Point2DType>(newPtCount, newPtCount);
    axom::copy((curveWithAxisPoints.data() + addFirst),
               m_sorCurve.data(),
               m_sorCurve.size() * sizeof(Point2DType));
    if(addFirst)
    {
      curveWithAxisPoints.front() = Point2DType({firstPt[0], 0.0});
    }
    if(addLast)
    {
      curveWithAxisPoints.back() = Point2DType({lastPt[0], 0.0});
    }
    return true;
  }
  return false;
}

/*
  Combine consecutive radial segments in SOR curve.
*/
void FSorClipper::combineRadialSegments(axom::Array<Point2DType>& sorCurve)
{
  int ptCount = sorCurve.size();
  if(ptCount < 3) { return; }

  constexpr double eps = 1e-14;

  // Compute in place.  Set sorCurve[j] to sorCurve[i] where
  // j <= i, skipping points joining consecutive radial segments.

  int j = 1;
  bool prevIsRadial = axom::utilities::isNearlyEqual(sorCurve[j][0] - sorCurve[j-1][0], eps);
  bool curIsRadial = false;
  for (int i = 2; i < ptCount; ++i)
  {
    curIsRadial = axom::utilities::isNearlyEqual(sorCurve[i][0] - sorCurve[i-1][0], eps);
    // If current and previous aren't consecutive radial segments,
    // copy pt to new point j.  Else, overwrite current point j.
    if (!(curIsRadial && prevIsRadial)) { ++j; }
    sorCurve[j] = sorCurve[i];
    prevIsRadial = curIsRadial;
  }
  sorCurve.resize(j + 1);
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
