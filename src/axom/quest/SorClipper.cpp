// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifdef AXOM_USE_CONDUIT
  #include "conduit_blueprint.hpp"
#endif

#include "axom/primal/operators/squared_distance.hpp"
#include "axom/quest/Discretize.hpp"
#include "axom/quest/SorClipper.hpp"
#include "axom/fmt.hpp"

#include <limits>

namespace axom
{
namespace quest
{

SorClipper::SorClipper(const klee::Geometry& kGeom, const std::string& name)
  : GeometryClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("Sor") : name)
  , m_maxRadius(0.0)
  , m_minRadius(std::numeric_limits<double>::max())
  , m_transformer()
{
  extractClipperInfo();

  // Combine internal and external rotations into m_transformer.
  m_transformer.addTranslation(m_sorOrigin.array());
  m_transformer.addRotation(Vector3DType({1,0,0}), m_sorDirection);
  m_transformer.addMatrix(m_extTrans);

  m_inverseTransformer = m_transformer.getInverse();

  const auto dCount = m_sorCurve.shape()[0];
  for(axom::IndexType i=0; i<dCount; ++i)
  {
    m_maxRadius = fmax(m_maxRadius, m_sorCurve(i, 1));
    m_minRadius = fmin(m_minRadius, m_sorCurve(i, 1));
  }

#ifdef AXOM_DEBUG
  axom::IndexType badZCount = 0;
  axom::IndexType badRadCount = m_sorCurve(0, 1) < 0;
  for(axom::IndexType i=1; i<dCount; ++i)
  {
    badZCount += m_sorCurve(i, 0) <= m_sorCurve(i-1, 0);
    badRadCount += m_sorCurve(i, 1) < 0;
  }
  SLIC_ERROR_IF(badZCount || badRadCount,
                axom::fmt::format("SorClipper '{}' has {} non-monotonically-increasing z-coordinates and {} negative radii", m_name, badZCount, badRadCount));
#endif


  computeRoughBlockings();
}

bool SorClipper::labelInOut(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  /*
    Planned implementation: (reverse) transform the mesh vertices to the
    r-z frame where the curve is defined as a 1D function.  It's easier
    to determine whether the point is in the sor that way.  I can check
    the SOR one conical section at a time.
  */
  AXOM_ANNOTATE_SCOPE("SorClipper::labelInOut");
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
void SorClipper::labelInOutImpl(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  SLIC_ERROR_IF(shapeeMesh.dimension() != 3, "SorClipper requires a 3D mesh.");

  int allocId = shapeeMesh.getAllocatorID();
  auto cellCount = shapeeMesh.getCellCount();

  auto inverseTransformer = m_inverseTransformer;
  auto transformer = m_transformer;

  axom::Array<BoundingBox2DType> bbOn(m_bbOn, allocId);
  axom::Array<BoundingBox2DType> bbUnder(m_bbUnder, allocId);
  auto bbOnView = bbOn.view();
  auto bbUnderView = bbUnder.view();

  auto cellLengths = shapeeMesh.getCellLengths();
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
    m_bbIn and m_bbsUnder.
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
        Point3DType pt3D({(vi & 1) ? boxMax[0] : boxMin[0],
                          (vi & 2) ? boxMax[1] : boxMin[1],
                          (vi & 4) ? boxMax[2] : boxMin[2]});
        Point3DType ptIt = inverseTransformer.getTransformed(pt3D);
        Point2DType ptRz({ptIt[0], sqrt(ptIt[1]*ptIt[1] + ptIt[2]*ptIt[2])});
        cellBbInRz.addPoint(ptRz);
      }
    });

  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
      const auto& cellBb = cellBbsInRzView[cellId];

      /*
        - If cellBb is close to any m_bbOn, label it ON.
        - Else if cellBb touches any m_bbUnder, label it IN.
          It cannot possibly be partially outside, because it
          doesn't cross the boundary or even touch any m_bbOn.
        - Else, label cellBb OUT.
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

  return;
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
void SorClipper::computeRoughBlockings()
{
  const axom::IndexType topI = m_sorCurve.shape()[0] - 1;
  const Point2DType basePt({m_sorCurve(0, 0), m_sorCurve(0, 1)});
  const Point2DType topPt({m_sorCurve(topI, 0), m_sorCurve(topI, 1)});

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
bool SorClipper::getGeometryAsOcts(quest::ShapeeMesh& shapeeMesh, axom::Array<OctahedronType>& octs)
{
  AXOM_ANNOTATE_BEGIN("SorClipper::getGeometryAsOcts");
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
  axom::ArrayView<Point2D> polyline((Point2D*)m_sorCurve.data(), m_sorCurve.shape()[0]);
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

  AXOM_ANNOTATE_END("SorClipper::getGeometryAsOcts");
  return true;
}

void SorClipper::extractClipperInfo()
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
      "***SorClipper: Discrete function must have an even number of values.  It has {}.",
      n));

  m_sorCurve.resize(axom::ArrayOptions::Uninitialized(), n / 2, 2);
  for(int i = 0; i < n; ++i)
  {
    m_sorCurve.flatIndex(i) = discreteFunctionArray[i];
  }

  m_levelOfRefinement = m_info.fetch_existing("levelOfRefinement").to_double();
}

}  // end namespace quest
}  // end namespace axom
