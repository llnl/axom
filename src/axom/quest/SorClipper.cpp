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

  const auto dCount = m_discreteFcn.shape()[0];
  for(axom::IndexType i=0; i<dCount; ++i)
  {
    m_maxRadius = fmax(m_maxRadius, m_discreteFcn(i, 1));
    m_minRadius = fmin(m_minRadius, m_discreteFcn(i, 1));
  }

#ifdef AXOM_DEBUG
  axom::IndexType badZCount = 0;
  axom::IndexType badRadCount = m_discreteFcn(0, 1) < 0;
  for(axom::IndexType i=1; i<dCount; ++i)
  {
    badZCount += m_discreteFcn(i, 0) <= m_discreteFcn(i-1, 0);
    badRadCount += m_discreteFcn(i, 1) < 0;
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

  /*
    Compute cell bounding boxes in rz plane, to be checked against
    m_bbIn and mbbsUnder.
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

  auto labelsView = labels.view();

  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
      const auto& cellBb = cellBbsInRz[cellId];

      // If cellBb is close to any m_bbOn, label it ON.
      // Else if cellBb touches any m_bbUnder, label it IN.
      // Else, label cellBb OUT.

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

#if 0
  Old stuff
  /*
    Generate a stack of cones to represent the SOR.
  */
  axom::Array<Cone3DType> cones(m_discreteFcn.size() - 1,
                                m_discreteFcn.size() - 1,
                                allocId);
  auto conesView = cones.view();

  axom::Array<double, 2> discreteFcn(m_discreteFcn, allocId);
  auto discreteFcnView = discreteFcn.view();

  axom::for_all<ExecSpace>(
    cones.size(),
    AXOM_LAMBDA(axom::IndexType coneId) {
      Point3DType basePt({discreteFcnView(coneId, 0), 0.0, 0.0});
      Point3DType topPt({discreteFcnView(coneId + 1, 0), 0.0, 0.0});
      basePt = transformer.getTransformed(basePt);
      conesView[coneId] =
        Cone3DType(discreteFcnView(coneId, 0), discreteFcnView(coneId, 1),
                   discreteFcnView(coneId + 1, 0), discreteFcnView(coneId, 1),
                   sorDirection,
                   sorOrigin);
    });

  /*
    Compute whether vertices are inside shape.
  */
  axom::Array<bool> vertIsInside {ArrayOptions::Uninitialized(), vertCount, vertCount, allocId};
  auto vertIsInsideView = vertIsInside.view();
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vX.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vY.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vZ.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vertIsInsideView.getAllocatorID()));

  double maxRadius = m_maxRadius;
  axom::for_all<ExecSpace>(
    vertCount,
    AXOM_LAMBDA(axom::IndexType vertId) {
      primal::Point3D vert {vX[vertId], vY[vertId], vZ[vertId]};
      inverseTransformer.transform(vert.array());
      if (vert[0] < conesView[0].getBaseZ() || vert[0] > conesView[conesView.size()-1].getBaseZ())
      {
        vertIsInsideView[vertId] = false;
      }
      else
      {
        double vertR = sqrt(vert[1]*vert[1] + vert[2]*vert[2]);
        if ( vertR > maxRadius )
        {
          vertIsInsideView[vertId] = false;
        }
        else
        {
          // Determine which cone's section vert is in and whether it is
          // under the rz curve.
          double vertZ = vert[0];
          axom::IndexType upperConeId = conesView.size() - 1;
          axom::IndexType lowerConeId = 0;
          while (upperConeId != lowerConeId) {
            axom::IndexType coneId = (upperConeId + lowerConeId)/2;
            if (vertZ < conesView[coneId].getBaseZ())
            {
              upperConeId = coneId;
            }
            else if (vertZ > conesView[coneId].getTopZ())
            {
              lowerConeId = coneId;
            }
          }
          SLIC_ASSERT(lowerConeId == upperConeId);
          const auto& cone = conesView[lowerConeId];
          double coneRadiusAtVertZ = cone.getRadiusAt(vertZ);
          vertIsInsideView[vertId] = vertR < coneRadiusAtVertZ;
        }
      }
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
#endif

  return;
}


/*
  Compute m_bbOn and m_bbUnder, the boxes that blocks of areas
  on and under the rz curve.  The blocking is rough but conservative.

  Currently, we have just 1 box under the curve, but we code for
  future arrays of boxes to block more efficiently.
  - m_bbOn has 3 boxes.  One for the base plane, one for the top
    plane and one for the m_discreteFcn curve.
  - m_bbUnder has only the boxes under m_discreteFcn.
*/
void SorClipper::computeRoughBlockings()
{
  const axom::IndexType topI = m_discreteFcn.shape()[0] - 1;
  const Point2DType basePt({m_discreteFcn(0, 0), m_discreteFcn(0, 1)});
  const Point2DType topPt({m_discreteFcn(topI, 0), m_discreteFcn(topI, 1)});

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
  axom::ArrayView<Point2D> polyline((Point2D*)m_discreteFcn.data(), m_discreteFcn.shape()[0]);
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

  m_discreteFcn.resize(axom::ArrayOptions::Uninitialized(), n / 2, 2);
  for(int i = 0; i < n; ++i)
  {
    m_discreteFcn.flatIndex(i) = discreteFunctionArray[i];
  }

  m_levelOfRefinement = m_info.fetch_existing("levelOfRefinement").to_double();
}

}  // end namespace quest
}  // end namespace axom
