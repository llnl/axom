// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifdef AXOM_USE_CONDUIT
  #include "conduit_blueprint.hpp"
#endif

#include "axom/mint/mesh/Mesh.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/quest/Discretize.hpp"
#include "axom/quest/SorClipper.hpp"
#include "axom/fmt.hpp"

namespace axom
{
namespace quest
{

SorClipper::SorClipper(const klee::Geometry& kGeom, const std::string& name)
  : GeometryClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("Sor") : name)
  , m_maxRadius(0.0)
  , m_transformer(m_transMat)
{
  extractClipperInfo();

  const auto dCount = m_discreteFcn.shape()[0];
  for(axom::IndexType i=0; i<dCount; ++i)
  {
    m_maxRadius = fmax(m_maxRadius, m_discreteFcn(i, 1));
  }

#ifdef AXOM_DEBUG
  int badRadCount = 0;
  for(axom::IndexType i=0; i<dCount; ++i)
  {
    badRadCount += m_discreteFcn(i, 1) < 0;
  }
  SLIC_ASSERT(badRadCount == 0);
#endif

  m_inverseTransformer = m_transformer.getInverse();
}

bool SorClipper::labelInOut(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  return false; // Work in progress.  Not checked.
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

  constexpr int NUM_VERTS_PER_CELL = 8;

  int allocId = shapeeMesh.getAllocatorID();
  auto cellCount = shapeeMesh.getCellCount();
  auto vertCount = shapeeMesh.getVertexCount();

  const auto& vertCoords = shapeeMesh.getVertexCoords3D();
  const auto& vX = vertCoords[0];
  const auto& vY = vertCoords[1];
  const auto& vZ = vertCoords[2];

  auto sorBase = m_sorBase;
  auto sorDirection = m_sorDirection;
  auto inverseTransformer = m_inverseTransformer;
  auto transformer = m_transformer;

  axom::Array<double, 2> discreteFcn(m_discreteFcn, allocId);
  auto discreteFcnView = discreteFcn.view();

  /*
    Generate a stack of cones to represent the SOR.
  */
  axom::Array<Cone3DType> cones(m_discreteFcn.size() - 1,
                                m_discreteFcn.size() - 1,
                                allocId);
  auto conesView = cones.view();

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
                   sorBase);
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

  return;
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

  // Rotate to the SOR axis direction and translate to the base location.
  // Then apply the external transformation.
  // TODO: Combine rotate, translate and m_transformer into a single Matrix
  //   and apply that.  It's faster.
  numerics::Matrix<double> rotate = sorAxisRotMatrix(m_sorDirection);
  const auto& translate = m_sorBase;
  auto octsView = octsOnHost.view();
  auto transformer = m_transformer;
  axom::for_all<axom::SEQ_EXEC>(
    octCount,
    AXOM_LAMBDA(axom::IndexType iOct) {
      OctahedronType& oct = octsView[iOct];
      for(int iVert = 0; iVert < OctType::NUM_VERTS; ++iVert)
      {
        Point3DType& newCoords = oct[iVert];
        Point3DType oldCoords = newCoords;
        numerics::matrix_vector_multiply(rotate, oldCoords.data(), newCoords.data());
        newCoords.array() += translate.array();
        transformer.transform(newCoords.array());
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

// Return a 3x3 matrix that rotates coordinates from the x-axis to the given direction.
numerics::Matrix<double> SorClipper::sorAxisRotMatrix(const Vector3DType& dir)
{
  // Note that the rotation matrix is not unique.
  static const Vector3DType x {1.0, 0.0, 0.0};
  Vector3DType a = dir.unitVector();
  Vector3DType u;  // Rotation vector, the cross product of x and a.
  numerics::cross_product(x.data(), a.data(), u.data());
  double sinT = u.norm();
  double cosT = numerics::dot_product(x.data(), a.data(), 3);
  double ccosT = 1 - cosT;

  // Degenerate case with angle near 0 or pi.
  if(utilities::isNearlyEqual(sinT, 0.0))
  {
    if(cosT > 0)
    {
      return numerics::Matrix<double>::identity(3);
    }
    else
    {
      // Give u a tiny component in any non-x direction
      // so we can rotate around it.
      u[1] = 1e-8;
    }
  }

  u = u.unitVector();
  numerics::Matrix<double> rot(3, 3, 0.0);
  rot(0, 0) = u[0] * u[0] * ccosT + cosT;
  rot(0, 1) = u[0] * u[1] * ccosT - u[2] * sinT;
  rot(0, 2) = u[0] * u[2] * ccosT + u[1] * sinT;
  rot(1, 0) = u[1] * u[0] * ccosT + u[2] * sinT;
  rot(1, 1) = u[1] * u[1] * ccosT + cosT;
  rot(1, 2) = u[1] * u[2] * ccosT - u[0] * sinT;
  rot(2, 0) = u[2] * u[0] * ccosT - u[1] * sinT;
  rot(2, 1) = u[2] * u[1] * ccosT + u[0] * sinT;
  rot(2, 2) = u[2] * u[2] * ccosT + cosT;

  return rot;
}

void SorClipper::extractClipperInfo()
{
  auto sorBaseArray = m_info.fetch_existing("sorBase").as_double_array();
  auto sorDirectionArray = m_info.fetch_existing("sorDirection").as_double_array();
  for(int d = 0; d < 3; ++d)
  {
    m_sorBase[d] = sorBaseArray[d];
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
