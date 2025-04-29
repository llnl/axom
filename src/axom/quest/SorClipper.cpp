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
  , m_cellCount(0)
{
  extractClipperInfo();
}

bool SorClipper::labelInOut(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  AXOM_UNUSED_VAR(shapeeMesh);
  AXOM_UNUSED_VAR(labels);
  return false;  // TODO: implement SorClipper::labelInOut
}

// TODO: Factor out execution-space-specific computations into a template method,
// instead of computing on host and copying to GPU.
bool SorClipper::getShapeAsOcts(quest::ShapeeMesh& shapeeMesh, axom::Array<OctahedronType>& octs)
{
  const int hostAllocId = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int allocId = shapeeMesh.getAllocatorId();

  if(octs.getAllocatorID() != allocId || octs.size() != m_cellCount)
  {
    octs = axom::Array<OctahedronType>(m_cellCount, m_cellCount, allocId);
  }

  axom::Array<OctahedronType> tmpOcts(0, 0, hostAllocId);

  // Host loop writes directly to octs if possible, else a temporary array.
  axom::Array<OctahedronType>& octsOnHost =
    axom::execution_space<axom::SEQ_EXEC>::usesAllocId(octs.getAllocatorID()) ? octs : tmpOcts;

  if(&octsOnHost == &tmpOcts)
  {
    tmpOcts.resize(m_cellCount, OctahedronType());
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
  numerics::Matrix<double> rotate = sorAxisRotMatrix(m_sorDirection);
  const auto& translate = m_sorBase;
  auto octsView = octsOnHost.view();
  axom::for_all<axom::SEQ_EXEC>(
    octCount,
    AXOM_LAMBDA(axom::IndexType iOct) {
      auto& oct = octsView[iOct];
      for(int iVert = 0; iVert < OctType::NUM_VERTS; ++iVert)
      {
        auto& newCoords = oct[iVert];
        auto oldCoords = newCoords;
        numerics::matrix_vector_multiply(rotate, oldCoords.data(), newCoords.data());
        newCoords.array() += translate.array();
      }
    });

  if(&octs != &octsOnHost)
  {
    octs.resize(octsOnHost.size());
    axom::copy(octs.data(), octsOnHost.data(), sizeof(OctahedronType) * octs.size());
  }

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
  m_info.print();

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
