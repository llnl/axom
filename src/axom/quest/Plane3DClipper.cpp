// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/quest/Plane3DClipper.hpp"

namespace axom
{
namespace quest
{

Plane3DClipper::Plane3DClipper(const klee::Geometry& kGeom, const std::string& name)
  : GeometryClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("Plane3D") : name)
{
  extractClipperInfo();
}

bool Plane3DClipper::labelInOut(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  AXOM_ANNOTATE_BEGIN("Plane3DClipper::labelInOut");
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
  AXOM_ANNOTATE_END("Plane3DClipper::labelInOut");
  return true;
}

bool Plane3DClipper::specializedClip(quest::ShapeeMesh& shapeeMesh,
                                     axom::ArrayView<double> ovlap,
                                     const axom::ArrayView<IndexType>& cellIds)
{
  AXOM_ANNOTATE_BEGIN("Plane3DClipper::specializedClip");
  switch(shapeeMesh.getRuntimePolicy())
  {
  case axom::runtime_policy::Policy::seq:
    specializedClipImpl<axom::SEQ_EXEC>(shapeeMesh, ovlap, cellIds);
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    specializedClipImpl<axom::OMP_EXEC>(shapeeMesh, ovlap, cellIds);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case axom::runtime_policy::Policy::cuda:
    specializedClipImpl<axom::CUDA_EXEC<256>>(shapeeMesh, ovlap, cellIds);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    specializedClipImpl<axom::HIP_EXEC<256>>(shapeeMesh, ovlap, cellIds);
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
  AXOM_ANNOTATE_END("Plane3DClipper::specializedClip");
  return true;
}

template <typename ExecSpace>
void Plane3DClipper::labelInOutImpl(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  SLIC_ERROR_IF(shapeeMesh.dimension() != 3, "Plane3DClipper requires a 3D mesh.");

  constexpr int NUM_VERTS_PER_CELL = 8;

  int allocId = shapeeMesh.getAllocatorID();
  auto cellCount = shapeeMesh.getCellCount();
  auto vertCount = shapeeMesh.getVertexCount();

  const auto& vertCoords = shapeeMesh.getVertexCoords3D();
  const auto& vX = vertCoords[0];
  const auto& vY = vertCoords[1];
  const auto& vZ = vertCoords[2];

  /*
    Compute whether vertices are inside shape.
  */
  axom::Array<bool> vertIsInside {ArrayOptions::Uninitialized(), vertCount, vertCount, allocId};
  auto vertIsInsideView = vertIsInside.view();
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vX.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vY.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vZ.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vertIsInsideView.getAllocatorID()));

  auto plane = m_plane;
  axom::for_all<ExecSpace>(
    vertCount,
    AXOM_LAMBDA(axom::IndexType vertId) {
      primal::Point3D vert {vX[vertId], vY[vertId], vZ[vertId]};
      double signedDist = plane.signedDistance(vert);
      vertIsInsideView[vertId] = signedDist > 0;
    });

  if(labels.size() < cellCount || labels.getAllocatorID() != shapeeMesh.getAllocatorID())
  {
    labels = axom::Array<LabelType>(ArrayOptions::Uninitialized(), cellCount, cellCount, allocId);
  }

  /*
    Label cell by whether it has vertices inside, outside or both.
  */
  axom::ArrayView<const axom::IndexType, 2> connView = shapeeMesh.getConnectivity();
  SLIC_ASSERT(connView.shape()[1] == NUM_VERTS_PER_CELL);

  auto labelsView = labels.view();

  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
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
      labelsView[cellId] = !hasOut ? 0 : !hasIn ? 2 : 1;
    });

  return;
}

template <typename ExecSpace>
void Plane3DClipper::specializedClipImpl(quest::ShapeeMesh& shapeeMesh,
                                         axom::ArrayView<double>& ovlap,
                                         const axom::ArrayView<IndexType>& cellIds)
{
  constexpr int NUM_VERTS_PER_CELL = 8;
  constexpr int NUM_TETS_PER_HEX = 24;
  constexpr double EPS = 1e-10;

  const axom::ArrayView<const axom::IndexType, 2> connView = shapeeMesh.getConnectivity();
  SLIC_ASSERT(connView.shape()[1] == NUM_VERTS_PER_CELL);

  auto& vertCoords = shapeeMesh.getVertexCoords3D();
  const auto& x = vertCoords[0];
  const auto& y = vertCoords[1];
  const auto& z = vertCoords[2];

  using TetsInHex = axom::StackArray<primal::Tetrahedron<double, 3>, NUM_TETS_PER_HEX>;

  auto plane = m_plane;
  axom::for_all<ExecSpace>(
    cellIds.size(),
    AXOM_LAMBDA(axom::IndexType i) {
      axom::IndexType cellId = cellIds[i];

      auto cellVerts = connView[cellId];
      primal::Hexahedron<double, 3> hex;
      for(int vi = 0; vi < NUM_VERTS_PER_CELL; ++vi)
      {
        axom::IndexType vertId = cellVerts[vi];
        auto& vCoords = hex[vi];
        vCoords[0] = x[vertId];
        vCoords[1] = y[vertId];
        vCoords[2] = z[vertId];
      }

      TetsInHex tetsInHex;
      hex.triangulate(tetsInHex);

      double vol = 0.0;
      for(int ti = 0; ti < NUM_TETS_PER_HEX; ++ti)
      {
        const auto& tet = tetsInHex[ti];
        primal::Polyhedron<double, 3> overlap = primal::clip(tet, plane, EPS);
        vol += overlap.volume();
      }
      ovlap[cellId] = vol;
    });
}

void Plane3DClipper::extractClipperInfo()
{
  const auto normal = m_info.fetch_existing("normal").as_double_array();
  const double offset = m_info.fetch_existing("offset").as_double();
  Vector3DType nVec;
  for(int d = 0; d < 3; ++d)
  {
    nVec[d] = normal[d];
  }
  m_plane = Plane3DType(nVec, offset);
}

}  // end namespace quest
}  // end namespace axom
