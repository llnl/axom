// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/quest/TetClipper.hpp"

namespace axom
{
namespace quest
{

TetClipper::TetClipper(const klee::Geometry& kGeom, const std::string& name)
  : GeometryClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("Tet") : name)
  , m_transformer(m_extTrans)
{
  extractClipperInfo();

  for(int i = 0; i < TetrahedronType::NUM_VERTS; ++i)
  {
    m_tet[i] = m_transformer.getTransformed(m_tetBeforeTrans[i]);
  }

  for(int i = 0; i < TetrahedronType::NUM_VERTS; ++i)
  {
    m_bb.addPoint(m_tet[i]);
  }

  for(int iPlane = 0; iPlane < 4; ++iPlane)
  {
    const Point3DType& a = m_tet[iPlane % 4];
    const Point3DType& b = m_tet[(iPlane + 1) % 4];
    const Point3DType& c = m_tet[(iPlane + 2) % 4];
    m_planes[iPlane] = axom::primal::make_plane(a, b, c);
    // For tet points ordered by right hand rule, odd planes
    // face outside.  Flip them to face inside.
    if(iPlane % 2 == 1) m_planes[iPlane].flip();
  }
  m_facets = m_tet.facets();

}

bool TetClipper::labelInOut(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  AXOM_ANNOTATE_BEGIN("TetClipper::labelInOut");
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
  AXOM_ANNOTATE_END("TetClipper::labelInOut");
  return true;
}

template <typename ExecSpace>
void TetClipper::labelInOutImpl(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  SLIC_ERROR_IF(shapeeMesh.dimension() != 3, "TetClipper requires a 3D mesh.");

  int allocId = shapeeMesh.getAllocatorID();
  auto vertCount = shapeeMesh.getVertexCount();

  const auto& vertCoords = shapeeMesh.getVertexCoords3D();
  const auto& vX = vertCoords[0];
  const auto& vY = vertCoords[1];
  const auto& vZ = vertCoords[2];
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vX.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vY.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vZ.getAllocatorID()));

  /*
    Compute signed distances to mesh vertices.
  */
  axom::Array<double> signedDists{ArrayOptions::Uninitialized(), vertCount, vertCount, allocId};
  auto signedDistsView = signedDists.view();
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(signedDists.getAllocatorID()));

  const auto tet = m_tet;
  const auto bb = m_bb;
  auto planes = m_planes;
  axom::for_all<ExecSpace>(
    vertCount,
    AXOM_LAMBDA(axom::IndexType vertId) {
      primal::Point3D vert {vX[vertId], vY[vertId], vZ[vertId]};

      auto& signedDist = signedDistsView[vertId];

      if(bb.contains(vert))
      {
        bool isInsideTet = true;
        double minDist = 0.0;
        for(const auto& plane : planes)
        {
          auto signedDistToPlane = plane.signedDistance(vert);
          isInsideTet &= signedDistToPlane >= 0;
          minDist = axom::utilities::min(minDist, signedDistToPlane);
        }

        if(isInsideTet)
        {
          SLIC_ASSERT(minDist >= 0);
          signedDist = minDist;
        }
        else
        {
          signedDist = sqrt(axom::primal::squared_distance(vert, tet));
        }
      }
      else
      {
        signedDist = sqrt(axom::primal::squared_distance(vert, tet));
      }
    });

  vertexSignedDistToLabel<ExecSpace>(shapeeMesh, signedDists.view(), labels);

  return;
}

template <typename ExecSpace>
void TetClipper::vertexSignedDistToLabel(quest::ShapeeMesh& shapeeMesh,
                                         axom::ArrayView<const double> vertexSignedDists,
                                         axom::Array<LabelType>& labels)
{
  /*
    Label cell by whether it has vertices inside, outside or both.
    The tet may intersect cell between its vertices,
    so we classify the cells conservatively.
    - If cell has any vertex less than a small distance outside the
      tet, the cell is classified as having a vertex inside.
    - We don't need to do the same for vertices a small distance
      inside the tet, because the tet is convex.  No cell with all
      vertices inside the tet can intersect the tet.
  */
  auto cellLengths = shapeeMesh.getCellLengths();
  const double lenFactor = 0.5;

  axom::IndexType cellCount = shapeeMesh.getCellCount();

  axom::ArrayView<const axom::IndexType, 2> connView = shapeeMesh.getCellNodeConnectivity();
  SLIC_ASSERT(connView.shape() ==
              (axom::StackArray<axom::IndexType, 2> {cellCount, HexahedronType::NUM_HEX_VERTS}));

  if(labels.size() < cellCount || labels.getAllocatorID() != shapeeMesh.getAllocatorID())
  {
    labels = axom::Array<LabelType>(ArrayOptions::Uninitialized(), cellCount, cellCount, shapeeMesh.getAllocatorID());
  }
  auto labelsView = labels.view();

  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
      LabelType& cellLabel = labelsView[cellId];
      auto cellVertIds = connView[cellId];
      const double proximityThreshold = cellLengths[cellId]*lenFactor;
      bool hasIn = vertexSignedDists[cellVertIds[0]] < proximityThreshold;
      bool hasOut = vertexSignedDists[cellVertIds[0]] > 0;
      for(int vi = 1; vi < HexahedronType::NUM_HEX_VERTS; ++vi)
      {
        int vertId = cellVertIds[vi];
        bool isIn = vertexSignedDists[vertId] < proximityThreshold;
        bool isOut = vertexSignedDists[vertId] > 0;
        hasIn |= isIn;
        hasOut |= isOut;
      }
      cellLabel = !hasOut ? LABEL_IN : !hasIn ? LABEL_OUT : LABEL_ON;
    });
}

bool TetClipper::getGeometryAsTets(quest::ShapeeMesh& shapeeMesh, axom::Array<TetrahedronType>& tets)
{
  AXOM_ANNOTATE_BEGIN("TetClipper::getGeometryAsTets");
  int allocId = shapeeMesh.getAllocatorID();
  if(tets.getAllocatorID() != allocId || tets.size() != 1)
  {
    tets = axom::Array<TetrahedronType>(1, 1, allocId);
  }
  axom::copy(tets.data(), &m_tet, sizeof(TetrahedronType));
  AXOM_ANNOTATE_END("TetClipper::getGeometryAsTets");
  return true;
}

void TetClipper::extractClipperInfo()
{
  const auto v0 = m_info.fetch_existing("v0").as_double_array();
  const auto v1 = m_info.fetch_existing("v1").as_double_array();
  const auto v2 = m_info.fetch_existing("v2").as_double_array();
  const auto v3 = m_info.fetch_existing("v3").as_double_array();
  for(int d = 0; d < 3; ++d)
  {
    m_tetBeforeTrans[0][d] = v0[d];
    m_tetBeforeTrans[1][d] = v1[d];
    m_tetBeforeTrans[2][d] = v2[d];
    m_tetBeforeTrans[3][d] = v3[d];
  }
}

}  // end namespace quest
}  // end namespace axom
