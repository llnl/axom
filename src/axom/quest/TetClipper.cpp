// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifdef AXOM_USE_CONDUIT
  #include "conduit_blueprint.hpp"
#endif

#include "axom/quest/TetClipper.hpp"
#include "axom/fmt.hpp"

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

  /*
    Compute whether mesh vertices are inside shape.
  */
  axom::Array<bool> vertIsInside {ArrayOptions::Uninitialized(), vertCount, vertCount, allocId};
  auto vertIsInsideView = vertIsInside.view();
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vX.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vY.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vZ.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vertIsInsideView.getAllocatorID()));

  auto bb = m_bb;
  auto planes = m_planes;
  axom::for_all<ExecSpace>(
    vertCount,
    AXOM_LAMBDA(axom::IndexType vertId) {
      primal::Point3D vert {vX[vertId], vY[vertId], vZ[vertId]};
      vertIsInsideView[vertId] = false;

      if(bb.contains(vert))
      {
        // Detailed check for verts inside the hex's bounding box.
        // If vert is in any tet, it's in the hex.
        // Be conservative: primal::ON_BOUNDARY is considered inside.
        // TODO: See if using Tetrahedron::physToBarycentric is faster.
        bool isInsideTet = true;
        for(int iPlane = 0; iPlane < 4; ++iPlane)
        {
          const auto& plane = planes[iPlane];
          isInsideTet &= plane.getOrientation(vert, 0.0) != primal::ON_NEGATIVE_SIDE;
        }
        if(isInsideTet)
        {
          vertIsInsideView[vertId] = true;
        }
      }
    });

  vertexInsideToCellLabel<ExecSpace>(shapeeMesh, vertIsInsideView, labels);

  /*
    Look for edge-tet intersections missed by vertex-only checks.
    These are small but real errors that can be reduced with mesh
    resolution, but cannot be eliminated without checking edges.
    With this correction there should be no clipping volume error
    except from floating-point round-off.
  */
  labelByEdges<ExecSpace>(shapeeMesh, labels);

  return;
}

/*
  Label cell outside if no vertex is in the bounding box.
  Otherwise, label it on boundary, because we don't know.
*/
template <typename ExecSpace>
void TetClipper::vertexInsideToCellLabel(
  quest::ShapeeMesh& shapeeMesh,
  axom::ArrayView<bool>& vertIsInside,
  axom::Array<LabelType>& labels)
{
  axom::ArrayView<const axom::IndexType, 2> connView = shapeeMesh.getConnectivity();
  SLIC_ASSERT(connView.shape() ==
              (axom::StackArray<axom::IndexType, 2> {shapeeMesh.getCellCount(), HexahedronType::NUM_HEX_VERTS}));

  if(labels.size() < shapeeMesh.getCellCount() || labels.getAllocatorID() != shapeeMesh.getAllocatorID())
  {
    labels = axom::Array<LabelType>(ArrayOptions::Uninitialized(), shapeeMesh.getCellCount(), shapeeMesh.getCellCount(), shapeeMesh.getAllocatorID());
  }
  auto labelsView = labels.view();

  axom::for_all<ExecSpace>(
    shapeeMesh.getCellCount(),
    AXOM_LAMBDA(axom::IndexType cellId) {
      auto cellVertIds = connView[cellId];
      bool hasIn = vertIsInside[cellVertIds[0]];
      bool hasOut = !hasIn;
      for(int vi = 0; vi < HexahedronType::NUM_HEX_VERTS; ++vi)
      {
        int vertId = cellVertIds[vi];
        bool isIn = vertIsInside[vertId];
        hasIn |= isIn;
        hasOut |= !isIn;
      }
      labelsView[cellId] = !hasOut ? LABEL_IN : !hasIn ? LABEL_OUT : LABEL_ON;
    });

  return;
}

template <typename ExecSpace>
void TetClipper::labelByEdges(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  /*
    If a cell's vertices are all outside, check if any edge crosses
    the geometry.  This check can find small overlaps that the
    vertex-only checks miss.
  */
  axom::ArrayView<const TetrahedronType> cellsAsTets = shapeeMesh.getCellsAsTets();
  axom::ArrayView<const BoundingBox3DType> cellBbs = shapeeMesh.getCellBoundingBoxes();
  auto labelsView = labels.view();
  constexpr int NUM_TETS_PER_HEX = primal::Hexahedron<double, 3>::NUM_TRIANGULATE;

  auto geomTet = m_tet;
  auto bb = m_bb;
  axom::for_all<ExecSpace>(
    labels.size(),
    AXOM_LAMBDA(axom::IndexType cellId) {
      LabelType& cellLabel = labelsView[cellId];
      const BoundingBox3DType cellBb = cellBbs[cellId];
      if (cellLabel == LABEL_OUT && bb.intersectsWith(cellBb))
      {
        const axom::IndexType tetIdxStart = cellId * NUM_TETS_PER_HEX;
        const axom::IndexType tetIdxEnd = (1 + cellId) * NUM_TETS_PER_HEX;
        for(axom::IndexType ti = tetIdxStart; ti < tetIdxEnd && cellLabel == LABEL_OUT; ++ti)
        {
          const TetrahedronType& cellTet = cellsAsTets[ti];
          if(axom::primal::intersect(cellTet, geomTet))
          {
            cellLabel = LABEL_ON;
          }
        }
      }
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
