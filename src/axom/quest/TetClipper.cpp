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
{
  extractClipperInfo();

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

  constexpr int NUM_VERTS_PER_CELL = 8;

  int allocId = shapeeMesh.getAllocatorId();
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

  if(labels.size() < cellCount || labels.getAllocatorID() != shapeeMesh.getAllocatorId())
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
      labelsView[cellId] = !hasOut ? LABEL_IN : !hasIn ? LABEL_OUT : LABEL_ON;
    });

  return;
}

bool TetClipper::getGeometryAsTets(quest::ShapeeMesh& shapeeMesh, axom::Array<TetrahedronType>& tets)
{
  AXOM_ANNOTATE_BEGIN("TetClipper::getGeometryAsTets");
  int allocId = shapeeMesh.getAllocatorId();
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
    m_tet[0][d] = v0[d];
    m_tet[1][d] = v1[d];
    m_tet[2][d] = v2[d];
    m_tet[3][d] = v3[d];
  }
}

}  // end namespace quest
}  // end namespace axom
