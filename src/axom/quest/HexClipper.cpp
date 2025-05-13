// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifdef AXOM_USE_CONDUIT
  #include "conduit_blueprint.hpp"
#endif

#include "axom/quest/Discretize.hpp"
#include "axom/quest/HexClipper.hpp"
#include "axom/fmt.hpp"
#include "axom/core/WhereMacro.hpp"

namespace axom
{
namespace quest
{

HexClipper::HexClipper(const klee::Geometry& kGeom, const std::string& name)
  : GeometryClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("Hex") : name)
  , m_transformer(m_transMat)
{
  extractClipperInfo();

  for(int i = 0; i < HexahedronType::NUM_HEX_VERTS; ++i)
  {
    m_hex[i] = m_transformer.getTransform(m_hexBeforeTrans[i]);
  }

  m_hex.triangulate(m_tets);

  for(int i = 0; i < HexahedronType::NUM_HEX_VERTS; ++i)
  {
    m_bb.addPoint(m_hex[i]);
  }
std::cout<<__WHERE<<m_hexBeforeTrans<<' '<<m_hex<<' '<<m_hex<<' '<<m_bb<<std::endl;

  for(int iTet = 0; iTet < HexahedronType::NUM_TRIANGULATE; ++iTet)
  {
    const auto& tet = m_tets[iTet];
    for(int iPlane = 0; iPlane < 4; ++iPlane)
    {
      const Point3DType& a = tet[iPlane % 4];
      const Point3DType& b = tet[(iPlane + 1) % 4];
      const Point3DType& c = tet[(iPlane + 2) % 4];
      m_planes[iTet][iPlane] = axom::primal::make_plane(a, b, c);
      // For tet points ordered by right hand rule, odd planes
      // face outside.  Flip them to face inside.
      if(iPlane % 2 == 1) m_planes[iTet][iPlane].flip();
    }
  }
}

bool HexClipper::labelInOut(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  AXOM_ANNOTATE_BEGIN("HexClipper::labelInOut");
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
  AXOM_ANNOTATE_END("HexClipper::labelInOut");
  return true;
}

template <typename ExecSpace>
void HexClipper::labelInOutImpl(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  SLIC_ERROR_IF(shapeeMesh.dimension() != 3, "HexClipper requires a 3D mesh.");

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
        for(int iTet = 0; iTet < HexahedronType::NUM_TRIANGULATE; ++iTet)
        {
          bool isInsideTet = true;
          for(int iPlane = 0; iPlane < 4; ++iPlane)
          {
            const auto& plane = planes[iTet][iPlane];
            isInsideTet &= plane.getOrientation(vert, 0.0) != primal::ON_NEGATIVE_SIDE;
          }
          if(isInsideTet)
          {
            vertIsInsideView[vertId] = true;
            break;
          }
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
      if(labelsView[cellId] == LABEL_OUT)
      {
        // Look for possible shallow intersection of cell edges may make with hex.
        // Will need to compute convex edges of the hex: 12 regular edges and
        // 6 diagonal edges.  The 6 alternate diagnoal edges are concave so they
        // don't need checking.  TODO:
      }
    });

  return;
}

bool HexClipper::getGeometryAsTets(quest::ShapeeMesh& shapeeMesh, axom::Array<TetrahedronType>& tets)
{
  AXOM_ANNOTATE_BEGIN("HexClipper::getGeometryAsTets");
  int allocId = shapeeMesh.getAllocatorId();
  if(tets.getAllocatorID() != allocId || tets.size() != HexahedronType::NUM_TRIANGULATE)
  {
    tets = axom::Array<TetrahedronType>(HexahedronType::NUM_TRIANGULATE,
                                        HexahedronType::NUM_TRIANGULATE,
                                        allocId);
  }
  axom::copy(tets.data(), m_tets.data(), m_tets.size() * sizeof(TetrahedronType));
  AXOM_ANNOTATE_END("HexClipper::getGeometryAsTets");
  return true;
}

void HexClipper::extractClipperInfo()
{
  const auto v0 = m_info.fetch_existing("v0").as_double_array();
  const auto v1 = m_info.fetch_existing("v1").as_double_array();
  const auto v2 = m_info.fetch_existing("v2").as_double_array();
  const auto v3 = m_info.fetch_existing("v3").as_double_array();
  const auto v4 = m_info.fetch_existing("v4").as_double_array();
  const auto v5 = m_info.fetch_existing("v5").as_double_array();
  const auto v6 = m_info.fetch_existing("v6").as_double_array();
  const auto v7 = m_info.fetch_existing("v7").as_double_array();
  for(int d = 0; d < 3; ++d)
  {
    m_hexBeforeTrans[0][d] = v0[d];
    m_hexBeforeTrans[1][d] = v1[d];
    m_hexBeforeTrans[2][d] = v2[d];
    m_hexBeforeTrans[3][d] = v3[d];
    m_hexBeforeTrans[4][d] = v4[d];
    m_hexBeforeTrans[5][d] = v5[d];
    m_hexBeforeTrans[6][d] = v6[d];
    m_hexBeforeTrans[7][d] = v7[d];
  }
}

}  // end namespace quest
}  // end namespace axom
