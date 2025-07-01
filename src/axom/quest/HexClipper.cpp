// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/quest/HexClipper.hpp"

namespace axom
{
namespace quest
{

HexClipper::HexClipper(const klee::Geometry& kGeom, const std::string& name)
  : GeometryClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("Hex") : name)
  , m_transformer(m_extTrans)
{
  extractClipperInfo();

  for(int i = 0; i < HexahedronType::NUM_HEX_VERTS; ++i)
  {
    m_hex[i] = m_transformer.getTransformed(m_hexBeforeTrans[i]);
  }

  m_hex.triangulate(m_tets);

  for(int i = 0; i < HexahedronType::NUM_HEX_VERTS; ++i)
  {
    m_hexBb.addPoint(m_hex[i]);
  }

  computeSurface();
}

bool HexClipper::labelInOut(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  AXOM_ANNOTATE_SCOPE("HexClipper::labelInOut");
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
void HexClipper::labelInOutImpl(quest::ShapeeMesh& shapeeMesh, axom::Array<LabelType>& labels)
{
  SLIC_ERROR_IF(shapeeMesh.dimension() != 3, "HexClipper requires a 3D mesh.");

  int allocId = shapeeMesh.getAllocatorID();
  auto cellCount = shapeeMesh.getCellCount();

  axom::ArrayView<const BoundingBox3DType> cellBbs = shapeeMesh.getCellBoundingBoxes();
  axom::ArrayView<const HexahedronType> cellsAsHexes = shapeeMesh.getCellsAsHexes();

  if(labels.size() < cellCount || labels.getAllocatorID() != shapeeMesh.getAllocatorID())
  {
    labels = axom::Array<LabelType>(ArrayOptions::Uninitialized(), cellCount, cellCount, allocId);
  }
  auto labelsView = labels.view();

  /*
    Label cell by eliminating where it can be w.r.t. the hex.
    The cell is conservatively represented by its bounding box, cellBb.
    - If cellBb doens't intersect hexBb, cell is OUT.
    - Else if cellBb intersects any of the hex's 24 surface triangles,
      the cell is ON the boundary.
    - Else if any cell vertex is in the hex (in any of the 24 tets),
      the cell is IN.  Otherwise OUT, since we've eliminated
      the possibility that it's ON the boundary.
  */
  const auto hexBb = m_hexBb;
  axom::Array<TetrahedronType> tets(m_tets.size(), m_tets.size(), allocId);
  axom::copy(tets.data(), m_tets.data(), m_tets.size() * sizeof(TetrahedronType));
  auto tetsView = tets.view();
  auto surfaceTriangles = m_surfaceTriangles;
  constexpr double eps = 1e-12;
  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
      auto& cellLabel = labelsView[cellId];
      auto& cellBb = cellBbs[cellId];
      // If bounding boxes don't intersect, nothing intersects.
      if(!hexBb.intersectsWith(cellBb))
      {
        cellLabel = LABEL_OUT;
        return;
      }
      // If cellBb intersects hex surface, there's a high chance cell does too.
      for(int ti = 0; ti < 24; ++ti)
      {
        const auto& surfTri = surfaceTriangles[ti];
        if(axom::primal::intersect(surfTri, cellBb))
        {
          cellLabel = LABEL_ON;
          return;
        }
      }
      // After eliminating possibility that cell is on surface,
      // the cell is either completely inside or completely out.
      // The cell is IN if any part of it is IN, so check an arbitrary vertex.
      // Note: Should the arbitrary vertex be some weird corner case, we could
      // use an alternative, like averaging two opposite corners, 0 and 6.
      const auto& cellAsHex = cellsAsHexes[cellId];
      const Point3DType& ptInCell(cellAsHex[0]);
      for(const auto& tet : tetsView)
      {
        if(tet.contains(ptInCell, eps))
        {
          cellLabel = LABEL_IN;
          return;
        }
      }
      cellLabel = LABEL_OUT;
    });

  return;
}

bool HexClipper::getGeometryAsTets(quest::ShapeeMesh& shapeeMesh, axom::Array<TetrahedronType>& tets)
{
  AXOM_ANNOTATE_SCOPE("HexClipper::getGeometryAsTets");
  int allocId = shapeeMesh.getAllocatorID();
  if(tets.getAllocatorID() != allocId || tets.size() != m_tets.size())
  {
    tets = axom::Array<TetrahedronType>(m_tets.size(), m_tets.size(), allocId);
  }
  axom::copy(tets.data(), m_tets.data(), m_tets.size() * sizeof(TetrahedronType));
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

void HexClipper::computeSurface()
{
  // Hex vertex shorthands
  // See the Hexahedron class documentation, especially the ASCII art.
  const auto& p = m_hex[0];
  const auto& q = m_hex[1];
  const auto& r = m_hex[2];
  const auto& s = m_hex[3];
  const auto& t = m_hex[4];
  const auto& u = m_hex[5];
  const auto& v = m_hex[6];
  const auto& w = m_hex[7];

  // 6 face centers, right handed, oriented inward.
  Point3DType pswt(axom::NumericArray<double, 3>{ p.array() + s.array() + w.array() + t.array() }/4);
  Point3DType quvr(axom::NumericArray<double, 3>{ q.array() + u.array() + v.array() + r.array() }/4);

  Point3DType ptuq(axom::NumericArray<double, 3>{ p.array() + t.array() + u.array() + q.array() }/4);
  Point3DType srvw(axom::NumericArray<double, 3>{ s.array() + r.array() + v.array() + w.array() }/4);

  Point3DType pqrs(axom::NumericArray<double, 3>{ p.array() + q.array() + r.array() + s.array() }/4);
  Point3DType twvu(axom::NumericArray<double, 3>{ t.array() + w.array() + v.array() + u.array() }/4);

  m_surfaceTriangles[0] = Triangle3DType(pswt, p, s);
  m_surfaceTriangles[1] = Triangle3DType(pswt, s, w);
  m_surfaceTriangles[2] = Triangle3DType(pswt, w, t);
  m_surfaceTriangles[3] = Triangle3DType(pswt, t, p);

  m_surfaceTriangles[4] = Triangle3DType(quvr, q, u);
  m_surfaceTriangles[5] = Triangle3DType(quvr, u, v);
  m_surfaceTriangles[6] = Triangle3DType(quvr, v, r);
  m_surfaceTriangles[7] = Triangle3DType(quvr, r, q);

  m_surfaceTriangles[8] = Triangle3DType(ptuq, p, t);
  m_surfaceTriangles[9] = Triangle3DType(ptuq, t, u);
  m_surfaceTriangles[10] = Triangle3DType(ptuq, u, q);
  m_surfaceTriangles[11] = Triangle3DType(ptuq, q, p);

  m_surfaceTriangles[12] = Triangle3DType(srvw, s, r);
  m_surfaceTriangles[13] = Triangle3DType(srvw, r, v);
  m_surfaceTriangles[14] = Triangle3DType(srvw, v, w);
  m_surfaceTriangles[15] = Triangle3DType(srvw, w, s);

  m_surfaceTriangles[16] = Triangle3DType(pqrs, p, q);
  m_surfaceTriangles[17] = Triangle3DType(pqrs, q, r);
  m_surfaceTriangles[18] = Triangle3DType(pqrs, r, s);
  m_surfaceTriangles[19] = Triangle3DType(pqrs, s, p);

  m_surfaceTriangles[20] = Triangle3DType(twvu, t, w);
  m_surfaceTriangles[21] = Triangle3DType(twvu, w, v);
  m_surfaceTriangles[22] = Triangle3DType(twvu, v, u);
  m_surfaceTriangles[23] = Triangle3DType(twvu, u, t);
}

}  // end namespace quest
}  // end namespace axom
