// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/quest/detail/clipping/TetClipper.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{

TetClipper::TetClipper(const klee::Geometry& kGeom, const std::string& name)
  : MeshClipperStrategy(kGeom)
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
    // face outside.  Make them face inside.
    if(iPlane % 2 == 1) m_planes[iPlane].flip();

    const Point3DType& apex = m_tet[(iPlane + 3) % 4];
    m_heights[iPlane] = m_planes[iPlane].signedDistance(apex);
    SLIC_ASSERT(m_heights[iPlane] >= 0);
  }
}

bool TetClipper::getGeometryAsTets(quest::experimental::ShapeMesh& shapeMesh, axom::Array<TetrahedronType>& tets)
{
  AXOM_ANNOTATE_SCOPE("TetClipper::getGeometryAsTets");
  int allocId = shapeMesh.getAllocatorID();
  if(tets.getAllocatorID() != allocId || tets.size() != 1)
  {
    tets = axom::Array<TetrahedronType>(1, 1, allocId);
  }
  axom::copy(tets.data(), &m_tet, sizeof(TetrahedronType));
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

}  // namespace experimental
}  // end namespace quest
}  // end namespace axom
