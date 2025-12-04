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
}

bool TetClipper::getGeometryAsTets(quest::experimental::ShapeMesh& shapeMesh,
                                   axom::Array<TetrahedronType>& tets)
{
  AXOM_ANNOTATE_SCOPE("TetClipper::getGeometryAsTets");
  int allocId = shapeMesh.getAllocatorID();
  if(tets.getAllocatorID() != allocId || tets.size() != 1)
  {
    tets = axom::Array<TetrahedronType>(1, 1, allocId);
  }
  // Copy tet into tets array, which may be in non-host memory.
  axom::copy(tets.data(), &m_tet, sizeof(TetrahedronType));
  return true;
}

void TetClipper::extractClipperInfo()
{
  const auto v0 = m_info.fetch_existing("v0").as_double_array();
  const auto v1 = m_info.fetch_existing("v1").as_double_array();
  const auto v2 = m_info.fetch_existing("v2").as_double_array();
  const auto v3 = m_info.fetch_existing("v3").as_double_array();
  SLIC_ASSERT(v0.number_of_elements() == 3);
  SLIC_ASSERT(v1.number_of_elements() == 3);
  SLIC_ASSERT(v2.number_of_elements() == 3);
  SLIC_ASSERT(v3.number_of_elements() == 3);
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
