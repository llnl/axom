// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_HEXCLIPPER_HPP
#define AXOM_QUEST_HEXCLIPPER_HPP

#include "axom/klee/Geometry.hpp"
#include "axom/quest/GeometryClipperStrategy.hpp"

namespace axom
{
namespace quest
{

/*!
  @brief GeometryClipper specialized for sphere geometries.
*/
class HexClipper : public GeometryClipperStrategy
{
public:
  /*!
    @brief Constructor.

    @param [in] kGeom Describes the shape to place
      into the mesh.
    @param [in] name To override the default strategy name
  */
  HexClipper(const klee::Geometry& kGeom, const std::string& name = "");

  virtual ~HexClipper() = default;

  const std::string& name() const override { return m_name; }

  bool labelInOut(quest::ShapeeMesh& shappeMesh, axom::Array<char>& label) override;

  bool getGeometryAsTets(quest::ShapeeMesh& shappeMesh, axom::Array<TetrahedronType>& tets) override;

#if !defined(__CUDACC__)
private:
#endif
  std::string m_name;

  HexahedronType m_hex;
  axom::primal::BoundingBox<double, 3> m_bb;
  axom::StackArray<TetrahedronType, HexahedronType::NUM_TRIANGULATE> m_tets;
  //!@brief 4 planes per tet, each oriented to the interior of the tet.
  axom::StackArray<axom::StackArray<Plane3DType, TetrahedronType::NUM_VERTS>, HexahedronType::NUM_TRIANGULATE>
    m_planes;

  template <typename ExecSpace>
  void labelInOutImpl(quest::ShapeeMesh& shapeeMesh, axom::Array<char>& label);

  void extractClipperInfo();
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_HEXCLIPPER_HPP
