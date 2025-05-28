// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_HEXCLIPPER_HPP
#define AXOM_QUEST_HEXCLIPPER_HPP

#include "axom/klee/Geometry.hpp"
#include "axom/quest/GeometryClipperStrategy.hpp"
#include "axom/primal/geometry/CoordinateTransformer.hpp"

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

  /*!
    If a mesh cell has all vertices outside the geometry, it labeled outside.
    This will miss cases where an edge of the cell passes through the geometry.

    TODO: Fix cell mislabled as outside.  For cells with all point
    outside the geometry, also check all its edges to see if any
    intersects the geometry.
  */
  bool labelInOut(quest::ShapeeMesh& shappeMesh, axom::Array<char>& label) override;

  bool getGeometryAsTets(quest::ShapeeMesh& shappeMesh, axom::Array<TetrahedronType>& tets) override;

#if !defined(__CUDACC__)
private:
#endif
  std::string m_name;

  //!@brief Hexahedron before transformation.
  HexahedronType m_hexBeforeTrans;

  //!@brief Hexahedron after transformation.
  HexahedronType m_hex;

  //!@brief Bounding box of m_hex.
  axom::primal::BoundingBox<double, 3> m_bb;

  //!@brief Tetrahedralized version of of m_hex.
  axom::StackArray<TetrahedronType, HexahedronType::NUM_TRIANGULATE> m_tets;

  //!@brief 4 planes per tet in m_tets, each oriented to the interior of the tet.
  axom::StackArray<axom::StackArray<Plane3DType, TetrahedronType::NUM_VERTS>, HexahedronType::NUM_TRIANGULATE>
    m_planes;

  axom::primal::CoordinateTransformer<double> m_transformer;

  template <typename ExecSpace>
  void labelInOutImpl(quest::ShapeeMesh& shapeeMesh, axom::Array<char>& label);

  void extractClipperInfo();
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_HEXCLIPPER_HPP
