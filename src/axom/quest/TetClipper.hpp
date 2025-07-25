// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_TETCLIPPER_HPP
#define AXOM_QUEST_TETCLIPPER_HPP

#include "axom/klee/Geometry.hpp"
#include "axom/quest/GeometryClipperStrategy.hpp"
#include "axom/primal/geometry/CoordinateTransformer.hpp"

namespace axom
{
namespace quest
{

/*!
  @brief Geometry clipping operations for sphere geometries.
*/
class TetClipper : public GeometryClipperStrategy
{
public:
  /*!
    @brief Constructor.

    @param [in] kGeom Describes the shape to place
      into the mesh.
    @param [in] name To override the default strategy name
  */
  TetClipper(const klee::Geometry& kGeom, const std::string& name = "");

  virtual ~TetClipper() = default;

  const std::string& name() const override { return m_name; }

  bool labelInOut(quest::ShapeeMesh& shappeMesh, axom::Array<char>& label) override;

  bool getGeometryAsTets(quest::ShapeeMesh& shappeMesh, axom::Array<TetrahedronType>& tets) override;

#if !defined(__CUDACC__)
private:
#endif
  std::string m_name;

  //!@brief Tetrahedron before transformation.
  TetrahedronType m_tetBeforeTrans;

  //!@brief Tetrahedron after transformation.
  TetrahedronType m_tet;

  axom::primal::BoundingBox<double, 3> m_bb;

  //!@brief 4 planes per Tet, each oriented to the interior of the tet.
  axom::StackArray<Plane3DType, 4> m_planes;

  //!@brief 4 triangular facets of Tet, each oriented to the interior of the tet.
  axom::StackArray<Triangle3DType, 4> m_facets;

  axom::primal::CoordinateTransformer<double> m_transformer;

  template <typename ExecSpace>
  void labelInOutImpl(quest::ShapeeMesh& shapeeMesh, axom::Array<char>& label);

  template <typename ExecSpace>
  void vertexSignedDistToLabel(quest::ShapeeMesh& shapeeMesh,
                               axom::ArrayView<const double> vertexSignedDists,
                               axom::Array<LabelType>& labels);

  // Extract clipper info from GeometryClipperStrategy::m_info.
  void extractClipperInfo();
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_TETCLIPPER_HPP
