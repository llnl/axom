// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_TETCLIPPER_HPP
#define AXOM_QUEST_TETCLIPPER_HPP

#include "axom/klee/Geometry.hpp"
#include "axom/quest/MeshClipperStrategy.hpp"
#include "axom/primal/geometry/CoordinateTransformer.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{

/*!
 * @brief Geometry clipping operations for sphere geometries.
*/
class TetClipper : public MeshClipperStrategy
{
public:
  /*!
   * @brief Constructor.
   *
   * @param [in] kGeom Describes the shape to place
   *   into the mesh.
   * @param [in] name To override the default strategy name
  */
  TetClipper(const klee::Geometry& kGeom, const std::string& name = "");

  virtual ~TetClipper() = default;

  const std::string& name() const override { return m_name; }

  bool labelCellsInOut(quest::experimental::ShapeeMesh& shappeMesh, axom::Array<char>& label) override;

  bool labelTetsInOut(quest::experimental::ShapeeMesh& shapeeMesh,
                      axom::ArrayView<const axom::IndexType> cellsOnBdry,
                      axom::Array<LabelType>& tetLabels) override;

  bool getGeometryAsTets(quest::experimental::ShapeeMesh& shappeMesh, axom::Array<TetrahedronType>& tets) override;

#if !defined(__CUDACC__)
private:
#endif
  std::string m_name;

  //!@brief Tetrahedron before transformation.
  TetrahedronType m_tetBeforeTrans;

  //!@brief Tetrahedron after transformation.
  TetrahedronType m_tet;

  axom::primal::BoundingBox<double, 3> m_bb;

  //!@brief 4 planes of the Tet, oriented to the interior of the tet.
  axom::StackArray<Plane3DType, 4> m_planes;

  //!@brief Height of the tet when resting on each facet.
  axom::StackArray<double, 4> m_heights;

  axom::primal::experimental::CoordinateTransformer<double> m_transformer;

  template <typename ExecSpace>
  void labelInOutImpl(quest::experimental::ShapeeMesh& shapeeMesh, axom::Array<char>& label);

  template <typename ExecSpace>
  void labelTetsInOutImpl(quest::experimental::ShapeeMesh& shapeeMesh,
                          axom::ArrayView<const axom::IndexType> cellsOnBdry,
                          axom::Array<LabelType>& tetLabels);

  // Extract clipper info from MeshClipperStrategy::m_info.
  void extractClipperInfo();
};

}  // namespace experimental
}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_TETCLIPPER_HPP
