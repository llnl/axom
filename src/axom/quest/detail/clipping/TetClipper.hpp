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
 * @brief Geometry clipping operations for a single tetrahedron geometry.
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
   *
   * \c kGeom.asHierarchy() must contain the following data:
   * - v0, v1, v2, v3: each contains a 3D coordinates of the
   *   tetrahedron vertices, in the order used by primal::Tetrahedron.
   *   The tet may be degenerate, but not inverted (negative volume).
   * - "fixOrientation": Whether to fix inverted tetrahedra
   *   instead of aborting.
   */
  TetClipper(const klee::Geometry& kGeom, const std::string& name = "");

  virtual ~TetClipper() = default;

  const std::string& name() const override { return m_name; }

  bool labelCellsInOut(quest::experimental::ShapeMesh& shapeMesh,
                       axom::Array<LabelType>& cellLabels) override;

  bool labelTetsInOut(quest::experimental::ShapeMesh& shapeMesh,
                      axom::ArrayView<const axom::IndexType> cellIds,
                      axom::Array<LabelType>& tetLabels) override;

  /*!
   * @copydoc MeshClipperStrategy::getGeometryAsTets()
   *
   * \c tets will have length one, because the geometry for this
   * class is a single tetrahedron.
   */
  bool getGeometryAsTets(quest::experimental::ShapeMesh& shapeMesh,
                         axom::Array<TetrahedronType>& tets) override;

#if !defined(__CUDACC__)
private:
#endif
  std::string m_name;

  //!@brief Tetrahedron before external transformation.
  TetrahedronType m_tetBeforeTrans;

  //!@brief Tetrahedron after external transformation.
  TetrahedronType m_tet;

  //!@brief External transformation.
  axom::primal::experimental::CoordinateTransformer<double> m_extTransformer;

  //!@brief Transformation of m_tet to unit tet.
  axom::primal::experimental::CoordinateTransformer<double> m_toUnitTet;

  template <typename ExecSpace>
  void labelCellsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                           axom::ArrayView<LabelType> cellLabel);

  template <typename ExecSpace>
  void labelTetsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                          axom::ArrayView<const axom::IndexType> cellsOnBdry,
                          axom::ArrayView<LabelType> tetLabels);

  // Extract clipper info from MeshClipperStrategy::m_info.
  void extractClipperInfo();
};

}  // namespace experimental
}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_TETCLIPPER_HPP
