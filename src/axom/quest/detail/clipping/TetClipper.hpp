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
   */
  TetClipper(const klee::Geometry& kGeom, const std::string& name = "");

  virtual ~TetClipper() = default;

  const std::string& name() const override { return m_name; }

  /*!
   * @copydoc MeshClipperStrategy::getGeometryAsTets()
   *
   * \c tets will have length one, because the geometry for this
   * class is a single tetrahedron.
   */
  bool getGeometryAsTets(quest::experimental::ShapeMesh& shappeMesh,
                         axom::Array<TetrahedronType>& tets) override;

#if !defined(__CUDACC__)
private:
#endif
  std::string m_name;

  //!@brief Tetrahedron before transformation.
  TetrahedronType m_tetBeforeTrans;

  //!@brief Tetrahedron after transformation.
  TetrahedronType m_tet;

  axom::primal::BoundingBox<double, 3> m_bb;

  axom::primal::experimental::CoordinateTransformer<double> m_transformer;

  // Extract clipper info from MeshClipperStrategy::m_info.
  void extractClipperInfo();
};

}  // namespace experimental
}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_TETCLIPPER_HPP
