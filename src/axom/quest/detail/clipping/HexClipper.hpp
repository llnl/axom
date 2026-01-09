// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_HEXCLIPPER_HPP
#define AXOM_QUEST_HEXCLIPPER_HPP

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
 * @brief Geometry clipping operations for hexahedron geometries.
*/
class HexClipper : public MeshClipperStrategy
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
   * - v0, v1, v2, ..., v7: each contains a 3D coordinates of the
   *   hexahedron vertices, in the order used by primal::Hexahedron.
   *   The hex may be degenerate, but when subdivided into tetrahedra,
   *   none of them may be inverted (have negative volume).
   */
  HexClipper(const klee::Geometry& kGeom, const std::string& name = "");

  virtual ~HexClipper() = default;

  const std::string& name() const override { return m_name; }

  bool labelCellsInOut(quest::experimental::ShapeMesh& shappeMesh,
                       axom::Array<LabelType>& label) override;

  bool labelTetsInOut(quest::experimental::ShapeMesh& shapeMesh,
                      axom::ArrayView<const axom::IndexType> cellIds,
                      axom::Array<LabelType>& tetLabels) override;

  bool getGeometryAsTets(quest::experimental::ShapeMesh& shappeMesh,
                         axom::Array<TetrahedronType>& tets) override;

#if !defined(__CUDACC__)
private:
#endif
  std::string m_name;

  //!@brief Hexahedron before transformation.
  HexahedronType m_hexBeforeTrans;

  //!@brief Hexahedron after transformation.
  HexahedronType m_hex;

  //!@brief External transformation.
  axom::primal::experimental::CoordinateTransformer<double> m_extTransformer;

  //!@brief Bounding box of m_hex.
  BoundingBox3DType m_hexBb;

  //!@brief Tetrahedralized version of of m_hex.
  axom::Array<TetrahedronType> m_tets;

  //!@brief Triangles on the discretized hex surface, oriented inward.
  axom::StackArray<Triangle3DType, 24> m_surfaceTriangles;

  template <typename ExecSpace>
  void labelCellsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                           axom::ArrayView<LabelType> label);

  template <typename ExecSpace>
  void labelTetsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                          axom::ArrayView<const axom::IndexType> cellIds,
                          axom::ArrayView<LabelType> tetLabels);

  /*!
   * @brief Compute whether a polyhedron is inside, outside or on the boundary
   * of a hexahedron.
   *
   * @param verts [in] The polyhedron (either hex or tet).
   * @param vertsBb [in] Bounding box of @c verts
   * @param hexBb [in] Bounding box of the hex
   * @param hexTets [in] The hex, decomposed into tets.
   * @param surfaceTriangles [in] A copy of m_surfaceTriangles,
   *   but in the appropriate memory space (host or device).
   *
   * This method should be fast.  It may label something as on
   * the boundary when it is inside or outside; this is a conservative
   * error and the way we use it doesn't lead to real errors.
   */
  template <typename Polyhedron>
  AXOM_HOST_DEVICE inline LabelType polyhedronToLabel(
    const Polyhedron& verts,
    const BoundingBox3DType& vertsBb,
    const BoundingBox3DType& hexBb,
    const axom::ArrayView<const TetrahedronType>& hexTets,
    const axom::StackArray<Triangle3DType, 24>& surfaceTriangles) const;

  void extractClipperInfo();

  void computeSurface();
};

}  // namespace experimental
}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_HEXCLIPPER_HPP
