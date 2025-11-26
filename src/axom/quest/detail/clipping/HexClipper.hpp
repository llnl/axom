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
 * @brief Geometry clipping operations for sphere geometries.
*/
class HexClipper : public MeshClipperStrategy
{
public:
  /*!
   * @brief Constructor.

   * @param [in] kGeom Describes the shape to place
   *   into the mesh.
   * @param [in] name To override the default strategy name
  */
  HexClipper(const klee::Geometry& kGeom, const std::string& name = "");

  virtual ~HexClipper() = default;

  const std::string& name() const override { return m_name; }

  /*!
   * If a mesh cell has all vertices outside the geometry, it labeled outside.
   * This will miss cases where an edge of the cell passes through the geometry.
  */
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

  //!@brief Bounding box of m_hex.
  BoundingBox3DType m_hexBb;

  //!@brief Tetrahedralized version of of m_hex.
  axom::Array<TetrahedronType> m_tets;

  //!@brief Triangles on the discretized hex surface, oriented inward.
  axom::StackArray<Triangle3DType, 24> m_surfaceTriangles;

  axom::primal::experimental::CoordinateTransformer<double> m_transformer;

  template <typename ExecSpace>
  void labelCellsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                           axom::ArrayView<LabelType> label);

  template <typename ExecSpace>
  void labelTetsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                          axom::ArrayView<const axom::IndexType> cellIds,
                          axom::ArrayView<LabelType> tetLabels);

  //!@brief Compute LabelType for a polyhedron (hex or tet in our case).
  template <typename Polyhedron>
  AXOM_HOST_DEVICE inline
  LabelType polyhedronToLabel(
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
