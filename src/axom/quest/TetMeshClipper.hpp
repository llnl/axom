// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_TETMESHCLIPPER_HPP
#define AXOM_QUEST_TETMESHCLIPPER_HPP

#include "axom/klee/Geometry.hpp"
#include "axom/quest/GeometryClipperStrategy.hpp"

namespace axom
{
namespace quest
{

/*!
  @brief GeometryClipper specialized for tetrahedral mesh geometries.
*/
class TetMeshClipper : public GeometryClipperStrategy
{
public:
  /*!
    @brief Constructor.

    @param [in] kGeom Describes the shape to place
      into the mesh.
    @param [in] name To override the default strategy name
  */
  TetMeshClipper(const klee::Geometry& kGeom, const std::string& name = "");

  virtual ~TetMeshClipper() = default;

  const std::string& name() const override { return m_name; }

  bool labelInOut(quest::ShapeeMesh& shappeMesh, axom::Array<char>& label) override;

  bool getShapeAsTets(quest::ShapeeMesh& shappeMesh, axom::Array<TetrahedronType>& tets) override;

#if !defined(__CUDACC__)
private:
#endif
  std::string m_name;

  //! @brief Topology to use in the Blueprint tet mesh.
  std::string m_topoName;

  //! @brief Tet mesh in Blueprint format.
  conduit::Node* m_bpMesh;

  //! @brief Number of cells in the tet mesh.
  axom::IndexType m_cellCount;

  template <typename ExecSpace>
  void labelInOutImpl(quest::ShapeeMesh& shapeeMesh, axom::Array<char>& label);

  // Extract clipper info from GeometryClipperStrategy::m_info.
  void extractClipperInfo();
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_TETMESHCLIPPER_HPP
