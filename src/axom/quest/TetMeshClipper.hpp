// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_TETMESHCLIPPER_HPP
#define AXOM_QUEST_TETMESHCLIPPER_HPP

#include "axom/klee/Geometry.hpp"
#include "axom/quest/GeometryClipperStrategy.hpp"
#include "axom/primal/geometry/CoordinateTransformer.hpp"

// Implementation requires Conduit.
#include "conduit_blueprint.hpp"

namespace axom
{
namespace quest
{

/*!
  @brief Geometry clipping operations for tetrahedral mesh geometries.

  @internal TODO: Implement load balancing.  The 1D array of shapee hexes
  should be load balanced for better performance.
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

  bool getGeometryAsTets(quest::ShapeeMesh& shappeMesh, axom::Array<TetrahedronType>& tets) override;

#if !defined(__CUDACC__)
private:
#endif
  std::string m_name;

  //! @brief Topology to use in the Blueprint tet mesh.
  std::string m_topoName;

  //! @brief Coordset to use in the Blueprint tet mesh.
  std::string m_coordsetName;

  //! @brief Tet mesh in Blueprint format.
  conduit::Node m_bpMesh;

  //! @brief Bounding box of the tet mesh.
  axom::primal::BoundingBox<double, 3> m_tetMeshBb;

  //! @brief Number of tets in the tet mesh.
  axom::IndexType m_tetCount;

  //! @brief Geometry as tetrahedra.
  axom::Array<TetrahedronType> m_tets;

  axom::primal::CoordinateTransformer<double> m_transformer;

  template <typename ExecSpace>
  void labelInOutImpl(quest::ShapeeMesh& shapeeMesh, axom::Array<char>& label);

  template <typename ExecSpace>
  void vertexInsideToCellLabel(
    quest::ShapeeMesh& shapeeMesh,
    axom::ArrayView<bool>& vertIsInside,
    axom::Array<LabelType>& labels);

  // Extract clipper info from GeometryClipperStrategy::m_info.
  void extractClipperInfo();

  void transformCoordset();

  void computeTets();

  //@{
  //!@name For computing surface of m_bpMesh.
  /*!
    @brief Entry point for computing geometry surface.

    This computation is independent of the shapee mesh, except that we
    need the policy and allocator id.
  */
  axom::Array<Triangle3DType> computeGeometrySurface(axom::runtime_policy::Policy policy, int allocId);
  template <typename ExecSpace>

  axom::Array<Triangle3DType> computeGeometrySurface(int allocId);

  //!@brief Add a polyhedral topology to an unstructured tet mesh.
  template <typename ExecSpace>
  void make_polyhedral_topology(const conduit::Node& inputBp,
                                conduit::Node& tetTopo,
                                conduit::Node& polyTopo);

  //!@brief Write out for debugging
  void writeTrianglesToVTK(
    const axom::Array<Triangle3DType>& triangles,
    const std::string& filename);
  //@}

  void copy_node_with_array_allocator(conduit::Node& src,
                                      conduit::Node& dst,
                                      conduit::index_t conduitAllocId);
  void copy_hierarchy_with_array_allocator(conduit::Node& hierarchy,
                                           const std::string& srcPath,
                                           const std::string& dstPath,
                                           int allocId);
  void copy_tetmesh_arrays_to(int allocId);
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_TETMESHCLIPPER_HPP
