// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_TETMESHCLIPPER_HPP
#define AXOM_QUEST_TETMESHCLIPPER_HPP

#include "axom/klee/Geometry.hpp"
#include "axom/quest/MeshClipperStrategy.hpp"
#include "axom/primal/geometry/CoordinateTransformer.hpp"
#include "axom/spin/BVH.hpp"

// Implementation requires Conduit.
#include "conduit_blueprint.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{

/*!
 * @brief Geometry clipping operations for tetrahedral mesh geometries.
 */
class TetMeshClipper : public MeshClipperStrategy
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
   * - "klee::Geometry:tetMesh": A blueprint tetrahedral mesh.
   *   The tetrahedra may be degenerate, but not inverted (negative volume).
   * - "topologyName": The mesh's blueprint topology name
   * - "fixOrientation": Whether to fix inverted tetrahedra
   *   instead of aborting.
   */
  TetMeshClipper(const klee::Geometry& kGeom, const std::string& name = "");

  virtual ~TetMeshClipper() = default;

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

  //! @brief Topology to use in the Blueprint tet mesh.
  std::string m_topoName;

  //! @brief Coordset to use in the Blueprint tet mesh.
  std::string m_coordsetName;

  //! @brief Tet mesh in Blueprint format.
  conduit::Node m_tetMesh;

  //! @brief Bounding box of the tet mesh.
  axom::primal::BoundingBox<double, 3> m_tetMeshBb;

  //! @brief Number of tets in the tet mesh.
  axom::IndexType m_tetCount;

  /*!
   * @brief Combined external transformation.
   *
   * (TetMesh has no internal transformation.)
   */
  axom::primal::experimental::CoordinateTransformer<double> m_transformer;

  template <typename ExecSpace>
  void labelCellsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                           axom::ArrayView<LabelType> label);

  template <typename ExecSpace>
  void labelTetsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                          axom::ArrayView<const axom::IndexType> cellIds,
                          axom::ArrayView<LabelType> tetLabels);

  template <typename ExecSpace>
  void vertexInsideToCellLabel(quest::experimental::ShapeMesh& shapeMesh,
                               axom::ArrayView<bool>& vertIsInside,
                               axom::Array<LabelType>& labels);

  //! @brief Compute rays from an interior point of each mesh hex, pointing away from tet mesh centroid.
  template <typename ExecSpace>
  void computeHexRays(quest::experimental::ShapeMesh& shapeMesh, axom::Array<Ray3DType>& hexRays);

  //! @brief Compute rays from an interior point of each mesh tet, pointing away from tet mesh centroid.
  template <typename ExecSpace>
  void computeTetRays(quest::experimental::ShapeMesh& shapeMesh,
                      axom::ArrayView<const axom::IndexType> cellIds,
                      axom::Array<Ray3DType>& tetRays,
                      axom::Array<BoundingBox3DType>& tetBbs);

  // Extract clipper info from MeshClipperStrategy::m_info.
  void extractClipperInfo();

  // Check validity of tetMesh for our purposes.
  bool isValidTetMesh(const conduit::Node& tetMesh, std::string& whyBad) const;

  /*!
   * @brief Add a transformed coordset to m_tetMesh.
   *
   * The transformed version @c m_tetMesh["coordsets"][m_coordsetName],
   * transformed through m_transformer.  It has the name
   * m_coordsetName + ".trans".
   */
  void transformCoordset();

  template <typename ExecSpace>
  void computeTets(axom::ArrayView<TetrahedronType> tetsView);

  //@{
  //!@name For computing surface of m_tetMesh.
  /*!
   * @brief Compute the surface triangles and their BVH.
   */
  template <typename ExecSpace>
  void computeSurfaceTrianglesAndBVH(int allocId,
                                     axom::Array<Triangle3DType>& surfTris,
                                     spin::BVH<3, ExecSpace, double>& bvh);

  /*!
   * @brief Compute the tet-mesh geometry surface as trianglular facets.
   */
  template <typename ExecSpace>
  axom::Array<Triangle3DType> computeGeometrySurface(int allocId);

  /*!
   * @brief Add a polyhedral topology to an unstructured tet mesh.
   * @param tetMesh Input unstructured tet mesh, single domain.
   * @param polyTopo Output unstructured polyhedral topology.
   */
  template <typename ExecSpace>
  void make_polyhedral_topology(conduit::Node& tetTopo, conduit::Node& polyTopo);

  //!@brief Write out for debugging
  void writeTrianglesToVTK(const axom::Array<Triangle3DType>& triangles, const std::string& filename);
  //@}

  void copy_node_with_array_allocator(conduit::Node& src,
                                      conduit::Node& dst,
                                      conduit::index_t conduitAllocId);
  void copy_hierarchy_with_array_allocator(conduit::Node& hierarchy,
                                           const std::string& srcPath,
                                           const std::string& dstPath,
                                           int allocId);
  void copy_topo_and_coords_to(int allocId);
};

}  // namespace experimental
}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_TETMESHCLIPPER_HPP
