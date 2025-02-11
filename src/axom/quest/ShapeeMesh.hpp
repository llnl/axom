// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_SHAPEEMESH_HPP
#define AXOM_QUEST_SHAPEEMESH_HPP

#ifndef AXOM_USE_CONDUIT
  #error "ShapeeMesh requires Conduit"
// TODO: Support MFEM and Sidre blueprint as well.
#endif

#include "axom/core.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"
#include "axom/primal/geometry/Hexahedron.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"

#include "conduit/conduit_node.hpp"
#include "conduit_blueprint.hpp"

namespace axom
{
namespace quest
{

/*!
  @brief Computational mesh and intermediate data typically used in shaping.

  The purpose of this class is to encapsulate mesh-dependent data and
  avoid redundant work.

  The mesh must have an unstructured 3D hex topology.  That is the only
  topology currently supported.  It can be extended to support 2D.

  TODO: Support sidre::Group blueprint and MFEM mesh.  First pass only
  supports Conduit blueprint.
*/
class ShapeeMesh
{
public:
  using RuntimePolicy = axom::runtime_policy::Policy;
  using TetrahedronType = primal::Tetrahedron<double, 3>;
  using HexahedronType = primal::Hexahedron<double, 3>;
  using BoundingBox3dType = primal::BoundingBox<double, 3>;

  /*!
    @brief Constructor

    @param [in] runtimePolicy
    @param [in] allocatorId Allocator id for internal memory.
    @param [in/out] bpMesh Blueprint mesh to shape into.
    @param [in] topo Name of the Blueprint topology.  If empty,
      use the first topology in @c bpMesh.

    It is an error if allocator id is not usable with the runtime policy.
    If @c allocatorId is axom::INVALID_ALLOCATOR_ID, the default
    allocator for the runtime policy will be used.

    Mesh array data are assumed to be accessible by the runtime policy,
    but they need not correspond to the allocator id.
  */
  ShapeeMesh(RuntimePolicy runtimePolicy,
             int allocatorId,
             conduit::Node& bpMesh,
             const std::string& topo = {});

  // TODO: Support other mesh forms: Blueprint Group, MFEM.

  //!@brief The mesh domain.
  conduit::Node& getMeshDomain(int domainLocalId);

  //!@brief Dimension of the mesh (2 or 3)
  int dimension() const { return m_dim; }

  //!@brief Number of cells in mesh.
  IndexType getCellCount() const { return m_cellCount; }

  //!@brief Number of vertices in mesh.
  IndexType getVertexCount() const { return m_vertexCount; }

  /*!
    @brief Expose runtime policy for related code to use.

    getAllocatorId() and getRuntimePolicy() are guaranteed to be
    compatible.
  */
  RuntimePolicy getRuntimePolicy() const { return m_runtimePolicy; }

  /*!
    @brief Expose allocator id for related code to use.

    getAllocatorId() and getRuntimePolicy() are guaranteed to be
    compatible.
  */
  int getAllocatorId() const { return m_allocId; }

  //@{
  //!@name Accessors to mesh data.
  //@}

  //@{
  //!@name Accessors to mesh-dependent intermediate data.
  // This data is dynamically generated as needed, and cached.
  axom::ArrayView<const TetrahedronType> getCellsAsTets();
  axom::ArrayView<const HexahedronType> getCellsAsHexes();
  axom::ArrayView<const double> getCellVolumes();
  axom::ArrayView<const BoundingBox3dType> getCellBoundingBoxes();
  axom::ArrayView<const IndexType, 2> getConnectivity();
  const axom::StackArray<axom::ArrayView<const double>, 3>& getVertexCoords3D() const
    { return m_vertCoordsViews3D; }
  //@}

private:
  const RuntimePolicy m_runtimePolicy;

  int m_allocId;

  //! @brief Mesh represented as Sidre Group
  const std::string m_bpTopo;

  //! @brief Mesh in an external Node, when provided as a Node.
  conduit::Node* m_bpNodeExt {nullptr};

  //! @brief Initial copy of mesh in an internal Node storage.
  conduit::Node m_bpNodeInt;

  //!@brief Dimension of mesh (2 or 3)
  int m_dim;

  //!@brief Number of cells in mesh.
  IndexType m_cellCount;

  //!@brief Number of vertices in mesh.
  IndexType m_vertexCount;

  //!@brief 3D Vertex coordinates as 1D ArrayViews.
  axom::StackArray<axom::ArrayView<const double>, 3> m_vertCoordsViews3D;

  //!@brief Vertex indices for each cell.
  axom::ArrayView<const axom::IndexType, 2> m_connectivity;

  //!@brief Mesh cells as an array of hexes.
  axom::Array<HexahedronType> m_cellsAsHexes;

  //!@brief m_cellsAsHexes, each split into 24 tets.
  axom::Array<TetrahedronType> m_cellsAsTets;

  //!@brief Bounding boxes for m_cellsAsHexes.
  axom::Array<BoundingBox3dType> m_hexBbs;

  //!@brief Volumes of hex cells.
  axom::Array<double> m_hexVolumes;

  void computeCellsAsHexes();
  void computeCellsAsTets();
  void computeHexVolumes();
  void computeHexBbs();
  void computeConnectivity();

#if defined(__CUDACC__)
public:
#endif

  template <typename ExecSpace>
  void computeCellsAsHexesImpl();

  template <typename ExecSpace>
  void computeCellsAsTetsImpl();

  template <typename ExecSpace>
  void computeHexVolumesImpl();

  template <typename ExecSpace>
  void computeHexBbsImpl();
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_SHAPEEMESH_HPP
