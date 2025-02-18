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

#ifdef AXOM_USE_SIDRE
#include "axom/sidre.hpp"
#endif

#include "conduit/conduit_node.hpp"
#include "conduit_blueprint.hpp"

namespace axom
{
namespace quest
{

/*!
  @brief Computational mesh and intermediate data typically used in shaping.

  The purpose of this class is to encapsulate mesh-dependent data and
  avoid redundant work.  Also to provide some convenience tools typically
  used in shaping.

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
    @brief Constructor with computational mesh in a conduit::Node.

    @param [in] runtimePolicy
    @param [in] allocatorId Allocator id for internal and scratch space.
    @param [in/out] bpMesh Blueprint mesh to shape into.
    @param [in] topoName Name of the Blueprint topology.  If empty,
      use the first topology in @c bpMesh.
    @param [in] matsetName Name of the Blueprint material set.
      If empty, use the first material set in @c bpMesh.

    It is an error if allocator id is not usable with the runtime policy.
    If @c allocatorId is axom::INVALID_ALLOCATOR_ID, the default
    allocator for the runtime policy will be used.

    Mesh array data are assumed to be accessible by the runtime policy,
    but they need not correspond to the allocator id.

    Because @c conduit::Node doesn't support application-specified
    allocator id for (only) arrays, the incoming @c bpNode must have
    all arrays pre-allocated in a space accessible by the runtime
    policy.  Any needed-but-missing space would lead to an exception.
  */
  ShapeeMesh(RuntimePolicy runtimePolicy,
             int allocatorId,
             conduit::Node& bpMesh,
             const std::string& topoName = {},
             const std::string& matsetName = {});

#ifdef AXOM_USE_SIDRE
  /*!
    @brief Constructor with computational mesh in a sidre::Group.

    @param [in] runtimePolicy
    @param [in] allocatorId Allocator id for internal and scratch space.
    @param [in/out] bpMesh Blueprint mesh to shape into.
    @param [in] topoName Name of the Blueprint topology.  If empty,
      use the first topology in @c bpMesh.
    @param [in] matsetName Name of the Blueprint material set.
      If empty, use the first material set in @c bpMesh.

    It is an error if allocator id is not usable with the runtime policy.
    If @c allocatorId is axom::INVALID_ALLOCATOR_ID, the default
    allocator for the runtime policy will be used.

    Mesh array data are assumed to be accessible by the runtime policy,
    but they need not correspond to the allocator id.
  */
  ShapeeMesh(RuntimePolicy runtimePolicy,
             int allocatorId,
             sidre::Group* bpMesh,
             const std::string& topoName = {},
             const std::string& matsetName = {});
#endif

  // TODO: Support other mesh forms: Blueprint Group, MFEM.

  /*!
    @brief Runtime policy set in constructor.

    getAllocatorId() and getRuntimePolicy() are guaranteed to be
    compatible.
  */
  RuntimePolicy getRuntimePolicy() const { return m_runtimePolicy; }

  /*!
    @brief Allocator id set in constructor.

    getAllocatorId() and getRuntimePolicy() are guaranteed to be
    compatible.
  */
  int getAllocatorId() const { return m_allocId; }

  /*
    !@brief Return computational mesh as a sidre::Group if it has
    that form, or nullptr otherwise.
  */
  sidre::Group* getMeshAsSidre() { return m_bpGrpExt; }

  /*!
    @brief Return computational mesh as a conduit::Node if it has
    that form, or nullptr otherwise.
  */
  conduit::Node* getMeshAsConduit() { return m_bpNodeExt; }

  //!@brief Dimension of the mesh (2 or 3)
  int dimension() const { return m_dim; }

  //!@brief Number of cells in mesh.
  IndexType getCellCount() const { return m_cellCount; }

  //!@brief Number of vertices in mesh.
  IndexType getVertexCount() const { return m_vertexCount; }

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
  {
    return m_vertCoordsViews3D;
  }
  //@}

  /*!
    @brief Check whether mesh meets requirements for shaping.
    @param whyNot [out] Diagnostic message if mesh is invalid.

    Note: The check is only performed on Blueprint meshes, not
    on MFEM meshes.
  */
  bool isValidForShaping(std::string& whyNot) const;

  //@{
  /*!
    @brief Create (Blueprint) matset in the mesh for a material.
    @param materialName [in]
    @param volumes [in] Cell-centered volumes
    @param isFraction [in] Whether @c volumes is actually
      volume fractions.

    @pre volumes.size() == getCellCount()
    @pre Mesh's matsets is multi-buffer, material-dominant
    form (see https://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html#material-sets).
    Currently not supporting other matset forms.
  */
  void setMatsetFromVolume(const std::string& materialName,
                           const axom::ArrayView<double>& volumes,
                           bool isFraction = false);

  void setFreeVolumeFractions(const std::string& freeName);
  //@}

private:
  const RuntimePolicy m_runtimePolicy;

  int m_allocId;

  //! @brief Mesh topology name.
  const std::string m_topoName;

  //! @brief Mesh matset name.
  const std::string m_matsetName;

#if defined(AXOM_USE_SIDRE)
  //! @brief External computational mesh, when provided as a Group.
  axom::sidre::Group* m_bpGrpExt;

  //! @brief Initial shallow copy of m_bpGrp in an internal Node storage.
  conduit::Node m_bpNodeInt;
#endif

  //! @brief External computational mesh, when provided as a Node.
  conduit::Node* m_bpNodeExt {nullptr};

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

  template <typename ExecSpace, typename T>
  void elementwiseDivideImpl(const T* numerator,
                             const T* denominator,
                             T* quotient,
                             axom::IndexType n);

  template <typename T>
  void fillNImpl(axom::ArrayView<T> a, const T& val) const;

  template <typename T>
  void elementwiseAddImpl(const axom::ArrayView<T> a,
                          const axom::ArrayView<T> b,
                          axom::ArrayView<T> result) const;

  template <typename T>
  void elementwiseComplementImpl(const axom::ArrayView<T> a,
                                 const T& val,
                                 axom::ArrayView<T> results) const;
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_SHAPEEMESH_HPP
