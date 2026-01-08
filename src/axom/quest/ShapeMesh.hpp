// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_SHAPEMESH_HPP
#define AXOM_QUEST_SHAPEMESH_HPP

#include "axom/config.hpp"

#ifndef AXOM_USE_CONDUIT
  #error "ShapeMesh requires Conduit"
#endif

#ifndef AXOM_USE_SIDRE
  #error "ShapeMesh requires sidre"
  // Note: We guard sidre use for mesh stored in sidre, but sidre::ConduitMemory
  // is required even when the mesh is stored in Conduit.  Hence the dependence
  // on sidre.
#endif

#include "axom/core.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"
#include "axom/primal/geometry/Hexahedron.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"

#include "axom/sidre.hpp"

#include "conduit/conduit_node.hpp"
#include "conduit_blueprint.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{

/*!
 * @brief Computational mesh and intermediate data typically used in shaping.
 *
 * This wrapper class:
 * - encapsulate mesh-dependent data
 * - provide a common interface to mesh data regardless of the mesh format
 * - avoid redundant work
 * - provides some convenience tools typically used in shaping
 *
 * The wrapped mesh must have an unstructured 3D hex topology.  That
 * is the only topology currently supported.  It can be extended to
 * support 2D.
 *
 * TODO: Support MFEM mesh.  First pass only supports blueprint.
*/
class ShapeMesh
{
public:
  using RuntimePolicy = axom::runtime_policy::Policy;
  using Point3DType = primal::Point<double, 3>;
  using TetrahedronType = primal::Tetrahedron<double, 3>;
  using HexahedronType = primal::Hexahedron<double, 3>;
  using Plane3DType = axom::primal::Plane<double, 3>;
  using BoundingBox3DType = primal::BoundingBox<double, 3>;

  /*!
   * @brief Number of tetrahedra that a hexahedron decomposes into
   * @see hexToTets()
   *
   * @internal Only values of 24 and 18 are valid.  18 is likely more
   * performant because it generates fewer tets.
   *
   * @internal The hexToTets() method branches on the value of
   * NUM_TETS_PER_HEX at a low level, but the branching should be
   * optimized out by the compiler.
   */
  static constexpr axom::IndexType NUM_TETS_PER_HEX = 18;

  //!@brief Number of vertices per cell.
  static constexpr axom::IndexType NUM_VERTS_PER_CELL_3D = 8;

  /*!
   * @brief Constructor with computational mesh in a conduit::Node.
   *
   * @param [in] runtimePolicy The run-time policy used for computations
   * @param [in] allocatorId Allocator id for internal and scratch space.
   *   It should be compatible with @c runtimePolicy.
   *   If axom::INVALID_ALLOCATOR_ID is specified, it will be
   *   replaced with the default allocator for @c runtimePolicy.
   *   For good performance, especially on GPUs, fast memory pools
   *   should be used.
   * @param [in/out] bpMesh Blueprint mesh to shape into.
   * @param [in] topoName Name of the Blueprint topology.  If empty,
   *   use the first topology in @c bpMesh.
   * @param [in] matsetName Name of the Blueprint material set.
   *   If empty, use the first material set in @c bpMesh.
   *
   * Incoming mesh array data are assumed to be accessible by the
   * runtime policy, or an error will result.  (However the data need
   * not correspond to the allocator id.)
   */
  ShapeMesh(RuntimePolicy runtimePolicy,
            int allocatorId,
            conduit::Node& bpMesh,
            const std::string& topoName = {},
            const std::string& matsetName = {});

#ifdef AXOM_USE_SIDRE
  /*!
   * @brief Constructor with computational mesh in a sidre::Group.

   * @param [in] runtimePolicy The run-time policy used for computations
   * @param [in] allocatorId Allocator id for internal and scratch space.
   *   It should be compatible with @c runtimePolicy.
   *   If axom::INVALID_ALLOCATOR_ID is specified, it will be
   *   replaced with the default allocator for @c runtimePolicy.
   *   For good performance, especially on GPUs, fast memory pools
   *   should be used.
   * @param [in/out] bpMesh Blueprint mesh to shape into.
   * @param [in] topoName Name of the Blueprint topology.  If empty,
   *   use the first topology in @c bpMesh.
   * @param [in] matsetName Name of the Blueprint material set.
   *   If empty, use the first material set in @c bpMesh.
   *
   * Incoming mesh array data are assumed to be accessible by the
   * runtime policy, or an error will result.  (However the data need
   * not correspond to the allocator id.)
  */
  ShapeMesh(RuntimePolicy runtimePolicy,
            int allocatorId,
            sidre::Group* bpMesh,
            const std::string& topoName = {},
            const std::string& matsetName = {});
#endif

  /*!
   * @brief Runtime policy set in constructor.
   */
  RuntimePolicy getRuntimePolicy() const { return m_runtimePolicy; }

  /*!
   * @brief Allocator id set in constructor.
   */
  int getAllocatorID() const { return m_allocId; }

  /*
   * !@brief Return computational mesh as a sidre::Group if it has
   * that form, or nullptr otherwise.
  */
  sidre::Group* getMeshAsSidre() { return m_bpGrpExt; }

  /*!
   * @brief Return computational mesh as a conduit::Node if it has
   * that form, or nullptr otherwise.
  */
  conduit::Node* getMeshAsConduit() { return m_bpNodeExt; }

  //!@brief Dimension of the mesh (2 or 3)
  int dimension() const { return m_dim; }

  //!@brief Number of cells in mesh.
  IndexType getCellCount() const { return m_cellCount; }

  //!@brief Number of vertices in mesh.
  IndexType getVertexCount() const { return m_vertexCount; }

  //!@brief Set the threshold to snapping vertex coordinates near
  // zero to zero.  Default threshold is 1e-10.
  void setZeroThreshold(double threshold) { m_zeroThreshold = threshold; }

  /*!
   * @brief Decompose a hexahedron into NUM_TETS_PER_HEX tetrahedra.
   * @param hex [in] The hexahedron
   * @param tets [out] Pointer to space for NUM_TETS_PER_HEX tetrahedra.
   *
   * To avoid ambiguity due to the choice of 2 diagonals for dividing
   * each hex face into 2 triangles, we introduce a face-centered
   * point at the average of the face vertices and decompose the face
   * into 4 triangles.
   *
   * It is expected that this method will be used in long inner loops,
   * so it is bare-bones for best performance.  Caller must ensure
   * tets points to at least NUM_TETS_PER_HEX objects. This method
   * neither checks the pointer nor reallocates the space.
   */
  AXOM_HOST_DEVICE inline static void hexToTets(const HexahedronType& hex, TetrahedronType* tets);

  //@{
  //!@name Accessors to mesh data.
  //@}

  //@{
  //!@name Accessors to mesh-dependent intermediate data.
  /*
   * This data is dynamically generated as needed, and cached for use
   * by multiple geometry clippers.  The idea is to eliminate redundant
   * code and computations.
   */
  /*!
   * @brief Tetrahedral version of mesh cells with cell \c i having tet ids in
   * [i*NUM_TETS_PER_HEX, (i+1)*NUM_TETS_PER_HEX).
   */
  axom::ArrayView<const TetrahedronType> getCellsAsTets();
  axom::ArrayView<const HexahedronType> getCellsAsHexes();
  //!@brief Get volume of mesh cells.
  axom::ArrayView<const double> getCellVolumes();
  //!@brief Get volumes of tets in getCellsAsTets().
  axom::ArrayView<const double> getTetVolumes();
  //!@brief Get characteristic lengths of mesh cells.
  axom::ArrayView<const double> getCellLengths();
  axom::ArrayView<const BoundingBox3DType> getCellBoundingBoxes();
  axom::ArrayView<const IndexType, 2> getCellNodeConnectivity();
  axom::ArrayView<const Point3DType> getVertexPoints();
  const axom::StackArray<axom::ArrayView<const double>, 3>& getVertexCoords3D() const
  {
    return m_vertCoordsViews3D;
  }
  /*!
   * @brief Precompute (and cache) mesh-dependent intermediate data
   * that may be used in clipping.
   *
   * Mesh-dependent data are computed as needed and cached, but this
   * computes them all at once.
  */
  void precomputeMeshData();
  //@}

  /*!
   * @brief Check whether mesh meets requirements for shaping.
   * @param whyNot [out] Diagnostic message if mesh is invalid.
   *
   * Requirements for the mesh are:
   * - Follow blueprint conventions.  See
   *   https://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html
   * - Be unstructured.
   * - Have hexahedral elements.
   */
  bool isValidForShaping(std::string& whyNot) const;

  //@{
  /*!
   * @brief Create (Blueprint) matset in the mesh for a material.
   * @param materialName [in] Name of material
   * @param volumes [in] Cell-centered volumes
   * @param isFraction [in] Whether @c volumes is actually
   *   volume fractions.
   *
   * @pre volumes.size() == getCellCount()
   * @pre Mesh's matsets is multi-buffer, material-dominant
   * form (see https://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html#material-sets).
   * Currently not supporting other matset forms.
   */
  void setMatsetFromVolume(const std::string& materialName,
                           const axom::ArrayView<double>& volumes,
                           bool isFraction = false);

  /*!
   * @brief Compute and set the free volume fraction.
   * @param freeName [in] Name of free material.
   */
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
  //! @brief External computational mesh, when mesh is provided as a Group.
  axom::sidre::Group* m_bpGrpExt;

  //! @brief Initial shallow copy of m_bpGrp in an internal Node storage.
  conduit::Node m_bpNodeInt;
#endif

  //! @brief External computational mesh, when mesh is provided as a Node.
  conduit::Node* m_bpNodeExt {nullptr};

  //!@brief Dimension of mesh (2 or 3)
  int m_dim;

  //!@brief Number of cells in mesh.
  IndexType m_cellCount;

  //!@brief Number of vertices in mesh.
  IndexType m_vertexCount;

  //!@brief Threshold for snapping vertex coordinates to zero.
  double m_zeroThreshold;

  //!@brief 3D Vertex coordinates as 1D ArrayViews.
  axom::StackArray<axom::ArrayView<const double>, 3> m_vertCoordsViews3D;

  //!@brief 3D Vertex coordinates as 1D Points.
  axom::Array<Point3DType> m_vertPoints3D;

  //!@brief Vertex indices for each cell.
  axom::ArrayView<const axom::IndexType, 2> m_connectivity;

  //!@brief Mesh cells as an array of hexes.
  axom::Array<HexahedronType> m_cellsAsHexes;

  //!@brief Volumes of hex cells.
  axom::Array<double> m_hexVolumes;

  //!@brief Volumes of cell tets.
  axom::Array<double> m_tetVolumes;

  //!@brief Characteristic lengths of cells.
  axom::Array<double> m_cellLengths;

  //!@brief Bounding boxes for m_cellsAsHexes.
  axom::Array<BoundingBox3DType> m_hexBbs;

  //!@brief Cells as NUM_TETS_PER_HEX*m_cellCount tets.
  axom::Array<TetrahedronType> m_cellsAsTets;

  void computeCellsAsHexes();
  void computeCellsAsTets();
  void computeHexVolumes();
  void computeTetVolumes();
  void computeHexBbs();
  void computeCellLengths();
  void computeVertPoints();
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
  void computeTetVolumesImpl();

  template <typename ExecSpace>
  void computeHexBbsImpl();

  template <typename ExecSpace>
  void computeCellLengthsImpl();

  template <typename ExecSpace>
  void computeVertPointsImpl();

  template <typename ExecSpace, typename T>
  void elementwiseDivideImpl(const T* numerator, const T* denominator, T* quotient, axom::IndexType n);

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

  /*!
   * @brief Get a Conduit hierarchy within another Conduit hierarchy.
   *
   * If specified path doesn't exist, create it.  If it does, verify
   * its conduit::DataType.
   *
   * If the node holds array data, verify that the data is compatible
   * with m_runtimePolicy.
   */
  conduit::Node& getMeshConduitPath(conduit::Node& node,
                                    const std::string& path,
                                    const conduit::DataType& dtype);
};

AXOM_HOST_DEVICE inline void ShapeMesh::hexToTets(const HexahedronType& hex, TetrahedronType* tets)
{
  AXOM_STATIC_ASSERT(NUM_TETS_PER_HEX == 24 || NUM_TETS_PER_HEX == 18);

  if(NUM_TETS_PER_HEX == 24)
  {
    hex.triangulate(tets);
  }
  else
  {
    /*
     * Six tets sharing the line segment between hex vertices 4 and 2.
     * Each tet also uses 2 of the remaining 6 hex vertices (any 2
     * that shares a hex edge).
     */
    tets[0][0] = hex[4];
    tets[0][1] = hex[2];
    tets[0][2] = hex[1];
    tets[0][3] = hex[0];

    tets[1][0] = hex[4];
    tets[1][1] = hex[2];
    tets[1][2] = hex[0];
    tets[1][3] = hex[3];

    tets[2][0] = hex[4];
    tets[2][1] = hex[2];
    tets[2][2] = hex[3];
    tets[2][3] = hex[7];

    tets[3][0] = hex[4];
    tets[3][1] = hex[2];
    tets[3][2] = hex[7];
    tets[3][3] = hex[6];

    tets[4][0] = hex[4];
    tets[4][1] = hex[2];
    tets[4][2] = hex[6];
    tets[4][3] = hex[5];

    tets[5][0] = hex[4];
    tets[5][1] = hex[2];
    tets[5][2] = hex[5];
    tets[5][3] = hex[1];

    // Centroids of the 6 hex faces.
    Point3DType mp0473 = Point3DType::midpoint(Point3DType::midpoint(hex[0], hex[4]),
                                               Point3DType::midpoint(hex[7], hex[3]));
    Point3DType mp1562 = Point3DType::midpoint(Point3DType::midpoint(hex[1], hex[5]),
                                               Point3DType::midpoint(hex[6], hex[2]));
    Point3DType mp0451 = Point3DType::midpoint(Point3DType::midpoint(hex[0], hex[4]),
                                               Point3DType::midpoint(hex[5], hex[1]));
    Point3DType mp3762 = Point3DType::midpoint(Point3DType::midpoint(hex[3], hex[7]),
                                               Point3DType::midpoint(hex[6], hex[2]));
    Point3DType mp0123 = Point3DType::midpoint(Point3DType::midpoint(hex[0], hex[1]),
                                               Point3DType::midpoint(hex[2], hex[3]));
    Point3DType mp4567 = Point3DType::midpoint(Point3DType::midpoint(hex[4], hex[5]),
                                               Point3DType::midpoint(hex[6], hex[7]));

    /*
     * Tets from the 6 hex faces (two per face).  If the face is
     * coplanar, its 2 tets are degenerate.
     */
    tets[6][0] = hex[4];
    tets[6][1] = hex[6];
    tets[6][2] = hex[7];
    tets[6][3] = mp4567;
    tets[7][0] = hex[4];
    tets[7][1] = hex[5];
    tets[7][2] = hex[6];
    tets[7][3] = mp4567;

    tets[8][0] = hex[0];
    tets[8][1] = hex[2];
    tets[8][2] = hex[3];
    tets[8][3] = mp0123;
    tets[9][0] = hex[0];
    tets[9][1] = hex[1];
    tets[9][2] = hex[2];
    tets[9][3] = mp0123;

    tets[10][0] = hex[4];
    tets[10][1] = hex[0];
    tets[10][2] = hex[3];
    tets[10][3] = mp0473;
    tets[11][0] = hex[4];
    tets[11][1] = hex[3];
    tets[11][2] = hex[7];
    tets[11][3] = mp0473;

    tets[12][0] = hex[5];
    tets[12][1] = hex[1];
    tets[12][2] = hex[2];
    tets[12][3] = mp1562;
    tets[13][0] = hex[5];
    tets[13][1] = hex[1];
    tets[13][2] = hex[2];
    tets[13][3] = mp1562;

    tets[14][0] = hex[4];
    tets[14][1] = hex[5];
    tets[14][2] = hex[1];
    tets[14][3] = mp0451;
    tets[15][0] = hex[4];
    tets[15][1] = hex[1];
    tets[15][2] = hex[0];
    tets[15][3] = mp0451;

    tets[16][0] = hex[7];
    tets[16][1] = hex[6];
    tets[16][2] = hex[2];
    tets[16][3] = mp3762;
    tets[17][0] = hex[7];
    tets[17][1] = hex[6];
    tets[17][2] = hex[3];
    tets[17][3] = mp3762;
  }
}

}  // namespace experimental
}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_SHAPEMESH_HPP
