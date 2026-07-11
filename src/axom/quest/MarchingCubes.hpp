// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/*!
 * @file MarchingCubes.hpp
 *
 * @brief Consists of classes implementing marching cubes algorithm to
 * compute isocontour from a scalar field in a blueprint mesh.
 */

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifdef AXOM_USE_CONDUIT

  // Axom includes
  #include "axom/core/execution/runtime_policy.hpp"
  #include "axom/mint/mesh/UnstructuredMesh.hpp"

  // Conduit includes
  #include "conduit_node.hpp"

  // C++ includes
  #include <string>

namespace axom::quest
{
namespace detail::marching_cubes
{
class MarchingCubesSingleDomain;
}  // namespace detail::marching_cubes

/*!
 * @brief Enum for the legacy marching cubes data-parallel implementation.
 *
 * Partial parallel implementation uses a non-parallizable loop and
 * processes less data.  It has been shown to work well on CPUs.
 * Full parallel implementation processes more data, but parallelizes
 * fully and has been shown to work well on GPUs.  byPolicy chooses
 * based on runtime policy.
 *
 * This setting controls only the legacy structured-mesh backend.
 * When MarchingCubes is configured to use the bump::extraction::CutField backend,
 * bump manages its own internal parallelism and this value is accepted only for API compatibility.
 */
enum class MarchingCubesDataParallelism
{
  byPolicy = 0,
  hybridParallel = 1,
  fullParallel = 2
};

/*!
 * @brief Enum controlling the meaning of the parent-cell ids reported for generated contour facets
 * (see MarchingCubes::getContourFacetParents and MarchingCubes::populateContourMesh).
 *
 * The legacy marching cubes implementation numbered parent cells by their flat
 * index in the same row- or column-major ordering as the input scalar function array
 * (i.e. following the function field's stride order).
 * The bump-backed implementation natively numbers cells by their Blueprint zone index,
 * which uses a fixed i-fastest ordering independent of how the field is stored in memory.
 * For structured input these two numberings coincide only when the field is stored i-fastest;
  *otherwise they differ by a stride-order permutation.
 *
 * This enum lets callers choose which numbering they receive:
 *  - \c blueprintZoneId (default): report the Blueprint zone index.  This is the natural,
 *    mesh-type-agnostic identifier and the only meaningful choice for unstructured input.
 *  - \c legacyFieldOrder: reproduce the legacy numbering (flat index in the function field's stride order).
 *     Provided so existing structured-mesh callers that depend on the historical meaning are unaffected.
 *    This option only applies to structured input; for unstructured input the Blueprint zone id is always used.
 */
enum class MarchingCubesParentCellIdMode
{
  blueprintZoneId = 0,
  legacyFieldOrder = 1
};

/*!
 * @brief Enum selecting the isosurface case-table / intersector robustness used by the bump backend
 *
 * The bump CutField backend determines per-cell topology with an intersector policy plus
 * VisIt-derived cut tables.  The default intersector (\c axom::bump::extraction::FieldIntersector)
 * classifies cell corners with a strict two-label test (corner value > isovalue),
 * evaluates edge crossings in single precision (its \c FieldType is \c float),
 * and uses a single fixed triangulation per case. Like the classic 1987 marching-cubes tables,
 * this resolves ambiguous (saddle) configurations consistently but not necessarily in a way
 * that matches the trilinear interpolant. It does not implement the +/-/0 (three-label) / asymptotic-decider
 * topology of Wenger's Isosurfaces or MC33.
 *
 * This enum is in anticipation of the more robust case that will be added soon and only applies to the
 * new bump-based backend:
 *  - \c standard (default): use bump's default intersector + tables.
 *    This is the only policy currently implemented.
 *  - \c robust: request a topologically-robust intersector/table set (double precision, +/-/0 aware).
 *    When bump provides such a policy it will be selected here with no further change to quest;
 *    until then, selecting \c robust behaves identically to \c standard (and may emit a one-time
 *    informational note), so callers can opt in now and benefit automatically once the robust policy lands.
 */
enum class MarchingCubesRobustnessPolicy
{
  standard = 0,
  robust = 1
};

/*!
 * @brief Class implementing marching cubes to compute a contour
 * mesh from a scalar function on an input mesh.
 *
 * This implementation is for the original 1987 algorithm:
 * Lorensen, William E.; Cline, Harvey E. (1 August 1987).
 * "Marching cubes: A high resolution 3D surface construction algorithm".
 * ACM SIGGRAPH Computer Graphics. 21 (4): 163-169
 *
 * Implementation is for 2D (marching squares) and 3D (marching cubes).
 *
 * The input mesh is a Conduit::Node following the Mesh Blueprint
 * convention.  The mesh must be in multi-domain format.
 *
 * Usage example:
 * @verbatim
 *   void foo( conduit::Node &meshNode,
 *             const std::string &topologyName,
 *             const std::string &functionName,
 *             double contourValue )
 *   {
 *     MarchingCubes mc(axom::runtime_policy::Policy::seq,
 *                      axom::getDefaultAllocatorID(),
 *                      MarchingCubesDataParallelism::byPolicy);
 *     mc.setMesh(meshNode, topologyName);
 *     mc.setFunctionField(functionName);
 *     mc.computeIsocontour(contourValue);
 *     axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>
 *       contourMesh(3, min::CellType::Triangle);
 *     mc.populateContourMesh( contourMesh, "cellIdField");
 *   }
 * @endverbatim
 *
 * To avoid confusion between the two meshes, we refer to the input
 * mesh with the scalar function as "parent" and the generated mesh
 * as the "contour".
 *
 * The output contour mesh format can be a mint::UnstructuredMesh or
 * Array data.  IDs of parent cell and domain that generated the
 * individual contour facets are provided.  Blueprint allows users to
 * specify ids for the domains.  If "state/domain_id" exists in the
 * domains, it is used as the domain id.  Otherwise, the domain's
 * iteration index within the multidomain mesh is used.
 *
 * Output arrays use the allocator id specified in the constructor.
 * However, the output mint mesh currently uses host data.  The data
 * output interfaces are interim and subject to change)
 */
class MarchingCubes
{
public:
  using RuntimePolicy = axom::runtime_policy::Policy;
  using DomainIdType = axom::IndexType;
  /*!
   * @brief Constructor sets up runtime preferences for the marching
   * cubes implementation.
   *
   * @param [in] runtimePolicy A value from RuntimePolicy.
   *             The simplest policy is RuntimePolicy::seq, which specifies
   *             running sequentially on the CPU.
   * @param [in] allocatorID Data allocator ID.  Choose something compatible
   *             with \c runtimePolicy.  See \c execution_space.
   * @param [in] dataParallelism Data parallel implementation choice for the legacy backend.
   *             The bump backend accepts but ignores this setting because
   *             bump manages its own internal parallelism.
  */
  MarchingCubes(RuntimePolicy runtimePolicy,
                int allocatorId,
                MarchingCubesDataParallelism dataParallelism);

  /*!
   * @brief Set the input mesh.
   * @param [in] bpMesh Blueprint multi-domain mesh containing scalar field.
   * @param [in] topologyName Name of Blueprint topology to use in \a bpMesh.
   * @param [in] maskField Cell-based std::int32_t mask field.  If provided,
   *             cells where this field evaluates to false are skipped.
   *
   * Array data in \a bpMesh must be accessible in the \a runtimePolicy
   * environment specified in the constructor.  It's an error if not,
   * e.g., using CPU memory with a GPU policy.
   * 
   * Some metadata from \a bpMesh may be cached.  Any change to it
   * after setMesh() leads to undefined behavior.
  */
  void setMesh(const conduit::Node& bpMesh,
               const std::string& topologyName,
               const std::string& maskField = {});

  /*!
   * @brief Set the field containing the nodal function.
   * @param [in] fcnField Name of node-based scalar function values.
  */
  void setFunctionField(const std::string& fcnField);

  /*!
   * @brief Set the mask value.
   * @param [in] maskVal mask value.  If a mask field is given in
   *   setMesh(), compute only for cells whose mask matches this value.
   *
   * The default vask value is 1 unless explicitly set by this method.
   * The mask value has no effect if a mask field is not specified.
  */
  void setMaskValue(int maskVal) { m_maskVal = maskVal; }

  /*!
   * @brief Set how parent-cell ids are numbered for generated contour facets.
   * @param [in] mode A value from MarchingCubesParentCellIdMode.
   *
   * The default is MarchingCubesParentCellIdMode::blueprintZoneId.
   * See that enum for the meaning of each mode and for why the two modes can differ
   * for structured input.  Has no effect unless parent-cell ids are requested
   * (via getContourFacetParents() or the cellIdField of populateContourMesh()).
   *
   * @note The legacyFieldOrder mode only affects structured input;
   *   unstructured input always reports the Blueprint zone id.
  */
  void setParentCellIdMode(MarchingCubesParentCellIdMode mode) { m_parentCellIdMode = mode; }

  /*!
   * @brief Select the bump::extraction::CutField backend 
   *   (vs. the legacy structured-only marching cubes kernel).
   * @param [in] useBump If true, isocontour extraction is delegated to bump,
   *   which additionally supports unstructured single-shape quad (2D) and hex
   *   (3D) meshes.  If false (default), the legacy kernel is used.
   *
   * Only available when Axom is configured with the bump component (AXOM_USE_BUMP).
   * requesting the bump backend otherwise has no effect.
   * The legacy backend supports only structured input.
   *
   * @note The MarchingCubesDataParallelism constructor argument is a legacy
   * backend scan-strategy selector.  The bump backend ignores it and relies on
   * bump's internal parallelism for the selected runtime policy.
   *
   * @note This is transitional: while the bump backend matures it is opt-in so
   * existing users are unaffected.  A future release is expected to make it the
   * default and retire the legacy kernel and its lookup tables.
  */
  void setUseBumpBackend(bool useBump) { m_useBumpBackend = useBump; }

  /*!
   * @brief Select the isosurface robustness policy for the bump backend.
   * @param [in] policy A value from MarchingCubesRobustnessPolicy.
   *
   * See MarchingCubesRobustnessPolicy for the meaning of each value.  The
   * default is MarchingCubesRobustnessPolicy::standard.  Has no effect on the
   * legacy backend, and selecting \c robust currently behaves as \c standard
   * until a robust bump intersector is available.
  */
  void setRobustnessPolicy(MarchingCubesRobustnessPolicy policy) { m_robustnessPolicy = policy; }

  /*!
   * @brief Computes the isocontour.
   * @param [in] contourVal isocontour value
   *
   * Each computeIsocontour call adds to previously computed contour mesh.
  */
  void computeIsocontour(double contourVal = 0.0);

  //!@brief Get number of cells (facets) in the generated contour mesh.
  axom::IndexType getContourCellCount() const { return m_facetCount; }
  //!@brief Get number of cells (facets) in the generated contour mesh.
  axom::IndexType getContourFacetCount() const { return m_facetCount; }

  //!@brief Get number of nodes in the generated contour mesh.
  axom::IndexType getContourNodeCount() const;

  ///@{
  //!@name Access to output contour mesh
  /*!
   * @brief Put generated contour in a mint::UnstructuredMesh.
   * @param mesh Output contour mesh
   * @param cellIdField Name of field to store the array of
   *   parent cells ids, numbered in the row- or column-major
   *   ordering of the nodal scalar function.
   *   If empty, the data is not provided.
   * @param domainIdField Name of field to store the
   *   parent domain ids.  The type of this data is \c DomainIdType.
   *   If omitted, the data is not provided.
   *
   *  If the fields aren't in the mesh, they will be created.
   *
   *  Important: mint::UnstructuredMesh only supports host memory, so
   *  regardless of the allocator ID, this method always deep-copies
   *  data to host memory.  To access the data without deep-copying, see
   *  the other output methods in this name group.
   *
   *  When the bump backend is enabled, its native 3D CutField output may contain polygonal surface elements.
   *  The adaptor triangulates those polygons (reusing bump's welded vertex coordinates).
  */
  void populateContourMesh(axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& mesh,
                           const std::string& cellIdField = {},
                           const std::string& domainIdField = {}) const;

  /*!
   * @brief Copy the richer bump-backed contour mesh into a Blueprint multi-domain mesh.
   * @param [out] bpMesh Output Blueprint multi-domain mesh.
   * @param triangulate If true, convert 3D polygonal surface elements into
   *   triangles while preserving bump's welded coordset.
   *
   * This accessor is available only for contours computed with the bump backend.
   * It preserves bump's native welded representation: line segments in 2D
   * and polygonal surface elements in 3D with Blueprint elements/{connectivity,sizes,offsets}.
   * When \a triangulate is true, 3D polygonal faces are triangulated in the returned Blueprint mesh.
   *
   * Array data in \a bpMesh is copied into the same memory space used by the MarchingCubes object.
   * If the contour was computed with a device policy, callers that need host-readable Blueprint data
   * should copy it to host.
  */
  void populateContourMeshBlueprint(conduit::Node& bpMesh, bool triangulate = false) const;

  /*!
   * @brief Return view of facet corner node indices (connectivity) Array.
   *
   * The array shape is (getContourCellCount(), <spatial dimension>), where
   * the second index is index of the facet corner.
   */
  axom::ArrayView<const axom::IndexType, 2> getContourFacetCorners() const
  {
    return m_facetNodeIds.view();
  }

  /*!
   * @brief Return view of node coordinates Array.
   *
   * The array shape is (getContourNodeCount(), <spatial dimension>), where
   * the second index is the spatial index.
   */
  axom::ArrayView<const double, 2> getContourNodeCoords() const { return m_facetNodeCoords.view(); }

  /*!
   *  @brief Return view of parent cell indices Array.
   *
   *  The buffer size is getContourCellCount().  The parent ID is the
   *  flat index of the cell in the parent domain (see MDMapping),
   *  not counting ghost cells, with row- or major-ordering same as that
   *  for the input scalar function array.
   */
  axom::ArrayView<const axom::IndexType> getContourFacetParents() const
  {
    return m_facetParentIds.view();
  }

  /*!
   *   @brief Return view of parent domain indices Array.
   *   @param allocatorID Allocator id for the output data.  If omitted,
   *          use the id set in the constructor.
   *   The buffer size is getContourCellCount().
   */
  axom::ArrayView<const axom::IndexType> getContourFacetDomainIds() const
  {
    return m_facetDomainIds.view();
  }

  /*!
   *  @brief Give caller posession of the contour data.
   *
   *  This efficiently gives the generated contour data to the caller,
   *  to stay in scope after the MarchingCubes object is deleted.
   *  @param [i] facetNodeIds Node ids for the node at the corners of each facet.
   *  @see getContourFacetCorners().
   *  @param [i] facetNodeCoords Coordinates of each facet node.
   *  @see getContourNodeCoords().
   *  @param [i] facetParentIds Parent cell id of each facet.
   *  @see getContourFacetParents().
   *  @param [i] facetDomainIds Domain id of each facet.
   *  @see getContourFacetDomainIds().
   * 
   *  @pre computeIsocontour() must have been called.
   *  @post outputs can no longer be accessed from object, as though clearOutput() has been called.
   */
  void relinquishContourData(axom::Array<axom::IndexType, 2>& facetNodeIds,
                             axom::Array<double, 2>& facetNodeCoords,
                             axom::Array<axom::IndexType, 1>& facetParentIds,
                             axom::Array<axom::IndexType>& facetDomainIds)
  {
    facetNodeIds.clear();
    facetNodeCoords.clear();
    facetParentIds.clear();
    facetDomainIds.clear();
    m_facetCount = 0;

    facetNodeIds.swap(m_facetNodeIds);
    facetNodeCoords.swap(m_facetNodeCoords);
    facetParentIds.swap(m_facetParentIds);
    facetDomainIds.swap(m_facetDomainIds);
  }

  /*!
   * @brief Give caller possession of the richer bump-backed Blueprint contour.
   * @param [out] bpMesh Output Blueprint multi-domain mesh.
   *
   * This moves the cached bump output nodes without deep-copying them.
   * It is available only for contours computed with the bump backend
   * and leaves this MarchingCubes object with no accessible contour output,
   * as though clearOutput() had been called.
  */
  void relinquishContourDataBlueprint(conduit::Node& bpMesh);
  ///@}

  //! @brief Clear the computed contour mesh.
  void clearOutput();

  // Allow single-domain code to share common scratch space.
  friend detail::marching_cubes::MarchingCubesSingleDomain;

  /*
    TODO: CrossingFlagType can be a boolean value but is wastefully
    stored in 32 bits because of a ROCM scan implementation that adds
    them in the input type without promoting them to our 32-bit output
    type.  When ROCM supports the promotion and RAJA uses it, we can
    change this type to something more efficient.
  */
  using CrossingFlagType = std::uint32_t;

private:
  RuntimePolicy m_runtimePolicy;
  int m_allocatorID = axom::INVALID_ALLOCATOR_ID;

  //! @brief Legacy backend data-parallel scan strategy, or byPolicy.
  MarchingCubesDataParallelism m_dataParallelism = MarchingCubesDataParallelism::byPolicy;

  //! @brief Number of domains.
  axom::IndexType m_domainCount;

  /*!
   * @brief Single-domain implementations.
   *
   * May be longer than m_domainCount (the real count).
  */
  axom::Array<std::shared_ptr<detail::marching_cubes::MarchingCubesSingleDomain>> m_singles;
  std::string m_topologyName;
  std::string m_fcnFieldName;
  std::string m_fcnPath;
  std::string m_maskFieldName;
  std::string m_maskPath;

  int m_maskVal = 1;

  //! @brief How to number parent-cell ids of generated facets.
  MarchingCubesParentCellIdMode m_parentCellIdMode = MarchingCubesParentCellIdMode::blueprintZoneId;

  //! @brief Whether to use the bump CutField backend (opt-in; default legacy).
  bool m_useBumpBackend = false;

  //! @brief Isosurface robustness policy for the bump backend (Phase 6 seam).
  MarchingCubesRobustnessPolicy m_robustnessPolicy = MarchingCubesRobustnessPolicy::standard;

  //! @brief First facet index from each parent domain.
  axom::Array<axom::IndexType> m_facetIndexOffsets;

  //! @brief Facet count over all parent domains.
  axom::IndexType m_facetCount = 0;

  ///@{
  //! @name Scratch space from m_allocatorID, shared among singles
  // Memory alloc is slow on CUDA, so this optimizes space AND time.
  axom::Array<std::uint16_t> m_caseIdsFlat;

  axom::Array<CrossingFlagType> m_crossingFlags;
  axom::Array<axom::IndexType> m_scannedFlags;
  axom::Array<axom::IndexType> m_facetIncrs;
  ///@}

  ///@{
  //!@name Generated contour mesh, shared with singles.

  axom::IndexType m_nodeCount {0};

  /*!
   * @brief Corners (index into m_facetNodeCoords) of generated facets.
   * @see allocateOutputBuffers().
  */
  axom::Array<axom::IndexType, 2> m_facetNodeIds;

  /*!
   * @brief Coordinates of generated surface mesh nodes.
   * @see allocateOutputBuffers().
  */
  axom::Array<double, 2> m_facetNodeCoords;

  axom::Array<axom::IndexType> m_nodeIndexOffsets;

  /*!
   * @brief Flat index of parent cell of facets.
   * @see allocateOutputBuffers().
  */
  axom::Array<IndexType, 1> m_facetParentIds;

  /// @brief Domain ids of facets
  axom::Array<IndexType, 1> m_facetDomainIds;
  ///@}

  //! @brief Allocate output buffers corresponding to runtime policy.
  void allocateOutputBuffers();
};

}  // namespace axom::quest

#endif  // AXOM_USE_CONDUIT
