// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/*!
 * @file MarchingCubes.hpp
 *
 * @brief Extracts 2D and 3D isocontours from scalar fields on Blueprint meshes.
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
 * @brief Selects the scan implementation for the legacy marching cubes data-parallel implementation.
 *
 * \c hybridParallel uses a serial loop but processes less data.
 * \c fullParallel processes more data but has no serial loop.
 * \c byPolicy chooses between them based on the runtime policy.
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
 * otherwise they differ by a stride-order permutation.
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
 * @brief Extracts a contour mesh from a scalar field.
 *
 * This implementation is for the original 1987 algorithm:
 * Lorensen, William E.; Cline, Harvey E. (1 August 1987).
 * "Marching cubes: A high resolution 3D surface construction algorithm".
 * ACM SIGGRAPH Computer Graphics. 21 (4): 163-169
 *
 * This class supports 2D ("marching squares") and 3D ("marching cubes") geometries.
 * The input is a single-domain or multi-domain Conduit Blueprint mesh.
 *
 * Usage example:
 * @verbatim
 *   void foo( conduit::Node &meshNode,
 *             const std::string &topologyName,
 *             const std::string &functionName,
 *             double contourValue )
 *   {
 *     axom::quest::MarchingCubes mc(axom::runtime_policy::Policy::seq,
 *                                   axom::getDefaultAllocatorID(),
 *                                   axom::quest::MarchingCubesDataParallelism::byPolicy);
 *     mc.setMesh(meshNode, topologyName);
 *     mc.setFunctionField(functionName);
 *     mc.computeIsocontour(contourValue);
 *     axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>
 *       contourMesh(3, axom::mint::CellType::TRIANGLE);
 *     mc.populateContourMesh(contourMesh, "cellIdField");
 *   }
 * @endverbatim
 *
 * The input is the parent mesh, and the generated output is the contour mesh.
 *
 * Output is available as arrays or an \c axom::mint::UnstructuredMesh.
 * The arrays identify the parent cell and domain of each facet.
 *
 * If a domain contains \c state/domain_id, that value becomes its domain id.
 * Otherwise, MarchingCubes uses the domain's iteration index.
 *
 * Output arrays use the allocator specified in the constructor.
 * The Mint output always uses host memory.
 */
class MarchingCubes
{
public:
  using RuntimePolicy = axom::runtime_policy::Policy;
  using DomainIdType = axom::IndexType;
  /*!
   * @brief Configure the execution policy, allocator, and scan strategy.
   *
   * @param [in] runtimePolicy A value from RuntimePolicy.
   *             The simplest policy is RuntimePolicy::seq, which specifies
   *             running sequentially on the CPU.
   * @param [in] allocatorId Data allocator ID. Choose one compatible with \c runtimePolicy.
   *             See \c execution_space.
   * @param [in] dataParallelism Data-parallel implementation choice for the legacy backend.
   *             The bump backend accepts but ignores this setting because
   *             bump manages its own internal parallelism.
  */
  MarchingCubes(RuntimePolicy runtimePolicy,
                int allocatorId,
                MarchingCubesDataParallelism dataParallelism);

  /*!
   * @brief Set the input mesh.
   * @param [in] bpMesh Blueprint single-domain or multi-domain mesh containing
   *             the topology and fields to use.
   * @param [in] topologyName Name of Blueprint topology to use in \a bpMesh.
   * @param [in] maskField Optional cell-based std::int32_t mask field.
   *             Cells whose values differ from the current mask value are skipped.
   *
   * Array data in \a bpMesh must be accessible to the \a runtimePolicy passed
   * to the constructor. For example, a GPU policy cannot use host-only memory.
   *
   * MarchingCubes retains references to data in \a bpMesh. Do not modify or destroy
   * that data before calling setMesh() again or destroying this object.
   */
  void setMesh(const conduit::Node& bpMesh,
               const std::string& topologyName,
               const std::string& maskField = {});

  /*!
   * @brief Select the nodal scalar field to contour.
   * @param [in] fcnField Name of the vertex-associated scalar field.
   */
  void setFunctionField(const std::string& fcnField);

  /*!
   * @brief Set the mask value.
   * @param [in] maskVal Mask value. If setMesh() received a mask field,
   *             compute only for cells whose mask matches this value.
   *
   * The default mask value is 1.
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
   * @brief Compute the isocontour.
   * @param [in] contourVal Isocontour value.
   *
   * Each call appends to the array output used by populateContourMesh().
   * Call clearOutput() first to replace prior results.
   */
  void computeIsocontour(double contourVal = 0.0);

  //! @brief Get the number of cells (facets) in the generated contour mesh.
  axom::IndexType getContourCellCount() const { return m_facetCount; }
  //! @brief Get the number of cells (facets) in the generated contour mesh.
  axom::IndexType getContourFacetCount() const { return m_facetCount; }

  //! @brief Get the number of nodes in the generated contour mesh.
  axom::IndexType getContourNodeCount() const;

  ///@{
  //!@name Access to output contour mesh
  /*!
   * @brief Copy the generated contour into a mint::UnstructuredMesh.
   * @param mesh Output contour mesh.
   * @param cellIdField Name of field to store the flat parent cell ids.
   *        If empty, the data is not provided.
   * @param domainIdField Name of field to store the parent domain ids.
   *        The type of this data is \c DomainIdType.
   *        If omitted, the data is not provided.
   *
   *  The method creates the requested fields when they do not exist.
   *
   *  mint::UnstructuredMesh supports only host memory, so this method always
   *  deep-copies data to the host. Use the array output methods to avoid that copy.
   *
   *  When the bump backend is enabled, its native 3D CutField output may contain polygonal
   *  surface elements.  The adaptor fan-triangulates those polygons when filling the legacy
   *  fixed-stride output consumed here, so this method still populates a triangle mesh in 3D.
   */
  void populateContourMesh(axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& mesh,
                           const std::string& cellIdField = {},
                           const std::string& domainIdField = {}) const;

  /*!
   * @brief Copy the richer bump-backed contour mesh into a Blueprint multi-domain mesh.
   * @param [out] bpMesh Output Blueprint multi-domain mesh.
   *
   * This accessor is available only for contours computed with the bump backend.
   * It preserves bump's native welded representation: line segments in 2D
   * and polygonal surface elements in 3D with Blueprint elements/{connectivity,sizes,offsets}.
   * The existing fixed-stride array accessors and populateContourMesh() still expose the
   * legacy un-welded triangle/segment soup; 3D polygonal faces are fan-triangulated there.
   *
   * Array data in \a bpMesh is copied into the same memory space used by the MarchingCubes object.
   * If the contour was computed with a device policy, callers that need host-readable Blueprint data
   * should copy it to host.
  */
  void populateContourMeshBlueprint(conduit::Node& bpMesh) const;

  /*!
   * @brief Return a view of the facet connectivity array.
   *
   * The array shape is (getContourCellCount(), <spatial dimension>),
   * where the second index identifies the facet corner.
   */
  axom::ArrayView<const axom::IndexType, 2> getContourFacetCorners() const
  {
    return m_facetNodeIds.view();
  }

  /*!
   * @brief Return a view of the node-coordinate array.
   *
   * The array shape is (getContourNodeCount(), <spatial dimension>),
   * where the second index is the coordinate axis.
   */
  axom::ArrayView<const double, 2> getContourNodeCoords() const { return m_facetNodeCoords.view(); }

  /*!
   *  @brief Return a view of the parent-cell index array.
   *
   *  The buffer size is getContourCellCount(). The parent ID is the flat cell
   *  index in the parent domain. For structured meshes, it excludes ghost cells
   *  and follows the scalar field's logical ordering.
   */
  axom::ArrayView<const axom::IndexType> getContourFacetParents() const
  {
    return m_facetParentIds.view();
  }

  /*!
   *  @brief Return a view of the parent-domain index array.
   *   The buffer size is getContourCellCount().
   */
  axom::ArrayView<const axom::IndexType> getContourFacetDomainIds() const
  {
    return m_facetDomainIds.view();
  }

  /*!
   *  @brief Transfer ownership of the contour arrays to the caller.
   *
   *  @param [out] facetNodeIds Node ids for the nodes at each facet corner.
   *  @see getContourFacetCorners().
   *  @param [out] facetNodeCoords Coordinates of each facet node.
   *  @see getContourNodeCoords().
   *  @param [out] facetParentIds Parent cell id of each facet.
   *  @see getContourFacetParents().
   *  @param [out] facetDomainIds Domain id of each facet.
   *  @see getContourFacetDomainIds().
   *
   *  @pre computeIsocontour() must have been called.
   *  @post The array accessors return empty views, as though clearOutput() had been called.
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
    m_nodeCount = 0;

    facetNodeIds.swap(m_facetNodeIds);
    facetNodeCoords.swap(m_facetNodeCoords);
    facetParentIds.swap(m_facetParentIds);
    facetDomainIds.swap(m_facetDomainIds);

    // The swaps left this object holding the caller's arrays, whose allocator
    // may differ from m_allocatorID. Recreate empty outputs in this object's
    // memory space so later computeIsocontour() kernels can write to them.
    const axom::StackArray<axom::IndexType, 2> twoZeros {0, 0};
    m_facetNodeIds = axom::Array<axom::IndexType, 2>(twoZeros, m_allocatorID);
    m_facetNodeCoords = axom::Array<double, 2>(twoZeros, m_allocatorID);
    m_facetParentIds = axom::Array<axom::IndexType>(0, 0, m_allocatorID);
    m_facetDomainIds = axom::Array<axom::IndexType>(0, 0, m_allocatorID);
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
  //! @brief Allocate output buffers corresponding to runtime policy.
  void allocateOutputBuffers();

private:
  RuntimePolicy m_runtimePolicy;
  int m_allocatorID {axom::INVALID_ALLOCATOR_ID};

  //! @brief Legacy backend data-parallel scan strategy, or byPolicy.
  MarchingCubesDataParallelism m_dataParallelism {MarchingCubesDataParallelism::byPolicy};

  //! @brief Number of domains.
  axom::IndexType m_domainCount {0};

  /*!
   * @brief Single-domain implementations.
   *
   * Workers are reused across setMesh() calls, so this array can be longer than the current domain count
   */
  axom::Array<std::shared_ptr<detail::marching_cubes::MarchingCubesSingleDomain>> m_singles;

  /*!
   * @brief Wrapper used when callers pass a single-domain Blueprint mesh.
   *
   * Single-domain workers cache references to this wrapper's child node,
   * so the wrapper must remain alive while the workers use it.
   */
  conduit::Node m_singleDomainMesh;

  std::string m_topologyName;
  std::string m_fcnFieldName;
  std::string m_fcnPath;
  std::string m_maskFieldName;
  std::string m_maskPath;

  int m_maskVal {1};

  //! @brief How to number parent-cell ids of generated facets.
  MarchingCubesParentCellIdMode m_parentCellIdMode {MarchingCubesParentCellIdMode::blueprintZoneId};

  //! @brief Whether to use the bump CutField backend (opt-in; default legacy).
  bool m_useBumpBackend {false};

  //! @brief Isosurface robustness policy for the bump backend
  MarchingCubesRobustnessPolicy m_robustnessPolicy {MarchingCubesRobustnessPolicy::standard};

  //! @brief First facet index from each parent domain.
  axom::Array<axom::IndexType> m_facetIndexOffsets;

  //! @brief Facet count over all parent domains.
  axom::IndexType m_facetCount = 0;

  ///@{
  //! @name Scratch arrays allocated with m_allocatorID and shared among workers
  // We reuse these arrays because device allocations are expensive.
  axom::Array<std::uint16_t> m_caseIdsFlat;

  axom::Array<CrossingFlagType> m_crossingFlags;
  axom::Array<axom::IndexType> m_scannedFlags;
  axom::Array<axom::IndexType> m_facetIncrs;
  ///@}

  ///@{
  //! @name Generated contour mesh, shared with single-domain workers.

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

  //! @brief First node index from each parent domain.
  axom::Array<axom::IndexType> m_nodeIndexOffsets;

  /*!
   * @brief Flat index of parent cell of facets.
   * @see allocateOutputBuffers().
  */
  axom::Array<IndexType, 1> m_facetParentIds;

  /// @brief Domain ids of facets.
  axom::Array<IndexType, 1> m_facetDomainIds;
  ///@}
};

}  // namespace axom::quest

#endif  // AXOM_USE_CONDUIT
