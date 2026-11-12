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
 *
 * @note The legacy backend is the default. See MarchingCubesBackend for supported topologies.
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
 * @brief Selects the legacy backend's scan implementation.
 *
 * \c hybridParallel uses a serial loop but processes less data.
 * \c fullParallel processes more data but has no serial loop.
 * \c byPolicy chooses between them based on the runtime policy.
 *
 * @note Bump ignores this setting and manages parallelism through the runtime policy.
 */
enum class MarchingCubesDataParallelism
{
  byPolicy = 0,
  hybridParallel = 1,
  fullParallel = 2
};

/*!
 * @brief Selects the MarchingCubes implementation.
 *
 * \c legacy supports structured meshes, including strided layouts.
 * \c bump uses Bump's CutField. It supports those structured layouts, uniform and
 * rectilinear topologies, and single-shape unstructured quads in 2D or hexes in 3D.
 * Only \c bump provides welded Blueprint contour output.
 */
enum class MarchingCubesBackend
{
  legacy = 0,
  bump = 1
};

/*!
 * @brief Extracts a contour mesh from a scalar field.
 *
 * The legacy backend implements the original 1987 algorithm:
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
 *     axom::quest::MarchingCubes mc(
 *       axom::runtime_policy::Policy::seq,
 *       axom::getDefaultAllocatorID(),
 *       axom::quest::MarchingCubesDataParallelism::byPolicy);
 *     mc.setMesh(meshNode, topologyName);
 *     mc.setFunctionField(functionName);
 *     mc.computeIsocontour(contourValue);
 *     axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>
 *       contourMesh(3, axom::mint::CellType::TRIANGLE);
 *     mc.populateContourMesh(contourMesh, "cellIdField");
 *   }
 * @endverbatim
 *
 * Output is available as arrays or an \c axom::mint::UnstructuredMesh.
 * The arrays identify the parent cell and domain of each facet.
 * The Bump backend can also return its welded Blueprint mesh.
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
   * @brief Configure the execution policy, allocator, and legacy scan strategy.
   *
   * @param [in] runtimePolicy Execution policy. RuntimePolicy::seq runs sequentially on the CPU.
   * @param [in] allocatorId Allocator ID compatible with \a runtimePolicy.
   *             See \c execution_space.
   * @param [in] dataParallelism Legacy scan strategy. Bump ignores this setting.
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
   * @brief Select the backend implementation.
   * @param [in] backend Backend to use. Defaults to MarchingCubesBackend::legacy.
   *
   * Call this before setMesh(). If the backend changes, call setMesh() again to
   * rebuild the domain workers. Otherwise extraction and Blueprint output report an error.
   */
  void setBackend(MarchingCubesBackend backend);

  //! @brief Return the selected backend.
  MarchingCubesBackend getBackend() const { return m_backend; }

  /*!
   * @brief Compute the isocontour.
   * @param [in] contourVal Isocontour value.
   *
   * Each call appends to the array output used by populateContourMesh().
   * Call clearOutput() first to replace prior results.
   * The Blueprint output methods hold one contour per domain,
   * so they require exactly one call since the last clearOutput().
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
   * @param [in,out] mesh Mesh to receive the contour nodes, cells and requested fields.
   * @param [in] cellIdField Field name for flat parent cell ids. Empty means no field.
   * @param [in] domainIdField Field name for parent domain ids of type \c DomainIdType.
   *        Empty means no field.
   *
   * This method creates missing fields.
   *
   * mint::UnstructuredMesh requires host memory, so this method copies data to the host.
   * Use the array or Blueprint output methods to avoid that copy.
   *
   * In 3D, this method fan-triangulates Bump's polygons and reuses its welded vertices.
   */
  void populateContourMesh(axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& mesh,
                           const std::string& cellIdField = {},
                           const std::string& domainIdField = {}) const;

  /*!
   * @brief Copy Bump's welded contour into a Blueprint multi-domain mesh.
   * @param [out] bpMesh Output Blueprint multi-domain mesh.
   * @param [in] triangulate If true, fan-triangulate 3D polygons. Coordinates remain welded.
   *
   * Requires the Bump backend. Output contains welded segments in 2D and polygons
   * in 3D unless \a triangulate is true. The topology includes connectivity, sizes and offsets.
   *
   * Arrays use the constructor's allocator. Copy device data to host before reading it there.
   *
   * @pre computeIsocontour() was called exactly once since the last clearOutput().
   *
   * Empty contours produce no output domain. Each output domain has \c state/domain_id.
   * If all contours are empty, \a bpMesh is empty.
   */
  void populateContourMeshBlueprint(conduit::Node& bpMesh, bool triangulate = false) const;

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
   *  @post The array accessors return empty views. Cached Bump Blueprint output remains
   *        available until clearOutput() or relinquishContourDataBlueprint() is called.
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
   * @brief Transfer ownership of Bump's welded Blueprint contour to the caller.
   * @param [out] bpMesh Output Blueprint multi-domain mesh.
   *
   * Requires the Bump backend. Moves cached output without copying its arrays,
   * then clears this object's contour output as clearOutput() would.
   * The precondition and empty-domain rules match populateContourMeshBlueprint().
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

  //! @brief Report an error if the domain workers do not use the selected backend.
  void checkBackendMatchesMesh() const;

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

  //! @brief Number of computeIsocontour() calls since the last clearOutput().
  axom::IndexType m_contourCount {0};

  //! @brief Selected backend implementation.
  MarchingCubesBackend m_backend {MarchingCubesBackend::legacy};

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
