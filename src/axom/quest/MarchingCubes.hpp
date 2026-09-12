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
 * @brief Selects the legacy backend's scan implementation.
 *
 * \c hybridParallel uses a serial loop but processes less data.
 * \c fullParallel processes more data but has no serial loop.
 * \c byPolicy chooses between them based on the runtime policy.
 *
 * @note This setting controls only the legacy structured-mesh backend.
 *       Bump manages its own parallelism for the selected runtime policy.
 */
enum class MarchingCubesDataParallelism
{
  byPolicy = 0,
  hybridParallel = 1,
  fullParallel = 2
};

/*!
 * @brief Specifies a Bump isosurface robustness policy.
 *
 * Both values currently select \c axom::bump::extraction::FieldIntersector.
 * It classifies a corner as inside when its value is greater than the
 * isovalue, computes edge crossings in single precision, and uses one fixed
 * triangulation per case. Each saddle case therefore has one fixed topology,
 * which may differ from the bilinear or trilinear interpolant. The intersector
 * does not distinguish negative, zero, and positive values or use an asymptotic decider.
 *
 * \c standard selects this implementation. \c robust currently behaves the same as \c standard.
 */
enum class MarchingCubesRobustnessPolicy
{
  standard = 0,
  robust = 1
};

/*!
 * @brief Extracts a contour mesh from a scalar field.
 *
 * The legacy backend implements the original 1987 algorithm:
 * Lorensen, William E.; Cline, Harvey E. (1 August 1987).
 * "Marching cubes: A high resolution 3D surface construction algorithm".
 * ACM SIGGRAPH Computer Graphics. 21 (4): 163-169
 *
 * The MarchingCubes class supports marching squares in 2D and marching cubes in 3D.
 *
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
 * The input mesh is called the parent mesh, and the generated mesh is called the contour mesh.
 *
 * Output is available as arrays or an \c axom::mint::UnstructuredMesh.
 * The arrays identify the parent cell and domain of each facet.
 * The Bump backend can also return its welded Blueprint mesh.
 *
 * If "state/domain_id" exists in the domains, it is used as the domain id.
 * Otherwise, the domain's iteration index within the multidomain mesh is used.
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
   * @param [in] runtimePolicy A value from RuntimePolicy.
   *             The simplest policy is RuntimePolicy::seq, which specifies
   *             running sequentially on the CPU.
   * @param [in] allocatorId Data allocator ID.  Choose something compatible
   *             with \c runtimePolicy.  See \c execution_space.
   * @param [in] dataParallelism Data-parallel implementation for the legacy backend.
   *             The Bump backend accepts but ignores this setting because
   *             Bump manages its own parallelism.
   */
  MarchingCubes(RuntimePolicy runtimePolicy,
                int allocatorId,
                MarchingCubesDataParallelism dataParallelism);

  /*!
   * @brief Set the input mesh.
   * @param [in] bpMesh Blueprint single-domain or multi-domain mesh containing a scalar field.
   * @param [in] topologyName Name of Blueprint topology to use in \a bpMesh.
   * @param [in] maskField Optional cell-based std::int32_t mask field.
   *             Cells whose values differ from the current mask value are skipped.
   *
   * Array data in \a bpMesh must be accessible in the \a runtimePolicy
   * environment specified in the constructor.  It's an error if not,
   * e.g., using CPU memory with a GPU policy.
   * 
   * MarchingCubes retains references to data in \a bpMesh. Do not modify or destroy
   * that data before calling setMesh() again or destroying this object.
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
   * @param [in] maskVal mask value.  If a mask field is given in setMesh(),
   *             compute only for cells whose mask matches this value.
   *
   * The default mask value is 1.
   * The mask value has no effect if a mask field is not specified.
   */
  void setMaskValue(int maskVal) { m_maskVal = maskVal; }

  /*!
   * @brief Enable or disable the \c bump::extraction::CutField backend.
   * @param [in] useBump If true, use Bump. If false, use the legacy backend.
   *
   * The legacy backend accepts structured meshes. Bump also accepts uniform
   * and rectilinear topologies and single-shape unstructured quad or hex meshes.
   * Enabling Bump requires \c AXOM_USE_BUMP.
   *
   * Call this method before setMesh(), which constructs backend-specific
   * workers for the input domains.
   *
   * @note The MarchingCubesDataParallelism constructor argument selects the
   *       legacy backend's scan strategy. Bump ignores it.
   */
  void setUseBumpBackend(bool useBump);

  /*!
   * @brief Select the isosurface robustness policy for the Bump backend.
   * @param [in] policy A value from MarchingCubesRobustnessPolicy.
   *
   * The default is \c MarchingCubesRobustnessPolicy::standard. The legacy
   * backend ignores this setting. Bump currently treats \c robust as \c standard.
   */
  void setRobustnessPolicy(MarchingCubesRobustnessPolicy policy) { m_robustnessPolicy = policy; }

  /*!
   * @brief Computes the isocontour.
   * @param [in] contourVal isocontour value
   *
   * Each call appends to the array output used by populateContourMesh().
   * Call clearOutput() first to replace prior results.
   * Bump's Blueprint output contains only the most recent call for each domain.
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
   * @param cellIdField Name of field to store the flat parent cell ids.
   *        If empty, the data is not provided.
   * @param domainIdField Name of field to store the parent domain ids.
   *        The type of this data is \c DomainIdType.
   *        If omitted, the data is not provided.
   *
   *  If the fields aren't in the mesh, they will be created.
   *
   *  Important: mint::UnstructuredMesh only supports host memory, so
   *  regardless of the allocator ID, this method always deep-copies
   *  data to host memory.  To access the data without deep-copying, see
   *  the other output methods in this name group.
   *
   *  Bump may produce polygonal faces in 3D. This method fan-triangulates
   *  those faces and reuses Bump's welded vertices.
   */
  void populateContourMesh(axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& mesh,
                           const std::string& cellIdField = {},
                           const std::string& domainIdField = {}) const;

  /*!
   * @brief Copy Bump's welded contour into a Blueprint multi-domain mesh.
   * @param [out] bpMesh Output Blueprint multi-domain mesh.
   * @param [in] triangulate If true, convert 3D polygonal surface elements
   *        into triangles while preserving Bump's welded coordset.
   *
   * This method requires the Bump backend. Without triangulation, it preserves
   * Bump's welded segments in 2D and polygonal faces in 3D. The topology stores
   * Blueprint connectivity, sizes, and offsets.
   *
   * Arrays in \a bpMesh use the allocator supplied to the constructor.
   * Callers must copy device data to host before reading it on the host.
   * The mesh contains only the most recent computeIsocontour() call for each input domain.
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
   *  The buffer size is getContourCellCount(). The parent ID is the flat cell
   *  index in the parent domain. For structured meshes, it excludes ghost
   *  cells and follows the scalar field's logical ordering.
   */
  axom::ArrayView<const axom::IndexType> getContourFacetParents() const
  {
    return m_facetParentIds.view();
  }

  /*!
   *   @brief Return view of parent domain indices Array.
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
   *  @post The array accessors return empty views. Cached Bump Blueprint output
   *        remains available until clearOutput() or
   *        relinquishContourDataBlueprint() is called.
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
  }

  /*!
   * @brief Transfer ownership of Bump's welded Blueprint contour to the caller.
   * @param [out] bpMesh Output Blueprint multi-domain mesh.
   *
   * This moves the cached Bump output nodes without deep-copying them.
   * It is available only for contours computed with the Bump backend
   * and leaves this MarchingCubes object with no accessible contour output,
   * as though clearOutput() had been called.
   * Only the most recent computeIsocontour() call for each domain is moved.
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
   * May be longer than m_domainCount (the real count).
   */
  axom::Array<std::shared_ptr<detail::marching_cubes::MarchingCubesSingleDomain>> m_singles;

  /*!
   * @brief Wrapper used when callers pass a single-domain Blueprint mesh.
   *
   * MarchingCubesSingleDomain caches references into the per-domain node, so
   * the synthetic multi-domain parent must outlive setMesh().
   */
  conduit::Node m_singleDomainMesh;

  std::string m_topologyName;
  std::string m_fcnFieldName;
  std::string m_fcnPath;
  std::string m_maskFieldName;
  std::string m_maskPath;

  int m_maskVal {1};

  //! @brief Whether to use the Bump CutField backend.
  bool m_useBumpBackend {false};

  //! @brief Isosurface robustness policy for the Bump backend.
  MarchingCubesRobustnessPolicy m_robustnessPolicy {MarchingCubesRobustnessPolicy::standard};

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
};

}  // namespace axom::quest

#endif  // AXOM_USE_CONDUIT
