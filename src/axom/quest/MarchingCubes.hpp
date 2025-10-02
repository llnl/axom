// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * @file MarchingCubes.hpp
 *
 * @brief Consists of classes implementing marching cubes algorithm to
 * compute isocontour from a scalar field in a blueprint mesh.
 */

#ifndef AXOM_QUEST_MARCHINGCUBES_H_
#define AXOM_QUEST_MARCHINGCUBES_H_

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

namespace axom
{
namespace quest
{
namespace detail
{
namespace marching_cubes
{
class MarchingCubesSingleDomain;
}  // namespace marching_cubes
}  // namespace detail

/*!
 * @brief Enum for implementation.
 *
 * Partial parallel implementation uses a non-parallizable loop and
 * processes less data.  It has been shown to work well on CPUs.
 * Full parallel implementation processes more data, but parallelizes
 * fully and has been shown to work well on GPUs.  byPolicy chooses
 * based on runtime policy.
 */
enum class MarchingCubesDataParallelism
{
  byPolicy = 0,
  hybridParallel = 1,
  fullParallel = 2
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
   * @param [in] dataParallelism Data parallel implementation choice.
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
  void setMesh(const conduit::Node &bpMesh,
               const std::string &topologyName,
               const std::string &maskField = {});

  /*!
   * @brief Set the field containing the nodal function.
   * @param [in] fcnField Name of node-based scalar function values.
  */
  void setFunctionField(const std::string &fcnField);

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
  */
  void populateContourMesh(axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> &mesh,
                           const std::string &cellIdField = {},
                           const std::string &domainIdField = {}) const;

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
  void relinquishContourData(axom::Array<axom::IndexType, 2> &facetNodeIds,
                             axom::Array<double, 2> &facetNodeCoords,
                             axom::Array<axom::IndexType, 1> &facetParentIds,
                             axom::Array<axom::IndexType> &facetDomainIds)
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

  //! @brief Choice of full or partial data-parallelism, or byPolicy.
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

}  // namespace quest
}  // namespace axom

#endif  // AXOM_USE_CONDUIT
#endif  // AXOM_QUEST_MARCHINGCUBES_H_
