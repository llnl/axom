// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/*!
 * \file MarchingCubesSingleDomain.hpp
 *
 * \brief Implements Marching Cubes for one Blueprint domain.
 */

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifndef AXOM_USE_CONDUIT
  #error "MarchingCubesSingleDomain.hpp requires conduit"
#endif

// Axom includes
#include "axom/core/execution/runtime_policy.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/quest/MarchingCubes.hpp"

// Conduit includes
#include "conduit_node.hpp"

// C++ includes
#include <string>

namespace axom::quest::detail::marching_cubes
{
template <int DIM, typename ExecSpace, typename SequentialLoopPolicy>
class MarchingCubesImpl;

/*!
 * @brief Applies Marching Cubes to one Blueprint domain.
 *
 * MarchingCubes uses this internal class for each local domain.
 *
 * \sa MarchingCubes
 */
class MarchingCubesSingleDomain
{
public:
  using RuntimePolicy = axom::runtime_policy::Policy;
  //! \brief Construct a single-domain worker for \a mc.
  MarchingCubesSingleDomain(MarchingCubes& mc);

  ~MarchingCubesSingleDomain() = default;

  /*!
   * @brief Set the Blueprint domain.
   * \param [in] dom Blueprint single-domain mesh containing the scalar field.
   * \param [in] topologyName Name of the Blueprint topology to use in \a dom.
   * \param [in] maskField Optional cell-based std::int32_t mask field.
   *             Cells whose values differ from the current mask value are skipped.
   *
   * Array data in \a dom must be accessible to the runtime policy passed to
   * the MarchingCubes constructor.
   *
   * This object retains references to data in \a dom. Do not modify or destroy
   * that data until setDomain() is called again or this object is destroyed.
   *
   * The legacy backend requires non-interleaved coordinates.
   * The Bump backend accepts any layout supported by its coordset views.
   */
  void setDomain(const conduit::Node& dom,
                 const std::string& topologyName,
                 const std::string& maskField);

  int spatialDimension() const { return m_ndim; }

  /*!
   * @brief Select the nodal scalar field to contour.
   * @param [in] fcnField Name of the vertex-associated scalar field.
  */
  void setFunctionField(const std::string& fcnField)
  {
    m_fcnFieldName = fcnField;
    m_fcnPath = "fields/" + fcnField;
    SLIC_ASSERT(m_dom->has_path(m_fcnPath));
    SLIC_ASSERT(m_dom->fetch_existing(m_fcnPath + "/association").as_string() == "vertex");
    SLIC_ASSERT(m_dom->has_path(m_fcnPath + "/values"));
    m_impl->setFunctionField(fcnField);
  }

  void setContourValue(double contourVal)
  {
    m_contourVal = contourVal;
    if(m_impl)
    {
      m_impl->setContourValue(m_contourVal);
    }
  }

  void setMaskValue(double maskVal)
  {
    m_maskVal = maskVal;
    if(m_impl)
    {
      m_impl->setMaskValue(m_maskVal);
    }
  }

  void setRobustnessPolicy(MarchingCubesRobustnessPolicy policy)
  {
    m_robustnessPolicy = policy;
    if(m_impl)
    {
      m_impl->setRobustnessPolicy(m_robustnessPolicy);
    }
  }

  // Methods trivially delegated to implementation.
  void markCrossings() { m_impl->markCrossings(); }
  void scanCrossings() { m_impl->scanCrossings(); }
  void computeFacets() { m_impl->computeFacets(); }

  //! @brief Return \c state/domain_id, or \a defaultId when the domain omits it.
  int32_t getDomainId(int32_t defaultId) const;

  //! @brief Return the number of cells in the generated contour mesh.
  axom::IndexType getContourCellCount() const { return m_impl->getContourCellCount(); }

  //! @brief Return the number of nodes in the generated contour mesh.
  axom::IndexType getContourNodeCount() const { return m_impl->getContourNodeCount(); }

  /*!
   * @brief Runtime interface for implementations templated on dimension and execution space.
   *
   * This interface lets \c m_impl hold the implementation chosen at runtime.
   */
  struct ImplBase
  {
    /*!
     * @brief Prepare internal data for operating on the given domain.
     *
     * Implementations use the compile-time dimension and execution space.
     */
    virtual void setDomain(const conduit::Node& dom,
                           const std::string& topologyName,
                           const std::string& maskPath = {}) = 0;

    virtual void setFunctionField(const std::string& fcnFieldName) = 0;
    virtual void setContourValue(double contourVal) = 0;
    virtual void setMaskValue(int maskVal) = 0;

    /*!
     * @brief Set the Bump isosurface robustness policy.
     *
     * The legacy implementation keeps this no-op default.
     */
    virtual void setRobustnessPolicy(MarchingCubesRobustnessPolicy) { }

    virtual void setDataParallelism(MarchingCubesDataParallelism dataPar) = 0;

    ///@{
    //! @name Distinct phases in contour generation.

    /*!
     * @brief Mark parent cells that cross the contour value.
     */
    virtual void markCrossings() = 0;

    //! @brief Determine output counts and offsets.
    virtual void scanCrossings() = 0;

    //! @brief Generate the contour data.
    virtual void computeFacets() = 0;
    ///@}

    ///@{
    //!@name Output methods

    //! @brief Return the number of generated contour facets.
    virtual axom::IndexType getContourCellCount() const = 0;

    //! @brief Return the number of generated contour nodes.
    virtual axom::IndexType getContourNodeCount() const = 0;

    //! @brief Whether this implementation has a Blueprint contour.
    virtual bool hasContourMeshBlueprint() const { return false; }

    /*!
     * @brief Copy the implementation's Blueprint contour.
     *
     * The legacy backend does not provide this representation; callers should
     * check hasContourMeshBlueprint() before invoking this method.
     */
    virtual void copyContourMeshBlueprint(conduit::Node& bpMesh, bool triangulate) const
    {
      AXOM_UNUSED_VAR(triangulate);
      bpMesh.reset();
    }

    /*!
     * @brief Move the implementation's Blueprint contour.
     *
     * The legacy backend does not provide this representation; callers should
     * check hasContourMeshBlueprint() before invoking this method.
     */
    virtual void relinquishContourMeshBlueprint(conduit::Node& bpMesh) { bpMesh.reset(); }
    ///@}

    void setOutputBuffers(axom::ArrayView<axom::IndexType, 2>& facetNodeIds,
                          axom::ArrayView<double, 2>& facetNodeCoords,
                          axom::ArrayView<axom::IndexType, 1>& facetParentIds,
                          axom::IndexType facetIndexOffset,
                          axom::IndexType nodeIndexOffset)
    {
      m_facetNodeIds = facetNodeIds;
      m_facetNodeCoords = facetNodeCoords;
      m_facetParentIds = facetParentIds;
      m_facetIndexOffset = facetIndexOffset;
      m_nodeIndexOffset = nodeIndexOffset;
    }

    virtual ~ImplBase() = default;

    virtual void clearDomain() = 0;

    MarchingCubesDataParallelism m_dataParallelism = MarchingCubesDataParallelism::byPolicy;

    double m_contourVal = 0.0;
    int m_maskVal = 1;
    axom::ArrayView<axom::IndexType, 2> m_facetNodeIds;
    axom::ArrayView<double, 2> m_facetNodeCoords;
    axom::ArrayView<IndexType> m_facetParentIds;
    axom::IndexType m_facetIndexOffset = -1;
    axom::IndexType m_nodeIndexOffset = -1;
  };

  ImplBase& getImpl() { return *m_impl; }
  const ImplBase& getImpl() const { return *m_impl; }

private:
  /*!
   * \brief Cache a Blueprint single-domain mesh.
   *
   * The implementation retains references to data in \a dom.
   */
  void setDomain(const conduit::Node& dom);

  //! @brief Create the backend implementation selected at runtime.
  std::unique_ptr<ImplBase> newMarchingCubesImpl();

private:
  //! @brief Owning multi-domain MarchingCubes object.
  MarchingCubes& m_mc;

  RuntimePolicy m_runtimePolicy;
  int m_allocatorID {axom::INVALID_ALLOCATOR_ID};

  //! @brief Choice of full or partial data-parallelism, or byPolicy.
  MarchingCubesDataParallelism m_dataParallelism {MarchingCubesDataParallelism::byPolicy};

  //! \brief Nonowning pointer to the input Blueprint domain.
  const conduit::Node* m_dom;
  int m_ndim;

  //! @brief Name of Blueprint topology in m_dom.
  std::string m_topologyName;

  std::string m_fcnFieldName;
  //! @brief Path to nodal scalar function in m_dom.
  std::string m_fcnPath;

  std::string m_maskFieldName;
  //! @brief Path to mask in m_dom.
  std::string m_maskPath;

  double m_contourVal {0.0};
  int m_maskVal {1};
  MarchingCubesRobustnessPolicy m_robustnessPolicy {MarchingCubesRobustnessPolicy::standard};

  std::unique_ptr<ImplBase> m_impl;
};

}  // namespace axom::quest::detail::marching_cubes
