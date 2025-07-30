// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file MarchingCubesSingleDomain.hpp
 *
 * \brief Consists of classes implementing marching cubes algorithm to
 * compute isocontour from a scalar field in a blueprint mesh.
 */

#ifndef AXOM_QUEST_MARCHINGCUBESSINGLEDOMAIN_H_
#define AXOM_QUEST_MARCHINGCUBESSINGLEDOMAIN_H_

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifdef AXOM_USE_CONDUIT

  // Axom includes
  #include "axom/core/execution/runtime_policy.hpp"
  #include "axom/mint/mesh/UnstructuredMesh.hpp"
  #include "axom/quest/MarchingCubes.hpp"

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
template <int DIM, typename ExecSpace, typename SequentialLoopPolicy>
class MarchingCubesImpl;

/*!
   \@brief Class implementing marching cubes algorithm for a single
   domain.

   This class is an internal detail for multi-domain implementation
   MarchinCubes class, and should not be used outside it.

   \sa MarchingCubes
*/
class MarchingCubesSingleDomain
{
public:
  using RuntimePolicy = axom::runtime_policy::Policy;
  /*!
   \brief Constructor for applying algorithm in a single domain.
  */
  MarchingCubesSingleDomain(MarchingCubes &mc);

  ~MarchingCubesSingleDomain() { }

  /*!
    @brief Intitialize object to a domain.
    \param [in] dom Blueprint single-domain mesh containing scalar field.
    \param [in] topologyName Name of Blueprint topology to use in \a dom
    \param [in] maskField Cell-based std::int32_t mask field.  If provided,
                cells where this field evaluates to false are skipped.

    Array data in \a dom must be accessible in the the \a
    runtimePolicy environment in the constructor.  It's an error if
    not, e.g., using CPU memory with a GPU policy.

    Some data from \a dom may be cached by the constructor.  Any
    change to it without re-initialization leads to undefined
    behavior.

    The mesh coordinates should be stored contiguously.  See
    conduit::blueprint::is_contiguous().  In the future, this
    requirement may be relaxed, possibly at the cost of a
    transformation and storage of the temporary contiguous layout.
  */
  void setDomain(const conduit::Node &dom,
                 const std::string &topologyName,
                 const std::string &maskfield);

  int spatialDimension() const { return m_ndim; }

  /*!
    @brief Specify the field containing the nodal scalar function
    in the input mesh.
    \param [in] fcnField Name of node-based scalar function values.
  */
  void setFunctionField(const std::string &fcnField)
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

  // Methods trivially delegated to implementation.
  void markCrossings() { m_impl->markCrossings(); }
  void scanCrossings() { m_impl->scanCrossings(); }
  void computeFacets() { m_impl->computeFacets(); }

  /*!
    @brief Get the Blueprint domain id specified in \a state/domain_id
    if it is provided, or use the given default if not provided.
  */
  int32_t getDomainId(int32_t defaultId) const;

  //!@brief Get number of cells in the generated contour mesh.
  axom::IndexType getContourCellCount() const { return m_impl->getContourCellCount(); }

  //!@brief Get number of nodes in the generated contour mesh.
  axom::IndexType getContourNodeCount() const { return m_ndim * getContourCellCount(); }

  /*!
    @brief Base class for implementations templated on dimension DIM
    and execution space ExecSpace.

    Implementation details templated on DIM and ExecSpace cannot
    be in MarchingCubesSingleDomain so should live in this class.

    This class allows m_impl to refer to any implementation used
    at runtime.
  */
  struct ImplBase
  {
    /*!
     * @brief Prepare internal data for operating on the given domain.
     *
     * Put in here codes that can't be in MarchingCubesSingleDomain
     * due to template use (DIM and ExecSpace).
     */
    virtual void setDomain(const conduit::Node &dom,
                           const std::string &topologyName,
                           const std::string &maskPath = {}) = 0;

    virtual void setFunctionField(const std::string &fcnFieldName) = 0;
    virtual void setContourValue(double contourVal) = 0;
    virtual void setMaskValue(int maskVal) = 0;

    virtual void setDataParallelism(MarchingCubesDataParallelism dataPar) = 0;

    ///@{
    //! @name Distinct phases in contour generation.

    /*!
     * @brief Compute the contour mesh.
     * @brief Mark parent cells that cross the contour value.
     */
    virtual void markCrossings() = 0;

    //!@brief Scan operations to determine counts and offsets.
    virtual void scanCrossings() = 0;

    //!@brief Compute contour data.
    virtual void computeFacets() = 0;
    ///@}

    ///@{
    //!@name Output methods

    //! @brief Return number of contour mesh facets generated.
    virtual axom::IndexType getContourCellCount() const = 0;
    ///@}

    void setOutputBuffers(axom::ArrayView<axom::IndexType, 2> &facetNodeIds,
                          axom::ArrayView<double, 2> &facetNodeCoords,
                          axom::ArrayView<axom::IndexType, 1> &facetParentIds,
                          axom::IndexType facetIndexOffset)
    {
      m_facetNodeIds = facetNodeIds;
      m_facetNodeCoords = facetNodeCoords;
      m_facetParentIds = facetParentIds;
      m_facetIndexOffset = facetIndexOffset;
    }

    virtual ~ImplBase() { }

    virtual void clearDomain() = 0;

    MarchingCubesDataParallelism m_dataParallelism = MarchingCubesDataParallelism::byPolicy;

    double m_contourVal = 0.0;
    int m_maskVal = 1;
    axom::ArrayView<axom::IndexType, 2> m_facetNodeIds;
    axom::ArrayView<double, 2> m_facetNodeCoords;
    axom::ArrayView<IndexType> m_facetParentIds;
    axom::IndexType m_facetIndexOffset = -1;
  };

  ImplBase &getImpl() { return *m_impl; }

private:
  //! @brief Multi-domain implementation this object is under.
  MarchingCubes &m_mc;

  RuntimePolicy m_runtimePolicy;
  int m_allocatorID = axom::INVALID_ALLOCATOR_ID;

  //! @brief Choice of full or partial data-parallelism, or byPolicy.
  MarchingCubesDataParallelism m_dataParallelism = MarchingCubesDataParallelism::byPolicy;

  //! \brief Computational mesh as a conduit::Node.
  const conduit::Node *m_dom;
  int m_ndim;

  //! @brief Name of Blueprint topology in m_dom.
  std::string m_topologyName;

  std::string m_fcnFieldName;
  //! @brief Path to nodal scalar function in m_dom.
  std::string m_fcnPath;

  std::string m_maskFieldName;
  //! @brief Path to mask in m_dom.
  std::string m_maskPath;

  double m_contourVal = 0.0;
  int m_maskVal = 1;

  std::unique_ptr<ImplBase> m_impl;

  /*!
   * \brief Set the blueprint single-domain mesh.
   *
   * Some data from \a dom may be cached.
   */
  void setDomain(const conduit::Node &dom);

  /// @brief Allocate MarchingCubesImpl object
  std::unique_ptr<ImplBase> newMarchingCubesImpl();

};  // class MarchingCubesSingleDomain

}  // end namespace marching_cubes
}  // end namespace detail
}  // namespace quest
}  // namespace axom

#endif  // AXOM_USE_CONDUIT
#endif  // AXOM_QUEST_MARCHINGCUBES_H_
