// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_CURVILINEARMESH_HPP_
#define MINT_CURVILINEARMESH_HPP_

#include "axom/config.hpp"
#include "axom/mint/mesh/StructuredMesh.hpp"
#include "axom/mint/config.hpp"
#include "axom/export/mint.h"

namespace axom
{
namespace mint
{
// Forward Declarations
class MeshCoordinates;

/*!
 * \class CurvilinearMesh
 *
 * \brief Provides the ability to represent and operate on structured,
 *  curvilinear meshes.
 *
 *  The CurvilinearMesh extends from the StructuredMesh base class and adds
 *  the ability to represent and operate on structured, curvilinear meshes.
 *
 *  A <em> CurvilinearMesh </em> divides the solution domain in cells that are
 *  topologically arranged on a logical, regular grid. Consequently, each node
 *  and cell in the mesh can be uniquely identified by a corresponding logical
 *  i-j-k index. However, the geometry, i.e., the coordinates of the nodes are
 *  explicitly defined and mapped on to the physical domain. For this reason,
 *  curvilinear meshes are also called <em> mapped </em> or <em> body-fitted
 *  </em> meshes.
 *
 *  A <em> CurvilinearMesh </em> object may be constructed using (a) a native
 *  constructor (b) an external constructor or (c) a Sidre constructor, when
 *  Mint is compiled with Sidre support.
 *
 * * <b> Native Constructor </b> <br />
 *
 *    When using native storage, the CurvilinearMesh object owns all associated
 *    memory. Once the object is deleted, all memory is returned to the system.
 *
 * * <b> External Constructor </b> <br />
 *
 *    A CurvilinearMesh object may also be constructed from externally supplied
 *    coordinate buffers. In this case, the calling code owns the memory and
 *    is responsible for properly deallocating the supplied buffers.
 *
 * * <b> Sidre Constructor </b> <br />
 *
 *    When the CurvilinearMesh is associated with a group within a Sidre
 *    hierarchy, Sidre owns all the memory. Once the object goes out-of-scope
 *    all mesh data remains persistent in Sidre.
 *
 * \note When using Sidre, the specified group must conform to the conventions
 *  defined by the <a href="http://llnl-conduit.readthedocs.io/en/latest/">
 *  computational mesh blueprint </a>
 *
 * \see StructuredMesh
 * \see Extent
 */
class CurvilinearMesh : public StructuredMesh
{
public:
  /*!
   * \brief Default constructor. Disabled.
   */
  CurvilinearMesh() = delete;

  /// \name Native Constructors
  /// @{

  /*!
   * \brief Constructs a CurvilinearMesh with the given dimensions.
   *
   * \param [in] Ni number of nodes in the i-direction.
   * \param [in] Nj number of nodes in the j-direction (for dimension >= 2).
   * \param [in] Nk number of nodes in the k-direction (for dimension==3).
   *
   * \pre Ni >= 1
   * \pre Nj >= 1 iff dimension >=2
   * \pre Nk >= 1 iff dimension == 3
   *
   * \post 1 <= getDimension() <= 3
   * \post getCoordinateArray( i ) != nullptr \f$ \forall i \f$
   * \post hasSidreGroup() == false
   * \post isExternal() == false
   */
  CurvilinearMesh(IndexType Ni, IndexType Nj = -1, IndexType Nk = -1);
  /// @}

  /// \name External Constructors
  /// @{

  /*!
   * \brief Constructs a CurvilinearMesh instance given the supplied, external
   *  coordinate buffers.
   *
   * \param [in] Ni number of nodes in the i-direction.
   * \param [in] x pointer to the x-coordinates
   * \param [in] Nj number of nodes in the j-direction (for dimension >= 2).
   * \param [in] y pointer to the y-coordinates (required only for 2D and 3D)
   * \param [in] Nk number of nodes in the k-direction (for dimension==3).
   * \param [in] z pointer to the z-coordinates (required only 3D)
   *
   * \note The calling code maintains ownership of the supplied coordinate
   *  buffers and the responsibility of properly deallocating them.
   *
   * \pre x != nullptr
   * \pre y != nullptr
   * \pre z != nullptr
   *
   * \post 1 <= getDimension() <= 3
   * \post getCoordinateArray( i ) != nullptr \f$ \forall i \f$
   * \post hasSidreGroup() == false.
   * \post isExternal() == true
   */
  CurvilinearMesh(IndexType Ni,
                  double* x,
                  IndexType Nj = -1,
                  double* y = nullptr,
                  IndexType Nk = -1,
                  double* z = nullptr);
  /// @}

#ifdef AXOM_MINT_USE_SIDRE

  /// \name Sidre Constructors
  /// @{

  /*!
   * \brief Create a curvilinear mesh instance from the given Sidre group that
   *  holds mesh data for a structured curvilinear mesh according to the
   *  computational mesh blueprint conventions.
   *
   * \param [in] group pointe to the root group within a Sidre hierarchy.
   * \param [in] topo the name of the topology for this mesh (optional).
   *
   * \note The supplied group is expected to contain a valid curvilinear
   *  structured mesh instance that is conforming to the conventions described
   *  in the conduit <a href="http://llnl-conduit.readthedocs.io/en/latest/">
   *  computational mesh blueprint </a>.
   *
   * \note If a topology name is not provided, the implementation will use the
   *  1st topology group under the parent "topologies" group.
   *
   * \note When using this constructor, all data is owned by Sidre. Once the
   *  mesh object goes out-of-scope, the data will remain persistent in Sidre.
   *
   * \pre group != nullptr
   * \pre blueprint::isValidRootGroup( group )
   *
   * \post hasSidreGroup() == true
   * \post isExternal() == false
   */
  /// @{
  explicit CurvilinearMesh(sidre::Group* group, const std::string& topo = "");
  /// @}

  /*!
   * \brief Create a curvilinear mesh instance, on an empty sidre::Group, with
   *  specified dimensions,  \f$ N_i, N_j, N_k f\$
   *
   * \param [in] group pointer to the Sidre group where to store the mesh
   * \param [in] topo the name of the associated topology (optional)
   * \param [in] coordset the name of the associated coordset group (optional)
   * \param [in] Ni number of nodes in the i-direction.
   * \param [in] Nj number of nodes in the j-direction (for dimension >= 2).
   * \param [in] Nk number of nodes in the k-direction (for dimension==3).
   *
   * \note If a topology and coordset name is not provided, internal defaults
   *  will be used by the implementation.
   *
   * \note When using this constructor, all data is owned by Sidre. Once the
   *  mesh object goes out-of-scope, the data will remain persistent in Sidre.
   *
   * \pre Ni >= 1
   * \pre Nj >= 1 iff dimension >=2
   * \pre Nk >= 1 iff dimension == 3
   *
   * \pre group != nullptr
   * \pre group->getNumViews()==0
   * \pre group->getNumGroups()==0
   *
   * \post hasSidreGroup() == true
   * \post isExternal() == false
   */
  /// @{

  CurvilinearMesh(sidre::Group* group,
                  const std::string& topo,
                  const std::string& coordset,
                  IndexType Ni,
                  IndexType Nj = -1,
                  IndexType Nk = -1);

  AXOM_MINT_EXPORT CurvilinearMesh(sidre::Group* group,
                                   IndexType Ni,
                                   IndexType Nj = -1,
                                   IndexType Nk = -1)
    : CurvilinearMesh(group, "", "", Ni, Nj, Nk)
  { }

    /// @}

/// @}
#endif

  /// \name Virtual methods
  /// @{

  /*!
   * \brief Destructor.
   */
  virtual ~CurvilinearMesh();

  /*!
   * \brief Returns true if the mesh points to external coordinate arrays.
   * \return status true iff external coordinate arrays were supplied.
   */
  virtual bool isExternal() const final override { return m_coordinates->isExternal(); }

  /// \name Nodes
  /// @{

  /*!
   * \brief Copy the coordinates of the given node into the provided buffer.
   *
   * \param [in] nodeID the ID of the node in question.
   * \param [in] coords the buffer to copy the coordinates into, of length at
   *  least getDimension().
   *
   * \pre 0 <= nodeID < getNumberOfNodes()
   * \pre coords != nullptr
   */
  virtual void getNode(IndexType nodeID, double* coords) const final override
  {
    m_coordinates->getCoordinates(nodeID, coords);
  }

  /*!
   * \brief Returns pointer to the requested mesh coordinate buffer.
   *
   * \param [in] dim the dimension of the requested coordinate buffer
   * \return ptr pointer to the coordinate buffer.
   *
   * \note The length of the returned buffer is getNumberOfNodes().
   *
   * \pre dim >= 0 && dim < dimension()
   * \pre dim == X_COORDINATE || dim == Y_COORDINATE || dim == Z_COORDINATE
   */
  /// @{

  virtual double* getCoordinateArray(int dim) final override
  {
    return m_coordinates->getCoordinateArray(dim);
  }

  virtual const double* getCoordinateArray(int dim) const final override
  {
    return m_coordinates->getCoordinateArray(dim);
  }

  /// @}

  /// @}

  /// @}

private:
  /*!
   * \brief Initializes the CurvilinearMesh
   * \note Helper method to be called from a constructor.
   */
  void initialize();

  MeshCoordinates* m_coordinates;

  DISABLE_COPY_AND_ASSIGNMENT(CurvilinearMesh);
  DISABLE_MOVE_AND_ASSIGNMENT(CurvilinearMesh);
};

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_CURVILINEARMESH_HPP_ */
