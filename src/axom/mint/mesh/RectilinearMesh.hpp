// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_RECTILINEARMESH_HPP_
#define MINT_RECTILINEARMESH_HPP_

#include "axom/config.hpp"
#include "axom/mint/mesh/StructuredMesh.hpp"
#include "axom/mint/config.hpp"
#include "axom/export/mint.h"

namespace axom
{
namespace mint
{
/*!
 * \class RectilinearMesh
 *
 * \brief Provides the ability to represent and operate on a structured
 *  rectilinear mesh.
 *
 * The RectilinearMesh object extends the StructuredMesh base class and adds
 * the ability to represent and operate on structured, rectilinear meshes.
 *
 * A <em> RectilinearMesh </em> divides the solution domain in nodes and cells
 * that are arranged on a regular lattice that is axis-aligned with the
 * Cartesian coordinate axes. Consequently, the cells are always rectangular.
 * Due to the regular topology, each node/cell on the mesh can be uniquely
 * identified by a corresponding logical i-j-k index. However, the nodal
 * coordinates along each axis are explicitly defined and may have variable
 * spacing in between them. Given a logical i-j-k index of a node, the physical
 * coordinates of the corresponding node can be evaluated by taking the
 * Cartesian product of the corresponding coordinate along each direction. For
 * that reason, <em> rectilinear </em> meshes are also often called product
 * meshes.
 *
 * A <em> RectilinearMesh </em> object may be constructed using (a) a native
 * constructor, (b) an external constructor, or (c) Sidre, when Mint is
 * compiled with Sidre support.
 *
 * * <b> Native Constructor </b>
 *
 *   When using native storage, the RectilinearMesh object owns all associated
 *   memory. Once the object is deleted, all memory is returned to the system.
 *
 * * <b> External Constructor </b>
 *
 *   A RectilinearMesh may also be constructed from externally supplied
 *   coordinate buffers. In this case, the callind code owns the memory and is
 *   responsible for properly deallocating the supplied buffers.
 *
 * * <b> Sidre Constructor </b>
 *
 *   A RectilinearMesh may also be constructed and be associated with a
 *   sidre::Group. In this case, all memory is managed by Sidre. Once the
 *   object goes out-of-scope the data remains persistent in Sidre. Last,
 *   the data is stored in Sidre according to the conventions described in the
 *   <a href="http://llnl-conduit.readthedocs.io/en/latest/"> computational
 *   mesh blueprint </a>.
 *
 * \see StructuredMesh
 * \see Extent
 */
class RectilinearMesh : public StructuredMesh
{
public:
  /*!
   * \brief Default constructor. Disabled.
   */
  RectilinearMesh() = delete;

  /// \name Native Constructors
  /// @{

  /*!
   * \brief Constructs a rectilinear mesh of specified dimensions.
   *
   * \param [in] Ni number of nodes in the i-direction
   * \param [in] Nj number of nodes in the j-direction (for dimension >= 2)
   * \param [in] Nk number of nodes in the k-direction (for dimension == 3)
   *
   * \pre Ni >= 1
   * \pre Nj >= 1 iff dimension >= 2
   * \pre Nk >= 1 iff dimension == 3
   *
   * \post 1 <= dimension <= 3
   * \post getCoordinateArray( i ) != nullptr \f$ \forall i \f$
   * \post hasSidreGroup() == false
   * \post isExternal() == false
   */
  RectilinearMesh(IndexType Ni, IndexType Nj = -1, IndexType Nk = -1);

  /// @}

  /// \name External Constructors
  /// @{

  /*!
   * \brief Constructs a rectilinear mesh with the specified logical extent and
   *  external coordinate buffers.
   *
   * \param [in] Ni number of nodes in the i-direction.
   * \param [in] x coordinates along the i-direction.
   * \param [in] Nj number of nodes in the j-direction (for dimension >= 2).
   * \param [in] y coordinates along the j-direction (required for 2D/3D).
   * \param [in] Nk number of nodes in the k-direction (for dimension == 3)
   * \param [in] z coordinates along the k-direction (required for 3D only)
   *
   * \note The calling code maintains ownership of the supplied coordinate
   *  buffers and the responsibility of properly deallocating them.
   *
   * \pre x != nullptr
   * \pre y != nullptr iff dimension >= 2
   * \pre z != nullptr iff dimension == 3
   *
   * \post 1 <= getDimension() <= 3
   * \post getCoordinateArray( i ) != nullptr \f$ \forall i \f$
   * \post hasSidreGroup() == false.
   * \post isExternal() == true
   */
  RectilinearMesh(IndexType Ni,
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
   * \brief Creates a RectilinearMesh from a given Sidre group that holds
   *  mesh data for a rectilinear mesh, according to the conventions described
   *  by the computational mesh blueprint.
   *
   * \param [in] group pointer to the Sidre group associated with the mesh
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
  explicit RectilinearMesh(sidre::Group* group, const std::string& topo = "");

  /*!
   * \brief Creates a RectilinearMesh, on an empty Sidre group, given the
   *  desired mesh dimensions, \f$ N_i, N_j, N_k f\$
   *
   * \param [in] group pointer to the Sidre group where to store the mesh.
   * \param [in] topo the name of the associated topology (optional)
   * \param [in] coordset the name of the associated coordset (optional)
   * \param [in] Ni number of nodes in the i-direction
   * \param [in] Nj number of nodes in the j-direction (for 2D or 3D)
   * \param [in] Nk number of nodes in the k-direction (for 3D only)
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
   *
   */
  /// @{

  RectilinearMesh(sidre::Group* group,
                  const std::string& topo,
                  const std::string& coordset,
                  IndexType Ni,
                  IndexType Nj = -1,
                  IndexType Nk = -1);

  AXOM_MINT_EXPORT RectilinearMesh(sidre::Group* group,
                                   IndexType Ni,
                                   IndexType Nj = -1,
                                   IndexType Nk = -1)
    : RectilinearMesh(group, "", "", Ni, Nj, Nk)
  { }

    /// @}

/// @}
#endif

  /// \name Virtual methods
  /// @{

  /*!
   * \brief Destructor.
   */
  virtual ~RectilinearMesh();

  /*!
   * \brief Returns true if coordinates of the mesh are externally supplied.
   * \return status true if the coordinates point to external buffers.
   */
  virtual bool isExternal() const final override { return m_coordinates[0]->isExternal(); }

  /*!
   * \brief Copy the coordinates of the given node into the provided buffer.
   *
   * \param [in] nodeID the ID of the node in question.
   * \param [in] coords the buffer to copy the coordinates into, of length at
   *  least getDimension().
   *
   * \note This method is provided primarily for convenience. It is a virtual
   *  method and should be used with caution. We advise to not call this method
   *  inside a tight loop to avoid the performance overhead of virtual calls.
   *  Instead, opt to use getCoordinateArray() method which can give direct
   *  memory access to the coordinates of a node.
   *
   * \pre 0 <= nodeID < getNumberOfNodes()
   * \pre coords != nullptr
   */
  virtual void getNode(IndexType nodeID, double* node) const final override;

  /*!
   * \brief Returns pointer to the nodal positions in the specified dimension.
   *
   * \param[in] dim the specified dimension
   *
   * \note the returned buffer is of length getNodeResolution( dim ).
   *
   * \pre dim >= 0 && dim < dimension()
   * \pre dim == X_COORDINATE || dim == Y_COORDINATE || dim == Z_COORDINATE
   */
  /// @{

  virtual double* getCoordinateArray(int dim) final override;
  virtual const double* getCoordinateArray(int dim) const final override;

  /// @}

  /// @}

private:
  /*!
   * \brief Initializes the RectilinearMesh
   * \note Helper method to be called from a constructor.
   */
  void initialize();

  /*!
   * \brief Allocates buffers for the coordinates along each coordinate axis.
   * \note Helper method to be called from constructors.
   */
  void allocateCoords();

#ifdef AXOM_MINT_USE_SIDRE

  /*!
   * \brief Allocates buffers for the coorindates on Sidre.
   * \note Helper method to be called from Sidre constructors.
   */
  void allocateCoordsOnSidre();

#endif

  axom::deprecated::MCArray<double>* m_coordinates[3] = {nullptr, nullptr, nullptr};

  DISABLE_COPY_AND_ASSIGNMENT(RectilinearMesh);
  DISABLE_MOVE_AND_ASSIGNMENT(RectilinearMesh);
};

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_RECTILINEARMESH_HPP_ */
