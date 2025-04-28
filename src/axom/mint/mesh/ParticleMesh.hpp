// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_PARTICLEMESH_HPP_
#define MINT_PARTICLEMESH_HPP_

#include "axom/core/Macros.hpp"  // for axom macros

#include "axom/mint/config.hpp"     // for mint compile-time definitions
#include "axom/mint/mesh/Mesh.hpp"  // for mint::Mesh base class

#include "axom/slic/interface/slic.hpp"  // for slic Macros

namespace axom
{
// Sidre Forward Declarations
namespace sidre
{
class Group;
}

namespace mint
{
// Mint Forward Declarations
class MeshCoordinates;

/*!
 * \class ParticleMesh
 *
 * \brief Provides the ability to store and operate on a set of particles.
 *
 *  The ParticleMesh class derives from the top-level Mesh base class and
 *  provides the ability to store and operate on a collection of particles and
 *  associated data. Each particle has a position and may have associated
 *  scalar, vector and tensor fields.
 *
 *  A ParticleMesh object may be constructed using (a) native storage, (b)
 *  external storage, or, (c) from a Sidre blueprint conforming hierarchy:
 *
 *  * <b> Native Storage </b> <br />
 *
 *    When using native storage, the ParticleMesh object owns all memory
 *    associated with the particle data. The storage can grow dynamically as
 *    needed by the application, i.e., adding more particles. Once the
 *    ParticleMesh object goes out-of-scope, all memory associated with it is
 *    returned to the system.
 *
 *  * <b> External Storage </b> <br />
 *
 *    A ParticleMesh may also be constructed from external, user-supplied
 *    buffers. In this case, all memory associated with the particle data
 *    is owned by the caller. Consequently, the number of particles and
 *    associated data cannot grow dynamically.
 *
 *  * <b> Sidre </b> <br />
 *
 *    A ParticleMesh may also be constructed from a Sidre hierarchy that is
 *    conforming to the <a href="http://llnl-conduit.readthedocs.io/en/latest/">
 *    mesh blueprint </a> conventions. In this case, all operations are
 *    supported, including dynamically adding new particles and growing the
 *    associated storage. However, Sidre owns all the memory. Once the
 *    ParticleMesh object goes out-of-scope, the data remains persistent in
 *    Sidre.
 *
 * \see mint::Mesh
 */
class ParticleMesh : public Mesh
{
public:
  /*!
   * \brief Default constructor. Disabled.
   */
  ParticleMesh() = delete;

  /// \name Native Storage Constructors
  /// @{

  /*!
   * \brief Constructs a ParticleMesh instance of specified dimension that
   *  holds the specified number of particles.
   *
   * \param [in] dimension the ambient dimension of the particle mesh.
   * \param [in] numParticles the number of particles in this
   * \param [in] capacity max particle capacity (optional)
   *
   * \pre 1 <= dimension <= 3
   * \pre numParticles >= 0
   *
   * \post getNumParticles() == numParticles
   * \post getNumParticles() <= capacity()
   * \post hasSidreGroup() == false
   */
  ParticleMesh(int dimension, IndexType numParticles, IndexType capacity = USE_DEFAULT);

  /// @}

  /// \name External Storage Constructors
  /// @{

  /*!
   * \brief Creates a ParticleMesh instance that points to the supplied external
   *  particle position buffers.
   *
   * \param [in] numParticles the number of particles in the supplied buffers.
   * \param [in] x pointer to the particle x-coordinate positions
   * \param [in] y pointer to the particle y-coordinate positions (optional)
   * \param [in] z pointer to the particle z-coordinate positions (optional)
   *
   * \note This constructor wraps the supplied particle position buffers.
   *  Consequently, the resulting ParticleMesh object does not own the memory
   *  associated with the particle positions.
   *
   * \warning All calls to shrink(), append(), resize() and reserve will fail
   *  on a ParticleMesh instance that is constructed using this constructor.
   *
   * \note The supplied buffers must have sufficient storage for numParticles
   *
   * \pre x != nullptr
   * \pre y != nullptr, if dimension==2 || dimension==3
   * \pre z != nullptr, if dimension==3
   * \pre numParticles >= 1
   *
   * \post 1 <= getDimension() <= 3
   * \post getNumParticles() == numParticles
   */
  ParticleMesh(IndexType numParticles, double* x, double* y = nullptr, double* z = nullptr);

  /// @}

#ifdef AXOM_MINT_USE_SIDRE
  /// \name Sidre Storage Constructors
  /// @{

  /*!
   * \brief Creates a ParticleMesh instance from a given Sidre group that holds
   *  particle mesh data that according to the conventions of the computational
   *  mesh blueprint.
   *
   * \param [in] group pointer to the blueprint root group in Sidre
   * \param [in] topo the name of the associated topology (optional)
   *
   * \note The supplied group is expected to contain particle mesh data
   *  that is valid and conforming to the conventions described in the mesh
   *  <a href="http://llnl-conduit.readthedocs.io/en/latest/"> blueprint </a>.
   *
   * \note If a topology name is not provided, the implementation will construct
   *  a mesh based on the 1st topology group under the parent "topologies"
   *  group.
   *
   * \note When using this constructor, all data is owned by Sidre. Once the
   *  ParticleMesh object goes out-of-scope, the data remains persistent in
   *  Sidre.
   *
   * \pre group != nullptr
   * \pre blueprint::isValidRootGroup( group )
   * \post hasSidreGroup() == true
   */
  explicit ParticleMesh(sidre::Group* group, const std::string& topo = "");

  /*!
   * \brief Creates a ParticleMesh object on an empty Sidre group.
   *
   * \param [in] dimension the ambient dimension of the particle mesh.
   * \param [in] numParticles the number of particles this instance holds
   * \param [in] group pointer to a group in Sidre where
   * \param [in] topo the name of the associated topology (optional)
   * \param [in] coordset the name of the coordset group in Sidre (optional)
   * \param [in] capacity max particle capacity (optional).
   *
   * \note If a topology and coordset name is not provided, internal defaults
   *  will be used by the implementation.
   *
   * \note When using this constructor, all data is owned by Sidre. Once the
   *  ParticleMesh object goes out-of-scope, the data remains persistent in
   *  Sidre.
   *
   * \pre 1 <= dimension <= 3
   * \pre group != nullptr
   * \pre group->getNumViews()==0
   * \pre group->getNumGroups()==0
   *
   * \post hasSidreGroup()==true
   */
  /// @{

  ParticleMesh(int dimension,
               IndexType numParticles,
               sidre::Group* group,
               const std::string& topo,
               const std::string& coordset,
               IndexType capacity = USE_DEFAULT);

  ParticleMesh(int dimension,
               IndexType numParticles,
               sidre::Group* group,
               IndexType capacity = USE_DEFAULT);
    /// @}

    /// @}

#endif /* AXOM_MINT_USE_SIDRE */

  /// \name Virtual methods
  /// @{

  /*!
   * \brief Destructor.
   */
  virtual ~ParticleMesh();

  /// \name Cells
  /// @{

  /*!
   * \brief Return the number of cells in the mesh.
   */
  virtual IndexType getNumberOfCells() const final override { return getNumberOfNodes(); }

  /*!
   * \brief Return the capacity for cells.
   */
  virtual IndexType getCellCapacity() const final override { return getNodeCapacity(); }

  virtual IndexType getNumberOfCellNodes(IndexType AXOM_UNUSED_PARAM(cellID) = 0) const final override
  {
    return 1;
  }

  virtual CellType getCellType(IndexType AXOM_UNUSED_PARAM(cellID) = 0) const final override
  {
    return VERTEX;
  }

  virtual IndexType getCellNodeIDs(IndexType cellID, IndexType* cell) const final override;

  /*!
   * \brief Return the number of faces associated with the given cell. For the
   *  ParticleMesh this is always zero.
   *
   * \param [in] cellID the ID of the cell in question.
   */
  virtual IndexType getNumberOfCellFaces(IndexType AXOM_UNUSED_PARAM(cellID) = 0) const final override
  {
    return 0;
  }

  /*!
   * \brief Populates the given buffer with the IDs of the faces of the given
   *  cell and returns the number of faces. Since the ParticleMesh has no faces
   *  this method errors out.
   * 
   * \param [in] cellID the ID of the cellID in question.
   * \param [out] faces buffer to populate with the face IDs. Must be of length
   *  at least getNumberOfCellFaces( cellID ).
   */
  virtual IndexType getCellFaceIDs(IndexType AXOM_UNUSED_PARAM(cellID),
                                   IndexType* AXOM_UNUSED_PARAM(faces)) const final override
  {
    SLIC_ERROR("ParticleMesh does not implement this method.");
    return 0;
  }

  /// @}

  /// \name Nodes
  /// @{

  /*!
   * \brief Return the number of nodes in the mesh.
   */
  virtual IndexType getNumberOfNodes() const final override { return m_positions->numNodes(); }

  /*!
   * \brief Return the capacity for nodes.
   */
  virtual IndexType getNodeCapacity() const final override { return m_positions->capacity(); }

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
  virtual void getNode(IndexType nodeID, double* node) const final override
  {
    m_positions->getCoordinates(nodeID, node);
  }

  /*!
   * \brief Returns pointer to the particle positions in the specified dimension
   *
   * \param[in] dim the specified dimension
   * \return coord pointer
   *
   * \pre 1 <= dim <= 3
   * \post coord != nullptr
   */
  /// @{

  virtual double* getCoordinateArray(int dim) final override
  {
    return m_positions->getCoordinateArray(dim);
  }

  virtual const double* getCoordinateArray(int dim) const final override
  {
    return m_positions->getCoordinateArray(dim);
  }

  /// @}

  /// @}

  /// \name Faces
  /// @{

  /*!
   * \brief Return the number of faces in the mesh.
   */
  virtual IndexType getNumberOfFaces() const final override { return 0; }

  /*!
   * \brief Return the type of the given face.
   *
   * \param [in] faceID the ID of the face in question.
   * 
   * \note The particle mesh does not have any faces so this call errors out.
   */
  virtual CellType getFaceType(IndexType AXOM_UNUSED_PARAM(faceID)) const final override
  {
    SLIC_ERROR("ParticleMesh does not implement this method.");
    return UNDEFINED_CELL;
  }

  /*!
   * \brief Return the number of nodes associated with the given face.
   *
   * \param [in] faceID the ID of the face in question.
   * 
   * \note The particle mesh does not have any faces so this call errors out.
   */
  virtual IndexType getNumberOfFaceNodes(IndexType AXOM_UNUSED_PARAM(faceID)) const final override
  {
    SLIC_ERROR("ParticleMesh does not implement this method.");
    return -1;
  }

  /*!
   * \brief Copy the IDs of the nodes that compose the given face into the
   *  provided buffer.
   *
   * \param [in] faceID the ID of the face in question.
   * \param [out] nodes the buffer into which the node IDs are copied, must
   *  be of length at least getNumberOfFaceNodes().
   *
   * \return The number of nodes for the given face, which is zero for the
   *  ParticleMesh.
   * 
   * \note The particle mesh does not have any faces so this call errors out.
   */
  virtual IndexType getFaceNodeIDs(IndexType AXOM_UNUSED_PARAM(faceID),
                                   IndexType* AXOM_UNUSED_PARAM(nodes)) const final override
  {
    SLIC_ERROR("ParticleMesh does not implement this method.");
    return -1;
  }

  /*!
   * \brief Copy the IDs of the cells adjacent to the given face into the
   *  provided indices.
   *
   * \param [in] faceID the ID of the face in question.
   * \param [out] cellIDOne the ID of the first cell.
   * \param [out] cellIDTwo the ID of the second cell.
   *
   * \note The particle mesh does not have any faces so this call errors out.
   */
  virtual void getFaceCellIDs(IndexType AXOM_UNUSED_PARAM(faceID),
                              IndexType& AXOM_UNUSED_PARAM(cellIDOne),
                              IndexType& AXOM_UNUSED_PARAM(cellIDTwo)) const final override
  {
    SLIC_ERROR("ParticleMesh does not implement this method.");
  }

  /// @}

  /// \name Edges
  /// @{

  /*!
   * \brief Return the number of edges in the mesh.
   */
  virtual IndexType getNumberOfEdges() const final override { return 0; }

  /*!
   * \brief Return the capacity for edges.
   */
  virtual IndexType getEdgeCapacity() const final override { return 0; }

  /// @}

  /*!
   * \brief Return true iff particle positions are stored in external arrays.
   * \return status true iff the particle positions point to external buffers.
   */
  virtual bool isExternal() const final override { return m_positions->isExternal(); }

  /// @}

  /// \name Attribute get/set Methods
  /// @{

  /// \name Nodes
  /// @{

  /*!
   * \brief Return the node resize ratio.
   */
  double getNodeResizeRatio() const { return m_positions->getResizeRatio(); }

  /*!
   * \brief Increase the number of particles this ParticleMesh instance can hold
   * \param [in] newSize the number of particles this instance will now hold
   * \post getNumParticles() == newSize
   */
  void resize(IndexType newSize);

  /*!
   * \brief Increase the max particle capacity of this ParticleMesh instance
   * \param [in] newCapacity
   */
  void reserve(IndexType newCapacity);

  /*!
   * \brief Shrinks the max particle capacity to the actual number of particles.
   *
   * \post getNumberOfNodes() == capacity()
   * \post f->getCapacity() == getNumberOfNodes() for all particle fields.
   */
  void shrink();

  /// @}

  /*!
   * \brief Return true iff the mesh holds no particles.
   */
  bool empty() const { return m_positions->empty(); }

  /*!
   * \brief Return true iff the particle positions are stored in sidre.
   */
  bool isInSidre() const { return m_positions->isInSidre(); }

  /// \name Data Access Methods
  /// @{

  /*!
   * \brief Appends a new particle to the ParticleMesh
   *
   * \param [in] x the x-coordinate of the particle position
   * \param [in] y the y-coordinate of the particle position (valid in 2-D)
   * \param [in] z the z-coordinate of the particle position (valid in 3-D)
   *
   * \post increments the number of particles by one.
   */
  /// @{

  void append(double x);
  void append(double x, double y);
  void append(double x, double y, double z);

  /// @}

  /// @}

private:
  /*!
   * \brief Helper method to initialize a ParticleMesh instance.
   * \note Called from the constructor.
   */
  void initialize();

  /*!
   * \brief Checks if the internal data array are consistent.
   * \return status true if the consistency checks pass, else, false.
   */
  bool checkConsistency();

  MeshCoordinates* m_positions;

  DISABLE_COPY_AND_ASSIGNMENT(ParticleMesh);
  DISABLE_MOVE_AND_ASSIGNMENT(ParticleMesh);
};

//------------------------------------------------------------------------------
// IN-LINE METHOD IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline IndexType ParticleMesh::getCellNodeIDs(IndexType cellID, IndexType* cell) const
{
  SLIC_ASSERT(cell != nullptr);
  SLIC_ASSERT(0 <= cellID && cellID <= getNumberOfCells());
  cell[0] = cellID;
  return 1;
}

//------------------------------------------------------------------------------
inline void ParticleMesh::append(double x)
{
  SLIC_ASSERT(m_positions != nullptr);
  SLIC_ERROR_IF(m_ndims != 1, "ParticleMesh::append(x) is only valid in 1-D");

  m_positions->append(x);
  m_mesh_fields[NODE_CENTERED]->resize(m_positions->numNodes());

  SLIC_ASSERT(checkConsistency());
}

//------------------------------------------------------------------------------
inline void ParticleMesh::append(double x, double y)
{
  SLIC_ASSERT(m_positions != nullptr);
  SLIC_ERROR_IF(m_ndims != 2, "ParticleMesh::append(x,y) is only valid in 2-D");

  m_positions->append(x, y);
  m_mesh_fields[NODE_CENTERED]->resize(m_positions->numNodes());

  SLIC_ASSERT(checkConsistency());
}

//------------------------------------------------------------------------------
inline void ParticleMesh::append(double x, double y, double z)
{
  SLIC_ASSERT(m_positions != nullptr);
  SLIC_ERROR_IF(m_ndims != 3, "ParticleMesh::append(x,y,z) is only valid in 3-D");

  m_positions->append(x, y, z);
  m_mesh_fields[NODE_CENTERED]->resize(m_positions->numNodes());

  SLIC_ASSERT(checkConsistency());
}

//------------------------------------------------------------------------------
inline void ParticleMesh::resize(IndexType newSize)
{
  SLIC_ASSERT(m_positions != nullptr);
  m_positions->resize(newSize);
  m_mesh_fields[NODE_CENTERED]->resize(newSize);
  SLIC_ASSERT(checkConsistency());
}

//------------------------------------------------------------------------------
inline void ParticleMesh::reserve(IndexType newCapacity)
{
  SLIC_ASSERT(m_positions != nullptr);
  m_positions->reserve(newCapacity);
  m_mesh_fields[NODE_CENTERED]->reserve(newCapacity);
}

//------------------------------------------------------------------------------
inline void ParticleMesh::shrink()
{
  SLIC_ASSERT(m_positions != nullptr);
  m_positions->shrink();
  m_mesh_fields[NODE_CENTERED]->shrink();
}

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_PARTICLEMESH_HPP_ */
