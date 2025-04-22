// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_MESH_HPP_
#define MINT_MESH_HPP_

#include "axom/core/Macros.hpp"  // for Axom macros

#include "axom/mint/mesh/CellTypes.hpp"         // for CellType enum
#include "axom/mint/config.hpp"                 // for mint compile-time type
#include "axom/mint/mesh/FieldAssociation.hpp"  // for FieldAssociation enum
#include "axom/mint/mesh/FieldData.hpp"         // for mint::FieldData
#include "axom/mint/mesh/MeshCoordinates.hpp"   // for mint::MeshCoordinates
#include "axom/mint/mesh/MeshTypes.hpp"         // for MeshType enum

#include "axom/slic/interface/slic.hpp"  // for SLIC macros

namespace axom
{
/* Forward declarations */
#ifdef AXOM_MINT_USE_SIDRE
namespace sidre
{
class Group;
}
#endif

namespace mint
{
/* Forward declarations */
class Field;
class Mesh;

/// \name Free Methods
/// @{

#ifdef AXOM_MINT_USE_SIDRE

/*!
 * \brief Creates a mesh instance from the given Sidre group.
 *
 * \param [in] group pointer to the root group of the mesh in Sidre.
 * \param [in] topo topology associated with the requested mesh (optional)
 *
 * \return m pointer to a mesh instance corresponding to the specified group.
 *
 * \note If a topology name is not provided, the implementation will construct
 *  a mesh based on the 1st topology group under the parent "topologies"
 *  group.
 *
 * \note Ownership of the resulting mesh object is passed to the caller.
 *
 * \note When using Mint with Sidre, Sidre maintains ownership of the data.
 *  Although the data can be modified via calls to Mint, groups and views
 *  cannot be deleted. The data remains persistent in Sidre once the mesh
 *  object goes out-of-scope.
 *
 * \pre  group != nullptr
 * \pre  blueprint::isValidRootGroup( group ) == true
 * \post m != nullptr
 */
Mesh* getMesh(sidre::Group* group, const std::string& topo = "");

#endif

/// @}

/*!
 * \class Mesh
 *
 * \brief Base class that defines the core API common to all Mesh types.
 *
 * A Mesh, \f$ \mathcal{M}(\Omega) \f$, provides an approximation of a physical
 * domain, \f$ \Omega \in \mathcal{R}^d \f$, where \f$ d \in [1,3] \f$ . The
 * Mesh is essentially a discrete representation of a problem and is used to
 * facilitate the analysis, e.g., FEA, CFD, etc. The solution domain is
 * approximated by dividing it into a finite number of <em> nodes </em> and
 * <em> cells </em> at which the variables of the underlying mathematical model
 * (i.e., a PDE) are then computed via a numerical method, such as,
 * Finite Difference (FD), Finite Volume (FV) or Finite Element (FE), chief
 * among them.
 *
 * There are a variety of mesh types. Mint supports the following mesh types:
 *
 *  * <b> Structured (Curvilinear) Mesh </b> <br />
 *
 *    A <em> structured mesh </em> divides the solution domain according to a
 *    logical grid where each node/cell of the mesh can be uniquely identified
 *    by a corresponding logical ijk-index. The nodes of the mesh are found at
 *    the intersection of the grid lines, but, are explicitly defined via a
 *    mapping onto the physical Cartesian coordinate system of the domain. For
 *    this reason, these types of meshes are also called <em> mapped </em> or
 *    <em> body-fitted </em> meshes.
 *
 *    However, the mesh topology (e.g., connectivity, neighbor information) is
 *    implicitly defined by the logical indexing scheme. For example,  a
 *    structured mesh is composed of <em> quadrilateral </em> cells in 2-D and
 *    <em> hexahedron </em> cells in 3-D. Given the logical index of a cell
 *    one can compute the node indices of the cell and neighbor information by
 *    performing simple shift operations on the associated index.
 *
 *  * <b> Rectilinear Mesh </b> <br />
 *
 *    A <em> rectilinear mesh </em>, also known as a product mesh, is similar
 *    to the <em> structured mesh </em> in that it is also defined according to
 *    a logical indexing scheme and has implicit topology.
 *
 *    However, the nodes and cells on a <em> rectilinear mesh </em> are arranged
 *    on a regular lattice that is axis-aligned with the Cartesian coordinate
 *    system. In contrast to the general <em> structured mesh </em>, the
 *    <em> rectilinear mesh </em> does not explicitly define **all** the nodes
 *    of the mesh. Instead, the nodes are only defined along each coordinate
 *    axis and may have variable spacing between nodes. Given a logical index,
 *    the corresponding physical position of a node can be evaluated by taking
 *    the cartesian product of the corresponding coordinate along each
 *    coordinate axis.
 *
 *  * <b> Uniform Mesh </b> <br />
 *
 *    A <em> uniform mesh </em>, also called a regular mesh, subdivides the
 *    domain in cells that have uniform spacing across each coordinate axis.
 *    Similar to the <em> structured mesh </em>, a <em> uniform mesh </em>
 *    adheres to the same logical indexing scheme and implicit topology
 *    representation. Moreover, The nodes and cells of a <em> uniform </em> mesh
 *    are arranged on a regular lattice as with the <em> rectilinear mesh </em>.
 *    However, both topology and geometry is implicitly defined on a
 *    <em> uniform mesh </em>. The geometry is solely defined by an origin,
 *    \f$ \hat{x_0} \f$ and spacing, \f$ \hat{h} \f$, along each axis. The
 *    coordinates of a node can be evaluated algebraically by the following:
 *    \f$ \hat{p} = \hat{x_0} + \hat{i} \times \hat{h} \f$, where \f$\hat{i}\f$
 *    is the logical ijk-index of the corresponding node.
 *
 *  * <b> Unstructured Mesh </b> <br />
 *
 *    An <em> unstructured mesh </em> stores both node and topology information
 *    explicitly. This allows the flexibility of discretizing the solution
 *    domain using a variety of cell types not just quadrilateral (in 2D) or
 *    hexahedral (in 3D) cells. Due to this added flexibility, the use of
 *    <em> unstructured meshes </em> is more common when dealing with complex
 *    geometries. However, <em> unstructured meshes </em> require additional
 *    storage and generally incur some performance penalty to store, create and
 *    access mesh topology information respectively.
 *
 *    Mint classifies <em> unstructured meshes </em> in two basic types based on
 *    the underlying mesh topology:
 *
 *    * <b> Single Cell Topology </b>
 *
 *      In this case, the <em> unstructured mesh </em> consists of a single cell
 *      type, e.g., a quad or triangle mesh in 2D, or, a hex or tet mesh in 3D.
 *      In this case  the underlying implementation is optimized for the
 *      specified cell type (specified in the constructor).
 *
 *    * <b> Mixed Cell Topology </b>
 *
 *      When <em> mixed cell topology </em> is specified, the <em> unstructured
 *      mesh </em> can be composed of any of the supported cell types, e.g.,
 *      a mesh consisting of both quads and triangles. This mode incurs
 *      additional overhead for storage and access to mesh topology information,
 *      since it requires indirection.
 *
 *    The list of supported cell types for an <em> unstructured mesh </em> is
 *    available in CellTypes.hpp
 *
 *  * <b> Particle Mesh </b> <br />
 *
 *    A <em> particle mesh </em> discretizes the solution domain using
 *    a set of discrete particle elements which correspond to the the nodes
 *    of the mesh. There is no ordering imposed on the particles and the
 *    coordinates of the particles are explicitly defined. A particle mesh
 *    has no connectivity information connecting the particles, which is why
 *    in literature methods using a particle discretization are also referred
 *    to as <em> meshless </em> or <em> meshfree </em> methods.
 *
 * The Mesh class provides the means to create, access and remove fields on
 * a mesh given the field name and its association. The field association
 * designates the corresponding mesh entity at which the field is stored, e.g.
 * whether the field is stored at the nodes or cell centers. A Field may be a
 * scalar quantity, e.g., pressure, or a vector field, such as, velocity.
 *
 * \warning When using Sidre, field names have to be unique. For example, if
 *  there exists a "pressure" node-centered field, there cannot be a
 *  corresponding cell-centered field.
 *
 * \note Mesh is a base class and cannot be instantiated directly
 *
 * \note Typically, the computational mesh can be defined across one or more
 *  <em> blocks </em>, e.g., for multi-block problems, where each block is
 *  subsequently decomposed into several <em> partitions </em> for parallel
 *  computation. A Mesh instance represents a <em> single </em> partition for
 *  a given block.
 *
 * \see mint::UnstructuredMesh
 * \see mint::StructuredMesh
 * \see mint::CurvilinearMesh
 * \see mint::RectilinearMesh
 * \see mint::UniformMesh
 * \see mint::Field
 * \see mint::FieldData
 * \see mint::MeshTypes
 */
class Mesh
{
public:
  /*!
   * \brief Default constructor. Disabled.
   */
  Mesh() = delete;

  /// \name Virtual methods
  /// @{

  /*!
   * \brief Destructor.
   */
  virtual ~Mesh();

  /// \name Cells
  /// @{

  /*!
   * \brief Returns the number of cells in this mesh instance.
   * \return N the number of cells
   * \post N >= 0
   */
  virtual IndexType getNumberOfCells() const = 0;

  /*!
   * \brief Returns the capacity for number of cell in this mesh instance.
   * \return N the cell capacity
   * \post N >= 0
   */
  virtual IndexType getCellCapacity() const { return getNumberOfCells(); }

  /*!
   * \brief Return the type of the given cell.
   *
   * \param [in] cellID the ID of the cell in question, this parameter is
   *  ignored if hasMixedCellTypes() == false.
   *
   * \pre 0 <= cellID < getNumberOfCells()
   */
  virtual CellType getCellType(IndexType cellID = 0) const = 0;

  /*!
   * \brief Return the number of nodes associated with the given cell.
   *
   * \param [in] cellID the ID of the cell in question, this parameter is
   *  ignored unless hasMixedCellTypes() == true.
   *
   * \pre 0 <= cellID < getNumberOfCells()
   */
  virtual IndexType getNumberOfCellNodes(IndexType cellID = 0) const = 0;

  /*!
   * \brief Copy the connectivity of the given cell into the provided buffer.
   *  The buffer must be of length at least getNumberOfCellNodes( cellID ).
   *
   * \param [in] cellID the ID of the cell in question.
   * \param [out] nodes the buffer into which the connectivity is copied, must
   *  be of length at least getNumberOfCellNodes( cellID ).
   *
   * \return The number of nodes for the given cell.
   *
   * \pre nodes != nullptr
   * \pre 0 <= cellID < getNumberOfCells()
   */
  virtual IndexType getCellNodeIDs(IndexType AXOM_UNUSED_PARAM(cellID),
                                   IndexType* AXOM_UNUSED_PARAM(nodes)) const = 0;

  /*!
   * \brief Return the number of faces associated with the given cell.
   *
   * \param [in] cellID the ID of the cell in question.
   */
  virtual IndexType getNumberOfCellFaces(IndexType AXOM_UNUSED_PARAM(cellID) = 0) const = 0;

  /*!
   * \brief Populates the given buffer with the IDs of the faces of the given
   *  cell and returns the number of faces.
   *
   * \param [in] cellID the ID of the cellID in question.
   * \param [out] faces buffer to populate with the face IDs. Must be of length
   *  at least getNumberOfCellFaces( cellID ).
   *
   * \pre faces != nullptr
   * \pre 0 <= cellID < getNumberOfCells()
   */
  virtual IndexType getCellFaceIDs(IndexType AXOM_UNUSED_PARAM(cellID),
                                   IndexType* AXOM_UNUSED_PARAM(faces)) const = 0;

  /// @}

  /// \name Nodes
  /// @{

  /*!
   * \brief Returns the number of nodes in this mesh instance.
   * \return N the number of nodes
   * \post N >= 0
   */
  virtual IndexType getNumberOfNodes() const = 0;

  /*!
   * \brief Returns the capacity for number of nodes in this mesh instance.
   * \return N the node capacity
   * \post N >= 0
   */
  virtual IndexType getNodeCapacity() const { return getNumberOfNodes(); }

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
  virtual void getNode(IndexType nodeID, double* node) const = 0;

  /*!
   * \brief Returns pointer to the requested mesh coordinate buffer.
   *
   * \param [in] dim the dimension of the requested coordinate buffer
   * \return ptr pointer to the coordinate buffer.
   *
   * \note if hasExplicitCoordinates() == true then the length of the returned
   *  buffer is getNumberOfNodes(). Otherwise the UniformMesh returns
   *  nullptr and the RectilinearMesh returns a pointer to the associated
   *  dimension scale which is of length
   *  static_cast< RectilinearMesh* >( this )->getNodeResolution().
   *
   * \pre dim >= 0 && dim < dimension()
   * \pre dim == X_COORDINATE || dim == Y_COORDINATE || dim == Z_COORDINATE
   */
  /// @{
  virtual double* getCoordinateArray(int dim) = 0;
  virtual const double* getCoordinateArray(int dim) const = 0;
  /// @}

  /// @}

  /// \name Faces
  /// @{

  /*!
   * \brief Returns the number of faces in this mesh instance.
   * \return N the number of faces
   * \post N >= 0
   */
  virtual IndexType getNumberOfFaces() const = 0;

  /*!
   * \brief Returns the capacity for number of faces in this mesh instance.
   * \return N the face capacity
   * \post N >= 0
   */
  virtual IndexType getFaceCapacity() const { return getNumberOfFaces(); }

  /*!
   * \brief Return the type of the given face.
   *
   * \param [in] faceID the ID of the face in question.
   */
  virtual CellType getFaceType(IndexType AXOM_UNUSED_PARAM(faceID)) const = 0;

  /*!
   * \brief Return the number of nodes associated with the given face.
   *
   * \param [in] faceID the ID of the face in question.
   */
  virtual IndexType getNumberOfFaceNodes(IndexType AXOM_UNUSED_PARAM(faceID)) const = 0;

  /*!
   * \brief Copy the IDs of the nodes that compose the given face into the
   *  provided buffer.
   *
   * \param [in] faceID the ID of the face in question.
   * \param [out] nodes the buffer into which the node IDs are copied, must
   *  be of length at least getNumberOfFaceNodes().
   *
   * \return The number of nodes for the given face.
   *
   * \pre nodes != nullptr
   * \pre 0 <= faceID < getNumberOfFaces()
   */
  virtual IndexType getFaceNodeIDs(IndexType AXOM_UNUSED_PARAM(faceID),
                                   IndexType* AXOM_UNUSED_PARAM(nodes)) const = 0;

  /*!
   * \brief Copy the IDs of the cells adjacent to the given face into the
   *  provided indices.
   *
   * \param [in] faceID the ID of the face in question.
   * \param [out] cellIDOne the ID of the first cell.
   * \param [out] cellIDTwo the ID of the second cell.
   *
   * \note If no cell exists (the face is external) then the ID will be set to
   * -1.
   *
   * \pre 0 <= faceID < getNumberOfFaces()
   */
  virtual void getFaceCellIDs(IndexType AXOM_UNUSED_PARAM(faceID),
                              IndexType& AXOM_UNUSED_PARAM(cellIDOne),
                              IndexType& AXOM_UNUSED_PARAM(cellIDTwo)) const = 0;

  /// @}

  /// \name Edges
  /// @{

  /*!
   * \brief Returns the number of edges in this mesh instance.
   * \return N the number of edges
   * \post N >= 0
   */
  virtual IndexType getNumberOfEdges() const = 0;

  /*!
   * \brief Returns the capacity for number of edges in this mesh instance.
   * \return N the edge capacity
   * \post N >= 0
   */
  virtual IndexType getEdgeCapacity() const { return getNumberOfEdges(); }

  /*!
   * \brief Returns true iff the mesh was constructed with external arrays.
   * \return status true if the mesh points to external buffers, else, false.
   */
  virtual bool isExternal() const = 0;

  /// @}

  /// @}

  /// \name Mesh Attribute get/set Methods
  /// @{

  /*!
   * \brief Returns the dimension for this mesh instance.
   * \return ndims the dimension of this mesh instance.
   * \post ndims >= 1 && ndims <= 3
   */
  inline int getDimension() const { return m_ndims; }

  /*!
   * \brief Returns the ID of this mesh instance.
   * \return Id the ID of the mesh.
   */
  inline int getBlockId() const { return m_block_idx; }

  /*!
   * \brief set the block ID of this mesh instance.
   *
   * \param [in] ID the new block ID.
   *
   * \post getBlockId() == ID
   */
  void setBlockId(int ID);

  /*!
   * \brief Returns the partition ID of this mesh instance.
   * \return partitionId the partition ID of the mesh.
   */
  inline int getPartitionId() const { return m_part_idx; }

  /*!
   * \brief set the partition ID of this mesh instance.
   *
   * \param [in] ID the new partition ID.
   *
   * \post getPartitionId() == ID
   */
  void setPartitionId(int ID);

  /*!
   * \brief Returns the mesh type of this mesh instance.
   * \return meshType the mesh type
   * \see MeshType
   */
  inline int getMeshType() const { return m_type; }

  /*!
   * \brief Checks if this mesh instance has explicit coordinates.
   * \return status true iff the mesh defines coordinates explicitly.
   */
  inline bool hasExplicitCoordinates() const { return m_explicit_coords; }

  /*!
   * \brief Checks if this mesh instance has explicit connectivity.
   * \return status true iff the mesh defines cell connectivity explicitly.
   */
  inline bool hasExplicitConnectivity() const { return m_explicit_connectivity; }

  /*!
   * \brief Checks if the mesh has mixed cell types, e.g., consisting of both
   *  triangle and quad elements or hex,pyramid,prisms and tets in 3-D.
   * \return status true iff the mesh has mixed cell types.
   */
  inline bool hasMixedCellTypes() const { return m_has_mixed_topology; }

  /*!
   * \brief Returns true if the mesh type is structured.
   * \return status true if the mesh type is structured, else, false.
   */
  inline bool isStructured() const
  {
    return ((m_type == STRUCTURED_CURVILINEAR_MESH) || (m_type == STRUCTURED_RECTILINEAR_MESH) ||
            (m_type == STRUCTURED_UNIFORM_MESH));
  }

  /*!
   * \brief Returns true if the mesh type is unstructured.
   * \return status true if the mesh type is unstructured, else, false.
   */
  inline bool isUnstructured() const { return (m_type == UNSTRUCTURED_MESH); }

  /*!
   * \brief Checks if this Mesh instance is associated with a Sidre Group.
   * \return status true if the Mesh is associated with a group in a Sidre
   *  hierarchy, else, false.
   */
  inline bool hasSidreGroup() const;

#ifdef AXOM_MINT_USE_SIDRE
  /*!
   * \brief Return a pointer to the sidre::Group associated with this Mesh
   *  instance or nullptr if none exists.
   */
  inline sidre::Group* getSidreGroup() { return m_group; }

  /*!
   * \brief Return a const pointer to the sidre::Group associated with this Mesh
   *  instance or nullptr if none exists.
   */
  inline const sidre::Group* getSidreGroup() const { return m_group; }

  /*!
   * \brief Return the name of the topology associated with this Mesh instance,
   *  the return value is undefined if the mesh is not in sidre.
   */
  inline const std::string& getTopologyName() const { return m_topology; }

  /*!
   * \brief Return the name of the coordset associated with this Mesh instance,
   *  the return value is undefined if the mesh is not in sidre.
   */
  inline const std::string& getCoordsetName() const { return m_coordset; }
#endif

  /// @}

  /// \name Methods to Create, Access & Remove Fields from a Mesh
  /// @{

  /*!
   * \brief Returns const pointer to the FieldData instance with the specified
   *  mesh field association, e.g., NODE_CENTERED, CELL_CENTERED, etc.
   *
   * \param [in] association the specified mesh field association
   * \return fd pointer to the requested FieldData instance
   *
   * \pre association >= 0 && association < NUM_FIELD_ASSOCIATION
   * \post fd != nullptr
   *
   * \see FieldAssociation
   * \see FieldData
   */
  inline const FieldData* getFieldData(int association) const;

  /*!
   * \brief Check if a field with the given name and association exists.
   *
   * \param [in] name the name of the field in query.
   * \param [in] association the field association (optional)
   *
   * \return status true if the field exists, else, false.
   *
   * \note If an association is not explicitly specified, the code will check
   *  if a field by the given name exists in any available centeering.
   *
   * \pre name.empty()==false
   * \pre association >= 0 && association < NUM_FIELD_ASSOCIATION
   *
   * \see FieldAssociation
   */
  inline bool hasField(const std::string& name, int association = ANY_CENTERING) const;

  /*!
   * \brief Creates a new field with the given name and specified mesh field
   *  association, e.g., NODE_CENTERED, CELL_CENTERED, etc.
   *
   * \param [in] name the name of the new field.
   * \param [in] association the mesh field association.
   * \param [in] num_components number of components of the field (optional).
   * \param [in] storeInSidre indicates whether to store the field in the
   *  corresponding Sidre group (optional).
   * \param [in] capacity
   *
   * \return ptr raw pointer to the data buffer of the new field.
   *
   * \note This method throws an error and aborts if any of the pre-conditions
   *  is not satisfied.
   *
   * \pre name.empty() == false
   * \pre hasField( name ) == false
   * \pre association >= 0 && association < NUM_FIELD_ASSOCIATION
   *
   * \post ptr != nullptr
   * \post hasField( name ) == true
   *
   * \see FieldAssociation
   */
  template <typename T>
  inline T* createField(const std::string& name,
                        int association,
                        IndexType num_components = 1,
                        bool storeInSidre = true);

  /*!
   * \brief Creates a new field from an external buffer that has the given name
   *  and specified mesh field association, e.g., NODE_CENTERED, CELL_CENTERED,
   *  etc.
   *
   * \param [in] name the name of the new field.
   * \param [in] association the mesh field association.
   * \param [in] data pointer to the external data buffer.
   * \param [in] num_components number of components of the field (optional).
   *
   * \return ptr raw pointer to the data buffer of the new field.
   *
   * \note This method throws an error and aborts if any of the pre-conditions
   *  is not satisfied.
   *
   * \pre name.empty() == false
   * \pre hasField( name ) == false
   * \pre data != nullptr
   * \pre association >= 0 && association < NUM_FIELD_ASSOCIATION
   *
   * \post ptr != nullptr
   * \post ptr == data
   * \post hasField( name ) == true
   *
   * \see FieldAssociation
   */
  template <typename T>
  inline T* createField(const std::string& name,
                        int association,
                        T* data,
                        IndexType num_components = 1,
                        IndexType capacity = USE_DEFAULT);

  /*!
   * \brief Removes the field with the given name and specified association.
   *
   * \param [in] name the name of the field to remove.
   * \param [in] association the mesh field association.
   *
   * \return status true if the field is removed successfully, else, false.
   *
   * \pre name.emtpy() == false
   * \pre association >= 0 && association < NUM_FIELD_ASSOCIATION
   *
   * \see FieldAssociation
   */
  inline bool removeField(const std::string& name, int association);

  /*!
   * \brief Returns pointer to buffer of the field with the given ane and
   *  specified mesh field association.
   *
   * \param [in] name the name of the requested field.
   * \param [in] association the mesh field association.
   * \param [out] num_components the number of components per tuple (optional).
   *
   * \return ptr raw pointer to the data buffer of the requested field.
   *
   * \pre name.empty() == false
   * \pre hasField( name )
   * \pre association >= 0 && association < NUM_FIELD_ASSOCIATION
   *
   * \see FieldAssociation
   */
  /// @{
  template <typename T>
  inline T* getFieldPtr(const std::string& name, int association, IndexType& num_components);

  template <typename T>
  inline T* getFieldPtr(const std::string& name, int association);

  template <typename T>
  inline const T* getFieldPtr(const std::string& name, int association, IndexType& num_components) const;

  template <typename T>
  inline const T* getFieldPtr(const std::string& name, int association) const;
  /// @}

  /// @}

protected:
  /// \name Protected Members
  /// @{

  int m_ndims;     /*! mesh dimension */
  int m_type;      /*! the type of the mesh */
  int m_block_idx; /*! the Block ID of the mesh */
  int m_part_idx;  /*! the partition ID of the mesh */

  bool m_explicit_coords;
  bool m_explicit_connectivity;
  bool m_has_mixed_topology;

  FieldData* m_mesh_fields[NUM_FIELD_ASSOCIATIONS];

#ifdef AXOM_MINT_USE_SIDRE
  sidre::Group* m_group;
  std::string m_topology;
  std::string m_coordset;
#endif

  /// @}

  /// \name Protected Constructors (used in derived classes )
  /// @{

  /*!
   * \brief Mesh Constructor.
   *
   * \param [in] ndims the number of dimensions
   * \param [in] type the mesh type.
   * \param [in] blockId the block ID for this mesh instance.
   * \param [in] partId the partition ID for this mesh instance.
   */
  Mesh(int ndims, int type);

#ifdef AXOM_MINT_USE_SIDRE

  /*!
   * \brief Constructor for use with a group that already has data.
   *
   * \param [in] group the sidre::Group to use.
   * \param [in] topo optional argument specifying the name of the topology
   *  associated with this Mesh instance.
   *
   * \note If a topology name is not provided, the implementation will construct
   *  a mesh based on the 1st topology group under the parent "topologies"
   *  group.
   *
   * \pre group != nullptr.
   * \pre blueprint::isValidRootGroup( group ) == true
   *
   * \see sidre::Group
   */
  Mesh(sidre::Group* group, const std::string& topo = "");

  /*!
   * \brief Constructor for use with an empty group.
   *
   * \param [in] ndims the number of dimensions
   * \param [in] type the mesh type.
   * \param [in] group the sidre::Group to use.
   * \param [in] topo the name of the associated topology group.
   * \param [in] coordset the name of the associated coordset group.
   *
   * \note If a topology and coordset name is not provided a default name is
   *  used by the implementation.
   *
   * \pre group != nullptr.
   * \pre group->getNumGroups() == 0
   * \pre group->getNumViews() == 0
   * \post blueprint::isValidRootGroup( group )
   *
   * \see sidre::Group
   */
  Mesh(int ndims, int type, sidre::Group* group, const std::string& topo, const std::string& coordset);

  /*!
   * \brief Helper method to return the associated coordset group.
   * \return coordset the associated coordset group.
   *
   * \pre  m_group != nullptr
   * \pre  blueprint::isValidRootGroup( m_group )
   * \post blueprint::isValidCoordsetGroup( coordset )
   */
  sidre::Group* getCoordsetGroup();

  /*!
   * \brief Helper method to return the associated topology group.
   * \return topology the associated topology group.
   *
   * \pre  m_group != nullptr
   * \pre  blueprint::isValidRootGroup( m_group )
   * \post blueprint::isValidTopologyGroup( topology )
   */
  sidre::Group* getTopologyGroup();

#endif

  /// @}

private:
  /*!
   * \brief Get the info corresponding to the given mesh field association.
   *
   * \param [in] association the mesh field association, e.g., NODE_CENTERED.
   * \param [out] num_tuples the number of tuples in the associated FieldData.
   * \param [out] capacity the capacity of the associated FieldData.
   */
  void getFieldInfo(int association, IndexType& num_tuples, IndexType& capacity) const;

  /*!
   * \brief Helper method to check if the mesh type is valid.
   * \return status true if the mesh type is valie, else, false.
   */
  inline bool validMeshType() const { return ((m_type >= 0) && (m_type < mint::NUM_MESH_TYPES)); }

  /*!
   * \brief Helper method to check if the mesh dimension is valid.
   * \return status true if the mesh dimension is valid, else, false.
   */
  inline bool validDimension() const { return (m_ndims >= 1 && m_ndims <= 3); }

  /*!
   * \brief Allocates the FieldData internal data-structures.
   * \note Helper method that is called from the constructor.
   */
  void allocateFieldData();

  /*!
   * \brief Deallocates the FieldData internal data-structures.
   * \note Helper method that is called by the destructor.
   */
  void deallocateFieldData();

  DISABLE_COPY_AND_ASSIGNMENT(Mesh);
  DISABLE_MOVE_AND_ASSIGNMENT(Mesh);
};

//------------------------------------------------------------------------------
//  IMPLEMENTATION OF TEMPLATE & IN-LINE METHODS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline bool Mesh::hasSidreGroup() const
{
#ifdef AXOM_MINT_USE_SIDRE
  return (m_group != nullptr);
#else
  return false;
#endif
}

//------------------------------------------------------------------------------
inline const FieldData* Mesh::getFieldData(int association) const
{
  SLIC_ERROR_IF(association < 0 || association >= NUM_FIELD_ASSOCIATIONS,
                "invalid field association [" << association << "]");
  SLIC_ERROR_IF(m_mesh_fields[association] == nullptr,
                "null field data object w/association [" << association << "]");
  SLIC_ERROR_IF(m_type == PARTICLE_MESH && association != NODE_CENTERED,
                "a particle mesh may only store node-centered fields");

  return m_mesh_fields[association];
}

//------------------------------------------------------------------------------
inline bool Mesh::hasField(const std::string& name, int association) const
{
  bool found = false;

  if(association == mint::ANY_CENTERING)
  {
    int N = (m_type == mint::PARTICLE_MESH) ? 1 : mint::NUM_FIELD_ASSOCIATIONS;
    for(int i = 0; !found && i < N; ++i)
    {
      const FieldData* fd = getFieldData(i);
      SLIC_ASSERT(fd != nullptr);
      found = fd->hasField(name);
    }
  }
  else
  {
    const FieldData* fd = getFieldData(association);
    SLIC_ASSERT(fd != nullptr);
    found = fd->hasField(name);
  }

  return (found);
}

//------------------------------------------------------------------------------
template <typename T>
inline T* Mesh::createField(const std::string& name,
                            int association,
                            IndexType num_components,
                            bool storeInSidre)
{
  SLIC_ERROR_IF(hasField(name), "a field with the same name already exists!");

  FieldData* fd = const_cast<FieldData*>(getFieldData(association));
  SLIC_ASSERT(fd != nullptr);

  IndexType num_tuples, capacity;
  getFieldInfo(association, num_tuples, capacity);
  T* ptr = fd->createField<T>(name, num_tuples, num_components, capacity, storeInSidre);
  if(num_tuples > 0)
  {
    SLIC_ASSERT(ptr != nullptr);
  }

  return (ptr);
}

//------------------------------------------------------------------------------
template <typename T>
inline T* Mesh::createField(const std::string& name,
                            int association,
                            T* data,
                            IndexType num_components,
                            IndexType capacity)
{
  SLIC_ERROR_IF(hasField(name), "a field with the same name already exists!");

  SLIC_ASSERT(data != nullptr);

  FieldData* fd = const_cast<FieldData*>(getFieldData(association));
  SLIC_ASSERT(fd != nullptr);

  IndexType num_tuples, dummy1;
  getFieldInfo(association, num_tuples, dummy1);
  T* ptr = fd->createField<T>(name, data, num_tuples, num_components, capacity);
  SLIC_ASSERT(ptr == data);

  return (ptr);
}

//------------------------------------------------------------------------------
inline bool Mesh::removeField(const std::string& name, int association)
{
  bool status = false;
  FieldData* fd = const_cast<FieldData*>(getFieldData(association));

  const bool hasField = fd->hasField(name);
  SLIC_WARNING_IF(!hasField, "field [" << name << "] does not exist!");

  if(hasField)
  {
    fd->removeField(name);
    status = true;
  }

  return (status);
}

//------------------------------------------------------------------------------
template <typename T>
inline T* Mesh::getFieldPtr(const std::string& name, int association)
{
  IndexType num_components = 0;
  return getFieldPtr<T>(name, association, num_components);
}

//------------------------------------------------------------------------------
template <typename T>
inline T* Mesh::getFieldPtr(const std::string& name, int association, IndexType& num_components)
{
  const Mesh* self = const_cast<const Mesh*>(this);
  const T* ptr = self->getFieldPtr<T>(name, association, num_components);
  return (const_cast<T*>(ptr));
}

//------------------------------------------------------------------------------
template <typename T>
inline const T* Mesh::getFieldPtr(const std::string& name, int association) const
{
  IndexType num_components = 0;
  return getFieldPtr<T>(name, association, num_components);
}

//------------------------------------------------------------------------------
template <typename T>
inline const T* Mesh::getFieldPtr(const std::string& name,
                                  int association,
                                  IndexType& num_components) const
{
  const FieldData* fd = getFieldData(association);
  SLIC_ASSERT(fd != nullptr);

  IndexType num_tuples = 0;
  const T* ptr = fd->getFieldPtr<T>(name, num_tuples, num_components);
  SLIC_ASSERT(ptr != nullptr);

  return (ptr);
}

//------------------------------------------------------------------------------
inline void Mesh::getFieldInfo(int association, IndexType& num_tuples, IndexType& capacity) const
{
  switch(association)
  {
  case NODE_CENTERED:
    num_tuples = getNumberOfNodes();
    capacity = getNodeCapacity();
    break;
  case CELL_CENTERED:
    num_tuples = getNumberOfCells();
    capacity = getCellCapacity();
    break;
  case FACE_CENTERED:
    num_tuples = getNumberOfFaces();
    capacity = getFaceCapacity();
    break;
  default:
    SLIC_ASSERT(association == EDGE_CENTERED);
    num_tuples = getNumberOfEdges();
    capacity = getEdgeCapacity();
    break;
  }  // END switch
}

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_MESH_HPP_ */
