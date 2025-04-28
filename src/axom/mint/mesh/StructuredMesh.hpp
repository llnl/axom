// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_STRUCTUREDMESH_HPP_
#define MINT_STRUCTUREDMESH_HPP_

#include "axom/core/Types.hpp"       // for axom types
#include "axom/core/Macros.hpp"      // for axom macros
#include "axom/core/StackArray.hpp"  // for StackArray

#include "axom/mint/config.hpp"          // for compile-time definitions
#include "axom/mint/mesh/CellTypes.hpp"  // for the CellTypes enum
#include "axom/mint/mesh/Mesh.hpp"       // for mint::Mesh base class

#include "axom/slic/interface/slic.hpp"  // for SLIC macros

#include <cstring>  // for std::memcpy

namespace axom
{
namespace mint
{
constexpr int I_DIRECTION = 0;
constexpr int J_DIRECTION = 1;
constexpr int K_DIRECTION = 2;

/*!
 * \class StructuredMesh
 *
 * \brief Base class that defines the core API common for structured mesh types.
 *
 *  The StructuredMesh class derives from the abstract Mesh base class and
 *  implements the core API for Structured meshes. Specifically, the
 *  StrucrturedMesh  defines and implements all the extent-based operations.
 *
 * \see UniformMesh
 * \see RectilinearMesh
 * \see CurvilinearMesh
 * \see Mesh
 */
class StructuredMesh : public Mesh
{
public:
  /*!
   * \brief Default constructor. Disabled.
   */
  StructuredMesh() = delete;

  /// \name Virtual methods
  /// @{

  /*!
   * \brief Destructor.
   */
  virtual ~StructuredMesh() { }

  /// \name Cells
  /// @{

  /*!
   * \brief Return the number of cells in the mesh.
   */
  virtual IndexType getNumberOfCells() const final override
  {
    IndexType n_cells = 1;
    for(int dim = 0; dim < m_ndims; ++dim)
    {
      n_cells *= getCellResolution(dim);
    }

    return n_cells;
  }

  /*!
   * \brief Return the type of cell this mesh holds. SEGMENT, QUAD, or HEX
   *  depending on the dimension.
   *
   * \param [in] cellID the ID of the cell in question, this parameter is
   *  ignored.
   */
  virtual CellType getCellType(IndexType AXOM_UNUSED_PARAM(cellID) = 0) const final override
  {
    return (m_ndims == 1) ? SEGMENT : (m_ndims == 2) ? QUAD : HEX;
  }

  /*!
   * \brief Return the number of nodes associated with the given cell.
   *
   * \param [in] cellID the ID of the cell in question, this parameter is
   *  ignored.
   */
  virtual IndexType getNumberOfCellNodes(IndexType AXOM_UNUSED_PARAM(cellID) = 0) const final override
  {
    return (m_ndims == 1) ? 2 : (m_ndims == 2) ? 4 : 8;
  }

  /*!
   * \brief Copy the node IDs of the given cell into the provided buffer.
   *  The buffer must be of length at least getNumberOfCellNodes().
   *
   * \param [in] cellID the ID of the cell in question.
   * \param [out] nodes the buffer into which the node IDs are copied, must
   *  be of length at least getNumberOfCellNodes().
   *
   * \return The number of nodes for the given cell.
   *
   * \pre nodes != nullptr
   * \pre 0 <= cellID < getNumberOfCells()
   */
  virtual IndexType getCellNodeIDs(IndexType cellID, IndexType* nodes) const final override;

  /*!
   * \brief Return the number of faces associated with the given cell.
   *
   * \param [in] cellID the ID of the cell in question, this parameter is
   *  ignored.
   */
  virtual IndexType getNumberOfCellFaces(IndexType AXOM_UNUSED_PARAM(cellID) = 0) const final override
  {
    CellType cell_type = getCellType();
    return getCellInfo(cell_type).num_faces;
  }

  /*!
   * \brief Returns the IDs of the faces of the cell at (i,j) or (i, j, k).
   *
   * \param [in] i logical index of the cell along the first dimension.
   * \param [in] j logical index of the cell along the second dimension.
   * \param [in] k logical index of the cell along the third dimension
   *  (optional).
   * \param [out] faces buffer to populate with the face IDs. Must be of length
   *  at least getNumberOfCellFaces().
   *
   * \note The faces are returned in the order of LOWER_I_FACE, UPPER_I_FACE,
   *  LOWER_J_FACE, UPPER_J_FACE and then LOWER_K_FACE, UPPER_K_FACE if 3D.
   *
   * \pre faces != nullptr
   * \pre 0 <= cellID < getNumberOfCells()
   */
  virtual IndexType getCellFaceIDs(IndexType cellID, IndexType* faces) const final override;

  /// @}

  /// \name Nodes
  /// @{

  /*!
   * \brief Return the number of nodes in the mesh.
   */
  virtual IndexType getNumberOfNodes() const final override
  {
    IndexType n_nodes = 1;
    for(int dim = 0; dim < m_ndims; ++dim)
    {
      n_nodes *= getNodeResolution(dim);
    }

    return n_nodes;
  }

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
  virtual void getNode(IndexType nodeID, double* node) const override = 0;

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
   *  static_cast< RectilinearMesh* >( this )->getNodeResolution( dim ).
   *
   * \pre dim >= 0 && dim < dimension()
   * \pre dim == X_COORDINATE || dim == Y_COORDINATE || dim == Z_COORDINATE
   */
  /// @{

  virtual double* getCoordinateArray(int dim) override = 0;
  virtual const double* getCoordinateArray(int dim) const override = 0;

  /// @}

  /// @}

  /// \name Faces
  /// @{

  /*!
   * \brief Return the number of faces in the mesh.
   */
  virtual IndexType getNumberOfFaces() const final override
  {
    return m_total_IJ_faces + getTotalNumFaces(2);
  }

  /*!
   * \brief Return the type of face this mesh holds. SEGMENT or QUAD depending
   * on the dimension.
   *
   * \param [in] faceID the ID of the face in question, this parameter is
   *  ignored.
   */
  virtual CellType getFaceType(IndexType AXOM_UNUSED_PARAM(faceID) = 0) const final override
  {
    return (m_ndims == 2) ? SEGMENT : (m_ndims == 3) ? QUAD : UNDEFINED_CELL;
  }

  /*!
   * \brief Return the number of nodes associated with the given face.
   *
   * \param [in] faceID the ID of the face in question, this parameter is
   *  ignored.
   */
  virtual IndexType getNumberOfFaceNodes(IndexType AXOM_UNUSED_PARAM(faceID) = 0) const final override
  {
    return (m_ndims == 2) ? 2 : (m_ndims == 3) ? 4 : 0;
  }

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
   * \pre 0 <= faceID < getNumberOfCells()
   */
  virtual IndexType getFaceNodeIDs(IndexType faceID, IndexType* nodes) const final override;

  /*!
   * \brief Copy the IDs of the cells adjacent to the given face into the
   *  provided indices.
   *
   * \param [in] faceID the ID of the face in question.
   * \param [out] cellIDOne the ID of the first cell.
   * \param [out] cellIDTwo the ID of the second cell.
   *
   * \note A face can be associated with one or two cells, depending on whether
   *  it is an external boundary face or interior face. By convention, if a face
   *  is an external boundary face, then only cellIDOne exists and cellIDTwo
   *  will be set to -1.
   *
   * \pre 0 <= faceID < getNumberOfCells()
   */
  virtual void getFaceCellIDs(IndexType faceID,
                              IndexType& cellIDOne,
                              IndexType& cellIDTwo) const final override;

  /// @}

  /// \name Edges
  /// @{

  /*!
   * \brief Return the number of edges in the mesh.
   */
  virtual IndexType getNumberOfEdges() const final override { return m_num_edges; }

  /// @}

  /*!
   * \brief Returns true iff the mesh was constructed with external arrays.
   * \return status true if the mesh points to external buffers, else, false.
   */
  virtual bool isExternal() const override { return false; }

  /// @}

  /// \name Data Accessor Methods
  /// @{

  /// \name Cells
  /// @{

  /*!
   * \brief Returns the cell node offsets array
   * \return offsets pointer to the cell-to-node offsets array.
   *
   * \post offsets != nullptr
   */
  const StackArray<IndexType, 8>& getCellNodeOffsetsArray() const { return m_cell_node_offsets; };

  /*!
   * \brief Copy the node IDs of the given cell into the provided buffer.
   *
   * \param [in] i logical index of the cell along the first dimension.
   * \param [in] j logical index of the cell along the second dimension.
   * \param [in] k logical index of the cell along the third dimension
   *(optional).
   *
   * \param [out] nodes pointer to buffer to populate with the node IDs.
   *
   * \return the number of nodes in the cell.
   */
  /// @{
  inline IndexType getCellNodeIDs(IndexType i, IndexType j, IndexType* nodes) const
  {
    return getCellNodeIDs(i, j, 0, nodes);
  }

  inline IndexType getCellNodeIDs(IndexType i, IndexType j, IndexType k, IndexType* nodes) const;
  /// @}

  /*!
   * \brief Returns the IDs of the faces of the cell at (i,j) or (i, j, k).
   *
   * \param [in] i logical index of the cell along the first dimension.
   * \param [in] j logical index of the cell along the second dimension.
   * \param [in] k logical index of the cell along the third dimension
   *  (optional).
   * \param [out] faces buffer to populate with the face IDs. Must be of length
   *  at least getNumberOfCellFaces().
   *
   * \note The faces are returned in the order of LOWER_I_FACE, UPPER_I_FACE,
   *  LOWER_J_FACE, UPPER_J_FACE and then LOWER_K_FACE, UPPER_K_FACE if 3D.
   *
   * \pre faces != nullptr
   */
  /// @{
  inline void getCellFaceIDs(IndexType i, IndexType j, IndexType* faces) const
  {
    const IndexType cellID = getCellLinearIndex(i, j);
    getCellFaceIDsInternal(cellID, j, faces);
  }

  inline void getCellFaceIDs(IndexType i, IndexType j, IndexType k, IndexType* faces) const
  {
    const IndexType cellID = getCellLinearIndex(i, j, k);
    getCellFaceIDsInternal(cellID, j, k, faces);
  }
  /// @}

  /// @}

  /// \name Faces
  /// @{

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
   * \note Each method is specialized for faces in the I, J, or K direction.
   *
   * \pre nodes != nullptr
   * \pre 0 <= faceID < getNumberOfCells()
   */
  /// @{
  inline IndexType getIFaceNodeIDs(IndexType faceID, IndexType* nodes) const;
  inline IndexType getJFaceNodeIDs(IndexType faceID, IndexType* nodes) const;
  inline IndexType getKFaceNodeIDs(IndexType faceID, IndexType* nodes) const;
  /// @}

  /*!
   * \brief Copy the IDs of the cells adjacent to the given face into the
   *  provided indices.
   *
   * \param [in] faceID the ID of the face in question.
   * \param [out] cellIDOne the ID of the first cell.
   * \param [out] cellIDTwo the ID of the second cell.
   *
   * \note A face can be associated with one or two cells, depending on whether
   *  it is an external boundary face or interior face. By convention, if a face
   *  is an external boundary face, then only cellIDOne exists and cellIDTwo
   *  will be set to -1.
   * 
   * \note Each method is specialized for faces in the I, J, or K direction.
   *
   * \pre 0 <= faceID < getNumberOfCells()
   */
  /// @{
  inline void getIFaceCellIDs(IndexType faceID, IndexType& cellIDOne, IndexType& cellIDTwo) const;
  inline void getJFaceCellIDs(IndexType faceID, IndexType& cellIDOne, IndexType& cellIDTwo) const;
  inline void getKFaceCellIDs(IndexType faceID, IndexType& cellIDOne, IndexType& cellIDTwo) const;
  /// @}

  /// @}

  /// @}

  /// \name Attribute Querying Methods
  /// @{

  /*!
   * \brief Returns the number of nodes along the given dimension.
   *
   * \param [in] dim the dimension to query.
   *
   * \pre 0 <= dim < 3
   */
  inline IndexType getNodeResolution(IndexType dim) const
  {
    SLIC_ASSERT(0 <= dim && dim < 3);
    return m_node_dims[dim];
  }

  /*!
   * \brief Get the global node extent of the mesh.
   *
   * \param [out] extent the buffer to copy the global node extent into, must be
   *  of length at least 6.
   */
  inline void getExtent(int64 extent[6]) const
  {
    SLIC_ASSERT(extent != nullptr);
    std::memcpy(extent, m_node_extent, sizeof(m_node_extent));
  }

  /*!
   * \brief Set the global node extent of the mesh.
   *
   * \param [in] ndims the number of dimensions to set.
   * \param [in] extent the values to set, of length at least 2 * ndims.
   *
   * \pre 0 <= dim < 3
   */
  void setExtent(int ndims, const int64* extent);

  /*!
   * \brief Returns the number of cells along the given dimension.
   *
   * \param [in] dim the dimension to query.
   *
   * \pre 0 <= dim < 3
   */
  inline IndexType getCellResolution(IndexType dim) const
  {
    SLIC_ASSERT(0 <= dim && dim < 3);
    return m_cell_dims[dim];
  }

  /*!
   * \brief Get the total number of I, J, or K faces.
   */
  inline IndexType getTotalNumFaces(int direction) const
  {
    SLIC_ASSERT(I_DIRECTION <= direction);
    SLIC_ASSERT(direction <= K_DIRECTION);
    return m_total_faces[direction];
  }

  /// @}

  /// \name Indexing Helper Methods
  /// @{

  /// \name Nodes
  /// @{

  /*!
   * \brief Returns stride to the second dimension of nodes.
   *
   * \note If the mesh is 1D the returned value is the largest number
   *  representable as an IndexType.
   */
  inline IndexType nodeJp() const { return m_node_jp; }

  /*!
   * \brief kp stride to the third dimension of nodes.
   *
   * \note If the mesh is 1D or 2D the returned value is the largest number
   *  representable as an IndexType.
   */
  inline IndexType nodeKp() const { return m_node_kp; }

  /*!
   * \brief Returns the linear index corresponding to the given logical node
   *  indices.
   *
   * \param [in] i logical node index of the first dimension.
   * \param [in] j logical node index of the second dimension.
   * \param [in] k logical node index of the third dimension (optional)
   *
   * \post 0 <= i < getNodeResolution( I_DIRECTION )
   * \post 0 <= j < getNodeResolution( J_DIRECTION )
   * \post 0 <= k < getNodeResolution( K_DIRECTION )
   */
  inline IndexType getNodeLinearIndex(IndexType i, IndexType j, IndexType k = 0) const
  {
    return i + j * nodeJp() + k * nodeKp();
  }

  /*!
   * \brief Given the 1D linear index of a node, this method computes the
   *  corresponding i-j-k grid index.
   *
   * \param [in] nodeID the local flat index.
   * \param [out] i the corresponding grid index along the I_DIRECTION.
   * \param [out] j the corresponding grid index along the J_DIRECTION.
   * \param [out] k the corresponding grid index along the K_DIRECTION
   *(optional).
   *
   * \note The first method is not valid for 3D meshes.
   *
   * \pre nodeID >= 0 && nodeID < getNumNodes()
   * \post 0 <= i < getNodeResolution( I_DIRECTION )
   * \post 0 <= j < getNodeResolution( J_DIRECTION )
   * \post 0 <= k < getNodeResolution( K_DIRECTION )
   */
  /// @{
  inline void getNodeGridIndex(IndexType nodeID, IndexType& i, IndexType& j) const
  {
    SLIC_ASSERT(m_ndims <= 2);
    j = nodeID / nodeJp();
    i = nodeID - j * nodeJp();
  }

  inline void getNodeGridIndex(IndexType nodeID, IndexType& i, IndexType& j, IndexType& k) const
  {
    k = nodeID / nodeKp();
    const IndexType temp = nodeID - k * nodeKp();
    j = temp / nodeJp();
    i = temp - j * nodeJp();
  }
  /// @}

  /// @}

  /// \name Cells
  /// @{

  /*!
   * \brief Returns stride to the second dimension of cells.
   * \note If the mesh is 1D the returned value is the largest number
   *  representable as an IndexType.
   * \post jp >= 0.
   */
  inline IndexType cellJp() const { return m_cell_jp; }

  /*!
   * \brief Returns stride to the third dimension of cells.
   * \note If the mesh is 1D or 2D the returned value is the largest number
   *  representable as an IndexType.
   * \post kp >= 0.
   */
  inline IndexType cellKp() const { return m_cell_kp; }

  /*!
   * \brief Returns the linear index corresponding to the given logical grid
   * cell indices.
   * \param [in] i logical cell index of the first dimension.
   * \param [in] j logical cell index of the second dimension.
   * \param [in] k logical cell index of the third dimension (optional)
   *
   * \pre 0 <= i < getCellResolution( I_DIRECTION )
   * \pre 0 <= j < getCellResolution( J_DIRECTION )
   * \pre 0 <= k < getCellResolution( K_DIRECTION )
   */
  inline IndexType getCellLinearIndex(IndexType i, IndexType j, IndexType k = 0) const
  {
    return i + j * cellJp() + k * cellKp();
  }

  /*!
   * \brief Given the 1D linear index of a cell, this method computes
   *  the corresponding i-j-k grid index.
   *
   * \param [in] cellID the ID of the cell in question.
   * \param [out] i the corresponding grid index along the I_DIRECTION.
   * \param [out] j the corresponding grid index along the J_DIRECTION.
   * \param [out] k the corresponding grid index along the K_DIRECTION.
   *
   * \note The first method is not valid for 3D meshes.
   *
   * \pre cellID >= 0 && cellID < getNumNodes()
   * \post 0 <= i < getCellResolution( I_DIRECTION )
   * \post 0 <= j < getCellResolution( J_DIRECTION )
   * \post 0 <= k < getCellResolution( K_DIRECTION )
   */
  /// @{
  inline void getCellGridIndex(IndexType cellID, IndexType& i, IndexType& j) const
  {
    SLIC_ASSERT(m_ndims <= 2);
    i = cellID % cellJp();
    j = cellID / cellJp();
  }

  inline void getCellGridIndex(IndexType cellID, IndexType& i, IndexType& j, IndexType& k) const
  {
    k = cellID / cellKp();
    const IndexType rest = cellID - k * cellKp();
    j = rest / cellJp();
    i = rest - j * cellJp();
  }
  /// @}

  /// @}

  /// \name Faces
  /// @{

  /*!
   * \brief Returns the linear index corresponding to the given logical grid
   * face indices.
   *
   * \param [in] dir the direction of the face in question.
   * \param [in] i logical face index of the first dimension.
   * \param [in] j logical face index of the second dimension.
   * \param [in] k logical face index of the third dimension (optional)
   */
  inline IndexType getFaceLinearIndex(int dir, IndexType i, IndexType j, IndexType k = 0) const;

  /*!
   * \brief Returns the linear index corresponding to the given logical face
   *  indices.
   *
   * \param [in] i logical face index of the first dimension.
   * \param [in] j logical face index of the second dimension.
   * \param [in] k logical face index of the third dimension (optional)
   *
   * \note Each method is valid only for the appropriate dimension of the mesh.
   * \note Each method is specialized for a particular direction.
   */
  /// @{
  inline IndexType getIFaceLinearIndex(IndexType i, IndexType j, IndexType k = 0) const
  {
    SLIC_ASSERT(m_ndims >= 2);
    return i + j * getNodeResolution(0) + k * m_num_I_faces_in_k_slice;
  }

  inline IndexType getJFaceLinearIndex(IndexType i, IndexType j, IndexType k = 0) const
  {
    SLIC_ASSERT(m_ndims >= 2);
    return getTotalNumFaces(0) + i + j * getCellResolution(0) + k * m_num_J_faces_in_k_slice;
  }

  inline IndexType getKFaceLinearIndex(IndexType i, IndexType j, IndexType k = 0) const
  {
    SLIC_ASSERT(m_ndims == 3);
    return m_total_IJ_faces + i + j * getCellResolution(0) + k * cellKp();
  }
  /// @}

  /*!
   * \brief Returns the grid indices corresponding to the given linear face
   *  index.
   *
   * \param [in] faceID the ID of the face in question.
   * \param [out] i logical face index of the first dimension.
   * \param [out] j logical face index of the second dimension.
   * \param [out] k logical face index of the third dimension (optional)
   *
   * \note Each method is valid only for the appropriate dimension of the mesh.
   * \note Each method is specialized for a particular direction.
   */
  /// @{
  inline void getIFaceGridIndex(IndexType faceID, IndexType& i, IndexType& j) const
  {
    SLIC_ASSERT(m_ndims == 2);
    SLIC_ASSERT(0 <= faceID && faceID < getTotalNumFaces(0));

    j = faceID / getNodeResolution(0);
    i = faceID - getNodeResolution(0) * j;
  }

  inline void getIFaceGridIndex(IndexType faceID, IndexType& i, IndexType& j, IndexType& k) const
  {
    SLIC_ASSERT(m_ndims >= 2);
    SLIC_ASSERT(0 <= faceID && faceID < getTotalNumFaces(0));

    k = getIFaceKIndex(faceID);
    const IndexType rest = m_num_I_faces_in_k_slice * k;
    j = (faceID - rest) / getNodeResolution(0);
    i = faceID - getNodeResolution(0) * j - rest;
  }

  inline void getJFaceGridIndex(IndexType faceID, IndexType& i, IndexType& j) const
  {
    SLIC_ASSERT(m_ndims == 2);
    const IndexType shiftedID = shiftJFaceID(faceID);
    SLIC_ASSERT(0 <= shiftedID && shiftedID < getTotalNumFaces(1));

    j = shiftedID / getCellResolution(0);
    i = shiftedID - getCellResolution(0) * j;
  }

  inline void getJFaceGridIndex(IndexType faceID, IndexType& i, IndexType& j, IndexType& k) const
  {
    SLIC_ASSERT(m_ndims >= 2);
    const IndexType shiftedID = shiftJFaceID(faceID);
    SLIC_ASSERT(0 <= shiftedID && shiftedID < getTotalNumFaces(1));

    k = getJFaceKIndex(shiftedID);
    j = getJFaceJIndex(shiftedID, k);
    i = shiftedID - j * getCellResolution(0) - k * m_num_J_faces_in_k_slice;
  }

  inline void getKFaceGridIndex(IndexType faceID, IndexType& i, IndexType& j, IndexType& k) const
  {
    SLIC_ASSERT(m_ndims == 3);
    const IndexType shiftedID = shiftKFaceID(faceID);
    SLIC_ASSERT(0 <= shiftedID && shiftedID < getTotalNumFaces(2));

    k = getKFaceKIndex(shiftedID);
    j = getKFaceJIndex(shiftedID, k);
    i = shiftedID - cellJp() * j - cellKp() * k;
  }
  /// @}

  /// @}

  /// @}

protected:
  /*!
   * \brief Constructs a structured mesh instance from the given extent.
   *
   * \param [in] meshType the mesh type
   * \param [in] Ni the number of nodes along the I direction.
   * \param [in] Nj the number of nodes along the J direction.
   * \param [in] Nk the number of nodes along the K direction.
   *
   * \note Nj and or Nk may be -1 to signify a 1D or 2D mesh.
   */
  StructuredMesh(int meshType, IndexType Ni, IndexType Nj, IndexType Nk);

#ifdef AXOM_MINT_USE_SIDRE

  /*!
   * \brief Constructs a structured mesh instance from the given group.
   *
   * \param [in] group pointer to the sidre::Group
   * \param [in] topo the topology name to use, an empty string may be supplied.
   *
   * \note If an empty string is supplied for the topology name, the code will
   *  use the 1st topology group under the parent "topologies" group.
   *
   * \pre group !=  nullptr
   * \pre blueprint::isValidRootGroup( group ) == true
   *
   * \note This constructor forwards this call to the parent Mesh class.
   *
   * \see Mesh( sidre::Group* group, const std::string& topo )
   */
  StructuredMesh(sidre::Group* group, const std::string& topo);

  /*!
   * \brief Constructs a structured mesh instance on the specified group.
   *
   * \param [in] meshType the mesh type.
   * \param [in] Ni the number of nodes along the I direction.
   * \param [in] Nj the number of nodes along the J direction.
   * \param [in] Nk the number of nodes along the K direction.
   * \param [in] group pointer to the group in the Sidre hierarchy.
   * \param [in] topo the topology name to use, may be an empty string.
   * \param [in] coordset the coordset name to use, may be an empty string.
   *
   * \note If an empty string is supplied for the topology and coordset name
   *  respectively, an internal default name will be provided by the
   *  implementation.
   *
   * \note Nj and or Nk may be -1 to signify a 1D or 2D mesh.
   *
   * \pre 1 <= dimension <= 3
   * \pre group != nullptr
   * \pre group->getNumGroups() == 0
   * \pre group->getNumViews() == 0
   *
   * \post blueprint::isValidRootGroup( group )
   *
   * \see Mesh( int ndims, int type, sidre::Group*,
   *            const std::string& topo, const std::string& coordset );
   */
  StructuredMesh(int meshType,
                 IndexType Ni,
                 IndexType Nj,
                 IndexType Nk,
                 sidre::Group* group,
                 const std::string& topo,
                 const std::string& coordset);

#endif

  /*!
   * \brief Initialize all the StructuredMesh members.
   */
  void structuredInit();

  StackArray<IndexType, 3> m_node_dims = {{0, 0, 0}};
  StackArray<int64, 6> m_node_extent = {{0, 0, 0, 0, 0, 0}};

  IndexType m_node_jp = 0;
  IndexType m_node_kp = 0;

  StackArray<IndexType, 3> m_cell_dims = {{0, 0, 0}};
  IndexType m_cell_jp = 0;
  IndexType m_cell_kp = 0;
  StackArray<IndexType, 8> m_cell_node_offsets = {{0, 0, 0, 0, 0, 0, 0, 0}};

  StackArray<IndexType, 3> m_total_faces = {{0, 0, 0}};
  IndexType m_total_IJ_faces = 0;
  IndexType m_num_I_faces_in_k_slice = 0;
  IndexType m_num_J_faces_in_k_slice = 0;

  IndexType m_num_edges = 0;

private:
  DISABLE_COPY_AND_ASSIGNMENT(StructuredMesh);
  DISABLE_MOVE_AND_ASSIGNMENT(StructuredMesh);

  /*!
   * \brief Copy the face IDs of the given cell into the provided buffer.
   *  The buffer must be of length at least getNumberOfCellFaces().
   *
   * \param [in] cellID the ID of the cell in question.
   * \param [in] j the grid j index of the cell in question.
   * \param [in] k the grid k index of the cell in question (optional).
   * \param [out] faces the buffer into which the face IDs are copied, must
   *  be of length at least getNumberOfCellFaces().
   *
   * \return The number of faces for the given cell.
   *
   * \note The faces are returned in the order of LOWER_I_FACE, UPPER_I_FACE,
   *  LOWER_J_FACE, UPPER_J_FACE and then LOWER_K_FACE, UPPER_K_FACE if 3D.
   *
   * \pre faces != nullptr
   * \pre 0 <= cellID < getNumberOfCells()
   */
  /// @{
  inline IndexType getCellFaceIDsInternal(IndexType cellID, IndexType j, IndexType* faces) const;

  inline IndexType getCellFaceIDsInternal(IndexType cellID,
                                          IndexType j,
                                          IndexType k,
                                          IndexType* faces) const;
  /// @}

  /*!
   * \brief Returns the shifted linear index of the given face ID.
   *
   * \param [in] faceID the ID of the face in question.
   */
  /// @{
  inline IndexType shiftJFaceID(IndexType faceID) const { return faceID - getTotalNumFaces(0); }

  inline IndexType shiftKFaceID(IndexType faceID) const
  {
    return faceID - getTotalNumFaces(0) - getTotalNumFaces(1);
  }
  /// @}

  /*!
   * \brief Return the K grid index of the given I direction face.
   *
   * \param [in] faceID the face ID of the I face in question.
   */
  inline IndexType getIFaceKIndex(IndexType faceID) const
  {
    return faceID / m_num_I_faces_in_k_slice;
  }

  /*!
   * \brief Return the K grid index of the given J direction face.
   *
   * \param [in] shiftedID the shifted face ID of the J face in question.
   */
  inline IndexType getJFaceKIndex(IndexType shiftedID) const
  {
    return shiftedID / m_num_J_faces_in_k_slice;
  }

  /*!
   * \brief Return the J grid index of the given J direction face.
   *
   * \param [in] shiftedID the shifted face ID of the J face in question.
   * \param [in] k the K grid index of the face.
   */
  inline IndexType getJFaceJIndex(IndexType shiftedID, IndexType k) const
  {
    return (shiftedID - m_num_J_faces_in_k_slice * k) / getCellResolution(0);
  }

  /*!
   * \brief Return the K grid index of the given K direction face.
   *
   * \param [in] shiftedID the shifted face ID of the K face in question.
   */
  inline IndexType getKFaceKIndex(IndexType shiftedID) const { return shiftedID / cellKp(); }

  /*!
   * \brief Return the J grid index of the given K direction face.
   *
   * \param [in] shiftedID the shifted face ID of the K face in question.
   * \param [in] k the K grid index of the face.
   */
  inline IndexType getKFaceJIndex(IndexType shiftedID, IndexType k) const
  {
    return (shiftedID - cellKp() * k) / cellJp();
  }
};

//------------------------------------------------------------------------------
//      In-lined Method Implementations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline IndexType StructuredMesh::getCellNodeIDs(IndexType cellID, IndexType* nodes) const
{
  SLIC_ASSERT(nodes != nullptr);
  SLIC_ASSERT(0 <= cellID && cellID < getNumberOfCells());

  // Calculate logical indices of the cell's first corner node.
  IndexType i, j, k;
  getCellGridIndex(cellID, i, j, k);
  return getCellNodeIDs(i, j, k, nodes);
}

//------------------------------------------------------------------------------
inline IndexType StructuredMesh::getCellFaceIDs(IndexType cellID, IndexType* faces) const
{
  SLIC_ASSERT(faces != nullptr);
  SLIC_ASSERT(0 <= cellID && cellID < getNumberOfCells());

  if(m_ndims == 2)
  {
    IndexType j = cellID / cellJp();
    return getCellFaceIDsInternal(cellID, j, faces);
  }
  else if(m_ndims == 3)
  {
    IndexType k = cellID / cellKp();
    IndexType j = (cellID - k * cellKp()) / cellJp();
    return getCellFaceIDsInternal(cellID, j, k, faces);
  }
  else
  {
    return 0;
  }
}

//------------------------------------------------------------------------------
inline IndexType StructuredMesh::getFaceNodeIDs(IndexType faceID, IndexType* nodes) const
{
  SLIC_ASSERT(0 <= faceID && faceID < getNumberOfFaces());

  /* Check which kind of face is being asked for */
  if(faceID < getTotalNumFaces(0))
  {
    /* It's a I_DIRECTION face */
    return getIFaceNodeIDs(faceID, nodes);
  }
  else if(faceID < m_total_IJ_faces)
  {
    /* It's a J_DIRECTION face */
    return getJFaceNodeIDs(faceID, nodes);
  }
  else
  {
    /* It's a K_DIRECTION face */
    return getKFaceNodeIDs(faceID, nodes);
  }
}

//------------------------------------------------------------------------------
inline void StructuredMesh::getFaceCellIDs(IndexType faceID,
                                           IndexType& cellIDOne,
                                           IndexType& cellIDTwo) const
{
  /* Check which kind of face is being asked for */
  if(faceID < getTotalNumFaces(0))
  {
    /* It's a I_DIRECTION face */
    return getIFaceCellIDs(faceID, cellIDOne, cellIDTwo);
  }
  else if(faceID < m_total_IJ_faces)
  {
    /* It's a J_DIRECTION face */
    return getJFaceCellIDs(faceID, cellIDOne, cellIDTwo);
  }
  else
  {
    /* It's a K_DIRECTION face */
    return getKFaceCellIDs(faceID, cellIDOne, cellIDTwo);
  }
}

//------------------------------------------------------------------------------
inline IndexType StructuredMesh::getCellNodeIDs(IndexType i,
                                                IndexType j,
                                                IndexType k,
                                                IndexType* nodes) const
{
  SLIC_ASSERT(nodes != nullptr);
  SLIC_ASSERT(0 <= i && i < getNodeResolution(0));
  SLIC_ASSERT(m_ndims < 2 || (0 <= j && j < getNodeResolution(1)));
  SLIC_ASSERT(m_ndims < 3 || (0 <= k && k < getNodeResolution(2)));

  /* Get the node index of the first node in the cell. */
  const IndexType n0 = getNodeLinearIndex(i, j, k);

  /* Fill in the node connectivity with the known offsets */
  const IndexType n_cell_nodes = getNumberOfCellNodes();
  for(IndexType ii = 0; ii < n_cell_nodes; ++ii)
  {
    nodes[ii] = n0 + m_cell_node_offsets[ii];
  }

  return n_cell_nodes;
}

//------------------------------------------------------------------------------
inline IndexType StructuredMesh::getIFaceNodeIDs(IndexType faceID, IndexType* nodes) const
{
  SLIC_ASSERT(0 <= faceID && faceID < getTotalNumFaces(0));

  if(m_ndims == 2)
  {
    nodes[0] = faceID;
    nodes[1] = nodes[0] + m_cell_node_offsets[3];
    return 2;
  }
  else if(m_ndims == 3)
  {
    SLIC_ASSERT(m_ndims == 3);
    const IndexType k = getIFaceKIndex(faceID);
    nodes[0] = faceID + getNodeResolution(0) * k;
    nodes[1] = nodes[0] + m_cell_node_offsets[4];
    nodes[2] = nodes[0] + m_cell_node_offsets[7];
    nodes[3] = nodes[0] + m_cell_node_offsets[3];
    return 4;
  }
  else
  {
    return 0;
  }
}

//------------------------------------------------------------------------------
inline IndexType StructuredMesh::getJFaceNodeIDs(IndexType faceID, IndexType* nodes) const
{
  const IndexType shiftedID = shiftJFaceID(faceID);
  SLIC_ASSERT(0 <= faceID && getTotalNumFaces(1));

  if(m_ndims == 2)
  {
    const IndexType j = getJFaceJIndex(shiftedID, 0);
    nodes[0] = shiftedID + j;
    nodes[1] = nodes[0] + 1;
    return 2;
  }
  else
  {
    SLIC_ASSERT(m_ndims == 3);
    const IndexType k = getJFaceKIndex(shiftedID);
    const IndexType j = getJFaceJIndex(shiftedID, k);
    nodes[0] = shiftedID + j + getNodeResolution(1) * k;
    nodes[1] = nodes[0] + 1;
    nodes[2] = nodes[0] + m_cell_node_offsets[5];
    nodes[3] = nodes[0] + m_cell_node_offsets[4];
    return 4;
  }
}

//------------------------------------------------------------------------------
inline IndexType StructuredMesh::getKFaceNodeIDs(IndexType faceID, IndexType* nodes) const
{
  SLIC_ASSERT(m_ndims == 3);
  const IndexType shiftedID = shiftKFaceID(faceID);
  SLIC_ASSERT(shiftedID >= 0 && shiftedID < getTotalNumFaces(2));

  const IndexType k = getKFaceKIndex(shiftedID);
  const IndexType j = getKFaceJIndex(shiftedID, k);
  nodes[0] = shiftedID + j + (getCellResolution(0) + getCellResolution(1) + 1) * k;
  nodes[1] = nodes[0] + 1;
  nodes[2] = nodes[0] + m_cell_node_offsets[2];
  nodes[3] = nodes[0] + m_cell_node_offsets[3];
  return 4;
}

//------------------------------------------------------------------------------
inline void StructuredMesh::getIFaceCellIDs(IndexType faceID,
                                            IndexType& cellIDOne,
                                            IndexType& cellIDTwo) const
{
  SLIC_ASSERT(m_ndims == 2 || m_ndims == 3);

  IndexType i, j, k;
  getIFaceGridIndex(faceID, i, j, k);

  cellIDOne = getCellLinearIndex(i - 1, j, k);
  cellIDTwo = getCellLinearIndex(i, j, k);
  if(i == 0)
  {
    cellIDOne = cellIDTwo;
    cellIDTwo = -1;
  }
  else if(i == getCellResolution(0))
  {
    cellIDTwo = -1;
  }
}

//-----------------------------------------------------------------------------
inline void StructuredMesh::getJFaceCellIDs(IndexType faceID,
                                            IndexType& cellIDOne,
                                            IndexType& cellIDTwo) const
{
  SLIC_ASSERT(m_ndims == 2 || m_ndims == 3);

  IndexType i, j, k;
  getJFaceGridIndex(faceID, i, j, k);

  cellIDOne = getCellLinearIndex(i, j - 1, k);
  cellIDTwo = getCellLinearIndex(i, j, k);
  if(j == 0)
  {
    cellIDOne = cellIDTwo;
    cellIDTwo = -1;
  }
  else if(j == getCellResolution(1))
  {
    cellIDTwo = -1;
  }
}

//------------------------------------------------------------------------------
inline void StructuredMesh::getKFaceCellIDs(IndexType faceID,
                                            IndexType& cellIDOne,
                                            IndexType& cellIDTwo) const
{
  SLIC_ASSERT(m_ndims == 3);

  IndexType i, j, k;
  getKFaceGridIndex(faceID, i, j, k);

  cellIDOne = getCellLinearIndex(i, j, k - 1);
  cellIDTwo = getCellLinearIndex(i, j, k);
  if(k == 0)
  {
    cellIDOne = cellIDTwo;
    cellIDTwo = -1;
  }
  else if(k == getCellResolution(2))
  {
    cellIDTwo = -1;
  }
}

//------------------------------------------------------------------------------
inline IndexType StructuredMesh::getFaceLinearIndex(int dir, IndexType i, IndexType j, IndexType k) const
{
  if(dir == I_DIRECTION)
  {
    return getIFaceLinearIndex(i, j, k);
  }
  else if(dir == J_DIRECTION)
  {
    return getJFaceLinearIndex(i, j, k);
  }
  else
  {
    SLIC_ASSERT(dir == K_DIRECTION);
    SLIC_ASSERT(m_ndims == 3);
    return getKFaceLinearIndex(i, j, k);
  }
}

//------------------------------------------------------------------------------
inline IndexType StructuredMesh::getCellFaceIDsInternal(IndexType cellID,
                                                        IndexType j,
                                                        IndexType* faces) const
{
  SLIC_ASSERT(m_ndims == 2);
  SLIC_ASSERT(0 <= cellID && cellID < getNumberOfCells());
  SLIC_ASSERT(0 <= j && j < getCellResolution(1));
  SLIC_ASSERT(faces != nullptr);

  /* The I_DIRECTION faces */
  faces[0] = cellID + j;
  faces[1] = faces[0] + 1;

  /* The J_DIRECTION faces */
  faces[2] = cellID + getTotalNumFaces(0);
  faces[3] = faces[2] + getCellResolution(0);

  return 4;
}

//------------------------------------------------------------------------------
inline IndexType StructuredMesh::getCellFaceIDsInternal(IndexType cellID,
                                                        IndexType j,
                                                        IndexType k,
                                                        IndexType* faces) const
{
  SLIC_ASSERT(m_ndims == 3);
  SLIC_ASSERT(0 <= cellID && cellID < getNumberOfCells());
  SLIC_ASSERT(0 <= j && j < getCellResolution(1));
  SLIC_ASSERT(0 <= k && k < getCellResolution(2));
  SLIC_ASSERT(faces != nullptr);

  /* The I_DIRECTION faces */
  faces[0] = cellID + j + getCellResolution(1) * k;
  faces[1] = faces[0] + 1;

  /* The J_DIRECTION faces */
  faces[2] = cellID + getTotalNumFaces(0) + getCellResolution(0) * k;
  faces[3] = faces[2] + getCellResolution(0);

  /* The K_DIRECTION faces */
  faces[4] = cellID + m_total_IJ_faces;
  faces[5] = faces[4] + cellKp();

  return 6;
}

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_STRUCTUREDMESH_HPP_ */
