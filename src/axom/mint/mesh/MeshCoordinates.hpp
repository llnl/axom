// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_MESHCOORDINATES_HPP_
#define MINT_MESHCOORDINATES_HPP_

#include "axom/core/Macros.hpp"  // for Axom macros and definitions
#include "axom/mint/deprecated/MCArray.hpp"

#ifdef AXOM_MINT_USE_SIDRE
  #include "axom/mint/deprecated/SidreMCArray.hpp"  // for sidre::deprecated::MCArray
#endif

#include "axom/slic/interface/slic.hpp"  // for slic logging macros

#include "axom/mint/config.hpp"  // for axom::IndexType

namespace axom
{
#ifdef AXOM_MINT_USE_SIDRE

// Sidre forward declarations
namespace sidre
{
class Group;
}

#endif

namespace mint
{
static constexpr int X_COORDINATE = 0;
static constexpr int Y_COORDINATE = 1;
static constexpr int Z_COORDINATE = 2;

/*!
 * \class MeshCoordinates
 *
 * \brief Provides functionality to store and operate on mesh node coordinates.
 *
 *  A MeshCoordinates object is generally tied to an associated Mesh instance
 *  and provides the means to store, access and modify the node coordinates of
 *  a mesh.
 *
 *  The MeshCoordinates object may be constructed using (a) native storage,
 *  (b) external storage, or, (c) Sidre:
 *
 *  * <b> Native Storage </b> <br />
 *
 *    When using native storage, the MeshCoordinates object owns all associated
 *    memory. The storage can dynamically grow as needed, e.g., when adding
 *    more nodes. Typically, extra space is allocated to minimize the number
 *    of re-allocations. At any given instance, the total node capacity
 *    can be queried by calling the capacity() function. The extra memory can
 *    be returned to the system by calling the shrink() method.
 *
 *    When all extra memory is exhausted, appending a new node triggers a
 *    re-allocation. The amount of extra space that is allocated is according
 *    to the <em> resize_ratio </em> parameter, which is set to 2.0 by default.
 *    The <em> resize_ratio </em> may be queried and set to a different value
 *    by the getResizeRatio() and setResizeRatio() functions respectively.
 *
 *    When the MeshCoordinates object goes out-of-scope, all memory associated
 *    with the given instance is returned to the system.
 *
 *  * <b> External Storage </b> <br />
 *
 *    A MeshCoordinates object may also be constructed from external,
 *    user-supplied buffers that store the mesh node coordinates. In this case,
 *    the memory is owned by the caller. The MeshCoordinates object just keeps
 *    pointers to the user-supplied buffers.
 *
 *    \warning Since the memory is not owned by the MeshCoordinates object
 *     when external buffers are supplied, the MeshCoordinates object cannot
 *     dynamically grow the storage. Consequently, the number of nodes the
 *     MeshCoordinates instance can hold is fixed. All calls to `shrink()` and
 *     `reserve()` will fail.
 *
 *    \warning Moreover, when the MeshCoordinates object goes out-of-scope, the
 *     associated buffers are not deleted. The caller owns the external data
 *     and has the responsibility of properly de-allocating the associated
 *     memory.
 *
 *  * <b> Sidre </b> <br />
 *
 *    A MeshCoordinates object may also be constructed from a sidre::Group which
 *    conforms to the <a href="http://llnl-conduit.readthedocs.io/en/latest/">
 *    mesh blueprint </a>.
 *
 *    A MeshCoordinates object that is bound to a particular sidre::Group
 *    supports all operations transparently including dynamically growing the
 *    storage to hold more nodes as needed, but, instead, Sidre owns the memory.
 *    All memory management operations are delegated to Sidre.
 *
 *    \warning Once the MeshCoordinates object goes out-of-scope, the data
 *     stays remains persistent in Sidre.
 *
 * \warning Reallocations tend to be costly operations in terms of performance.
 *  Use `reserve()` when the number of nodes is known a priori, or opt to
 *  use a constructor that takes an actual size and capacity when possible.
 *
 * \see axom::deprecated::MCArray
 * \see mint::Mesh
 * \see sidre::Group
 * \see sidre::deprecated::MCArray
 */
class MeshCoordinates
{
public:
  /// \name Native Storage Constructors
  /// @{

  /*!
   * \brief Creates a MeshCoordinates instance of specified dimension and
   *  given number of nodes. Optionally, an initial max capacity may be
   *  specified in the constructor to reserve additional space for storing
   *  nodes. If a capacity is not specified, an internal initial default
   *  capacity will be computed instead.
   *
   * \param [in] dimension the mesh dimension.
   * \param [in] numNodes the number of nodes this instance will hold initially
   * \param [in] capacity initial max capacity to reserve space for (optional).
   *
   * \pre 1 <= dimension() <= 3
   * \post numNodes() <= capacity()
   * \post numNodes() == numNodes
   * \post if capacity == USE_DEFAULT then
   *  capacity() == max(DEFAULT_CAPACITY, numNodes()*DEFAULT_RESIZE_RATIO)
   */
  explicit MeshCoordinates(int dimension, IndexType numNodes = 0, IndexType capacity = USE_DEFAULT);

  /// @}

  /// \name External Storage Constructors
  /// @{

  /*!
   * \brief Creates a MeshCoordinates object from the given coordinate buffers.
   *
   * \param [in] numNodes the number of nodes in the supplied array
   * \param [in] capacity the actual node capacity on the supplied buffers.
   * \param [in] x pointer to the x-coordinates
   * \param [in] y pointer to the y-coordinates, may be nullptr if 1D
   * \param [in] z pointer to the z-coordinates, may be nullptr if 2D or 1D
   *
   * \warning This constructor wraps the given coordinate buffers and does not
   *  own the data. Consequently, new nodes cannot be added to this
   *  MeshCoordinates instance when using this constructor.
   *
   * \warning All calls to shrink(), append(), resize() and reserve() will fail.
   *
   * \pre x != nullptr
   * \pre y != nullptr if dimension==2 || dimension==3
   * \pre z != nullptr if dimension==3
   * \pre numNodes >= 1
   *
   * \post 1 <= dimension() <= 3
   * \post empty() == false
   * \post numNodes() == numNodes
   * \post numNodes() <= capacity()
   */
  /// @{
  MeshCoordinates(IndexType numNodes,
                  IndexType capacity,
                  double* x,
                  double* y = nullptr,
                  double* z = nullptr);

  MeshCoordinates(IndexType numNodes, double* x, double* y = nullptr, double* z = nullptr);
  /// @}

  /// @}

  /// \name Sidre Constructors
  /// @{

#ifdef AXOM_MINT_USE_SIDRE

  /*!
   * \brief Creates a MeshCoordinates object from the given sidre::Group that
   *  holds mesh coordinate data according to the computational mesh blueprint
   *  conventions.
   *
   * \param [in] group the sidre::Group instance that holds mesh coordinates.
   *
   * \note The supplied sidre::Group must conform to the mesh blueprint
   *  computational mesh protocol for storing coordinates. For details on
   *  the mesh blueprint, see http://llnl-conduit.readthedocs.io/en/latest/.
   *
   * \pre group != nullptr
   *
   * \see sidre::Group
   */
  MeshCoordinates(sidre::Group* group);

  /*!
   * \brief Creates a MeshCoordinates object on an empty sidre::Group
   *
   * \param [in,out] group the sidre::Group to hold the mesh coordinates
   * \param [in] dimension the mesh dimension
   * \param [in] numNodes the number of nodes initially this instance will hold.
   * \param [in] capacity the initial max node capacity (optional)
   *
   * \note This constructor will create the appropriate groups and views under
   *  the given sidre::Group.
   *
   * \warning The user-supplied sidre::Group is expected to be empty
   *
   * \pre group != nullptr.
   * \pre group->getNumGroups()==0
   * \pre group->getNumViews()==0
   *
   * \post group->getNumGroups()==2
   * \post group->getNumViews()==dimension
   * \post numNodes() == numNodes
   * \post numNodes() <= capacity()
   *
   * \see sidre::Group
   */
  MeshCoordinates(sidre::Group* group,
                  int dimension,
                  IndexType numNodes = 0,
                  IndexType capacity = USE_DEFAULT);

#endif

  /// @}

  /*!
   * \brief Destructor, free's the allocated vectors.
   */
  ~MeshCoordinates();

  /// \name Attribute get/set Methods
  /// @{

  /*!
   * \brief Returns the dimension of this MeshCoordinates instance.
   * \return dim the dimension
   * \post 1 <= dim <= 3
   */
  int dimension() const { return m_ndims; };

  /*!
   * \brief Returns the number of nodes in this MeshCoordinates instance.
   * \return N the number of nodes in this MeshCoordinates instance.
   */
  IndexType numNodes() const { return m_coordinates[0]->size(); }

  /*!
   * \brief Changes the number of nodes to the specified size.
   *
   * \post numNodes()==size
   * \post numNodes() <= capacity()
   */
  void resize(IndexType size);

  /*!
   * \brief Get the maximum number of points that can currently be held.
   * \return N the capacity of m_coordinates.
   * \post N >= numNodes()
   */
  IndexType capacity() const { return m_coordinates[0]->capacity(); }

  /*!
   * \brief Changes the max capacity to the specified capacity.
   *
   * \warning The reserve() operation is invalid when a MeshCoordinates object
   *  is constructed using external buffers. Since, in this case the object does
   *  not own the associated memory, dynamic resizing is not allowed.
   *
   * \param [in] capacity the new max node capacity.
   *
   * \post capacity() == capacity
   * \post numNodes() <= capacity()
   */
  void reserve(IndexType capacity);

  /*!
   * \brief Returns all extra memory to the system.
   *
   * \warning The shrink() operation is invalid when a MeshCoordinates object
   *  is constructed using external buffers. Since, in this case the object does
   *  not own the associated memory, dynamic resizing is not allowed.
   *
   * \post numNodes() == capacity()
   */
  void shrink();

  /*!
   * \brief Returns the resize ratio by which the capacity will increase.
   * \return N the ratio by which the capacity will increase
   * \see Array::getResizeRatio()
   */
  double getResizeRatio() const { return m_coordinates[0]->getResizeRatio(); }

  /*!
   * \brief Sets the resize ratio by which the capacity will increase when
   *  a dynamically re-allocation is triggered.
   *
   * \param [in] ratio the ratio by which the capacity will increase
   * \see Array::setResizeRatio()
   */
  void setResizeRatio(double ratio);

  /*
   * \brief Checks if this MeshCoordinates instance is empty.
   * \return status true if numNodes()==0, else, false.
   */
  bool empty() const { return numNodes() == 0; }

  /*!
   * \brief Return true iff constructed via the external constructor.
   */
  bool isExternal() const;

  /*!
   * \brief Return true iff constructed via the sidre constructors.
   */
  bool isInSidre() const;

#ifdef AXOM_MINT_USE_SIDRE
  /*!
   * \brief Return the sidre group associated with this instance or nullptr
   *  if none exits.
   */
  sidre::Group* getSidreGroup() { return m_group; }
#endif

  /// @}

  /// \name Data Access Methods
  /// @{

  /*!
   * \brief Appends a new node to the MeshCoordinates instance
   *
   * \param [in] x the first coordinate to append.
   * \param [in] y the second coordinate to append.
   * \param [in] z the third coordinate to append.
   *
   * \note Each method is valid only for the appropriate dimension of the mesh.
   * \post the number of nodes is incremented by one.
   */
  /// @{

  IndexType append(double x);
  IndexType append(double x, double y);
  IndexType append(double x, double y, double z);

  /// @}

  /*!
   * \brief Appends multiple nodes to the MeshCoordinates instance
   *
   * \param [in] coords pointer to the nodes to append, of length
   *  n * getDimension().
   * \param [in] n the number of nodes to append.
   *
   * \note coords is assumed to be in the array of structs format, ie
   *  coords = {x0, y0, z0, x1, y1, z1, ..., xn, yn, zn}.
   *
   * \pre coords != nullptr
   * \pre n >= 0
   */
  void append(const double* coords, IndexType n = 1);

  /*!
   * \brief Appends new nodes to the MeshCoordinates instance
   *
   * \param [in] x array of the first coordinates to append, of length n.
   * \param [in] y array of the second coordinates to append, of length n.
   * \param [in] z array of the third coordinates to append, of length n.
   * \param [in] n the number of coordinates to append.
   *
   * \note The first method is only valid for 2D meshes while the second
   *  is only for 3D.
   *
   * \pre x != nullptr
   * \pre y != nullptr
   * \pre z != nullptr
   * \pre n >= 0
   * \post the number of nodes is incremented by n
   */
  /// @{

  void append(const double* x, const double* y, IndexType n);
  void append(const double* x, const double* y, const double* z, IndexType n);

  /// @}

  /*!
   * \brief Sets the coordinates for the given node.
   *
   * \param [in] nodeID the index of the node whose coordinates to set.
   * \param [in] x the new value of the first coordinate.
   * \param [in] y the new value of the second coordinate.
   * \param [in] z the new value of the third coordinate.
   *
   * \note Each method is valid only for the appropriate dimension of the mesh.
   * \post idx >=0 && idx == numNodes()-1
   */
  /// @{

  void set(IndexType nodeID, double x);
  void set(IndexType nodeID, double x, double y);
  void set(IndexType nodeID, double x, double y, double z);

  /// @}

  /*!
   * \brief Insert a node to the MeshCoordinates instance.
   *
   * \param [in] nodeID the position to insert at.
   * \param [in] x the value of the first coordinate to insert.
   * \param [in] y the value of the second coordinate to insert.
   * \param [in] z the value of the third coordinate to insert.
   *
   * \note Each method is valid only for the appropriate dimension of the mesh.
   * \pre 0 <= nodeID <= getNumberOfNodes
   * \post the number of nodes is incremented by 1
   */
  /// @{

  void insert(IndexType nodeID, double x);
  void insert(IndexType nodeID, double x, double y);
  void insert(IndexType nodeID, double x, double y, double z);

  /// @}

  /*!
   * \brief Inserts multiple nodes to the MeshCoordinates instance
   *
   * \param [in] coords pointer to the nodes to insert, of length
   *  n * getDimension().
   * \param [in] n the number of nodes to append.
   *
   * \note coords is assumed to be in the array of structs format, ie
   *  coords = {x0, y0, z0, x1, y1, z1, ..., xn, yn, zn}.
   *
   * \pre 0 <= nodeID <= getNumberOfNodes
   * \pre coords != nullptr
   * \pre n >= 0
   * \post the number of nodes is incremented by n
   */
  void insert(IndexType nodeID, const double* coords, IndexType n = 1);

  /*!
   * \brief Insert multiple nodes to the MeshCoordinates instance.
   *
   * \param [in] nodeID the position to insert at.
   * \param [in] x the array of the first coordinates to insert.
   * \param [in] y the array of the second coordinates to insert.
   * \param [in] z the array of the third coordinates to insert.
   * \param [in] n the number of nodes to insert.
   *
   * \pre getDimension() == 3
   * \pre 0 <= nodeID <= getNumberOfNodes
   * x != nullptr
   * y != nullptr
   * z != nullptr
   * \pre n >= 0
   * \post the number of nodes is incremented by n
   */
  /// @{

  void insert(IndexType nodeID, const double* x, const double* y, IndexType n);
  void insert(IndexType nodeID, const double* x, const double* y, const double* z, IndexType n);

  /// @}

  /*!
   * \brief Returns the coordinate of a point at the given dimension.
   *
   * \param [in] nodeID the index of the point in query.
   * \param [in] dim the dimension in query.
   *
   * \return coord the coordinate of the point.
   *
   * \pre dim < m_ndims
   * \pre (nodeID >= 0) && (nodeID < getNumberOfPoints())
   */
  double getCoordinate(IndexType nodeID, int dim) const;

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
  void getCoordinates(IndexType nodeID, double* coords) const;

  /*!
   * \brief Returns pointer to the coordinate array at the requested dimension.
   * \param [in] dim the requested dimension.
   * \return coord_array pointer to the coordinate array
   *
   * \pre dim < m_ndims
   * \post coord_array != nullptr.
   */
  /// @{

  double* getCoordinateArray(int dim)
  {
    SLIC_ERROR_IF(!indexInRange(dim, 0, m_ndims - 1),
                  "invalid request for coordinate array along dimension [" << dim << "]"
                                                                           << "ndims=" << m_ndims);

    SLIC_ASSERT(m_coordinates[dim] != nullptr);
    return m_coordinates[dim]->getData();
  }

  const double* getCoordinateArray(int dim) const
  {
    SLIC_ERROR_IF(!indexInRange(dim, 0, m_ndims - 1),
                  "invalid request for coordinate array along dimension [" << dim << "]"
                                                                           << "ndims=" << m_ndims);

    SLIC_ASSERT(m_coordinates[dim] != nullptr);
    return m_coordinates[dim]->getData();
  }

  /// @}

  /// @}

private:
  /// \name Private helper methods
  /// @{

  /*!
   * \brief Checks if the given index is within the specified range.
   *
   * \param [in] i the index to check
   * \param [in] A lower bound
   * \param [in] B upper bound
   *
   * \return status true \f$ \iff i \in [A,B] \f$, otherwise false
   */
  bool indexInRange(int i, int A, int B) const { return (i >= A && i <= B); }

  /*!
   * \brief Helper method to check to validate the supplied dimension.
   * \note Used primarily to do pre-condition checks
   * \return status true if the dimension is invalid, otherwise, false.
   */
  bool invalidDimension() const { return (m_ndims < 1 || m_ndims > 3); }

  /*!
   * \brief Checks if the given node index is within the valid range.
   * \param [in] idx the node index to check
   * \return status true if the index is valid, false, otherwise.
   */
  bool validIndex(IndexType idx) const
  {
    return indexInRange(static_cast<int>(idx), 0, numNodes() - 1);
  }

  /*!
   * \brief Helper method to initialize the internal array data-structures.
   *
   * \param [in] numNodes the number of nodes to store
   * \param [in] maxCapacity initial max capacity of the arrays
   */
  void initialize(IndexType numNodes, IndexType maxCapacity);

  /*!
   * \brief Helper method to check if the internal arrays are consistent.
   * \return status true if the arrays are consistent, else, false.
   */
  bool consistencyCheck() const;

  /// @}

#ifdef AXOM_MINT_USE_SIDRE
  sidre::Group* m_group;
#endif
  int m_ndims;
  axom::deprecated::MCArray<double>* m_coordinates[3] = {nullptr, nullptr, nullptr};

  DISABLE_COPY_AND_ASSIGNMENT(MeshCoordinates);
  DISABLE_MOVE_AND_ASSIGNMENT(MeshCoordinates);
};

//------------------------------------------------------------------------------
// In-lined method implementations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline bool MeshCoordinates::isExternal() const
{
  bool is_external = m_coordinates[0]->isExternal();
  bool consistent = true;
  for(int i = 1; i < m_ndims; ++i)
  {
    consistent &= m_coordinates[i]->isExternal() == is_external;
  }

  SLIC_WARNING_IF(!consistent, "External state not consistent.");
  return is_external;
}

//------------------------------------------------------------------------------
inline bool MeshCoordinates::isInSidre() const
{
#ifdef AXOM_MINT_USE_SIDRE
  bool is_in_sidre = m_group != nullptr;
  bool consistent = true;
  for(int i = 0; i < m_ndims; ++i)
  {
    consistent &= m_coordinates[i]->isInSidre() == is_in_sidre;
  }

  SLIC_WARNING_IF(!consistent, "Sidre state not consistent.");
  return m_coordinates[0]->isInSidre();
#else
  return false;
#endif
}

//------------------------------------------------------------------------------
inline IndexType MeshCoordinates::append(double x)
{
  SLIC_ASSERT(dimension() == 1);
  SLIC_ASSERT(m_coordinates[0] != nullptr);

  IndexType idx = numNodes();
  m_coordinates[0]->append(x);

  SLIC_ASSERT(idx == numNodes() - 1);
  SLIC_ASSERT(validIndex(idx));
  SLIC_ASSERT(consistencyCheck());

  return idx;
}

//------------------------------------------------------------------------------
inline IndexType MeshCoordinates::append(double x, double y)
{
  SLIC_ASSERT(dimension() == 2);
  SLIC_ASSERT(m_coordinates[0] != nullptr);
  SLIC_ASSERT(m_coordinates[1] != nullptr);

  IndexType idx = numNodes();
  m_coordinates[0]->append(x);
  m_coordinates[1]->append(y);

  SLIC_ASSERT(idx == numNodes() - 1);
  SLIC_ASSERT(validIndex(idx));
  SLIC_ASSERT(consistencyCheck());

  return idx;
}

//------------------------------------------------------------------------------
inline IndexType MeshCoordinates::append(double x, double y, double z)
{
  SLIC_ASSERT(dimension() == 3);
  SLIC_ASSERT(m_coordinates[0] != nullptr);
  SLIC_ASSERT(m_coordinates[1] != nullptr);
  SLIC_ASSERT(m_coordinates[2] != nullptr);

  IndexType idx = numNodes();
  m_coordinates[0]->append(x);
  m_coordinates[1]->append(y);
  m_coordinates[2]->append(z);

  SLIC_ASSERT(idx == numNodes() - 1);
  SLIC_ASSERT(validIndex(idx));
  SLIC_ASSERT(consistencyCheck());

  return idx;
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::append(const double* coords, IndexType n)
{
  SLIC_ASSERT(coords != nullptr);
  SLIC_ASSERT(n >= 0);

  if(m_ndims == 1)
  {
    SLIC_ASSERT(m_coordinates[0] != nullptr);
    m_coordinates[0]->append(coords, n);
    return;
  }

  IndexType original_size = numNodes();
  double* coord_arrays[3];
  for(int dim = 0; dim < m_ndims; ++dim)
  {
    SLIC_ASSERT(m_coordinates[dim] != nullptr);
    m_coordinates[dim]->resize(original_size + n);
    coord_arrays[dim] = getCoordinateArray(dim);
  }

  for(IndexType i = 0; i < n; ++i)
  {
    for(int dim = 0; dim < m_ndims; ++dim)
    {
      coord_arrays[dim][original_size + i] = coords[m_ndims * i + dim];
    }
  }

  SLIC_ASSERT(consistencyCheck());
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::append(const double* x, const double* y, IndexType n)
{
  SLIC_ASSERT(m_ndims == 2);
  SLIC_ASSERT(m_coordinates[0] != nullptr);
  SLIC_ASSERT(m_coordinates[1] != nullptr);
  SLIC_ASSERT(x != nullptr);
  SLIC_ASSERT(y != nullptr);
  SLIC_ASSERT(n >= 0);

  m_coordinates[0]->append(x, n);
  m_coordinates[1]->append(y, n);

  SLIC_ASSERT(consistencyCheck());
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::append(const double* x, const double* y, const double* z, IndexType n)
{
  SLIC_ASSERT(m_ndims == 3);
  SLIC_ASSERT(m_coordinates[0] != nullptr);
  SLIC_ASSERT(m_coordinates[1] != nullptr);
  SLIC_ASSERT(m_coordinates[2] != nullptr);
  SLIC_ASSERT(x != nullptr);
  SLIC_ASSERT(y != nullptr);
  SLIC_ASSERT(z != nullptr);
  SLIC_ASSERT(n >= 0);

  m_coordinates[0]->append(x, n);
  m_coordinates[1]->append(y, n);
  m_coordinates[2]->append(z, n);

  SLIC_ASSERT(consistencyCheck());
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::set(IndexType nodeID, double x)
{
  SLIC_ASSERT(dimension() == 1);
  SLIC_ASSERT(validIndex(nodeID));
  SLIC_ASSERT(m_coordinates[0] != nullptr);

  double* data = m_coordinates[0]->getData();
  SLIC_ASSERT(data != nullptr);
  data[nodeID] = x;

  SLIC_ASSERT(consistencyCheck());
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::set(IndexType nodeID, double x, double y)
{
  SLIC_ASSERT(dimension() == 2);
  SLIC_ASSERT(validIndex(nodeID));
  SLIC_ASSERT(m_coordinates[0] != nullptr);
  SLIC_ASSERT(m_coordinates[1] != nullptr);

  double* data = m_coordinates[0]->getData();
  SLIC_ASSERT(data != nullptr);
  data[nodeID] = x;

  data = m_coordinates[1]->getData();
  SLIC_ASSERT(data != nullptr);
  data[nodeID] = y;

  SLIC_ASSERT(consistencyCheck());
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::set(IndexType nodeID, double x, double y, double z)
{
  SLIC_ASSERT(dimension() == 3);
  SLIC_ASSERT(validIndex(nodeID));
  SLIC_ASSERT(m_coordinates[0] != nullptr);
  SLIC_ASSERT(m_coordinates[1] != nullptr);
  SLIC_ASSERT(m_coordinates[2] != nullptr);

  double* data = m_coordinates[0]->getData();
  SLIC_ASSERT(data != nullptr);
  data[nodeID] = x;

  data = m_coordinates[1]->getData();
  SLIC_ASSERT(data != nullptr);
  data[nodeID] = y;

  data = m_coordinates[2]->getData();
  SLIC_ASSERT(data != nullptr);
  data[nodeID] = z;

  SLIC_ASSERT(consistencyCheck());
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::insert(IndexType nodeID, double x)
{
  SLIC_ASSERT(dimension() == 1);
  SLIC_ASSERT(m_coordinates[0] != nullptr);
  SLIC_ASSERT(0 <= nodeID && nodeID <= numNodes());

  m_coordinates[0]->insert(x, nodeID);

  SLIC_ASSERT(consistencyCheck());
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::insert(IndexType nodeID, double x, double y)
{
  SLIC_ASSERT(dimension() == 2);
  SLIC_ASSERT(m_coordinates[0] != nullptr);
  SLIC_ASSERT(m_coordinates[1] != nullptr);
  SLIC_ASSERT(0 <= nodeID && nodeID <= numNodes());

  m_coordinates[0]->insert(x, nodeID);
  m_coordinates[1]->insert(y, nodeID);

  SLIC_ASSERT(consistencyCheck());
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::insert(IndexType nodeID, double x, double y, double z)
{
  SLIC_ASSERT(dimension() == 3);
  SLIC_ASSERT(m_coordinates[0] != nullptr);
  SLIC_ASSERT(m_coordinates[1] != nullptr);
  SLIC_ASSERT(m_coordinates[2] != nullptr);
  SLIC_ASSERT(0 <= nodeID && nodeID <= numNodes());

  m_coordinates[0]->insert(x, nodeID);
  m_coordinates[1]->insert(y, nodeID);
  m_coordinates[2]->insert(z, nodeID);

  SLIC_ASSERT(consistencyCheck());
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::insert(IndexType nodeID, const double* coords, IndexType n)
{
  SLIC_ASSERT(coords != nullptr);
  SLIC_ASSERT(n >= 0);
  SLIC_ASSERT(0 <= nodeID && nodeID <= numNodes());

  if(m_ndims == 1)
  {
    SLIC_ASSERT(m_coordinates[0] != nullptr);
    m_coordinates[0]->insert(coords, n, nodeID);
    return;
  }

  double* coord_arrays[3];
  for(int dim = 0; dim < m_ndims; ++dim)
  {
    SLIC_ASSERT(m_coordinates[dim] != nullptr);
    m_coordinates[dim]->emplace(n, nodeID);
    coord_arrays[dim] = getCoordinateArray(dim);
  }

  for(IndexType i = 0; i < n; ++i)
  {
    for(int dim = 0; dim < m_ndims; ++dim)
    {
      coord_arrays[dim][nodeID + i] = coords[m_ndims * i + dim];
    }
  }

  SLIC_ASSERT(consistencyCheck());
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::insert(IndexType nodeID, const double* x, const double* y, IndexType n)
{
  SLIC_ASSERT(m_ndims == 2);
  SLIC_ASSERT(m_coordinates[0] != nullptr);
  SLIC_ASSERT(m_coordinates[1] != nullptr);
  SLIC_ASSERT(x != nullptr);
  SLIC_ASSERT(y != nullptr);
  SLIC_ASSERT(0 <= nodeID && nodeID <= numNodes());
  SLIC_ASSERT(n >= 0);

  m_coordinates[0]->insert(x, n, nodeID);
  m_coordinates[1]->insert(y, n, nodeID);

  SLIC_ASSERT(consistencyCheck());
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::insert(IndexType nodeID,
                                    const double* x,
                                    const double* y,
                                    const double* z,
                                    IndexType n)
{
  SLIC_ASSERT(m_ndims == 3);
  SLIC_ASSERT(m_coordinates[0] != nullptr);
  SLIC_ASSERT(m_coordinates[1] != nullptr);
  SLIC_ASSERT(m_coordinates[2] != nullptr);
  SLIC_ASSERT(x != nullptr);
  SLIC_ASSERT(y != nullptr);
  SLIC_ASSERT(z != nullptr);
  SLIC_ASSERT(0 <= nodeID && nodeID <= numNodes());
  SLIC_ASSERT(n >= 0);

  m_coordinates[0]->insert(x, n, nodeID);
  m_coordinates[1]->insert(y, n, nodeID);
  m_coordinates[2]->insert(z, n, nodeID);

  SLIC_ASSERT(consistencyCheck());
}

//------------------------------------------------------------------------------
inline double MeshCoordinates::getCoordinate(IndexType nodeID, int dim) const
{
  SLIC_ASSERT(dim < m_ndims);
  SLIC_ASSERT((nodeID >= 0) && (nodeID < numNodes()));
  SLIC_ASSERT(m_coordinates[dim] != nullptr);

  return (*m_coordinates[dim])(nodeID);
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::getCoordinates(IndexType nodeID, double* coords) const
{
  SLIC_ASSERT(coords != nullptr);
  for(int dim = 0; dim < m_ndims; ++dim)
  {
    coords[dim] = getCoordinate(nodeID, dim);
  }
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::reserve(IndexType capacity)
{
  SLIC_ERROR_IF(isExternal() && capacity > this->capacity(),
                "cannot exceed initial capacity of external buffer!");

  for(int dim = 0; dim < m_ndims; ++dim)
  {
    SLIC_ASSERT(m_coordinates[dim] != nullptr);
    m_coordinates[dim]->reserve(capacity);
  }

  SLIC_ASSERT(consistencyCheck());
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::resize(IndexType size)
{
  for(int dim = 0; dim < m_ndims; ++dim)
  {
    SLIC_ASSERT(m_coordinates[dim] != nullptr);
    m_coordinates[dim]->resize(size);
  }

  SLIC_ASSERT(consistencyCheck());
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::setResizeRatio(double ratio)
{
  for(int dim = 0; dim < m_ndims; ++dim)
  {
    SLIC_ASSERT(m_coordinates[dim] != nullptr);
    m_coordinates[dim]->setResizeRatio(ratio);
  }

  SLIC_ASSERT(consistencyCheck());
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::shrink()
{
  for(int dim = 0; dim < m_ndims; ++dim)
  {
    SLIC_ASSERT(m_coordinates[dim] != nullptr);
    m_coordinates[dim]->shrink();
  }

  SLIC_ASSERT(consistencyCheck());
}

//------------------------------------------------------------------------------
// in-line  private method implementations
//------------------------------------------------------------------------------
inline void MeshCoordinates::initialize(IndexType numNodes, IndexType maxCapacity)
{
  SLIC_ERROR_IF((numNodes < 0), "supplied numNodes must be positive!");

  SLIC_ASSERT(numNodes <= maxCapacity);
  SLIC_ASSERT(!invalidDimension());

  for(int i = 0; i < m_ndims; ++i)
  {
    m_coordinates[i] = new axom::deprecated::MCArray<double>(numNodes, 1, maxCapacity);
  }

  SLIC_ASSERT(consistencyCheck());
}

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_MESHCOORDINATES_HPP_ */
