// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file View.hpp
 *
 * \brief   Header file containing definition of View class.
 *
 ******************************************************************************
 */

#ifndef SIDRE_VIEW_HPP_
#define SIDRE_VIEW_HPP_

// Standard C++ headers
#include <string>
#include <set>
#include <typeinfo>

// Other axom headers
#include "axom/config.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/core/Array.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"
#include "axom/slic.hpp"

// Sidre headers
#include "axom/sidre/core/ConduitMemory.hpp"
#include "axom/sidre/core/SidreTypes.hpp"
#include "axom/sidre/core/AttrValues.hpp"

namespace axom
{
namespace sidre
{
// Helper macro for defining a prepend string for sidre::View log messages
// We are using it to add the pathName() of the view
#ifndef SIDRE_VIEW_LOG_PREPEND
  #define SIDRE_VIEW_LOG_PREPEND "[View: '" << this->getPathName() << "'] "
#endif

class Buffer;
class Group;
class DataStore;
class Attribute;

/*!
 * \class View
 *
 * \brief A View object describes data, which may be
 *        owned by the view object (e.g., via an attached Buffer) or
 *        owned externally.
 *
 * The View class has the following properties:
 *
 *    - View objects can only be created via the Group interface,
 *      not constructed directly. A View object is owned by the Group
 *      object that creates it. A View object owned by a Group object
 *      that is a descendant of some ancestor Group is a descendant
 *      View of the ancestor Group.
 *    - A View object has a unique name (string) within the Group
 *      that owns it.
 *    - A View holds a pointer to the Group that created it and which
 *      owns it.
 *    - A View object can describe and provide access to data in one of
 *      four ways:
 *        * A view can describe (a subset of) data owned by an existing
 *          Buffer. In this case, the data can be (re)allocated or
 *          deallocated by the view if and only if it is the only view
 *          attached to the buffer.
 *        * A view can describe and allocate data using semantics similar
 *          to Buffer data description and allocation. In this case, no
 *          other view is allowed to (re)allocate or deallocate the data held
 *          by the associated data buffer.
 *        * A view can describe data associated with a pointer to an
 *          "external" data object. In this case, the view cannot (re)allocate
 *          or deallocate the data. However, all other view operations are
 *          essentially the same as the previous two cases.
 *        * It can hold a pointer to an undescribed (i.e., "opaque") data
 *          object. In this case, the view knows nothing about the type or
 *          structure of the data; it is essentially just a handle to the data.
 *    - For any View object that is "external" or associated with a
 *      Buffer, the data description of the view may be specified, or
 *      changed, by calling one of the apply() methods.
 *
 */
class View
{
public:
  //
  // Friend declaration to constrain usage via controlled access to
  // private members.
  //
  friend class Group;
  friend class Buffer;

  ///@{
  //!  @name View query and accessor methods

  /*!
   * \brief Return index of View within owning Group.
   *
   * If View is detached, return sidre::InvalidIndex.
   */
  IndexType getIndex() const { return m_index; }

  /*!
   * \brief Return const reference to name of View.
   *
   * \sa getPath(), getPathName()
   */
  const std::string& getName() const { return m_name; }

  /*!
   * \brief Return path of View's owning Group object.
   *
   * \sa getName(), getPathName()
   */
  std::string getPath() const;

  /*!
   * \brief Return full path of View object, including its name.
   *
   * If a DataStore contains a Group tree structure a/b/c/d/e, with
   * group d owning a view v, the following results are expected:
   *
   * Method Call      | Result
   * -----------------|----------
   * v->getName()     | v
   * v->getPath()     | a/b/c/d
   * v->getPathName() | a/b/c/d/v
   *
   * \sa getName(), getPath(), Group::getPathName()
   */
  std::string getPathName() const;

  /*!
   * \brief Return pointer to non-const Group that owns View object.
   */
  Group* getOwningGroup() { return m_owning_group; }

  /*!
   * \brief Return pointer to const Group that owns View object.
   */
  const Group* getOwningGroup() const { return m_owning_group; }

  /*!
   * \brief Return true if view has a an associated Buffer object
   *        (e.g., view is not opaque, it has been created with a buffer,
   *         it has been allocated, etc.); false otherwise.
   */
  bool hasBuffer() const { return m_data_buffer != nullptr; }

  /*!
   * \brief Return pointer to non-const Buffer associated with View.
   */
  Buffer* getBuffer() { return m_data_buffer; }

  /*!
   * \brief Return pointer to const Buffer associated with View.
   */
  const Buffer* getBuffer() const { return m_data_buffer; }

  /*!
   * \brief Return true if view holds external data; false otherwise.
   */
  bool isExternal() const { return m_state == EXTERNAL; }

  /*!
   * \brief Return true if view is described and refers to a buffer
   * that has been allocated.
   */
  bool isAllocated() const;

  /*!
   * \brief Return true if data description (schema) has been applied to data
   *        in buffer associated with view; false otherwise.
   */
  bool isApplied() const { return m_is_applied; }

  /*!
   * \brief Return true if data description exists.  It may/may not have been
   * applied to the data yet.  ( Check isApplied() for that. )
   */
  bool isDescribed() const { return !m_schema.dtype().is_empty(); }

  /*!
   * \brief Return true if view is empty.
   */
  bool isEmpty() const { return m_state == EMPTY; }

  /*!
   * \brief Convenience function that returns true if view is opaque
   *        (i.e., has access to data but has no knowledge of the data
   *        type or structure); false otherwise.
   */
  bool isOpaque() const { return m_state == EXTERNAL && !isApplied(); }

  /*!
   * \brief Return true if view contains a tuple or scalar value.
   */
  bool isTuple() const { return m_state == TUPLE; }

  /*!
   * \brief Return true if view contains a scalar value.
   */
  bool isScalar() const { return m_state == TUPLE && getNumElements() == 1; }

  /*!
   * \brief Return true if view contains a string value.
   */
  bool isString() const { return m_state == STRING; }

  /*!
   * \brief Return type of data for this View object.
   *        Return NO_TYPE_ID for an undescribed view.
   */
  TypeID getTypeID() const
  {
    if(isDescribed())
    {
      return static_cast<TypeID>(m_schema.dtype().id());
    }
    else
    {
      return NO_TYPE_ID;
    }
  }

  /*!
   * \brief Return total number of bytes associated with this View object.
   *
   * \attention This is the total bytes described by the view; they may not
   *            yet be allocated.
   */
  IndexType getTotalBytes() const { return m_schema.total_strided_bytes(); }

  /*!
   * \brief Return total number of elements described by this View object.
   *
   * \attention This is the number of elements described by the view;
   *            they may not yet be allocated.
   */
  IndexType getNumElements() const { return m_schema.dtype().number_of_elements(); }

  /*!
   * \brief Return number of bytes per element in the described view.
   *
   * \attention This is the number of bytes per element described by the view
   *            which may not yet be allocated.
   */
  IndexType getBytesPerElement() const { return m_schema.dtype().element_bytes(); }

  /*!
   * \brief Return the offset in number of elements for the data described by
   *  this View object.
   *
   * \warning The code currently assumes that offsets into a view are given in
   *  terms of whole elements.  It is an assertion error if this is not the
   *  case.  If you have a different use case, please talk to the Sidre team.
   *
   * \note View::getData() and View::getArray() already account for the offset
   *        and return a pointer to the first element in the array:
   *        View::getVoidPtr() does not account for the offset.
   *
   * \attention This function is based on the view description.  It does not
   * imply that the data is allocated.
   *
   * \return The offset, in terms of the number of elements, from the described
   *  array to the first element.
   */
  IndexType getOffset() const;

  /*!
   * \brief Return the stride in number of elements for the data described by
   *  this View object.
   *
   * \warning The code currently assumes that strides into a view are given in
   *  terms of whole elements.  It is an assertion error if this is not the
   *  case.  If you have a different use case, please talk to the Sidre team.
   *
   * \attention This function is based on the view description.  It does not
   * imply that the data is allocated.
   *
   * \return The stride, in terms of the number of elements, between elements in
   *  the described array.
   */
  IndexType getStride() const;

  /*!
   * \brief Return dimensionality of this View's data.
   *
   * \sa getShape()
   */
  int getNumDimensions() const { return static_cast<int>(m_shape.size()); }

  /*!
   * \brief Return number of dimensions in data view and fill in shape
   *        information of this data view object.
   *
   *  ndims - maximum number of dimensions to return.
   *  shape - user supplied buffer assumed to be ndims long.
   *
   *  Return the number of dimensions of the view.
   *  Return -1 if shape is too short to hold all dimensions.
   */
  int getShape(int ndims, IndexType* shape) const;

  /*!
   * \brief Return const reference to schema describing data.
   */
  const Schema& getSchema() const { return m_node.schema(); }

  /*!
   * \brief Return non-const reference to Conduit node holding data.
   */
  Node& getNode() { return m_node; }

  /*!
   * \brief Return const reference to Conduit node holding data.
   */
  const Node& getNode() const { return m_node; }

  /*!
   * \brief Returns boolean telling whether two Views have equivalent
   *  internal description, in terms of name, datatype, and current state of the
   *  object. Values of the data are not checked.
   */
  bool isEquivalentTo(const View* other) const;

  /*!
   * \brief Returns true if both Views are either associated with a buffer or
   * external, they span the same number of bytes and have unit stride.
   */
  bool isUpdateableFrom(const View* other) const;

  ///@}

  ///@{
  //!  @name View allocation methods

  /*!
   * \brief Allocate data for a view, previously described.
   *
   * \note Allocation from a view is allowed only if it is the only
   *       view associated with its buffer (when it has one), or the view
   *       is not external, not a string view, or not a scalar view.
   *       If none of these condition is true, this method does nothing.
   *
   * \return pointer to this View object.
   */
  View* allocate(int allocID = INVALID_ALLOCATOR_ID);

  /*!
   * \brief Allocate data for view given type and number of elements.
   *
   * \note The allocate() method (above) describes conditions where View
   *       allocation is allowed.  If the conditions are not met,
   *       type is NO_TYPE_ID, or num_elems < 0, this method does nothing.
   *
   * \return pointer to this View object.
   */
  View* allocate(TypeID type, IndexType num_elems, int allocID = INVALID_ALLOCATOR_ID);

  /*!
   * \brief Allocate data for view given type and shape.
   *
   * \note The allocate() method (above) describes conditions where View
   *       allocation is allowed.  If the conditions are not met,
   *       type is NO_TYPE_ID, or ndims < 0, or shape is nullptr,
   *       or any element of shape < 0, this method does nothing.
   *
   * If shape is nullptr, this is a no-op.
   *
   * \return pointer to this View object.
   */
  View* allocate(TypeID type, int ndims, const IndexType* shape, int allocID);

  /*!
   * \brief Allocate data for view described by a Conduit data type object.
   *
   * \note The allocate() method describes conditions where view
   *       allocation is allowed. If the conditions are not met,
   *       this method does nothing.
   *
   * \return pointer to this View object.
   */
  View* allocate(const DataType& dtype, int allocID = INVALID_ALLOCATOR_ID);

  /*!
   * \brief  Reallocate data for view to given number of elements (type
   *         stays the same).
   *
   * \note Reallocation from a view is only allowed under that same conditions
   *       for the allocate() method. If the conditions are not met
   *       or num_elems < 0 this method does nothing.
   *
   * \return pointer to this View object.
   */
  View* reallocate(IndexType num_elems);

  /*!
   * \brief Reallocate data for this View as specified by Conduit data type object.
   *
   * \note Reallocation from a view is allowed under the conditions
   *       described by the allocate() method. If the conditions are not met
   *       or dtype is undefined, this method does nothing.
   *
   * \note The given data type object must match the view type, if it is
   *       defined. If not, the method does nothing.
   *
   * \return pointer to this View object.
   */
  View* reallocate(const DataType& dtype);

  /*!
   * \brief  Deallocate data for view.
   *
   * \note Deallocation from a view is only allowed under the conditions
   *       described by the allocate() method. If the conditions are not met
   *       or a Buffer is not attached this method does nothing.
   *
   * \return pointer to this View object.
   */
  View* deallocate();

  ///@}

  /*!
   * \brief Reshape the array without changing its size.
   *
   * The View state must be either BUFFER or EXTERNAL.
   * The new shape must have the same size as the current.
   * If either conditions are not met, a warning is issued
   * and this is a no-op.
   *
   * \pre Old and new shapes must be the same size.
   * \pre !isEmpty() && (isExternal() || hasBuffer())
   * \post getNumDimensions() == ndims
   * \post getNumElements() is unchanged.
   */
  View* reshapeArray(int ndims, const IndexType* shape);

  /*!
   * \brief Attach Buffer object to data view.
   *
   * If the view has no description, then the buffer's description
   * is copied into the view.
   *
   * Note that, in general, the view cannot be used to access data in
   * buffer until one of the apply() methods is called. However, if
   * the view has a valid data description with a total number of bytes
   * that is <= number of bytes held in the buffer, then apply() will
   * be called internally.
   *
   * If data view already has a buffer, or it is an external view,
   * a scalar view, or a string view, this method does nothing.
   *
   * If data view already has a buffer and buff is nullptr, the attached
   * buffer will be detached. After the view is detached from the
   * buffer, if the buffer has no views attached to it, then it will
   * be destroyed.
   *
   * \return pointer to this View object.
   */
  View* attachBuffer(Buffer* buff);

  /*!
   * \brief Describe the data view and attach Buffer object.
   *
   * \return pointer to this View object.
   */
  View* attachBuffer(TypeID type, IndexType num_elems, Buffer* buff)
  {
    describe(type, num_elems);
    attachBuffer(buff);
    return this;
  }

  /*!
   * \brief Describe the data view and attach Buffer object.
   *
   * If shape is nullptr, this is a no-op.
   *
   * \return pointer to this View object.
   */
  View* attachBuffer(TypeID type, int ndims, const IndexType* shape, Buffer* buff)
  {
    SLIC_CHECK_MSG(shape != nullptr,
                   SIDRE_VIEW_LOG_PREPEND << "Could not allocate: specified shape is nullptr.");
    if(shape == nullptr)
    {
      return this;
    }
    describe(type, ndims, shape);
    attachBuffer(buff);
    return this;
  }

  /*!
   * \brief Detach this view from its Buffer.
   *
   * If the view has no buffer, the method does nothing.
   *
   * \return pointer to detached buffer.
   */
  Buffer* detachBuffer();

  /*!
   * \brief Clear data and metadata from a View.
   *
   * The view will be EMPTY. There will be no description
   * or data associated with the View.
   */
  void clear();

  ///@{
  //!  @name Methods to apply View description to data.

  /*!
   * \brief Apply view description to data.
   *
   * If view holds a scalar, tuple or a string, the method does nothing.
   *
   * \return pointer to this View object.
   */
  View* apply();

  /*!
   * \brief Apply data description defined by number of elements, and
   *        optionally offset and stride to data view (type remains the same).
   *
   * \note The units for offset and stride are in number of elements, which
   *       is different than the Conduit DataType usage below where offset
   *       and stride are in number of bytes.
   *
   * \attention If view has been previously described (or applied), this
   *            operation will apply the new data description to the view.
   *
   * If view holds a scalar or a string, is external and does not have a
   * sufficient data description to get type information, or given number
   * of elements < 0, or offset < 0, the method does nothing.
   *
   * \return pointer to this View object.
   */
  View* apply(IndexType num_elems, IndexType offset = 0, IndexType stride = 1);

  /*!
   * \brief Apply data description defined by type and number of elements, and
   *        optionally offset and stride to data view.
   *
   * \note The units for offset and stride are in number of elements, which
   *       is different than the Conduit DataType usage below where offset
   *       and stride are in number of bytes.
   *
   * \attention If view has been previously described (or applied), this
   *            operation will apply the new data description to the view.
   *
   * If view holds a scalar or a string, or type is NO_TYPE_ID,
   * or given number of elements < 0, or offset < 0, the method does nothing.
   *
   * \return pointer to this View object.
   */
  View* apply(TypeID type, IndexType num_elems, IndexType offset = 0, IndexType stride = 1);

  /*!
   * \brief Apply data description defined by type and shape information
   *        to data view.
   *
   * \note The units for the shape are in number of elements.
   *
   * \attention If view has been previously described (or applied), this
   *            operation will apply the new data description to the view.
   *
   * If view holds a scalar or a string, or type is NO_TYPE_ID,
   * or given number of dimensions < 0, or pointer to shape is null,
   * the method does nothing.
   *
   * \return pointer to this View object.
   */
  View* apply(TypeID type, int ndims, const IndexType* shape);

  /*!
   * \brief Apply data description of given Conduit data type to data view.
   *
   * If view holds a scalar or a string, the method does nothing.
   *
   * \return pointer to this View object.
   */
  View* apply(const DataType& dtype);

  ///@}

  ///@{
  //!  @name Methods to set data in the view (scalar, string, tuple or external data).

  /*!
   * \brief Set the view to hold the given tuple.
   *
   * \return pointer to this View object.
   */
  template <typename ScalarType>
  View* setTuple(const axom::ArrayView<ScalarType>& values, int allocID = INVALID_ALLOCATOR_ID)
  {
    // If this view already contains a tuple, issue a warning if the user is
    // changing the underlying type ( ie: integer -> float ).
#if defined(AXOM_DEBUG)
    if(m_state == TUPLE)
    {
      DataTypeId arg_id = detail::SidreTT<ScalarType>::id;
      SLIC_CHECK_MSG(arg_id == m_node.dtype().id(),
                     SIDRE_VIEW_LOG_PREPEND << "You are setting a scalar value which has changed "
                                            << " the underlying data type. "
                                            << "Old type: " << m_node.dtype().name() << ", "
                                            << "new type: " << DataType::id_to_name(arg_id) << ".");
    }
#endif

    // Note: most of these calls that set the view class members are
    //       unnecessary if the view already holds a tuple.
    if(m_state == EMPTY || m_state == TUPLE)
    {
      auto conduitAllocId = ConduitMemory::axomAllocIdToConduit(getValidTupleAllocatorId(allocID));
      m_node.set_allocator(conduitAllocId);
      std::vector<ScalarType> tmpValues(values.begin(), values.end());
      m_node.set(tmpValues);
      m_schema.set(m_node.schema());
      m_state = TUPLE;
      m_is_applied = true;
      describeShape();
    }
    else
    {
      SLIC_CHECK_MSG(m_state == EMPTY || m_state == TUPLE,
                     SIDRE_VIEW_LOG_PREPEND << "Unable to set scalar value on view "
                                            << " with state: " << getStateStringName(m_state));
    }
    return this;
  }

  /*!
   * \brief Set the view to hold the given scalar.
   *
   * \return pointer to this View object.
   */
  template <typename ScalarType>
  View* setScalar(ScalarType value, int allocID = INVALID_ALLOCATOR_ID)
  {
    // If this view already contains a scalar, issue a warning if the user is
    // changing the underlying type ( ie: integer -> float ).
#if defined(AXOM_DEBUG)
    if(m_state == TUPLE)
    {
      DataTypeId arg_id = detail::SidreTT<ScalarType>::id;
      SLIC_CHECK_MSG(arg_id == m_node.dtype().id(),
                     SIDRE_VIEW_LOG_PREPEND << "You are setting a scalar value which has changed "
                                            << " the underlying data type. "
                                            << "Old type: " << m_node.dtype().name() << ", "
                                            << "new type: " << DataType::id_to_name(arg_id) << ".");
    }
#endif

    // Note: most of these calls that set the view class members are
    //       unnecessary if the view already holds a scalar.
    if(m_state == EMPTY || m_state == TUPLE)
    {
      auto conduitAllocId = ConduitMemory::axomAllocIdToConduit(getValidTupleAllocatorId(allocID));
      m_node.set_allocator(conduitAllocId);
      m_node.set(value);
      m_schema.set(m_node.schema());
      m_state = TUPLE;
      m_is_applied = true;
      describeShape();
    }
    else
    {
      SLIC_CHECK_MSG(m_state == EMPTY || m_state == TUPLE,
                     SIDRE_VIEW_LOG_PREPEND << "Unable to set scalar value on view "
                                            << " with state: " << getStateStringName(m_state));
    }
    return this;
  }

  /*!
   * \brief Set the view to hold the given scalar.
   *
   * \return pointer to this View object.
   */
  View* setScalar(Node& value, int allocID = INVALID_ALLOCATOR_ID)
  {
    // If this view already contains a scalar, issue a warning if the user is
    // changing the underlying type ( ie: integer -> float ).
#if defined(AXOM_DEBUG)
    if(m_state == TUPLE)
    {
      SLIC_CHECK_MSG(value.dtype().id() == m_node.dtype().id(),
                     SIDRE_VIEW_LOG_PREPEND
                       << "Setting a scalar value in view  which has changed "
                       << "the underlying data type."
                       << "Old type: " << m_node.dtype().name() << ", "
                       << "New type: " << DataType::id_to_name(value.dtype().id()) << ".");
    }
#endif

    // Note: most of these calls that set the view class members are
    //       unnecessary if the view already holds a scalar.
    if(m_state == EMPTY || m_state == TUPLE)
    {
      auto conduitAllocId = ConduitMemory::axomAllocIdToConduit(getValidAllocatorId(allocID));
      m_node.set_allocator(conduitAllocId);
      m_node.set(value);
      m_schema.set(m_node.schema());
      m_state = TUPLE;
      m_is_applied = true;
      describeShape();
    }
    else
    {
      SLIC_CHECK_MSG(m_state == EMPTY || m_state == TUPLE,
                     SIDRE_VIEW_LOG_PREPEND << "Unable to set scalar value on view with state: "
                                            << getStateStringName(m_state));
    }
    return this;
  }

  /*
 * \brief set the View to hold an array in a Buffer
 *
 * This takes a Node that holds an array of data and sets up the View to
 * hold a copy of that data in an attached Buffer.
 *
 * The prerequisites are that the Node must have numerical conduit data type
 * (it returns true for a call to conduit::DataType::is_number()), and the
 * state of this View must be either EMPTY or BUFFER when entering.  If the
 * prerequisites are not met, a warning will be issued and this View will be
 * unchanged.
 *
 * If this View has the state BUFFER when this method is called, it will
 * become detached from its previous Buffer.
 *
 * \return pointer to this View object
 */
  View* importArrayNode(const Node& array);

  //
  // RDH -- Add an overload of the following that takes a const char *.
  //
  /*!
 * \brief Set the view to hold the given string.
 *
 * \return pointer to this View object.
 */
  View* setString(const std::string& value, int allocID = INVALID_ALLOCATOR_ID)
  {
    // Note: most of these calls that set the view class members are
    //       unnecessary if the view already holds a string.
    if(m_state == EMPTY || m_state == STRING)
    {
      m_state = STRING;
      auto conduitAllocId = ConduitMemory::axomAllocIdToConduit(getValidAllocatorId(allocID));
      m_node.set_allocator(conduitAllocId);
      m_node.set_string(value);
      m_schema.set(m_node.schema());
      m_is_applied = true;
      describeShape();
    }
    else
    {
      SLIC_CHECK_MSG(m_state == EMPTY || m_state == STRING,
                     SIDRE_VIEW_LOG_PREPEND << "Unable to set string value on view with state: "
                                            << getStateStringName(m_state));
    }
    return this;
  };

  /*!
   * \brief Set view to hold external data.
   *
   * Data is undescribed (i.e., view is opaque) until an apply methods
   * is called on the view.
   *
   * If external_ptr is nullptr, the view will be EMPTY.
   * Any existing description is unchanged.
   *
   * \return pointer to this View object.
   */
  View* setExternalDataPtr(void* external_ptr);

  /*!
   * \brief Set view to hold described external data.
   *
   * If external_ptr is nullptr, the view will be EMPTY.
   *
   * \return pointer to this View object.
   */
  View* setExternalDataPtr(TypeID type, IndexType num_elems, void* external_ptr)
  {
    describe(type, num_elems);
    setExternalDataPtr(external_ptr);
    return this;
  }

  /*!
   * \brief Set view to hold described external data.
   *
   * If external_ptr is nullptr, the view will be EMPTY.
   *
   * If shape is nullptr, this is a no-op.
   *
   * \return pointer to this View object.
   */
  View* setExternalDataPtr(TypeID type, int ndims, const IndexType* shape, void* external_ptr)
  {
    SLIC_CHECK_MSG(
      shape != nullptr,
      SIDRE_VIEW_LOG_PREPEND << "Could not set external data ptr: specified shape is nullptr.");
    if(shape == nullptr)
    {
      return this;
    }

    describe(type, ndims, shape);
    setExternalDataPtr(external_ptr);
    return this;
  }

  ///@}

  /*!
 * \brief Update the data in this View with the data in other
 * if isUpdateableFrom( other ). Otherwise nothing is done.
 *
 * \return pointer to this View object.
 */
  View* updateFrom(const View* other);

  ///@{
  //! @name Methods to retrieve data in a view.

  /*!
   * \brief Return a pointer or conduit array object to the view's array data.
   *
   * Return value depends on variable type caller assigns it to. For example,
   * if view holds an integer array, the following usage is possible:
   *
   *      int* a = view->getArray();      // Get array as int pointer
   *      int_array a = view->getArray(); // Get array as Conduit array struct.
   *
   * \note The returned pointer accounts for the View's offset, so getArray()[0]
   *       always points to the first element in the array.
   */
  Node::Value getArray()
  {
    //TODO add check that view holds array data.  Will be added in later commit.
    //If debug, should trigger assert.  If release, issue warning.
    return getData();
  }

  /// \overload
  Node::ConstValue getArray() const { return getData(); }

  /*!
   * \brief Returns a pointer to the string contained in the view.
   *
   *  If the view is not a STRING, then nullptr is returned.
   */
  const char* getString() const
  {
    if(m_state == STRING)
    {
      return m_node.as_char8_str();
    }
    else
    {
      return nullptr;
    }
  }

  /*!
   * \brief Returns a copy of the scalar value contained in the view.
   */
  Node::ConstValue getScalar() const
  {
    SLIC_CHECK_MSG(isScalar(),
                   SIDRE_VIEW_LOG_PREPEND << "View::getScalar() called on non-tuple view.");
    return getData();
  }

  /*!
   * \brief Return data held by view and cast it to any compatible type
   *  allowed by Conduit (return type depends on type caller assigns it to).
   *
   *  If view does not contain allocated data, an empty Node::Value will be
   *  returned.
   *
   *  \note The return value already accounts for the View's offset
   *   (when present), so, if the View is an array, getData()[0] already points
   *   to the first element
   */
  /// @{
  Node::Value getData()
  {
    SLIC_CHECK_MSG(
      isAllocated(),
      SIDRE_VIEW_LOG_PREPEND << "No view data present, memory has not been allocated.");
    SLIC_CHECK_MSG(isDescribed(), SIDRE_VIEW_LOG_PREPEND << "View data description not present.");

    // this will return a default value
    return m_node.value();
  }

  Node::ConstValue getData() const
  {
    SLIC_CHECK_MSG(
      isAllocated(),
      SIDRE_VIEW_LOG_PREPEND << "No view data present, memory has not been allocated.");
    SLIC_CHECK_MSG(isDescribed(), SIDRE_VIEW_LOG_PREPEND "View data description not present.");

    // this will return a default value
    return m_node.value();
  }
  /// @}

  /*!
   * \brief Lightweight templated wrapper around getData() that can be used when
   *  you are calling getData(), but not assigning the return type.
   *
   * \sa getData()
   */
  template <typename DataType>
  DataType getData()
  {
    DataType data = m_node.value();
    return data;
  }

  /*!
   * \brief Returns a void pointer to the view's data
   *
   * \note This function returns the base pointer that was used to set up the
   *  view. It does not account for any offsets or strides in the View's
   *  description.
   *
   * To access the first data element, you will need to cast to the appropriate
   * type and add the offset.  E.g. if the underlying data is an array of
   * integers you can access the first element as follows:
   *
   *      void* vptr = view->getVoidPtr();
   *      int*  iptr = static_cast<int*>(vptr) + view->getOffset();
   *
   *
   * \sa getData(), getArray()
   */
  void* getVoidPtr() const;

  ///@}

  ///@{
  //!  @name View print methods.

  /*!
   * \brief Print JSON description of data view to stdout.
   */
  void print() const;

  /*!
   * \brief Print JSON description of data view to an ostream.
   */
  void print(std::ostream& os) const;

  /*!
   * \brief Print data in a way that won't crash for non-host data.
   *
   * If data is not host-accessible, print the pointer and a comment.
   */
  void hostPrint(std::ostream& os = std::cout) const;

  ///@}

  /*!
   * \brief Copy data view description to given Conduit node.
   */
  void copyToConduitNode(Node& n) const;

  /*!
   * \brief Copy data view native layout to given Conduit node.
   *
   * The native layout is a Conduit Node hierarchy that maps the Conduit Node
   * data externally to the Sidre View data so that it can be filled in from the
   * data in the file (independent of file format) and can be accessed as a
   * Conduit tree.
   */
  void createNativeLayout(Node& n) const;

  /*!
   * \brief Deep copy View into the given conduit::Node.
   * \param dst [in/out] Destination
   * \param allocId [in] If not equal to INVALID_ALLOCATOR_ID,
   *   use this allocator for the destination.
   */
  void deepCopyToConduit(Node& dst, int allocId = INVALID_ALLOCATOR_ID) const;

  /*!
   * \brief Copy metadata of the View to the given Conduit node
   *
   * The View's datatype schema and state information will be copied into the
   * node, but not the data held by the View.
   */
  void copyMetadataToNode(Node& n) const;

  /*!
   * \brief Change the name of this View.
   *
   * The name of this view is changed to the new name.  This also changes
   * the name for this view held by the owning group.
   *
   * Warnings will occur and the name will not be changed under these
   * conditions:  If the new name is an empty string, if the new name
   * contains a path delimiter (usually '/'), or if the new name is
   * identical to a name that is already held by the parent for another
   * Group or View object.
   *
   * \param new_name    The new name for this view.
   *
   * \return            Success or failure of rename.
   */
  bool rename(const std::string& new_name);

  ///@{
  //!  @name Attribute Value query and accessor methods

  /*!
   * \brief Get an Attribute.
   *
   * \param [in] idx Attribute identifier.
   *
   * \return the specified Attribute, or nullptr if it doesn't exist.
   */
  Attribute* getAttribute(IndexType idx);

  /*!
   * \overload
   */
  const Attribute* getAttribute(IndexType idx) const;

  /*!
   * \overload
   */
  Attribute* getAttribute(const std::string& name);

  /*!
   * \overload
   */
  const Attribute* getAttribute(const std::string& name) const;

  /*!
   * \brief Return true if the attribute has been explicitly set; else false.
   */
  bool hasAttributeValue(IndexType idx) const
  {
    const Attribute* attr = getAttribute(idx);
    return m_attr_values.hasValue(attr);
  }

  /*!
   * \brief Return true if the attribute has been explicitly set; else false.
   */
  bool hasAttributeValue(const std::string& name) const
  {
    const Attribute* attr = getAttribute(name);
    return m_attr_values.hasValue(attr);
  }

  /*!
   * \brief Return true if the attribute has been explicitly set; else false.
   */
  bool hasAttributeValue(const Attribute* attr) const
  {
    SLIC_CHECK_MSG(attr != nullptr,
                   SIDRE_VIEW_LOG_PREPEND << "hasAttributeValue: called with a null Attribute");

    return m_attr_values.hasValue(attr);
  }

  /*!
   * \brief Set Attribute to its default value.
   * \param [in] idx Attribute identifier
   *
   * This causes hasAttributeValue to return false for the attribute.
   */
  bool setAttributeToDefault(IndexType idx)
  {
    const Attribute* attr = getAttribute(idx);
    return m_attr_values.setToDefault(attr);
  }

  /*!
   * \overload
   */
  bool setAttributeToDefault(const std::string& name)
  {
    const Attribute* attr = getAttribute(name);
    return m_attr_values.setToDefault(attr);
  }

  /*!
   * \overload
   */
  bool setAttributeToDefault(const Attribute* attr)
  {
    SLIC_CHECK_MSG(attr != nullptr,
                   SIDRE_VIEW_LOG_PREPEND << "getAttributeToDefault: called with a null Attribute");

    return m_attr_values.setToDefault(attr);
  }

  /*!
   * \brief Set Attribute to a scalar value.
   * \param [in] idx Attribute identifier
   * \param [in] value Value of Attribute
   *
   * If no such Attribute exists, this is a no-op.
   */
  template <typename ScalarType>
  bool setAttributeScalar(IndexType idx, ScalarType value)
  {
    const Attribute* attr = getAttribute(idx);
    if(attr == nullptr)
    {
      return false;
    }

    return m_attr_values.setScalar(attr, value);
  }

  /*!
   * \overload
   */
  template <typename ScalarType>
  bool setAttributeScalar(const std::string& name, ScalarType value)
  {
    const Attribute* attr = getAttribute(name);
    if(attr == nullptr)
    {
      return false;
    }

    return m_attr_values.setScalar(attr, value);
  }

  /*!
   * \overload
   */
  template <typename ScalarType>
  bool setAttributeScalar(const Attribute* attr, ScalarType value)
  {
    if(attr == nullptr)
    {
      SLIC_CHECK_MSG(attr != nullptr,
                     SIDRE_VIEW_LOG_PREPEND << "setAttributeScalar: called with a null Attribute");
      return false;
    }

    return m_attr_values.setScalar(attr, value);
  }

  /*!
   * \brief Set Attribute to a string value.
   * \param [in] indx Attribute identifier.
   * \param [in] value Value of attribute.
   *
   * If no such Attribute exists, this is a no-op.
   */
  bool setAttributeString(IndexType indx, const std::string& value);

  /*!
   * \overload
   */
  bool setAttributeString(const std::string& name, const std::string& value);

  /*!
   * \overload
   */
  bool setAttributeString(const Attribute* attr, const std::string& value);

  /*!
   * \brief Return scalar Attribute value.
   * \param [in] idx Attribute identifier.
   *
   * If no such Attribute exists, return empty scalar value.
   */
  Node::ConstValue getAttributeScalar(IndexType idx) const
  {
    const Attribute* attr = getAttribute(idx);
    if(attr == nullptr)
    {
      return m_attr_values.getEmptyNodeRef().value();
    }

    return m_attr_values.getScalar(attr);
  }

  /*!
   * \overload
   */
  Node::ConstValue getAttributeScalar(const std::string& name) const
  {
    const Attribute* attr = getAttribute(name);
    if(attr == nullptr)
    {
      return m_attr_values.getEmptyNodeRef().value();
    }

    return m_attr_values.getScalar(attr);
  }

  /*!
   * \overload
   */
  Node::ConstValue getAttributeScalar(const Attribute* attr) const
  {
    if(attr == nullptr)
    {
      SLIC_CHECK_MSG(attr != nullptr,
                     SIDRE_VIEW_LOG_PREPEND << "getScalar: called with a null Attribute");
      return m_attr_values.getEmptyNodeRef().value();
    }

    return m_attr_values.getScalar(attr);
  }

  /*!
   * \brief Lightweight templated wrapper around getAttributeScalar()
   *  that can be used when you are calling getAttributeScalar(), but not
   *  assigning the return type.
   *
   * \sa getAttributeScalar()
   */
  template <typename DataType>
  DataType getAttributeScalar(IndexType idx)
  {
    const Attribute* attr = getAttribute(idx);
    const Node& node = m_attr_values.getValueNodeRef(attr);
    DataType data = node.value();
    return data;
  }

  /*!
   * \brief Lightweight templated wrapper around getAttributeScalar()
   *  that can be used when you are calling getAttributeScalar(), but not
   *  assigning the return type.
   *
   * \sa getAttributeScalar()
   */
  template <typename DataType>
  DataType getAttributeScalar(const std::string& name)
  {
    const Attribute* attr = getAttribute(name);
    const Node& node = m_attr_values.getValueNodeRef(attr);
    DataType data = node.value();
    return data;
  }

  /*!
   * \brief Lightweight templated wrapper around getAttributeScalar()
   *  that can be used when you are calling getAttributeScalar(), but not
   *  assigning the return type.
   *
   * \sa getAttributeScalar()
   */
  template <typename DataType>
  DataType getAttributeScalar(const Attribute* attr)
  {
    SLIC_CHECK_MSG(attr != nullptr,
                   SIDRE_VIEW_LOG_PREPEND << "getAttributeScalar: called with a null Attribute");

    const Node& node = m_attr_values.getValueNodeRef(attr);
    DataType data = node.value();
    return data;
  }

  /*!
   * \brief Return a string attribute.
   * \param [in] idx Attribute identifer.
   *
   * If the value has not been explicitly set, return the current default.
   *
   * \return The Attribute name, or nullptr if no such Attribute.
   */
  const char* getAttributeString(IndexType idx) const;

  /*!
   * \overload
   */
  const char* getAttributeString(const std::string& name) const;

  /*!
   * \overload
   */
  const char* getAttributeString(const Attribute* attr) const;

  /*!
   * \brief Return reference to an Attribute node.
   * \param [in] idx Attribute identifier
   *
   * If the value has not been explicitly set, return the current default.
   *
   * \return The Attribute Node, or an empty Node if no such Attribute.
   */
  const Node& getAttributeNodeRef(IndexType idx) const
  {
    const Attribute* attr = getAttribute(idx);
    return m_attr_values.getValueNodeRef(attr);
  }

  /*!
   * \overload
   */
  const Node& getAttributeNodeRef(const std::string& name) const
  {
    const Attribute* attr = getAttribute(name);
    return m_attr_values.getValueNodeRef(attr);
  }

  /*!
   * \overload
   */
  const Node& getAttributeNodeRef(const Attribute* attr) const
  {
    SLIC_CHECK_MSG(attr != nullptr,
                   SIDRE_VIEW_LOG_PREPEND << "getAttributeNodeRef: called with a null Attribute");

    return m_attr_values.getValueNodeRef(attr);
  }

  /*!
   * \brief Return first valid Attribute index for a set Attribute in
   *        View object (i.e., smallest index over all Attributes).
   *
   * sidre::InvalidIndex is returned if View has no Attributes set.
   */
  IndexType getFirstValidAttrValueIndex() const
  {
    return m_attr_values.getFirstValidAttrValueIndex();
  }

  /*!
   * \brief Return next valid Attribute index for a set Attribute in
   *        View object after given index (i.e., smallest index over
   *        all Attribute indices larger than given one).
   *
   * sidre::InvalidIndex is returned if there is no valid index greater
   * than given one.
   * getNextAttrValueIndex(InvalidIndex) returns InvalidIndex.
   */
  IndexType getNextValidAttrValueIndex(IndexType idx) const
  {
    return m_attr_values.getNextValidAttrValueIndex(idx);
  }

  ///@}

private:
  DISABLE_DEFAULT_CTOR(View);
  DISABLE_MOVE_AND_ASSIGNMENT(View);

  ///@{
  //!  @name Private View ctor and dtor
  //!        (callable only by Group and View methods).

  /*!
   *  \brief Private ctor that creates a View with given name
   *         which has no data associated with it.
   */
  View(const std::string& name);

  /*!
   * \brief Private copy ctor.
   */
  View(const View& source);

  /*!
   * \brief Private dtor.
   */
  ~View();

  ///@}

  ///@{
  //!  @name Private View declaration methods.
  //!        (callable only by Group and View methods).

  /*!
   * \brief Describe a data view with given type and number of elements.
   *
   *
   * \attention If view has been previously described, this operation will
   *            re-describe the view. To have the new description take effect,
   *            the apply() method must be called.
   *
   * If given type of NO_TYPE_ID, or number of elements < 0, or view is opaque,
   * method does nothing.
   */
  void describe(TypeID type, IndexType num_elems);

  /*!
   * \brief Describe a data view with given type, number of dimensions, and
   *        number of elements per dimension.
   *
   * \attention If view has been previously described, this operation will
   *            re-describe the view. To have the new description take effect,
   *            the apply() method must be called.
   *
   * If given type of NO_TYPE_ID, or number of dimensions or total
   * number of elements < 0, or view is opaque, method does nothing.
   *
   * If shape is nullptr, this is a no-op.
   */
  void describe(TypeID type, int ndims, const IndexType* shape);

  /*!
   * \brief Declare a data view with a Conduit data type object.
   *
   * \attention If view has been previously described, this operation will
   *            re-describe the view. To have the new description take effect,
   *            the apply() method must be called.
   *
   * If view is opaque, the method does nothing.
   */
  void describe(const DataType& dtype);

  /*!
   * \brief Set the shape to be a one dimension with the described number of
   * elements.
   */
  void describeShape();

  /*!
   * \brief Set the shape to be a ndims dimensions with shape.
   *
   * If shape is nullptr, this is a no-op.
   */
  void describeShape(int ndims, const IndexType* shape);

  /*!
   * \brief Private method to remove user provided description.
   */
  void undescribe()
  {
    m_schema.reset();
    m_shape.clear();
  }

  /*!
   * \brief Copy contents of this View contents into an undescribed EMPTY View.
   *
   * For TUPLE and STRING the data is copied; EXTERNAL,
   * data pointer is copied; BUFFER attaches the buffer.
   */
  void copyView(View* copy) const;

  /*!
   * \brief Deep copy this View into an undescribed EMPTY View.
   *
   * For TUPLE and STRING the data is copied and the state is preserved.
   * For BUFFER and EXTERNAL, the data described by this View is copied into a
   * new Buffer that is of the size needed to hold the copied data. Any
   * parts of the source Buffer or external array that are not seen due
   * to offsets and strides in the description will not be copied. The copied
   * View will have BUFFER state with zero offset and a stride of one.
   */
  void deepCopyView(View* copy,
                    int arrayAllocId = INVALID_ALLOCATOR_ID,
                    int tupleAllocID = INVALID_ALLOCATOR_ID) const;

  /*!
   * \brief Add view description and references to it's data to a conduit tree.
   */
  void exportTo(conduit::Node& data_holder, std::set<IndexType>& buffer_indices) const;

  /*!
   * \brief Restore a view's description and data from a conduit tree.
   * This does not include a view's buffer data, that is done in the buffer
   */
  void importFrom(conduit::Node& data_holder, const std::map<IndexType, IndexType>& buffer_id_map);

  /*!
   * \brief Add view's description to a conduit tree.
   */
  void exportDescription(conduit::Node& data_holder) const;

  /*!
   * \brief Restore a view's description from a conduit tree.
   */
  void importDescription(conduit::Node& data_holder);

  /*!
   * \brief Add view's attributes to a conduit tree.
   */
  void exportAttribute(conduit::Node& data_holder) const;

  /*!
   * \brief Restore a view's attributes from a conduit tree.
   */
  void importAttribute(conduit::Node& data_holder);

  /*!
   *  \brief Private method to remove any applied description;
   *         but preserves user provided description.
   *
   *  Note: The description is stored in m_schema.
   */
  void unapply()
  {
    m_node.reset();
    m_is_applied = false;
  }

  /*!
   *  \brief Private method to reset a view from BUFFER to EMPTY.
   *
   *         Used by Buffer when detaching from a view.
   */
  void setBufferViewToEmpty()
  {
    m_data_buffer = nullptr;
    m_state = EMPTY;
    unapply();
  }

  ///@}

  ///@{
  //!  @name Private methods that indicate when certain view operations are valid.

  /*!
   *  \brief Private method returns true if data allocation on view is a
   *         valid operation; else false
   */
  bool isAllocateValid() const;

  /*!
   *  \brief Private method returns true if apply is a valid operation on
   *         view; else false
   */
  bool isApplyValid() const;

  ///@}

  ///
  /// Enum with constants that identify the state of a view.
  ///
  /// Note that these states are mutually-exclusive. These constants
  /// combined with the boolean m_is_applied uniquely identify the view
  /// state, or how it was created and defined.
  ///
  enum State
  {
    EMPTY,     // View created with name only :
               //    has no data or data description
    BUFFER,    // View has a buffer attached, either via call to
               // View::attachBuffer(), View::allocate(), etc. :
               //    applied may be true or false
    EXTERNAL,  // View holds pointer to external data (no buffer) :
               //    applied may be true or false
    TUPLE,     // View holds tuple (including scalar) data (via setScalar() or setTuple()):
               //    applied is true
    STRING     // View holds string data (via setString()):
               //    applied is true
  };

  /*!
   *  \brief Private method returns string name of given view state enum value.
   */
  static char const* getStateStringName(State state);

  /*!
   *  \brief Private method returns state enum value give a state name.
   */
  State getStateId(const std::string& name) const;

  /*!
   * If allocID == INVALID_ALLOCATOR_ID, return the default allocator id,
   * which depends on the View's data semantic and owning Group.
   *
   * Note: The data semantic depends on m_state, so that must be set first.
   */
  int getValidAllocatorId(int allocId);

  int getValidTupleAllocatorId(int allocId);

  int getValidArrayAllocatorId(int allocId);

  /*!
   * \brief Return whether view data is accessible on the host CPU,
   * as determined by Axom's memory management.
   *
   * If data is external and not allocated through Umpire or malloc,
   * the behavior is undefined.
   *
   * By convention, this returns a false if data is not allocated,
   * because we expect null pointers to be correctly checked before use.
   */
  bool isHostAccessible() const;

  //!@brief Print as a single line, in a way that won't crash for non-host data.
  template <typename T>
  void hostPrintScalar(std::ostream& os = std::cout) const
  {
    if(isHostAccessible())
    {
      os << ' ' << T(m_node.value());
    }
    else
    {
      os << ' ' << getVoidPtr() << " # non-host " << typeid(T).name() << " data";
    }
  }

  //!@brief Print on host, as a single line.
  template <typename T>
  void hostPrintArray(std::ostream& os = std::cout) const
  {
    constexpr IndexType maxPrintItems = 10;
    if(isHostAccessible())
    {
      os << " [";
      auto start = (T*)(getVoidPtr());
      auto end = (T*)(getVoidPtr()) + getNumElements();
      if(getNumElements() <= maxPrintItems)
      {
        for(auto i = start; i < end; ++i)
        {
          os << *i;
          if(i != end - 1) os << ", ";
        }
      }
      else
      {
        auto a = start + maxPrintItems / 2;
        auto b = end - maxPrintItems / 2;
        for(auto i = start; i < a; ++i)
        {
          os << *i << ", ";
        }
        os << "..., ";
        for(auto i = b; i < end; ++i)
        {
          os << *i;
          if(i != end - 1) os << ", ";
        }
      }
      os << "]";
    }
    else
    {
      axom::Array<axom::IndexType> shape(getNumDimensions());
      getShape(static_cast<int>(shape.size()), shape.data());
      os << ' ' << getVoidPtr() << " # non-host " << typeid(T).name() << " array of (" << shape[0];
      for(axom::IndexType i = 1; i < shape.size(); ++i)
      {
        os << " x " << shape[i];
      }
      os << ") elements";
    }
  }

  /// Name of this View object.
  std::string m_name;

  /// Index of this View object within m_owning_group.
  IndexType m_index;

  /// Group object that owns this View object.
  Group* m_owning_group;

  /// Buffer associated with this View object.
  Buffer* m_data_buffer;

  /// Data description (schema) that describes the view's data.
  Schema m_schema;

  /// Conduit node used to access the data in this View.
  Node m_node;

  /// Shape information
  std::vector<IndexType> m_shape;

  /// Pointer to external memory
  void* m_external_ptr;

  /// State of view.
  State m_state;

  /// Has data description been applied to the view's data?
  bool m_is_applied;

  /// Attribute Values
  AttrValues m_attr_values;
};

} /* end namespace sidre */
} /* end namespace axom */

#endif /* SIDRE_VIEW_HPP_ */
