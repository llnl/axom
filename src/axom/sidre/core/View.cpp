// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Associated header file
#include "View.hpp"

// Sidre component headers
#include "Buffer.hpp"
#include "Group.hpp"
#include "DataStore.hpp"
#include "Attribute.hpp"
#include "SidreTypes.hpp"

#include "axom/core/execution/execution_space.hpp"
#include "axom/core/Macros.hpp"
#include "axom/sidre/core/ConduitMemory.hpp"

namespace axom
{
namespace sidre
{
/*
 *************************************************************************
 *
 * Return path of View's owning Group object.
 * Needs to be in the .cpp file because Group methods aren't
 * accessible in .hpp file.
 *
 *************************************************************************
 */
std::string View::getPath() const { return getOwningGroup()->getPathName(); }

/*
 *************************************************************************
 *
 * Return full path of View object, including its name.
 * Needs to be in the .cpp file because Group methods aren't
 * accessible in .hpp file.
 *
 *************************************************************************
 */
std::string View::getPathName() const
{
  const auto path = getPath();

  if(path.length() < 1)
  {
    return getName();
  }

  return path + getOwningGroup()->getPathDelimiter() + getName();
}

/*
 *************************************************************************
 *
 * Allocate data for view, previously described.
 * The state may transition from EMPTY to BUFFER;
 * otherwise, the state must already be BUFFER.
 *
 *************************************************************************
 */
View* View::allocate(int allocID)
{
  if(isAllocateValid())
  {
    if(m_state == EMPTY)
    {
      SLIC_ASSERT_MSG(m_data_buffer == nullptr,
                      SIDRE_VIEW_LOG_PREPEND << "State was EMPTY, but data buffer was not null.");
      m_data_buffer = m_owning_group->getDataStore()->createBuffer();
      m_data_buffer->attachToView(this);
      m_state = BUFFER;
    }

    TypeID type = static_cast<TypeID>(m_schema.dtype().id());
    IndexType num_elems = m_schema.dtype().number_of_elements();
    allocID = getValidAllocatorId(allocID);
    m_data_buffer->allocate(type, num_elems, allocID);
    apply();
  }

  return this;
}

/*
 *************************************************************************
 *
 * Allocate array data for view with type and number of elements.
 *
 *************************************************************************
 */
View* View::allocate(TypeID type, IndexType num_elems, int allocID)
{
  if(type == NO_TYPE_ID || num_elems < 0)
  {
    SLIC_CHECK_MSG(type != NO_TYPE_ID,
                   SIDRE_VIEW_LOG_PREPEND << "Could not allocate: Data type was 'NO_TYPE_ID'.");
    SLIC_CHECK_MSG(
      num_elems >= 0,
      SIDRE_VIEW_LOG_PREPEND << "Could not allocate: num_elems cannot be less than zero.");

    return this;
  }

  describe(type, num_elems);
  allocID = getValidArrayAllocatorId(allocID);
  allocate(allocID);

  return this;
}

/*
 *************************************************************************
 *
 * Allocate array data for view with type and shape.
 *
 *************************************************************************
 */
View* View::allocate(TypeID type, int ndims, const IndexType* shape, int allocID)
{
  SLIC_CHECK_MSG(ndims > 0, SIDRE_VIEW_LOG_PREPEND << "Could not allocate: ndim is non-positive.");

  SLIC_CHECK_MSG(type != NO_TYPE_ID,
                 SIDRE_VIEW_LOG_PREPEND << "Could not allocate: Data type was 'NO_TYPE_ID'.");

  SLIC_CHECK_MSG(shape != nullptr,
                 SIDRE_VIEW_LOG_PREPEND << "Could not allocate: specified shape is nullptr.");

  IndexType num_elems = 1;
  if(shape != nullptr)
  {
    for(int d = 0; d < ndims; ++d)
    {
      SLIC_CHECK_MSG(shape[d] > 0,
                     SIDRE_VIEW_LOG_PREPEND << "Could not allocate: shape is non-positive.");
      num_elems *= shape[d];
      if(num_elems <= 0)
      {
        break;
      }
    }
  }

  if(ndims > 0 && shape != nullptr && num_elems > 0 && type != NO_TYPE_ID)
  {
    describe(type, ndims, shape);
    allocID = getValidArrayAllocatorId(allocID);
    allocate(allocID);
  }

  return this;
}

/*
 *************************************************************************
 *
 * Allocate array data for view described by a Conduit data type object.
 *
 *************************************************************************
 */
View* View::allocate(const DataType& dtype, int allocID)
{
  if(dtype.is_empty())
  {
    SLIC_CHECK_MSG(!dtype.is_empty(),
                   SIDRE_VIEW_LOG_PREPEND << "Unable to allocate View with empty data type.");
    return this;
  }

  describe(dtype);
  allocID = getValidArrayAllocatorId(allocID);
  allocate(allocID);

  return this;
}

/*
 *************************************************************************
 *
 * Reallocate data for view to given number of elements.
 * This function requires that the view is already described.
 * The state may transition from EMPTY to BUFFER;
 * otherwise, the state must already be BUFFER.
 *
 *************************************************************************
 */
View* View::reallocate(IndexType num_elems)
{
  TypeID vtype = static_cast<TypeID>(m_schema.dtype().id());

  if(num_elems < 0)
  {
    SLIC_CHECK_MSG(false, SIDRE_VIEW_LOG_PREPEND << "Unable to reallocate, num_elems must be >= 0");
  }
  else if(isAllocateValid())
  {
    if(m_state == EMPTY)
    {
      allocate(vtype, num_elems);
    }
    else if(m_data_buffer->isAllocated())  //  XXX if ( isAllocated() )
    {
      describe(vtype, num_elems);
      m_data_buffer->reallocate(num_elems);
      apply();
    }
    else
    {
      allocate(vtype, num_elems);
    }
  }

  return this;
}

/*
 *************************************************************************
 *
 * Deallocate data for view.
 *
 *************************************************************************
 */
View* View::deallocate()
{
  if(!isAllocateValid())
  {
    SLIC_CHECK_MSG(isAllocateValid(),
                   SIDRE_VIEW_LOG_PREPEND << "View's state " << getStateStringName(m_state)
                                          << " does not allow data deallocation");
    return this;
  }

  if(hasBuffer())
  {
    m_data_buffer->deallocate();
  }

  return this;
}

/*
 *************************************************************************
 *
 * Reallocate data for view using a Conduit data type object.
 *
 *************************************************************************
 */
View* View::reallocate(const DataType& dtype)
{
  // If we don't have an allocated buffer, we can just call allocate.
  if(!isAllocated())
  {
    return allocate(dtype);
  }

  TypeID type = static_cast<TypeID>(dtype.id());
  TypeID view_type = static_cast<TypeID>(m_schema.dtype().id());

  if(dtype.is_empty() || !isAllocateValid() || type != view_type)
  {
    SLIC_CHECK_MSG(!dtype.is_empty(),
                   SIDRE_VIEW_LOG_PREPEND << "Unable to re-allocate View with empty data type.");
    SLIC_CHECK_MSG(isAllocateValid(),
                   SIDRE_VIEW_LOG_PREPEND << "View's state " << getStateStringName(m_state)
                                          << " does not allow data re-allocation");
    SLIC_CHECK_MSG(type == view_type,
                   SIDRE_VIEW_LOG_PREPEND << "Attempting to re-allocate view with different type.");
    return this;
  }

  describe(dtype);
  IndexType num_elems = dtype.number_of_elements();
  m_data_buffer->reallocate(num_elems);
  apply();

  return this;
}

/*
 *************************************************************************
 *
 * Reshape an array View.
 *
 *************************************************************************
 */
View* View::reshapeArray(int ndims, const IndexType* shape)
{
  if(m_state != BUFFER && m_state != EXTERNAL)
  {
    SLIC_WARNING(SIDRE_VIEW_LOG_PREPEND
                 << "View can only reshape array states (BUFFER or EXTERNAL).");
    return this;
  }

  if(shape == nullptr)
  {
    SLIC_WARNING(SIDRE_VIEW_LOG_PREPEND << "Invalid shape.");
    return this;
  }

  IndexType newSize = shape[0];
  for(int d = 1; d < ndims; ++d)
  {
    newSize *= shape[d];
  }
  if(newSize != getNumElements())
  {
    SLIC_WARNING(SIDRE_VIEW_LOG_PREPEND << "View reshape must not change the number of elements."
                                        << " Current num elements: " << getNumElements()
                                        << "; Desired size: " << newSize);
    return this;
  }

  // If View was applied before reshape, then reapply it.
  const bool is_applied = m_is_applied;
  describe(getTypeID(), ndims, shape);
  if(is_applied)
  {
    apply();
  }

  return this;
}

/*
 *************************************************************************
 *
 * Attach/detach buffer to view.
 *
 *************************************************************************
 */
View* View::attachBuffer(Buffer* buff)
{
  if(m_state == BUFFER && buff == nullptr)
  {
    Buffer* old_buffer = detachBuffer();
    if(old_buffer->getNumViews() == 0)
    {
      getOwningGroup()->getDataStore()->destroyBuffer(old_buffer);
    }
    unapply();
  }
  else if(m_state == EMPTY && buff != nullptr)
  {
    m_data_buffer = buff;
    buff->attachToView(this);
    m_state = BUFFER;
    SLIC_ASSERT(m_is_applied == false);

    // If view is described and the buffer is allocated, then call apply.
    if(isDescribed() && m_data_buffer->isAllocated())
    {
      apply();
    }
  }

  return this;
}

/*
 *************************************************************************
 *
 * Detach buffer from view.
 *
 *************************************************************************
 */
Buffer* View::detachBuffer()
{
  Buffer* buff = nullptr;

  if(m_state == BUFFER)
  {
    buff = m_data_buffer;
    m_data_buffer->detachFromView(this);
  }

  return buff;
}

/*
 *************************************************************************
 *
 * Clear data and metadata from a View.
 *
 *************************************************************************
 */
void View::clear()
{
  switch(m_state)
  {
  case EMPTY:
    undescribe();
    break;
  case BUFFER:
    attachBuffer(nullptr);
    undescribe();
    break;
  case EXTERNAL:
    setExternalDataPtr(nullptr);
    undescribe();
    break;
  case STRING:
  case TUPLE:
    unapply();
    undescribe();
    break;
  default:
    SLIC_ASSERT_MSG(
      false,
      SIDRE_VIEW_LOG_PREPEND << "View is in unexpected state: " << getStateStringName(m_state));
  }

  m_attr_values.clear();
  m_state = EMPTY;
}

/*
 *************************************************************************
 *
 * Apply data description to data.
 *
 *************************************************************************
 */
View* View::apply()
{
  if(!isApplyValid())
  {
    SLIC_CHECK_MSG(isApplyValid(),
                   SIDRE_VIEW_LOG_PREPEND << "View's state, '" << getStateStringName(m_state)
                                          << "', does not allow apply operation");
    return this;
  }

  void* data_pointer = nullptr;

  if(hasBuffer())
  {
    data_pointer = m_data_buffer->getVoidPtr();
  }
  else
  {
    SLIC_ASSERT(m_state == EXTERNAL);
    data_pointer = m_external_ptr;
  }

  m_node.set_external(m_schema, data_pointer);
  m_is_applied = true;

  return this;
}

/*
 *************************************************************************
 *
 * Apply given # elems, offset, stride description to data view.
 *
 *************************************************************************
 */
View* View::apply(IndexType num_elems, IndexType offset, IndexType stride)
{
  if(num_elems < 0)
  {
    SLIC_CHECK_MSG(num_elems >= 0,
                   SIDRE_VIEW_LOG_PREPEND << "Could not apply -- num_elems was less than zero.");
    return this;
  }

  DataType dtype(m_schema.dtype());
  if(dtype.is_empty())
  {
    dtype = conduit::DataType::default_dtype(m_data_buffer->getTypeID());
  }

  const size_t bytes_per_elem = dtype.element_bytes();

  dtype.set_number_of_elements(num_elems);
  dtype.set_offset(offset * bytes_per_elem);
  dtype.set_stride(stride * bytes_per_elem);

  describe(dtype);

  apply();

  return this;
}

/*
 *************************************************************************
 *
 * Apply given type, # elems, offset, stride description to data view.
 *
 *************************************************************************
 */
View* View::apply(TypeID type, IndexType num_elems, IndexType offset, IndexType stride)
{
  if(type == NO_TYPE_ID || num_elems < 0)
  {
    SLIC_CHECK_MSG(type != NO_TYPE_ID, SIDRE_VIEW_LOG_PREPEND << "Could not apply -- invalid type.");

    SLIC_CHECK_MSG(num_elems >= 0,
                   SIDRE_VIEW_LOG_PREPEND << "Could not apply -- num_elems was less than zero.");

    return this;
  }

  DataType dtype = conduit::DataType::default_dtype(type);

  const size_t bytes_per_elem = dtype.element_bytes();

  dtype.set_number_of_elements(num_elems);
  dtype.set_offset(offset * bytes_per_elem);
  dtype.set_stride(stride * bytes_per_elem);

  describe(dtype);
  apply();

  return this;
}

/*
 *************************************************************************
 *
 * Apply given type, number of dimensions and shape to data view.
 * If ndims is 1 then do not save in m_shape.
 *
 *************************************************************************
 */
View* View::apply(TypeID type, int ndims, const IndexType* shape)
{
  if(type == NO_TYPE_ID || ndims < 1 || shape == nullptr)
  {
    SLIC_CHECK_MSG(type != NO_TYPE_ID, SIDRE_VIEW_LOG_PREPEND << "Could not apply -- invalid type.");
    SLIC_CHECK_MSG(ndims >= 1,
                   SIDRE_VIEW_LOG_PREPEND << "Could not apply -- ndims was less than one.");
    SLIC_CHECK_MSG(shape != nullptr, SIDRE_VIEW_LOG_PREPEND << "Could not apply -- shape was null.");

    return this;
  }

  describe(type, ndims, shape);
  apply();

  return this;
}

/*
 *************************************************************************
 *
 * Apply a data type description to data view.
 *
 *************************************************************************
 */
View* View::apply(const DataType& dtype)
{
  if(dtype.is_empty())
  {
    SLIC_CHECK_MSG(!dtype.is_empty(),
                   SIDRE_VIEW_LOG_PREPEND << "Unable to apply description, data type is empty.");
    return this;
  }

  describe(dtype);
  apply();

  return this;
}

/*
 *************************************************************************
 *
 * Get void pointer to any data held by the view.
 *
 *************************************************************************
 */
void* View::getVoidPtr() const
{
  void* rv = nullptr;

  switch(m_state)
  {
  case EMPTY:
    break;
  case EXTERNAL:
    if(isApplied())
    {
      rv = const_cast<void*>(m_node.data_ptr());
    }
    else
    {
      rv = m_external_ptr;  // Opaque
    }
    break;
  case BUFFER:
    if(isApplied())
    {
      rv = const_cast<void*>(m_node.data_ptr());
    }
    else
    {
      SLIC_CHECK_MSG(false, SIDRE_VIEW_LOG_PREPEND << "View has no applied data.");
    }
    break;
  case STRING:
  case TUPLE:
    rv = const_cast<void*>(m_node.data_ptr());
    break;
  default:
    SLIC_ASSERT_MSG(
      false,
      SIDRE_VIEW_LOG_PREPEND << "View is in unexpected state: " << getStateStringName(m_state));
  }

  return rv;
}

/*
 *************************************************************************
 *
 * Set data view to hold external data.
 *
 *************************************************************************
 */
View* View::setExternalDataPtr(void* external_ptr)
{
  if(m_state == EMPTY || m_state == EXTERNAL)
  {
    if(external_ptr == nullptr)
    {
      unapply();
      m_external_ptr = nullptr;
      m_state = EMPTY;
    }
    else
    {
      m_external_ptr = external_ptr;
      m_state = EXTERNAL;

      if(isDescribed())
      {
        apply();
      }
    }
  }
  else
  {
    SLIC_CHECK_MSG(m_state == EMPTY || m_state == EXTERNAL,
                   SIDRE_VIEW_LOG_PREPEND << "Calling setExternalDataPtr on View with "
                                          << getStateStringName(m_state) << " data is not allowed.");
  }

  return this;
}

/*
 *************************************************************************
 *
 * Update the data in this View with the data from other
 *
 *************************************************************************
 */
View* View::updateFrom(const View* other)
{
  if(!isUpdateableFrom(other))
  {
    SLIC_WARNING(SIDRE_VIEW_LOG_PREPEND << "View '" << getPathName() << "' is not updateable "
                                        << "from View '" << other->getPathName() << "'");
    return this;
  }

  SLIC_WARNING_IF(getTypeID() != other->getTypeID(),
                  SIDRE_VIEW_LOG_PREPEND << "Updating View " << getPathName() << " with type "
                                         << getTypeID() << " from View " << other->getPathName()
                                         << " with type " << other->getTypeID());

  char* dst = static_cast<char*>(getVoidPtr());
  dst += getOffset() * getBytesPerElement();

  char* src = static_cast<char*>(other->getVoidPtr());
  src += other->getOffset() * other->getBytesPerElement();

  copy(dst, src, getTotalBytes());

  return this;
}

/*
 *************************************************************************
 *
 * Return true if view contains allocated data. This could mean a buffer
 * with allocated data, or a scalar value, or a string.
 *
 * Note: Most of our isXXX functions are implemented in the header.
 * This one is in not, because we are only forward declaring the buffer
 * class in the view header.
 *************************************************************************
 */
bool View::isAllocated() const
{
  bool rv = false;

  switch(m_state)
  {
  case EMPTY:
    break;
  case BUFFER:
    // false if buffer is not allocated or view is not described
    rv = isDescribed() && m_data_buffer->isAllocated();
    break;
  case EXTERNAL:
  case STRING:
  case TUPLE:
    rv = true;
    break;
  default:
    SLIC_ASSERT_MSG(
      false,
      SIDRE_VIEW_LOG_PREPEND << "View is in unexpected state: " << getStateStringName(m_state));
  }

  return rv;
}

/*
 *************************************************************************
 *
 * Return number of dimensions and fill in shape information.
 *
 *************************************************************************
 */
int View::getShape(int ndims, IndexType* shape) const
{
  const int viewDims = getNumDimensions();
  if(ndims < viewDims)
  {
    SLIC_WARNING(SIDRE_VIEW_LOG_PREPEND << "getShape(): passed in array is not large enough for "
                                           "view's shape. Array only has space for "
                                        << ndims << " dims, but view has " << viewDims);
    return -1;
  }

  for(int i = 0; i < viewDims; ++i)
  {
    shape[i] = m_shape[i];
  }

  // Fill the rest of the array with zeros (when ndims > viewDims)
  if(ndims > viewDims)
  {
    for(int i = viewDims; i < ndims; ++i)
    {
      shape[i] = 0;
    }
  }

  return viewDims;
}

/*
 *************************************************************************
 *
 * Return offset from description in terms of number of elements (0 if not
 * described)
 *
 *************************************************************************
 */
IndexType View::getOffset() const
{
  IndexType offset = 0;

  if(isDescribed())
  {
    offset = m_schema.dtype().offset();

    const auto bytes_per_elem = getBytesPerElement();
    if(bytes_per_elem != 0)
    {
      SLIC_ERROR_IF(offset % bytes_per_elem != 0,
                    SIDRE_VIEW_LOG_PREPEND
                      << "Error calculating offset. "
                      << "Sidre assumes that offsets are given as integral number "
                      << "of elements into the array.  In this View, the offset was " << offset
                      << " bytes and each element is " << bytes_per_elem
                      << " bytes. If you have a need for "
                      << "non-integral offsets, please contact the Sidre team");

      offset /= bytes_per_elem;
    }
  }

  return offset;
}

/*
 *************************************************************************
 *
 * Return stride from description in terms of number of elements (1 if not
 * described)
 *
 *************************************************************************
 */
IndexType View::getStride() const
{
  IndexType stride = 1;

  if(isDescribed())
  {
    stride = m_schema.dtype().stride();

    const auto bytes_per_elem = getBytesPerElement();
    if(bytes_per_elem != 0)
    {
      SLIC_ERROR_IF(stride % bytes_per_elem != 0,
                    SIDRE_VIEW_LOG_PREPEND
                      << "Error caclulating stride. "
                      << "Sidre assumes that strides are given as integral number "
                      << "of elements into the array. In this View, the stride was " << stride
                      << " bytes and each element is " << bytes_per_elem
                      << " bytes. If you have a need for "
                      << "non-integral strides, please contact the Sidre team");

      stride /= bytes_per_elem;
    }
  }

  return stride;
}

/*
 *************************************************************************
 *
 * Test equivalence of two Views
 *
 *************************************************************************
 */
bool View::isEquivalentTo(const View* other) const
{
  //add isAllocated() if it can be declared const
  return (getName() == other->getName()) && (getTypeID() == other->getTypeID()) &&
    (isApplied() == other->isApplied()) && (hasBuffer() == other->hasBuffer()) &&
    (getTotalBytes() == other->getTotalBytes());
}

/*
 *************************************************************************
 *
 * Test whether the two Views can update each other
 *
 *************************************************************************
 */
bool View::isUpdateableFrom(const View* other) const
{
  const bool valid_state = (m_state == BUFFER) || (m_state == EXTERNAL);
  const bool other_valid_state = (other->m_state == BUFFER) || (other->m_state == EXTERNAL);
  const bool same_length = (getTotalBytes() == other->getTotalBytes());
  const bool unit_stride = (getStride() == 1) && (other->getStride() == 1);

  return valid_state && other_valid_state && same_length && unit_stride;
}

/*
 *************************************************************************
 *
 * Print JSON description of data view to stdout.
 *
 *************************************************************************
 */
void View::print() const { print(std::cout); }

/*
 *************************************************************************
 *
 * Print JSON description of data view to an ostream.
 *
 *************************************************************************
 */
void View::print(std::ostream& os) const
{
  Node n;
  copyToConduitNode(n);
  n.to_json_stream(os);
}

/*
  Print data in a way that won't crash when data is not host-accessible.

  If data is not host-accessible, print the pointer and a comment.
*/
void View::hostPrint(std::ostream& os) const
{
  if(isString())
  {
    if(isHostAccessible())
    {
      os << ' ' << '"' << getString() << '"';
    }
    else
    {
      os << ' ' << getVoidPtr() << " # non-host string data";
    }
  }
  else if(isScalar())
  {
    switch(getTypeID())
    {
    case detail::SidreTT<axom::float32>::id:
      hostPrintScalar<axom::float32>();
      break;
    case detail::SidreTT<axom::float64>::id:
      hostPrintScalar<axom::float64>();
      break;
    case detail::SidreTT<std::int8_t>::id:
      hostPrintScalar<std::int8_t>();
      break;
    case detail::SidreTT<std::int16_t>::id:
      hostPrintScalar<std::int16_t>();
      break;
    case detail::SidreTT<std::int32_t>::id:
      hostPrintScalar<std::int32_t>();
      break;
    case detail::SidreTT<std::int64_t>::id:
      hostPrintScalar<std::int64_t>();
      break;
    case detail::SidreTT<std::uint8_t>::id:
      hostPrintScalar<std::uint8_t>();
      break;
    case detail::SidreTT<std::uint16_t>::id:
      hostPrintScalar<std::uint16_t>();
      break;
    case detail::SidreTT<std::uint32_t>::id:
      hostPrintScalar<std::int32_t>();
      break;
    case detail::SidreTT<std::uint64_t>::id:
      hostPrintScalar<std::int64_t>();
      break;
    default:
      os << ' ' << getVoidPtr() << " # non-host unknown scalar data";
    }
  }
  else if(hasBuffer() || (isExternal() && !isOpaque()))
  {
    switch(getTypeID())
    {
    case detail::SidreTT<axom::float32>::id:
      hostPrintArray<axom::float32>();
      break;
    case detail::SidreTT<axom::float64>::id:
      hostPrintArray<axom::float64>();
      break;
    case detail::SidreTT<std::int8_t>::id:
      hostPrintArray<std::int8_t>();
      break;
    case detail::SidreTT<std::int16_t>::id:
      hostPrintArray<std::int16_t>();
      break;
    case detail::SidreTT<std::int32_t>::id:
      hostPrintArray<std::int32_t>();
      break;
    case detail::SidreTT<std::int64_t>::id:
      hostPrintArray<std::int64_t>();
      break;
    case detail::SidreTT<std::uint8_t>::id:
      hostPrintArray<std::uint8_t>();
      break;
    case detail::SidreTT<std::uint16_t>::id:
      hostPrintArray<std::uint16_t>();
      break;
    case detail::SidreTT<std::uint32_t>::id:
      hostPrintArray<std::int32_t>();
      break;
    case detail::SidreTT<std::uint64_t>::id:
      hostPrintArray<std::int64_t>();
      break;
    default:
      os << ' ' << getVoidPtr() << " # " << (isHostAccessible() ? "" : "non-")
         << "host data of unrecognized type id " << getTypeID();
    }
  }
  else if(isOpaque())
  {
    if(isHostAccessible())
    {
      os << ' ' << getVoidPtr() << " # opaque host data";
    }
    else
    {
      os << ' ' << getVoidPtr() << " # opaque non-host data";
    }
  }
}

/*
 *************************************************************************
 *
 * Copy data view description to given Conduit node.
 *
 *************************************************************************
 */
void View::copyToConduitNode(Node& n) const
{
  n["name"] = m_name;
  n["schema"] = m_schema.to_json();
  n["value"] = m_node.to_json();
  n["state"] = getStateStringName(m_state);
  n["is_applied"] = m_is_applied;
}

/*
 *************************************************************************
 *
 * Copy data view native layout to given Conduit node.
 *
 *************************************************************************
 */
void View::createNativeLayout(Node& n) const
{
  // see ATK-726 - Handle undescribed and unallocated views in Sidre's
  // createNativeLayout()
  // TODO: Need to handle cases where the view is not described
  // TODO: Need to handle cases where the view is not allocated
  // TODO: Need to handle cases where the view is not applied

  // Note: We are using conduit's pointer rather than the View pointer
  //    since the conduit pointer handles offsetting
  // Note: const_cast the pointer to satisfy conduit's interface
  void* data_ptr = const_cast<void*>(m_node.data_ptr());
  n.set_external(m_node.schema(), data_ptr);
}

/*
 *************************************************************************
 *
 * Deep-copy data to given Conduit node.
 *
 *************************************************************************
 */
void View::deepCopyToConduit(Node& dst, int allocId) const
{
  // see ATK-726 - Handle undescribed and unallocated views in Sidre's
  // createNativeLayout()
  // TODO: Need to handle cases where the view is not described
  // TODO: Need to handle cases where the view is not allocated
  // TODO: Need to handle cases where the view is not applied

  const conduit::DataType& srcDtype = m_node.dtype();
  dst.set(srcDtype);
  if(allocId != INVALID_ALLOCATOR_ID)
  {
    dst.set_allocator(ConduitMemory::axomAllocIdToConduit(allocId));
  }
  if(isAllocated())
  {
    // Using set_node to set dst: would reset dst's allocator id to
    // the Conduit default (not what we want) if dst is an object or
    // list.  Fortunately, Sidre never uses a Conduit node in those
    // modes.
#ifdef AXOM_DEBUG
    const auto oldAllocatorId = dst.allocator();
#endif
    dst.set_node(m_node);
    SLIC_ASSERT(dst.allocator() == oldAllocatorId);
  }
}

/*
 *************************************************************************
 *
 * Copy the metadata for the View.
 *
 *************************************************************************
 */
void View::copyMetadataToNode(Node& n) const
{
  n["state"] = getStateStringName(m_state);
  n["schema"] = m_schema.to_json();
  n["is_applied"] = m_is_applied;
}

/*
 *************************************************************************
 *
 * PRIVATE ctor for View not associated with any data.
 *
 *************************************************************************
 */
View::View(const std::string& name)
  : m_name(name)
  , m_index(InvalidIndex)
  , m_owning_group(nullptr)
  , m_data_buffer(nullptr)
  , m_schema()
  , m_node()
  , m_shape()
  , m_external_ptr(nullptr)
  , m_state(EMPTY)
  , m_is_applied(false)
{ }

/*
 *************************************************************************
 *
 * PRIVATE dtor.
 *
 *************************************************************************
 */
View::~View()
{
  if(m_data_buffer != nullptr)
  {
    m_data_buffer->detachFromView(this);
  }
}

/*
 *************************************************************************
 *
 * PRIVATE method to describe data view with type and number of elements.
 *         Caller has already checked arguments.
 *
 *************************************************************************
 */
void View::describe(TypeID type, IndexType num_elems)
{
  DataType dtype = conduit::DataType::default_dtype(type);
  dtype.set_number_of_elements(num_elems);
  m_schema.set(dtype);
  describeShape();
  m_is_applied = false;
}

/*
 *************************************************************************
 *
 * PRIVATE method to describe data view with type, number of dimensions,
 *         and number of elements per dimension.
 *         Caller has already checked arguments.
 *
 *************************************************************************
 */
void View::describe(TypeID type, int ndims, const IndexType* shape)
{
  SLIC_CHECK_MSG(shape != nullptr,
                 SIDRE_VIEW_LOG_PREPEND << "Could not allocate: specified shape is nullptr.");
  if(shape == nullptr)
  {
    return;
  }

  IndexType num_elems = 0;
  if(ndims > 0)
  {
    num_elems = shape[0];
    for(int i = 1; i < ndims; i++)
    {
      num_elems *= shape[i];
    }
  }

  describe(type, num_elems);
  describeShape(ndims, shape);
}

/*
 *************************************************************************
 *
 * PRIVATE method to describe data view with a Conduit data type object.
 *         Caller has already checked arguments.
 *
 *************************************************************************
 */
void View::describe(const DataType& dtype)
{
  m_schema.set(dtype);
  describeShape();
  m_is_applied = false;
}

/*
 *************************************************************************
 *
 * PRIVATE method set shape to described length.
 * This is called after describe to set the shape.
 *
 *************************************************************************
 */
void View::describeShape()
{
  m_shape.clear();
  m_shape.push_back(m_schema.dtype().number_of_elements());
}

/*
 *************************************************************************
 *
 * PRIVATE method set shape from user input.
 *
 *************************************************************************
 */
void View::describeShape(int ndims, const IndexType* shape)
{
  SLIC_CHECK_MSG(shape != nullptr,
                 SIDRE_VIEW_LOG_PREPEND << "Could not allocate: specified shape is nullptr.");
  if(shape == nullptr)
  {
    return;
  }

  m_shape.clear();
  for(int i = 0; i < ndims; i++)
  {
    m_shape.push_back(shape[i]);
  }
}

/*
 *************************************************************************
 *
 * PRIVATE method copy the contents of this into a undescribed EMPTY view.
 *
 *************************************************************************
 */
void View::copyView(View* copy) const
{
  SLIC_ASSERT_MSG(copy->m_state == EMPTY && !copy->isDescribed(),
                  SIDRE_VIEW_LOG_PREPEND << "copyView can only copy into undescribed view "
                                         << "with empty state.");

  if(isDescribed())
  {
    copy->describe(m_schema.dtype());
  }

  switch(m_state)
  {
  case EMPTY:
    // Nothing more to do
    break;
  case STRING:
  case TUPLE:
    copy->m_node = m_node;
    copy->m_state = m_state;
    copy->m_is_applied = true;
    break;
  case EXTERNAL:
    copy->setExternalDataPtr(m_external_ptr);
    break;
  case BUFFER:
    copy->attachBuffer(m_data_buffer);
    break;
  default:
    SLIC_ASSERT_MSG(
      false,
      SIDRE_VIEW_LOG_PREPEND << "View is in unexpected state: " << getStateStringName(m_state));
  }
}

/*
 *************************************************************************
 *
 * PRIVATE method copy the contents of this into a undescribed EMPTY view.
 *
 *************************************************************************
 */
void View::deepCopyView(View* copy, int arrayAllocId, int tupleAllocId) const
{
  SLIC_ASSERT_MSG(copy->m_state == EMPTY && !copy->isDescribed(),
                  SIDRE_VIEW_LOG_PREPEND << "deepCopyView can only copy into undescribed view "
                                         << "with empty state.");

  if(isDescribed())
  {
    if(hasBuffer() || m_state == EXTERNAL)
    {
      copy->allocate(getTypeID(), getNumDimensions(), m_shape.data(), arrayAllocId);
    }
    else
    {
      copy->describe(getTypeID(), getNumDimensions(), m_shape.data());
    }
  }

  switch(m_state)
  {
  case EMPTY:
    // Nothing more to do
    break;
  case STRING:
  case TUPLE:
    copy->m_node.set_allocator(ConduitMemory::axomAllocIdToConduit(tupleAllocId));
    copy->m_node.set_node(m_node);
    copy->m_state = m_state;
    copy->m_is_applied = true;
    break;
  case EXTERNAL:
    if(!copy->isAllocated())
    {
      copy->allocate(arrayAllocId);
    }
    if(isApplied())
    {
      copy->apply();
    }
    {
      IndexType stride = getStride();
      IndexType src_offset = getOffset();
      IndexType dst_offset = copy->getOffset();
      IndexType num_bytes = getBytesPerElement();
      IndexType j = 0;
      for(IndexType i = 0; i < getNumElements(); ++i)
      {
        char* copy_dst = static_cast<char*>(copy->getVoidPtr()) + (dst_offset + i) * num_bytes;
        const char* copy_src = static_cast<const char*>(getVoidPtr()) + (src_offset + j) * num_bytes;
        axom::copy(copy_dst, copy_src, num_bytes);

        j += stride;
      }
    }
    break;
  case BUFFER:
    if(isAllocated() && !copy->isAllocated())
    {
      copy->allocate(arrayAllocId);
    }
    if(isApplied())
    {
      copy->apply();
    }
    if(hasBuffer())
    {
      IndexType stride = getStride();
      IndexType src_offset = getOffset();
      IndexType dst_offset = copy->getOffset();
      IndexType num_bytes = getBytesPerElement();
      IndexType j = 0;
      for(IndexType i = 0; i < getNumElements(); ++i)
      {
        char* copy_dst = static_cast<char*>(copy->getVoidPtr()) + (dst_offset + i) * num_bytes;
        const char* copy_src = static_cast<const char*>(getVoidPtr()) + (src_offset + j) * num_bytes;
        axom::copy(copy_dst, copy_src, num_bytes);

        j += stride;
      }
    }
    break;
  default:
    SLIC_ASSERT_MSG(
      false,
      SIDRE_VIEW_LOG_PREPEND << "View is in unexpected state: " << getStateStringName(m_state));
  }
}

/*
 *************************************************************************
 *
 * PRIVATE method returns true if view can allocate data; else false.
 *
 * This method does not need to emit the view state as part of it's
 * checking.  The caller functions are already printing out the view
 * state if this function returns false.
 *************************************************************************
 */
bool View::isAllocateValid() const
{
  bool rv = false;

  switch(m_state)
  {
  case EMPTY:
    rv = isDescribed();
    break;
  case STRING:
  case TUPLE:
  case EXTERNAL:
    SLIC_CHECK_MSG(false,
                   SIDRE_VIEW_LOG_PREPEND << "Allocate is not valid for view in '"
                                          << getStateStringName(m_state) << "' state.");
    break;
  case BUFFER:
    rv = isDescribed() && m_data_buffer->getNumViews() == 1;
    break;
  default:
    SLIC_ASSERT_MSG(
      false,
      SIDRE_VIEW_LOG_PREPEND << "View is in unexpected state: " << getStateStringName(m_state));
  }

  return rv;
}

/*
 *************************************************************************
 *
 * PRIVATE method returns true if apply is a valid operation on view;
 * else false.
 *
 * For an EXTERNAL view, assume user provided m_external_ptr and
 * description are consistent. This includes m_external_ptr == NULL.
 *
 *************************************************************************
 */
bool View::isApplyValid() const
{
  bool rv = false;

  if(!isDescribed())
  {
    SLIC_CHECK_MSG(
      false,
      SIDRE_VIEW_LOG_PREPEND << "Apply is not valid. View does not have a description.");
    return rv;
  }

  switch(m_state)
  {
  case EMPTY:
  case STRING:
  case TUPLE:
    SLIC_CHECK_MSG(false,
                   SIDRE_VIEW_LOG_PREPEND << "Apply is not valid for View with state '"
                                          << getStateStringName(m_state) << "'.'");
    break;
  case EXTERNAL:
    SLIC_ASSERT(m_external_ptr != nullptr);
    rv = isDescribed();
    break;
  case BUFFER:
    rv = 0 <= getTotalBytes() && getTotalBytes() <= m_data_buffer->getTotalBytes();
    SLIC_CHECK_MSG(0 <= getTotalBytes(),
                   SIDRE_VIEW_LOG_PREPEND << "Apply is not valid on data with zero length.");
    SLIC_CHECK_MSG(getTotalBytes() <= m_data_buffer->getTotalBytes(),
                   SIDRE_VIEW_LOG_PREPEND << "Apply is not valid. "
                                          << "View's datatype length exceeds bytes in buffer.");
    break;
  default:
    SLIC_ASSERT_MSG(
      false,
      SIDRE_VIEW_LOG_PREPEND << "View is in unexpected state: " << getStateStringName(m_state));
  }

  return rv;
}

/*
 *************************************************************************
 *
 * PRIVATE method returns string name of given view state enum value.
 *
 *************************************************************************
 */
char const* View::getStateStringName(State state)
{
  char const* ret_string = nullptr;

  switch(state)
  {
  case EMPTY:
    ret_string = "EMPTY";
    break;
  case BUFFER:
    ret_string = "BUFFER";
    break;
  case EXTERNAL:
    ret_string = "EXTERNAL";
    break;
  case TUPLE:
    ret_string = "TUPLE";
    break;
  case STRING:
    ret_string = "STRING";
    break;
  default:
    ret_string = "UNKNOWN";
  }

  return ret_string;
}

/*
 *************************************************************************
 *
 * PRIVATE method returns state enum value when given string with a
 * state name.
 *
 *************************************************************************
 */
View::State View::getStateId(const std::string& name) const
{
  State res = EMPTY;
  if(name == "EMPTY")
  {
    res = EMPTY;
  }
  else if(name == "BUFFER")
  {
    res = BUFFER;
  }
  else if(name == "EXTERNAL")
  {
    res = EXTERNAL;
  }
  else if(name == "TUPLE")
  {
    res = TUPLE;
  }
  else if(name == "STRING")
  {
    res = STRING;
  }
  else if(name == "UNKNOWN")
  {
    res = EMPTY;
  }

  return res;
}

/*
 * Return whether view data is accessible on the host CPU.
 */
bool View::isHostAccessible() const
{
  bool rval = false;
  void* dataPtr = getVoidPtr();
  if(dataPtr != nullptr)
  {
    int allocId = axom::getAllocatorIDFromPointer(dataPtr);
    rval = axom::execution_space<axom::SEQ_EXEC>::usesAllocId(allocId);
  }
  return rval;
}

/*
 *************************************************************************
 *
 * PRIVATE method to copy view data to given Conduit node using
 * given set of ids to maintain correct association of data buffers
 * to data views.
 *
 *************************************************************************
 */
void View::exportTo(conduit::Node& data_holder, std::set<IndexType>& buffer_indices) const
{
  data_holder["state"] = getStateStringName(m_state);
  exportAttribute(data_holder);

  switch(m_state)
  {
  case EMPTY:
    if(isDescribed())
    {
      exportDescription(data_holder);
    }
    break;
  case BUFFER:
  {
    IndexType buffer_id = getBuffer()->getIndex();
    data_holder["buffer_id"] = buffer_id;
    if(isDescribed())
    {
      exportDescription(data_holder);
    }
    data_holder["is_applied"] = static_cast<unsigned char>(m_is_applied);
    buffer_indices.insert(buffer_id);
    break;
  }
  case EXTERNAL:
    if(isDescribed())
    {
      exportDescription(data_holder);
    }
    else
    {
      // If there is no description, make it an EMPTY view
      data_holder["state"] = getStateStringName(EMPTY);
    }
    break;
  case TUPLE:
  case STRING:
    data_holder["value"] = getNode();
    break;
  default:
    SLIC_ASSERT_MSG(
      false,
      SIDRE_VIEW_LOG_PREPEND << "View is in unexpected state: " << getStateStringName(m_state));
  }
}

/*
 *************************************************************************
 * TODO
 *
 *************************************************************************
 */
void View::importFrom(conduit::Node& data_holder, const std::map<IndexType, IndexType>& buffer_id_map)
{
  m_state = getStateId(data_holder["state"].as_string());
  importAttribute(data_holder);

  switch(m_state)
  {
  case EMPTY:
    importDescription(data_holder);
    break;
  case BUFFER:
  {
    // If view has a buffer, the easiest way to restore it is to use a series of
    // API calls.
    // Start from scratch
    m_state = EMPTY;

    IndexType old_buffer_id = data_holder["buffer_id"].to_int64();
    bool is_applied = data_holder["is_applied"].as_unsigned_char();

    SLIC_ASSERT_MSG(buffer_id_map.find(old_buffer_id) != buffer_id_map.end(),
                    SIDRE_VIEW_LOG_PREPEND << "Buffer id map is old."
                                           << "New id entry for buffer " << old_buffer_id);

    Buffer* buffer = m_owning_group->getDataStore()->getBuffer(buffer_id_map.at(old_buffer_id));

    importDescription(data_holder);
    attachBuffer(buffer);
    if(is_applied)
    {
      apply();
    }
    break;
  }
  case EXTERNAL:
    importDescription(data_holder);
    break;
  case TUPLE:
  case STRING:
    m_node = data_holder["value"];
    m_schema.set(m_node.schema());
    m_is_applied = true;
    break;
  default:
    SLIC_ASSERT_MSG(
      false,
      SIDRE_VIEW_LOG_PREPEND << "View is in unexpected state: " << getStateStringName(m_state));
  }
}

/*
 *************************************************************************
 *
 * Import Node holding an array into a View with an attached Buffer.
 *
 *************************************************************************
 */
View* View::importArrayNode(const Node& array)
{
  conduit::DataType array_dtype = array.dtype();

  if(array_dtype.is_number())
  {
    if(m_state == BUFFER)
    {
      setBufferViewToEmpty();
    }
    if(m_state == EMPTY)
    {
      Buffer* buff = m_owning_group->getDataStore()->createBuffer();

      conduit::index_t num_ele = array_dtype.number_of_elements();
      conduit::index_t ele_bytes = DataType::default_bytes(array_dtype.id());

      int allocID = m_owning_group->getDefaultArrayAllocatorID();
      buff->allocate((TypeID)array_dtype.id(), num_ele, allocID);

      // copy the data in a way that matches
      // to compact representation of the buffer
      conduit::uint8* data_ptr = (conduit::uint8*)buff->getVoidPtr();
      for(conduit::index_t i = 0; i < num_ele; i++)
      {
        axom::copy(data_ptr, array.element_ptr(i), ele_bytes);
        data_ptr += ele_bytes;
      }

      attachBuffer(buff);

      // it is important to not use the data type directly
      // it could contain offsets that are no longer
      // valid our new buffer
      apply((TypeID)array_dtype.id(), array_dtype.number_of_elements());
    }
    else
    {
      SLIC_CHECK_MSG(m_state == EMPTY,
                     SIDRE_VIEW_LOG_PREPEND << "Unable to import array Node to View with state: "
                                            << getStateStringName(m_state));
    }
  }
  else
  {
    SLIC_CHECK_MSG(
      array_dtype.is_number(),
      SIDRE_VIEW_LOG_PREPEND << "Unable to import array from Node of type: " << array_dtype.name());
  }

  return this;
}

/*
 *************************************************************************
 *
 * PRIVATE method to save view's description to a conduit tree.
 * The shape information is only written if there is more than
 * one dimension.
 *
 *************************************************************************
 */
void View::exportDescription(conduit::Node& data_holder) const
{
  data_holder["schema"] = m_schema.to_json();
  if(getNumDimensions() > 1)
  {
    data_holder["shape"].set(m_shape);
  }
}

/*
 *************************************************************************
 *
 * PRIVATE method to restore a view's description from a conduit tree.
 *
 *************************************************************************
 */
void View::importDescription(conduit::Node& data_holder)
{
  if(data_holder.has_path("schema"))
  {
    conduit::Schema schema(data_holder["schema"].as_string());
    describe(schema.dtype());
    if(data_holder.has_path("shape"))
    {
      Node& n = data_holder["shape"];
      IndexType* shape = n.value();
      const auto ndims = static_cast<int>(n.dtype().number_of_elements());
      describeShape(ndims, shape);
    }
  }
}

/*
 *************************************************************************
 *
 * PRIVATE method to save view's attributes to a conduit tree.
 * Only add "attribute" Node if the View has any attributes.
 *
 *************************************************************************
 */
void View::exportAttribute(conduit::Node& data_holder) const
{
  IndexType aidx = getFirstValidAttrValueIndex();

  if(aidx == InvalidIndex)
  {
    return;
  }

  Node& node = data_holder["attribute"];
  node.set(DataType::object());

  while(indexIsValid(aidx))
  {
    const Attribute* attr = getAttribute(aidx);

    node[attr->getName()] = getAttributeNodeRef(attr);

    aidx = getNextValidAttrValueIndex(aidx);
  }
}

/*
 *************************************************************************
 *
 * PRIVATE method to restore a view's attributes from a conduit tree.
 *
 *************************************************************************
 */
void View::importAttribute(conduit::Node& data_holder)
{
  if(data_holder.has_path("attribute"))
  {
    conduit::NodeIterator attrs_itr = data_holder["attribute"].children();
    while(attrs_itr.has_next())
    {
      Node& n_attr = attrs_itr.next();
      std::string attr_name = attrs_itr.name();

      Attribute* attr = getAttribute(attr_name);
      if(attr != nullptr)
      {
        m_attr_values.setNode(attr, n_attr);
      }
    }
  }
}

/*
 *************************************************************************
 *
 * Rename this View with a new string name.
 *
 *************************************************************************
 */
bool View::rename(const std::string& new_name)
{
  bool do_rename = true;
  if(new_name != m_name)
  {
    Group* parent = getOwningGroup();
    SLIC_CHECK(parent != nullptr);

    if(new_name.empty())
    {
      SLIC_WARNING(SIDRE_VIEW_LOG_PREPEND << "Cannot rename View to an empty string.");
      do_rename = false;
    }
    else if(new_name.find(parent->getPathDelimiter()) != std::string::npos)
    {
      SLIC_WARNING(SIDRE_VIEW_LOG_PREPEND << "Cannot rename View " << getPathName()
                                          << " to path name '" << new_name
                                          << "'. Only strings without path delimiters can "
                                          << "be passed into the rename method.");
      do_rename = false;
    }
    else if(parent->hasGroup(new_name) || parent->hasView(new_name))
    {
      SLIC_WARNING(SIDRE_VIEW_LOG_PREPEND << "Parent group '" << parent->getPathName()
                                          << "' already has a child object named " << new_name << ". "
                                          << "View " << getPathName() << " will not be renamed.");
      do_rename = false;
    }
    else
    {
      View* detached_view = parent->detachView(m_name);
      SLIC_CHECK(detached_view == this);

      m_name = new_name;

      View* attached_view = parent->attachView(detached_view);
      AXOM_UNUSED_VAR(attached_view);
      SLIC_CHECK(attached_view == this);
    }
  }

  return do_rename;
}

/*
 *************************************************************************
 *
 * Return pointer to Attribute from Attribute index.
 *
 *************************************************************************
 */
Attribute* View::getAttribute(IndexType idx)
{
  Attribute* attr = getOwningGroup()->getDataStore()->getAttribute(idx);
  return attr;
}

/*
 *************************************************************************
 *
 * Return pointer to Attribute from Attribute index.
 *
 *************************************************************************
 */
const Attribute* View::getAttribute(IndexType idx) const
{
  const Attribute* attr = getOwningGroup()->getDataStore()->getAttribute(idx);
  return attr;
}

/*
 *************************************************************************
 *
 * Return pointer to Attribute from Attribute name.
 *
 *************************************************************************
 */
Attribute* View::getAttribute(const std::string& name)
{
  Attribute* attr = getOwningGroup()->getDataStore()->getAttribute(name);
  return attr;
}

/*
 *************************************************************************
 *
 * Return pointer to Attribute from Attribute name.
 *
 *************************************************************************
 */
const Attribute* View::getAttribute(const std::string& name) const
{
  const Attribute* attr = getOwningGroup()->getDataStore()->getAttribute(name);
  return attr;
}

/*
 *************************************************************************
 *
 * Set Attribute for a View from Attribute name.
 *
 *************************************************************************
 */
bool View::setAttributeString(IndexType idx, const std::string& value)
{
  const Attribute* attr = getAttribute(idx);

  if(attr == nullptr)
  {
    return false;
  }

  return m_attr_values.setString(attr, value);
}

/*
 *************************************************************************
 *
 * Set Attribute for a View from Attribute name.
 *
 *************************************************************************
 */
bool View::setAttributeString(const std::string& name, const std::string& value)
{
  const Attribute* attr = getAttribute(name);

  if(attr == nullptr)
  {
    return false;
  }

  return m_attr_values.setString(attr, value);
}

/*
 *************************************************************************
 *
 * Set Attribute for a View from Attribute pointer.
 *
 *************************************************************************
 */
bool View::setAttributeString(const Attribute* attr, const std::string& value)
{
  if(attr == nullptr)
  {
    SLIC_CHECK_MSG(attr != nullptr,
                   SIDRE_VIEW_LOG_PREPEND << "setAttributeString: called with a null Attribute");
    return false;
  }

  return m_attr_values.setString(attr, value);
}

/*
 *************************************************************************
 *
 * Return a string attribute from the Attribute index.
 *
 * If the value has not been explicitly set, return the current default.
 *
 *************************************************************************
 */
const char* View::getAttributeString(IndexType idx) const
{
  const Attribute* attr = getAttribute(idx);

  if(attr == nullptr)
  {
    return nullptr;
  }

  return m_attr_values.getString(attr);
}

/*
 *************************************************************************
 *
 * Return a string attribute from the Attribute name.
 *
 * If the value has not been explicitly set, return the current default.
 *
 *************************************************************************
 */
const char* View::getAttributeString(const std::string& name) const
{
  const Attribute* attr = getAttribute(name);

  if(attr == nullptr)
  {
    return nullptr;
  }

  return m_attr_values.getString(attr);
}

/*
 *************************************************************************
 *
 * Return a string attribute from the Attribute pointer.
 *
 * If the value has not been explicitly set, return the current default.
 *
 *************************************************************************
 */
const char* View::getAttributeString(const Attribute* attr) const
{
  if(attr == nullptr)
  {
    SLIC_CHECK_MSG(attr != nullptr,
                   SIDRE_VIEW_LOG_PREPEND << "getAttributeString: called with a null Attribute");
    return nullptr;
  }

  return m_attr_values.getString(attr);
}

int View::getValidArrayAllocatorId(int allocId)
{
  return getOwningGroup()->getValidArrayAllocatorId(allocId);
}

int View::getValidTupleAllocatorId(int allocId)
{
  return getOwningGroup()->getValidTupleAllocatorId(allocId);
}

int View::getValidAllocatorId(int allocId)
{
  if(allocId != axom::INVALID_ALLOCATOR_ID)
  {
    return allocId;
  }
  if(m_state == BUFFER || m_state == EXTERNAL)
  {
    return getOwningGroup()->getValidArrayAllocatorId(allocId);
  }
  if(m_state == TUPLE || m_state == STRING)
  {
    return getOwningGroup()->getValidTupleAllocatorId(allocId);
  }
  SLIC_ASSERT_MSG(
    false,
    "Axom internal error: Cannot determine semantic valid allocator id");  // Should never get here.
  return axom::INVALID_ALLOCATOR_ID;
}

} /* end namespace sidre */
} /* end namespace axom */
