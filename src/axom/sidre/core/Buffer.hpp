// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file Buffer.hpp
 *
 * \brief   Header file containing definition of Buffer class.
 *
 ******************************************************************************
 */

#ifndef SIDRE_BUFFER_HPP_
#define SIDRE_BUFFER_HPP_

// Standard C++ headers
#include <set>

// Other axom headers
#include "axom/core/memory_management.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"
#include "axom/slic/interface/slic.hpp"

// Sidre project headers
#include "SidreTypes.hpp"

namespace axom
{
namespace sidre
{
class DataStore;
class Group;
class View;

/*!
 * \class Buffer
 *
 * \brief Buffer is a container that describes and holds data in memory.
 *
 * The Buffer class has the following properties:
 *
 *    - Buffer objects can only be created using the DataStore
 *      createBuffer() methods. The Buffer ctor is private.
 *    - A Buffer object has a unique identifier within a DataStore,
 *      which is assigned by the DataStore when the Buffer is created.
 *    - The data owned by a Buffer is unique to that Buffer
 *      object; i.e.,  Buffers do not share their data with other Buffers.
 *    - Typical usage is to describe the data a Buffer will hold and then
 *      allocate it by calling one of the Buffer allocate or
 *      reallocate methods.
 *    - A Buffer object maintains a collection of Views that
 *      refer to its data. These references are created when a Buffer
 *      object is attached to a View.
 */
class Buffer
{
public:
  /*!
   * Friend declarations to constrain usage via controlled access to
   * private members.
   */
  friend class DataStore;
  friend class Group;
  friend class View;

  ///@{
  //!  @name Basic query and accessor methods

  /*!
   * \brief Return the unique index of this Buffer object.
   */
  IndexType getIndex() const { return m_index; }

  /*!
   * \brief Return number of Views this Buffer is attached to.
   */
  IndexType getNumViews() const
  {
    // need error checking for this conversion
    return static_cast<IndexType>(m_views.size());
  }

  ///@}

  ///@{
  //!  @name Methods to query and access Buffer data

  /*!
   * \brief Return void-pointer to data held by Buffer.
   */
  void* getVoidPtr() { return m_node.data_ptr(); }

  /*!
   * \brief Return data held by Buffer (return type is type caller assigns
   *        return value to).
   *
   *        Note that if Buffer is not allocated, an empty Conduit
   *        Node::Value is returned.
   */
  Node::Value getData()
  {
    if(!isAllocated())
    {
      SLIC_CHECK_MSG(isAllocated(), "Buffer data is not allocated.");
      return Node().value();
    }

    return m_node.value();
  }

  /*!
   * \brief Return type of data owned by this Buffer object.
   */
  TypeID getTypeID() const { return static_cast<TypeID>(m_node.dtype().id()); }

  /*!
   * \brief Return total number of data elements (of its type) owned by
   *        this Buffer object.
   */
  IndexType getNumElements() const { return m_node.dtype().number_of_elements(); }

  /*!
   * \brief Return total number of bytes of data owned by this Buffer object.
   */
  IndexType getTotalBytes() const { return m_node.dtype().strided_bytes(); }

  /*!
   * \brief Return the number of bytes per element owned by this Buffer object.
   */
  IndexType getBytesPerElement() const { return m_node.dtype().element_bytes(); }

  /*!
   * \brief Return true if Buffer has been (re)allocated with length >= 0, else
   *  false.
   */
  bool isAllocated() const { return (m_node.data_ptr() != nullptr); }

  /*!
   * \brief Return true if data description exists.
   */
  bool isDescribed() const { return !m_node.dtype().is_empty(); }

  ///@}

  ///@{
  //!  @name Data description and allocation methods

  /*!
   * \brief Describe a Buffer with data given data type and number of elements.
   *
   * To use the Buffer, the data must be allocated by calling allocate().
   *
   * If Buffer is already allocated or given number of elements is < 0,
   * method is a no-op.
   *
   * \return pointer to this Buffer object.
   */
  Buffer* describe(TypeID type, IndexType num_elems);

  /*!
   * \brief Allocate data for a Buffer.
   *
   * If the Buffer is not described or already allocated the method is a no-op.
   *
   * \return pointer to this Buffer object.
   */
  Buffer* allocate(int allocID = INVALID_ALLOCATOR_ID);

  /*!
   * \brief Allocate Buffer with data type and number of elements.
   *
   * This is equivalent to: buff->describe(type, num_elems)->allocate().
   *
   * Method is a no-op under the same conditions as either of those methods.
   *
   * \return pointer to this Buffer object.
   */
  Buffer* allocate(TypeID type, IndexType num_elems, int allocID = INVALID_ALLOCATOR_ID);

  /*!
   * \brief Reallocate data to given number of elements.
   *
   * This is equivalent to: buff->describe(type, num_elems)->allocate()
   * if the Buffer is not allocated.
   *
   * If given number of elements < 0, or the Buffer is not already described
   * with type information, this method is a no-op.
   *
   * \return pointer to this Buffer object.
   */
  Buffer* reallocate(IndexType num_elems);

  /*!
   * \brief Deallocate data in a Buffer.
   *
   * If Buffer is attached to Views, it will remain attached to those Views
   * and the descriptions will remain intact. However, the View descriptions
   * will be 'un-applied' since there is no data to apply them to.
   *
   * If the Buffer is subsequently redescribed and/or re-allocated, the
   * associated Views may need to be re-described. They will need to be
   * re-applied if the Views will be used to access the Buffer data.
   *
   * If the Buffer is not allocated, method is a no-op.
   *
   * \return pointer to this Buffer object.
   */
  Buffer* deallocate();

  ///@}

  /*!
   * \brief Copy given number of bytes of data from src into Buffer.
   *
   * If nbytes < 0 or nbytes > getTotalBytes(), or src is null,
   * method is a no-op.
   *
   * \return pointer to this Buffer object.
   */
  Buffer* copyBytesIntoBuffer(void* src, IndexType nbytes);

  /*!
   * \brief Copy Buffer description to a Conduit node.
   */
  void copyToConduitNode(Node& n) const;

  /*!
   * \brief Print JSON description of Buffer to std::cout.
   */
  void print() const;

  /*!
   * \brief Print JSON description of Buffer to given output stream.
   */
  void print(std::ostream& os) const;

  /*!
   * \brief Exports Buffer's state to a Conduit node.
   */
  void exportTo(conduit::Node& data_holder);

  /*!
   * \brief Exports Buffer's metadata to a Conduit node.
   */
  void exportMetadata(conduit::Node& data_holder);

  /*!
   * \brief Import Buffer's state from a Conduit node.
   */
  void importFrom(conduit::Node& data_holder);

private:
  DISABLE_DEFAULT_CTOR(Buffer);
  DISABLE_MOVE_AND_ASSIGNMENT(Buffer);

  /*!
   *  \brief Private ctor assigns id generated by DataStore (must be
   *         unique among Buffers in DataStore.
   */
  Buffer(IndexType uid);

  /*!
   * \brief Private copy ctor.
   */
  Buffer(const Buffer& source);

  /*!
   * \brief Private dtor.
   */
  ~Buffer();

  /*!
   * \brief Private method to attach Buffer to View.
   *
   * Note: If View's Buffer pointer does not match 'this', method is a no-op.
   */
  void attachToView(View* view);

  /*!
   * \brief Private method to detach Buffer from View.
   *
   * Note: If View's Buffer pointer does not match 'this', method is a no-op.
   */
  void detachFromView(View* view);

  /*!
   * \brief Private method to detach Buffer from all Views it is attached to.
   */
  void detachFromAllViews();

  /*!
   * \brief Private method to allocate num_bytes bytes of data using the given
   * allocator and return a void-pointer to the allocation.
   */
  void* allocateBytes(IndexType num_bytes, int allocID = INVALID_ALLOCATOR_ID);

  /*!
   * \brief Private method to delete data referenced by pointer.
   */
  void releaseBytes(void* ptr);

  /// Buffer's unique index within DataStore object that created it.
  IndexType m_index;

  /// Container of Views attached to this Buffer.
  std::set<View*> m_views;

  /// Conduit Node that holds Buffer data.
  Node m_node;
};

} /* end namespace sidre */
} /* end namespace axom */

#endif /* SIDRE_BUFFER_HPP_ */
