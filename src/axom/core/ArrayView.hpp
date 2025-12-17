// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_ARRAYVIEW_HPP_
#define AXOM_ARRAYVIEW_HPP_

#include "axom/core/memory_management.hpp"  // for memory allocation functions
#include "axom/core/ArrayBase.hpp"
#include "axom/core/ArrayIteratorBase.hpp"

namespace axom
{
// Forward declare the templated classes and operator function(s)
template <typename T, int DIM, MemorySpace SPACE>
class ArrayView;

namespace detail
{
// Static information to pass to ArrayBase
template <typename T, int DIM, MemorySpace SPACE>
struct ArrayTraits<ArrayView<T, DIM, SPACE>>
{
  constexpr static bool is_view = true;
};

}  // namespace detail

/// \name ArrayView to wrap a pointer and provide indexing semantics
/// @{

/*!
 * \class ArrayView
 *
 * \brief Provides a view over a generic array container.
 * 
 * The ArrayView expresses a non-owning relationship over a pointer
 *
 * \tparam T the type of the values to hold.
 * \tparam DIM The dimension of the array.
 *
 */
template <typename T, int DIM = 1, MemorySpace SPACE = MemorySpace::Dynamic>
class ArrayView : public ArrayBase<T, DIM, ArrayView<T, DIM, SPACE>>
{
public:
  using value_type = T;
  static constexpr int dimension = DIM;
  static constexpr MemorySpace space = SPACE;
  using ArrayViewIterator = ArrayIteratorBase<const ArrayView<T, DIM, SPACE>, T>;

  /// \brief Default constructor
  AXOM_HOST_DEVICE ArrayView()
#ifndef AXOM_DEVICE_CODE
    : m_allocator_id(axom::detail::getAllocatorID<SPACE>())
#endif
  { }

  /*!
   * \brief Generic constructor for an ArrayView of arbitrary dimension with external data
   *
   *  This constructor assumes that the data is laid out in a dense row-major
   *  layout; that is, elements are laid out contiguously with respect to the
   *  last (DIM-1) dimension.
   *
   * \param [in] data the external data this ArrayView will wrap.
   * \param [in] args The parameter pack containing the "shape" of the ArrayView
   *
   * \pre sizeof...(Args) == DIM
   */
  template <typename... Args,
            typename Enable = typename std::enable_if<
              sizeof...(Args) == DIM && detail::all_types_are_integral<Args...>::value>::type>
  AXOM_HOST_DEVICE ArrayView(T* data, Args... args);

  /*!
   * \brief Generic constructor for an ArrayView of arbitrary dimension with external data
   *
   *  This constructor assumes that the data is laid out in a row-major layout,
   *  but accepts a min_stride parameter for the stride between consecutive
   *  elements in the last dimension.
   *
   * \param [in] data the external data this ArrayView will wrap.
   * \param [in] shape Array size in each dimension.
   * \param [in] min_stride Minimum stride between consecutive items in the
   *  last dimension
   */
  AXOM_HOST_DEVICE ArrayView(T* data,
                             const StackArray<IndexType, DIM>& shape,
                             IndexType min_stride = 1);

  /*!
   * \brief Generic constructor for an ArrayView of arbitrary dimension with external data
   *
   * \param [in] data the external data this ArrayView will wrap.
   * \param [in] shape Array size in each dimension.
   * \param [in] strides Array strides for each dimension.
   */
  AXOM_HOST_DEVICE ArrayView(T* data,
                             const StackArray<IndexType, DIM>& shape,
                             const StackArray<IndexType, DIM>& stride);

  /*! 
   * \brief Constructor for transferring between memory spaces
   * 
   * \param [in] other The array in a different memory space to copy from
   * 
   * \note The parameter is non-const because \a other can be modified through the constructed View
   * 
   * \note This constructor is left implicit to allow for convenient function calls that convert
   * from \p Array -> \p ArrayView or from dynamic memory spaces to an \p ArrayView of explicitly specified
   * space.
   */
  template <typename OtherArrayType>
  AXOM_HOST_DEVICE ArrayView(ArrayBase<T, DIM, OtherArrayType>& other);
  /// \overload
  template <typename OtherArrayType>
  AXOM_HOST_DEVICE ArrayView(
    const ArrayBase<typename std::remove_const<T>::type, DIM, OtherArrayType>& other);

  /*!
   * \brief Return the number of elements stored in the data array.
   */
  inline AXOM_HOST_DEVICE IndexType size() const { return m_num_elements; }

  /*!
   * \brief Returns true iff the ArrayView stores no elements.
   */
  AXOM_HOST_DEVICE
  inline bool empty() const { return m_num_elements == 0; }

  /*!
   * \brief Returns an ArrayViewIterator to the first element of the Array
   */
  AXOM_HOST_DEVICE
  ArrayViewIterator begin() const { return ArrayViewIterator(0, this); }

  /*!
   * \brief Returns an ArrayViewIterator to the element following the last
   *  element of the Array.
   */
  AXOM_HOST_DEVICE
  ArrayViewIterator end() const { return ArrayViewIterator(size(), this); }

  /*!
   * \brief Return a pointer to the array of data.
   */
  /// @{

  AXOM_HOST_DEVICE inline T* data() const { return m_data; }

  /// @}

  /*!
   * \brief Get the ID for the umpire allocator
   */
  AXOM_HOST_DEVICE int getAllocatorID() const { return m_allocator_id; }

  /*!
   * \brief Returns an ArrayView that is a subspan of the original range of
   *  elements.
   *
   * \param [in] offset The index where the subspan should begin.
   * \param [in] count The number of elements to include in the subspan, or -1
   *  to take all elements after offset (default).
   *
   * \return An ArrayView that spans the indices [offset, offset + count),
   *  or [offset, num_elements) if count < 0.
   *
   * \pre offset + count <= m_num_elements if count < 0
   */
  template <int UDIM = DIM, typename Enable = typename std::enable_if<UDIM == 1>::type>
  AXOM_HOST_DEVICE ArrayView subspan(IndexType offset, IndexType count = -1) const
  {
    assert(offset >= 0 && offset < m_num_elements);
    if(count >= 0)
    {
      assert(offset + count <= m_num_elements);
    }

    ArrayView slice = *this;
    slice.m_data += offset;
    if(count < 0)
    {
      slice.m_num_elements -= offset;
    }
    else
    {
      slice.m_num_elements = count;
    }
    return slice;
  }

  /*!
   * \brief Returns an ArrayView that is a subspan of the original range of
   *  elements.
   *
   * \param [in] offset The index in each dimension where the subspan should
   *  begin.
   * \param [in] count The number of elements in each dimension to include in
   *  the subspan. -1 in any dimension indicates to take all elements after
   *  the offset for that domain.
   *
   * \return An ArrayView that spans in each dimension the indices:
   *   * [offsets[i], offsets[i] + counts[i])
   *   * or [offsets[i], shape[i]) if count < 0.
   *
   * \pre offsets[i] + counts[i] < shape[i] for 0 <= i < DIM
   */
  AXOM_HOST_DEVICE ArrayView subspan(const StackArray<IndexType, DIM>& offsets,
                                     const StackArray<IndexType, DIM>& counts);

  /*!
   * \brief Fill the ArrayView with a value.
   *
   * \param value The value to be used for filling the ArrayView.
   */
  void fill(const T& value);

  /*!
   * \brief Set a range of elements to a given value.
   *
   * \param [in] value the value to set to.
   * \param [in] n the number of elements to write.
   * \param [in] pos the position at which to begin the fill operation.
   *
   * \note The size is unchanged by calls to fill.
   *
   * \pre pos + n <= m_num_elements.
   */
  void fill(const T& value, IndexType n, IndexType pos);

  /*!
   * \brief Modify the values of existing elements.
   *
   * \param [in] elements the new elements to use in the set operation.
   * \param [in] n the number of elements in the set operation.
   * \param [in] pos the position at which to begin replacing values.
   *
   * \note It's assumed that elements is of length n.
   * \note The size is unchanged by calls to set.
   *
   * \pre pos + n <= m_num_elements.
   */
  void set(const T* elements, IndexType n, IndexType pos);

  /*!
   * \brief Set the array view contents.
   *
   * \param [in] count The new number of elements.
   * \param [in] value The value to store in the elements.
   *
   * \note It's assumed that it is safe to store \a count elements.
   * \note The size is set to \a count by calls to assign.
   */
  void assign(axom::IndexType count, const T& value);

  /*!
   * \brief Replaces contents with copies of objects in the range [first, last).
   *
   * \param [in] first The iterator that begins the range used for assignment.
   * \param [in] last The iterator that ends the range used for assignment.
   *                  The value at this iterator is not assigned into to the view.
   *
   * \post Size of ArrayView is changed to the number of items in the iterator
   *       range, and values values referenced by the iterator range are copied
   *       into the Array.
   */
  template <class InputIt>
  void assign(InputIt first, InputIt last);

  /*!
   * \brief Set the array contents using an initializer list.
   *
   * \param [in] elems An initializer list containing the new array values.
   *
   * \post The ArrayView contains copies of the initializer list elements.
   */
  void assign(std::initializer_list<T> elems);

private:
  T* m_data = nullptr;
  /// \brief The full number of elements in the array
  ///  i.e., 3 for a 1D Array of size 3, 9 for a 3x3 2D array, etc
  IndexType m_num_elements = 0;
  /// \brief The allocator ID for the memory space in which m_data was allocated
  int m_allocator_id;
};

/// \brief Helper alias for multi-component arrays
template <typename T>
using MCArrayView = ArrayView<T, 2>;

//------------------------------------------------------------------------------
//                            ArrayView IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
template <typename... Args, typename Enable>
AXOM_HOST_DEVICE ArrayView<T, DIM, SPACE>::ArrayView(T* data, Args... args)
  : ArrayView(data, StackArray<IndexType, DIM> {{static_cast<IndexType>(args)...}})
{
  static_assert(sizeof...(Args) == DIM, "Array size must match number of dimensions");
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
AXOM_HOST_DEVICE ArrayView<T, DIM, SPACE>::ArrayView(T* data,
                                                     const StackArray<IndexType, DIM>& shape,
                                                     IndexType min_stride)
  : ArrayBase<T, DIM, ArrayView<T, DIM, SPACE>>(shape, min_stride)
  , m_data(data)
#ifndef AXOM_DEVICE_CODE
  , m_allocator_id(axom::detail::getAllocatorID<SPACE>())
#endif
{
#if defined(AXOM_DEVICE_CODE) && defined(AXOM_USE_UMPIRE)
  static_assert((SPACE != MemorySpace::Constant) || std::is_const<T>::value,
                "T must be const if memory space is Constant memory");
#endif
  // Intel hits internal compiler error when casting as part of function call
  if(data != nullptr)
  {
    m_num_elements = detail::packProduct(shape.m_data);
  }
  else
  {
    m_num_elements = 0;
  }

#if !defined(AXOM_DEVICE_CODE) && defined(AXOM_USE_UMPIRE)
  // If we have Umpire, we can try and see what space the pointer is allocated in
  // Probably not worth checking this if SPACE != Dynamic, we *could* error out
  // if e.g., the user gives a host pointer to ArrayView<T, DIM, Device>, but even
  // Thrust doesn't guard against this.

  // FIXME: Is it worth trying to get rid of this at compile time?
  // (using a workaround since we don't have "if constexpr")
  if(SPACE == MemorySpace::Dynamic)
  {
    auto& rm = umpire::ResourceManager::getInstance();

    using NonConstT = typename std::remove_const<T>::type;
    // TODO: There's no reason these Umpire methods should take a non-const pointer.
    if(rm.hasAllocator(const_cast<NonConstT*>(data)))
    {
      auto alloc = rm.getAllocator(const_cast<NonConstT*>(data));
      m_allocator_id = alloc.getId();
    }
  }
#endif
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
AXOM_HOST_DEVICE ArrayView<T, DIM, SPACE>::ArrayView(T* data,
                                                     const StackArray<IndexType, DIM>& shape,
                                                     const StackArray<IndexType, DIM>& stride)
  : ArrayBase<T, DIM, ArrayView<T, DIM, SPACE>>(shape, stride)
  , m_data(data)
#ifndef AXOM_DEVICE_CODE
  , m_allocator_id(axom::detail::getAllocatorID<SPACE>())
#endif
{
#if defined(AXOM_DEVICE_CODE) && defined(AXOM_USE_UMPIRE)
  static_assert((SPACE != MemorySpace::Constant) || std::is_const<T>::value,
                "T must be const if memory space is Constant memory");
#endif
  // Intel hits internal compiler error when casting as part of function call
  if(data != nullptr)
  {
    m_num_elements = detail::packProduct(shape.m_data);
  }
  else
  {
    m_num_elements = 0;
  }

#if !defined(AXOM_DEVICE_CODE) && defined(AXOM_USE_UMPIRE)
  // If we have Umpire, we can try and see what space the pointer is allocated in
  // Probably not worth checking this if SPACE != Dynamic, we *could* error out
  // if e.g., the user gives a host pointer to ArrayView<T, DIM, Device>, but even
  // Thrust doesn't guard against this.

  // FIXME: Is it worth trying to get rid of this at compile time?
  // (using a workaround since we don't have "if constexpr")
  if(SPACE == MemorySpace::Dynamic)
  {
    auto& rm = umpire::ResourceManager::getInstance();

    using NonConstT = typename std::remove_const<T>::type;
    // TODO: There's no reason these Umpire methods should take a non-const pointer.
    if(rm.hasAllocator(const_cast<NonConstT*>(data)))
    {
      auto alloc = rm.getAllocator(const_cast<NonConstT*>(data));
      m_allocator_id = alloc.getId();
    }
  }
#endif
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
template <typename OtherArrayType>
AXOM_HOST_DEVICE ArrayView<T, DIM, SPACE>::ArrayView(ArrayBase<T, DIM, OtherArrayType>& other)
  : ArrayBase<T, DIM, ArrayView<T, DIM, SPACE>>(other)
  , m_data(static_cast<OtherArrayType&>(other).data())
  , m_num_elements(static_cast<OtherArrayType&>(other).size())
  , m_allocator_id(static_cast<OtherArrayType&>(other).getAllocatorID())
{
#if !defined(AXOM_DEVICE_CODE) && defined(AXOM_DEBUG)
  // If it's not dynamic, the allocator ID from the argument array has to match the template param.
  // If that's not the case then things have gone horribly wrong somewhere.
  if(SPACE != MemorySpace::Dynamic && SPACE != axom::detail::getAllocatorSpace(m_allocator_id))
  {
    std::cerr << "Input argument allocator does not match the explicitly "
                 "provided memory space\n";
    utilities::processAbort();
  }
#endif
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
template <typename OtherArrayType>
AXOM_HOST_DEVICE ArrayView<T, DIM, SPACE>::ArrayView(
  const ArrayBase<typename std::remove_const<T>::type, DIM, OtherArrayType>& other)
  : ArrayBase<T, DIM, ArrayView<T, DIM, SPACE>>(other)
  , m_data(static_cast<const OtherArrayType&>(other).data())
  , m_num_elements(static_cast<const OtherArrayType&>(other).size())
  , m_allocator_id(static_cast<const OtherArrayType&>(other).getAllocatorID())
{
  static_assert(std::is_const<T>::value || detail::ArrayTraits<OtherArrayType>::is_view,
                "Cannot create an ArrayView of non-const type from a const Array");
#if !defined(AXOM_DEVICE_CODE) && defined(AXOM_DEBUG)
  // If it's not dynamic, the allocator ID from the argument array has to match the template param.
  // If that's not the case then things have gone horribly wrong somewhere.
  if(SPACE != MemorySpace::Dynamic && SPACE != axom::detail::getAllocatorSpace(m_allocator_id))
  {
    std::cerr << "Input argument allocator does not match the explicitly "
                 "provided memory space\n";
    utilities::processAbort();
  }
#endif
}

//------------------------------------------------------------------------------
AXOM_SUPPRESS_HD_WARN
template <typename T, int DIM, MemorySpace SPACE>
AXOM_HOST_DEVICE ArrayView<T, DIM, SPACE> ArrayView<T, DIM, SPACE>::subspan(
  const StackArray<IndexType, DIM>& offsets,
  const StackArray<IndexType, DIM>& counts)
{
#ifdef AXOM_DEBUG
  for(int dim = 0; dim < DIM; dim++)
  {
    assert(offsets[dim] >= 0 && offsets[dim] < this->shape()[dim]);
    if(counts[dim] >= 0)
    {
      assert(offsets[dim] + counts[dim] <= this->shape()[dim]);
    }
  }
#endif
  // Compute flat offset into existing data.
  IndexType offset = numerics::dot_product(offsets.m_data, this->strides().m_data, DIM);

  // Setup new shape array.
  StackArray<IndexType, DIM> newShape;
  for(int dim = 0; dim < DIM; dim++)
  {
    if(counts[dim] < 0)
    {
      newShape[dim] = this->shape()[dim] - offsets[dim];
    }
    else
    {
      newShape[dim] = counts[dim];
    }
  }

  ArrayView slice = *this;
  slice.m_data += offset;
  slice.m_num_elements = detail::packProduct(newShape.m_data);
  slice.setShapeAndStride(newShape, this->strides());
  return slice;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
void ArrayView<T, DIM, SPACE>::fill(const T& value)
{
  using OpHelper = detail::ArrayOps<T>;
  const bool executeOnGPU = axom::isDeviceAllocator(m_allocator_id);
  OpHelper {m_allocator_id, executeOnGPU}.destroy(m_data, 0, m_num_elements);
  OpHelper {m_allocator_id, executeOnGPU}.fill(m_data, 0, m_num_elements, value);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
void ArrayView<T, DIM, SPACE>::fill(const T& value, IndexType n, IndexType pos)
{
  assert(pos >= 0);
  assert(pos + n <= m_num_elements);

  using OpHelper = detail::ArrayOps<T>;
  const bool executeOnGPU = axom::isDeviceAllocator(m_allocator_id);
  OpHelper {m_allocator_id, executeOnGPU}.destroy(m_data, pos, n);
  OpHelper {m_allocator_id, executeOnGPU}.fill(m_data, pos, n, value);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
void ArrayView<T, DIM, SPACE>::set(const T* elements, IndexType n, IndexType pos)
{
  assert(pos >= 0);
  assert(pos + n <= m_num_elements);

  using OpHelper = detail::ArrayOps<T>;
  const bool executeOnGPU = axom::isDeviceAllocator(m_allocator_id);
  OpHelper {m_allocator_id, executeOnGPU}.destroy(m_data, pos, n);
  OpHelper {m_allocator_id, executeOnGPU}.fill_range(m_data, pos, n, elements, MemorySpace::Dynamic);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
void ArrayView<T, DIM, SPACE>::assign(axom::IndexType count, const T& value)
{
  assert(count >= 0);

  using OpHelper = detail::ArrayOps<T>;
  const bool executeOnGPU = axom::isDeviceAllocator(m_allocator_id);
  OpHelper {m_allocator_id, executeOnGPU}.destroy(m_data, 0, m_num_elements);
  OpHelper {m_allocator_id, executeOnGPU}.fill(m_data, 0, count, value);
  m_num_elements = count;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
template <class InputIt>
void ArrayView<T, DIM, SPACE>::assign(InputIt first, InputIt last)
{
  using OpHelper = detail::ArrayOps<T>;
  const bool executeOnGPU = axom::isDeviceAllocator(m_allocator_id);
  OpHelper {m_allocator_id, executeOnGPU}.destroy(m_data, 0, m_num_elements);
  m_num_elements = 0;
  for(auto it = first; it != last; it++)
  {
    OpHelper {m_allocator_id, executeOnGPU}.emplace(m_data, m_num_elements++, *it);
  }
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
void ArrayView<T, DIM, SPACE>::assign(std::initializer_list<T> elems)
{
  assign(elems.begin(), elems.end());
}

} /* namespace axom */

#endif /* AXOM_ARRAYVIEW_HPP_ */
