// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_ARRAY_HPP_
#define AXOM_ARRAY_HPP_

#include "axom/config.hpp"
#include "axom/core/MDMapping.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/execution/atomics.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/core/Types.hpp"
#include "axom/core/ArrayBase.hpp"
#include "axom/core/ArrayIteratorBase.hpp"
#include "axom/core/ArrayView.hpp"

// C/C++ includes
#include <algorithm>
#include <iostream>

namespace axom
{
namespace ArrayOptions
{
/// \brief A "tag type" for constructing an Array without initializing its memory
struct Uninitialized
{ };

}  // namespace ArrayOptions

// Forward declare the templated classes and operator function(s)
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
class Array;

namespace detail
{
// Static information to pass to ArrayBase
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
struct ArrayTraits<Array<T, DIM, SPACE, StoragePolicy>>
{
  constexpr static bool is_view = false;
  constexpr static int dimension = DIM;
};

/*!
 * \class DefaultStoragePolicy
 *
 * \brief Default storage policy for axom::Array.
 *  Uses Umpire to reallocate buffers.
 */
template <typename T>
struct DefaultStoragePolicy
{
  /*!
   * \brief Callback to report changes in shape of data in Array.
   *
   * \param [in] shape the current dimensions of the array
   */
  template <int Dims>
  void onShapeUpdate(StackArray<IndexType, Dims> AXOM_UNUSED_PARAM(shape))
  { }

  /*!
   * \brief Reallocates a buffer.
   *
   * \param [in] old_data pointer to the old buffer
   * \param [in] old_capacity the capacity of the currently allocated buffer
   * \param [in] allocator_id the allocator ID to use
   * \param [in] new_capacity the capacity to allocate
   * \param [in] nontrivial_move a callback to move elements that aren't
   *  trivially copyable
   *
   * \return a pointer to the new buffer with moved elements
   */
  template <typename Func>
  T* reallocate(T* old_data,
                int AXOM_UNUSED_PARAM(old_capacity),
                int allocator_id,
                int new_capacity,
                Func&& nontrivial_move)
  {
    // Create a new block of memory, and move the elements over.
    T* new_data = axom::allocate<T>(new_capacity, allocator_id);
    nontrivial_move(new_data);

    // Destroy the original array.
    axom::deallocate(old_data);

    return new_data;
  }

  /*!
   * \brief Frees a buffer.
   */
  void deallocate(T* data)
  {
    if(data != nullptr)
    {
      axom::deallocate(data);
    }
  }
};

}  // namespace detail

/*!
 * \class Array
 *
 * \brief Provides a generic multidimensional array container.
 *
 *  The Array class provides a generic multidimensional array 
 *  container with dynamic reallocation and insertion.  The dimensionality
 *  of the array must be known at compile time but the extents in each dimension
 *  are dynamic and can be changed at runtime.  Array elements are stored
 *  contiguously.
 *
 *  \note For a multi-component array container, where each element
 *  is a tuple of 1 or more components, Axom provides the MCArray alias, which
 *  corresponds to Array<T, 2>.
 *
 *  The Array class mirrors std::vector, with future support for GPUs
 *  in-development.  The class's multidimensional array functionality roughly
 *  mirrors the multidimensional array support provided by numpy's ndarray.
 *
 *  \see https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html
 *
 *  Some interfaces accomodate data ordering specifications.
 *  Unless otherwise specified, data storage defaults to row-major.
 *
 *  This class is meant to be a drop-in replacement for std::vector.
 *  However, it differs in its memory management and construction semantics.
 *  Specifically, we do not require axom::Array to initialize/construct
 *  its memory at allocation time and we use axom's memory_management
 *  and allocator ID abstractions rather than std::allocator.
 *
 *  Array always retains exclusive ownership of its data and is responsible for
 *  freeing its memory.
 *
 *  \see ArrayView for non-owning views of one- or multi-dimensional data
 *  Depending on which constructor is used, the Array object can have two
 *  different underlying storage types:
 *
 * \tparam T the type of the values to hold.
 * \tparam DIM The dimension of the array.
 * \tparam SPACE The memory space of the array.
 * 
 * \pre T must be CopyAssignable and Erasable
 * \see https://en.cppreference.com/w/cpp/named_req
 *
 * \pre When Array is allocated on the device, T must be relocatable, i.e.
 *  moving the object to a new index and destroying the original should be
 *  equivalent to a memcpy.
 * \see https://github.com/facebook/folly/blob/main/folly/docs/FBVector.md#object-relocation
 *
 */
template <typename T,
          int DIM = 1,
          MemorySpace SPACE = MemorySpace::Dynamic,
          typename StoragePolicy = detail::DefaultStoragePolicy<T>>
class Array : public ArrayBase<T, DIM, Array<T, DIM, SPACE, StoragePolicy>>, protected StoragePolicy
{
public:
  static constexpr double DEFAULT_RESIZE_RATIO = 2.0;
  static constexpr IndexType MIN_DEFAULT_CAPACITY = 32;
  using value_type = T;
  static constexpr MemorySpace space = SPACE;
  using ArrayIterator = ArrayIteratorBase<Array<T, DIM, SPACE, StoragePolicy>, T>;
  using ConstArrayIterator = ArrayIteratorBase<const Array<T, DIM, SPACE, StoragePolicy>, const T>;

  using ArrayViewType = ArrayView<T, DIM, SPACE>;
  using ConstArrayViewType = ArrayView<const T, DIM, SPACE>;

private:
  using OpHelper = detail::ArrayOps<T>;

public:
  /// \name Native Storage Array Constructors
  /// @{

  /*! 
   * \brief Default constructor. Constructs an Array instance with no elements
   *  and default allocator ID. 
   *
   */
  Array();

  /*!
   * \brief Constructs a 1D Array instance with the given number of elements.
   *
   * \param [in] num_elements the number of elements the Array holds.
   * \param [in] capacity the number of elements to allocate space for.
   * \param [in] allocator_id the ID of the allocator to use (optional)
   *
   * \note If no capacity or capacity less than num_elements is specified
   *  then it will default to at least num_elements * DEFAULT_RESIZE_RATIO.
   * \note a capacity is specified for the number of elements to store in the
   *  array and does not correspond to the actual bytesize.
   * \note The option to select a capacity is only available for 1-dimensional Arrays
   * 
   * \note Some overloads have an `ArrayOptions::Uninitialized` first parameter.
   * These are intended for cases where the array data should not be initialized
   * when memory is allocated, e.g. if the code is known to initialize the data
   *
   * \pre num_elements >= 0
   *
   * \post capacity() >= size()
   * \post size() == num_elements
   * \post getResizeRatio() == DEFAULT_RESIZE_RATIO
   */
  template <IndexType SFINAE_DIM = DIM,
            MemorySpace SFINAE_SPACE = SPACE,
            typename std::enable_if<SFINAE_DIM == 1>::type* = nullptr>
  Array(IndexType num_elements,
        IndexType capacity = 0,
        int allocator_id = axom::detail::getAllocatorID<SPACE>());

  /// \overload
  template <IndexType SFINAE_DIM = DIM,
            MemorySpace SFINAE_SPACE = SPACE,
            typename std::enable_if<SFINAE_DIM == 1>::type* = nullptr>
  Array(ArrayOptions::Uninitialized,
        IndexType num_elements,
        IndexType capacity = 0,
        int allocator_id = axom::detail::getAllocatorID<SPACE>());

  Array(const axom::StackArray<axom::IndexType, DIM>& shape,
        int allocator_id = axom::detail::getAllocatorID<SPACE>());

  /*!
    \brief Construct Array with row- or column-major data ordering.

    \pre rowOrColumn must be either ArrayStrideOrder::ROW or
    ArrayStrideOrder::COLUMN (or, if DIM is 1, ArrayStrideOrder::BOTH).
  */
  Array(const axom::StackArray<axom::IndexType, DIM>& shape,
        axom::ArrayStrideOrder rowOrColumn,
        int allocator_id = axom::detail::getAllocatorID<SPACE>());

  /*!
    \brief Construct Array with data ordering specifications.

    Example of 3D Array where j is slowest and i is fastest:
        Array<3, int> ar(
          shape,
          axom::StackArray<axom::IndexType, DIM>{2, 0, 1});

    \pre slowestDirs must be a permutation of [0 ... DIM-1]
  */
  template <typename DirType = axom::IndexType>
  Array(const axom::StackArray<axom::IndexType, DIM>& shape,
        const axom::StackArray<DirType, DIM>& slowestDirs,
        int allocator_id = axom::detail::getAllocatorID<SPACE>());

  /*!
   * \brief Generic constructor for an Array of arbitrary dimension
   *
   * \param [in] args The parameter pack containing the "shape" of the Array
   * \see https://numpy.org/doc/stable/reference/generated/numpy.empty.html#numpy.empty
   *
   * \pre sizeof...(Args) == DIM
   *
   * \post capacity() >= size()
   * \post size() == num_elements
   * \post getResizeRatio() == DEFAULT_RESIZE_RATIO
   */
  template <typename... Args,
            typename Enable = typename std::enable_if<
              sizeof...(Args) == DIM && detail::all_types_are_integral<Args...>::value>::type>
  Array(Args... args);

  /// \overload
  template <typename... Args,
            typename Enable = typename std::enable_if<
              sizeof...(Args) == DIM && detail::all_types_are_integral<Args...>::value>::type>
  Array(ArrayOptions::Uninitialized, Args... args);

  /*!
   * \brief Initializer list constructor for a one-dimensional Array
   *
   * \param [in] elems The elements to initialize the array with
   * \param [in] allocator_id the ID of the allocator to use (optional)
   */
  template <int UDIM = DIM, typename Enable = typename std::enable_if<UDIM == 1>::type>
  Array(std::initializer_list<T> elems, int allocator_id = axom::detail::getAllocatorID<SPACE>());

  /*! 
   * \brief Copy constructor for an Array instance 
   */
  AXOM_HOST_DEVICE Array(const Array& other);

  /*! 
   * \brief Move constructor for an Array instance 
   */
  Array(Array&& other) noexcept;

  /*!
   * \brief Constructor for transferring between memory spaces
   *
   * \param [in] other The array in a different memory space to copy from
   *
   * \note The new Array will be constructed with the allocator ID of the
   *  copied-from array, if compatible with the memory space of the target
   *  array. Otherwise, the new Array will be constructed with the default
   *  allocator for the memory space.
   *
   *  An Array specified with the default Dynamic memory space will always
   *  propagate the allocator ID from the source array.
   */
  template <typename OtherArrayType>
  Array(const ArrayBase<T, DIM, OtherArrayType>& other);

  /// \overload
  template <typename OtherArrayType>
  Array(const ArrayBase<const T, DIM, OtherArrayType>& other);

  /*!
   * \brief Constructor for transferring between memory spaces, with a user-
   *  specified allocator
   *
   * \param [in] other The array in a different memory space to copy from
   * \param [in] allocator_id the ID of the allocator to use
   *
   * \note The specified allocator ID must be compatible with the memory space
   *  of the new Array.
   */
  template <typename OtherArrayType>
  Array(const ArrayBase<T, DIM, OtherArrayType>& other, int allocator_id);

  /// \overload
  template <typename OtherArrayType>
  Array(const ArrayBase<const T, DIM, OtherArrayType>& other, int allocator_id);

  /// @}

  /// \name Array copy and move operators
  /// @{

  /*! 
   * \brief Copy assignment operator for Array
   * 
   * \pre T must be TriviallyCopyable
   */
  Array& operator=(const Array& other)
  {
    if(this != &other)
    {
      this->clear();
      static_cast<ArrayBase<T, DIM, Array<T, DIM, SPACE, StoragePolicy>>&>(*this) = other;
      m_allocator_id = other.m_allocator_id;
      m_executeOnGPU = axom::isDeviceAllocator(m_allocator_id);
      m_arrayOps = other.m_arrayOps;
      m_resize_ratio = other.m_resize_ratio;
      setCapacity(other.capacity());
      // Use fill_range to ensure that copy constructors are invoked for each element
      MemorySpace srcSpace = SPACE;
      if(srcSpace == MemorySpace::Dynamic)
      {
        srcSpace = axom::detail::getAllocatorSpace(other.m_allocator_id);
      }
      m_arrayOps.fill_range(m_data, 0, other.size(), other.data(), srcSpace);
      updateNumElements(other.size());
    }

    return *this;
  }

  /*! 
   * \brief Move assignment operator for Array
   */
  Array& operator=(Array&& other) noexcept
  {
    if(this != &other)
    {
      this->clear();
      StoragePolicy::deallocate(m_data);
      static_cast<ArrayBase<T, DIM, Array>&>(*this) = std::move(other);

      m_data = other.m_data;
      m_num_elements = other.m_num_elements;
      m_capacity = other.m_capacity;
      m_resize_ratio = other.m_resize_ratio;
      m_allocator_id = other.m_allocator_id;
      m_executeOnGPU = axom::isDeviceAllocator(m_allocator_id);
      m_arrayOps = other.m_arrayOps;

      other.m_data = nullptr;
      other.m_num_elements = 0;
      other.m_capacity = 0;
      other.m_resize_ratio = DEFAULT_RESIZE_RATIO;
      other.m_allocator_id = INVALID_ALLOCATOR_ID;
    }

    return *this;
  }

  /*!
   * \brief Initializer list assignment operator for Array.
   *
   * \param [in] elems the elements to set the array to.
   */
  template <int UDIM = DIM, typename Enable = typename std::enable_if<UDIM == 1>::type>
  Array& operator=(std::initializer_list<T> elems)
  {
    clear();
    insert(0, elems.size(), elems.begin());
    return *this;
  }

  /// @}

  /*!
   * Destructor. Frees the associated buffer.
   */
  ~Array();

  /// \name Array element access operators
  /// @{

  // TODO: Implement View class for the case where sizeof...(Args) < DIM (i.e., where the indexing results in a nonscalar)

  /*!
   * \brief Return a pointer to the array of data.
   */
  /// @{

  AXOM_HOST_DEVICE inline T* data() { return m_data; }
  AXOM_HOST_DEVICE inline const T* data() const { return m_data; }

  /// @}

  /*!
    @brief Convert 1D Array into a StackArray.
  */
  template <int LENGTH1D, typename TT = T, int TDIM = DIM>
  AXOM_HOST_DEVICE inline typename std::enable_if<TDIM == 1, axom::StackArray<TT, LENGTH1D>>::type
  to_stack_array() const
  {
    axom::StackArray<TT, LENGTH1D> rval;
    IndexType copyCount = LENGTH1D <= m_num_elements ? LENGTH1D : m_num_elements;
    for(IndexType i = 0; i < copyCount; ++i)
    {
      rval[i] = m_data[i];
    }
    return rval;
  }

  /// @}

  /// \name Array methods to modify the data.
  /// @{

  /*!
   * \brief Set all the values of the array.
   *
   * \param [in] value the value to set to.
   */
  void fill(const T& value);

  /*!
   * \brief Set a range of elements to a given value.
   *
   * \param [in] value the value to set to.
   * \param [in] n the number of elements to write.
   * \param [in] pos the position at which to begin writing.
   *
   * \note The size is unchanged by calls to fill.
   *
   * \pre pos + n <= m_num_elements.
   */
  void fill(const T& value, IndexType n, IndexType pos);

  /*!
   * \brief Modify the values of existing elements.
   *
   * \param [in] elements the new elements to write.
   * \param [in] n the number of elements to write.
   * \param [in] pos the position at which to begin writing.
   *
   * \note It's assumed that elements is of length n.
   * \note The size is unchanged by calls to set.
   *
   * \pre pos + n <= m_num_elements.
   */
  void set(const T* elements, IndexType n, IndexType pos);

  /*!
   * \brief Set the array contents.
   *
   * \param [in] count The new number of elements.
   * \param [in] value The value to store in the elements.
   *
   * \post Size of Array is \a count, all elements contain \a value.
   */
  void assign(axom::IndexType count, const T& value);

  /*!
   * \brief Replaces contents with copies of objects in [first, last).
   *
   * \param [in] first The iterator for the first value to use in a container.
   * \param [in] last The last iterator to use in a container.
   *
   * \post Size of Array is changed to the number of items in the range
   *       designated by the iterators (last-first), and values referenced by
   *       the iterator range are copied into the Array.
   */
  template <class InputIt>
  void assign(InputIt first, InputIt last);

  /*!
   * \brief Set the array contents using an initializer list.
   *
   * \param [in] elems An initializer list containing the new array values.
   *
   * \post The Array contains copies of the initializer list elements.
   */
  void assign(std::initializer_list<T> elems);

  /*!
   * \brief Clears the contents of the array
   * 
   * \post size of Array is 0
   * \post capacity is unchanged
   */
  void clear();

  /*!
   * \brief Insert an element into the array at the given position.
   *
   * \param [in] pos the position at which to insert.
   * \param [in] value the element value to insert.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note The size increases by 1.
   *
   */
  void insert(IndexType pos, const T& value);

  /*!
   * \brief Insert an element into the array at the value before pos.
   *
   * \param [in] pos the ArrayIterator before which value will be inserted.
   * \param [in] value the element value to insert.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note The size increases by 1.
   *
   * \return ArrayIterator to inserted value
   */
  ArrayIterator insert(ArrayIterator pos, const T& value);

  /*!
   * \brief Insert elements into the array at the given position.
   *
   * \param [in] pos the position at which to begin the insertion.
   * \param [in] n the number of elements to insert.
   * \param [in] values the element values to insert.
   *
   * \note It's assumed that elements is of length n.
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note The size increases by n.
   *
   * \pre pos <= m_num_elements.
   */
  void insert(IndexType pos, IndexType n, const T* values);

  /*!
   * \brief Insert elements into the array at the value before pos.
   *
   * \param [in] pos the ArrayIterator before which value will be inserted.
   * \param [in] n the number of elements to insert.
   * \param [in] values the element values to insert.
   *
   * \note It's assumed that elements is of length n.
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note The size increases by n.
   *
   * \pre pos <= end()
   *
   * \return ArrayIterator to first element inserted (pos if n == 0)
   */
  ArrayIterator insert(ArrayIterator pos, IndexType n, const T* values);

  /*!
   * \brief Insert n copies of element into the array at the given position.
   *
   * \param [in] pos the position at which to begin the insertion.
   * \param [in] n the number of elements to insert.
   * \param [in] value the element value to insert.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note The size increases by n.
   * \note This method is used to create space for elements in the middle of
   *  the array.
   *
   * \pre pos <= m_num_elements.
   */
  void insert(IndexType pos, IndexType n, const T& value);

  /*!
   * \brief Insert n copies of element into the array at the value before pos.
   *
   * \param [in] pos the ArrayIterator before which value will be inserted.
   * \param [in] n the number of elements to insert.
   * \param [in] value the element value to insert.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note The size increases by n.
   * \note This method is used to create space for elements in the middle of
   *  the array.
   *
   * \pre pos <= end()
   *
   * \return ArrayIterator to first element inserted (pos if n == 0)
   */
  ArrayIterator insert(ArrayIterator pos, IndexType n, const T& value);

  /*!
   * \brief Inserts an Array to the end of the calling object
   *
   * \param [in] other The Array to append
   *
   * \pre The shapes of the calling Array and @a other are the same
   * (excluding the leading dimension), i.e., shape()[1:] == other.shape()[1:]
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  template <MemorySpace OtherSpace>
  void insert(IndexType pos, ArrayView<const T, DIM, OtherSpace> other);

  /// \overload
  template <MemorySpace OtherSpace>
  void insert(IndexType pos, ArrayView<T, DIM, OtherSpace> other)
  {
    insert(pos, ArrayView<const T, DIM, OtherSpace>(other));
  }

  /*!
   * \brief Appends an Array to the end of the calling object
   *
   * \param [in] other The Array to append
   * \tparam OtherArrayType The underlying type of the other array
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  template <MemorySpace OtherSpace>
  void append(ArrayView<const T, DIM, OtherSpace> other)
  {
    insert(size(), other);
  }

  /// \overload
  template <MemorySpace OtherSpace>
  void append(ArrayView<T, DIM, OtherSpace> other)
  {
    insert(size(), other);
  }

  /*!
   * \brief Erases an element from the Array 
   *
   * \param [in] pos the ArrayIterator to the element in the Array
   *
   * \return An ArrayIterator following the last element removed.
   */
  ArrayIterator erase(ArrayIterator pos);

  /*!
   * \brief Erases elements in the range [first, last) from the Array
   *
   * \param [in] first the ArrayIterator to the beginning of the range.
   * \param [in] last the ArrayIterator to end of range.
   *
   * \return An ArrayIterator following the last element removed. 
   */
  ArrayIterator erase(ArrayIterator first, ArrayIterator last);

  /*!
   * \brief Inserts new element into Array at the given position.
   *
   * \param [in] pos the position to insert element at.
   * \param [in] args the arguments to forward to constructor of the element.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note The size increases by 1.
   *
   * \pre T must be MoveAssignable
   */
  template <typename... Args>
  void emplace(IndexType pos, Args&&... args);

  /*!
   * \brief Inserts new element into Array before pos.
   *
   * \param [in] pos the ArrayIterator to insert element before.
   * \param [in] args the arguments to forward to constructor of the element.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note The size increases by 1.
   *
   * \pre T must be MoveAssignable
   *
   * \return An ArrayIterator to the emplaced element.
   */
  template <typename... Args>
  ArrayIterator emplace(ArrayIterator pos, Args&&... args);

  /*!
   * \brief Push a value to the back of the array.
   *
   * \param [in] value the value to be added to the back.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   *
   * \pre DIM == 1
   */
  void push_back(const T& value);

  /*!
   * \brief Push a value to the back of the array.
   *
   * \param [in] value the value to move to the back.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   *
   * \pre DIM == 1
   */
  void push_back(T&& value);

  /*!
   * \brief Inserts new element at the end of the Array.
   *
   * \param [in] args the arguments to forward to constructor of the element.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note The size increases by 1.
   *
   * \pre DIM == 1
   */
  template <typename... Args>
  void emplace_back(Args&&... args);

  /*!
   * \brief Push a value to the back of the array.
   *
   * \param [in] value the value to move to the back.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note If used in a device kernel, the number of push_backs must not exceed
   *  the capacity, since device-side reallocations aren't supported.
   * \note Array must be allocated in unified memory if calling on the device.
   *
   * \pre DIM == 1
   */
  /// @{
  AXOM_HOST_DEVICE void push_back_device(const T& value);
  AXOM_HOST_DEVICE void push_back_device(T&& value);
  /// @}

  /*!
   * \brief Inserts new element at the end of the Array.
   *
   * \param [in] args the arguments to forward to constructor of the element.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note The size increases by 1.
   * \note If used in a device kernel, the number of push_backs must not exceed
   *  the capacity, since device-side reallocations aren't supported.
   * \note Array must be allocated in unified memory if calling on the device.
   *
   * \pre DIM == 1
   */
  template <typename... Args>
  AXOM_HOST_DEVICE void emplace_back_device(Args&&... args);

  /// @}

  /// \name Array methods to query and set attributes
  /// @{

  /*!
   * \brief Return the number of elements allocated for the data array.
   */
  IndexType capacity() const { return m_capacity; }

  /*!
   * \brief Increase the capacity. Does nothing if the new capacity is less
   *  than the current capacity.
   *
   * \param [in] capacity the new number of elements to allocate.
   */
  void reserve(IndexType capacity)
  {
    if(capacity > m_capacity)
    {
      setCapacity(capacity);
    }
  }

  /*!
   * \brief Returns an ArrayIterator to the first element of the Array
   */
  ArrayIterator begin() { return ArrayIterator(0, this); }

  /// \overload
  ConstArrayIterator begin() const { return ConstArrayIterator(0, this); }

  /*!
   * \brief Returns an ArrayIterator to the element following the last
   *  element of the Array.
   */
  ArrayIterator end() { return ArrayIterator(size(), this); }

  /// \overload
  ConstArrayIterator end() const { return ConstArrayIterator(size(), this); }

  /*!
   * \brief Returns a reference to the first element in the Array.
   *
   * \pre array.empty() == false
   */
  T& front()
  {
    assert(!empty());
    return *begin();
  }

  /// \overload
  const T& front() const
  {
    assert(!empty());
    return *begin();
  }

  /*!
   * \brief Returns a reference to the last element in the Array.
   *
   * \pre array.size() > 0
   */
  T& back()
  {
    assert(!empty());
    return *(end() - 1);
  }

  /// \overload
  const T& back() const
  {
    assert(!empty());
    return *(end() - 1);
  }

  /*!
   * \brief Shrink the capacity to be equal to the size.
   */
  void shrink() { setCapacity(m_num_elements); }

  /*!
   * \brief Returns true iff the Array stores no elements.
   *
   * \note If the Array is empty the capacity can still be greater than zero.
   */
  AXOM_HOST_DEVICE inline bool empty() const { return m_num_elements == 0; }

  /*!
   * \brief Return the number of elements stored in the data array.
   */
  AXOM_HOST_DEVICE inline IndexType size() const { return m_num_elements; }

  /*!
   * \brief Update the number of elements stored in the data array.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  template <typename... Args,
            typename Enable = typename std::enable_if<
              sizeof...(Args) == DIM && detail::all_types_are_integral<Args...>::value>::type>
  void resize(Args... args)
  {
    static_assert(std::is_default_constructible<T>::value,
                  "Cannot call Array<T>::resize() when T is non-trivially-"
                  "constructible. Use Array<T>::reserve() and emplace_back()"
                  "instead.");
    const StackArray<IndexType, DIM> dims {{static_cast<IndexType>(args)...}};
    resizeImpl(dims, true);
  }

  /// \overload
  template <typename... Args,
            typename Enable = typename std::enable_if<
              sizeof...(Args) == DIM && detail::all_types_are_integral<Args...>::value>::type>
  void resize(ArrayOptions::Uninitialized, Args... args)
  {
    const StackArray<IndexType, DIM> dims {{static_cast<IndexType>(args)...}};
    resizeImpl(dims, false);
  }

  template <int Dims = DIM, typename Enable = std::enable_if_t<Dims == 1>>
  void resize(IndexType size, const T& value)
  {
    resizeImpl({{size}}, true, &value);
  }

  void resize(const StackArray<IndexType, DIM>& size, const T& value)
  {
    resizeImpl(size, true, &value);
  }

  /*!
   * \brief Exchanges the contents of this Array with the other.
   */
  void swap(Array& other);

  /*!
   * \brief Get the ratio by which the capacity increases upon dynamic resize.
   */
  double getResizeRatio() const { return m_resize_ratio; }

  /*!
   * \brief Set the ratio by which the capacity increases upon dynamic resize.
   *
   * \param [in] ratio the new resize ratio.
   */
  void setResizeRatio(double ratio) { m_resize_ratio = ratio; }

  /*!
   * \brief Get the ID for the umpire allocator
   */
  AXOM_HOST_DEVICE int getAllocatorID() const { return m_allocator_id; }

  /*!
   * \brief Sets the preferred space where operations on this array should be
   *  performed.
   *
   *  This option only has an effect for memory which is both accessible on the
   *  CPU and the GPU. For CUDA this is the Unified and Pinned memory spaces,
   *  while for HIP this is the Unified, Pinned, and Device memory spaces.
   */
  void setDevicePreference(bool on_device) { m_executeOnGPU = on_device; }

  /*!
   * \brief Returns a view of the array
   * \sa ArrayView
   */
  ArrayViewType view() { return ArrayViewType(*this); }
  /// \overload
  ConstArrayViewType view() const { return ConstArrayViewType(*this); }

  /// @}

protected:
  /*!
   * \brief Initialize an Array instance with the given number of elements.
   *
   * \param [in] num_elements the number of elements the Array holds.
   * \param [in] capacity the number of elements to allocate space for.
   * \param [in] should_default_construct whether to create default-constructed
   *  objects in the region [0, num_elements). Defaults to true.
   *
   * \note If no capacity or capacity less than num_elements is specified
   *  then it will default to at least num_elements * DEFAULT_RESIZE_RATIO.
   * \note a capacity is specified for the number of elements to store in the
   *  array and does not correspond to the actual bytesize.
   *
   * \pre num_elements >= 0
   *
   * \post capacity() >= size()
   * \post size() == num_elements
   * \post getResizeRatio() == DEFAULT_RESIZE_RATIO
   */
  void initialize(IndexType num_elements, IndexType capacity, bool should_default_construct = true);

  /*!
   * \brief Helper function for initializing an Array instance with an existing
   *  range of elements.
   *
   * \param [in] data pointer to the existing array of elements
   * \param [in] num_elements the number of elements in the existing array
   * \param [in] data_space the memory space in which data has been allocated
   * \param [in] user_provided_allocator true if the Array's allocator ID was
   *  provided by the user
   */
  void initialize_from_other(const T* data,
                             IndexType num_elements,
                             MemorySpace data_space,
                             bool user_provided_allocator);

  /*!
   * \brief Updates the number of elements stored in the data array.
   *
   * \param [in] dims the number of elements to allocate in each dimension
   * \param [in] construct_with_values if true, sets new elements in the array
   *             to a specified value
   * \param [in] value pointer to the value to fill new elements in the array
   *             with. If null, will default-construct elements in place.
   */
  void resizeImpl(const StackArray<IndexType, DIM>& dims,
                  bool construct_with_values,
                  const T* value = nullptr);

  /*!
   * \brief Make space for a subsequent insertion into the array.
   *
   * \param [in] n the number of elements to insert.
   * \param [in] pos the position at which to begin the insertion.
   *
   * \return a pointer to the beginning of the insertion space.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  T* reserveForInsert(IndexType n, IndexType pos);

  /*!
   * \brief Make space for a subsequent insertion into the array.
   *
   * \note This version supports concurrent GPU insertions.
   * \note Reallocation is not supported on the GPU.
   */
  AXOM_HOST_DEVICE IndexType reserveForPushBack();

  /*!
   * \brief Update the number of elements.
   *
   * \param [in] new_num_elements the new number of elements.
   */
  void updateNumElements(IndexType new_num_elements);

  /*!
   * \brief Set the number of elements allocated for the data array.
   *
   * \param [in] capacity the new number of elements to allocate.
   */
  void setCapacity(IndexType new_capacity);

  /*!
   * \brief Reallocates the data array when the size exceeds the capacity.
   *
   * \param [in] new_num_elements the number of elements which exceeds the
   *  current capacity.
   */
  void dynamicRealloc(IndexType new_num_elements);

  T* m_data = nullptr;
  /// \brief The full number of elements in the array
  ///  i.e., 3 for a 1D Array of size 3, 9 for a 3x3 2D array, etc
  IndexType m_num_elements = 0;
  IndexType m_capacity = 0;
  double m_resize_ratio = DEFAULT_RESIZE_RATIO;
  int m_allocator_id = INVALID_ALLOCATOR_ID;
  bool m_executeOnGPU = false;
  OpHelper m_arrayOps;
};

/// \brief Helper alias for multi-component arrays
template <typename T>
using MCArray = Array<T, 2>;

//------------------------------------------------------------------------------
//                            Array IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
Array<T, DIM, SPACE, StoragePolicy>::Array()
  : m_allocator_id(axom::detail::getAllocatorID<SPACE>())
  , m_executeOnGPU(axom::isDeviceAllocator(m_allocator_id))
  , m_arrayOps(m_allocator_id, m_executeOnGPU)
{ }

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
Array<T, DIM, SPACE, StoragePolicy>::Array(const axom::StackArray<axom::IndexType, DIM>& shape,
                                           int allocator_id)
  : ArrayBase<T, DIM, Array<T, DIM, SPACE, StoragePolicy>>(shape)
  , m_allocator_id(allocator_id)
  , m_arrayOps(m_allocator_id, m_executeOnGPU)
{
  initialize(detail::packProduct(shape.m_data), detail::packProduct(shape.m_data), false);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
Array<T, DIM, SPACE, StoragePolicy>::Array(const axom::StackArray<axom::IndexType, DIM>& shape,
                                           axom::ArrayStrideOrder rowOrColumn,
                                           int allocator_id)
  : ArrayBase<T, DIM, Array<T, DIM, SPACE, StoragePolicy>>(shape,
                                                           MDMapping<DIM> {shape, rowOrColumn, 1})
  , m_allocator_id(allocator_id)
  , m_arrayOps(m_allocator_id, m_executeOnGPU)
{
  assert(rowOrColumn == axom::ArrayStrideOrder::ROW || rowOrColumn == axom::ArrayStrideOrder::COLUMN ||
         (DIM == 1 && rowOrColumn == axom::ArrayStrideOrder::BOTH));
  initialize(detail::packProduct(shape.m_data), detail::packProduct(shape.m_data), false);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
template <typename DirType>
Array<T, DIM, SPACE, StoragePolicy>::Array(const axom::StackArray<axom::IndexType, DIM>& shape,
                                           const axom::StackArray<DirType, DIM>& slowestDirs,
                                           int allocator_id)
  : ArrayBase<T, DIM, Array<T, DIM, SPACE, StoragePolicy>>(shape, {shape, slowestDirs, 1})
  , m_allocator_id(allocator_id)
  , m_arrayOps(m_allocator_id, m_executeOnGPU)
{
  initialize(detail::packProduct(shape.m_data), detail::packProduct(shape.m_data), false);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
template <typename... Args, typename Enable>
Array<T, DIM, SPACE, StoragePolicy>::Array(Args... args)
  : ArrayBase<T, DIM, Array<T, DIM, SPACE, StoragePolicy>>(
      StackArray<IndexType, DIM> {{static_cast<IndexType>(args)...}})
  , m_allocator_id(axom::detail::getAllocatorID<SPACE>())
  , m_arrayOps(m_allocator_id, m_executeOnGPU)
{
  static_assert(sizeof...(Args) == DIM, "Array size must match number of dimensions");
  // Intel hits internal compiler error when casting as part of function call
  const IndexType tmp_args[] = {static_cast<IndexType>(args)...};
  assert(detail::allNonNegative(tmp_args));
  initialize(detail::packProduct(tmp_args), 0);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
template <typename... Args, typename Enable>
Array<T, DIM, SPACE, StoragePolicy>::Array(ArrayOptions::Uninitialized, Args... args)
  : ArrayBase<T, DIM, Array<T, DIM, SPACE, StoragePolicy>>(
      StackArray<IndexType, DIM> {{static_cast<IndexType>(args)...}})
  , m_allocator_id(axom::detail::getAllocatorID<SPACE>())
  , m_arrayOps(m_allocator_id, m_executeOnGPU)
{
  static_assert(sizeof...(Args) == DIM, "Array size must match number of dimensions");
  // Intel hits internal compiler error when casting as part of function call
  const IndexType tmp_args[] = {static_cast<IndexType>(args)...};
  assert(detail::allNonNegative(tmp_args));
  initialize(detail::packProduct(tmp_args), 0, false);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
template <IndexType SFINAE_DIM, MemorySpace SFINAE_SPACE, typename std::enable_if<SFINAE_DIM == 1>::type*>
Array<T, DIM, SPACE, StoragePolicy>::Array(IndexType num_elements, IndexType capacity, int allocator_id)
  : m_allocator_id(allocator_id)
  , m_arrayOps(m_allocator_id, m_executeOnGPU)
{
  // If a memory space has been explicitly set for the Array object, check that
  // the space of the user-provided allocator matches the explicit space.
  if(SPACE != MemorySpace::Dynamic && SPACE != axom::detail::getAllocatorSpace(m_allocator_id))
  {
#ifdef AXOM_DEBUG
    std::cerr << "Incorrect allocator ID was provided for an Array object with "
                 "explicit memory space - using default for space\n";
#endif
    m_allocator_id = axom::detail::getAllocatorID<SPACE>();
  }
  m_arrayOps = OpHelper {m_allocator_id, m_executeOnGPU};
  initialize(num_elements, capacity);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
template <IndexType SFINAE_DIM, MemorySpace SFINAE_SPACE, typename std::enable_if<SFINAE_DIM == 1>::type*>
Array<T, DIM, SPACE, StoragePolicy>::Array(ArrayOptions::Uninitialized,
                                           IndexType num_elements,
                                           IndexType capacity,
                                           int allocator_id)
  : m_allocator_id(allocator_id)
  , m_arrayOps(m_allocator_id, m_executeOnGPU)
{
  // If a memory space has been explicitly set for the Array object, check that
  // the space of the user-provided allocator matches the explicit space.
  if(SPACE != MemorySpace::Dynamic && SPACE != axom::detail::getAllocatorSpace(m_allocator_id))
  {
#ifdef AXOM_DEBUG
    std::cerr << "Incorrect allocator ID was provided for an Array object with "
                 "explicit memory space - using default for space\n";
#endif
    m_allocator_id = axom::detail::getAllocatorID<SPACE>();
  }
  m_arrayOps = OpHelper {m_allocator_id, m_executeOnGPU};
  initialize(num_elements, capacity, false);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
template <int UDIM, typename Enable>
Array<T, DIM, SPACE, StoragePolicy>::Array(std::initializer_list<T> elems, int allocator_id)
  : m_allocator_id(allocator_id)
  , m_arrayOps(m_allocator_id, m_executeOnGPU)
{
  initialize_from_other(elems.begin(), elems.size(), MemorySpace::Dynamic, true);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
AXOM_HOST_DEVICE Array<T, DIM, SPACE, StoragePolicy>::Array(const Array& other)
  : ArrayBase<T, DIM, Array<T, DIM, SPACE, StoragePolicy>>(
      static_cast<const ArrayBase<T, DIM, Array<T, DIM, SPACE, StoragePolicy>>&>(other))
  , m_allocator_id(other.m_allocator_id)
  , m_arrayOps(other.m_arrayOps)
{
#if defined(AXOM_DEVICE_CODE)
  #if defined(AXOM_DEBUG)
  printf(
    "axom::Array: cannot copy-construct on the device.\n"
    "This is usually the result of capturing an array by-value in a lambda. "
    "Use axom::ArrayView for value captures instead.\n");
  #endif
  #if defined(__CUDA_ARCH__)
  asm("trap;");
  #endif
  #if defined(__HIP_DEVICE_COMPILE__)
  abort();
  #endif
#else
  this->setCapacity(other.capacity());
  m_executeOnGPU = axom::isDeviceAllocator(m_allocator_id);
  m_arrayOps = OpHelper {m_allocator_id, m_executeOnGPU};
  // Use fill_range to ensure that copy constructors are invoked for each
  // element.
  MemorySpace srcSpace = SPACE;
  if(srcSpace == MemorySpace::Dynamic)
  {
    srcSpace = axom::detail::getAllocatorSpace(other.m_allocator_id);
  }
  m_arrayOps.fill_range(m_data, 0, other.size(), other.data(), srcSpace);
  this->updateNumElements(other.size());
#endif
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
Array<T, DIM, SPACE, StoragePolicy>::Array(Array&& other) noexcept
  : ArrayBase<T, DIM, Array<T, DIM, SPACE, StoragePolicy>>(
      static_cast<ArrayBase<T, DIM, Array<T, DIM, SPACE, StoragePolicy>>&&>(std::move(other)))
  , m_resize_ratio(0.0)
  , m_arrayOps(other.m_arrayOps)
{
  m_data = other.m_data;
  m_num_elements = other.m_num_elements;
  m_capacity = other.m_capacity;
  m_resize_ratio = other.m_resize_ratio;
  m_allocator_id = other.m_allocator_id;
  m_executeOnGPU = axom::isDeviceAllocator(m_allocator_id);
  m_arrayOps = OpHelper {m_allocator_id, m_executeOnGPU};

  other.m_data = nullptr;
  other.m_num_elements = 0;
  other.m_capacity = 0;
  other.m_resize_ratio = DEFAULT_RESIZE_RATIO;
  other.m_allocator_id = INVALID_ALLOCATOR_ID;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
template <typename OtherArrayType>
Array<T, DIM, SPACE, StoragePolicy>::Array(const ArrayBase<T, DIM, OtherArrayType>& other)
  : ArrayBase<T, DIM, Array<T, DIM, SPACE, StoragePolicy>>(other)
  , m_allocator_id(static_cast<const OtherArrayType&>(other).getAllocatorID())
  , m_arrayOps(m_allocator_id, m_executeOnGPU)
{
  initialize_from_other(static_cast<const OtherArrayType&>(other).data(),
                        static_cast<const OtherArrayType&>(other).size(),
                        axom::detail::getAllocatorSpace(m_allocator_id),
                        false);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
template <typename OtherArrayType>
Array<T, DIM, SPACE, StoragePolicy>::Array(const ArrayBase<const T, DIM, OtherArrayType>& other)
  : ArrayBase<T, DIM, Array<T, DIM, SPACE, StoragePolicy>>(other)
  , m_allocator_id(static_cast<const OtherArrayType&>(other).getAllocatorID())
  , m_arrayOps(m_allocator_id, m_executeOnGPU)
{
  initialize_from_other(static_cast<const OtherArrayType&>(other).data(),
                        static_cast<const OtherArrayType&>(other).size(),
                        axom::detail::getAllocatorSpace(m_allocator_id),
                        false);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
template <typename OtherArrayType>
Array<T, DIM, SPACE, StoragePolicy>::Array(const ArrayBase<T, DIM, OtherArrayType>& other,
                                           int allocatorId)
  : ArrayBase<T, DIM, Array<T, DIM, SPACE, StoragePolicy>>(other)
  , m_allocator_id(allocatorId)
  , m_arrayOps(m_allocator_id, m_executeOnGPU)
{
  int src_allocator = static_cast<const OtherArrayType&>(other).getAllocatorID();

  initialize_from_other(static_cast<const OtherArrayType&>(other).data(),
                        static_cast<const OtherArrayType&>(other).size(),
                        axom::detail::getAllocatorSpace(src_allocator),
                        true);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
template <typename OtherArrayType>
Array<T, DIM, SPACE, StoragePolicy>::Array(const ArrayBase<const T, DIM, OtherArrayType>& other,
                                           int allocatorId)
  : ArrayBase<T, DIM, Array<T, DIM, SPACE, StoragePolicy>>(other)
  , m_allocator_id(allocatorId)
  , m_arrayOps(m_allocator_id, m_executeOnGPU)
{
  int src_allocator = static_cast<const OtherArrayType&>(other).getAllocatorID();

  initialize_from_other(static_cast<const OtherArrayType&>(other).data(),
                        static_cast<const OtherArrayType&>(other).size(),
                        axom::detail::getAllocatorSpace(src_allocator),
                        true);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
Array<T, DIM, SPACE, StoragePolicy>::~Array()
{
  clear();
  StoragePolicy::deallocate(m_data);

  m_data = nullptr;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline void Array<T, DIM, SPACE, StoragePolicy>::fill(const T& value)
{
  m_arrayOps.destroy(m_data, 0, m_num_elements);
  m_arrayOps.fill(m_data, 0, m_num_elements, value);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline void Array<T, DIM, SPACE, StoragePolicy>::fill(const T& value, IndexType n, IndexType pos)
{
  assert(pos >= 0);
  assert(pos + n <= m_num_elements);

  m_arrayOps.destroy(m_data, pos, n);
  m_arrayOps.fill(m_data, pos, n, value);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline void Array<T, DIM, SPACE, StoragePolicy>::set(const T* elements, IndexType n, IndexType pos)
{
  assert(elements != nullptr);
  assert(pos >= 0);
  assert(pos + n <= m_num_elements);

  m_arrayOps.destroy(m_data, pos, n);
  m_arrayOps.fill_range(m_data, pos, n, elements, MemorySpace::Dynamic);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline void Array<T, DIM, SPACE, StoragePolicy>::assign(axom::IndexType count, const T& value)
{
  assert(count >= 0);
  resize(count, value);
  m_arrayOps.destroy(m_data, 0, count);
  m_arrayOps.fill(m_data, 0, count, value);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
template <class InputIt>
inline void Array<T, DIM, SPACE, StoragePolicy>::assign(InputIt first, InputIt last)
{
  Array<T, DIM, axom::MemorySpace::Dynamic> tmp;
  for(auto it = first; it != last; it++)
  {
    tmp.push_back(*it);
  }
  initialize_from_other(tmp.data(), tmp.size(), MemorySpace::Dynamic, true);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline void Array<T, DIM, SPACE, StoragePolicy>::assign(std::initializer_list<T> elems)
{
  resize(elems.size());
  m_arrayOps.destroy(m_data, 0, elems.size());
  m_arrayOps.fill_range(m_data, 0, elems.size(), elems.begin(), MemorySpace::Dynamic);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline void Array<T, DIM, SPACE, StoragePolicy>::clear()
{
  if(m_num_elements > 0)
  {
    m_arrayOps.destroy(m_data, 0, m_num_elements);

    updateNumElements(0);
  }
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline void Array<T, DIM, SPACE, StoragePolicy>::insert(IndexType pos, const T& value)
{
  static_assert(DIM == 1, "Insertion not supported for multidimensional Arrays");
  reserveForInsert(1, pos);

  m_arrayOps.emplace(m_data, pos, value);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline typename Array<T, DIM, SPACE, StoragePolicy>::ArrayIterator
Array<T, DIM, SPACE, StoragePolicy>::insert(Array<T, DIM, SPACE, StoragePolicy>::ArrayIterator pos,
                                            const T& value)
{
  static_assert(DIM == 1, "Insertion not supported for multidimensional Arrays");
  assert(pos >= begin() && pos <= end());
  insert(pos - begin(), value);
  return pos;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline void Array<T, DIM, SPACE, StoragePolicy>::insert(IndexType pos, IndexType n, const T* values)
{
  assert(values != nullptr);
  reserveForInsert(n, pos);
  m_arrayOps.fill_range(m_data, pos, n, values, MemorySpace::Dynamic);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline typename Array<T, DIM, SPACE, StoragePolicy>::ArrayIterator
Array<T, DIM, SPACE, StoragePolicy>::insert(Array<T, DIM, SPACE, StoragePolicy>::ArrayIterator pos,
                                            IndexType n,
                                            const T* values)
{
  static_assert(DIM == 1, "Insertion not supported for multidimensional Arrays");
  assert(pos >= begin() && pos <= end());
  insert(pos - begin(), n, values);
  return pos;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline void Array<T, DIM, SPACE, StoragePolicy>::insert(IndexType pos, IndexType n, const T& value)
{
  static_assert(DIM == 1, "Insertion not supported for multidimensional Arrays");
  reserveForInsert(n, pos);
  m_arrayOps.fill(m_data, pos, n, value);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline typename Array<T, DIM, SPACE, StoragePolicy>::ArrayIterator
Array<T, DIM, SPACE, StoragePolicy>::insert(Array<T, DIM, SPACE, StoragePolicy>::ArrayIterator pos,
                                            IndexType n,
                                            const T& value)
{
  static_assert(DIM == 1, "Insertion not supported for multidimensional Arrays");
  assert(pos >= begin() && pos <= end());
  insert(pos - begin(), n, value);
  return pos;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
template <MemorySpace OtherSpace>
inline void Array<T, DIM, SPACE, StoragePolicy>::insert(IndexType pos,
                                                        ArrayView<const T, DIM, OtherSpace> other)
{
  // First update the dimensions
  this->updateShapeOnInsert(other.shape());
  // Then add the raw data to the buffer
  insert(pos, other.size(), other.data());
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline typename Array<T, DIM, SPACE, StoragePolicy>::ArrayIterator
Array<T, DIM, SPACE, StoragePolicy>::erase(Array<T, DIM, SPACE, StoragePolicy>::ArrayIterator pos)
{
  assert(pos >= begin() && pos < end());

  IndexType posIdx = pos - begin();

  // Destroy element at posIdx and shift elements over by 1
  m_arrayOps.destroy(m_data, posIdx, 1);
  m_arrayOps.move(m_data, posIdx + 1, m_num_elements, posIdx);
  updateNumElements(m_num_elements - 1);

  return ArrayIterator(posIdx, this);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline typename Array<T, DIM, SPACE, StoragePolicy>::ArrayIterator
Array<T, DIM, SPACE, StoragePolicy>::erase(Array<T, DIM, SPACE, StoragePolicy>::ArrayIterator first,
                                           Array<T, DIM, SPACE, StoragePolicy>::ArrayIterator last)
{
  assert(first >= begin() && first < end());
  assert(last >= first && last <= end());

  // Empty range, return last
  if(first == last)
  {
    return last;
  }

  // Erase [first,last) elements
  IndexType firstIdx = first - begin();
  IndexType lastIdx = last - begin();
  IndexType nelems = last - first;
  m_arrayOps.destroy(m_data, firstIdx, nelems);

  // Shift [last, end) elements over
  m_arrayOps.move(m_data, lastIdx, m_num_elements, firstIdx);

  IndexType count = lastIdx - firstIdx;
  updateNumElements(m_num_elements - count);
  return ArrayIterator(firstIdx, this);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
template <typename... Args>
inline void Array<T, DIM, SPACE, StoragePolicy>::emplace(IndexType pos, Args&&... args)
{
  reserveForInsert(1, pos);
  m_arrayOps.emplace(m_data, pos, std::forward<Args>(args)...);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
template <typename... Args>
inline typename Array<T, DIM, SPACE, StoragePolicy>::ArrayIterator
Array<T, DIM, SPACE, StoragePolicy>::emplace(Array<T, DIM, SPACE, StoragePolicy>::ArrayIterator pos,
                                             Args&&... args)
{
  assert(pos >= begin() && pos <= end());
  emplace(pos - begin(), std::forward<Args>(args)...);
  return pos;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline void Array<T, DIM, SPACE, StoragePolicy>::push_back(const T& value)
{
  static_assert(DIM == 1, "push_back is only supported for 1D arrays");
  emplace_back(value);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline void Array<T, DIM, SPACE, StoragePolicy>::push_back(T&& value)
{
  static_assert(DIM == 1, "push_back is only supported for 1D arrays");
  emplace_back(std::move(value));
}

//------------------------------------------------------------------------------
AXOM_SUPPRESS_HD_WARN
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
AXOM_HOST_DEVICE inline void Array<T, DIM, SPACE, StoragePolicy>::push_back_device(const T& value)
{
  static_assert(DIM == 1, "push_back_device is only supported for 1D arrays");
  emplace_back_device(value);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
AXOM_HOST_DEVICE inline void Array<T, DIM, SPACE, StoragePolicy>::push_back_device(T&& value)
{
  static_assert(DIM == 1, "push_back_device is only supported for 1D arrays");
  emplace_back_device(std::move(value));
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
template <typename... Args>
inline void Array<T, DIM, SPACE, StoragePolicy>::emplace_back(Args&&... args)
{
  static_assert(DIM == 1, "emplace_back is only supported for 1D arrays");
  IndexType insertIndex = reserveForPushBack();
  m_arrayOps.emplace(m_data, insertIndex, std::forward<Args>(args)...);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
template <typename... Args>
AXOM_HOST_DEVICE inline void Array<T, DIM, SPACE, StoragePolicy>::emplace_back_device(Args&&... args)
{
  static_assert(DIM == 1, "emplace_back is only supported for 1D arrays");
  IndexType insertIndex = reserveForPushBack();
#ifdef AXOM_DEVICE_CODE
  // Construct in-place in uninitialized memory.
  new(m_data + insertIndex) T(std::forward<Args>(args)...);
#else
  m_arrayOps.emplace(m_data, insertIndex, std::forward<Args>(args)...);
#endif
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline void Array<T, DIM, SPACE, StoragePolicy>::resizeImpl(const StackArray<IndexType, DIM>& dims,
                                                            bool construct_with_values,
                                                            const T* value)
{
  assert(detail::allNonNegative(dims.m_data));
  const auto new_num_elements = detail::packProduct(dims.m_data);

  static_cast<ArrayBase<T, DIM, Array<T, DIM, SPACE, StoragePolicy>>&>(*this) =
    ArrayBase<T, DIM, Array<T, DIM, SPACE, StoragePolicy>> {dims};

  const IndexType prev_num_elements = m_num_elements;

  if(new_num_elements > m_capacity)
  {
    dynamicRealloc(new_num_elements);
  }

  if(prev_num_elements < new_num_elements && construct_with_values)
  {
    if(value)
    {
      // Copy-construct new elements with value
      m_arrayOps.fill(m_data, prev_num_elements, new_num_elements - prev_num_elements, *value);
    }
    else
    {
      // Default-initialize the new elements
      m_arrayOps.init(m_data, prev_num_elements, new_num_elements - prev_num_elements);
    }
  }
  else if(prev_num_elements > new_num_elements)
  {
    // Destroy any elements above new_num_elements
    m_arrayOps.destroy(m_data, new_num_elements, prev_num_elements - new_num_elements);
  }

  updateNumElements(new_num_elements);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline void Array<T, DIM, SPACE, StoragePolicy>::swap(Array<T, DIM, SPACE, StoragePolicy>& other)
{
  ArrayBase<T, DIM, Array<T, DIM, SPACE, StoragePolicy>>::swap(other);
  axom::utilities::swap(m_data, other.m_data);
  axom::utilities::swap(m_num_elements, other.m_num_elements);
  axom::utilities::swap(m_capacity, other.m_capacity);
  axom::utilities::swap(m_resize_ratio, other.m_resize_ratio);
  axom::utilities::swap(m_allocator_id, other.m_allocator_id);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline void Array<T, DIM, SPACE, StoragePolicy>::initialize(IndexType num_elements,
                                                            IndexType capacity,
                                                            bool default_construct)
{
  assert(num_elements >= 0);

  if(capacity < 0 || num_elements > capacity)
  {
    capacity = 0;
  }

  if(capacity == 0)
  {
    capacity = (num_elements > MIN_DEFAULT_CAPACITY) ? num_elements : MIN_DEFAULT_CAPACITY;
  }
  m_executeOnGPU = axom::isDeviceAllocator(m_allocator_id);
  setCapacity(capacity);
  if(default_construct)
  {
    m_arrayOps.init(m_data, 0, num_elements);
  }
  updateNumElements(num_elements);

  // quick checks
  assert(m_data != nullptr);
  assert(m_num_elements >= 0);
  assert(m_capacity >= m_num_elements);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline void Array<T, DIM, SPACE, StoragePolicy>::initialize_from_other(
  const T* other_data,
  IndexType num_elements,
  MemorySpace other_data_space,
  bool AXOM_DEBUG_PARAM(user_provided_allocator))
{
  // If a memory space has been explicitly set for the Array object, check that
  // the space of the user-provided allocator matches the explicit space.
  if(SPACE != MemorySpace::Dynamic && SPACE != axom::detail::getAllocatorSpace(m_allocator_id))
  {
#ifdef AXOM_DEBUG
    if(user_provided_allocator)
    {
      std::cerr << "Incorrect allocator ID was provided for an Array object "
                   "with explicit memory space - using default for space\n";
    }
#endif
    m_allocator_id = axom::detail::getAllocatorID<SPACE>();
  }
  m_executeOnGPU = axom::isDeviceAllocator(m_allocator_id);
  m_arrayOps = OpHelper {m_allocator_id, m_executeOnGPU};
  this->setCapacity(num_elements);
  // Use fill_range to ensure that copy constructors are invoked for each
  // element.
  m_arrayOps.fill_range(m_data, 0, num_elements, other_data, other_data_space);
  this->updateNumElements(num_elements);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline T* Array<T, DIM, SPACE, StoragePolicy>::reserveForInsert(IndexType n, IndexType pos)
{
  assert(n >= 0);
  assert(pos >= 0);
  assert(pos <= m_num_elements);

  if(n == 0)
  {
    return m_data + pos;
  }

  IndexType new_size = m_num_elements + n;
  if(new_size > m_capacity)
  {
    dynamicRealloc(new_size);
  }

  m_arrayOps.move(m_data, pos, m_num_elements, pos + n);

  updateNumElements(new_size);
  return m_data + pos;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
AXOM_HOST_DEVICE inline IndexType Array<T, DIM, SPACE, StoragePolicy>::reserveForPushBack()
{
#ifndef AXOM_DEVICE_CODE
  if AXOM_UNLIKELY(m_num_elements >= m_capacity)
  {
    dynamicRealloc(m_num_elements + 1);
  }
  IndexType end = m_num_elements;
  updateNumElements(m_num_elements + 1);
  return end;
#else
  // Device path: supports insertion while m_num_elements < m_capacity
  // Does not support insertions which require reallocating the underlying
  // buffer.
  IndexType new_pos = axom::atomicAdd<axom::auto_atomic>(&m_num_elements, IndexType {1});
  if(new_pos >= m_capacity)
  {
  #ifdef AXOM_DEBUG
    printf(
      "Array::reserveForInsert: size() exceeded capacity() when inserting "
      "on the device.\n");
  #endif
  #ifdef AXOM_USE_CUDA
    assert(false);
  #elif defined(AXOM_USE_HIP)
    abort();
  #endif
  }
  return new_pos;
#endif
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline void Array<T, DIM, SPACE, StoragePolicy>::updateNumElements(IndexType new_num_elements)
{
  assert(new_num_elements >= 0);
  assert(new_num_elements <= m_capacity);

  m_num_elements = new_num_elements;
  // Needed for Sidre array
  StoragePolicy::onShapeUpdate(this->shape());
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline void Array<T, DIM, SPACE, StoragePolicy>::setCapacity(IndexType new_capacity)
{
  assert(new_capacity >= 0);

  if(new_capacity < m_num_elements)
  {
    // No need to default-initialize here because this is only
    // when the array is being shrunk
    updateNumElements(new_capacity);
  }
  T* new_data =
    StoragePolicy::reallocate(m_data, m_num_elements, m_allocator_id, new_capacity, [this](T* new_data) {
      // Call helper method to move underlying elements if T is non-trivial.
      m_arrayOps.realloc_move(new_data, static_cast<IndexType>(this->m_num_elements), this->m_data);
    });

  if(new_data)
  {
    this->m_data = new_data;
    this->m_capacity = new_capacity;
  }

  // Needed for Sidre array
  StoragePolicy::onShapeUpdate(this->shape());

  assert(m_data != nullptr || m_capacity <= 0);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
inline void Array<T, DIM, SPACE, StoragePolicy>::dynamicRealloc(IndexType new_num_elements)
{
  assert(m_resize_ratio >= 1.0);

  // Using resize strategy from LLVM libc++ (vector::__recommend()):
  //   new_capacity = max(capacity() * resize_ratio, new_num_elements)
  IndexType new_capacity =
    axom::utilities::max<IndexType>(this->capacity() * m_resize_ratio + 0.5, new_num_elements);
  const IndexType block_size = this->blockSize();
  const IndexType remainder = new_capacity % block_size;
  if(remainder != 0)
  {
    new_capacity += block_size - remainder;
  }

  if(m_resize_ratio < 1.0)
  {
    std::cerr << "ERROR: resize ratio must be greater than 1.0.\n";
    std::cerr << "Set a valid resize ratio via calling setResizeRatio() with "
              << "an appropriate value.\n";

    utilities::processAbort();
  }

  setCapacity(new_capacity);

  assert(m_data != nullptr || m_capacity <= 0);
}

} /* namespace axom */

#endif /* AXOM_ARRAY_HPP_ */
