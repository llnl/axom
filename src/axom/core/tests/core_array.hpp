// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/core/execution/for_all.hpp"

#include "gtest/gtest.h"

#include <algorithm>

namespace
{
/*!
 * \brief Calculate the new capacity for an Array given an increase in the size.
 * \param [in] v, the Array in question.
 * \param [in] increase, the amount the size will increase by
 * \return the new capacity.
 */
template <typename T>
axom::IndexType calc_new_capacity(axom::Array<T>& v, axom::IndexType increase)
{
  axom::IndexType new_num_elements = v.size() + increase;
  if(new_num_elements > v.capacity())
  {
    axom::IndexType capacity_expanded = v.capacity() * v.getResizeRatio() + 0.5;
    return axom::utilities::max<axom::IndexType>(capacity_expanded, new_num_elements);
  }

  return v.capacity();
}

/*!
 * \brief Check if two Arrays are copies. Does not check the resize ratio.
 * \param [in] lhs, the first Array to compare.
 * \param [in] rhs, the second Array to compare.
 * \return the new capacity.
 */
template <typename T>
void check_copy(const axom::Array<T>& lhs, const axom::Array<T>& rhs)
{
  EXPECT_EQ(lhs.size(), rhs.size());
  EXPECT_EQ(lhs.capacity(), rhs.capacity());

  const T* lhs_data = lhs.data();
  const T* rhs_data = rhs.data();
  EXPECT_EQ(lhs_data, rhs_data);
}

/*!
 * \brief Check that the storage of an Array is working properly.
 * \param [in] v the Array to check.
 */
template <typename T>
void check_storage(axom::Array<T>& v)
{
  EXPECT_TRUE(v.empty());
  EXPECT_EQ(v.size(), 0);

  axom::IndexType capacity = v.capacity();
  const T* data_ptr = v.data();

  /* Push back up to half the capacity. */
  for(T i = 0; i < capacity / 2; ++i)
  {
    v.push_back(i);
  }

  /* Check the array metadata. */
  EXPECT_TRUE(!v.empty());
  EXPECT_EQ(v.size(), capacity / 2);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.data(), data_ptr);

  /* Push back up to the full capacity. */
  for(T i = capacity / 2; i < capacity; ++i)
  {
    v.push_back(i);
  }

  /* Check the array metadata. */
  EXPECT_TRUE(!v.empty());
  EXPECT_EQ(v.size(), capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.data(), data_ptr);

  /* Check the array data using the [] operator and the raw pointer. */
  for(axom::IndexType i = 0; i < capacity; ++i)
  {
    EXPECT_EQ(v[i], i);
    EXPECT_EQ(data_ptr[i], i);
  }

  /* Set the array data to new values using the [] operator. */
  for(axom::IndexType i = 0; i < capacity; ++i)
  {
    v[i] = i - 5 * i + 7;
  }

  /* Check the array data using the [] operator and the raw pointer. */
  for(axom::IndexType i = 0; i < capacity; ++i)
  {
    EXPECT_EQ(v[i], i - 5 * i + 7);
    EXPECT_EQ(data_ptr[i], i - 5 * i + 7);
  }

  /* Check the array metadata. */
  EXPECT_TRUE(!v.empty());
  EXPECT_EQ(v.size(), capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.data(), data_ptr);
}

/*!
 * \brief Check that the fill method is working properly.
 * \param [in] v the Array to check.
 */
template <typename T, int DIM, axom::MemorySpace SPACE>
void check_fill(axom::Array<T, DIM, SPACE>& v)
{
  constexpr T MAGIC_NUM_0 = 55;
  constexpr T MAGIC_NUM_1 = 6834;
  const axom::IndexType capacity = v.capacity();
  const axom::IndexType size = v.size();
  const double ratio = v.getResizeRatio();
  const T* const data_ptr = v.data();

  /* Fill the Array with MAGIC_NUM_0. */
  v.fill(MAGIC_NUM_0);

  /* Check the meta data. */
  EXPECT_EQ(capacity, v.capacity());
  EXPECT_EQ(size, v.size());
  EXPECT_EQ(ratio, v.getResizeRatio());
  EXPECT_EQ(data_ptr, v.data());

  // To check entries, we copy data to a dynamic array
  int host_alloc_id = axom::getDefaultAllocatorID();
  axom::Array<T, DIM> v_host(v, host_alloc_id);

  /* Check that the entries are all MAGIC_NUM_0. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    EXPECT_EQ(v_host[i], MAGIC_NUM_0);
  }

  /* Fill the Array with MAGIC_NUM_1. */
  v.fill(MAGIC_NUM_1);

  /* Check the meta data. */
  EXPECT_EQ(capacity, v.capacity());
  EXPECT_EQ(size, v.size());
  EXPECT_EQ(ratio, v.getResizeRatio());
  EXPECT_EQ(data_ptr, v.data());

  v_host = axom::Array<T, DIM>(v, host_alloc_id);

  /* Check that the entries are all MAGIC_NUM_1. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    EXPECT_EQ(v_host[i], MAGIC_NUM_1);
  }

  //--------------------------------------------------------
  // Test view fill. All elements contain MAGIC_NUM_1 now.

  /* Fill the Array View with MAGIC_NUM_0. */
  v.view().fill(MAGIC_NUM_0);

  v_host = axom::Array<T, DIM>(v, host_alloc_id);

  /* Check that the entries are all MAGIC_NUM_0. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    EXPECT_EQ(v_host[i], MAGIC_NUM_0);
  }
}

/*!
 * \brief Check that the set method is working properly.
 * \param [in] v the Array to check.
 */
template <typename T>
void check_set(axom::Array<T>& v)
{
  constexpr T ZERO = 0;
  const axom::IndexType capacity = v.capacity();
  const axom::IndexType size = v.size();
  const double ratio = v.getResizeRatio();
  const T* const data_ptr = v.data();

  /* Allocate a buffer half the size of the array. Fill it up with sequential values. */
  const axom::IndexType buffer_size = size / 2;
  T* buffer = axom::allocate<T>(buffer_size);
  for(axom::IndexType i = 0; i < buffer_size; ++i)
  {
    buffer[i] = i;
  }

  /* Set all the values in the array to zero. */
  v.fill(ZERO);

  /* Set the first half of the elements in the array to the sequential values in buffer. */
  v.set(buffer, buffer_size, 0);

  /* Check the array metadata. */
  EXPECT_EQ(capacity, v.capacity());
  EXPECT_EQ(size, v.size());
  EXPECT_EQ(ratio, v.getResizeRatio());
  EXPECT_EQ(data_ptr, v.data());

  /* Check that the first half of the elements in the array are equivalent to those in buffer. */
  for(axom::IndexType i = 0; i < buffer_size; ++i)
  {
    EXPECT_EQ(v[i], buffer[i]);
  }

  /* Check that the second half of the elements in the array are all zero. */
  for(axom::IndexType i = buffer_size; i < size; ++i)
  {
    EXPECT_EQ(v[i], ZERO);
  }

  /* Reset the values in buffer to the next sequential values. */
  for(axom::IndexType i = 0; i < buffer_size; ++i)
  {
    buffer[i] = i + buffer_size;
  }

  /* Set the second half of the elements in the array to the new sequential values in buffer. */
  v.set(buffer, buffer_size, buffer_size);

  /* Check the array metadata. */
  EXPECT_EQ(capacity, v.capacity());
  EXPECT_EQ(size, v.size());
  EXPECT_EQ(ratio, v.getResizeRatio());
  EXPECT_EQ(data_ptr, v.data());

  /* Check that all the elements in the array now hold sequential values. */
  for(axom::IndexType i = 0; i < 2 * buffer_size; ++i)
  {
    EXPECT_EQ(v[i], i);
  }

  // Check setting through the view.
  v.fill(ZERO);
  v.view().set(buffer, buffer_size, 0);
  for(axom::IndexType i = 0; i < buffer_size; ++i)
  {
    EXPECT_EQ(v[i], buffer[i]);
  }

  axom::deallocate(buffer);
}

/*!
 * \brief Check that the assign method is working properly.
 * \param [in] v the Array to check.
 */
template <typename ArrayType>
void check_assign(ArrayType& v)
{
  using T = typename ArrayType::value_type;
  constexpr int DIM = 1;
  constexpr T MAGIC_NUM_0 = 55;
  const axom::IndexType size = v.size();

  // Fill the Array with MAGIC_NUM_0. This may be filling device memory.
  v.assign(size, MAGIC_NUM_0);
  EXPECT_EQ(size, v.size());

  // Copy data onto host and compare.
  axom::Array<T, DIM, axom::MemorySpace::Dynamic> vHost(size, size);
  axom::copy(vHost.data(), v.data(), size * sizeof(T));
  for(axom::IndexType i = 0; i < size; i++)
  {
    EXPECT_EQ(vHost[i], MAGIC_NUM_0);
  }

  // Try the iterator method using source data on host.
  std::set<T> s;
  for(axom::IndexType i = 0; i < size; i++)
  {
    s.insert(static_cast<T>(i));
  }
  v.assign(s.begin(), s.end());
  axom::copy(vHost.data(), v.data(), v.size() * sizeof(T));
  for(axom::IndexType i = 0; i < size; i++)
  {
    EXPECT_EQ(vHost[i], static_cast<T>(i));
  }

  // Try initializer list using source data on host.
  v.assign({T {0}, T {11}, T {22}, T {33}, T {44}, T {55}, T {66}, T {77}, T {88}, T {99}});
  const axom::IndexType smallSize = 10;
  EXPECT_EQ(smallSize, v.size());
  axom::copy(vHost.data(), v.data(), v.size() * sizeof(T));
  for(axom::IndexType i = 0; i < smallSize; i++)
  {
    EXPECT_EQ(vHost[i], static_cast<T>(i * 11));
  }
}

/*!
 * \brief Check that the resizing of an array is working properly.
 * \param [in] v the Array to check.
 */
template <typename T>
void check_resize(axom::Array<T>& v)
{
  /* Resize the array up to the capacity */
  axom::IndexType capacity = v.capacity();
  v.resize(capacity);
  axom::IndexType size = capacity;

  /* Check that the size equals the capacity. */
  EXPECT_EQ(v.size(), v.capacity());

  /* Set the existing data in v */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    v[i] = static_cast<T>(i - 5 * i + 7);
  }

  /* Push back a new element, should resize. */
  axom::IndexType old_capacity = capacity;
  capacity = calc_new_capacity(v, 1);
  v.push_back(size - 5 * size + 7);
  size++;

  /* Check that it resized properly */
  EXPECT_GT(capacity, old_capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check that the data is still intact. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    EXPECT_EQ(v[i], i - 5 * i + 7);
  }

  /* Push back 1000 elements */
  const axom::IndexType n_elements = 1000;

  T* values = axom::allocate<T>(n_elements);
  for(axom::IndexType i = 0; i < n_elements; ++i)
  {
    axom::IndexType i_real = i + size;
    values[i] = i_real - 5 * i_real + 7;
  }

  /* Push back the new elements. */
  capacity = calc_new_capacity(v, n_elements);
  v.insert(size, n_elements, values);
  size += n_elements;

  /* Check that size and capacity are as expected. */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check that the data is still intact. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    EXPECT_EQ(v[i], i - 5 * i + 7);
  }

  /* Reduce the size down to 500 elements */
  T* data_address = v.data();
  size = 500;
  v.resize(size);

  /* Check the metadata. */
  EXPECT_EQ(v.size(), size);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.data(), data_address);

  /* Check the data. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    EXPECT_EQ(v[i], i - 5 * i + 7);
  }

  /* Shrink the vector */
  capacity = size;
  v.shrink();

  /* Check that the capacity and size are as expected. */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check that the data is intact. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    EXPECT_EQ(v[i], i - 5 * i + 7);
  }

  /* Push back a new element, should resize. */
  old_capacity = capacity;
  capacity = calc_new_capacity(v, 1);
  v.push_back(size - 5 * size + 7);
  size++;

  /* Check the new size and capacity. */
  EXPECT_GT(capacity, old_capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check that the data is intact. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    EXPECT_EQ(v[i], i - 5 * i + 7);
  }

  /* Reset the data */
  T* data_ptr = v.data();
  for(axom::IndexType i = 0; i < size; ++i)
  {
    data_ptr[i] = i;
  }

  /* Push back a bunch of elements to fill in up to the capacity. Resize should
   * not occur. */
  old_capacity = capacity;
  for(axom::IndexType i = size; i < old_capacity; ++i)
  {
    v.push_back(i);
    size++;
    EXPECT_EQ(v.capacity(), old_capacity);
    EXPECT_EQ(v.size(), size);
    EXPECT_EQ(v.data(), data_ptr);
  }

  EXPECT_EQ(v.size(), old_capacity);

  /* Push back a final element that should trigger a resize. */
  capacity = calc_new_capacity(v, old_capacity - size + 1);
  v.push_back(size);
  size++;

  /* Check the new capacity and size. */
  EXPECT_GT(capacity, old_capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check the data. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    EXPECT_EQ(v[i], i);
  }

  axom::deallocate(values);
  values = nullptr;
}

/*!
 * \brief Check that the insertion into an Array is working properly.
 * \param [in] v the Array to check.
 */
template <typename T>
void check_insert(axom::Array<T>& v)
{
  /* Resize the array up to the capacity */
  axom::IndexType capacity = v.capacity();
  v.resize(capacity);
  axom::IndexType size = capacity;

  EXPECT_EQ(v.size(), v.capacity());

  /* Set the existing data in v */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    v[i] = i - 5 * i + 7;
  }

  /* Insert a new element, should resize. */
  axom::IndexType old_capacity = capacity;
  capacity = calc_new_capacity(v, 1);
  v.insert(v.size(), 1, size - 5 * size + 7);
  size++;

  /* Check that it resized properly */
  EXPECT_GT(capacity, old_capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);
  for(axom::IndexType i = 0; i < size; ++i)
  {
    EXPECT_EQ(v[i], i - 5 * i + 7);
  }

  /* Insert 1000 elements */
  const axom::IndexType n_elements = 1000;
  T* values = axom::allocate<T>(n_elements);
  for(axom::IndexType i = 0; i < n_elements; ++i)
  {
    axom::IndexType i_real = i + size;
    values[i] = i_real - 5 * i_real + 7;
  }

  capacity = calc_new_capacity(v, n_elements);
  v.insert(size, n_elements, values);
  size += n_elements;

  /* Check that it resizes properly */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);
  for(axom::IndexType i = 0; i < size; ++i)
  {
    EXPECT_EQ(v[i], i - 5 * i + 7);
  }

  capacity = size;
  v.shrink();
  axom::IndexType n_insert_front = 100;

  /* Reset the data */
  T* data_ptr = v.data();
  for(axom::IndexType i = 0; i < size; ++i)
  {
    data_ptr[i] = i + n_insert_front;
  }

  /* Insert into the front of the array. */
  for(axom::IndexType i = n_insert_front - 1; i >= 0; i--)
  {
    capacity = calc_new_capacity(v, 1);
    v.insert(0, 1, i);
    size++;
  }

  /* Check that the insertion worked as expected */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);
  for(axom::IndexType i = 0; i < size; ++i)
  {
    EXPECT_EQ(v[i], i);
  }

  axom::deallocate(values);
  values = nullptr;
}

/*!
 * \brief Check that the insertion into an Array is working properly for iterators.
 * \param [in] v the Array to check.
 */
template <typename T>
void check_insert_iterator(axom::Array<T>& v)
{
  /* Resize the array up to the capacity */
  axom::IndexType capacity = v.capacity();
  v.resize(capacity);
  axom::IndexType size = capacity;

  EXPECT_EQ(v.size(), v.capacity());

  /* Set the existing data in v */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    v[i] = i - 5 * i + 7;
  }

  /* Insert a new element, should resize. */
  axom::IndexType old_capacity = capacity;
  capacity = calc_new_capacity(v, 1);
  auto ret = v.insert(v.end(), 1, size - 5 * size + 7);
  size++;

  /* Check that it resized properly */
  EXPECT_GT(capacity, old_capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);
  EXPECT_EQ(ret, v.end() - 1);
  for(axom::IndexType i = 0; i < size; ++i)
  {
    EXPECT_EQ(v[i], i - 5 * i + 7);
  }

  /* Insert 1000 elements */
  const axom::IndexType n_elements = 1000;
  T* values = axom::allocate<T>(n_elements);
  for(axom::IndexType i = 0; i < n_elements; ++i)
  {
    axom::IndexType i_real = i + size;
    values[i] = i_real - 5 * i_real + 7;
  }

  capacity = calc_new_capacity(v, n_elements);
  auto ret2 = v.insert(v.end(), n_elements, values);

  EXPECT_EQ(ret2, v.begin() + size);
  size += n_elements;

  /* Check that it resizes properly */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);
  for(axom::IndexType i = 0; i < size; ++i)
  {
    EXPECT_EQ(v[i], i - 5 * i + 7);
  }

  capacity = size;
  v.shrink();
  axom::IndexType n_insert_front = 100;

  /* Reset the data */
  T* data_ptr = v.data();
  for(axom::IndexType i = 0; i < size; ++i)
  {
    data_ptr[i] = i + n_insert_front;
  }

  /* Insert into the front of the Array. */
  for(axom::IndexType i = n_insert_front - 1; i >= 0; i--)
  {
    capacity = calc_new_capacity(v, 1);
    auto ret3 = v.insert(v.begin(), 1, i);
    EXPECT_EQ(ret3, v.begin());
    EXPECT_EQ(i, v.front());
    size++;
  }

  /* Check that the insertion worked as expected */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);
  for(axom::IndexType i = 0; i < size; ++i)
  {
    EXPECT_EQ(v[i], i);
  }

  axom::deallocate(values);
  values = nullptr;
}

/*!
 * \brief Check that emplace() into an Array is working properly.
 * \param [in] v the Array to check.
 */
template <typename T>
void check_emplace(axom::Array<T>& v)
{
  /* Resize the array up to the capacity */
  axom::IndexType capacity = v.capacity();
  v.resize(capacity);
  axom::IndexType size = capacity;

  EXPECT_EQ(v.size(), v.capacity());

  /* Set the existing data in v */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    v[i] = i - 5 * i + 7;
  }

  /* Emplace a new element, should resize. */
  axom::IndexType old_capacity = capacity;
  capacity = calc_new_capacity(v, 1);
  auto ret = v.emplace(v.end(), size - 5 * size + 7);
  size++;

  /* Check that it resized properly */
  EXPECT_GT(capacity, old_capacity);
  EXPECT_EQ(v.size(), size);
  EXPECT_EQ(ret, v.end() - 1);
  for(axom::IndexType i = 0; i < size; ++i)
  {
    EXPECT_EQ(v[i], i - 5 * i + 7);
  }

  /* Emplace_back 1000 elements */
  const axom::IndexType n_elements = 1000;
  for(axom::IndexType i = 0; i < n_elements; ++i)
  {
    axom::IndexType i_real = i + size;
    v.emplace_back(i_real - 5 * i_real + 7);
  }

  size += n_elements;

  /* Check that it resizes properly */
  EXPECT_EQ(v.size(), size);
  for(axom::IndexType i = 0; i < size; ++i)
  {
    EXPECT_EQ(v[i], i - 5 * i + 7);
  }

  capacity = size;
  v.shrink();
  axom::IndexType n_insert_front = 100;

  /* Reset the data */
  T* data_ptr = v.data();
  for(axom::IndexType i = 0; i < size; ++i)
  {
    data_ptr[i] = i + n_insert_front;
  }

  /* Emplace into the front of the Array. */
  for(axom::IndexType i = n_insert_front - 1; i >= 0; i--)
  {
    capacity = calc_new_capacity(v, 1);
    auto ret3 = v.emplace(v.begin(), i);
    EXPECT_EQ(ret3, v.begin());
    EXPECT_EQ(i, v.front());
    size++;
  }

  /* Check that the emplace worked as expected */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);
  for(axom::IndexType i = 0; i < size; ++i)
  {
    EXPECT_EQ(v[i], i);
  }
}

// Clumsy implementation of https://numpy.org/doc/stable/reference/generated/numpy.zeros_like.html
// for testing purposes
template <typename T, int DIM>
axom::Array<T, DIM> zeros_like(const axom::Array<T, DIM>& other)
{
  axom::Array<T, DIM> zeros(other);
  zeros.fill(T {});
  return zeros;
}

template <typename T, int DIM>
void check_swap(axom::Array<T, DIM>& v)
{
  auto v_two = zeros_like(v);

  /* Push 0...size elements */
  for(int i = 0; i < v.size(); i++)
  {
    v.flatIndex(i) = i;
    v_two.flatIndex(i) = -i;
  }

  /* Create copies */
  axom::Array<T, DIM> v_copy(v);
  axom::Array<T, DIM> v_two_copy(v_two);

  EXPECT_EQ(v, v_copy);
  EXPECT_EQ(v_two, v_two_copy);

  /* Swap */
  v.swap(v_two);

  EXPECT_NE(v, v_two);
  EXPECT_NE(v, v_copy);

  /* Swap back */
  v.swap(v_two);

  EXPECT_EQ(v, v_copy);
  EXPECT_EQ(v_two, v_two_copy);
}

template <typename T, int DIM, axom::MemorySpace SPACE>
void check_alloc(axom::Array<T, DIM, SPACE>& v, const int id)
{
  // Verify allocation
  EXPECT_EQ(v.getAllocatorID(), id);
  EXPECT_EQ(v.size(), v.capacity());

// Use introspection to verify pointer if available
#if defined(AXOM_USE_UMPIRE)
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  auto v_allocator = rm.getAllocator(v.data());
  EXPECT_EQ(v_allocator.getId(), id);
#endif
}

/*!
 * \brief Check an external array for defects.
 * \param [in] v the external array to check.
 */
template <typename T>
void check_external_view(axom::ArrayView<T>& v)
{
  const axom::IndexType size = v.size();
  const axom::IndexType num_values = size;
  T* const data_ptr = v.data();

  /* Set the elements in the array. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    v[i] = i;
  }

  /* Check the elements using the raw pointer. */
  for(axom::IndexType i = 0; i < num_values; ++i)
  {
    EXPECT_EQ(data_ptr[i], i);
  }

  /* Set the elements using the raw pointer. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    data_ptr[i] = i * i;
  }

  /* Check the elements using the () operator. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    EXPECT_EQ(v[i], i * i);
  }

  EXPECT_EQ(size, v.size());
  EXPECT_EQ(data_ptr, v.data());
}

#if defined(AXOM_USE_GPU) && defined(AXOM_GPUCC) && defined(AXOM_USE_UMPIRE)

template <typename T>
__global__ void assign_raw(T* data, int N)
{
  for(int i = 0; i < N; i++)
  {
    data[i] = i;
  }
}

template <typename T, int DIM, axom::MemorySpace SPACE>
__global__ void assign_view(axom::ArrayView<T, DIM, SPACE> view)
{
  for(int i = 0; i < view.size(); i++)
  {
    view[i] = i * 2;
  }
}

/*!
 * \brief Check that an array can be modified/accessed from device code
 * \param [in] v the array to check.
 */
template <typename T, int DIM, axom::MemorySpace SPACE>
void check_device(axom::Array<T, DIM, SPACE>& v)
{
  const axom::IndexType size = v.size();
  // Then assign to it via a raw device pointer
  assign_raw<<<1, 1>>>(v.data(), size);

  // Check the contents of the array by assigning to a Dynamic array
  // The Umpire allocator should be Host, so we can access it from the CPU
  int host_alloc_id = axom::getDefaultAllocatorID();
  axom::Array<T, 1> check_raw_array_dynamic(v, host_alloc_id);
  EXPECT_EQ(check_raw_array_dynamic.size(), size);
  for(int i = 0; i < check_raw_array_dynamic.size(); i++)
  {
    EXPECT_EQ(check_raw_array_dynamic[i], i);
  }

  // Then check the contents by assigning to an explicitly Host Array
  axom::Array<T, 1, axom::MemorySpace::Host> check_raw_array_host = v;
  EXPECT_EQ(check_raw_array_host.size(), size);
  for(int i = 0; i < check_raw_array_host.size(); i++)
  {
    EXPECT_EQ(check_raw_array_host[i], i);
  }

  // Then modify the underlying data via a view
  axom::ArrayView<T, DIM, SPACE> view(v);
  assign_view<<<1, 1>>>(view);

  // Check the contents of the array by assigning to a Dynamic array
  // The Umpire allocator should be Host, so we can access it from the CPU
  axom::Array<T, 1> check_view_array_dynamic(view, host_alloc_id);
  EXPECT_EQ(check_view_array_dynamic.size(), size);
  for(int i = 0; i < check_view_array_dynamic.size(); i++)
  {
    EXPECT_EQ(check_view_array_dynamic[i], i * 2);
  }

  // Then check the contents by assigning to an explicitly Host array
  axom::Array<T, 1, axom::MemorySpace::Host> check_view_array_host = view;
  EXPECT_EQ(check_view_array_host.size(), size);
  for(int i = 0; i < check_view_array_host.size(); i++)
  {
    EXPECT_EQ(check_view_array_host[i], i * 2);
  }
}

template <typename T>
__global__ void assign_raw_2d(T* data, int M, int N)
{
  for(int i = 0; i < M; i++)
  {
    for(int j = 0; j < N; j++)
    {
      data[i * N + j] = i * i + j;
    }
  }
}

template <typename T, axom::MemorySpace SPACE>
__global__ void assign_view_2d(axom::ArrayView<T, 2, SPACE> view)
{
  for(int i = 0; i < view.shape()[0]; i++)
  {
    for(int j = 0; j < view.shape()[1]; j++)
    {
      view(i, j) = j * j + i;
    }
  }
}

/*!
 * \brief Check that a 2D array can be modified/accessed from device code
 * \param [in] v the array to check.
 */
template <typename T, axom::MemorySpace SPACE>
void check_device_2D(axom::Array<T, 2, SPACE>& v)
{
  const axom::IndexType size = v.size();
  const axom::IndexType M = v.shape()[0];
  const axom::IndexType N = v.shape()[1];
  // Then assign to it via a raw device pointer
  assign_raw_2d<<<1, 1>>>(v.data(), M, N);

  // Check the contents of the array by assigning to a Dynamic array
  // The Umpire allocator should be Host, so we can access it from the CPU
  int host_alloc_id = axom::getDefaultAllocatorID();
  axom::Array<T, 2> check_raw_array_dynamic(v, host_alloc_id);
  EXPECT_EQ(check_raw_array_dynamic.size(), size);
  EXPECT_EQ(check_raw_array_dynamic.shape(), v.shape());

  for(int i = 0; i < M; i++)
  {
    for(int j = 0; j < N; j++)
    {
      EXPECT_EQ(check_raw_array_dynamic(i, j), i * i + j);
    }
  }

  // Then check the contents by assigning to an explicitly Host array
  axom::Array<T, 2, axom::MemorySpace::Host> check_raw_array_host = v;
  EXPECT_EQ(check_raw_array_host.size(), size);
  EXPECT_EQ(check_raw_array_host.shape(), v.shape());

  for(int i = 0; i < M; i++)
  {
    for(int j = 0; j < N; j++)
    {
      EXPECT_EQ(check_raw_array_host(i, j), i * i + j);
    }
  }

  // Then modify the underlying data via a view
  axom::ArrayView<T, 2, SPACE> view(v);
  assign_view_2d<<<1, 1>>>(view);

  // Check the contents of the array by assigning to a Dynamic array
  // The Umpire allocator should be Host, so we can access it from the CPU
  axom::Array<T, 2> check_view_array_dynamic(view, host_alloc_id);
  EXPECT_EQ(check_view_array_dynamic.size(), size);
  EXPECT_EQ(check_view_array_dynamic.shape(), v.shape());

  for(int i = 0; i < M; i++)
  {
    for(int j = 0; j < N; j++)
    {
      EXPECT_EQ(check_view_array_dynamic(i, j), j * j + i);
    }
  }

  // Then check the contents by assigning to an explicitly Host array
  axom::Array<T, 2, axom::MemorySpace::Host> check_view_array_host = view;
  EXPECT_EQ(check_view_array_host.size(), size);
  EXPECT_EQ(check_view_array_host.shape(), v.shape());

  for(int i = 0; i < M; i++)
  {
    for(int j = 0; j < N; j++)
    {
      EXPECT_EQ(check_view_array_host(i, j), j * j + i);
    }
  }
}

#endif  // defined(AXOM_USE_GPU) && defined(AXOM_GPUCC) && defined(AXOM_USE_UMPIRE)

}  // end anonymous namespace

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

TEST(core_array, checkSlowestDirsConstructor)
{
  axom::Array<double, 1> a1({3}, {{0}});
  EXPECT_TRUE(a1.size() == 3);
  EXPECT_TRUE(a1.mapping().fastestStrideLength() == 1);
  EXPECT_TRUE(a1.mapping().slowestDirs()[0] == 0);
  EXPECT_TRUE(a1.mapping().strides()[0] == 1);

  axom::Array<double, 2> a2({3, 4}, {1, 0});
  EXPECT_TRUE(a2.size() == 12);
  EXPECT_TRUE(a2.mapping().fastestStrideLength() == 1);
  EXPECT_TRUE(a2.mapping().slowestDirs()[0] == 1);
  EXPECT_TRUE(a2.mapping().slowestDirs()[1] == 0);
  EXPECT_TRUE(a2.mapping().strides()[0] == 1);
  EXPECT_TRUE(a2.mapping().strides()[1] == 3);

  axom::Array<double, 3> a3({3, 4, 5}, {2, 0, 1});
  EXPECT_TRUE(a3.size() == 60);
  EXPECT_TRUE(a3.mapping().fastestStrideLength() == 1);
  EXPECT_TRUE(a3.mapping().slowestDirs()[0] == 2);
  EXPECT_TRUE(a3.mapping().slowestDirs()[1] == 0);
  EXPECT_TRUE(a3.mapping().slowestDirs()[2] == 1);
  EXPECT_TRUE(a3.mapping().strides()[0] == 4);
  EXPECT_TRUE(a3.mapping().strides()[1] == 1);
  EXPECT_TRUE(a3.mapping().strides()[2] == 12);
}

//------------------------------------------------------------------------------
TEST(core_array, checkStorage)
{
  constexpr axom::IndexType ZERO = 0;

  for(axom::IndexType capacity = 2; capacity < 512; capacity *= 2)
  {
    axom::Array<int> v_int(ZERO, capacity);
    ::check_storage(v_int);

    axom::Array<double> v_double(ZERO, capacity);
    ::check_storage(v_double);
  }
}

//------------------------------------------------------------------------------
TEST(core_array, checkFill)
{
  for(axom::IndexType capacity = 2; capacity < 512; capacity *= 2)
  {
    axom::IndexType size = capacity / 2;
    axom::Array<int> v_int(size, capacity);
    ::check_fill(v_int);

    axom::Array<double> v_double(size, capacity);
    ::check_fill(v_double);
  }
}

//------------------------------------------------------------------------------
#if defined(AXOM_USE_GPU) && defined(AXOM_GPUCC) && defined(AXOM_USE_UMPIRE)
TEST(core_array, checkFillDevice)
{
  for(axom::IndexType capacity = 2; capacity < 512; capacity *= 2)
  {
    axom::IndexType size = capacity / 2;
    axom::Array<int, 1, axom::MemorySpace::Device> v_int(size, capacity);
    ::check_fill(v_int);

    axom::Array<double, 1, axom::MemorySpace::Device> v_double(size, capacity);
    ::check_fill(v_double);
  }
}
#endif

//------------------------------------------------------------------------------
TEST(core_array, checkSet)
{
  for(axom::IndexType size = 2; size < 512; size *= 2)
  {
    axom::Array<int> v_int(size);
    ::check_set(v_int);

    axom::Array<double> v_double(size);
    ::check_set(v_double);
  }
}

//------------------------------------------------------------------------------
TEST(core_array, checkAssign)
{
  const axom::IndexType size = 100;
  axom::Array<int> v_int(size);
  ::check_assign(v_int);

  axom::Array<double> v_double(size);
  ::check_assign(v_double);
}

//------------------------------------------------------------------------------
TEST(core_array, checkAssignView)
{
  const axom::IndexType size = 100;
  axom::Array<int> v_int(size);
  auto vv_int = v_int.view();
  ::check_assign(vv_int);

  axom::Array<double> v_double(size);
  auto vv_double = v_double.view();
  ::check_assign(vv_double);
}

//------------------------------------------------------------------------------
#if defined(AXOM_USE_GPU) && defined(AXOM_GPUCC) && defined(AXOM_USE_UMPIRE)
TEST(core_array, checkAssignDevice)
{
  // Check Array::assign methods when using device memory.
  const axom::IndexType size = 100, capacity = 100;
  axom::Array<int, 1, axom::MemorySpace::Device> v_int(size, capacity);
  ::check_assign(v_int);

  axom::Array<double, 1, axom::MemorySpace::Device> v_double(size, capacity);
  ::check_assign(v_double);
}

TEST(core_array, checkAssignViewDevice)
{
  // Check Array::assign methods when using device memory.
  const axom::IndexType size = 100, capacity = 100;
  axom::Array<int, 1, axom::MemorySpace::Device> v_int(size, capacity);
  auto vv_int = v_int.view();
  ::check_assign(vv_int);

  axom::Array<double, 1, axom::MemorySpace::Device> v_double(size, capacity);
  auto vv_double = v_double.view();
  ::check_assign(vv_double);
}
#endif

//------------------------------------------------------------------------------
TEST(core_array, checkResize)
{
  constexpr axom::IndexType ZERO = 0;

  for(double ratio = 1.0; ratio <= 2.0; ratio += 0.5)
  {
    for(axom::IndexType capacity = 2; capacity <= 512; capacity *= 2)
    {
      axom::Array<int> v_int(ZERO, capacity);
      v_int.setResizeRatio(ratio);
      ::check_resize(v_int);

      axom::Array<double> v_double(ZERO, capacity);
      v_double.setResizeRatio(ratio);
      ::check_resize(v_double);
    }
  }
}

//------------------------------------------------------------------------------
TEST(core_array_DeathTest, checkResize)
{
  constexpr axom::IndexType ZERO = 0;
  axom::IndexType size = 100;

  /* Resizing isn't allowed with a ratio less than 1.0. */
  axom::Array<int> v_int(ZERO, size);
  v_int.setResizeRatio(0.99);
  EXPECT_DEATH_IF_SUPPORTED(::check_resize(v_int), "");
}

//------------------------------------------------------------------------------
TEST(core_array, checkInsert)
{
  constexpr axom::IndexType ZERO = 0;

  for(double ratio = 1.0; ratio <= 2.0; ratio += 0.5)
  {
    for(axom::IndexType capacity = 2; capacity <= 512; capacity *= 2)
    {
      axom::Array<int> v_int(ZERO, capacity);
      v_int.setResizeRatio(ratio);
      ::check_insert(v_int);

      axom::Array<double> v_double(ZERO, capacity);
      v_double.setResizeRatio(ratio);
      ::check_insert(v_double);
    }
  }
}

//------------------------------------------------------------------------------
TEST(core_array, checkInsertIterator)
{
  constexpr axom::IndexType ZERO = 0;

  for(double ratio = 1.0; ratio <= 2.0; ratio += 0.5)
  {
    for(axom::IndexType capacity = 10; capacity <= 512; capacity *= 2)
    {
      axom::Array<int> v_int(ZERO, capacity);
      v_int.setResizeRatio(ratio);
      ::check_insert_iterator(v_int);

      axom::Array<double> v_double(ZERO, capacity);
      v_double.setResizeRatio(ratio);
      ::check_insert_iterator(v_double);
    }
  }
}

//------------------------------------------------------------------------------
TEST(core_array, checkEmplace)
{
  constexpr axom::IndexType ZERO = 0;

  for(double ratio = 1.0; ratio <= 2.0; ratio += 0.5)
  {
    for(axom::IndexType capacity = 10; capacity <= 512; capacity *= 2)
    {
      axom::Array<int> v_int(ZERO, capacity);
      v_int.setResizeRatio(ratio);
      ::check_emplace(v_int);

      axom::Array<double> v_double(ZERO, capacity);
      v_double.setResizeRatio(ratio);
      ::check_emplace(v_double);
    }
  }
}

//------------------------------------------------------------------------------
TEST(core_array, checkSwap)
{
  for(axom::IndexType size = 10; size <= 512; size *= 2)
  {
    axom::Array<int> v_int(size);
    ::check_swap(v_int);

    axom::Array<double> v_double(size);
    ::check_swap(v_double);

    axom::Array<int, 2> v_int_2d(size, size);
    ::check_swap(v_int_2d);

    axom::Array<double, 2> v_double_2d(size, size);
    ::check_swap(v_double);
  }
}

//------------------------------------------------------------------------------
TEST(core_array, checkAlloc)
{
  std::vector<int> memory_locations {
#if defined(AXOM_USE_UMPIRE)
    axom::getUmpireResourceAllocatorID(umpire::resource::Host)
  #if defined(UMPIRE_ENABLE_DEVICE)
      ,
    axom::getUmpireResourceAllocatorID(umpire::resource::Device)
  #endif
  #if defined(UMPIRE_ENABLE_UM)
      ,
    axom::getUmpireResourceAllocatorID(umpire::resource::Unified)
  #endif
  #if defined(UMPIRE_ENABLE_CONST)
      ,
    axom::getUmpireResourceAllocatorID(umpire::resource::Constant)
  #endif
  #if defined(UMPIRE_ENABLE_PINNED)
      ,
    axom::getUmpireResourceAllocatorID(umpire::resource::Pinned)
  #endif
#endif
  };

  for(double ratio = 1.0; ratio <= 2.0; ratio += 0.5)
  {
    for(axom::IndexType capacity = 4; capacity <= 512; capacity *= 2)
    {
      // First use the dynamic option
      for(int id : memory_locations)
      {
        axom::Array<int, 1, axom::MemorySpace::Dynamic> v_int(capacity, capacity, id);
        ::check_alloc(v_int, id);

        axom::Array<double, 1, axom::MemorySpace::Dynamic> v_double(capacity, capacity, id);
        ::check_alloc(v_double, id);
      }
// Then, if Umpire is available, we can use the space as an explicit template parameter
#ifdef AXOM_USE_UMPIRE
  #ifdef UMPIRE_ENABLE_DEVICE
      axom::Array<int, 1, axom::MemorySpace::Device> v_int_device(capacity, capacity);
      ::check_alloc(v_int_device, axom::getUmpireResourceAllocatorID(umpire::resource::Device));
      axom::Array<double, 1, axom::MemorySpace::Device> v_double_device(capacity, capacity);
      ::check_alloc(v_double_device, axom::getUmpireResourceAllocatorID(umpire::resource::Device));
  #endif
  #ifdef UMPIRE_ENABLE_UM
      axom::Array<int, 1, axom::MemorySpace::Unified> v_int_unified(capacity, capacity);
      ::check_alloc(v_int_unified, axom::getUmpireResourceAllocatorID(umpire::resource::Unified));
      axom::Array<double, 1, axom::MemorySpace::Unified> v_double_unified(capacity, capacity);
      ::check_alloc(v_double_unified, axom::getUmpireResourceAllocatorID(umpire::resource::Unified));
  #endif
  #ifdef UMPIRE_ENABLE_CONST
      axom::Array<int, 1, axom::MemorySpace::Constant> v_int_const(capacity, capacity);
      ::check_alloc(v_int_const, axom::getUmpireResourceAllocatorID(umpire::resource::Constant));
      axom::Array<double, 1, axom::MemorySpace::Constant> v_double_const(capacity, capacity);
      ::check_alloc(v_double_const, axom::getUmpireResourceAllocatorID(umpire::resource::Constant));
  #endif
  #ifdef UMPIRE_ENABLE_PINNED
      axom::Array<int, 1, axom::MemorySpace::Pinned> v_int_pinned(capacity, capacity);
      ::check_alloc(v_int_pinned, axom::getUmpireResourceAllocatorID(umpire::resource::Pinned));
      axom::Array<double, 1, axom::MemorySpace::Pinned> v_double_pinned(capacity, capacity);
      ::check_alloc(v_double_pinned, axom::getUmpireResourceAllocatorID(umpire::resource::Pinned));
  #endif
#endif
    }
  }
}

//------------------------------------------------------------------------------
TEST(core_array, checkExternalView)
{
  constexpr double MAGIC_NUM = 5683578.8;
  constexpr axom::IndexType MAX_SIZE = 256;
  constexpr axom::IndexType MAX_VALUES = MAX_SIZE;
  union DataBuffer
  {
    int ints[MAX_SIZE];
    double doubles[MAX_SIZE];
  };

  DataBuffer buffer;
  std::fill_n(buffer.doubles, MAX_VALUES, MAGIC_NUM);

  for(axom::IndexType size = 16; size <= MAX_SIZE; size *= 2)
  {
    axom::ArrayView<int> v_int_view(buffer.ints, size);
    EXPECT_EQ(v_int_view.data(), buffer.ints);
    ::check_external_view(v_int_view);

    axom::ArrayView<double> v_double_view(buffer.doubles, size);
    EXPECT_EQ(v_double_view.data(), buffer.doubles);
    ::check_external_view(v_double_view);

    /* Set v_double's data to MAGIC_NUM */
    std::fill_n(v_double_view.data(), size, MAGIC_NUM);

    /* Check that the data still exists in the buffer */
    for(axom::IndexType i = 0; i < MAX_VALUES; ++i)
    {
      EXPECT_EQ(buffer.doubles[i], MAGIC_NUM);
    }
  }
}

//------------------------------------------------------------------------------
TEST(core_array, checkIterator)
{
  constexpr int SIZE = 1000;
  axom::Array<int> v_int(SIZE);

  /* Push 0...999 elements */
  for(int i = 0; i < SIZE; i++)
  {
    v_int[i] = i;
  }

  EXPECT_EQ(*v_int.begin(), 0);
  EXPECT_EQ(v_int.front(), 0);
  EXPECT_EQ(*(v_int.end() - 1), SIZE - 1);
  EXPECT_EQ(v_int.back(), SIZE - 1);
  EXPECT_EQ(v_int.size(), SIZE);

  /* Erase nothing */
  auto ret1 = v_int.erase(v_int.begin() + SIZE / 2, v_int.begin() + SIZE / 2);

  EXPECT_EQ(ret1, v_int.begin() + SIZE / 2);
  EXPECT_EQ(v_int.size(), SIZE);

  /* Erase half the elements */
  auto ret2 = v_int.erase(v_int.begin(), v_int.begin() + SIZE / 2);

  EXPECT_EQ(ret2, v_int.begin());
  EXPECT_EQ(*v_int.begin(), SIZE / 2);
  EXPECT_EQ(v_int.front(), SIZE / 2);
  EXPECT_EQ(*(v_int.end() - 1), SIZE - 1);
  EXPECT_EQ(v_int.back(), SIZE - 1);
  EXPECT_EQ(v_int.size(), SIZE / 2);

  /* Erase first, last elements */
  auto ret3 = v_int.erase(v_int.begin());

  EXPECT_EQ(ret3, v_int.begin());
  EXPECT_EQ(*v_int.begin(), SIZE / 2 + 1);
  EXPECT_EQ(v_int.front(), SIZE / 2 + 1);

  auto ret4 = v_int.erase(v_int.end() - 1);

  EXPECT_EQ(ret4, v_int.end());
  EXPECT_EQ(*(v_int.end() - 1), SIZE - 2);
  EXPECT_EQ(v_int.back(), SIZE - 2);

  /* Clear the rest of the array */
  v_int.clear();
  EXPECT_EQ(v_int.size(), 0);
}

//------------------------------------------------------------------------------
#if defined(AXOM_USE_GPU) && defined(AXOM_GPUCC) && defined(AXOM_USE_UMPIRE)
void checkIteratorDeviceImpl()
{
  constexpr int SIZE = 1000;
  axom::Array<int, 1, axom::MemorySpace::Host> v_int_host(SIZE);
  axom::Array<int, 1, axom::MemorySpace::Device> v_int(SIZE);

  /* Push 0...999 elements */
  for(int i = 0; i < SIZE; i++)
  {
    v_int_host[i] = i;
  }
  v_int = v_int_host;

  EXPECT_EQ(*v_int_host.begin(), 0);
  EXPECT_EQ(*(v_int_host.end() - 1), SIZE - 1);
  EXPECT_EQ(v_int.size(), SIZE);

  /* Erase nothing */
  auto ret1 = v_int.erase(v_int.begin() + SIZE / 2, v_int.begin() + SIZE / 2);

  EXPECT_EQ(ret1, v_int.begin() + SIZE / 2);
  EXPECT_EQ(v_int.size(), SIZE);

  /* Erase half the elements */
  auto ret2 = v_int.erase(v_int.begin(), v_int.begin() + SIZE / 2);

  EXPECT_EQ(ret2, v_int.begin());
  EXPECT_EQ(v_int.size(), SIZE / 2);
  v_int_host = v_int;
  EXPECT_EQ(*v_int_host.begin(), SIZE / 2);
  EXPECT_EQ(*(v_int_host.end() - 1), SIZE - 1);

  /* Erase first, last elements */
  auto ret3 = v_int.erase(v_int.begin());

  EXPECT_EQ(ret3, v_int.begin());
  v_int_host = v_int;
  EXPECT_EQ(*v_int_host.begin(), SIZE / 2 + 1);

  auto ret4 = v_int.erase(v_int.end() - 1);

  EXPECT_EQ(ret4, v_int.end());
  v_int_host = v_int;
  EXPECT_EQ(*(v_int_host.end() - 1), SIZE - 2);

  /* Clear the rest of the array */
  v_int.clear();
  EXPECT_EQ(v_int.size(), 0);
}

TEST(core_array, checkIteratorDevice) { checkIteratorDeviceImpl(); }
#endif

//------------------------------------------------------------------------------
TEST(core_array, check_move_copy)
{
  constexpr int MAGIC_INT = 255;
  constexpr double MAGIC_DOUBLE = 5683578.8;

  for(axom::IndexType capacity = 2; capacity < 512; capacity *= 2)
  {
    axom::IndexType size = capacity;

    /* Check copy and move semantics for array of ints */
    axom::Array<int> v_int(size, capacity);
    v_int.fill(MAGIC_INT);

    axom::Array<int> v_int_copy_ctor(v_int);
    axom::Array<int> v_int_copy_assign;
    v_int_copy_assign = v_int;
    EXPECT_EQ(v_int, v_int_copy_ctor);
    EXPECT_EQ(v_int, v_int_copy_assign);

    axom::Array<int> v_int_move_assign;
    v_int_move_assign = std::move(v_int_copy_assign);
    axom::Array<int> v_int_move_ctor = std::move(v_int_copy_ctor);
    EXPECT_EQ(v_int, v_int_move_assign);
    EXPECT_EQ(v_int, v_int_move_ctor);
    EXPECT_EQ(v_int_copy_assign.data(), nullptr);
    EXPECT_EQ(v_int_copy_ctor.data(), nullptr);

    /* Check copy and move semantics for array of doubles */
    axom::Array<double> v_double(size, capacity);
    v_double.fill(MAGIC_DOUBLE);

    axom::Array<double> v_double_copy_ctor(v_double);
    axom::Array<double> v_double_copy_assign;
    v_double_copy_assign = v_double;
    EXPECT_EQ(v_double, v_double_copy_ctor);
    EXPECT_EQ(v_double, v_double_copy_assign);

    axom::Array<double> v_double_move_assign;
    v_double_move_assign = std::move(v_double_copy_assign);
    axom::Array<double> v_double_move_ctor = std::move(v_double_copy_ctor);
    EXPECT_EQ(v_double, v_double_move_assign);
    EXPECT_EQ(v_double, v_double_move_ctor);
    EXPECT_EQ(v_double_copy_assign.data(), nullptr);
    EXPECT_EQ(v_double_copy_ctor.data(), nullptr);
  }
}

//------------------------------------------------------------------------------
TEST(core_array, check_move_copy_multidimensional)
{
  constexpr int MAGIC_INT = 255;
  constexpr double MAGIC_DOUBLE = 5683578.8;

  for(axom::IndexType capacity = 2; capacity < 512; capacity *= 2)
  {
    axom::IndexType size = capacity;

    /* Check copy and move semantics for array of ints */
    axom::Array<int, 2> v_int(size, size);
    v_int.fill(MAGIC_INT);

    axom::Array<int, 2> v_int_copy_ctor(v_int);
    axom::Array<int, 2> v_int_copy_assign;
    v_int_copy_assign = v_int;
    EXPECT_EQ(v_int, v_int_copy_ctor);
    EXPECT_EQ(v_int, v_int_copy_assign);
    // operator== already checks the shape, but we should also check the strides
    EXPECT_EQ(v_int.strides(), v_int_copy_ctor.strides());
    EXPECT_EQ(v_int.strides(), v_int_copy_assign.strides());

    axom::Array<int, 2> v_int_move_assign;
    v_int_move_assign = std::move(v_int_copy_assign);
    axom::Array<int, 2> v_int_move_ctor = std::move(v_int_copy_ctor);
    EXPECT_EQ(v_int, v_int_move_assign);
    EXPECT_EQ(v_int, v_int_move_ctor);
    EXPECT_EQ(v_int.strides(), v_int_move_assign.strides());
    EXPECT_EQ(v_int.strides(), v_int_move_ctor.strides());
    EXPECT_EQ(v_int_copy_assign.data(), nullptr);
    EXPECT_EQ(v_int_copy_ctor.data(), nullptr);

    /* Check copy and move semantics for array of doubles */
    axom::Array<double, 2> v_double(size, size);
    v_double.fill(MAGIC_DOUBLE);

    axom::Array<double, 2> v_double_copy_ctor(v_double);
    axom::Array<double, 2> v_double_copy_assign;
    v_double_copy_assign = v_double;
    EXPECT_EQ(v_double, v_double_copy_ctor);
    EXPECT_EQ(v_double, v_double_copy_assign);
    // operator== already checks the shape, but we should also check the strides
    EXPECT_EQ(v_double.strides(), v_double_copy_ctor.strides());
    EXPECT_EQ(v_double.strides(), v_double_copy_assign.strides());

    axom::Array<double, 2> v_double_move_assign;
    v_double_move_assign = std::move(v_double_copy_assign);
    axom::Array<double, 2> v_double_move_ctor = std::move(v_double_copy_ctor);
    EXPECT_EQ(v_double, v_double_move_assign);
    EXPECT_EQ(v_double, v_double_move_ctor);
    EXPECT_EQ(v_double.strides(), v_double_move_assign.strides());
    EXPECT_EQ(v_double.strides(), v_double_move_ctor.strides());
    EXPECT_EQ(v_double_copy_assign.data(), nullptr);
    EXPECT_EQ(v_double_copy_ctor.data(), nullptr);
  }
}

//------------------------------------------------------------------------------
TEST(core_array, check_view_move_copy)
{
  constexpr int MAGIC_INT = 255;
  constexpr double MAGIC_DOUBLE = 5683578.8;

  for(axom::IndexType capacity = 2; capacity < 512; capacity *= 2)
  {
    axom::IndexType size = capacity;

    /* Check copy and move semantics for ArrayView of ints */
    std::vector<int> ints(size, MAGIC_INT);
    axom::ArrayView<int> v_int_view(ints.data(), size);

    axom::ArrayView<int> v_int_view_copy_ctor(v_int_view);
    axom::ArrayView<int> v_int_view_copy_assign;
    v_int_view_copy_assign = v_int_view;
    EXPECT_EQ(v_int_view, v_int_view_copy_ctor);
    EXPECT_EQ(v_int_view, v_int_view_copy_assign);

    /* Check copy and move semantics for ArrayView of doubles */
    std::vector<double> doubles(size, MAGIC_DOUBLE);
    axom::ArrayView<double> v_double_view(doubles.data(), size);

    axom::ArrayView<double> v_double_view_copy_ctor(v_double_view);
    axom::ArrayView<double> v_double_view_copy_assign;
    v_double_view_copy_assign = v_double_view;
    EXPECT_EQ(v_double_view, v_double_view_copy_ctor);
    EXPECT_EQ(v_double_view, v_double_view_copy_assign);
  }
}

//------------------------------------------------------------------------------
TEST(core_array, check_multidimensional)
{
  constexpr int MAGIC_INT = 255;
  constexpr double MAGIC_DOUBLE = 5683578.8;

  // First test multidimensional int arrays
  axom::Array<int, 2> v_int;
  v_int.resize(2, 2);
  v_int.fill(MAGIC_INT);
  // Make sure the number of elements and contents are correct
  EXPECT_EQ(v_int.size(), 2 * 2);
  axom::StackArray<axom::IndexType, 2> expected_shape = {2, 2};
  EXPECT_EQ(v_int.shape(), expected_shape);
  for(const auto val : v_int)
  {
    EXPECT_EQ(val, MAGIC_INT);
  }
  // Then assign different values to each element
  v_int(0, 0) = 1;
  v_int(0, 1) = 2;
  v_int(1, 0) = 3;
  v_int(1, 1) = 4;

  // FIXME: Should we add a std::initializer_list ctor?
  axom::Array<int> v_int_flat(4);
  v_int_flat[0] = 1;
  v_int_flat[1] = 2;
  v_int_flat[2] = 3;
  v_int_flat[3] = 4;
  axom::StackArray<axom::IndexType, 1> expected_flat_shape = {4};
  EXPECT_EQ(v_int_flat.shape(), expected_flat_shape);

  for(int i = 0; i < v_int_flat.size(); i++)
  {
    // For a multidim array, flatIndex(i) is a "flat" index into the raw data
    EXPECT_EQ(v_int.flatIndex(i), v_int_flat[i]);
  }

  axom::Array<double, 3> v_double(4, 3, 2);
  v_double.fill(MAGIC_DOUBLE);
  EXPECT_EQ(v_double.size(), 4 * 3 * 2);
  axom::StackArray<axom::IndexType, 3> expected_double_shape = {4, 3, 2};
  EXPECT_EQ(v_double.shape(), expected_double_shape);
  for(const auto val : v_double)
  {
    EXPECT_EQ(val, MAGIC_DOUBLE);
  }

  axom::Array<double> v_double_flat(4 * 3 * 2);
  int double_flat_idx = 0;
  for(int i = 0; i < v_double.shape()[0]; i++)
  {
    for(int j = 0; j < v_double.shape()[1]; j++)
    {
      for(int k = 0; k < v_double.shape()[2]; k++)
      {
        v_double(i, j, k) = (i * i) + (5 * j) + k;
        v_double_flat[double_flat_idx++] = (i * i) + (5 * j) + k;
      }
    }
  }

  for(int i = 0; i < v_double.size(); i++)
  {
    // For a multidim array, flatIndex(i) is a "flat" index into the raw data
    EXPECT_EQ(v_double.flatIndex(i), v_double_flat[i]);
  }

  for(int i = 0; i < v_double.shape()[0]; i++)
  {
    axom::Array<double, 2> v_double_subarray_2d = v_double[i];
    for(int j = 0; j < v_double.shape()[1]; j++)
    {
      axom::Array<double> v_double_subarray_1d = v_double(i, j);
      for(int k = 0; k < v_double.shape()[2]; k++)
      {
        EXPECT_EQ(v_double(i, j, k), v_double[i][j][k]);
        EXPECT_EQ(v_double(i, j, k), v_double_subarray_2d(j, k));
        EXPECT_EQ(v_double(i, j, k), v_double_subarray_1d[k]);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(core_array, check_multidimensional_view)
{
  constexpr int MAGIC_INT = 255;
  constexpr double MAGIC_DOUBLE = 5683578.8;

  // First test multidimensional int arrays
  int v_int_arr[] = {MAGIC_INT, MAGIC_INT, MAGIC_INT, MAGIC_INT};
  axom::ArrayView<int, 2> v_int_view(v_int_arr, 2, 2);
  // Make sure the number of elements and contents are correct
  EXPECT_EQ(v_int_view.size(), 2 * 2);
  axom::StackArray<axom::IndexType, 2> expected_shape = {2, 2};
  EXPECT_EQ(v_int_view.shape(), expected_shape);
  for(const auto val : v_int_view)
  {
    EXPECT_EQ(val, MAGIC_INT);
  }
  // Then assign different values to each element
  v_int_view(0, 0) = 1;
  v_int_view(0, 1) = 2;
  v_int_view(1, 0) = 3;
  v_int_view(1, 1) = 4;

  // FIXME: Should we add a std::initializer_list ctor?
  int v_int_flat_arr[] = {1, 2, 3, 4};
  axom::ArrayView<int> v_int_flat_view(v_int_flat_arr, 4);
  axom::StackArray<axom::IndexType, 1> expected_flat_shape = {4};
  EXPECT_EQ(v_int_flat_view.shape(), expected_flat_shape);

  for(int i = 0; i < v_int_flat_view.size(); i++)
  {
    // For a multidim array, flatIndex(i) is a "flat" index into the raw data
    EXPECT_EQ(v_int_view.flatIndex(i), v_int_flat_view[i]);
  }

  double v_double_arr[4 * 3 * 2];
  std::fill_n(v_double_arr, 4 * 3 * 2, MAGIC_DOUBLE);
  axom::ArrayView<double, 3> v_double_view(v_double_arr, 4, 3, 2);
  EXPECT_EQ(v_double_view.size(), 4 * 3 * 2);
  axom::StackArray<axom::IndexType, 3> expected_double_shape = {4, 3, 2};
  EXPECT_EQ(v_double_view.shape(), expected_double_shape);
  for(const auto val : v_double_view)
  {
    EXPECT_EQ(val, MAGIC_DOUBLE);
  }

  double v_double_flat_arr[4 * 3 * 2];
  axom::ArrayView<double> v_double_flat_view(v_double_flat_arr, 4 * 3 * 2);
  int double_flat_idx = 0;
  for(int i = 0; i < v_double_view.shape()[0]; i++)
  {
    for(int j = 0; j < v_double_view.shape()[1]; j++)
    {
      for(int k = 0; k < v_double_view.shape()[2]; k++)
      {
        v_double_view(i, j, k) = (i * i) + (5 * j) + k;
        v_double_flat_view[double_flat_idx++] = (i * i) + (5 * j) + k;
      }
    }
  }

  for(int i = 0; i < v_double_view.size(); i++)
  {
    // For a multidim array, flatIndex(i) is a "flat" index into the raw data
    EXPECT_EQ(v_double_view.flatIndex(i), v_double_flat_view[i]);
  }

  for(int i = 0; i < v_double_view.shape()[0]; i++)
  {
    axom::Array<double, 2> v_double_subarray_2d = v_double_view[i];
    for(int j = 0; j < v_double_view.shape()[1]; j++)
    {
      axom::Array<double> v_double_subarray_1d = v_double_view(i, j);
      for(int k = 0; k < v_double_view.shape()[2]; k++)
      {
        EXPECT_EQ(v_double_view(i, j, k), v_double_view[i][j][k]);
        EXPECT_EQ(v_double_view(i, j, k), v_double_subarray_2d(j, k));
        EXPECT_EQ(v_double_view(i, j, k), v_double_subarray_1d[k]);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(core_array, check_multidimensional_view_subspan)
{
  constexpr double GHOST_DOUBLE = -1.0;

  // Create a block of 2 * 3 * 4 doubles wrapped by a layer of -1s, representing
  // a layer of "ghost zones."
  double v_double_arr[4 * 5 * 6];
  std::fill_n(v_double_arr, 4 * 5 * 6, GHOST_DOUBLE);

  // For the non-ghost values, fill them with unique values.
  constexpr int I_STRIDE = 6 * 5;
  constexpr int J_STRIDE = 6;
  constexpr int K_STRIDE = 1;
  for(int i = 1; i < 3; i++)
  {
    for(int j = 1; j < 4; j++)
    {
      for(int k = 1; k < 5; k++)
      {
        double i_comp = (i - 1) * 0.1;
        double j_comp = (j - 1) * 1.0;
        double k_comp = (k - 1) * 10.0;
        double value = i_comp + j_comp + k_comp;
        v_double_arr[i * I_STRIDE + j * J_STRIDE + k] = value;
      }
    }
  }

  axom::ArrayView<double, 3> double_view(v_double_arr, {4, 5, 6});

  EXPECT_EQ(double_view.size(), 4 * 5 * 6);
  EXPECT_EQ(double_view.strides()[0], I_STRIDE);
  EXPECT_EQ(double_view.strides()[1], J_STRIDE);
  EXPECT_EQ(double_view.strides()[2], K_STRIDE);

  // Construct a subspan view starting at (1, 1, 1) and spanning 2 x 3 x 4 elements.
  axom::ArrayView<double, 3> double_view_real = double_view.subspan({1, 1, 1}, {2, 3, 4});

  EXPECT_EQ(double_view_real.size(), 2 * 3 * 4);
  EXPECT_EQ(double_view_real.shape()[0], 2);
  EXPECT_EQ(double_view_real.shape()[1], 3);
  EXPECT_EQ(double_view_real.shape()[2], 4);
  EXPECT_EQ(double_view_real.strides(), double_view.strides());

  // Iterate over logical indexes and check for validity.
  for(int i = 0; i < 2; i++)
  {
    for(int j = 0; j < 3; j++)
    {
      for(int k = 0; k < 4; k++)
      {
        double i_comp = i * 0.1;
        double j_comp = j * 1.0;
        double k_comp = k * 10.0;
        double value = i_comp + j_comp + k_comp;
        EXPECT_EQ(double_view_real(i, j, k), value);
        EXPECT_EQ(double_view_real[i][j][k], value);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(core_array, check_multidimensional_view_subspan_colmaj)
{
  constexpr double GHOST_DOUBLE = -1.0;

  // Create a block of 2 * 3 * 4 doubles wrapped by a layer of -1s, representing
  // a layer of "ghost zones."
  double v_double_arr[4 * 5 * 6];
  std::fill_n(v_double_arr, 4 * 5 * 6, GHOST_DOUBLE);

  // For the non-ghost values, fill them with unique values.
  constexpr int I_STRIDE = 1;
  constexpr int J_STRIDE = 4;
  constexpr int K_STRIDE = 4 * 5;
  for(int i = 1; i < 3; i++)
  {
    for(int j = 1; j < 4; j++)
    {
      for(int k = 1; k < 5; k++)
      {
        double i_comp = (i - 1) * 0.1;
        double j_comp = (j - 1) * 1.0;
        double k_comp = (k - 1) * 10.0;
        double value = i_comp + j_comp + k_comp;
        v_double_arr[i * I_STRIDE + j * J_STRIDE + k * K_STRIDE] = value;
      }
    }
  }

  axom::ArrayView<double, 3> double_view(v_double_arr, {4, 5, 6}, {I_STRIDE, J_STRIDE, K_STRIDE});

  EXPECT_EQ(double_view.size(), 4 * 5 * 6);
  EXPECT_EQ(double_view.minStride(), 1);
  EXPECT_EQ(double_view.strides()[0], I_STRIDE);
  EXPECT_EQ(double_view.strides()[1], J_STRIDE);
  EXPECT_EQ(double_view.strides()[2], K_STRIDE);

  // Construct a subspan view starting at (1, 1, 1) and spanning 2 x 3 x 4 elements
  axom::ArrayView<double, 3> double_view_real = double_view.subspan({1, 1, 1}, {2, 3, 4});

  EXPECT_EQ(double_view_real.size(), 2 * 3 * 4);
  EXPECT_EQ(double_view_real.shape()[0], 2);
  EXPECT_EQ(double_view_real.shape()[1], 3);
  EXPECT_EQ(double_view_real.shape()[2], 4);
  EXPECT_EQ(double_view_real.strides(), double_view.strides());

  // Iterate over logical indexes and check for validity.
  for(int i = 0; i < 2; i++)
  {
    for(int j = 0; j < 3; j++)
    {
      for(int k = 0; k < 4; k++)
      {
        double i_comp = i * 0.1;
        double j_comp = j * 1.0;
        double k_comp = k * 10.0;
        double value = i_comp + j_comp + k_comp;
        EXPECT_EQ(double_view_real(i, j, k), value);
        EXPECT_EQ(double_view_real[i][j][k], value);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(core_array, check_multidimensional_view_subspan_stride)
{
  constexpr double GHOST_DOUBLE = -1.0;

  // Create a block of 2 * 3 * 4 doubles wrapped by a layer of -1s, representing
  // a layer of "ghost zones."
  double v_double_arr[4 * 5 * 6];
  std::fill_n(v_double_arr, 4 * 5 * 6, GHOST_DOUBLE);

  // For the non-ghost values, fill them with unique values.
  constexpr int I_STRIDE = 5 * 6;
  constexpr int J_STRIDE = 1;
  constexpr int K_STRIDE = 5;
  for(int i = 1; i < 3; i++)
  {
    for(int j = 1; j < 4; j++)
    {
      for(int k = 1; k < 5; k++)
      {
        double i_comp = (i - 1) * 0.1;
        double j_comp = (j - 1) * 1.0;
        double k_comp = (k - 1) * 10.0;
        double value = i_comp + j_comp + k_comp;
        v_double_arr[i * I_STRIDE + j * J_STRIDE + k * K_STRIDE] = value;
      }
    }
  }

  axom::ArrayView<double, 3> double_view(v_double_arr, {4, 5, 6}, {I_STRIDE, J_STRIDE, K_STRIDE});

  EXPECT_EQ(double_view.size(), 4 * 5 * 6);
  EXPECT_EQ(double_view.minStride(), 1);
  EXPECT_EQ(double_view.strides()[0], I_STRIDE);
  EXPECT_EQ(double_view.strides()[1], J_STRIDE);
  EXPECT_EQ(double_view.strides()[2], K_STRIDE);

  // Construct a subspan view starting at (1, 1, 1) and spanning 2 x 3 x 4 elements.
  axom::ArrayView<double, 3> double_view_real = double_view.subspan({1, 1, 1}, {2, 3, 4});

  EXPECT_EQ(double_view_real.size(), 2 * 3 * 4);
  EXPECT_EQ(double_view_real.shape()[0], 2);
  EXPECT_EQ(double_view_real.shape()[1], 3);
  EXPECT_EQ(double_view_real.shape()[2], 4);
  EXPECT_EQ(double_view_real.strides(), double_view.strides());

  // Iterate over logical indexes and check for validity.
  for(int i = 0; i < 2; i++)
  {
    for(int j = 0; j < 3; j++)
    {
      for(int k = 0; k < 4; k++)
      {
        double i_comp = i * 0.1;
        double j_comp = j * 1.0;
        double k_comp = k * 10.0;
        double value = i_comp + j_comp + k_comp;
        EXPECT_EQ(double_view_real(i, j, k), value);
        EXPECT_EQ(double_view_real[i][j][k], value);
      }
    }
  }
}

struct Elem
{
  double earth, wind, fire;
};
//------------------------------------------------------------------------------
TEST(core_array, check_multidimensional_view_spacing)
{
  // NUM_COMPS is the spacing between consecutive values of the
  // same type (earth, wind or fire).
  constexpr int NUM_COMPS = sizeof(Elem) / sizeof(double);
  constexpr int NUM_DIMS = 3;

  // Initialize a multidimensional, multicomponent array.
  axom::StackArray<axom::IndexType, NUM_DIMS> shape {3, 2, 4};
  axom::Array<Elem, 3> mdElemArray(shape[0], shape[1], shape[2]);
  for(int i = 0; i < shape[0]; ++i)
  {
    for(int j = 0; j < shape[1]; ++j)
    {
      for(int k = 0; k < shape[2]; ++k)
      {
        int pvalue = 1000 * i + 100 * j + 10 * k;
        mdElemArray(i, j, k).earth = pvalue + .1;
        mdElemArray(i, j, k).wind = pvalue + .2;
        mdElemArray(i, j, k).fire = pvalue + .3;
      }
    }
  }

  // Verify views of the 3 individual components in mdElemarray.
  axom::IndexType spacing = NUM_COMPS;
  axom::ArrayView<double, 3> earthOnly(&mdElemArray(0, 0, 0).earth, shape, spacing);
  axom::ArrayView<double, 3> windOnly(&mdElemArray(0, 0, 0).wind, shape, spacing);
  axom::ArrayView<double, 3> fireOnly(&mdElemArray(0, 0, 0).fire, shape, spacing);

  for(int i = 0; i < shape[0]; ++i)
  {
    for(int j = 0; j < shape[1]; ++j)
    {
      for(int k = 0; k < shape[2]; ++k)
      {
        int pvalue = 1000 * i + 100 * j + 10 * k;
        EXPECT_EQ(earthOnly(i, j, k), pvalue + .1);
        EXPECT_EQ(windOnly(i, j, k), pvalue + .2);
        EXPECT_EQ(fireOnly(i, j, k), pvalue + .3);

        EXPECT_EQ(earthOnly[i][j][k], pvalue + .1);
        EXPECT_EQ(windOnly[i][j][k], pvalue + .2);
        EXPECT_EQ(fireOnly[i][j][k], pvalue + .3);

        EXPECT_EQ(&earthOnly(i, j, k), &mdElemArray(i, j, k).earth);
        EXPECT_EQ(&windOnly(i, j, k), &mdElemArray(i, j, k).wind);
        EXPECT_EQ(&fireOnly(i, j, k), &mdElemArray(i, j, k).fire);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(core_array, checkDevice)
{
#if !defined(AXOM_USE_GPU) || !defined(AXOM_GPUCC) || !defined(AXOM_USE_UMPIRE) || \
  !defined(UMPIRE_ENABLE_DEVICE)
  GTEST_SKIP() << "CUDA or HIP is not available, skipping tests that use Array "
                  "in device code";
#else
  for(axom::IndexType capacity = 2; capacity < 512; capacity *= 2)
  {
    // Allocate a Dynamic array in Device memory
    axom::Array<int, 1, axom::MemorySpace::Dynamic> v_int_dynamic(
      capacity,
      capacity,
      axom::getUmpireResourceAllocatorID(umpire::resource::Device));

    ::check_device(v_int_dynamic);

    axom::Array<double, 1, axom::MemorySpace::Dynamic> v_double_dynamic(
      capacity,
      capacity,
      axom::getUmpireResourceAllocatorID(umpire::resource::Device));

    ::check_device(v_double_dynamic);

    // Then allocate an explicitly Device array
    axom::Array<int, 1, axom::MemorySpace::Device> v_int_device(capacity, capacity);
    ::check_device(v_int_device);

    axom::Array<double, 1, axom::MemorySpace::Device> v_double_device(capacity, capacity);
    ::check_device(v_double_device);
  }
#endif
}

//------------------------------------------------------------------------------
TEST(core_array, checkDevice2D)
{
#if !defined(AXOM_USE_GPU) || !defined(AXOM_GPUCC) || !defined(AXOM_USE_UMPIRE) || \
  !defined(UMPIRE_ENABLE_DEVICE)
  GTEST_SKIP() << "CUDA or HIP is not available, skipping tests that use Array "
                  "in device code";
#else
  for(axom::IndexType capacity = 2; capacity < 512; capacity *= 2)
  {
    // Allocate an explicitly Device array
    axom::Array<int, 2, axom::MemorySpace::Device> v_int_device(capacity, capacity);
    ::check_device_2D(v_int_device);

    axom::Array<double, 2, axom::MemorySpace::Device> v_double_device(capacity, capacity);
    ::check_device_2D(v_double_device);
  }
#endif
}

struct HasDefault
{
  int member = 255;

  bool operator==(const HasDefault& other) const { return (member == other.member); }

  friend std::ostream& operator<<(std::ostream& os, const HasDefault& hd)
  {
    return os << hd.member;
  }
};

//------------------------------------------------------------------------------
TEST(core_array, checkDefaultInitialization)
{
  constexpr int MAGIC_INT = 255;
  for(axom::IndexType capacity = 2; capacity < 512; capacity *= 2)
  {
    // Make sure the array gets zero-initialized
    axom::Array<int> v_int(capacity);
    for(const auto ele : v_int)
    {
      EXPECT_EQ(ele, 0);
    }

    axom::Array<int, 2> v_int_2d(capacity, capacity);
    for(const auto ele : v_int_2d)
    {
      EXPECT_EQ(ele, 0);
    }

    axom::Array<HasDefault> v_has_default(capacity);
    for(const auto& ele : v_has_default)
    {
      EXPECT_EQ(ele.member, MAGIC_INT);
    }

    axom::Array<HasDefault, 2> v_has_default_2d(capacity, capacity);
    for(const auto& ele : v_has_default_2d)
    {
      EXPECT_EQ(ele.member, MAGIC_INT);
    }
  }
}

//------------------------------------------------------------------------------
TEST(core_array, checkDefaultInitializationDevice)
{
#if !defined(AXOM_USE_GPU) || !defined(AXOM_GPUCC) || !defined(AXOM_USE_UMPIRE) || \
  !defined(UMPIRE_ENABLE_DEVICE)
  GTEST_SKIP() << "CUDA or HIP is not available, skipping tests that use Array "
                  "in device code";
#else
  constexpr int MAGIC_INT = 255;
  for(axom::IndexType capacity = 2; capacity < 512; capacity *= 2)
  {
    // Allocate an explicitly Device array of ints (zero-initialized)
    axom::Array<int, 1, axom::MemorySpace::Device> v_int(capacity);

    // Then copy it to the host
    axom::Array<int, 1, axom::MemorySpace::Host> v_int_host(v_int);

    for(const auto ele : v_int_host)
    {
      EXPECT_EQ(ele, 0);
    }

    // Allocate an explicitly Device array of a default-constructible type
    axom::Array<HasDefault, 1, axom::MemorySpace::Device> v_has_default_device(capacity);

    // Then copy it to the host
    axom::Array<HasDefault, 1, axom::MemorySpace::Host> v_has_default_host(v_has_default_device);

    for(const auto& ele : v_has_default_host)
    {
      EXPECT_EQ(ele.member, MAGIC_INT);
    }
  }
#endif
}

//------------------------------------------------------------------------------

/**
 * A struct that fails if the default constructor is called.
 * It is in support of the checkUninitialized test for axom::Array
 */
struct FailsOnConstruction
{
  int member = 255;

  FailsOnConstruction() { EXPECT_TRUE(false); }

  bool operator==(const FailsOnConstruction& other) const { return member == other.member; }
};

TEST(core_array, checkUninitialized)
{
  for(axom::IndexType capacity = 2; capacity < 512; capacity *= 2)
  {
    // Test uninitialized functionality with 1D Array using FailsOnConstruction type
    {
      axom::Array<FailsOnConstruction> arr(axom::ArrayOptions::Uninitialized {}, capacity);

      EXPECT_LE(capacity, arr.capacity());
      EXPECT_EQ(capacity, arr.size());

      arr.resize(axom::ArrayOptions::Uninitialized {}, capacity * 2);

      EXPECT_LE(capacity * 2, arr.capacity());
      EXPECT_EQ(capacity * 2, arr.size());
    }

    // Test default 1D Array with trivially copyable HasDefault type
    // Note: this test will not fail if HasDefault is copied, but will at least
    // check that the code compiles and that we can copy and assign these arrays
    {
      axom::Array<HasDefault> arr(axom::ArrayOptions::Uninitialized {}, capacity);

      EXPECT_LE(capacity, arr.capacity());
      EXPECT_EQ(capacity, arr.size());

      axom::Array<HasDefault> copied(arr);
      EXPECT_LE(arr.capacity(), copied.capacity());
      EXPECT_EQ(arr.size(), copied.size());
      EXPECT_EQ(arr, copied);

      axom::Array<HasDefault> assigned;
      assigned = arr;
      EXPECT_LE(arr.capacity(), assigned.capacity());
      EXPECT_EQ(arr.size(), assigned.size());
      EXPECT_EQ(arr, assigned);
    }

    // Tests uninitialized with 2D Array
    {
      axom::Array<FailsOnConstruction, 2> arr(axom::ArrayOptions::Uninitialized {}, capacity, capacity);

      EXPECT_LE(capacity * capacity, arr.capacity());
      EXPECT_EQ(capacity * capacity, arr.size());

      arr.resize(axom::ArrayOptions::Uninitialized {}, capacity * 2, capacity * 2);

      EXPECT_LE(capacity * capacity * 4, arr.capacity());
      EXPECT_EQ(capacity * capacity * 4, arr.size());
    }

    // Tests uninitialized with 1D Array with user-supplied allocator
    {
      axom::Array<FailsOnConstruction> arr(axom::ArrayOptions::Uninitialized {},
                                           capacity,
                                           capacity,
                                           axom::getDefaultAllocatorID());

      EXPECT_EQ(capacity, arr.capacity());
      EXPECT_EQ(capacity, arr.size());
      EXPECT_EQ(axom::getDefaultAllocatorID(), arr.getAllocatorID());

      arr.resize(axom::ArrayOptions::Uninitialized {}, capacity * 2);

      EXPECT_LE(capacity * 2, arr.capacity());
      EXPECT_EQ(capacity * 2, arr.size());
      EXPECT_EQ(axom::getDefaultAllocatorID(), arr.getAllocatorID());
    }
  }
}

//------------------------------------------------------------------------------
TEST(core_array, checkConstConversion)
{
  constexpr axom::IndexType size = 10;
  axom::Array<int> v_int(size);
  for(int i = 0; i < v_int.size(); i++)
  {
    v_int[i] = i;
  }
  const axom::Array<int>& v_int_cref = v_int;
  axom::ArrayView<const int> v_int_view = v_int_cref;
  EXPECT_EQ(v_int, v_int_view);
  axom::ArrayView<const int> v_int_view_copy = v_int_view;
  EXPECT_EQ(v_int, v_int_view_copy);
  // axom::ArrayView<int> v_int_view_copy_non_const = v_int_view; // Fails to compile as it should
  // v_int_view[1] = 12; // Fails to compile as it should

  // Check begin() const and end() const
  int idx = 0;
  for(const auto& ele : v_int_cref)
  {
    EXPECT_EQ(ele, idx++);
  }

  axom::Array<double, 2> v_double_2d(size, size);
  for(int i = 0; i < size; i++)
  {
    for(int j = 0; j < size; j++)
    {
      v_double_2d(i, j) = 7.1 * i + j * 2.3;
    }
  }
  const axom::Array<double, 2>& v_double_2d_cref = v_double_2d;
  axom::ArrayView<const double, 2> v_double_2d_view = v_double_2d_cref;
  EXPECT_EQ(v_double_2d, v_double_2d_view);
  axom::ArrayView<const double, 2> v_double_2d_view_copy = v_double_2d_view;
  EXPECT_EQ(v_double_2d, v_double_2d_view_copy);
  // v_double_2d_view_copy(0,0) = 1.1; // Fails to compile as it should
}

//------------------------------------------------------------------------------

TEST(core_array, checkVariadicCtors)
{
  int i = 5;
  size_t s = 6;

  axom::Array<int, 2> arr1(i, i);
  axom::Array<int, 2> arr2(i, s);
  axom::Array<int, 2> arr3(s, i);
  axom::Array<int, 2> arr4(s, s);

  axom::Array<int, 3> arr5(i, i, i);
  axom::Array<int, 3> arr6(i, i, s);
  axom::Array<int, 3> arr7(i, s, i);
  axom::Array<int, 3> arr8(i, s, s);
  axom::Array<int, 3> arr9(s, i, i);
  axom::Array<int, 3> arr10(s, i, s);
  axom::Array<int, 3> arr11(s, s, i);
  axom::Array<int, 3> arr12(s, s, s);
}

//------------------------------------------------------------------------------

TEST(core_array, check_subspan_range)
{
  int m = 10;
  int n = 3;
  int l = 5;
  axom::Array<int> arr(m);

  axom::ArrayView<int> arrv1(arr);
  EXPECT_EQ(arrv1.size(), arr.size());

  axom::ArrayView<int> arrv2 = arrv1.subspan(n);
  EXPECT_EQ(arrv2.size() + n, arrv1.size());
  EXPECT_EQ(&arrv2[0], &arr[n]);
  EXPECT_EQ(&arrv2[arrv2.size() - 1], &arr[arr.size() - 1]);

  axom::ArrayView<int> arrv3 = arrv1.subspan(n, l);
  EXPECT_EQ(arrv3.size(), l);
  EXPECT_EQ(&arrv3[0], &arr[n]);
  EXPECT_EQ(&arrv3[arrv3.size() - 1], &arr[n + l - 1]);
}

//------------------------------------------------------------------------------

template <typename DataType>
void test_resize_with_stackarray(DataType value)
{
  const int I_DIMS = 3;
  const int J_DIMS = 5;
  const int K_DIMS = 7;
  axom::Array<DataType, 2> arr2;

  axom::StackArray<axom::IndexType, 2> dims2 = {{I_DIMS, J_DIMS}};
  arr2.resize(dims2, value);
  EXPECT_EQ(arr2.size(), I_DIMS * J_DIMS);
  EXPECT_EQ(arr2.shape()[0], I_DIMS);
  EXPECT_EQ(arr2.shape()[1], J_DIMS);
  for(int i = 0; i < I_DIMS; i++)
  {
    for(int j = 0; j < J_DIMS; j++)
    {
      EXPECT_EQ(arr2[i][j], value);
    }
  }

  axom::Array<DataType, 3> arr3;

  axom::StackArray<axom::IndexType, 3> dims3 = {{I_DIMS, J_DIMS, K_DIMS}};
  arr3.resize(dims3, value);
  EXPECT_EQ(arr3.size(), I_DIMS * J_DIMS * K_DIMS);
  EXPECT_EQ(arr3.shape()[0], I_DIMS);
  EXPECT_EQ(arr3.shape()[1], J_DIMS);
  EXPECT_EQ(arr3.shape()[2], K_DIMS);
  for(int i = 0; i < I_DIMS; i++)
  {
    for(int j = 0; j < J_DIMS; j++)
    {
      for(int k = 0; k < K_DIMS; k++)
      {
        EXPECT_EQ(arr3[i][j][k], value);
      }
    }
  }
}

TEST(core_array, resize_stackarray)
{
  test_resize_with_stackarray<bool>(false);
  test_resize_with_stackarray<int>(-1);
}

//------------------------------------------------------------------------------
constexpr static int NONTRIVIAL_RELOC_MAGIC = 123;
struct NonTriviallyRelocatable
{
  NonTriviallyRelocatable() : m_member(NONTRIVIAL_RELOC_MAGIC), m_localMemberPointer(&m_member) { }

  NonTriviallyRelocatable(const NonTriviallyRelocatable& other)
    : m_member(other.m_member)
    , m_localMemberPointer(&m_member)
  { }

  NonTriviallyRelocatable& operator=(const NonTriviallyRelocatable& other)
  {
    if(this != &other)
    {
      m_member = other.m_member;
    }
    return *this;
  }

  ~NonTriviallyRelocatable() = default;

  int m_member;
  int* m_localMemberPointer;
};

TEST(core_array, reserve_nontrivial_reloc)
{
  const int NUM_ELEMS = 1024;
  axom::Array<NonTriviallyRelocatable> array(NUM_ELEMS, NUM_ELEMS);

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    // Check initial values.
    EXPECT_EQ(array[i].m_member, NONTRIVIAL_RELOC_MAGIC);
    EXPECT_EQ(&(array[i].m_member), array[i].m_localMemberPointer);
  }

  // Reallocation should work for non-trivially relocatable types.
  array.reserve(NUM_ELEMS * 4);

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    EXPECT_EQ(array[i].m_member, NONTRIVIAL_RELOC_MAGIC);
    EXPECT_EQ(&(array[i].m_member), array[i].m_localMemberPointer);
  }
}

TEST(core_array, reserve_nontrivial_reloc_2)
{
  const int NUM_ELEMS = 1024;
  axom::Array<NonTriviallyRelocatable> array(NUM_ELEMS, NUM_ELEMS);

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    // Check initial values.
    EXPECT_EQ(array[i].m_member, NONTRIVIAL_RELOC_MAGIC);
    EXPECT_EQ(&(array[i].m_member), array[i].m_localMemberPointer);
  }

  EXPECT_EQ(array.capacity(), NUM_ELEMS);

  // Emplace to trigger a resize.
  array.emplace_back(NonTriviallyRelocatable {});
  EXPECT_GE(array.capacity(), NUM_ELEMS);

  for(int i = 0; i < NUM_ELEMS + 1; i++)
  {
    EXPECT_EQ(array[i].m_member, NONTRIVIAL_RELOC_MAGIC);
    EXPECT_EQ(&(array[i].m_member), array[i].m_localMemberPointer);
  }
}

#if defined(AXOM_USE_GPU) && defined(AXOM_USE_UMPIRE)
TEST(core_array, reserve_nontrivial_reloc_um)
{
  const int NUM_ELEMS = 1024;
  axom::Array<NonTriviallyRelocatable, 1, axom::MemorySpace::Unified> array(NUM_ELEMS, NUM_ELEMS);

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    // Check initial values.
    EXPECT_EQ(array[i].m_member, NONTRIVIAL_RELOC_MAGIC);
    EXPECT_EQ(&(array[i].m_member), array[i].m_localMemberPointer);
  }

  // Reallocation should work for non-trivially relocatable types.
  array.reserve(NUM_ELEMS * 4);

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    EXPECT_EQ(array[i].m_member, NONTRIVIAL_RELOC_MAGIC);
    EXPECT_EQ(&(array[i].m_member), array[i].m_localMemberPointer);
  }
}
#endif

class AllocatingDefaultInit
{
public:
  AllocatingDefaultInit() { m_value.resize(32); }

private:
  axom::Array<int> m_value;
};

TEST(core_array, array_ctors_leak)
{
  using Array = axom::Array<AllocatingDefaultInit>;

  // Initialize test array
  constexpr int SZ1 = 20;
  Array data(SZ1);

  // axom::Array copy constructor
  {
    Array arr_copy {data};

    EXPECT_EQ(SZ1, arr_copy.size());
  }

  // axom::Array copy assignment operator
  {
    Array arr_copy {SZ1 + 10};
    arr_copy = data;

    EXPECT_EQ(SZ1, arr_copy.size());
  }

  // axom::Array move constructor
  {
    Array tmp_move_from {data};
    Array arr_move {std::move(tmp_move_from)};

    EXPECT_EQ(SZ1, arr_move.size());
  }

  // axom::Array move assignment operator
  {
    Array tmp_move_from {data};
    Array arr_move {SZ1 + 10};
    arr_move = std::move(tmp_move_from);

    EXPECT_EQ(SZ1, arr_move.size());
  }
}
