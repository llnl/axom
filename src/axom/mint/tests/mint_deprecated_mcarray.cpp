// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/mint/deprecated/MCArray.hpp"

#include "gtest/gtest.h"

// C/C++ includes
#include <algorithm>

namespace axom
{
const char IGNORE_OUTPUT[] = ".*";

namespace internal
{
/*!
 * \brief Calculate the new capacity for an MCArray given an increase in the
 *  size.
 * \param [in] v, the MCArray in question.
 * \param [in] increase, the amount the size will increase by
 * \return the new capacity.
 */
template <typename T>
IndexType calc_new_capacity(deprecated::MCArray<T>& v, IndexType increase)
{
  IndexType new_num_tuples = v.size() + increase;
  if(new_num_tuples > v.capacity())
  {
    return static_cast<IndexType>(new_num_tuples * v.getResizeRatio() + 0.5);
  }

  return v.capacity();
}

/*!
 * \brief Check if two MCArrays are equivalent. Does not check the resize ratio.
 * \param [in] lhs, the first MCArray to compare.
 * \param [in] rhs, the second MCArray to compare.
 * \return the new capacity.
 */
template <typename T>
void check_copy(const deprecated::MCArray<T>& lhs, const deprecated::MCArray<T>& rhs)
{
  EXPECT_EQ(lhs.size(), rhs.size());
  EXPECT_EQ(lhs.numComponents(), rhs.numComponents());
  EXPECT_EQ(lhs.capacity(), rhs.capacity());

  const T* lhs_data = lhs.getData();
  const T* rhs_data = rhs.getData();
  EXPECT_EQ(lhs_data, rhs_data);
}

/*!
 * \brief Check that the storage of an MCArray is working properly.
 * \param [in] v the MCArray to check.
 */
template <typename T>
void check_storage(deprecated::MCArray<T>& v)
{
  EXPECT_TRUE(v.empty());
  EXPECT_EQ(v.size(), 0);

  IndexType capacity = v.capacity();
  IndexType num_components = v.numComponents();
  const T* data_ptr = v.getData();

  /* Append up to half the capacity. */
  if(num_components == 1)
  {
    for(T i = 0; i < capacity / 2; ++i)
    {
      v.append(i);
    }
  }
  else
  {
    T* tuple = new T[num_components];
    for(IndexType i = 0; i < capacity / 2; ++i)
    {
      for(IndexType j = 0; j < num_components; ++j)
      {
        tuple[j] = i * num_components + j;
      }
      v.append(tuple, 1);
    }
    delete[] tuple;
    tuple = nullptr;
  }

  /* Check the MCArray metadata. */
  EXPECT_TRUE(!v.empty());
  EXPECT_EQ(v.size(), capacity / 2);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.getData(), data_ptr);

  /* Append up to the full capacity. */
  if(num_components == 1)
  {
    for(T i = capacity / 2; i < capacity; ++i)
    {
      v.append(i);
    }
  }
  else
  {
    T* tuple = allocate<T>(num_components);
    for(IndexType i = capacity / 2; i < capacity; ++i)
    {
      for(IndexType j = 0; j < num_components; ++j)
      {
        tuple[j] = i * num_components + j;
      }
      v.append(tuple, 1);
    }
    deallocate(tuple);
    tuple = nullptr;
  }

  /* Check the MCArray metadata. */
  EXPECT_TRUE(!v.empty());
  EXPECT_EQ(v.size(), capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.getData(), data_ptr);

  /* Check the MCArray data using the () and [] operators and the raw pointer. */
  for(IndexType i = 0; i < capacity; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * num_components + j);
      EXPECT_EQ(v[i * num_components + j], i * num_components + j);
      EXPECT_EQ(data_ptr[i * num_components + j], i * num_components + j);
    }
  }

  /* Set the MCArray data to new values using the () operator. */
  for(IndexType i = 0; i < capacity; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      v(i, j) = i * j - 5 * i + 7 * j;
    }
  }

  /* Check the MCArray data using the () and [] operators and the raw pointer. */
  for(IndexType i = 0; i < capacity; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
      EXPECT_EQ(v[i * num_components + j], i * j - 5 * i + 7 * j);
      EXPECT_EQ(data_ptr[i * num_components + j], i * j - 5 * i + 7 * j);
    }
  }

  /* Check the MCArray metadata. */
  EXPECT_TRUE(!v.empty());
  EXPECT_EQ(v.size(), capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.getData(), data_ptr);
}

/*!
 * \brief Check that the fill method is working properly.
 * \param [in] v the MCArray to check.
 */
template <typename T>
void check_fill(deprecated::MCArray<T>& v)
{
  constexpr T MAGIC_NUM_0 = 55;
  constexpr T MAGIC_NUM_1 = 6834;
  const IndexType capacity = v.capacity();
  const IndexType size = v.size();
  const IndexType num_components = v.numComponents();
  const double ratio = v.getResizeRatio();
  const T* const data_ptr = v.getData();

  /* Fill the MCArray with MAGIC_NUM_0. */
  v.fill(MAGIC_NUM_0);

  /* Check the meta data. */
  EXPECT_EQ(capacity, v.capacity());
  EXPECT_EQ(size, v.size());
  EXPECT_EQ(num_components, v.numComponents());
  EXPECT_EQ(ratio, v.getResizeRatio());
  EXPECT_EQ(data_ptr, v.getData());

  /* Check that the entries are all MAGIC_NUM_0. */
  for(IndexType i = 0; i < size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM_0);
    }
  }

  /* Fill the MCArray with MAGIC_NUM_1. */
  v.fill(MAGIC_NUM_1);

  /* Check the meta data. */
  EXPECT_EQ(capacity, v.capacity());
  EXPECT_EQ(size, v.size());
  EXPECT_EQ(num_components, v.numComponents());
  EXPECT_EQ(ratio, v.getResizeRatio());
  EXPECT_EQ(data_ptr, v.getData());

  /* Check that the entries are all MAGIC_NUM_1. */
  for(IndexType i = 0; i < size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM_1);
    }
  }
}

/*!
 * \brief Check that the set method is working properly.
 * \param [in] v the MCArray to check.
 */
template <typename T>
void check_set(deprecated::MCArray<T>& v)
{
  constexpr T ZERO_VAL = 0;
  const IndexType capacity = v.capacity();
  const IndexType size = v.size();
  const IndexType num_components = v.numComponents();
  const double ratio = v.getResizeRatio();
  const T* const data_ptr = v.getData();

  /* Allocate a buffer half the size of the MCArray. Fill it up with sequential
   * values. */
  const IndexType buffer_size = size / 2;
  T* buffer = allocate<T>(buffer_size * num_components);
  for(IndexType i = 0; i < buffer_size * num_components; ++i)
  {
    buffer[i] = i;
  }

  /* Set all the values in the MCArray to zero. */
  v.fill(ZERO_VAL);

  /* Set the first half of the tuples in the MCArray to the sequential values in
   * buffer. */
  v.set(buffer, buffer_size, 0);

  /* Check the MCArray metadata. */
  EXPECT_EQ(capacity, v.capacity());
  EXPECT_EQ(size, v.size());
  EXPECT_EQ(num_components, v.numComponents());
  EXPECT_EQ(ratio, v.getResizeRatio());
  EXPECT_EQ(data_ptr, v.getData());

  /* Check that the first half of the tuples in the MCArray are equivalent to
   * those in buffer. */
  for(IndexType i = 0; i < buffer_size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), buffer[i * num_components + j]);
    }
  }

  /* Check that the second half of the tuples in the MCArray are all zero. */
  for(IndexType i = buffer_size; i < size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), ZERO_VAL);
    }
  }

  /* Reset the values in buffer to the next sequential values. */
  for(IndexType i = 0; i < buffer_size * num_components; ++i)
  {
    buffer[i] = i + buffer_size * num_components;
  }

  /* Set the second half of the tuples in the MCArray to the new sequential values
   * in buffer. */
  v.set(buffer, buffer_size, buffer_size);

  /* Check the MCArray metadata. */
  EXPECT_EQ(capacity, v.capacity());
  EXPECT_EQ(size, v.size());
  EXPECT_EQ(num_components, v.numComponents());
  EXPECT_EQ(ratio, v.getResizeRatio());
  EXPECT_EQ(data_ptr, v.getData());

  /* Check that all the tuples in the MCArray now hold sequential values. */
  for(IndexType i = 0; i < 2 * buffer_size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * num_components + j);
    }
  }

  deallocate(buffer);
}

/*!
 * \brief Check that the resizing of an MCArray is working properly.
 * \param [in] v the MCArray to check.
 */
template <typename T>
void check_resize(deprecated::MCArray<T>& v)
{
  /* Resize the MCArray up to the capacity */
  IndexType capacity = v.capacity();
  v.resize(capacity);
  IndexType size = capacity;
  IndexType num_components = v.numComponents();

  /* Check that the size equals the capacity. */
  EXPECT_EQ(v.size(), v.capacity());

  /* Set the existing data in v */
  for(IndexType i = 0; i < size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      v(i, j) = static_cast<T>(i * j - 5 * i + 7 * j);
    }
  }

  /* Append a new tuple, should resize. */
  IndexType old_capacity = capacity;
  capacity = calc_new_capacity(v, 1);
  T* tuple = allocate<T>(num_components);
  for(IndexType j = 0; j < num_components; ++j)
  {
    tuple[j] = size * j - 5 * size + 7 * j;
  }
  v.append(tuple, 1);
  size++;

  /* Check that it resized properly */
  EXPECT_GT(capacity, old_capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check that the data is still intact. */
  for(IndexType i = 0; i < size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
    }
  }

  /* Prepare 1000 tuples to be appended. */
  const IndexType n_tuples = 1000;
  T* values = allocate<T>(n_tuples * num_components);
  for(IndexType i = 0; i < n_tuples; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      IndexType i_real = i + size;
      values[i * num_components + j] = i_real * j - 5 * i_real + 7 * j;
    }
  }

  /* Append the new tuples. */
  capacity = calc_new_capacity(v, n_tuples);
  v.append(values, n_tuples);
  size += n_tuples;

  /* Check that size and capacity are as expected. */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check that the data is still intact. */
  for(IndexType i = 0; i < size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
    }
  }

  /* Reduce the size down to 500 tuples */
  T* data_address = v.getData();
  size = 500;
  v.resize(size);

  /* Check the metadata. */
  EXPECT_EQ(v.size(), size);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.getData(), data_address);

  /* Check the data. */
  for(IndexType i = 0; i < size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
    }
  }

  /* Shrink the vector */
  capacity = size;
  v.shrink();

  /* Check that the capacity and size are as expected. */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check that the data is intact. */
  for(IndexType i = 0; i < size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
    }
  }

  /* Append a new tuple, should resize. */
  old_capacity = capacity;
  capacity = calc_new_capacity(v, 1);
  for(IndexType j = 0; j < num_components; ++j)
  {
    tuple[j] = size * j - 5 * size + 7 * j;
  }
  v.append(tuple, 1);
  size++;

  /* Check the new size and capacity. */
  EXPECT_GT(capacity, old_capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check that the data is intact. */
  for(IndexType i = 0; i < size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
    }
  }

  /* Reset the data */
  T* data_ptr = v.getData();
  for(IndexType i = 0; i < size * num_components; ++i)
  {
    data_ptr[i] = i;
  }

  /* Append a bunch of tuples to fill in up to the capacity. Resize should
   * not occur. */
  old_capacity = capacity;
  for(IndexType i = size; i < old_capacity; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      tuple[j] = i * num_components + j;
    }

    v.append(tuple, 1);
    size++;
    EXPECT_EQ(v.capacity(), old_capacity);
    EXPECT_EQ(v.size(), size);
    EXPECT_EQ(v.getData(), data_ptr);
  }

  EXPECT_EQ(v.size(), old_capacity);

  /* Append a final tuple that should trigger a resize. */
  for(IndexType j = 0; j < num_components; ++j)
  {
    tuple[j] = size * num_components + j;
  }

  capacity = calc_new_capacity(v, old_capacity - size + 1);
  v.append(tuple, 1);
  size++;

  /* Check the new capacity and size. */
  EXPECT_GT(capacity, old_capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check the data. */
  for(IndexType i = 0; i < size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * num_components + j);
    }
  }

  deallocate(tuple);
  tuple = nullptr;

  deallocate(values);
  values = nullptr;
}

/*!
 * \brief Check that the insertion into an MCArray is working properly.
 * \param [in] v the MCArray to check.
 */
template <typename T>
void check_insert(deprecated::MCArray<T>& v)
{
  /* Resize the MCArray up to the capacity */
  IndexType capacity = v.capacity();
  v.resize(capacity);
  IndexType size = capacity;
  IndexType num_components = v.numComponents();

  EXPECT_EQ(v.size(), v.capacity());

  /* Set the existing data in v */
  for(IndexType i = 0; i < size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      v(i, j) = i * j - 5 * i + 7 * j;
    }
  }

  /* Append a new tuple, should resize. */
  IndexType old_capacity = capacity;
  capacity = calc_new_capacity(v, 1);
  T* tuple = allocate<T>(num_components);
  for(IndexType j = 0; j < num_components; ++j)
  {
    tuple[j] = size * j - 5 * size + 7 * j;
  }
  v.insert(tuple, 1, v.size());
  size++;

  /* Check that it resized properly */
  EXPECT_GT(capacity, old_capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);
  for(IndexType i = 0; i < size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
    }
  }

  /* Append 1000 tuples */
  const IndexType n_tuples = 1000;
  T* values = allocate<T>(n_tuples * num_components);
  for(IndexType i = 0; i < n_tuples; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      IndexType i_real = i + size;
      values[i * num_components + j] = i_real * j - 5 * i_real + 7 * j;
    }
  }

  capacity = calc_new_capacity(v, n_tuples);
  v.insert(values, n_tuples, size);
  size += n_tuples;

  /* Check that it resizes properly */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);
  for(IndexType i = 0; i < size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
    }
  }

  capacity = size;
  v.shrink();
  IndexType n_insert_front = 100;

  /* Reset the data */
  T* data_ptr = v.getData();
  for(IndexType i = 0; i < size * num_components; ++i)
  {
    data_ptr[i] = i + num_components * n_insert_front;
  }

  /* Insert into the front of the MCArray. */
  for(IndexType i = n_insert_front - 1; i >= 0; i--)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      tuple[j] = i * num_components + j;
    }
    capacity = calc_new_capacity(v, 1);
    v.insert(tuple, 1, 0);
    size++;
  }

  /* Check that the insertion worked as expected */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);
  for(IndexType i = 0; i < size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * num_components + j);
    }
  }

  deallocate(tuple);
  tuple = nullptr;

  deallocate(values);
  values = nullptr;
}

/*!
 * \brief Check that the emplace method is working properly.
 * \param [in] v the MCArray to check.
 */
template <typename T>
void check_emplace(deprecated::MCArray<T>& v)
{
  constexpr T MAGIC_NUM = 52706;

  /* Resize the MCArray up to the capacity */
  IndexType capacity = v.capacity();
  v.resize(capacity);
  IndexType size = capacity;
  const IndexType num_components = v.numComponents();

  EXPECT_EQ(v.size(), v.capacity());

  /* Set the existing data in v */
  for(IndexType i = 0; i < size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      v(i, j) = i * num_components + j;
    }
  }

  /* Emplace 1 tuple at the end of the MCArray with the default value,
   * should resize. */
  capacity = calc_new_capacity(v, 1);
  v.emplace(1, v.size());

  /* Emplace 9 tuples at the end of the MCArray with MAGIC_NUM. */
  capacity = calc_new_capacity(v, 9);
  v.emplace(9, v.size(), MAGIC_NUM);
  size += 10;

  /* Check that it resized properly */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check the data */
  for(IndexType i = 0; i < size - 10; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * num_components + j);
    }
  }

  for(IndexType j = 0; j < num_components; ++j)
  {
    EXPECT_EQ(v(size - 10, j), T());
  }

  for(IndexType i = size - 9; i < size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM);
    }
  }

  /* Emplace 10 tuples at the beginning of the MCArray with MAGIC_NUM. */
  capacity = calc_new_capacity(v, 9);
  v.emplace(9, 0, MAGIC_NUM);

  /* Emplace 1 tuple at the beginning of the MCArray with the default value. */
  capacity = calc_new_capacity(v, 1);
  v.emplace(1, 0);
  size += 10;

  /* Check that it resized properly */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check the beginning */
  for(IndexType j = 0; j < num_components; ++j)
  {
    EXPECT_EQ(v(0, j), T());
  }

  for(IndexType i = 1; i < 10; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM);
    }
  }

  /* Check the middle */
  for(IndexType i = 10; i < size - 10; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), (i - 10) * num_components + j);
    }
  }

  /* Check the end */
  for(IndexType j = 0; j < num_components; ++j)
  {
    EXPECT_EQ(v(size - 10, j), T());
  }

  for(IndexType i = size - 9; i < size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM);
    }
  }

  /* Emplace 9 tuples in the middle of the MCArray with MAGIC_NUM. */
  IndexType middle = size / 2;
  capacity = calc_new_capacity(v, 9);
  v.emplace(9, middle, MAGIC_NUM);

  /* Emplace 1 tuple in the middle of the MCArray with the default value. */
  capacity = calc_new_capacity(v, 1);
  v.emplace(1, middle);
  size += 10;

  /* Check that it resized properly */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check the beginning */
  for(IndexType j = 0; j < num_components; ++j)
  {
    EXPECT_EQ(v(0, j), T());
  }

  for(IndexType i = 1; i < 10; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM);
    }
  }

  /* Check the first section */
  for(IndexType i = 10; i < middle; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), (i - 10) * num_components + j);
    }
  }

  /* Check the middle */
  for(IndexType j = 0; j < num_components; ++j)
  {
    EXPECT_EQ(v(middle, j), T());
  }

  for(IndexType i = middle + 1; i < middle + 10; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM);
    }
  }

  /* Check the second section */
  for(IndexType i = middle + 10; i < size - 10; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), (i - 20) * num_components + j);
    }
  }

  /* Check the end */
  for(IndexType j = 0; j < num_components; ++j)
  {
    EXPECT_EQ(v(size - 10, j), T());
  }

  for(IndexType i = size - 9; i < size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM);
    }
  }
}

/*!
 * \brief Check an external MCArray for defects.
 * \param [in] v the external MCArray to check.
 */
template <typename T>
void check_external(deprecated::MCArray<T>& v)
{
  ASSERT_TRUE(v.isExternal());

  /* Check that the MCArray is full. */
  ASSERT_EQ(v.size(), v.capacity());

  const IndexType size = v.size();
  const IndexType num_components = v.numComponents();
  const IndexType num_values = size * num_components;
  T* const data_ptr = v.getData();

  /* Set the tuples in the MCArray. */
  for(IndexType i = 0; i < size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      v(i, j) = i * num_components + j;
    }
  }

  /* Check the tuples using the raw pointer. */
  for(IndexType i = 0; i < num_values; ++i)
  {
    EXPECT_EQ(data_ptr[i], i);
  }

  /* Set the tuples using the raw pointer. */
  for(IndexType i = 0; i < size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      data_ptr[i * num_components + j] = i * j - i - j;
    }
  }

  /* Check the tuples using the () operator. */
  for(IndexType i = 0; i < size; ++i)
  {
    for(IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - i - j);
    }
  }

  EXPECT_EQ(size, v.size());
  EXPECT_EQ(size, v.capacity());
  EXPECT_EQ(data_ptr, v.getData());

  /* Since the MCArray is full all of the following calls should require a
   * reallocation and cause a fatal error. */
  T* tuple = allocate<T>(num_components);
  EXPECT_DEATH_IF_SUPPORTED(v.append(tuple, 1), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(v.insert(tuple, 1, 0), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(v.reserve(size + 1), IGNORE_OUTPUT);

  deallocate(tuple);
  tuple = nullptr;
}

} /* end namespace internal */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
TEST(core_mcarray, checkStorage)
{
  constexpr IndexType ZERO_VAL = 0;

  for(IndexType capacity = 2; capacity < 512; capacity *= 2)
  {
    for(IndexType n_components = 1; n_components <= 4; n_components++)
    {
      deprecated::MCArray<int> v_int(ZERO_VAL, n_components, capacity);
      internal::check_storage(v_int);

      deprecated::MCArray<double> v_double(ZERO_VAL, n_components, capacity);
      internal::check_storage(v_double);
    }
  }
}

//------------------------------------------------------------------------------
TEST(core_mcarray, checkFill)
{
  for(IndexType capacity = 2; capacity < 512; capacity *= 2)
  {
    IndexType size = capacity / 2;
    for(IndexType n_components = 1; n_components <= 4; n_components++)
    {
      deprecated::MCArray<int> v_int(size, n_components, capacity);
      internal::check_fill(v_int);

      deprecated::MCArray<double> v_double(size, n_components, capacity);
      internal::check_fill(v_double);
    }
  }
}

//------------------------------------------------------------------------------
TEST(core_mcarray, checkSet)
{
  for(IndexType size = 2; size < 512; size *= 2)
  {
    for(IndexType n_components = 1; n_components <= 4; n_components++)
    {
      deprecated::MCArray<int> v_int(size, n_components);
      internal::check_set(v_int);

      deprecated::MCArray<double> v_double(size, n_components);
      internal::check_set(v_double);
    }
  }
}

//------------------------------------------------------------------------------
TEST(core_mcarray, checkResize)
{
  constexpr IndexType ZERO_VAL = 0;

  for(double ratio = 1.0; ratio <= 2.0; ratio += 0.5)
  {
    for(IndexType capacity = 2; capacity <= 512; capacity *= 2)
    {
      for(IndexType n_components = 1; n_components <= 4; n_components++)
      {
        deprecated::MCArray<int> v_int(ZERO_VAL, n_components, capacity);
        v_int.setResizeRatio(ratio);
        internal::check_resize(v_int);

        deprecated::MCArray<double> v_double(ZERO_VAL, n_components, capacity);
        v_double.setResizeRatio(ratio);
        internal::check_resize(v_double);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(core_mcarray_DeathTest, deathtest_checkResize)
{
  /* Resizing isn't allowed with a ratio less than 1.0. */
  constexpr IndexType ZERO = 0;
  deprecated::MCArray<int> v_int(ZERO, 1, 100);
  v_int.setResizeRatio(0.99);
  EXPECT_DEATH_IF_SUPPORTED(internal::check_resize(v_int), IGNORE_OUTPUT);
}

//------------------------------------------------------------------------------
TEST(core_mcarray, checkInsert)
{
  constexpr IndexType ZERO_VAL = 0;

  for(double ratio = 1.0; ratio <= 2.0; ratio += 0.5)
  {
    for(IndexType capacity = 2; capacity <= 512; capacity *= 2)
    {
      for(IndexType n_components = 1; n_components <= 3; n_components++)
      {
        deprecated::MCArray<int> v_int(ZERO_VAL, n_components, capacity);
        v_int.setResizeRatio(ratio);
        internal::check_insert(v_int);

        deprecated::MCArray<double> v_double(ZERO_VAL, n_components, capacity);
        v_double.setResizeRatio(ratio);
        internal::check_insert(v_double);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(core_mcarray, checkEmplace)
{
  constexpr IndexType ZERO_VAL = 0;

  for(double ratio = 1.0; ratio <= 2.0; ratio += 0.5)
  {
    for(IndexType capacity = 2; capacity <= 512; capacity *= 2)
    {
      for(IndexType n_components = 1; n_components <= 3; n_components++)
      {
        deprecated::MCArray<int> v_int(ZERO_VAL, n_components, capacity);
        v_int.setResizeRatio(ratio);
        internal::check_emplace(v_int);

        deprecated::MCArray<double> v_double(ZERO_VAL, n_components, capacity);
        v_double.setResizeRatio(ratio);
        internal::check_emplace(v_double);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(core_mcarray_DeathTest, checkExternal)
{
  constexpr double MAGIC_NUM = 5683578.8;
  constexpr IndexType MAX_SIZE = 256;
  constexpr IndexType MAX_COMPONENTS = 3;
  constexpr IndexType MAX_VALUES = MAX_SIZE * MAX_COMPONENTS;
  union DataBuffer
  {
    int ints[MAX_SIZE * MAX_COMPONENTS];
    double doubles[MAX_SIZE * MAX_COMPONENTS];
  };

  DataBuffer buffer;
  std::fill_n(buffer.doubles, MAX_VALUES, MAGIC_NUM);

  for(IndexType size = 16; size <= MAX_SIZE; size *= 2)
  {
    for(IndexType n_comp = 1; n_comp <= MAX_COMPONENTS; n_comp++)
    {
      deprecated::MCArray<int> v_int(buffer.ints, size, n_comp);
      EXPECT_EQ(v_int.getData(), buffer.ints);
      internal::check_external(v_int);

      deprecated::MCArray<double> v_double(buffer.doubles, size, n_comp);
      EXPECT_EQ(v_double.getData(), buffer.doubles);
      internal::check_external(v_double);

      /* Set v_double's data to MAGIC_NUM */
      v_double.fill(MAGIC_NUM);
    }

    /* Check that the data still exists in the buffer */
    for(IndexType i = 0; i < MAX_VALUES; ++i)
    {
      EXPECT_EQ(buffer.doubles[i], MAGIC_NUM);
    }
  }
}

} /* end namespace axom */

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
