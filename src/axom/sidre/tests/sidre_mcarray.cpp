// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/sidre.hpp"

#include "gtest/gtest.h"

// C/C++ includes
#include <algorithm> /* for std::fill_n */

namespace axom
{
namespace sidre
{
const char IGNORE_OUTPUT[] = ".*";

namespace internal
{
/*!
 * \brief Calculate the new capacity for an MCArray given an increase in the
 *  size.
 * \param [in] v, the MCArray in question.
 * \param [in] increase, the number of tuples the size will increase by
 * \return the new capacity, in number of tuples
 */
template <typename T>
axom::IndexType calc_new_capacity(MCArray<T>& v, axom::IndexType increase)
{
  const auto num_components = v.shape()[1];
  axom::IndexType new_num_tuples = (v.size() / num_components) + increase;
  if((new_num_tuples * num_components) > v.capacity())
  {
    axom::IndexType new_capacity = v.capacity() * v.getResizeRatio() + 0.5;
    axom::IndexType remainder = new_capacity % num_components;
    if(remainder > 0)
    {
      new_capacity += num_components - remainder;
    }
    return axom::utilities::max<axom::IndexType>(new_capacity / num_components, new_num_tuples);
  }

  return v.capacity() / num_components;
}

/*!
 * \brief Check if two MCArrays are equivalent. Does not check the resize ratio.
 * \param [in] lhs, the first MCArray to compare.
 * \param [in] rhs, the second MCArray to compare.
 * \return the new capacity.
 */
template <typename T>
void check_copy(const MCArray<T>& lhs, const MCArray<T>& rhs)
{
  EXPECT_EQ(lhs.size(), rhs.size());
  EXPECT_EQ(lhs.shape()[1], rhs.shape()[1]);
  EXPECT_EQ(lhs.capacity(), rhs.capacity());

  for(int i = 0; i < lhs.shape()[0]; i++)
  {
    for(int j = 0; j < lhs.shape()[1]; j++)
    {
      EXPECT_EQ(lhs(i, j), rhs(i, j));
    }
  }
}

/*!
 * \brief Check if two MCArrays are referentially equivalent. Does not check the resize ratio.
 * \param [in] lhs, the first MCArray to compare.
 * \param [in] rhs, the second MCArray to compare.
 * \return the new capacity.
 */
template <typename T>
void check_referential_copy(const MCArray<T>& lhs, const MCArray<T>& rhs)
{
  EXPECT_EQ(lhs.size(), rhs.size());
  EXPECT_EQ(lhs.shape()[1], rhs.shape()[1]);
  EXPECT_EQ(lhs.capacity(), rhs.capacity());

  const T* lhs_data = lhs.data();
  const T* rhs_data = rhs.data();
  EXPECT_EQ(lhs_data, rhs_data);
}

/*!
 * \brief Check that the storage of an MCArray is working properly.
 * \param [in] v the MCArray to check.
 */
template <typename T>
void check_storage(MCArray<T>& v)
{
  EXPECT_TRUE(v.empty());
  EXPECT_EQ(v.size(), 0);

  axom::IndexType num_components = v.shape()[1];
  axom::IndexType num_tuples_capacity = v.capacity() / num_components;
  // FIXME: Messy way of making sure the number of tuples is divisible by two
  if(num_tuples_capacity % 2 != 0)
  {
    v.reserve(v.capacity() + num_components);
    num_tuples_capacity++;
  }
  axom::IndexType capacity = v.capacity();
  const T* data_ptr = v.data();

  /* Append up to half the capacity. */
  {
    T* tuple = new T[num_components];
    for(axom::IndexType i = 0; i < num_tuples_capacity / 2; ++i)
    {
      for(axom::IndexType j = 0; j < num_components; ++j)
      {
        tuple[j] = i * num_components + j;
      }
      v.append(axom::MCArrayView<T>(tuple, 1, num_components));
    }
    delete[] tuple;
    tuple = nullptr;
  }

  /* Check the MCArray metadata. */
  EXPECT_TRUE(!v.empty());
  EXPECT_EQ(v.size(), capacity / 2);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.data(), data_ptr);

  /* Append up to the full capacity. */
  {
    T* tuple = new T[num_components];
    for(axom::IndexType i = num_tuples_capacity / 2; i < num_tuples_capacity; ++i)
    {
      for(axom::IndexType j = 0; j < num_components; ++j)
      {
        tuple[j] = i * num_components + j;
      }
      v.append(axom::MCArrayView<T>(tuple, 1, num_components));
    }
    delete[] tuple;
    tuple = nullptr;
  }

  /* Check the MCArray metadata. */
  EXPECT_TRUE(!v.empty());
  EXPECT_EQ(v.size(), capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.data(), data_ptr);

  /* Check the MCArray data using the () and [] operators and the raw pointer. */
  for(axom::IndexType i = 0; i < num_tuples_capacity; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * num_components + j);
      EXPECT_EQ(v.flatIndex(i * num_components + j), i * num_components + j);
      EXPECT_EQ(data_ptr[i * num_components + j], i * num_components + j);
    }
  }

  /* Set the MCArray data to new values using the () operator. */
  for(axom::IndexType i = 0; i < num_tuples_capacity; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      v(i, j) = i * j - 5 * i + 7 * j;
    }
  }

  /* Check the MCArray data using the () and [] operators and the raw pointer. */
  for(axom::IndexType i = 0; i < num_tuples_capacity; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
      EXPECT_EQ(v.flatIndex(i * num_components + j), i * j - 5 * i + 7 * j);
      EXPECT_EQ(data_ptr[i * num_components + j], i * j - 5 * i + 7 * j);
    }
  }

  /* Check the MCArray metadata. */
  EXPECT_TRUE(!v.empty());
  EXPECT_EQ(v.size(), capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.data(), data_ptr);
}

/*!
 * \brief Check that the fill method is working properly.
 * \param [in] v the MCArray to check.
 */
template <typename T>
void check_fill(MCArray<T>& v)
{
  constexpr T MAGIC_NUM_0 = 55;
  constexpr T MAGIC_NUM_1 = 6834;
  const axom::IndexType capacity = v.capacity();
  const axom::IndexType size = v.size();
  const axom::IndexType num_components = v.shape()[1];
  const double ratio = v.getResizeRatio();
  const T* const data_ptr = v.data();

  /* Fill the MCArray with MAGIC_NUM_0. */
  v.fill(MAGIC_NUM_0);

  /* Check the meta data. */
  EXPECT_EQ(capacity, v.capacity());
  EXPECT_EQ(size, v.size());
  EXPECT_EQ(num_components, v.shape()[1]);
  EXPECT_EQ(ratio, v.getResizeRatio());
  EXPECT_EQ(data_ptr, v.data());

  /* Check that the entries are all MAGIC_NUM_0. */
  // NOTE: size now refers to the number of Ts and not the number of tuples
  for(axom::IndexType i = 0; i < size / num_components; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM_0);
    }
  }

  /* Fill the MCArray with MAGIC_NUM_1. */
  v.fill(MAGIC_NUM_1);

  /* Check the meta data. */
  EXPECT_EQ(capacity, v.capacity());
  EXPECT_EQ(size, v.size());
  EXPECT_EQ(num_components, v.shape()[1]);
  EXPECT_EQ(ratio, v.getResizeRatio());
  EXPECT_EQ(data_ptr, v.data());

  /* Check that the entries are all MAGIC_NUM_1. */
  for(axom::IndexType i = 0; i < size / num_components; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
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
void check_set(MCArray<T>& v)
{
  constexpr T ZERO = 0;
  const axom::IndexType capacity = v.capacity();
  const axom::IndexType num_components = v.shape()[1];
  const axom::IndexType num_tuples = v.size() / num_components;
  const double ratio = v.getResizeRatio();
  const T* const data_ptr = v.data();

  /* Allocate a buffer half the size of the MCArray. Fill it up with sequential
   * values. */
  // buffer_size is the number of tuples
  const axom::IndexType buffer_size = num_tuples / 2;
  T* buffer = new T[buffer_size * num_components];
  for(axom::IndexType i = 0; i < buffer_size * num_components; ++i)
  {
    buffer[i] = i;
  }

  /* Set all the values in the MCArray to zero. */
  v.fill(ZERO);

  /* Set the first half of the tuples in the MCArray to the sequential values in
   * buffer. */
  v.set(buffer, buffer_size * num_components, 0);

  /* Check the MCArray metadata. */
  EXPECT_EQ(capacity, v.capacity());
  EXPECT_EQ(num_tuples, v.size() / num_components);
  EXPECT_EQ(num_components, v.shape()[1]);
  EXPECT_EQ(ratio, v.getResizeRatio());
  EXPECT_EQ(data_ptr, v.data());

  /* Check that the first half of the tuples in the MCArray are equivalent to
   * those in buffer. */
  for(axom::IndexType i = 0; i < buffer_size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), buffer[i * num_components + j]);
    }
  }

  /* Check that the second half of the tuples in the MCArray are all zero. */
  for(axom::IndexType i = buffer_size; i < num_tuples; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), ZERO);
    }
  }

  /* Reset the values in buffer to the next sequential values. */
  for(axom::IndexType i = 0; i < buffer_size * num_components; ++i)
  {
    buffer[i] = i + (buffer_size * num_components);
  }

  /* Set the second half of the tuples in the MCArray to the new sequential values
   * in buffer. */
  v.set(buffer, buffer_size * num_components, buffer_size * num_components);

  /* Check the MCArray metadata. */
  EXPECT_EQ(capacity, v.capacity());
  EXPECT_EQ(num_tuples, v.size() / num_components);
  EXPECT_EQ(num_components, v.shape()[1]);
  EXPECT_EQ(ratio, v.getResizeRatio());
  EXPECT_EQ(data_ptr, v.data());

  /* Check that all the tuples in the MCArray now hold sequential values. */
  for(axom::IndexType i = 0; i < 2 * buffer_size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * num_components + j);
    }
  }

  delete[] buffer;
}

/*!
 * \brief Check that the resizing of an MCArray is working properly.
 * \param [in] v the MCArray to check.
 */
template <typename T>
void check_resize(MCArray<T>& v)
{
  /* Resize the array up to the capacity */
  axom::IndexType capacity = v.capacity();
  axom::IndexType num_components = v.shape()[1];
  v.resize(capacity / num_components, num_components);
  axom::IndexType size = capacity;

  /* Check that the size equals the capacity. */
  EXPECT_EQ(v.size(), v.capacity());

  /* Set the existing data in v */
  for(axom::IndexType i = 0; i < size / num_components; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      v(i, j) = i * j - 5 * i + 7 * j;
    }
  }

  /* Append a new tuple, should resize. */
  axom::IndexType old_capacity = capacity;
  capacity = calc_new_capacity(v, 1) * num_components;
  T* tuple = new T[num_components];
  for(axom::IndexType j = 0; j < num_components; ++j)
  {
    tuple[j] = (size / num_components) * j - 5 * (size / num_components) + 7 * j;
  }
  v.append(axom::MCArrayView<T>(tuple, 1, num_components));
  size += num_components;

  /* Check that it resized properly */
  EXPECT_GT(capacity, old_capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check that the data is still intact. */
  for(axom::IndexType i = 0; i < size / num_components; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
    }
  }

  /* Prepare 1000 tuples to be appended. */
  const axom::IndexType n_tuples = 1000;
  T* values = new T[n_tuples * num_components];
  for(axom::IndexType i = 0; i < n_tuples; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      axom::IndexType i_real = i + (size / num_components);
      values[i * num_components + j] = i_real * j - 5 * i_real + 7 * j;
    }
  }

  /* Append the new tuples. */
  capacity = calc_new_capacity(v, n_tuples) * num_components;
  v.append(axom::MCArrayView<T>(values, n_tuples, num_components));
  size += n_tuples * num_components;

  /* Check that size and capacity are as expected. */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check that the data is still intact. */
  for(axom::IndexType i = 0; i < size / num_components; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
    }
  }

  /* Reduce the size down to 500 tuples */
  T* data_address = v.data();
  size = 500 * num_components;
  v.resize(size / num_components, num_components);

  /* Check the metadata. */
  EXPECT_EQ(v.size(), size);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.data(), data_address);

  /* Check the data. */
  for(axom::IndexType i = 0; i < size / num_components; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
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
  for(axom::IndexType i = 0; i < size / num_components; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
    }
  }

  /* Append a new tuple, should resize. */
  old_capacity = capacity;
  capacity = calc_new_capacity(v, 1) * num_components;
  for(axom::IndexType j = 0; j < num_components; ++j)
  {
    tuple[j] = (size / num_components) * j - 5 * (size / num_components) + 7 * j;
  }
  v.append(axom::MCArrayView<T>(tuple, 1, num_components));
  size += num_components;

  /* Check the new size and capacity. */
  EXPECT_GT(capacity, old_capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check that the data is intact. */
  for(axom::IndexType i = 0; i < size / num_components; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
    }
  }

  /* Reset the data */
  T* data_ptr = v.data();
  for(axom::IndexType i = 0; i < size; ++i)
  {
    data_ptr[i] = i;
  }

  /* Append a bunch of tuples to fill in up to the capacity. Resize should
   * not occur. */
  old_capacity = capacity;
  for(axom::IndexType i = size / num_components; i < old_capacity / num_components; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      tuple[j] = i * num_components + j;
    }

    v.append(axom::MCArrayView<T>(tuple, 1, num_components));
    size += num_components;
    EXPECT_EQ(v.capacity(), old_capacity);
    EXPECT_EQ(v.size(), size);
    EXPECT_EQ(v.data(), data_ptr);
  }

  EXPECT_EQ(v.size(), old_capacity);

  /* Append a final tuple that should trigger a resize. */
  for(axom::IndexType j = 0; j < num_components; ++j)
  {
    tuple[j] = size + j;
  }

  capacity = calc_new_capacity(v, old_capacity - size + 1) * num_components;
  v.append(axom::MCArrayView<T>(tuple, 1, num_components));
  size += num_components;

  /* Check the new capacity and size. */
  EXPECT_GT(capacity, old_capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check the data. */
  for(axom::IndexType i = 0; i < size / num_components; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * num_components + j);
    }
  }

  delete[] tuple;
  tuple = nullptr;

  delete[] values;
  values = nullptr;
}

/*!
 * \brief Check that the insertion into an MCArray is working properly.
 * \param [in] v the MCArray to check.
 */
template <typename T>
void check_insert(MCArray<T>& v)
{
  /* Resize the MCArray up to the capacity */
  axom::IndexType capacity = v.capacity();
  axom::IndexType num_components = v.shape()[1];
  v.resize(capacity / num_components, num_components);
  axom::IndexType size = capacity;

  EXPECT_EQ(v.size(), v.capacity());

  /* Set the existing data in v */
  for(axom::IndexType i = 0; i < size / num_components; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      v(i, j) = i * j - 5 * i + 7 * j;
    }
  }

  /* Append a new tuple, should resize. */
  axom::IndexType old_capacity = capacity;
  capacity = calc_new_capacity(v, 1) * num_components;
  T* tuple = new T[num_components];
  for(axom::IndexType j = 0; j < num_components; ++j)
  {
    tuple[j] = (size / num_components) * j - 5 * (size / num_components) + 7 * j;
  }
  v.append(axom::MCArrayView<T>(tuple, 1, num_components));
  size += num_components;

  /* Check that it resized properly */
  EXPECT_GT(capacity, old_capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);
  for(axom::IndexType i = 0; i < size / num_components; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
    }
  }

  /* Append 1000 tuples */
  const axom::IndexType n_tuples = 1000;
  T* values = new T[n_tuples * num_components];
  for(axom::IndexType i = 0; i < n_tuples; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      axom::IndexType i_real = i + (size / num_components);
      values[i * num_components + j] = i_real * j - 5 * i_real + 7 * j;
    }
  }

  capacity = calc_new_capacity(v, n_tuples) * num_components;
  v.append(axom::MCArrayView<T>(values, n_tuples, num_components));
  size += n_tuples * num_components;

  /* Check that it resizes properly */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);
  for(axom::IndexType i = 0; i < size / num_components; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
    }
  }
  // std::cout << "capcity 1 is " << v.capacity() << "\n";
  capacity = size;
  v.shrink();
  axom::IndexType n_insert_front = 2;
  // std::cout << "capcity 1-1 is " << v.capacity() << "\n";

  /* Reset the data */
  T* data_ptr = v.data();
  for(axom::IndexType i = 0; i < size; ++i)
  {
    data_ptr[i] = i + num_components * n_insert_front;
  }
  // std::cout << "capcity 2 is " << v.capacity() << "\n";

  /* Insert into the front of the MCArray. */
  for(axom::IndexType i = n_insert_front - 1; i >= 0; i--)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      tuple[j] = i * num_components + j;
    }
    capacity = calc_new_capacity(v, 1) * num_components;
    v.insert(0, axom::MCArrayView<T>(tuple, 1, num_components));
    size += num_components;
  }
  // std::cout << "capcity 3 is " << v.capacity() << "\n";

  /* Check that the insertion worked as expected */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);
  for(axom::IndexType i = 0; i < size / num_components; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * num_components + j);
    }
  }

  delete[] tuple;
  tuple = nullptr;

  delete[] values;
  values = nullptr;
}

/*!
 * \brief Check that the emplace method is working properly.
 * \param [in] v the MCArray to check.
 */
template <typename T>
void check_emplace(MCArray<T>& v)
{
  constexpr T MAGIC_NUM = 52706;

  /* Resize the MCArray up to the capacity */
  axom::IndexType capacity = v.capacity();
  axom::IndexType num_components = v.shape()[1];
  v.resize(capacity / num_components, num_components);
  axom::IndexType size = capacity;

  EXPECT_EQ(v.size(), v.capacity());

  /* Set the existing data in v */
  for(axom::IndexType i = 0; i < size / num_components; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      v(i, j) = i * num_components + j;
    }
  }

  /* Emplace 1 tuple at the end of the MCArray with the default value,
   * should resize. */
  capacity = calc_new_capacity(v, 1) * num_components;
  v.emplace(v.size(), T {});

  /* Emplace 9 tuples at the end of the MCArray with MAGIC_NUM. */
  for(axom::IndexType i = 0; i < 9; i++)
  {
    capacity = calc_new_capacity(v, 1) * num_components;
    v.emplace(v.size(), MAGIC_NUM);
  }
  size += 10 * num_components;

  /* Check that it resized properly */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check the data */
  for(axom::IndexType i = 0; i < (size / num_components) - 10; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * num_components + j);
    }
  }

  for(axom::IndexType j = 0; j < num_components; ++j)
  {
    EXPECT_EQ(v((size / num_components) - 10, j), T());
  }

  for(axom::IndexType i = (size / num_components) - 9; i < size / num_components; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM);
    }
  }

  /* Emplace 9 tuples at the beginning of the MCArray with MAGIC_NUM. */
  for(axom::IndexType i = 0; i < 9; i++)
  {
    capacity = calc_new_capacity(v, 1) * num_components;
    v.emplace(0, MAGIC_NUM);
  }

  /* Emplace 1 tuple at the beginning of the MCArray with the default value. */
  capacity = calc_new_capacity(v, 1) * num_components;
  v.emplace(0, T {});
  size += 10 * num_components;

  /* Check that it resized properly */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check the beginning */
  for(axom::IndexType j = 0; j < num_components; ++j)
  {
    EXPECT_EQ(v(0, j), T());
  }

  for(axom::IndexType i = 1; i < 10; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM);
    }
  }

  /* Check the middle */
  for(axom::IndexType i = 10; i < (size / num_components) - 10; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), (i - 10) * num_components + j);
    }
  }

  /* Check the end */
  for(axom::IndexType j = 0; j < num_components; ++j)
  {
    EXPECT_EQ(v(size - 10, j), T());
  }

  for(axom::IndexType i = (size / num_components) - 9; i < size / num_components; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM);
    }
  }

  /* Emplace 9 tuples in the middle of the MCArray with MAGIC_NUM. */
  axom::IndexType middle = (size / num_components) / 2;
  for(axom::IndexType i = 0; i < 9; i++)
  {
    capacity = calc_new_capacity(v, 1) * num_components;
    v.emplace(middle, MAGIC_NUM);
  }

  /* Emplace 1 tuple in the middle of the MCArray with the default value. */
  capacity = calc_new_capacity(v, 1) * num_components;
  v.emplace(middle, T {});
  size += 10 * num_components;

  /* Check that it resized properly */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check the beginning */
  for(axom::IndexType j = 0; j < num_components; ++j)
  {
    EXPECT_EQ(v(0, j), T());
  }

  for(axom::IndexType i = 1; i < 10; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM);
    }
  }

  /* Check the first section */
  for(axom::IndexType i = 10; i < middle; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), (i - 10) * num_components + j);
    }
  }

  /* Check the middle */
  for(axom::IndexType j = 0; j < num_components; ++j)
  {
    EXPECT_EQ(v(middle, j), T());
  }

  for(axom::IndexType i = middle + 1; i < middle + 10; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM);
    }
  }

  /* Check the second section */
  for(axom::IndexType i = middle + 10; i < size - 10; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), (i - 20) * num_components + j);
    }
  }

  /* Check the end */
  for(axom::IndexType j = 0; j < num_components; ++j)
  {
    EXPECT_EQ(v(size - 10, j), T());
  }

  for(axom::IndexType i = size - 9; i < size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM);
    }
  }
}

/*!
 * \brief Check that the move constructors/assignment operators are working properly.
 * \param [in] v the MCArray to check.
 */
template <typename T>
void check_move(MCArray<T>& v)
{
  /* Resize the MCArray up to the capacity */
  axom::IndexType capacity = v.capacity();
  axom::IndexType num_components = v.shape()[1];
  v.resize(capacity / num_components, num_components);
  axom::IndexType size = capacity;

  /* Set the existing data in v */
  for(axom::IndexType i = 0; i < size / num_components; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      v(i, j) = i * j - 5 * i + 7 * j;
    }
  }

  /* Create a copy. */
  MCArray<T> cpy(const_cast<sidre::View*>(v.getView()));

  // Move from the copy (call the move constructor)
  MCArray<T> moved(std::move(cpy));
  EXPECT_EQ(cpy.getView(), nullptr);

  // Check that the moved-to object is the same as the original
  check_copy(v, moved);
  EXPECT_EQ(v.getView(), moved.getView());
}

/*!
 * \brief Make a copy of an MCArray through sidre and check it for defects.
 * \param [in] v the MCArray to copy.
 */
template <typename T>
void check_sidre(MCArray<T>& v)
{
  // These should be equivalent
  ASSERT_NE(dynamic_cast<sidre::MCArray<T>*>(&v), nullptr);
  // ASSERT_TRUE(v.isInSidre());
  /* Create a copy. */
  MCArray<T> cpy(const_cast<sidre::View*>(v.getView()));
  cpy.setResizeRatio(v.getResizeRatio());

  /* Check that the copy holds the same data. */
  check_referential_copy(v, cpy);

  /* Resize the copy and check that it functions correctly. */
  cpy.resize(0, cpy.shape()[1]);
  check_storage(cpy);
  check_insert(cpy);
}

} /* end namespace internal */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
TEST(sidre_core_mcarray, checkStorage)
{
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr axom::IndexType ZERO = 0;

  for(axom::IndexType capacity = 2; capacity < 512; capacity *= 2)
  {
    for(axom::IndexType n_components = 1; n_components <= 4; n_components++)
    {
      MCArray<int> v_int_sidre(root->createView("int"), ZERO, n_components, capacity * n_components);
      internal::check_storage(v_int_sidre);

      MCArray<double> v_double_sidre(root->createView("double"),
                                     ZERO,
                                     n_components,
                                     capacity * n_components);
      internal::check_storage(v_double_sidre);

      root->destroyViewsAndData();
    }
  }
}

//------------------------------------------------------------------------------
TEST(sidre_core_mcarray, checkFill)
{
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  for(axom::IndexType capacity = 2; capacity < 512; capacity *= 2)
  {
    axom::IndexType size = capacity / 2;
    for(axom::IndexType n_components = 1; n_components <= 4; n_components++)
    {
      MCArray<int> v_int_sidre(root->createView("int"), size, n_components, capacity * n_components);
      internal::check_fill(v_int_sidre);

      MCArray<double> v_double_sidre(root->createView("double"),
                                     size,
                                     n_components,
                                     capacity * n_components);
      internal::check_fill(v_double_sidre);

      root->destroyViewsAndData();
    }
  }
}

//------------------------------------------------------------------------------
TEST(sidre_core_mcarray, checkSet)
{
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  for(axom::IndexType capacity = 2; capacity < 512; capacity *= 2)
  {
    axom::IndexType size = capacity / 2;
    for(axom::IndexType n_components = 1; n_components <= 4; n_components++)
    {
      MCArray<int> v_int_sidre(root->createView("int"), size, n_components, capacity * n_components);
      internal::check_set(v_int_sidre);

      MCArray<double> v_double_sidre(root->createView("double"),
                                     size,
                                     n_components,
                                     capacity * n_components);
      internal::check_set(v_double_sidre);

      root->destroyViewsAndData();
    }
  }
}

//------------------------------------------------------------------------------
TEST(sidre_core_mcarray, checkResize)
{
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr axom::IndexType ZERO = 0;

  for(double ratio = 1.0; ratio <= 2.0; ratio += 0.5)
  {
    for(axom::IndexType capacity = 2; capacity <= 512; capacity *= 2)
    {
      for(axom::IndexType n_components = 1; n_components <= 4; n_components++)
      {
        MCArray<int> v_int_sidre(root->createView("int"), ZERO, n_components, capacity * n_components);
        v_int_sidre.setResizeRatio(ratio);
        internal::check_resize(v_int_sidre);

        MCArray<double> v_double_sidre(root->createView("double"),
                                       ZERO,
                                       n_components,
                                       capacity * n_components);
        v_double_sidre.setResizeRatio(ratio);
        internal::check_resize(v_double_sidre);

        root->destroyViewsAndData();
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(sidre_core_mcarray_DeathTest, checkResize)
{
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  /* Resizing isn't allowed with a ratio less than 1.0. */
  MCArray<int> v_int(root->createView("int"), internal::ZERO, 1, 100);
  v_int.setResizeRatio(0.99);
  EXPECT_DEATH_IF_SUPPORTED(internal::check_resize(v_int), IGNORE_OUTPUT);
}

//------------------------------------------------------------------------------
TEST(sidre_core_mcarray, checkInsert)
{
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr axom::IndexType ZERO = 0;

  for(double ratio = 1.0; ratio <= 2.0; ratio += 0.5)
  {
    for(axom::IndexType capacity = 2; capacity <= 512; capacity *= 2)
    {
      for(axom::IndexType n_components = 1; n_components <= 3; n_components++)
      {
        MCArray<int> v_int_sidre(root->createView("int"), ZERO, n_components, capacity * n_components);
        v_int_sidre.setResizeRatio(ratio);
        internal::check_insert(v_int_sidre);

        MCArray<double> v_double_sidre(root->createView("double"),
                                       ZERO,
                                       n_components,
                                       capacity * n_components);
        v_double_sidre.setResizeRatio(ratio);
        internal::check_insert(v_double_sidre);

        root->destroyViewsAndData();
      }
    }
  }
}

//------------------------------------------------------------------------------
// NOTE: The functionality provided by the old MCArray::emplace() is completely different
// than emplace() as defined on e.g. std::vector
// This test has been disabled for that reason
// TEST(sidre_core_mcarray, checkEmplace)
// {
//   sidre::DataStore ds;
//   sidre::Group* root = ds.getRoot();

//   constexpr axom::IndexType ZERO = 0;

//   for(double ratio = 1.0; ratio <= 2.0; ratio += 0.5)
//   {
//     for(axom::IndexType capacity = 2; capacity <= 512; capacity *= 2)
//     {
//       for(axom::IndexType n_components = 1; n_components <= 3; n_components++)
//       {
//         MCArray<int> v_int_sidre(root->createView("int"),
//                                  ZERO,
//                                  n_components,
//                                  capacity);
//         v_int_sidre.setResizeRatio(ratio);
//         internal::check_emplace(v_int_sidre);

//         MCArray<double> v_double_sidre(root->createView("double"),
//                                        ZERO,
//                                        n_components,
//                                        capacity);
//         v_double_sidre.setResizeRatio(ratio);
//         internal::check_emplace(v_double_sidre);

//         root->destroyViewsAndData();
//       }
//     }
//   }
// }

//------------------------------------------------------------------------------
TEST(sidre_core_mcarray, checkMove)
{
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr axom::IndexType ZERO = 0;

  for(double ratio = 1.0; ratio <= 2.0; ratio += 0.5)
  {
    for(axom::IndexType capacity = 2; capacity <= 512; capacity *= 2)
    {
      for(axom::IndexType n_components = 1; n_components <= 3; n_components++)
      {
        MCArray<int> v_int_sidre(root->createView("int"), ZERO, n_components, capacity * n_components);
        v_int_sidre.setResizeRatio(ratio);
        internal::check_move(v_int_sidre);

        MCArray<double> v_double_sidre(root->createView("double"),
                                       ZERO,
                                       n_components,
                                       capacity * n_components);
        v_double_sidre.setResizeRatio(ratio);
        internal::check_move(v_double_sidre);

        root->destroyViewsAndData();
      }
    }
  }
}

/* Sidre specific tests */

//------------------------------------------------------------------------------
TEST(sidre_core_mcarray, checkSidre)
{
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr axom::IndexType ZERO = 0;

  for(double ratio = 1.0; ratio <= 2.0; ratio += 0.5)
  {
    for(axom::IndexType capacity = 2; capacity <= 512; capacity *= 2)
    {
      for(axom::IndexType n_components = 1; n_components <= 3; n_components++)
      {
        MCArray<int> v_int(root->createView("int"), ZERO, n_components, capacity * n_components);
        v_int.setResizeRatio(ratio);
        internal::check_storage(v_int);
        internal::check_sidre(v_int);

        MCArray<double> v_double(root->createView("double"),
                                 ZERO,
                                 n_components,
                                 capacity * n_components);
        v_double.setResizeRatio(ratio);
        internal::check_storage(v_double);
        internal::check_sidre(v_double);

        root->destroyViewsAndData();
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(sidre_core_mcarray, checkSidrePermanence)
{
  constexpr double MAGIC_NUM = 5683578.8;

  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr axom::IndexType ZERO = 0;

  for(double ratio = 1.0; ratio <= 2.0; ratio += 0.5)
  {
    for(axom::IndexType capacity = 2; capacity <= 512; capacity *= 2)
    {
      for(axom::IndexType n_components = 1; n_components <= 3; n_components++)
      {
        const double* MCArray_data_ptr;
        axom::IndexType num_values;
        { /* Begin scope */
          MCArray<double> v(root->createView("double"), ZERO, n_components, capacity * n_components);
          MCArray_data_ptr = v.data();
          num_values = v.size();
          v.setResizeRatio(ratio);
          internal::check_storage(v);

          /* Set v's data to MAGIC_NUM */
          for(axom::IndexType i = 0; i < v.size() / n_components; ++i)
          {
            for(axom::IndexType j = 0; j < n_components; ++j)
            {
              v(i, j) = MAGIC_NUM;
            }
          }
        } /* End scope, v has been deallocated. */

        /* Check that the data still exists in sidre */
        sidre::View* view = root->getView("double");
        const double* view_data_ptr = static_cast<const double*>(view->getVoidPtr());
        EXPECT_EQ(view_data_ptr, MCArray_data_ptr);
        EXPECT_EQ(view->getNumDimensions(), 2);

        IndexType dims[2];
        view->getShape(2, dims);
        EXPECT_EQ(dims[0], capacity);
        EXPECT_EQ(dims[1], n_components);

        for(axom::IndexType i = 0; i < num_values; ++i)
        {
          EXPECT_EQ(view_data_ptr[i], MAGIC_NUM);
        }

        root->destroyViewsAndData();
      }
    }
  }
}

} /* end namespace sidre */
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
