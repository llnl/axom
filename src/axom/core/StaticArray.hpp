// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_STATICARRAY_HPP_
#define AXOM_STATICARRAY_HPP_

#include "axom/config.hpp"           // for compile-time defines
#include "axom/core/Macros.hpp"      // for axom macros
#include "axom/core/StackArray.hpp"  // for StackArray

namespace axom
{
/*!
 * \accelerated
 * \class StaticArray
 *
 * \brief This class extends StackArray with some std::vector-like convenience methods.
 *
 *
 * \tparam T the type of the values to hold.
 * \tparam N the number of values in the array.
 *
 * \note Type \a T must be default-constructible on device for device
 *       execution.
 */
template <typename T, int N>
class StaticArray : public StackArray<T, N>
{
public:
  /*!
   * \brief Constructor
   */
  AXOM_HOST_DEVICE StaticArray() : StackArray<T, N>(), m_size(0)
  {
#if defined(AXOM_DEBUG)
    for(axom::IndexType i = 0; i < m_size; i++)
    {
      StackArray<T, N>::m_data[i] = T {};
    }
#endif
  }

  /*!
   * \brief Copy Constructor
   * \param obj The object to be copied.
   */
  AXOM_HOST_DEVICE StaticArray(const StaticArray &obj) : StackArray<T, N>(obj), m_size(obj.m_size)
  {
    for(axom::IndexType i = 0; i < obj.m_size; i++)
    {
      StackArray<T, N>::m_data[i] = obj.StackArray<T, N>::m_data[i];
    }
  }

  /*!
   * \brief Destructor.
   */
  AXOM_HOST_DEVICE ~StaticArray() { }

  /*!
   * \brief Copy assignment operator.
   * \param obj The object to be copied.
   */
  AXOM_HOST_DEVICE StaticArray operator=(const StaticArray &obj)
  {
    for(axom::IndexType i = 0; i < obj.m_size; i++)
    {
      StackArray<T, N>::m_data[i] = obj.StackArray<T, N>::m_data[i];
    }
    m_size = obj.m_size;
    return *this;
  }

  /*!
   * \brief Returns the capacity of the static array
   *
   * \return The capacity of the static array
   */
  AXOM_HOST_DEVICE
  constexpr axom::IndexType capacity() const { return static_cast<axom::IndexType>(N); }

  /*!
   * \brief Returns the size of the static array
   *
   * \return The size of the static array
   */
  AXOM_HOST_DEVICE
  axom::IndexType size() const { return m_size; }

  /*!
   * \brief Pushes an object to the back of the static array
   *
   * \param [in] obj the object to be added to the back.
   *
   * \note The number of push_backs must not exceed N,
   *       the max number of values in the array.
   *
   * \note If the static array is full, push_back
   *       will not modify the static array.
   */
  AXOM_HOST_DEVICE
  void push_back(const T &obj)
  {
    assert(m_size < capacity());
    if(m_size < capacity())
    {
      StackArray<T, N>::m_data[m_size++] = obj;
    }
  }

  /*!
   * \brief Pops the last element off the list.
   */
  AXOM_HOST_DEVICE
  void pop_back() { m_size = (m_size > 0) ? (m_size - 1) : 0; }

  /*!
   * \brief Clears the data from the static array
   */
  AXOM_HOST_DEVICE
  void clear()
  {
    m_size = 0;
#if defined(AXOM_DEBUG)
    for(int i = 0; i < N; i++)
    {
      StackArray<T, N>::m_data[i] = T {};
    }
#endif
  }

  /**
   * \brief Determines whether the container is empty.
   * \return True if empty, false otherwise.
   */
  AXOM_HOST_DEVICE
  bool empty() const { return m_size == 0; }

  /*!
   * \brief Fills the container with the supplied value.
   *
   * \param fill_value The fill value.
   */
  AXOM_HOST_DEVICE
  void fill(const T &fill_value)
  {
    for(T &datum : StackArray<T, N>::m_data)
    {
      datum = fill_value;
    }
  }

private:
  axom::IndexType m_size {0};
};

} /* namespace axom */

#endif /* AXOM_STATICARRAY_HPP_ */
