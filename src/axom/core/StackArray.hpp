// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_STACKARRAY_HPP_
#define AXOM_STACKARRAY_HPP_

#include "axom/config.hpp"       // for compile-time defines
#include "axom/core/Macros.hpp"  // for axom macros
#include "axom/core/Types.hpp"   // for axom types

#include <iostream>

namespace axom
{
/*!
 * \accelerated
 * \class StackArray
 *
 * \brief Provides a wrapper for a compile time sized array, similar to
 *  std::array. This class is needed because NVCC doesn't capture standard
 *  stack arrays in device lambdas. Furthermore we can't use std::array because
 *  it is not host-device decorated.
 *
 * \tparam T the type of the values to hold.
 * \tparam N the number of values in the array.
 */
template <typename T, int N>
struct StackArray
{
  using value_type = T;

  /*!
   * \brief Return size of the array.
   */
  AXOM_HOST_DEVICE
  constexpr static int size() { return N; }

  /*!
   * \brief Accessor, returns a reference to the value at the given index.
   *
   * \param [in] i the index to access.
   */
  /// @{

  AXOM_HOST_DEVICE
  T& operator[](IndexType i) noexcept { return m_data[i]; }

  AXOM_HOST_DEVICE
  constexpr const T& operator[](IndexType i) const noexcept { return m_data[i]; }

  /// @}

  /*!
   * \brief User defined conversion to a raw pointer.
   */
  /// @{

  AXOM_HOST_DEVICE operator T*() noexcept { return &m_data[0]; }

  AXOM_HOST_DEVICE
  constexpr operator const T*() const noexcept { return &m_data[0]; }

  AXOM_HOST_DEVICE T* data() noexcept { return &m_data[0]; }
  AXOM_HOST_DEVICE const T* data() const noexcept { return &m_data[0]; }

  /// @}

  /*!
   * \brief Begin/end iterators
   */
  /// @{

  AXOM_HOST_DEVICE T* begin() noexcept { return &m_data[0]; }
  AXOM_HOST_DEVICE const T* begin() const noexcept { return &m_data[0]; }

  AXOM_HOST_DEVICE T* end() noexcept { return &m_data[0] + N; }
  AXOM_HOST_DEVICE const T* end() const noexcept { return &m_data[0] + N; }

  /// @}

  T m_data[N];
};

/*!
 * \brief Equality comparison operator for StackArray
 *
 * \param [in] lhs left StackArray to compare
 * \param [in] rhs right StackArray to compare
 * \return true if the StackArrays have the same element values
 */
template <typename T, int N>
AXOM_HOST_DEVICE bool operator==(const StackArray<T, N>& lhs, const StackArray<T, N>& rhs)
{
  for(int i = 0; i < N; ++i)
  {
    if(lhs[i] != rhs[i])
    {
      return false;
    }
  }
  return true;
}

/*!
 * \brief Inequality comparison operator for StackArray
 *
 * \param [in] lhs left StackArray to compare
 * \param [in] rhs right StackArray to compare
 * \return true if the StackArrays have different element values
 */
template <typename T, int N>
AXOM_HOST_DEVICE bool operator!=(const StackArray<T, N>& lhs, const StackArray<T, N>& rhs)
{
  return !(lhs == rhs);
}

/*!
 * \brief Less than operator for StackArray
 *
 * \param [in] lhs left StackArray to compare
 * \param [in] rhs right StackArray to compare
 * \return true if \a lhs is lexicographically less than \a rhs, false otherwise
 * \note It is only valid to call this function when values of type \a T are comparable, 
 * e.g. when T has an operator<()
 */
template <typename T, int N>
AXOM_HOST_DEVICE bool operator<(const StackArray<T, N>& lhs, const StackArray<T, N>& rhs)
{
  for(int i = 0; i < N; ++i)
  {
    if(lhs[i] < rhs[i])
    {
      return true;
    }
    else if(lhs[i] > rhs[i])
    {
      return false;
    }
  }
  return false;
}

/**
 * \brief Print the StackArray to a stream.
 * \param os The stream to use.
 * \param obj The StackArray to print.
 * \return The input stream.
 */
template <typename T, int N>
std::ostream& operator<<(std::ostream& os, const StackArray<T, N>& obj)
{
  os << "(";
  for(int i = 0; i < N; i++)
  {
    if(i > 0)
    {
      os << ", ";
    }
    os << obj.m_data[i];
  }
  os << ")";
  return os;
}

} /* namespace axom */

#endif /* AXOM_STACKARRAY_HPP_ */
