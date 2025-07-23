// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_NUMERIC_ARRAY_HPP_
#define AXOM_PRIMAL_NUMERIC_ARRAY_HPP_

#include "axom/core/Macros.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include <cassert>

// C/C++ includes
#include <algorithm>
#include <ostream>
#include <initializer_list>
#include "axom/fmt.hpp"

namespace axom
{
// Forward declare the templated classes and operator functions
template <typename T, int SIZE>
class NumericArray;

/// \name Forward Declared Overloaded Operators
///@{

/*!
 * \brief Checks if two numeric arrays are component-wise equal.
 * \param [in] lhs numeric array instance on the left-hand side.
 * \param [in] rhs numeric array instance on the right-hand side.
 * \return status true if lhs==rhs, otherwise, false.
 */
template <typename T, int SIZE>
AXOM_HOST_DEVICE bool operator==(const NumericArray<T, SIZE>& lhs, const NumericArray<T, SIZE>& rhs);

/*!
 * \brief Checks if two numeric arrays are *not* component-wise equal.
 * \param [in] lhs numeric array instance on the left-hand side.
 * \param [in] rhs numeric array instance on the right-hand side.
 * \return status true if lhs!=rhs, otherwise, false.
 */
template <typename T, int SIZE>
AXOM_HOST_DEVICE bool operator!=(const NumericArray<T, SIZE>& lhs, const NumericArray<T, SIZE>& rhs);

/*!
 * \brief Performs component-wise addition of two numeric arrays.
 * \param [in] lhs numeric array instance on the left-hand side.
 * \param [in] rhs numeric array instance on the right-hand side.
 * \return C resulting numeric array from the component-wise addition.
 */
template <typename T, int SIZE>
AXOM_HOST_DEVICE NumericArray<T, SIZE> operator+(const NumericArray<T, SIZE>& lhs,
                                                 const NumericArray<T, SIZE>& rhs);

/*!
 * \brief Performs component-wise subtraction of two numeric arrays.
 * \param [in] lhs numeric array instance on the left-hand side.
 * \param [in] rhs numeric array instance on the right-hand side.
 * \result C resulting numeric array from component-wise subtraction.
 */
template <typename T, int SIZE>
AXOM_HOST_DEVICE NumericArray<T, SIZE> operator-(const NumericArray<T, SIZE>& lhs,
                                                 const NumericArray<T, SIZE>& rhs);

/*!
 * \brief Unary negation of a numeric array instance.
 * \param [in] arr numeric array instance on the left-hand side.
 * \result C resulting numeric array from unary negation.
 */
template <typename T, int SIZE>
AXOM_HOST_DEVICE NumericArray<T, SIZE> operator-(const NumericArray<T, SIZE>& arr);

/*!
 * \brief Scalar multiplication a numeric array; Scalar on rhs.
 * \param [in] arr numeric array instance.
 * \param [in] scalar user-supplied scalar.
 * \return C resutling numeric array, \f$ \ni: C_i = scalar*arr_i, \forall i\f$
 */
template <typename T, int SIZE>
AXOM_HOST_DEVICE NumericArray<T, SIZE> operator*(const NumericArray<T, SIZE>& arr, double scalar);

/*!
 * \brief Scalar multiplication a numeric array; Scalar on lhs.
 * \param [in] scalar user-supplied scalar.
 * \param [in] arr numeric array instance.
 * \return C resulting numeric array, \f$ \ni: C_i = scalar*arr_i, \forall i\f$
 */
template <typename T, int SIZE>
AXOM_HOST_DEVICE NumericArray<T, SIZE> operator*(double scalar, const NumericArray<T, SIZE>& arr);

/*!
 * \brief Component-wise multiplication of NumericArrays
 * \param [in] lhs numeric array instance on the left-hand side.
 * \param [in] rhs numeric array instance on the right-hand side.
 * \return C resulting numeric array, \f$ \ni: C_i = lhs_i * rhs_i, \forall i\f$
 */
template <typename T, int SIZE>
AXOM_HOST_DEVICE NumericArray<T, SIZE> operator*(const NumericArray<T, SIZE>& lhs,
                                                 const NumericArray<T, SIZE>& rhs);

/*!
 * \brief Component-wise division of NumericArrays
 * \param [in] lhs numeric array instance on the left-hand side.
 * \param [in] rhs numeric array instance on the right-hand side.
 * \return C resulting numeric array, \f$ \ni: C_i = lhs_i / rhs_i, \forall i\f$
 * \pre \f$ rhs_i != 0.0, \forall i \f$
 */
template <typename T, int SIZE>
AXOM_HOST_DEVICE NumericArray<T, SIZE> operator/(const NumericArray<T, SIZE>& lhs,
                                                 const NumericArray<T, SIZE>& rhs);

/*!
 * \brief Scalar division of NumericArray; Scalar on rhs
 * \param [in] arr numeric array instance
 * \param [in] scalar user-supplied scalar
 * \return C resulting numeric array, \f$ \ni: C_i = arr_i/scalar, \forall i\f$
 * \pre scalar != 0.0
 */
template <typename T, int SIZE>
AXOM_HOST_DEVICE NumericArray<T, SIZE> operator/(const NumericArray<T, SIZE>& arr, double scalar);

/*!
 * \brief Coordinate-wise absolute value on the NumericArray
 * \param [in] arr numeric array instance
 * \pre std::abs is defined for template type T
 * \return A NumericArray whose coordinates are the absolute value of arr
 */
template <typename T, int SIZE>
AXOM_HOST_DEVICE NumericArray<T, SIZE> abs(const NumericArray<T, SIZE>& arr);

/*!
 * \brief Overloaded output operator for numeric arrays
 * \param [in] os C++ output stream
 * \param [in] arr numeric array instance.
 */
template <typename T, int SIZE>
std::ostream& operator<<(std::ostream& os, const NumericArray<T, SIZE>& arr);

///@}

/**
 * \brief Type trait to avoid outputting chars when a value is expected
 *  This avoids unintentionally outputting system beeps
 */
template <typename T>
struct NonChar
{
  using type = T;  // The non-char type to return
};

template <>
struct NonChar<char>
{
  /** A non-char signed type to which we can cast a char for output */
  using type = int;
};

template <>
struct NonChar<unsigned char>
{
  /** A non-char unsigned type to which we can cast a char for output */
  using type = unsigned int;
};

/*!
 * \accelerated
 * \class NumericArray
 *
 * \brief A simple statically sized array of data with component-wise operators.
 *
 * \tparam T the numeric type of the elements in the array, e.g., float, double.
 * \tparam SIZE the size of the array
 */
template <typename T, int SIZE>
class NumericArray  // NOLINT
{
public:
  /*!
   * \brief Fill the first sz coordinates with val and zeros the rest
   * \param [in] val The value to set the coordinates to. Defaults to zero.
   * \param [in] sz The number of components to set to val.
   * The rest will be set to zero.  Defaults is SIZE.
   * If sz is greater than SIZE, we set all coordinates to val
   */
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE explicit NumericArray(T val = T(), int sz = SIZE);

  /*!
   * \brief Creates a numeric array from the first sz values of the input array.
   * \param [in] vals An array containing at least sz values
   * \param [in] sz number of coordinates. Defaults to SIZE.
   * \note If sz is greater than SIZE, we only take the first SIZE values.
   */
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE explicit NumericArray(const T* vals, int sz = SIZE);

  /*!
   * \brief Creates a numeric array from an initializer list
   * \param [in] values an initializer list containing the values of the
   * array. If the size is not the same as the size of this array, this
   * behaves the same way as the constructor which takes a pointer and size.
   */
  AXOM_HOST_DEVICE
  NumericArray(std::initializer_list<T> values)
    : NumericArray {values.begin(), static_cast<int>(values.size())}
  { }

  /*!
   * \brief Returns the dimension of this numeric array instance.
   * \return d the dimension (size) of the array
   * \post d >= 1.
   */
  static int size() { return SIZE; };

  /*!
   * \brief Access operator for individual components.
   * \param [in] i the component index to access
   * \return \f$ p_i \f$ the value at the given component index.
   * \pre \f$  0 \le i < SIZE \f$
   */
  AXOM_HOST_DEVICE const T& operator[](int i) const;

  AXOM_HOST_DEVICE T& operator[](int i);

  /*!
   * \brief Returns a pointer to the underlying data.
   */
  AXOM_HOST_DEVICE const T* data() const;

  AXOM_HOST_DEVICE T* data();

  /*!
   * \brief Copy the coordinate data to the provided array
   * \param [in] arr The array to which we are copying.
   * \pre The user needs to make sure that the provided array has been allocated
   * and has sufficient space for SIZE coordinates.
   */
  AXOM_HOST_DEVICE void to_array(T* arr) const;

  /*!
   * \brief Simple formatted print of a numeric array instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const;

  /*!
   * \brief Component-wise addition assignment operator.
   * \param [in] arr the array to add.
   * Adds the numeric array arr to this instance (component-wise).
   * \return A reference to the NumericArray instance after addition.
   */
  AXOM_HOST_DEVICE NumericArray& operator+=(const NumericArray& arr);

  /*!
   * \brief Component-wise subtraction assignment operator.
   * \param [in] arr the array to subtract.
   * Subtracts the numeric array arr from this instance (component-wise).
   * \return A reference to the NumericArray instance after subtraction.
   */
  AXOM_HOST_DEVICE NumericArray& operator-=(const NumericArray& arr);

  /*!
   * \brief Scalar multiplication on the NumericArray instance.
   * \param [in] scalar the scalar value with which to multiply.
   * Each element of the numeric array is multiplied by scalar
   * \return A reference to the NumericArray instance after scalar
   * multiplication.
   */
  AXOM_HOST_DEVICE NumericArray& operator*=(double scalar);

  /*!
   * \brief Scalar division on the NumericArray instance.
   * \param [in] scalar the scalar value with which to divide .
   * \pre scalar != 0
   * Each element of the numeric array is divided by scalar
   * \return A reference to the NumericArray instance after scalar division.
   */
  AXOM_HOST_DEVICE NumericArray& operator/=(double scalar);

  /*!
   * \brief Component-wise multiplication assignment operator.
   * \param [in] arr the array to multiply (component-wise).
   * Multiplies the numeric array arr with this instance (component-wise).
   * \return A reference to the NumericArray instance after cwise
   * multiplication.
   */
  AXOM_HOST_DEVICE NumericArray& operator*=(const NumericArray& arr);

  /*!
   * \brief Component-wise division assignment operator.
   * \param [in] arr the array to divide (component-wise).
   * Divides the numeric array arr with this instance (component-wise).
   * \pre forall i, arr[i] != 0
   * \return A reference to the NumericArray instance after cwise division.
   */
  AXOM_HOST_DEVICE NumericArray& operator/=(const NumericArray& arr);

  /*!
   * \brief Ensures that the highest value of the coordinates is at most
   *  upperVal.
   *
   * \param [in] upperVal The highest possible value
   * \post forall i, arr[i] <= upperVal
   * \return A reference to the NumericArray instance after clamping upper
   */
  NumericArray& clampUpper(const T& upperVal);

  /*!
   * \brief Ensures that the lowest value of the coordinates is at least
   *  lowerVal.
   *
   * \param [in] lowerVal The lowest possible value
   *
   * \post forall i, arr[i] >= lowerVal
   *
   * \return A reference to the NumericArray instance after clamping lower
   */
  NumericArray& clampLower(const T& lowerVal);

  /*!
   * \brief Ensures that each coordinate's value is in range
   *  [lowerVal,upperVal].
   *
   * \param [in] lowerVal The lowest possible value
   * \param [in] upperVal The highest possible value
   *
   * \pre lowerVal <= upperVal
   * \post forall i, lowerVal <= arr[i] <= upperVal
   *
   * \return A reference to the NumericArray instance after clamping
   */
  NumericArray& clamp(const T& lowerVal, const T& upperVal);

  /*!
   * \brief Find the max component.
   * \return The value of the largest component.
   */
  AXOM_HOST_DEVICE T max() const;

  /*!
   * \brief Find the min component.
   * \return The value of the smallest component.
   */
  AXOM_HOST_DEVICE T min() const;

  /*!
   * \brief Find the index of the max component.
   * \return The index of the largest component ( \f$ 0 \le ret < SIZE \f$)
   */
  AXOM_HOST_DEVICE int argMax() const;

  /*!
   * \brief Find the index of the min component.
   * \return The index of the smallest component ( \f$ 0 \le ret < SIZE \f$)
   */
  AXOM_HOST_DEVICE int argMin() const;

  /// \brief Computes the sum of the components
  AXOM_HOST_DEVICE T sum() const;

private:
  AXOM_HOST_DEVICE void verifyIndex(int idx) const
  {
    AXOM_UNUSED_VAR(idx);
    assert(idx >= 0 && idx < SIZE);
  }

protected:
  T m_components[SIZE];  /// The encapsulated array
};

}  // namespace axom

//------------------------------------------------------------------------------
//  NumericArray implementation
//------------------------------------------------------------------------------

namespace axom
{
//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE NumericArray<T, SIZE>::NumericArray(T val, int sz)
{
  // NOTE (KW): This should be a static assert in the class
  assert(SIZE >= 1);

  // Fill first nvals coordinates with val ( 0 <= nvals <= SIZE )
  const int nvals = axom::utilities::clampVal(sz, 0, SIZE);
  for(int i = 0; i < nvals; i++)
  {
    m_components[i] = val;
  }

  // Fill any remaining coordinates with zero
  for(int j = nvals; j < SIZE; j++)
  {
    m_components[j] = T();
  }
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE NumericArray<T, SIZE>::NumericArray(const T* vals, int sz)
{
  assert(SIZE >= 1);

  const int nvals = axom::utilities::clampVal(sz, 0, SIZE);

  // Copy first nvals coordinates from vals array ( 0 <= nvals <= SIZE )
  for(int i = 0; i < nvals; i++)
  {
    m_components[i] = vals[i];
  }

  // Fill any remaining coordinates with zero
  for(int j = nvals; j < SIZE; j++)
  {
    m_components[j] = T();
  }
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline T& NumericArray<T, SIZE>::operator[](int i)
{
  verifyIndex(i);
  return m_components[i];
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline const T& NumericArray<T, SIZE>::operator[](int i) const
{
  verifyIndex(i);
  return m_components[i];
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline const T* NumericArray<T, SIZE>::data() const
{
  return m_components;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline T* NumericArray<T, SIZE>::data()
{
  return m_components;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE void NumericArray<T, SIZE>::to_array(T* arr) const
{
  assert(arr != nullptr);
  for(int dim = 0; dim < SIZE; ++dim)
  {
    arr[dim] = m_components[dim];
  }
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
std::ostream& NumericArray<T, SIZE>::print(std::ostream& os) const
{
  os << "[ ";
  for(int dim = 0; dim < SIZE - 1; ++dim)
  {
    os << static_cast<typename NonChar<T>::type>(m_components[dim]) << " ";
  }

  os << static_cast<typename NonChar<T>::type>(m_components[SIZE - 1]) << "]";

  return os;
}

//------------------------------------------------------------------------------
// Member function arithmetic operators (component-wise)
//------------------------------------------------------------------------------

template <typename T, int SIZE>
AXOM_HOST_DEVICE inline NumericArray<T, SIZE>& NumericArray<T, SIZE>::operator*=(double scalar)
{
  for(int i = 0; i < SIZE; ++i)
  {
    m_components[i] = static_cast<T>(m_components[i] * scalar);
  }

  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline NumericArray<T, SIZE>& NumericArray<T, SIZE>::operator/=(double scalar)
{
  assert(scalar != 0.);
  return operator*=(1. / scalar);
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline NumericArray<T, SIZE>& NumericArray<T, SIZE>::operator*=(
  const NumericArray<T, SIZE>& v)
{
  for(int i = 0; i < SIZE; ++i)
  {
    m_components[i] *= v[i];
  }

  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline NumericArray<T, SIZE>& NumericArray<T, SIZE>::operator/=(
  const NumericArray<T, SIZE>& v)
{
  for(int i = 0; i < SIZE; ++i)
  {
    assert(v[i] != 0.);
    m_components[i] /= v[i];
  }

  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline NumericArray<T, SIZE>& NumericArray<T, SIZE>::operator+=(
  const NumericArray<T, SIZE>& v)
{
  for(int i = 0; i < SIZE; ++i)
  {
    m_components[i] += v[i];
  }

  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline NumericArray<T, SIZE>& NumericArray<T, SIZE>::operator-=(
  const NumericArray<T, SIZE>& v)
{
  for(int i = 0; i < SIZE; ++i)
  {
    m_components[i] -= v[i];
  }

  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
inline NumericArray<T, SIZE>& NumericArray<T, SIZE>::clamp(const T& lowerVal, const T& upperVal)
{
  assert(lowerVal <= upperVal);

  for(int i = 0; i < SIZE; ++i)
  {
    m_components[i] = axom::utilities::clampVal(m_components[i], lowerVal, upperVal);
  }

  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
inline NumericArray<T, SIZE>& NumericArray<T, SIZE>::clampLower(const T& lowerVal)
{
  for(int i = 0; i < SIZE; ++i)
  {
    m_components[i] = std::max(m_components[i], lowerVal);
  }

  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
inline NumericArray<T, SIZE>& NumericArray<T, SIZE>::clampUpper(const T& upperVal)
{
  for(int i = 0; i < SIZE; ++i)
  {
    m_components[i] = std::min(m_components[i], upperVal);
  }

  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline T NumericArray<T, SIZE>::max() const
{
  T result = this->m_components[0];
  for(int i = 1; i < SIZE; ++i)
  {
    T tmp = m_components[i];

    if(tmp > result)
    {
      result = tmp;
    }
  }

  return result;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline T NumericArray<T, SIZE>::min() const
{
  T result = this->m_components[0];
  for(int i = 1; i < SIZE; ++i)
  {
    T tmp = this->m_components[i];

    if(tmp < result)
    {
      result = tmp;
    }
  }

  return result;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline int NumericArray<T, SIZE>::argMax() const
{
  int idx = 0;
  for(int i = 1; i < SIZE; ++i)
  {
    if(m_components[i] > m_components[idx])
    {
      idx = i;
    }
  }

  return idx;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline int NumericArray<T, SIZE>::argMin() const
{
  int idx = 0;
  for(int i = 1; i < SIZE; ++i)
  {
    if(m_components[i] < m_components[idx])
    {
      idx = i;
    }
  }

  return idx;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline T NumericArray<T, SIZE>::sum() const
{
  T result {};
  for(int i = 0; i < SIZE; ++i)
  {
    result += this->m_components[i];
  }

  return result;
}

//------------------------------------------------------------------------------
/// Free functions implementing comparison and arithmetic operators
//------------------------------------------------------------------------------

template <typename T, int SIZE>
AXOM_HOST_DEVICE bool operator==(const NumericArray<T, SIZE>& lhs, const NumericArray<T, SIZE>& rhs)
{
  for(int dim = 0; dim < SIZE; ++dim)
  {
    if(lhs[dim] != rhs[dim])
    {
      return false;
    }
  }

  return true;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE bool operator!=(const NumericArray<T, SIZE>& lhs, const NumericArray<T, SIZE>& rhs)
{
  return !(lhs == rhs);
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
std::ostream& operator<<(std::ostream& os, const NumericArray<T, SIZE>& arr)
{
  arr.print(os);
  return os;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline NumericArray<T, SIZE> operator*(const NumericArray<T, SIZE>& arr,
                                                        double scalar)
{
  NumericArray<T, SIZE> result(arr);
  result *= scalar;
  return result;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline NumericArray<T, SIZE> operator*(double scalar,
                                                        const NumericArray<T, SIZE>& arr)
{
  NumericArray<T, SIZE> result(arr);
  result *= scalar;
  return result;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline NumericArray<T, SIZE> operator+(const NumericArray<T, SIZE>& lhs,
                                                        const NumericArray<T, SIZE>& rhs)
{
  NumericArray<T, SIZE> result(lhs);
  result += rhs;
  return result;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline NumericArray<T, SIZE> operator*(const NumericArray<T, SIZE>& lhs,
                                                        const NumericArray<T, SIZE>& rhs)
{
  NumericArray<T, SIZE> result(lhs);
  result *= rhs;
  return result;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline NumericArray<T, SIZE> operator/(const NumericArray<T, SIZE>& lhs,
                                                        const NumericArray<T, SIZE>& rhs)
{
  NumericArray<T, SIZE> result(lhs);
  result /= rhs;
  return result;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline NumericArray<T, SIZE> operator/(const NumericArray<T, SIZE>& arr,
                                                        double scalar)
{
  NumericArray<T, SIZE> result(arr);
  result /= scalar;
  return result;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline NumericArray<T, SIZE> operator-(const NumericArray<T, SIZE>& lhs,
                                                        const NumericArray<T, SIZE>& rhs)
{
  NumericArray<T, SIZE> result(lhs);
  result -= rhs;
  return result;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline NumericArray<T, SIZE> operator-(const NumericArray<T, SIZE>& arr)
{
  NumericArray<T, SIZE> result;
  result -= arr;
  return result;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline NumericArray<T, SIZE> abs(const NumericArray<T, SIZE>& arr)
{
  NumericArray<T, SIZE> result(arr);

  for(int i = 0; i < SIZE; ++i)
  {
    result[i] = axom::utilities::abs(result[i]);
  }

  return result;
}

}  // namespace axom

/// Overload to format a axom::NumericArray using fmt
template <typename T, int NDIMS>
struct axom::fmt::formatter<axom::NumericArray<T, NDIMS>> : ostream_formatter
{ };

#endif  // AXOM_PRIMAL_NUMERIC_ARRAY_HPP_
