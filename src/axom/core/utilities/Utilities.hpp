// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *
 * \file Utilities.hpp
 *
 * \brief Header file containing utility functions.
 *
 */

#ifndef AXOM_UTILITIES_HPP_
#define AXOM_UTILITIES_HPP_

#include "axom/config.hpp"  // for compile-time definitions
#include "axom/core/Types.hpp"
#include "axom/core/Macros.hpp"  // for AXOM_STATIC_ASSERT

#include <cassert>  // for assert()
#include <cmath>    // for log2()

#include <random>       // for random  number generator
#include <type_traits>  // for std::is_floating_point()

namespace axom
{
namespace utilities
{
/*!
 * \brief Gracefully aborts the application
 */
void processAbort();

/*!
 * \brief Returns the absolute value of x.
 * \accelerated
 * \param [in] x value whose absolute value is computed.
 * \return abs(x) the absolute value of x.
 */
template <typename T>
inline AXOM_HOST_DEVICE T abs(const T& x)
{
  return (x < T(0)) ? -x : x;
}

/*!
 * \brief Returns the largest integer less than or equal to x.
 * \accelerated
 * \param [in] x value whose floor value is computed.
 * \return floor(x) the largest integer less than or equal to x.
 */
template <typename T>
inline AXOM_HOST_DEVICE T floor(const T& x)
{
  return ::floor(x);
}

/*!
 * \brief Returns the smallest integer greater than or equal to x.
 * \accelerated
 * \param [in] x value whose ceil value is computed.
 * \return ceil(x) the smallest integer greater than or equal to x.
 */
template <typename T>
inline AXOM_HOST_DEVICE T ceil(const T& x)
{
  return ::ceil(x);
}

/*!
 * \brief Returns the max value of x and y.
 * \accelerated
 * \param [in] x the first value to check.
 * \param [in] y the second value to check.
 * \return max(x, y) the max value of x and y.
 */
template <typename T>
inline AXOM_HOST_DEVICE const T& max(const T& x, const T& y)
{
  return (y < x) ? x : y;
}

/*!
 * \brief Returns the min value of x and y.
 * \accelerated
 * \param [in] x the first value to check.
 * \param [in] y the second value to check.
 * \return min(x, y) the min value of x and y.
 */
template <typename T>
inline AXOM_HOST_DEVICE const T& min(const T& x, const T& y)
{
  return (y < x) ? y : x;
}

/*!
 * \brief Swaps the values of a, b.
 * \accelerated
 * \param [in,out] a 1st object to swap.
 * \param [in,out] b 2nd object to swap.
 */
template <typename T>
inline AXOM_HOST_DEVICE void swap(T& a, T& b)
{
  T tmp = a;
  a = b;
  b = tmp;
}

/*!
 * \brief Returns the fused-multiply-add of a, b, and c
 * \accelerated
 * \param [in] a,b,c the values
 * \return a*b+c computed using an fma operation
 */
template <typename T>
inline AXOM_HOST_DEVICE T fma(const T& a, const T& b, const T& c)
{
#if defined(AXOM_DEVICE_CODE)
  return ::fma(a, b, c);
#else
  return std::fma(a, b, c);
#endif
}

/*!
 * \brief Returns the base 2 logarithm of the input.
 * \param [in] val The input value
 */
template <typename T>
inline T log2(T val)
{
  return static_cast<T>(std::log2(val));
}

/*!
 * \brief Linearly interpolates between two values
 * \param [in] val0 The first value
 * \param [in] val2 The second value
 * \param [in] t The interpolation parameter.
 * \return The interpolated value
 */
template <typename T>
inline AXOM_HOST_DEVICE T lerp(T v0, T v1, T t)
{
  constexpr T one = T(1);
  return (one - t) * v0 + t * v1;
}

/*!
 * \brief Clamps an input value to a given range.
 * \accelerated
 * \param [in] val  The value to clamp.
 * \param [in] lower The lower range.
 * \param [in] upper The upper range.
 * \return The clamped value.
 * \pre lower <= upper
 * \post lower <= returned value <= upper.
 */
template <typename T>
inline AXOM_HOST_DEVICE T clampVal(T val, T lower, T upper)
{
  return (val < lower) ? lower : (val > upper) ? upper : val;
}

/*!
 * \brief Clamps the upper range on an input value
 * \accelerated
 * \param [in] val The value to clamp
 * \param [in] upper The upper range
 * \return upper if val > upper, else val
 * \post returned value is less than or equal to upper
 */
template <typename T>
inline AXOM_HOST_DEVICE T clampUpper(T val, T upper)
{
  return val > upper ? upper : val;
}

/*!
 * \brief Clamps the lower range on an input value
 * \accelerated
 * \param [in] val The value to clamp
 * \param [in] lower The lower range
 * \return lower if val < lower, else val
 * \post returned value is greater than or equal to lower
 */
template <typename T>
inline AXOM_HOST_DEVICE T clampLower(T val, T lower)
{
  return val < lower ? lower : val;
}

/*!
 * \brief Computes the binomial coefficient `n choose k`
 *
 * \return \f$ {n\choose k}  = n! / (k! * (n-k)!)\f$
 * when \f$ n \ge k \ge 0 \f$, 0 otherwise.
 */
int binomialCoefficient(int n, int k);

/*!
 * \brief Returns a random real number within the specified interval
 *
 * \param [in] a the interval's lower bound
 * \param [in] b the interval's upper bound
 * \return r the random real number within the specified interval.
 *
 * \tparam T a built-in floating point type, e.g., double, float, long double.
 *
 * \note Consecutive calls to this method will generate a non-deterministic
 *  sequence of numbers.
 *
 * \pre a < b
 * \post a <= r < b
 */
template <typename T>
inline T random_real(const T& a, const T& b)
{
  AXOM_STATIC_ASSERT(std::is_floating_point<T>::value);
  assert((a < b) && "invalid bounds, a < b");

  static std::random_device rd;
  static std::mt19937_64 mt(rd());
  static std::uniform_real_distribution<T> dist(0.0, 1.0);

  T temp = dist(mt);
  return temp * (b - a) + a;
}

/*!
 * \brief Returns a random real number within the specified interval given
 *  the bounds of the interval and a seed value for the underlying random
 *  number generator.
 *
 * \param [in] a the interval's lower bound
 * \param [in] b the interval's upper bound
 * \param [in] seed user-supplied seed for the random number generator
 * \return r the random real number within the specified interval.
 *
 * \tparam T a built-in floating type, e.g., double, float, long double.
 *
 * \note Consecutive calls to this method will generate a deterministic
 *  sequence of numbers.
 *
 * \pre a < b
 * \post a <= r < b
 */
template <typename T>
inline T random_real(const T& a, const T& b, unsigned int seed)
{
  AXOM_STATIC_ASSERT(std::is_floating_point<T>::value);
  assert((a < b) && "invalid bounds, a < b");

  static std::mt19937_64 mt(seed);
  static std::uniform_real_distribution<T> dist(0.0, 1.0);

  double temp = dist(mt);
  return temp * (b - a) + a;
}

/*!
 * \brief Tests the endianness of the system.
 * \return True, if the system is little endian, false otherwise.
 */
inline bool isLittleEndian()
{
  // Note: Endianness test adapted from: https://stackoverflow.com/a/2103095

  enum
  {
    O32_LITTLE_ENDIAN = 0x03020100ul,
    O32_BIG_ENDIAN = 0x00010203ul,
    O32_PDP_ENDIAN = 0x01000302ul
  };

  const union
  {
    std::uint8_t raw[4];
    std::uint32_t value;
  } host_order = {{0, 1, 2, 3}};

  return host_order.value == O32_LITTLE_ENDIAN;
}

/*!
 * \brief Reverses the bytes of the input value.
 * \param [in] val The input value
 * \return The value with swapped bytes
 * \note Assumes endianness is either little or big (not PDP).
 * \pre T is a native arithmetic type (i.e. integral or floating point).
 * \pre sizeof(T) must be 1, 2, 4, or 8 bytes.
 */
template <typename T>
constexpr T byteswap(T val) noexcept
{
  constexpr int NBYTES = sizeof(T);

  AXOM_STATIC_ASSERT_MSG(NBYTES == 1 || NBYTES == 2 || NBYTES == 4 || NBYTES == 8,
                         "byteswap only valid for types of size 1, 2, 4 or 8 bytes.");

  AXOM_STATIC_ASSERT_MSG(std::is_arithmetic<T>::value,
                         "byteswap only valid for native arithmetic types");

  union
  {
    std::uint8_t raw[NBYTES];
    T val;
  } swp;

  std::uint8_t* src = reinterpret_cast<std::uint8_t*>(&val);

  // Reverse the bytes
  for(int i = 0; i < NBYTES; ++i)
  {
    swp.raw[i] = src[(NBYTES - 1) - i];
  }

  return swp.val;
}

/*!
 * \brief Fuzzy comparison of two real valued quantities.
 * \accelerated
 * \param [in] a The first real valued quantities we are comparing.
 * \param [in] b The second real valued quantities we are comparing.
 * \param [in] thresh The threshold of the fuzzy comparison.  Default is 1.0e-8.
 * \return True if the absolute value of the difference is less than thresh and
 *  false otherwise.
 */
template <typename RealType>
inline AXOM_HOST_DEVICE bool isNearlyEqual(RealType a, RealType b, RealType thresh = 1.0e-8)
{
  return abs(a - b) <= thresh;
}

/*!
 * \brief Fuzzy comparison of two real valued quantities.
 * \accelerated
 * \param [in] a The first real valued quantities we are comparing.
 * \param [in] b The second real valued quantities we are comparing.
 * \param [in] relThresh The relative threshold of the fuzzy comparison.
 *  Default is 1.0e-6.
 * \param [in] absThresh The absolute threshold of the fuzzy comparison.
 *  Default is 1.0e-8.
 * \return True if the absolute value of the difference is less than the sum of
 *  absThresh and the relative difference (relThresh times the absolute max of
 *  a and b).
 */
template <typename RealType>
inline AXOM_HOST_DEVICE bool isNearlyEqualRelative(RealType a,
                                                   RealType b,
                                                   RealType relThresh = 1.0e-6,
                                                   RealType absThresh = 1.0e-8)
{
  RealType maxFabs = max(abs(a), abs(b));
  return abs(a - b) <= (maxFabs * relThresh + absThresh);

  // Equation from Real-Time Collision Detection book --
  // http://realtimecollisiondetection.net/pubs/Tolerances/
  // Note: If we use this, we must update the doxygen
  // return abs(a-b) <= max(absThresh, relThresh * maxFabs );
}

/*!
 * \brief Sign of a value of any type that supports comparison and
 * negation operators.
 * \return 0 for v in [-eps, eps], -1 for v < eps and +1 for v > eps
 */
template <typename T>
inline int sign_of(const T& v, const T& eps = {0})
{
  assert(eps >= 0);
  return v > eps ? 1 : v < -eps ? -1 : 0;
}

/*!
 * \brief Insertion sort of an array.
 * \accelerated
 * \param [in] array The array to sort.
 * \param [in] n The number of entries in the array.
 * \param [in] cmp The comparator to use for comparing elements; "less than"
 *  by default.
 */
AXOM_SUPPRESS_HD_WARN
template <typename DataType, typename Predicate = std::less<DataType>>
inline AXOM_HOST_DEVICE void insertionSort(DataType* array, IndexType n, Predicate cmp = {})
{
  for(int i = 1; i < n; i++)
  {
    int j = i;
    // Keep swapping elements until we're not out-of-order.
    while(j > 0 && cmp(array[j], array[j - 1]))
    {
      axom::utilities::swap(array[j], array[j - 1]);
      j--;
    }
  }
}

/*!
 * \brief Compares std::vector< T > in lexicographic order
 *
 * \note T must support > and < operators.
 */
template <typename T>
class LexiComparator
{
public:
  bool operator()(const std::vector<T>& v1, const std::vector<T>& v2) const
  {
    size_t s1 = v1.size();
    size_t s2 = v2.size();
    size_t size = std::min(s1, s2);
    size_t i = 0;

    while(i < size)
    {
      if(v1[i] < v2[i])
      {
        return true;
      }
      else if(v1[i] > v2[i])
      {
        return false;
      }

      ++i;
    }

    if(s1 < s2)
    {
      return true;
    }

    return false;
  }
};

//------------------------------------------------------------------------------
/*!
 * \brief Use binary search to find the index of the \a value in the supplied
 *        sorted container.
 *
 * \tparam ContainerT A container or view of T. It must support size(), operator[].
 * \tparam T The type of elements we're handling.
 *
 * \param[in] cont A view or container that contains the sorted search data values.
 * \param[in] value The search value.
 *
 * \return The index where value was located in view or -1 if not found. If the
 *         container contained repeats of \a value, the first index is returned.
 */
template <typename ContainerT, typename T>
AXOM_HOST_DEVICE std::int32_t binary_search(const ContainerT& cont, T value)
{
  std::int32_t index = -1;
  std::int32_t left = 0;
  std::int32_t right = cont.size() - 1;
  while(left <= right)
  {
    std::int32_t m = (left + right) / 2;
    if(cont[m] < value)
      left = m + 1;
    else if(cont[m] > value)
      right = m - 1;
    else
    {
      index = m;
      // Keep searching to the left to find the first value.
      right = m - 1;
    }
  }

  return index;
}

//------------------------------------------------------------------------------
/*!
 * \brief Hash a stream of bytes into a uint64_t hash value.
 *
 * \param[in] data The bytes to be hashed.
 * \param[in] length The number of bytes to be hashed.
 *
 * \return A uint64_t hash for the byte stream.
 *
 * \note The byte stream is hashed using a Jenkins-hash algorithm forwards and
 *       backwards and the two results are merged into a uint64_t. The length is
 *       also part of the hash to guard against a lot of repeated values in the
 *       byte stream hashing to the same thing.
 *
 * \note We make this function inline since it is not a template and we want to
 *       use it in both host and device code.
 */
AXOM_HOST_DEVICE
inline std::uint64_t hash_bytes(const std::uint8_t* data, std::uint32_t length)
{
  std::uint32_t hash = 0;

  // Build the length into the hash.
  const auto ldata = reinterpret_cast<const std::uint8_t*>(&length);
  for(int e = 0; e < 4; e++)
  {
    hash += ldata[e];
    hash += hash << 10;
    hash ^= hash >> 6;
  }

  std::uint32_t hashr = hash;
  for(std::uint32_t i = 0; i < length; i++)
  {
    hash += data[i];
    hash += hash << 10;
    hash ^= hash >> 6;

    hashr += data[length - 1 - i];
    hashr += hashr << 10;
    hashr ^= hashr >> 6;
  }
  hash += hash << 3;
  hash ^= hash >> 11;
  hash += hash << 15;

  hashr += hashr << 3;
  hashr ^= hashr >> 11;
  hashr += hashr << 15;

  return (static_cast<std::uint64_t>(hash) << 32) | hashr;
}

}  // namespace utilities
}  // namespace axom

#endif  // AXOM_UTILITIES_HPP_
