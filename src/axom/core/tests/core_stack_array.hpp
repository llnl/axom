// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/StackArray.hpp"
#include "gtest/gtest.h"
#include <string>

namespace axom
{
namespace internal
{
template <typename T, int N, typename LAMBDA>
void check(const StackArray<T, N>& arr, LAMBDA&& getValue)
{
  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(arr[i], getValue(i));
  }

  const T* ptr = arr;
  EXPECT_EQ(ptr, &arr[0]);

  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(ptr[i], getValue(i));
  }
}

template <typename T, int N, typename LAMBDA>
StackArray<T, N> test_storage(LAMBDA&& getValue)
{
  StackArray<T, N> arr;

  for(int i = 0; i < N; ++i)
  {
    arr[i] = getValue(i);
  }

  check(arr, std::forward<LAMBDA>(getValue));
  return arr;
}

template <typename T, int N, typename LAMBDA>
void test_copy(LAMBDA&& getValue)
{
  StackArray<T, N> arr;

  for(int i = 0; i < N; ++i)
  {
    arr[i] = getValue(i);
  }

  StackArray<T, N> copy(arr);
  EXPECT_NE(&arr[0], &copy[0]);

  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(arr[i], copy[i]);
  }

  check(copy, std::forward<LAMBDA>(getValue));
}

template <typename T, typename LAMBDA>
void test_list_initilization(LAMBDA&& getValue)
{
  const StackArray<T, 5> arr = {{getValue(0), getValue(1), getValue(2), getValue(3), getValue(4)}};

  check(arr, std::forward<LAMBDA>(getValue));
}

template <typename T, int N, typename LAMBDA>
void test_less_than(LAMBDA&& getValue)
{
  StackArray<T, N> arr1, arr2;

  for(int i = 0; i < N; ++i)
  {
    arr1[i] = getValue(i);
    arr2[i] = getValue(i);
  }

  EXPECT_EQ(arr1, arr2);

  // check that operator< returns false when the arrays are equal
  EXPECT_FALSE(arr1 < arr2);
  EXPECT_FALSE(arr2 < arr1);

  // Modify a few values and check for less than
  for(int index : {0, 1, N / 2, N - 1})
  {
    EXPECT_EQ(arr1, arr2);

    T val = arr2[index];

    arr2[index] = val + val;
    EXPECT_TRUE(arr1 < arr2);
    EXPECT_FALSE(arr2 < arr1);
    EXPECT_LT(arr1, arr2);

    arr2[index] = val;
  }
}

} /* namespace internal */

struct Tensor
{
  double x, y, z;

  Tensor() = default;
  Tensor(const Tensor& other) = default;

  explicit Tensor(double val) : x(3 * val), y(3 * val + 1), z(3 * val + 2) { }

  Tensor& operator=(const Tensor& other)
  {
    x = other.x;
    y = other.y;
    z = other.z;
    return *this;
  }

  bool operator==(const Tensor& other) const
  {
#ifndef WIN32
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wfloat-equal"
#endif

    return x == other.x && y == other.y && z == other.z;

#ifndef WIN32
  #pragma GCC diagnostic pop
#endif
  }
};

std::ostream& operator<<(std::ostream& stream, const Tensor& t)
{
  stream << "(" << t.x << ", " << t.y << ", " << t.z << ")";
  return stream;
}

TEST(core_stack_array, storage)
{
  constexpr int N = 100; /* Number of values to store */
  internal::test_storage<int, N>([](int i) { return i; });

  internal::test_storage<Tensor, N>([](int i) { return Tensor(i); });

  internal::test_storage<std::string, N>([](int i) { return std::to_string(i); });
}

TEST(core_stack_array, copy)
{
  constexpr int N = 100; /* Number of values to store */
  internal::test_copy<int, N>([](int i) { return i; });

  internal::test_copy<Tensor, N>([](int i) { return Tensor(i); });

  internal::test_copy<std::string, N>([](int i) { return std::to_string(i); });
}

TEST(core_stack_array, list_initilization)
{
  internal::test_list_initilization<int>([](int i) { return i; });

  internal::test_list_initilization<Tensor>([](int i) { return Tensor(i); });

  internal::test_list_initilization<std::string>([](int i) { return std::to_string(i); });
}

TEST(core_stack_array, less_than)
{
  // Note: Using `i+1` instead of `i` in lambdas
  // due to how `test_less_than()` is implemented

  constexpr int N = 100;  // Number of values to store

  internal::test_less_than<int, N>([](int i) { return i + 1; });

  internal::test_less_than<double, N>([](int i) { return double(i + 1); });

  // The following would produce a compiler error since Tensor does not have an operator<() !
  //internal::test_less_than<Tensor, N>([](int i) { return Tensor(i + 1); });

#if !defined(AXOM_USE_GPU) || !defined(AXOM_GPUCC)
  internal::test_less_than<std::string, N>([](int i) { return std::to_string(i + 1); });
#endif
}

} /* namespace axom */
