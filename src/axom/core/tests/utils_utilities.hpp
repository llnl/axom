// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/utilities/Utilities.hpp"

// C/C++ includes
#include <type_traits>
#include <vector>

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
template <typename T>
void check_random_real(int offset)
{
  AXOM_STATIC_ASSERT(std::is_floating_point<T>::value);

  constexpr T min = 0;
  constexpr T max = 1;

  const T curr_min = min - static_cast<T>(offset);
  const T curr_max = max - static_cast<T>(offset);

  constexpr int MAX_ITERS = 100;
  for(int i = 0; i < MAX_ITERS; ++i)
  {
    T val = axom::utilities::random_real(curr_min, curr_max);
    EXPECT_GE(val, curr_min);
    EXPECT_LE(val, curr_max);
  }
}

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(utils_utilities, log2)
{
  std::cout << "Testing log2 functions." << std::endl;

  // Test integer log2 value of type int
  {
    int val = 64;
    int exp = 6;
    EXPECT_EQ(exp, axom::utilities::log2(val));
  }

  // Test non-integer log2 value of type int
  {
    int val = 72;  // not a power of 2
    int exp = 6;
    EXPECT_EQ(exp, axom::utilities::log2(val));
  }

  // Test integer log2 value of type double
  {
    double val = 16.;
    double exp = 4.;
    EXPECT_EQ(exp, axom::utilities::log2(val));
  }

  // Test non-integer log2 value of type double
  {
    double val = 20.;  // not a power of 2
    double exp = 4.3219281;
    EXPECT_NEAR(exp, axom::utilities::log2(val), 1e-5);
  }
}

//------------------------------------------------------------------------------
TEST(utils_utilities, random_real)
{
  std::cout << "Testing random_real functions (non-deterministic)." << std::endl;

  for(int offset = 0; offset < 10; ++offset)
  {
    check_random_real<float>(offset);
    check_random_real<double>(offset);
    check_random_real<long double>(offset);
  }
}

//------------------------------------------------------------------------------
TEST(utils_utilities, random_real_with_seed)
{
  std::cout << "Testing random_real functions (deterministic)." << std::endl;
  constexpr unsigned int seed = 123456789;

  constexpr double a = -5.0;
  constexpr double b = 5.0;

  for(int i = 0; i < 1000; ++i)
  {
    const double real = axom::utilities::random_real(a, b, seed);
    EXPECT_GE(real, a);
    EXPECT_LE(real, b);
  }
}

//------------------------------------------------------------------------------
TEST(utils_utilities, minmax)
{
  std::cout << "Testing min and max functions." << std::endl;

  // Test simple min, max comparisons on ints
  {
    int a = 5;
    int b = 7;

    int temp_min = axom::utilities::min(a, b);
    int temp_max = axom::utilities::max(a, b);

    EXPECT_EQ(a, temp_min);
    EXPECT_EQ(b, temp_max);
  }

  // Test simple min, max comparisons on doubles
  {
    double a = 5.2;
    double b = -1.7;

    double temp_min = axom::utilities::min(a, b);
    double temp_max = axom::utilities::max(a, b);

    EXPECT_EQ(b, temp_min);
    EXPECT_EQ(a, temp_max);
  }
}

//------------------------------------------------------------------------------
TEST(utils_utilities, floor_ceil)
{
  std::cout << "Testing floor and ceil functions." << std::endl;

  {
    double val = 5.2;
    EXPECT_EQ(5.0, axom::utilities::floor(val));
    EXPECT_EQ(6.0, axom::utilities::ceil(val));
  }

  {
    double val = -5.2;
    EXPECT_EQ(-6.0, axom::utilities::floor(val));
    EXPECT_EQ(-5.0, axom::utilities::ceil(val));
  }

  {
    double val = 5.0;
    EXPECT_EQ(5.0, axom::utilities::floor(val));
    EXPECT_EQ(5.0, axom::utilities::ceil(val));
  }

  {
    double val = -5.0;
    EXPECT_EQ(-5.0, axom::utilities::floor(val));
    EXPECT_EQ(-5.0, axom::utilities::ceil(val));
  }
}

//------------------------------------------------------------------------------
TEST(utils_utilities, binomial_coefficient)
{
  std::cout << "Testing binomial coefficient function." << std::endl;

  // test n less than zero
  {
    const int n = -1;
    const int exp = 0;
    for(int k = -1; k < 10; ++k)
    {
      auto binom_k_n = axom::utilities::binomialCoefficient(n, k);
      EXPECT_EQ(exp, binom_k_n);
    }
  }

  // test n := 0
  {
    const int n = 0;

    EXPECT_EQ(1, axom::utilities::binomialCoefficient(n, 0));

    EXPECT_EQ(0, axom::utilities::binomialCoefficient(n, -1));
    EXPECT_EQ(0, axom::utilities::binomialCoefficient(n, 1));
  }

  // test n := 1
  {
    const int n = 1;

    EXPECT_EQ(0, axom::utilities::binomialCoefficient(n, -1));

    EXPECT_EQ(1, axom::utilities::binomialCoefficient(n, 0));
    EXPECT_EQ(1, axom::utilities::binomialCoefficient(n, 1));

    EXPECT_EQ(0, axom::utilities::binomialCoefficient(n, 2));
  }

  // test n := 2
  {
    const int n = 2;

    EXPECT_EQ(1, axom::utilities::binomialCoefficient(n, 0));
    EXPECT_EQ(2, axom::utilities::binomialCoefficient(n, 1));
    EXPECT_EQ(1, axom::utilities::binomialCoefficient(n, 2));
  }

  // test n := 3
  {
    const int n = 3;

    EXPECT_EQ(1, axom::utilities::binomialCoefficient(n, 0));
    EXPECT_EQ(3, axom::utilities::binomialCoefficient(n, 1));
    EXPECT_EQ(3, axom::utilities::binomialCoefficient(n, 2));
    EXPECT_EQ(1, axom::utilities::binomialCoefficient(n, 3));
  }

  // test n := 4
  {
    const int n = 4;

    EXPECT_EQ(1, axom::utilities::binomialCoefficient(n, 0));
    EXPECT_EQ(4, axom::utilities::binomialCoefficient(n, 1));
    EXPECT_EQ(6, axom::utilities::binomialCoefficient(n, 2));
    EXPECT_EQ(4, axom::utilities::binomialCoefficient(n, 3));
    EXPECT_EQ(1, axom::utilities::binomialCoefficient(n, 4));
  }

  // test recurrence relation  nCk = (n-1)C(k-1) + (n-1)C(k)
  {
    for(int n = 1; n < 10; ++n)
    {
      for(int k = 1; k <= n; ++k)
      {
        auto binom_n_k = axom::utilities::binomialCoefficient(n, k);
        auto binom_n1_k1 = axom::utilities::binomialCoefficient(n - 1, k - 1);
        auto binom_n1_k = axom::utilities::binomialCoefficient(n - 1, k);

        EXPECT_EQ(binom_n_k, binom_n1_k1 + binom_n1_k);
      }
    }
  }
}

TEST(utils_utilities, lerp)
{
  std::cout << "Testing lerp function." << std::endl;

  {
    double A = 0.;
    double B = 1.;

    EXPECT_DOUBLE_EQ(0., axom::utilities::lerp(A, B, 0.));
    EXPECT_DOUBLE_EQ(1., axom::utilities::lerp(A, B, 1.));
    EXPECT_DOUBLE_EQ(.5, axom::utilities::lerp(A, B, .5));

    for(double t = -2.0; t < 2.0; t += 0.05)
    {
      EXPECT_DOUBLE_EQ(t, axom::utilities::lerp(A, B, t));
    }
  }

  // Test endpoint interpolation
  {
    const double lower = -.23;  // Arbitrary end points
    const double upper = 5.73;
    for(int i = 0; i < 100; ++i)
    {
      double A = axom::utilities::random_real(lower, upper);
      double B = axom::utilities::random_real(lower, upper);

      EXPECT_DOUBLE_EQ(A, axom::utilities::lerp(A, B, 0.));
      EXPECT_DOUBLE_EQ(B, axom::utilities::lerp(B, A, 0.));
      EXPECT_DOUBLE_EQ(B, axom::utilities::lerp(A, B, 1.));
      EXPECT_DOUBLE_EQ(A, axom::utilities::lerp(B, A, 1.));
    }
  }

  // Compute using different form
  {
    const double lower = -.23;  // Arbitrary end points
    const double upper = 5.73;
    for(int i = 0; i < 100; ++i)
    {
      double A = axom::utilities::random_real(lower, upper);
      double B = axom::utilities::random_real(lower, upper);

      // Test interpolation and also extrapolation beyond endpoints.
      double t = axom::utilities::random_real(-1.5, 1.5);

      double exp = A + (B - A) * t;
      EXPECT_NEAR(exp, axom::utilities::lerp(A, B, t), 1e-12);
    }
  }
}

//------------------------------------------------------------------------------
template <typename T>
axom::Array<T> make_array(int index, int repeats, int size)
{
  axom::Array<T> arr;
  arr.reserve(size + repeats - 1);
  for(int i = 0; i < size; i++)
  {
    int n = (i == index) ? repeats : 1;
    for(int j = 0; j < n; j++)
    {
      arr.push_back(i);
    }
  }
  return arr;
}

TEST(core_Utilities, binary_search)
{
  std::cout << "Testing binary_search function." << std::endl;
  const int SIZE = 100;
  const int REPEATS = 4;

  for(int size = 10; size < SIZE; size++)
  {
    for(int repeats = 1; repeats < REPEATS; repeats++)
    {
      for(int index = 0; index < size; index++)
      {
        {
          const int searchValue = index;
          const int notFoundValue = SIZE + 1;

          const auto arr = make_array<int>(index, repeats, size);
          const auto foundIndex = axom::utilities::binary_search(arr, searchValue);
          const auto notFoundIndex = axom::utilities::binary_search(arr, notFoundValue);
          EXPECT_EQ(foundIndex, index);
          EXPECT_EQ(notFoundIndex, -1);
        }

        {
          const double searchValue = index;
          const double notFoundValue = SIZE + 1;

          const auto arr = make_array<double>(index, repeats, size);
          const auto foundIndex = axom::utilities::binary_search(arr, searchValue);
          const auto notFoundIndex = axom::utilities::binary_search(arr, notFoundValue);
          EXPECT_EQ(foundIndex, index);
          EXPECT_EQ(notFoundIndex, -1);
        }
      }
    }
  }
}
