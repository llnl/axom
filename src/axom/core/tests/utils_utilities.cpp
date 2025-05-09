// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/utilities/Utilities.hpp"

TEST(core_Utilities, log2)
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

TEST(core_Utilities, random_real)
{
  std::cout << "Testing random_real functions (non-deterministic)." << std::endl;

  int min = 0;
  int max = 1;
  for(int offset = 0; offset < 10; ++offset)
  {
    int cur_min = min - offset;
    int cur_max = max + offset;
    for(int i = 0; i < 100; ++i)
    {
      float f_val = axom::utilities::random_real<float>(cur_min, cur_max);
      EXPECT_GE(f_val, cur_min);
      EXPECT_LT(f_val, cur_max);

      double d_val = axom::utilities::random_real<double>(cur_min, cur_max);
      EXPECT_GE(d_val, cur_min);
      EXPECT_LT(d_val, cur_max);

      long double ld_val = axom::utilities::random_real<long double>(cur_min, cur_max);
      EXPECT_GE(ld_val, cur_min);
      EXPECT_LT(ld_val, cur_max);
    }
  }
}

TEST(core_Utilities, random_real_with_seed)
{
  std::cout << "Testing random_real functions (deterministic)." << std::endl;
  constexpr unsigned int seed = 123456789;

  constexpr double a = -5.0;
  constexpr double b = 5.0;

  const double expected_reals[5] = {-1.5112829544380526,
                                    -2.3311429024686219,
                                    -3.6335370551231403,
                                    -4.714431326610093,
                                    3.6893326916732878};

  for(int i = 0; i < 5; ++i)
  {
    const double real = axom::utilities::random_real(a, b, seed);
    EXPECT_DOUBLE_EQ(real, expected_reals[i]);
    EXPECT_GE(real, a);
    EXPECT_LT(real, b);
  }
}

TEST(core_Utilities, minmax)
{
  std::cout << "Testing min and max functions." << std::endl;

  // Test simple min, max comparisons on ints
  {
    int a = 5;
    int b = 7;

    EXPECT_EQ(a, axom::utilities::min(a, b));
    EXPECT_EQ(b, axom::utilities::max(a, b));
  }

  // Test simple min, max comparisons on doubles
  {
    double a = 5.2;
    double b = -1.7;

    EXPECT_EQ(b, axom::utilities::min(a, b));
    EXPECT_EQ(a, axom::utilities::max(a, b));
  }
}

TEST(core_Utilities, lerp)
{
  std::cout << "Testing linear interpolation (lerp) function." << std::endl;

  double f0 = 50.0;
  double f1 = 100.0;

  // Test end points
  {
    EXPECT_DOUBLE_EQ(f0, axom::utilities::lerp(f0, f1, 0.));
    EXPECT_DOUBLE_EQ(f1, axom::utilities::lerp(f1, f0, 0.));

    EXPECT_DOUBLE_EQ(f1, axom::utilities::lerp(f0, f1, 1.));
    EXPECT_DOUBLE_EQ(f0, axom::utilities::lerp(f1, f0, 1.));
  }

  // Test midpoint
  {
    double t = 0.5;
    double exp = 75.;
    EXPECT_DOUBLE_EQ(exp, axom::utilities::lerp(f0, f1, t));
  }

  // Another test
  {
    double t = 0.66;
    double exp = 83.;
    EXPECT_DOUBLE_EQ(exp, axom::utilities::lerp(f0, f1, t));
  }
}
