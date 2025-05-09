// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/core/numerics/polynomial_solvers.hpp"
#include "axom/core/utilities/Utilities.hpp"  // for isNearlyEqual()

// Google Test include
#include "gtest/gtest.h"

int count_mismatches(double* standard, double* test, int n, double thresh = 1.0e-8);

int count_mismatches(double* standard, double* test, int n, double thresh)
{
  int mcount = 0;

  for(int i = 0; i < n; ++i)
  {
    if(!axom::utilities::isNearlyEqual(standard[i], test[i], thresh))
    {
      ++mcount;
    }
  }

  return mcount;
}

TEST(numerics_polynomial_solvers, solve_linear)
{
  double coeff[2];
  double roots[1];
  double expected[1];
  int n = 1;
  int rc;

  // In these tests, we are solving ax + b = 0, so
  // coeff[1] = a, coeff[0] = b.
  {
    SCOPED_TRACE("Line 1 through origin.");
    coeff[0] = 0;
    coeff[1] = 1;
    roots[0] = 0;
    expected[0] = 0;
    rc = axom::numerics::solve_linear(coeff, roots, n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(count_mismatches(expected, roots, 1), 0);
  }

  {
    SCOPED_TRACE("Line 2 through origin.");
    coeff[0] = 0;
    coeff[1] = 18;
    roots[0] = 0;
    expected[0] = 0;
    rc = axom::numerics::solve_linear(coeff, roots, n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(count_mismatches(expected, roots, 1), 0);
  }

  {
    SCOPED_TRACE("Off origin 1");
    coeff[0] = -1;
    coeff[1] = 0.5;
    roots[0] = 0;
    expected[0] = 2;
    rc = axom::numerics::solve_linear(coeff, roots, n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(count_mismatches(expected, roots, 1), 0);
  }

  {
    SCOPED_TRACE("Off origin 2");
    coeff[0] = 0.5;
    coeff[1] = -1;
    roots[0] = 0;
    expected[0] = 0.5;
    rc = axom::numerics::solve_linear(coeff, roots, n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(count_mismatches(expected, roots, 1), 0);
  }

  {
    SCOPED_TRACE("X-axis");
    coeff[0] = 0;
    coeff[1] = 0;
    roots[0] = 0;
    expected[0] = 0;
    rc = axom::numerics::solve_linear(coeff, roots, n);
    // rc == 0 because there are real solutions
    EXPECT_EQ(rc, 0);
    // n == -1 because there are infinitely many solutions
    EXPECT_EQ(n, -1);
    EXPECT_EQ(count_mismatches(expected, roots, 1), 0);
  }
}

TEST(numerics_polynomial_solvers, solve_quadratic)
{
  double coeff[3];
  double roots[2];
  double expected[2];
  int n = 2;
  int rc;

  // In these tests, we are solving ax^2 + bx + c = 0, so
  // coeff[2] = a, coeff[1] = b, coeff[0] = c.
  {
    // y = (x + 2.3)(x + 2.3)
    SCOPED_TRACE("Double root at x = -2.3");
    coeff[0] = 5.29;
    coeff[1] = 4.6;
    coeff[2] = 1;
    roots[0] = 0;
    roots[1] = 0;
    expected[0] = -2.3;
    expected[1] = -2.3;
    rc = axom::numerics::solve_quadratic(coeff, roots, n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(count_mismatches(expected, roots, 2), 0);
  }

  {
    // y = (-x + 1.5)(x - 1.5)
    SCOPED_TRACE("Double root at x = 1.5 (opening down)");
    coeff[0] = -2.25;
    coeff[1] = 3;
    coeff[2] = -1;
    roots[0] = 0;
    roots[1] = 0;
    expected[0] = 1.5;
    expected[1] = 1.5;
    rc = axom::numerics::solve_quadratic(coeff, roots, n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(count_mismatches(expected, roots, 2), 0);
  }

  {
    // y = 3.2(x + 0.7)(x - 2)
    SCOPED_TRACE("Roots at -0.7 and 2");
    coeff[0] = -4.48;
    coeff[1] = -4.16;
    coeff[2] = 3.2;
    roots[0] = 0;
    roots[1] = 0;
    expected[0] = 2;
    expected[1] = -0.7;
    rc = axom::numerics::solve_quadratic(coeff, roots, n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 2);
    EXPECT_EQ(count_mismatches(expected, roots, 2), 0);
  }

  {
    // y = 0.1x^2 + 0.2x + 6
    SCOPED_TRACE("No real roots (opening up)");
    coeff[0] = 6;
    coeff[1] = 0.2;
    coeff[2] = 0.1;
    roots[0] = 0;
    roots[1] = 0;
    expected[0] = 0;
    expected[1] = 0;
    rc = axom::numerics::solve_quadratic(coeff, roots, n);
    EXPECT_EQ(rc, -1);
    EXPECT_EQ(n, 0);
    EXPECT_EQ(count_mismatches(expected, roots, 2), 0);
  }

  {
    // y = -5x^2 + 0.2x -20
    SCOPED_TRACE("No real roots (opening up)");
    coeff[0] = 6;
    coeff[1] = 0.2;
    coeff[2] = 0.1;
    roots[0] = 0;
    roots[1] = 0;
    expected[0] = 0;
    expected[1] = 0;
    rc = axom::numerics::solve_quadratic(coeff, roots, n);
    EXPECT_EQ(rc, -1);
    EXPECT_EQ(n, 0);
    EXPECT_EQ(count_mismatches(expected, roots, 2), 0);
  }
}

TEST(numerics_polynomial_solvers, solve_cubic)
{
  double coeff[4];
  double roots[3];
  double expected[3];
  int n = 3;
  int rc;

  // In these tests, we are solving ax^3 + bx^2 + cx + d = 0, so
  // coeff[3] = a, coeff[2] = b, coeff[1] = c, coeff[0] = d.
  {
    // y = (x - 1.2)^3 = 1.0 x^3 - 3.6 x^2 + 4.32 x - 1.728
    SCOPED_TRACE("Triple root at x = 1.2 (NOTE: LOOSE TOLERANCE HERE)");
    coeff[0] = -1.728;
    coeff[1] = 4.32;
    coeff[2] = -3.6;
    coeff[3] = 1;
    roots[0] = 0;
    roots[1] = 0;
    roots[2] = 0;
    expected[0] = 1.2;
    expected[1] = 1.2;
    expected[2] = 1.2;
    rc = axom::numerics::solve_cubic(coeff, roots, n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(count_mismatches(expected, roots, 3, 5.0e-5), 0);
  }

  {
    // y = -x^3 - x + 10
    SCOPED_TRACE("Single real root at x = 2");
    coeff[0] = 10;
    coeff[1] = -1;
    coeff[2] = 0;
    coeff[3] = -1;
    roots[0] = 0;
    roots[1] = 0;
    roots[2] = 0;
    expected[0] = 2;
    expected[1] = 0;
    expected[2] = 0;
    rc = axom::numerics::solve_cubic(coeff, roots, n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(count_mismatches(expected, roots, 3), 0);
  }

  {
    // y = x^3 + x^2 - 8*x - 12
    SCOPED_TRACE("Two real roots, at x = 3, twice at x = -2");
    coeff[0] = -12;
    coeff[1] = -8;
    coeff[2] = 1;
    coeff[3] = 1;
    roots[0] = 0;
    roots[1] = 0;
    roots[2] = 0;
    expected[0] = 3;
    expected[1] = -2;
    expected[2] = -2;
    rc = axom::numerics::solve_cubic(coeff, roots, n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 2);
    EXPECT_EQ(count_mismatches(expected, roots, 3), 0);
  }

  {
    // y = (x + 0.8)(-x + 1)(x - 8) = -3x^3 + 24.6x^2 - 2.4x - 19.2
    SCOPED_TRACE("Three real roots, at x = -0.8, 1, 8");
    coeff[0] = -19.2;
    coeff[1] = -2.4;
    coeff[2] = 24.6;
    coeff[3] = -3;
    roots[0] = 0;
    roots[1] = 0;
    roots[2] = 0;
    expected[0] = 8;
    expected[1] = 1;
    expected[2] = -0.8;
    rc = axom::numerics::solve_cubic(coeff, roots, n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 3);
    EXPECT_EQ(count_mismatches(expected, roots, 3), 0);
  }

  {
    // y = 4.3(x + 38)(x + 1)(x - 0.001)
    //   = 4.3x^3 + 167.6957x^2 + 163.2323x - 0.1634
    SCOPED_TRACE("Three real roots, at x = -38, -1, 0.001");
    coeff[0] = -0.1634;
    coeff[1] = 163.2323;
    coeff[2] = 167.6957;
    coeff[3] = 4.3;
    roots[0] = 0;
    roots[1] = 0;
    roots[2] = 0;
    expected[0] = 0.001;
    expected[1] = -1;
    expected[2] = -38;
    rc = axom::numerics::solve_cubic(coeff, roots, n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 3);
    EXPECT_EQ(count_mismatches(expected, roots, 3), 0);
  }
}
