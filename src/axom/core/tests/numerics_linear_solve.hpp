// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/numerics/Matrix.hpp"
#include "axom/core/numerics/linear_solve.hpp"

TEST(numerics_linear_solve, linear_solve_with_identity_matrix)
{
  const int N = 25;

  for(int i = 1; i < N; ++i)
  {
    axom::numerics::Matrix<double> A = axom::numerics::Matrix<double>::identity(i);

    // form right-hand side
    double* b = new double[i];
    for(int j = 0; j < i; ++j)
    {
      b[j] = j;
    }

    // allocate solution vector
    double* x = new double[i];

    // solve
    int rc = axom::numerics::linear_solve(A, b, x);
    EXPECT_EQ(0, rc);

    // check answer
    for(int j = 0; j < i; ++j)
    {
      EXPECT_DOUBLE_EQ(j, x[j]);
    }

    delete[] b;
    delete[] x;
  }
}

//------------------------------------------------------------------------------
TEST(numerics_linear_solve, linear_solve2x2)
{
  const int N = 2;
  axom::numerics::Matrix<double> A(N, N);

  // clang-format off
  A(0, 0) = 2; A(0, 1) = 1;
  A(1, 0) = 3; A(1, 1) = 4;
  // clang-format on

  double b[2] = {1, 14};
  double x[2];

  int rc = axom::numerics::linear_solve(A, b, x);
  EXPECT_EQ(0, rc);

  EXPECT_DOUBLE_EQ(-2, x[0]);
  EXPECT_DOUBLE_EQ(5, x[1]);
}

//------------------------------------------------------------------------------
TEST(numerics_linear_solve, linear_solve3x3)
{
  const int N = 3;
  axom::numerics::Matrix<double> A(N, N);

  // clang-format off
  A(0, 0) = 4; A(0, 1) = 5; A(0, 2) =-2;
  A(1, 0) = 7; A(1, 1) =-1; A(1, 2) = 2;
  A(2, 0) = 3; A(2, 1) = 1; A(2, 2) = 4;
  // clang-format on

  double b[3] = {-14, 42, 28};
  double x[3];

  int rc = axom::numerics::linear_solve(A, b, x);
  EXPECT_EQ(0, rc);

  EXPECT_DOUBLE_EQ(4, x[0]);
  EXPECT_DOUBLE_EQ(-4, x[1]);
  EXPECT_DOUBLE_EQ(5, x[2]);
}

//------------------------------------------------------------------------------
TEST(numerics_linear_solve, linear_solve_singular)
{
  const int N = 3;
  axom::numerics::Matrix<double> A = axom::numerics::Matrix<double>::zeros(N, N);
  double b[3] = {-14, 42, 28};
  double x[3];

  int rc = axom::numerics::linear_solve(A, b, x);
  EXPECT_NE(0, rc);
}
