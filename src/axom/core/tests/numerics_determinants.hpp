// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/numerics/Matrix.hpp"
#include "axom/core/numerics/Determinants.hpp"

TEST(numerics_determinants, determinant_of_In)
{
  const int N = 25;

  for(int i = 2; i < N; ++i)
  {
    axom::numerics::Matrix<double> In = axom::numerics::Matrix<double>::identity(i);
    double det = axom::numerics::determinant(In);
    EXPECT_DOUBLE_EQ(1.0, det);
  }
}

//------------------------------------------------------------------------------
TEST(numerics_determinants, determinant5x5)
{
  const int N = 5;
  const double EPS = 1e-11;

  axom::numerics::Matrix<double> A(N, N);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(0, 2) = 4;
  A(0, 3) = 3;
  A(0, 4) = 0;
  A(1, 0) = 2;
  A(1, 1) = 1;
  A(1, 2) = -1;
  A(1, 3) = 1;
  A(1, 4) = 3;
  A(2, 0) = 4;
  A(2, 1) = -1;
  A(2, 2) = -2;
  A(2, 3) = 5;
  A(2, 4) = 1;
  A(3, 0) = 7;
  A(3, 1) = 3;
  A(3, 2) = 6;
  A(3, 3) = 2;
  A(3, 4) = 1;
  A(4, 0) = 1;
  A(4, 1) = 0;
  A(4, 2) = -1;
  A(4, 3) = 1;
  A(4, 4) = 1;

  double computed_det = axom::numerics::determinant(A);
  double expected_det = -34.0;
  EXPECT_NEAR(expected_det, computed_det, EPS);
}
