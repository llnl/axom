// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/numerics/LU.hpp"
#include "axom/core/numerics/Matrix.hpp"
#include "axom/core/numerics/matvecops.hpp"

using IndexType = axom::IndexType;

TEST(numerics_lu, lu_decompose_non_square_matrix)
{
  int pivots[3];
  const int M = 2;
  const int N = 3;

  EXPECT_TRUE(M != N);
  axom::numerics::Matrix<double> A = axom::numerics::Matrix<double>::ones(M, N);
  int rc = axom::numerics::lu_decompose(A, pivots);
  EXPECT_EQ(axom::numerics::LU_NONSQUARE_MATRIX, rc);
}

//------------------------------------------------------------------------------
TEST(numerics_lu, lu_decompose_singular_matrix)
{
  int pivots[3];
  const int N = 3;

  axom::numerics::Matrix<double> A = axom::numerics::Matrix<double>::zeros(N, N);
  int rc = axom::numerics::lu_decompose(A, pivots);
  EXPECT_EQ(axom::numerics::LU_SINGULAR_MATRIX, rc);
}

//------------------------------------------------------------------------------
TEST(numerics_lu, lu_decompose)
{
  int pivots[3];
  const int N = 3;

  // initial matrix
  axom::numerics::Matrix<double> A(N, N);
  // clang-format off
  A(0, 0) = 1; A(0, 1) = -2; A(0, 2) =  3;
  A(1, 0) = 2; A(1, 1) = -5; A(1, 2) = 12;
  A(2, 0) = 0; A(2, 1) =  2; A(2, 2) =-10;
  // clang-format on

  // decompose matrix
  axom::numerics::Matrix<double> A_decomposed = A;
  int rc = axom::numerics::lu_decompose(A_decomposed, pivots);
  EXPECT_EQ(axom::numerics::LU_SUCCESS, rc);

  // extract lower/upper triangular matrix
  axom::numerics::Matrix<double> L = axom::numerics::lower_triangular(A_decomposed);
  axom::numerics::Matrix<double> U = axom::numerics::upper_triangular(A_decomposed);

  // construct permutation matrix on lhs
  axom::numerics::Matrix<double> P = axom::numerics::Matrix<double>::identity(N);
  for(int i = 0; i < N; ++i)
  {
    if(pivots[i] != i)
    {
      P.swapRows(i, pivots[i]);
    }
  }

  axom::numerics::Matrix<double> PA(N, N);
  axom::numerics::matrix_multiply(P, A, PA);

  axom::numerics::Matrix<double> LU(N, N);
  axom::numerics::matrix_multiply(L, U, LU);

  for(int i = 0; i < N; ++i)
  {
    for(int j = 0; j < N; ++j)
    {
      EXPECT_EQ(PA(i, j), LU(i, j));
    }  // END for all colummns
  }  // END for all rows
}
//------------------------------------------------------------------------------
TEST(numerics_lu, lu_solve)
{
  int pivots[3];
  const int N = 3;

  // initial matrix
  axom::numerics::Matrix<double> A(N, N);

  // clang-format off
  A(0, 0) = 1; A(0, 1) = 2; A(0, 2) =  4;
  A(1, 0) = 3; A(1, 1) = 8; A(1, 2) = 14;
  A(2, 0) = 2; A(2, 1) = 6; A(2, 2) = 13;
  // clang-format on

  double b[3] = {3, 13, 4};
  double x[3];
  axom::numerics::Matrix<double> A_decomposed = A;
  int rc = axom::numerics::lu_decompose(A_decomposed, pivots);
  EXPECT_EQ(axom::numerics::LU_SUCCESS, rc);

  rc = axom::numerics::lu_solve(A_decomposed, pivots, b, x);
  EXPECT_EQ(axom::numerics::LU_SUCCESS, rc);
  EXPECT_DOUBLE_EQ(3, x[0]);
  EXPECT_DOUBLE_EQ(4, x[1]);
  EXPECT_DOUBLE_EQ(-2, x[2]);

  double rhs[3];
  axom::numerics::matrix_vector_multiply(A, x, rhs);

  for(IndexType i = 0; i < N; ++i)
  {
    EXPECT_DOUBLE_EQ(b[i], rhs[i]);
  }
}
