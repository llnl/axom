// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/numerics/matvecops.hpp"
#include "axom/core/utilities/Utilities.hpp"

namespace
{
/*!
 * \brief Checks if the given two vectors are equal.
 *
 * \param [in] u rhs vector to compare
 * \param [in] v lhs vector to compare
 * \param [in] N the size of the vector
 */
void expect_vector_eq(const double* u, const double* v, int N)
{
  for(int i = 0; i < N; ++i)
  {
    EXPECT_DOUBLE_EQ(u[i], v[i]);
  }
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// UNIT TESTS FOR MATRIX OPS
//------------------------------------------------------------------------------
TEST(numerics_matvecops, check_invalid_matrix_addition)
{
  axom::numerics::Matrix<int> A(2, 2);
  axom::numerics::Matrix<int> B(3, 3);
  axom::numerics::Matrix<int> C(3, 3);
  axom::numerics::Matrix<int> C2(3, 3);

  EXPECT_FALSE(axom::numerics::matrix_add(A, B, C));
  EXPECT_EQ(C.getNumColumns(), 1);
  EXPECT_EQ(C.getNumRows(), 1);
  EXPECT_EQ(C(0, 0), 0);

  EXPECT_FALSE(axom::numerics::matrix_add(B, C2, A));
  EXPECT_EQ(A.getNumColumns(), 1);
  EXPECT_EQ(A.getNumRows(), 1);
  EXPECT_EQ(A(0, 0), 0);
}

//------------------------------------------------------------------------------
TEST(numerics_matvecops, matrix_add)
{
  constexpr int N = 3;
  axom::numerics::Matrix<int> I1 = axom::numerics::Matrix<int>::identity(N);
  axom::numerics::Matrix<int> I2 = axom::numerics::Matrix<int>::identity(N);

  axom::numerics::Matrix<int> A(N, N);
  EXPECT_TRUE(axom::numerics::matrix_add(I1, I2, A));

  for(int i = 0; i < N; ++i)
  {
    for(int j = 0; j < N; ++j)
    {
      int expected = (i == j) ? 2 : 0;
      EXPECT_EQ(expected, A(i, j));
    }  // END for all columns
  }  // END for all rows
}

//------------------------------------------------------------------------------
TEST(numerics_matvecops, check_invalid_matrix_subtraction)
{
  axom::numerics::Matrix<int> A(2, 2);
  axom::numerics::Matrix<int> B(3, 3);
  axom::numerics::Matrix<int> C(3, 3);
  axom::numerics::Matrix<int> C2(3, 3);

  EXPECT_FALSE(axom::numerics::matrix_subtract(A, B, C));
  EXPECT_EQ(C.getNumColumns(), 1);
  EXPECT_EQ(C.getNumRows(), 1);
  EXPECT_EQ(C(0, 0), 0);

  EXPECT_FALSE(axom::numerics::matrix_subtract(B, C2, A));
  EXPECT_EQ(A.getNumColumns(), 1);
  EXPECT_EQ(A.getNumRows(), 1);
  EXPECT_EQ(A(0, 0), 0);
}

//------------------------------------------------------------------------------
TEST(numerics_matvecops, matrix_subtract)
{
  constexpr int ZERO = 0;
  constexpr int N = 3;

  axom::numerics::Matrix<int> I1 = axom::numerics::Matrix<int>::identity(N);
  axom::numerics::Matrix<int> I2 = axom::numerics::Matrix<int>::identity(N);

  axom::numerics::Matrix<int> A(N, N);
  EXPECT_TRUE(axom::numerics::matrix_subtract(I1, I2, A));

  for(int i = 0; i < N; ++i)
  {
    for(int j = 0; j < N; ++j)
    {
      EXPECT_EQ(A(i, j), ZERO);
    }  // END for all columns
  }  // END for all rows
}

//------------------------------------------------------------------------------
TEST(numerics_matvecops, check_invalid_matrix_multiplication)
{
  axom::numerics::Matrix<int> A(2, 3);
  axom::numerics::Matrix<int> B(3, 2);
  axom::numerics::Matrix<int> C(3, 3);

  EXPECT_FALSE(axom::numerics::matrix_multiply(A, B, C));
  EXPECT_EQ(C.getNumColumns(), 1);
  EXPECT_EQ(C.getNumRows(), 1);
  EXPECT_EQ(C(0, 0), 0);

  EXPECT_FALSE(axom::numerics::matrix_multiply(A, C, B));
  EXPECT_EQ(B.getNumColumns(), 1);
  EXPECT_EQ(B.getNumRows(), 1);
  EXPECT_EQ(B(0, 0), 0);
}

//------------------------------------------------------------------------------
TEST(numerics_matvecops, matrix_matrix_multiply)
{
  // Setup test matrix A
  axom::numerics::Matrix<int> A(2, 3);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(0, 2) = 3;
  A(1, 0) = 4;
  A(1, 1) = 5;
  A(1, 2) = 6;

  // Setup test matrix B
  axom::numerics::Matrix<int> B(3, 2);
  B(0, 0) = 1;
  B(0, 1) = 2;
  B(1, 0) = 3;
  B(1, 1) = 4;
  B(2, 0) = 5;
  B(2, 1) = 6;

  axom::numerics::Matrix<int> C(2, 2);
  EXPECT_TRUE(axom::numerics::matrix_multiply(A, B, C));

  EXPECT_EQ(22, C(0, 0));
  EXPECT_EQ(28, C(0, 1));
  EXPECT_EQ(49, C(1, 0));
  EXPECT_EQ(64, C(1, 1));

  axom::numerics::Matrix<int> I2 = axom::numerics::Matrix<int>::identity(2);
  axom::numerics::Matrix<int> C2(2, 2);
  EXPECT_TRUE(axom::numerics::matrix_multiply(C, I2, C2));

  EXPECT_EQ(22, C2(0, 0));
  EXPECT_EQ(28, C2(0, 1));
  EXPECT_EQ(49, C2(1, 0));
  EXPECT_EQ(64, C2(1, 1));
}

//------------------------------------------------------------------------------
TEST(numerics_matvecops, matrix_scalar_multiply)
{
  constexpr int N = 3;
  constexpr int SCALAR = 2;

  axom::numerics::Matrix<int> A = axom::numerics::Matrix<int>::identity(N);
  axom::numerics::matrix_scalar_multiply(A, SCALAR);

  for(int i = 0; i < N; ++i)
  {
    for(int j = 0; j < N; ++j)
    {
      int expected = (i == j) ? SCALAR : 0;
      EXPECT_EQ(expected, A(i, j));
    }
  }
}

//------------------------------------------------------------------------------
TEST(numerics_matvecops, matrix_vector_multiply)
{
  constexpr int N = 2;
  axom::numerics::Matrix<int> A(N, N);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;

  int res[2];
  int v[2] = {1, 2};
  int expected[2] = {5, 11};

  axom::numerics::matrix_vector_multiply(A, v, res);

  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(expected[i], res[i]);
  }
}

//------------------------------------------------------------------------------
TEST(numerics_matvecops, check_invalid_matrix_transpose)
{
  axom::numerics::Matrix<int> A(2, 3);
  axom::numerics::Matrix<int> M(10, 10);

  EXPECT_FALSE(axom::numerics::matrix_transpose(A, M));
  EXPECT_EQ(M.getNumRows(), 1);
  EXPECT_EQ(M.getNumColumns(), 1);
  EXPECT_EQ(M(0, 0), 0);
}

//------------------------------------------------------------------------------
TEST(numerics_matvecops, matrix_transpose)
{
  // Setup matrix: (1 2 3)
  //               (4 5 6)
  axom::numerics::Matrix<int> A(2, 3);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(0, 2) = 3;
  A(1, 0) = 4;
  A(1, 1) = 5;
  A(1, 2) = 6;

  // Matrix transpose
  axom::numerics::Matrix<int> M(A.getNumColumns(), A.getNumRows());
  EXPECT_TRUE(axom::numerics::matrix_transpose(A, M));

  EXPECT_EQ(1, M(0, 0));
  EXPECT_EQ(2, M(1, 0));
  EXPECT_EQ(3, M(2, 0));

  EXPECT_EQ(4, M(0, 1));
  EXPECT_EQ(5, M(1, 1));
  EXPECT_EQ(6, M(2, 1));
}

//------------------------------------------------------------------------------
TEST(numerics_matvecops, matrix_norm)
{
  constexpr int N = 3;
  constexpr bool ZERO_COPY = true;
  constexpr double EXPECTED_P1_NORM = 19;
  constexpr double EXPECTED_INFTY_NORM = 15;
  constexpr double EXPECTED_FROBENIOUS_NORM = 14.387494569938159;

  // STEP 0: construct test matrix
  // clang-format off
  double data[9 /* 3 x 3 */] = {
    -3, 2, 0,     // column 1
     5, 6, 2,     // column 2
     7, 4, 8      // column 3
  };
  // clang-format on

  axom::numerics::Matrix<double> A(N, N, data, ZERO_COPY);

  // STEP 1: test p1-norm
  double p1norm = axom::numerics::matrix_norm(A, axom::numerics::P1_NORM);
  EXPECT_DOUBLE_EQ(p1norm, EXPECTED_P1_NORM);

  // STEP 2: test infinity-norm
  double inftynorm = axom::numerics::matrix_norm(A, axom::numerics::INF_NORM);
  EXPECT_DOUBLE_EQ(inftynorm, EXPECTED_INFTY_NORM);

  // STEP 3: test frobenius norm
  double frobnorm = axom::numerics::matrix_norm(A, axom::numerics::FROBENIUS_NORM);
  EXPECT_DOUBLE_EQ(frobnorm, EXPECTED_FROBENIOUS_NORM);
}

//------------------------------------------------------------------------------
// UNIT TESTS FOR VECTOR OPS
//------------------------------------------------------------------------------
TEST(numerics_matvecops, linspace)
{
  double vc[6];  // storage for the computed vector

  const double v1[6] = {-2, -1, 0, 1, 2, 3};
  const double v2[6] = {3, 2, 1, 0, -1, -2};
  const double v3[3] = {0.0, 0.5, 1.0};

  const double x0 = -2.0;
  const double x1 = 3.0;

  EXPECT_FALSE(axom::numerics::linspace(x0, x1, vc, 0));
  EXPECT_FALSE(axom::numerics::linspace(x0, x1, vc, 1));

  EXPECT_TRUE(axom::numerics::linspace(x0, x1, vc, 6));
  expect_vector_eq(vc, v1, 6);

  EXPECT_TRUE(axom::numerics::linspace(x1, x0, vc, 6));
  expect_vector_eq(vc, v2, 6);

  EXPECT_TRUE(axom::numerics::linspace(0.0, 1.0, vc, 3));
  expect_vector_eq(vc, v3, 3);
}

//------------------------------------------------------------------------------
TEST(numerics_matvecops, vector_cross_product)
{
  const int NDIMS = 3;

  const double e1[3] = {1.0, 0.0, 0.0};
  const double e2[3] = {0.0, 1.0, 0.0};
  const double e3[3] = {0.0, 0.0, 1.0};
  const double me3[3] = {0.0, 0.0, -1.0};

  const double u[3] = {2.0, 1.0, -1.0};
  const double v[3] = {-3.0, 4.0, 1.0};
  const double u_x_v[3] = {5.0, 1.0, 11.0};

  double w[3];

  axom::numerics::cross_product(e1, e2, w);
  expect_vector_eq(w, e3, NDIMS);

  axom::numerics::cross_product(e2, e1, w);
  expect_vector_eq(w, me3, NDIMS);

  axom::numerics::cross_product(e3, e1, w);
  expect_vector_eq(w, e2, NDIMS);

  axom::numerics::cross_product(e2, e3, w);
  expect_vector_eq(w, e1, NDIMS);

  axom::numerics::cross_product(u, v, w);
  expect_vector_eq(w, u_x_v, NDIMS);
}

//------------------------------------------------------------------------------
TEST(numerics_matvecops, vector_dot_product)
{
  const int dim = 10;
  const double EPS = 1E-4;

  double* u = new double[dim]();  // 0 initialization
  double* v = new double[dim]();  // 0 initialization

  u[0] = 1.;
  v[1] = 1.;

  EXPECT_NEAR(axom::numerics::dot_product<double>(u, v, dim), 0., EPS);

  u[1] = 1.;
  v[0] = 1.;

  EXPECT_NEAR(axom::numerics::dot_product<double>(u, v, dim), 2., EPS);

  u[5] = 10.;
  v[5] = -10.;

  EXPECT_NEAR(axom::numerics::dot_product<double>(u, v, dim), -98., EPS);

  // free up allocated memory
  delete[] u;
  delete[] v;
}

//------------------------------------------------------------------------------
TEST(numerics_matvecops, vector_make_orthogonal)
{
  const int dim = 3;
  const double ONE_THIRD = 0.3333333333;
  const double EPS = 1E-4;

  double* u = new double[dim]();  // 0 initialization
  double* v = new double[dim]();  // 0 initialization

  u[0] = 1.;
  u[1] = 1.;
  u[2] = -1.;

  v[0] = 1.;
  v[1] = 1.;
  v[2] = 1.;

  axom::numerics::make_orthogonal<double>(u, v, dim);

  EXPECT_NEAR(u[0], 1. - ONE_THIRD, EPS);
  EXPECT_NEAR(u[1], 1. - ONE_THIRD, EPS);
  EXPECT_NEAR(u[2], -1. - ONE_THIRD, EPS);

  axom::numerics::make_orthogonal(u, u, dim);

  EXPECT_NEAR(u[0], 0., EPS);
  EXPECT_NEAR(u[1], 0., EPS);
  EXPECT_NEAR(u[2], 0., EPS);

  // free up allocated memory
  delete[] u;
  delete[] v;
}

//------------------------------------------------------------------------------
TEST(numerics_matvecops, vector_orthonormalize)
{
  const int dim = 3;
  const int size = 2;
  const double EPS = 1E-4;
  const double ONE_OVER_SQRT_TWO = 0.70710678;

  double* basis = new double[dim * dim]();  // 0 initialization

  double* u = basis;
  double* v = basis + dim;

  u[0] = 1.;
  u[1] = 1.;
  v[0] = 2.;

  EXPECT_TRUE(axom::numerics::orthonormalize<double>(basis, size, dim));

  EXPECT_NEAR(u[0], ONE_OVER_SQRT_TWO, EPS);
  EXPECT_NEAR(u[1], ONE_OVER_SQRT_TWO, EPS);
  EXPECT_NEAR(u[2], 0., EPS);

  EXPECT_NEAR(v[0], ONE_OVER_SQRT_TWO, EPS);
  EXPECT_NEAR(v[1], -ONE_OVER_SQRT_TWO, EPS);
  EXPECT_NEAR(v[2], 0., EPS);

  // doesn't span
  EXPECT_FALSE(axom::numerics::orthonormalize<double>(basis, dim, dim));

  // free up allocated memory
  delete[] basis;
}

//------------------------------------------------------------------------------
TEST(numerics_matvecops, vector_normalize)
{
  const int dim = 3;
  const double EPS = 1E-4;
  const double ONE_OVER_SQRT_THREE = 0.57735026;
  const double ONE_OVER_SQRT_TEN = 0.316227766;

  double* u = new double[dim]();  // 0 initialization
  double* v = new double[dim]();  // 0 initialization

  u[0] = 1.;
  u[1] = 1.;
  u[2] = 1.;

  v[0] = 0.;
  v[1] = 3.;
  v[2] = 1.;

  EXPECT_TRUE(axom::numerics::normalize<double>(u, dim));

  EXPECT_NEAR(u[0], ONE_OVER_SQRT_THREE, EPS);
  EXPECT_NEAR(u[1], ONE_OVER_SQRT_THREE, EPS);
  EXPECT_NEAR(u[2], ONE_OVER_SQRT_THREE, EPS);

  EXPECT_TRUE(axom::numerics::normalize<double>(v, dim));

  EXPECT_NEAR(v[0], 0., EPS);
  EXPECT_NEAR(v[1], 3. * ONE_OVER_SQRT_TEN, EPS);
  EXPECT_NEAR(v[2], ONE_OVER_SQRT_TEN, EPS);

  // now make u = (0, 0, 0)
  axom::numerics::make_orthogonal(u, u, dim);

  // normalize should fail
  EXPECT_FALSE(axom::numerics::normalize<double>(u, dim));

  delete[] u;
  delete[] v;
}
