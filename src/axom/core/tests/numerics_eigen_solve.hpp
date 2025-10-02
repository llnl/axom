// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/numerics/Matrix.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/core/numerics/eigen_solve.hpp"

TEST(numerics_eigen_solve, eigen_solve_with_diagonal_matrix)
{
  const int N = 10;
  const double EPS = 1E-4;

  double* u = new double[N * N];
  double* lambdas = new double[N];

  axom::numerics::Matrix<double> A(N, N);

  // make a diagonal matrix with descending values on the diagonal
  for(int i = 0; i < N; i++)
  {
    A(i, i) = N - i;
  }

  EXPECT_EQ(1, axom::numerics::eigen_solve(A, N, u, lambdas));

  // eigenvalues should be N, N - 1, ... etc. vectors should be standard basis
  // vectors in standard order
  for(int i = 0; i < N; i++)
  {
    EXPECT_NEAR(lambdas[i], static_cast<double>(N - i), EPS);
    for(int j = 0; j < N; j++)
    {
      if(j == i)
      {
        EXPECT_NEAR(u[i * N + j] * u[i * N + j], 1., EPS) << "At entry j == i == " << j;
      }
      else
      {
        EXPECT_NEAR(u[i * N + j], 0, EPS) << "At entry j = " << j << ",  i = " << i;
      }
    }
  }

  delete[] u;
  delete[] lambdas;
}

TEST(numerics_eigen_solve, eigen_solve_with_partial_diagonal)
{
  const int N = 10;
  const double EPS = 1E-4;

  double* u = new double[N * N];
  double* lambdas = new double[N];

  axom::numerics::Matrix<double> A(N, N);

  // same thing as above test, but zero out the diagonal halfway through
  for(int i = 0; i < N / 2; i++)
  {
    A(i, i) = N - i;
  }

  EXPECT_EQ(1, axom::numerics::eigen_solve(A, N, u, lambdas));

  // now eigenvals should be N, N - 1, ... N/2 + 1, and 0's for the rest
  // vecs are same as before for first half. After that any basis for the
  // remaining subspace would be valid; we'll verify the last half of
  // eigenvecs are orthonormal
  for(int i = 0; i < N / 2; i++)
  {
    EXPECT_NEAR(lambdas[i], static_cast<double>(N - i), EPS);
    for(int j = 0; j < N; j++)
    {
      if(j == i)
      {
        EXPECT_NEAR(u[i * N + j] * u[i * N + j], 1., EPS);
      }
      else
      {
        EXPECT_NEAR(u[i * N + j], 0., EPS);
      }
    }
  }

  for(int i = N / 2; i < N; i++)
  {
    EXPECT_NEAR(lambdas[i], 0., EPS);
    for(int j = N / 2; j < N; j++)
    {
      double dot = axom::numerics::dot_product<double>(u + i * N, u + j * N, N);
      if(j == i)  // check it's unit norm
      {
        EXPECT_NEAR(dot, 1., EPS);
      }
      else  // check they're orthogonal
      {
        EXPECT_NEAR(dot, 0., EPS);
      }
    }
  }

  delete[] u;
  delete[] lambdas;
}

TEST(numerics_eigen_solve, eigen_solve_with_two_by_two)
{
  const int N = 2;
  const double EPS = 1E-4;

  double* u = new double[N * N];
  double* lambdas = new double[N];

  axom::numerics::Matrix<double> A(N, N);

  A(0, 0) = 3.;
  A(0, 1) = 1.;
  A(1, 0) = 1.;
  A(1, 1) = 3.;

  int rc = axom::numerics::eigen_solve(A, N, u, lambdas);
  EXPECT_EQ(1, rc);

  // check lambdas are correct
  EXPECT_NEAR(lambdas[0], 4., EPS);
  EXPECT_NEAR(lambdas[1], 2., EPS);

  // check it has correct eigen vecs; we're doing this weird squaring thing
  // because they could be in one of two orders; namely (1/sqrt(2), 1/sqrt(2))
  // followed by (1/sqrt(2), -1/sqrt(2)) or some rearrangement or multiplied
  // by -1
  EXPECT_NEAR(u[0] * u[0], 0.5, EPS);
  EXPECT_NEAR(u[1] * u[1], 0.5, EPS);
  EXPECT_NEAR(u[2] * u[2], 0.5, EPS);
  EXPECT_NEAR(u[3] * u[3], 0.5, EPS);

  delete[] u;
  delete[] lambdas;
}

TEST(numerics_eigen_solve, eigen_solve_with_three_by_three)
{
  const int N = 3;
  const double EPS = 1E-4;

  double* u = new double[N * N];
  double* lambdas = new double[N];

  axom::numerics::Matrix<double> A(N, N);

  A(0, 0) = 5.;
  A(0, 1) = 1.;
  A(1, 0) = 1.;
  A(1, 1) = 5.;
  A(2, 2) = 10.;

  int rc = axom::numerics::eigen_solve(A, N, u, lambdas);
  EXPECT_EQ(1, rc);

  // note we are comparing squares here, since there is some ambiguity
  // in the eigenvectors (i.e. they are known up to a sign)
  EXPECT_NEAR(u[0 + 0], 0., EPS);
  EXPECT_NEAR(u[0 + 1], 0., EPS);
  EXPECT_NEAR(u[0 + 2] * u[0 + 2], 1., EPS);

  EXPECT_NEAR(u[1 * 3 + 0] * u[1 * 3 + 0], 0.5, EPS);
  EXPECT_NEAR(u[1 * 3 + 1] * u[1 * 3 + 1], 0.5, EPS);
  EXPECT_NEAR(u[1 * 3 + 2], 0., EPS);

  EXPECT_NEAR(u[2 * 3 + 0] * u[2 * 3 + 0], 0.5, EPS);
  EXPECT_NEAR(u[2 * 3 + 1] * u[2 * 3 + 1], 0.5, EPS);
  EXPECT_NEAR(u[2 * 3 + 2], 0., EPS);

  EXPECT_NEAR(lambdas[0], 10., EPS);
  EXPECT_NEAR(lambdas[1], 6., EPS);
  EXPECT_NEAR(lambdas[2], 4., EPS);

  delete[] u;
  delete[] lambdas;
}
