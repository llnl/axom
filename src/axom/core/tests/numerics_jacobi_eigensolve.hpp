// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/core/numerics/Matrix.hpp"             // for numerics::Matrix
#include "axom/core/utilities/Utilities.hpp"         // random_real()/isNearlyEqual
#include "axom/core/numerics/jacobi_eigensolve.hpp"  // for jacobi_eigensolve()
#include "axom/core/numerics/matvecops.hpp"          // for matrix operators

// gtest includes
#include "gtest/gtest.h"

// C/C++ includes
#include <sstream>  // for std::ostringstream

using IndexType = axom::IndexType;

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
/*!
 * \brief Generates a random symmetric matrix consisting of values within
 *  the specified interval.
 *
 * \param [in] A the matrix to fill
 * \param [in] lo lower interval bound
 * \param [in] hi upper interval bound
 *
 * \pre A.isSquare()
 * \pre lo < hi
 */
template <typename T>
void random_symmetric_matrix_init(axom::numerics::Matrix<T>& A, T lo, T hi)
{
  EXPECT_TRUE(A.isSquare());

  constexpr unsigned int seed = 123456789;
  const int N = A.getNumRows();

  for(IndexType i = 0; i < N; ++i)
  {
    for(IndexType j = 0; j <= i; ++j)
    {
      A(i, j) = axom::utilities::random_real(lo, hi, seed);
      A(j, i) = A(i, j);
    }  // END for all j
  }  // END for all i
}

//------------------------------------------------------------------------------

/*!
 * \brief Correctness check for the eigenvalues and eigenvectors
 *
 * \param [in] A the input matrix
 * \param [in] V matrix consisting of the eigenvectors column-wise
 * \param [in] lambdas corresponding vector of eigenvalues
 * \param [in] do_gtest_checks toggles gtest checks. Default is set to true.
 *
 * \return status true if correct
 *
 * \pre A.isSquare()
 * \pre V.isSquare()
 * \pre lambdas != nullptr
 */
template <typename T>
bool check_eigen_decomposition(const axom::numerics::Matrix<T>& A,
                               const axom::numerics::Matrix<T>& V,
                               const T* lambdas,
                               bool do_gtest_checks = true)
{
  // sanity checks
  EXPECT_TRUE(A.isSquare());
  EXPECT_TRUE(V.isSquare());
  EXPECT_TRUE(lambdas != nullptr);

  bool status = true;
  constexpr double TOL = 1.e-12;

  const int N = A.getNumRows();

  // compute A*V - V*lambda
  axom::numerics::Matrix<double> tmp(N, N);
  axom::numerics::Matrix<double> test(N, N);
  axom::numerics::matrix_multiply(A, V, tmp);

  for(IndexType i = 0; i < N; ++i)
  {
    for(IndexType j = 0; j < N; ++j)
    {
      test(i, j) = tmp(i, j) - V(i, j) * lambdas[j];
      status = status && axom::utilities::isNearlyEqual(test(i, j), 0.0, TOL);
      if(do_gtest_checks)
      {
        EXPECT_TRUE(status);
      }
    }
  }

  double p1norm = axom::numerics::matrix_norm(test, axom::numerics::P1_NORM);
  double inftynorm = axom::numerics::matrix_norm(test, axom::numerics::INF_NORM);
  double frobnorm = axom::numerics::matrix_norm(test, axom::numerics::FROBENIUS_NORM);

  status = status && axom::utilities::isNearlyEqual(p1norm, 0.0, TOL);
  if(do_gtest_checks)
  {
    EXPECT_TRUE(status);
  }

  status = status && axom::utilities::isNearlyEqual(inftynorm, 0.0, TOL);
  if(do_gtest_checks)
  {
    EXPECT_TRUE(status);
  }

  status = status && axom::utilities::isNearlyEqual(frobnorm, 0.0, TOL);
  if(do_gtest_checks)
  {
    EXPECT_TRUE(status);
  }

  return status;
}

}  // namespace

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(numerics_jacobi_eigensolve, check_non_convergence)
{
  constexpr int N = 40;
  constexpr int MAX_ITERS = 1;  // set to 1 so that it does not converge

  double lambdas[40];
  axom::numerics::Matrix<double> V(N, N);

  // create 40 x 40 random symmetric matrix A
  axom::numerics::Matrix<double> A(N, N);
  random_symmetric_matrix_init(A, 0.0, 1.0);

  int numIters = 0;
  int rc = axom::numerics::jacobi_eigensolve(A, V, lambdas, MAX_ITERS, &numIters);
  EXPECT_EQ(numIters, MAX_ITERS);
  EXPECT_EQ(rc, axom::numerics::JACOBI_EIGENSOLVE_FAILURE);

  bool converged = check_eigen_decomposition(A, V, lambdas, false);
  EXPECT_FALSE(converged);
}

//------------------------------------------------------------------------------
TEST(numerics_jacobi_eigensolve, random_symmetric_matrix)
{
  constexpr int MAX_SIZE = 41;

  // 40x40 test. Generate random symmetric matrix A
  for(int k = 2; k < MAX_SIZE; ++k)
  {
    const int N = k;

    std::ostringstream oss;
    oss << "checking [" << N << " x " << N << "] matrix";
    SCOPED_TRACE(oss.str());

    axom::numerics::Matrix<double> A_test(N, N);
    random_symmetric_matrix_init(A_test, 0.0, 1.0);

    double* lambdas = new double[N];
    axom::numerics::Matrix<double> V(N, N);

    int numIterations = 0;
    int rc = axom::numerics::jacobi_eigensolve(A_test,
                                               V,
                                               lambdas,
                                               axom::numerics::JACOBI_DEFAULT_MAX_ITERATIONS,
                                               &numIterations,
                                               axom::numerics::JACOBI_DEFAULT_TOLERANCE);

    EXPECT_EQ(rc, axom::numerics::JACOBI_EIGENSOLVE_SUCCESS);
    EXPECT_TRUE(numIterations > 0);
    EXPECT_TRUE(numIterations < 16);

    // compute A*V - V*lambda
    bool converged = check_eigen_decomposition(A_test, V, lambdas);
    EXPECT_TRUE(converged);

    delete[] lambdas;
    lambdas = nullptr;

  }  // END for all matrix sizes
}
