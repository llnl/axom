// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! \file core_numerics.cpp
 *  \brief This example code is a demonstration of the Axom Core numerics.
 */

/* This example code contains snippets used in the Core Sphinx documentation.
  * They begin and end with comments such as
  *
  * timer_start
  * timer_end
  *
  * each prepended with an underscore.
  */

// Axom includes
#include "axom/core/Macros.hpp"
#include "axom/core/numerics/eigen_solve.hpp"
#include "axom/core/numerics/jacobi_eigensolve.hpp"
#include "axom/core/numerics/linear_solve.hpp"
#include "axom/core/numerics/LU.hpp"
#include "axom/core/numerics/Matrix.hpp"
#include "axom/core/numerics/matvecops.hpp"
#include "axom/core/numerics/polynomial_solvers.hpp"
#include "axom/core/utilities/Timer.hpp"
#include "axom/core/utilities/Utilities.hpp"

#ifdef WIN32
  #include "windows.h"
void sleep(int numSeconds)
{
  int numMilliSecs = numSeconds * 1000;
  Sleep(numMilliSecs);
}
#else
  #include <unistd.h>  // for sleep()
#endif

// C/C++ includes
#include <iostream>
#include <vector>

void demoVectorOps()
{
  // _vecops_start
  namespace numerics = axom::numerics;

  // First and second vectors
  double u[] = {4., 1., 0.};
  double v[] = {1., 2., 3.};
  double w[] = {0., 0., 0.};

  std::cout << std::endl;

  std::cout << "Originally, u and v are" << std::endl
            << "u = [" << u[0] << ", " << u[1] << ", " << u[2] << "]" << std::endl
            << "v = [" << v[0] << ", " << v[1] << ", " << v[2] << "]" << std::endl;

  std::cout << std::endl;

  // Calculate dot and cross products
  double dotprod = numerics::dot_product(u, v, 3);
  numerics::cross_product(u, v, w);

  std::cout << "The dot product is " << dotprod << " and the cross product is" << std::endl
            << "[" << w[0] << ", " << w[1] << ", " << w[2] << "]" << std::endl;

  std::cout << std::endl;

  // Make u orthogonal to v, then normalize v
  numerics::make_orthogonal(u, v, 3);
  numerics::normalize(v, 3);

  std::cout << "Now orthogonal u and normalized v are" << std::endl
            << "u = [" << u[0] << ", " << u[1] << ", " << u[2] << "]" << std::endl
            << "v = [" << v[0] << ", " << v[1] << ", " << v[2] << "]" << std::endl;

  std::cout << std::endl;

  // Fill a linear space
  const int lincount = 45;
  double s[lincount];
  numerics::linspace(1., 9., s, lincount);

  // Find the real roots of a cubic equation.
  // (x + 2)(x - 1)(2x - 3) = 0 = 2x^3 - x^2 - 7x + 6 has real roots at
  // x = -2, x = 1, x = 1.5.
  double coeff[] = {6., -7., -1., 2.};
  double roots[3];
  int numRoots;
  int result = numerics::solve_cubic(coeff, roots, numRoots);

  std::cout << "Root-finding returned " << result
            << " (should be 0, success)."
               "  Found "
            << numRoots << " roots (should be 3)" << std::endl
            << "at x = " << roots[0] << ", " << roots[1] << ", " << roots[2]
            << " (should be x = -2, 1, 1.5 in arbitrary order)." << std::endl;
  // _vecops_end
}

void display_eigs(double* eigvec, double* eigval, int nrows, int i)
{
  // display eigen vectors and values
  std::cout << "Eigenvalue " << i << " = " << eigval[i] << " Eigenvector " << i << " = [";
  for(int j = 0; j < nrows; ++j)
  {
    if(j > 0)
    {
      std::cout << ", ";
    }
    std::cout << eigvec[i * nrows + j];
  }
  std::cout << "]" << std::endl;
}

void display_eigs(axom::numerics::Matrix<double>& eigvec, double* eigval, int i)
{
  double* p_eigvec = eigvec.data();
  int nrows = eigvec.getNumRows();
  display_eigs(p_eigvec, eigval, nrows, i);
}

void demoMatrix()
{
  // _matctor_start
  namespace numerics = axom::numerics;

  // Here's a 3X3 matrix of double values, initialized from an array.
  const int nrows = 3;
  const int ncols = 3;
  double val[9] = {0.6, 2.4, 1.1, 2.4, 0.6, -.1, 1.1, -.1, 0.6};
  numerics::Matrix<double> A(nrows, ncols, val, true);

  // We'll make a 3X3 identity matrix.
  // The third argument specifies the value to fill the matrix.
  numerics::Matrix<double> m(nrows, ncols, 0.);
  m.fillDiagonal(1.);
  // _matctor_end

  std::cout << std::endl;

  // _matops_start
  std::cout << "Originally, the matrix A = " << std::endl << A << std::endl;

  // Multiply, add matrices
  numerics::Matrix<double> result(nrows, ncols, 0.);
  numerics::matrix_add(A, m, result);
  std::cout << "A + identity matrix = " << std::endl << result << std::endl;
  numerics::matrix_scalar_multiply(m, 2.);
  numerics::matrix_multiply(A, m, result);
  std::cout << "A * 2*(identity matrix) = " << std::endl << result << std::endl;

  double x1[3] = {1., 2., -.5};
  double b1[3];
  std::cout << "Vector x1 = [" << x1[0] << ", " << x1[1] << ", " << x1[2] << "]" << std::endl;
  numerics::matrix_vector_multiply(A, x1, b1);
  std::cout << "A * x1 = [" << b1[0] << ", " << b1[1] << ", " << b1[2] << "]" << std::endl;

  std::cout << std::endl;

  // Calculate determinant
  std::cout << "Determinant of A = " << numerics::determinant(A) << std::endl;

  // Get lower, upper triangle.
  // By default the diagonal entries are copied from A, but you can get the
  // identity vector main diagonal entries by passing true as the second
  // argument.
  numerics::Matrix<double> ltri = lower_triangular(A);
  numerics::Matrix<double> utri = upper_triangular(A, true);
  std::cout << "A's lower triangle = " << std::endl << ltri << std::endl;
  std::cout << "A's upper triangle (with 1s in the main diagonal) = " << std::endl
            << utri << std::endl;

  // Get a column from the matrix.
  double* col1 = A.getColumn(1);
  std::cout << "A's column 1 is [" << col1[0] << ", " << col1[1] << ", " << col1[2] << "]"
            << std::endl;
  // _matops_end

  std::cout << std::endl;

  // _eigs_start
  // Solve for eigenvectors and values using the power method
  // The power method calls rand(), so we need to initialize it with srand().
  std::srand(std::time(nullptr));
  double eigvec[nrows * ncols];
  double eigval[nrows];
  int res = numerics::eigen_solve(A, nrows, eigvec, eigval);
  std::cout << "Tried to find " << nrows
            << " eigenvectors and values for"
               " matrix "
            << std::endl
            << A << std::endl
            << "and the result code was " << res << " (1 = success)." << std::endl
            << std::endl;
  if(res > 0)
  {
    for(int i = 0; i < nrows; ++i)
    {
      display_eigs(eigvec, eigval, nrows, i);
    }
  }

  std::cout << std::endl;

  // Solve for eigenvectors and values using the Jacobi method.
  numerics::Matrix<double> evecs(nrows, ncols);
  res = numerics::jacobi_eigensolve(A, evecs, eigval);
  std::cout << "Using the Jacobi method, tried to find eigenvectors and "
               "eigenvalues of matrix "
            << std::endl
            << A << std::endl
            << "and the result code was " << res << " (" << numerics::JACOBI_EIGENSOLVE_SUCCESS
            << " = success)." << std::endl
            << std::endl;
  if(res == numerics::JACOBI_EIGENSOLVE_SUCCESS)
  {
    for(int i = 0; i < nrows; ++i)
    {
      display_eigs(evecs, eigval, i);
    }
  }
  // _eigs_end

  std::cout << std::endl;

  // _solve_start
  {
    // Solve a linear system Ax = b
    numerics::Matrix<double> A(nrows, ncols);
    A(0, 0) = 1;
    A(0, 1) = 2;
    A(0, 2) = 4;
    A(1, 0) = 3;
    A(1, 1) = 8;
    A(1, 2) = 14;
    A(2, 0) = 2;
    A(2, 1) = 6;
    A(2, 2) = 13;
    double b[3] = {3., 13., 4.};
    double x[3];

    int rc = numerics::linear_solve(A, b, x);

    std::cout << "Solved for x in the linear system Ax = b, where" << std::endl
              << "A = " << std::endl
              << A << " and b = [" << b[0] << ", " << b[1] << ", " << b[2] << "]." << std::endl
              << std::endl
              << "Result code is " << rc << " (0 = success)" << std::endl;
    if(rc == 0)
    {
      std::cout << "Found x = [" << x[0] << ", " << x[1] << ", " << x[2] << "]" << std::endl;
    }
  }

  std::cout << std::endl;

  {
    // Solve a linear system Ax = b using LU decomposition and back-substitution
    numerics::Matrix<double> A(nrows, ncols);
    A(0, 0) = 1;
    A(0, 1) = 2;
    A(0, 2) = 4;
    A(1, 0) = 3;
    A(1, 1) = 8;
    A(1, 2) = 14;
    A(2, 0) = 2;
    A(2, 1) = 6;
    A(2, 2) = 13;
    double b[3] = {3., 13., 4.};
    double x[3];
    int pivots[3];

    int rc = numerics::lu_decompose(A, pivots);

    std::cout << "Decomposed to " << std::endl
              << A << " with pivots [" << pivots[0] << ", " << pivots[1] << ", " << pivots[2] << "]"
              << " with result " << rc << " (" << numerics::LU_SUCCESS << " is success)"
              << std::endl;

    rc = numerics::lu_solve(A, pivots, b, x);
    if(rc == numerics::LU_SUCCESS)
    {
      std::cout << "Found x = [" << x[0] << ", " << x[1] << ", " << x[2] << "]" << std::endl;
    }
  }
  // _solve_end
}

int main(int AXOM_UNUSED_PARAM(argc), char** AXOM_UNUSED_PARAM(argv))
{
  // _timer_start
  axom::utilities::Timer t;

  t.start();

  demoVectorOps();
  demoMatrix();

  t.stop();

  std::cout << "The tests took " << t.elapsedTimeInMilliSec() << " ms." << std::endl;
  // _timer_end

  return 0;
}
