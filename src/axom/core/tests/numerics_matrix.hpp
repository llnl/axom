// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/numerics/Matrix.hpp"

// C/C++ includes
#include <sstream>  // for std::ostringstream

using IndexType = axom::IndexType;

//-----------------------------------------------------------------------------
// HELPER ROUTINES
//-----------------------------------------------------------------------------
namespace
{
void testConstAccess(const axom::numerics::Matrix<double>& A)
{
  const int MROWS = 10;
  const int NCOLS = 10;
  EXPECT_EQ(MROWS, A.getNumRows());
  EXPECT_EQ(NCOLS, A.getNumColumns());

  for(int irow = 0; irow < MROWS; ++irow)
  {
    for(int jcol = 0; jcol < NCOLS; ++jcol)
    {
      double expval = static_cast<double>(irow * NCOLS + jcol);
      EXPECT_EQ(expval, A(irow, jcol));
    }  // END for all columns
  }  // END for all rows
}

//------------------------------------------------------------------------------
void testCopyConstructor(axom::numerics::Matrix<double> A)
{
  const int MROWS = 10;
  const int NCOLS = 10;
  EXPECT_EQ(MROWS, A.getNumRows());
  EXPECT_EQ(NCOLS, A.getNumColumns());

  for(int irow = 0; irow < MROWS; ++irow)
  {
    for(int jcol = 0; jcol < NCOLS; ++jcol)
    {
      double expval = static_cast<double>(irow * NCOLS + jcol);
      EXPECT_EQ(expval, A(irow, jcol));
    }  // END for all columns
  }  // END for all rows
}

//------------------------------------------------------------------------------
void testExternalBufferPassByValue(axom::numerics::Matrix<int> A)
{
  const int FILL_VAL = 42;

  EXPECT_FALSE(A.usesExternalBuffer());
  A.fill(FILL_VAL);

  const int nrows = A.getNumRows();
  const int ncols = A.getNumColumns();

  for(int i = 0; i < nrows; ++i)
  {
    for(int j = 0; j < ncols; ++j)
    {
      EXPECT_EQ(FILL_VAL, A(i, j));
    }  // END for all columns
  }  // END for all rows
}

} /* end unnamed namespace */

//-----------------------------------------------------------------------------
// TEST ROUTINES
//------------------------------------------------------------------------------
TEST(numerics_matrix, basic_constructor)
{
  const int MROWS = 5;
  const int NCOLS = 10;

  axom::numerics::Matrix<double> A(MROWS, NCOLS);
  EXPECT_EQ(MROWS, A.getNumRows());
  EXPECT_EQ(NCOLS, A.getNumColumns());
  EXPECT_FALSE(A.isSquare());

  for(int irow = 0; irow < MROWS; ++irow)
  {
    for(int jcol = 0; jcol < NCOLS; ++jcol)
    {
      EXPECT_DOUBLE_EQ(0.0, A(irow, jcol));
    }
  }
}

//------------------------------------------------------------------------------
TEST(numerics_matrix, basic_constructor_with_value)
{
  const int MROWS = 5;
  const int NCOLS = 10;
  const double FILL_VAL = 2.5;

  axom::numerics::Matrix<double> A(MROWS, NCOLS, FILL_VAL);
  EXPECT_EQ(MROWS, A.getNumRows());
  EXPECT_EQ(NCOLS, A.getNumColumns());
  EXPECT_FALSE(A.isSquare());

  for(int irow = 0; irow < MROWS; ++irow)
  {
    for(int jcol = 0; jcol < NCOLS; ++jcol)
    {
      EXPECT_DOUBLE_EQ(FILL_VAL, A(irow, jcol));
    }
  }
}

//------------------------------------------------------------------------------
TEST(numerics_matrix, array_constructor)
{
  const int MROWS = 2;
  const int NCOLS = 2;
  int data[4] = {1, 2, 3, 4};

  axom::numerics::Matrix<int> A(MROWS, NCOLS, data);
  EXPECT_EQ(MROWS, A.getNumRows());
  EXPECT_EQ(NCOLS, A.getNumColumns());

  int idx = 0;
  for(int j = 0; j < NCOLS; ++j)
  {
    for(int i = 0; i < MROWS; ++i)
    {
      ++idx;
      EXPECT_EQ(idx, A(i, j));
    }
  }
}

//------------------------------------------------------------------------------
TEST(numerics_matrix, array_constructor_with_external_buffer)
{
  const int MROWS = 4;
  const int NCOLS = 2;
  int data[8] = {1, 1, 1, 1, 2, 2, 2, 2};
  // BEGIN SCOPE
  {
    axom::numerics::Matrix<int> A(MROWS, NCOLS, data, true);
    EXPECT_TRUE(A.usesExternalBuffer());

    for(int i = 0; i < NCOLS; ++i)
    {
      int* col = A.getColumn(i);
      for(int j = 0; j < MROWS; ++j)
      {
        EXPECT_EQ(i + 1, col[j]);
      }  // END for all rows
    }  // END for all cols

    A.swapColumns(0, 1);
  }
  // END SCOPE

  // ensure the data is still there when the matrix goes out-of-scope
  for(int i = 0; i < 8; ++i)
  {
    int expected = (i < 4) ? 2 : 1;
    EXPECT_EQ(expected, data[i]);
  }

  // test assignment
  // BEGIN SCOPE
  {
    axom::numerics::Matrix<int> A(MROWS, NCOLS, data, true);
    EXPECT_TRUE(A.usesExternalBuffer());
    A.swapColumns(0, 1);

    // deep copy A into B
    axom::numerics::Matrix<int> B = A;
    EXPECT_FALSE(B.usesExternalBuffer());

    const int nrows = B.getNumRows();
    const int ncols = B.getNumColumns();
    EXPECT_EQ(nrows, A.getNumRows());
    EXPECT_EQ(ncols, A.getNumColumns());

    for(int i = 0; i < nrows; ++i)
    {
      for(int j = 0; j < ncols; ++j)
      {
        EXPECT_EQ(A(i, j), B(i, j));
      }
    }

    const int FILL_VAL = 3;
    B.fill(FILL_VAL);
    for(int i = 0; i < nrows; ++i)
    {
      for(int j = 0; j < ncols; ++j)
      {
        EXPECT_EQ(FILL_VAL, B(i, j));
        EXPECT_FALSE(A(i, j) == B(i, j));
      }
    }

    for(int i = 0; i < NCOLS; ++i)
    {
      int* col = A.getColumn(i);
      for(int j = 0; j < MROWS; ++j)
      {
        EXPECT_EQ(i + 1, col[j]);
      }  // END for all rows
    }  // END for all cols
  }
  // END SCOPE

  // test copy constructor
  // BEGIN SCOPE
  {
    axom::numerics::Matrix<int> A(MROWS, NCOLS, data, true);
    EXPECT_TRUE(A.usesExternalBuffer());

    testExternalBufferPassByValue(A);

    for(int i = 0; i < NCOLS; ++i)
    {
      int* col = A.getColumn(i);
      for(int j = 0; j < MROWS; ++j)
      {
        EXPECT_EQ(i + 1, col[j]);
      }  // END for all rows
    }  // END for all cols
  }
  // END SCOPE
}

//------------------------------------------------------------------------------
TEST(numerics_matrix, is_square)
{
  const int MROWS = 10;
  const int NCOLS = 10;

  axom::numerics::Matrix<double> A(MROWS, NCOLS);
  EXPECT_TRUE(A.isSquare());
}

//------------------------------------------------------------------------------
TEST(numerics_matrix, equality)
{
  using Mtx = axom::numerics::Matrix<double>;
  constexpr double val = 2.2;

  {
    Mtx a, b;
    EXPECT_EQ(a, b);
    EXPECT_EQ(b, a);
  }

  {
    Mtx a, b(5, 7);
    EXPECT_NE(a, b);

    EXPECT_NE(b, a);
  }

  {
    Mtx a(5, 7, val), b(5, 7, val);
    EXPECT_EQ(a, b);
    EXPECT_EQ(b, a);

    a(3, 4) = 1;
    EXPECT_NE(a, b);
    EXPECT_NE(b, a);
  }

  {
    Mtx a = Mtx::identity(3);
    Mtx b = Mtx::identity(4);
    EXPECT_NE(a, b);
    EXPECT_NE(b, a);
  }

  {
    Mtx a = Mtx::identity(3);
    Mtx b = Mtx::identity(3);
    EXPECT_EQ(a, b);
    EXPECT_EQ(b, a);
  }

  {
    Mtx a = Mtx::zeros(3, 3);
    Mtx b = Mtx::ones(3, 3);
    EXPECT_NE(a, b);
    EXPECT_NE(b, a);
  }

  {
    Mtx a = Mtx::identity(3);
    Mtx b = Mtx::zeros(3, 3);
    b(0, 0) = 1.;
    b(1, 1) = 1.;
    b(2, 2) = 1.;

    EXPECT_EQ(a, b);
    EXPECT_EQ(b, a);
  }

  {
    Mtx a(7, 5);
    Mtx b(5, 7);
    EXPECT_NE(a, b);
    EXPECT_NE(b, a);

    a(3, 4) = val;
    b(4, 3) = val;
    EXPECT_NE(a, b);
    EXPECT_NE(b, a);

    Mtx a_tr(5, 7);
    axom::numerics::matrix_transpose(a, a_tr);
    EXPECT_EQ(a_tr, b);
    EXPECT_EQ(b, a_tr);
  }
}

//------------------------------------------------------------------------------
TEST(numerics_matrix, random_access_operators)
{
  const int MROWS = 10;
  const int NCOLS = 10;

  axom::numerics::Matrix<double> A(MROWS, NCOLS);

  for(int irow = 0; irow < MROWS; ++irow)
  {
    for(int jcol = 0; jcol < NCOLS; ++jcol)
    {
      double val = static_cast<double>(irow * NCOLS + jcol);
      A(irow, jcol) = val;
      EXPECT_EQ(val, A(irow, jcol));
    }  // END for all columns
  }  // END for all rows

  testConstAccess(A);
}

//------------------------------------------------------------------------------
TEST(numerics_matrix, copy_constructor)
{
  const int MROWS = 10;
  const int NCOLS = 10;

  axom::numerics::Matrix<double> A(MROWS, NCOLS);

  for(int irow = 0; irow < MROWS; ++irow)
  {
    for(int jcol = 0; jcol < NCOLS; ++jcol)
    {
      double val = static_cast<double>(irow * NCOLS + jcol);
      A(irow, jcol) = val;
      EXPECT_EQ(val, A(irow, jcol));
    }  // END for all columns
  }  // END for all rows

  testCopyConstructor(A);
}

//------------------------------------------------------------------------------
TEST(numerics_matrix, assignment)
{
  const int MROWS = 3;
  const int NCOLS = 3;

  axom::numerics::Matrix<int> A(MROWS, NCOLS);
  A.fill(3);

  axom::numerics::Matrix<int> B(2, 2);
  B.fill(1);

  B = A;

  EXPECT_EQ(MROWS, B.getNumRows());
  EXPECT_EQ(NCOLS, B.getNumColumns());

  for(int i = 0; i < MROWS; ++i)
  {
    for(int j = 0; j < NCOLS; ++j)
    {
      EXPECT_EQ(A(i, j), B(i, j));
    }
  }
}

//------------------------------------------------------------------------------
TEST(numerics_matrix, getColumn)
{
  const int N = 3;
  axom::numerics::Matrix<int> M = axom::numerics::Matrix<int>::identity(N);
  EXPECT_EQ(N, M.getNumRows());
  EXPECT_EQ(N, M.getNumColumns());

  int expected[] = {1, 0, 0, 0, 1, 0, 0, 0, 1};

  for(int j = 0; j < N; ++j)
  {
    const int* column = M.getColumn(j);
    for(int i = 0; i < N; ++i)
    {
      EXPECT_EQ(expected[j * N + i], column[i]);
    }  // END for all i
  }  // END for all j
}

//------------------------------------------------------------------------------
TEST(numerics_matrix, getRow)
{
  const int N = 3;
  axom::numerics::Matrix<int> A(N, N);

  int row_sums[] = {0, 0, 0};

  for(IndexType i = 0; i < N; ++i)
  {
    A.fillRow(i, i + 1);
    row_sums[i] = N * (i + 1);
  }

  IndexType p = 0;
  IndexType size = 0;
  for(IndexType i = 0; i < N; ++i)
  {
    int* row = A.getRow(i, p, size);
    for(IndexType j = 0; j < size; j += p)
    {
      EXPECT_EQ(i + 1, row[j]);
    }

    double sum = 0.0;
    for(IndexType j = 0; j < size; j += p)
    {
      sum += row[j];
    }

    EXPECT_EQ(row_sums[i], sum);
  }
}

//------------------------------------------------------------------------------
TEST(numerics_matrix, getDiagonal)
{
  const int N = 3;
  axom::numerics::Matrix<int> M = axom::numerics::Matrix<int>::identity(N);
  EXPECT_EQ(N, M.getNumRows());
  EXPECT_EQ(N, M.getNumColumns());
  EXPECT_EQ(N, M.getDiagonalSize());

  const int EXPECTED_SUM = N;

  int* diagonal = new int[N];
  M.getDiagonal(diagonal);

  int sum = 0;
  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(1, diagonal[i]);
    sum += diagonal[i];
  }

  EXPECT_EQ(EXPECTED_SUM, sum);
  delete[] diagonal;

  IndexType p = 0;
  IndexType size = 0;
  const int* diag = M.getDiagonal(p, size);
  EXPECT_EQ(N + 1, p);
  EXPECT_EQ(M.getDiagonalSize() * N, size);

  int sum2 = 0;
  for(IndexType i = 0; i < size; i += p)
  {
    EXPECT_EQ(1, diag[i]);
    sum2 += diag[i];
  }
  EXPECT_EQ(EXPECTED_SUM, sum2);
}

//------------------------------------------------------------------------------
TEST(numerics_matrix, fillDiagonal)
{
  const int N = 3;
  axom::numerics::Matrix<int> M = axom::numerics::Matrix<int>::identity(N);
  EXPECT_EQ(N, M.getNumRows());
  EXPECT_EQ(N, M.getNumColumns());
  EXPECT_EQ(N, M.getDiagonalSize());

  M.fillDiagonal(3);

  int* diagonal = new int[N];
  M.getDiagonal(diagonal);

  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(3, diagonal[i]);
  }

  delete[] diagonal;
}

//------------------------------------------------------------------------------
TEST(numerics_matrix, fillRow)
{
  const int FILL_VAL = 3;
  const int TARGET_ROW = 1;

  const int N = 3;
  axom::numerics::Matrix<int> M = axom::numerics::Matrix<int>::identity(N);
  EXPECT_EQ(N, M.getNumRows());
  EXPECT_EQ(N, M.getNumColumns());

  M.fillRow(TARGET_ROW, FILL_VAL);

  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(FILL_VAL, M(TARGET_ROW, i));
  }
}

//------------------------------------------------------------------------------
TEST(numerics_matrix, fillColumn)
{
  const int FILL_VAL = 3;
  const int TARGET_COL = 1;

  const int N = 3;
  axom::numerics::Matrix<int> M = axom::numerics::Matrix<int>::identity(N);
  EXPECT_EQ(N, M.getNumRows());
  EXPECT_EQ(N, M.getNumColumns());

  M.fillColumn(TARGET_COL, FILL_VAL);

  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(FILL_VAL, M(i, TARGET_COL));
  }
}

//------------------------------------------------------------------------------
TEST(numerics_matrix, fill)
{
  const int FILL_VAL = 3;

  const int N = 3;
  axom::numerics::Matrix<int> M = axom::numerics::Matrix<int>::identity(N);
  EXPECT_EQ(N, M.getNumRows());
  EXPECT_EQ(N, M.getNumColumns());

  M.fill(FILL_VAL);

  for(int i = 0; i < N; ++i)
  {
    for(int j = 0; j < N; ++j)
    {
      EXPECT_EQ(FILL_VAL, M(i, j));
    }
  }
}

//------------------------------------------------------------------------------
TEST(numerics_matrix, swapRows)
{
  const int M = 2;
  const int N = 3;
  axom::numerics::Matrix<int> A = axom::numerics::Matrix<int>::zeros(M, N);

  EXPECT_EQ(2, M);
  int FILL_VAL[2] = {3, 9};

  A.fillRow(0, FILL_VAL[0]);
  A.fillRow(1, FILL_VAL[1]);

  A.swapRows(0, 1);

  for(IndexType i = 0; i < N; ++i)
  {
    int* column = A.getColumn(i);
    for(IndexType j = 0; j < M; ++j)
    {
      EXPECT_EQ(FILL_VAL[(j + 1) % M], column[j]);
    }
  }
}

//------------------------------------------------------------------------------
TEST(numerics_matrix, swapColumns)
{
  const int M = 3;
  const int N = 4;

  // setup a test matrix
  axom::numerics::Matrix<int> A(M, N);
  for(IndexType i = 0; i < N; ++i)
  {
    A.fillColumn(i, i + 1);
  }

  int first_column = 0;
  int last_column = N - 1;
  A.swapColumns(first_column, last_column);

  int* first_column_data = A.getColumn(first_column);
  int* last_column_data = A.getColumn(last_column);
  for(IndexType i = 0; i < M; ++i)
  {
    EXPECT_EQ(last_column + 1, first_column_data[i]);
    EXPECT_EQ(first_column + 1, last_column_data[i]);
  }
}

//------------------------------------------------------------------------------
TEST(numerics_matrix, output_stream)
{
  axom::numerics::Matrix<int> A(2, 3);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(0, 2) = 3;
  A(1, 0) = 4;
  A(1, 1) = 5;
  A(1, 2) = 6;

  std::ostringstream expected_stream;
  expected_stream << "[ 1 2 3 ]\n"
                  << "[ 4 5 6 ]\n";

  std::ostringstream actual_stream;
  actual_stream << A;

  EXPECT_EQ(expected_stream.str(), actual_stream.str());
}

//------------------------------------------------------------------------------
TEST(numerics_matrix, is_identity)
{
  const int N = 3;
  axom::numerics::Matrix<int> M = axom::numerics::Matrix<int>::zeros(N, N);
  axom::numerics::Matrix<int> I = axom::numerics::Matrix<int>::identity(N);

  EXPECT_FALSE(M.isIdentity());
  EXPECT_TRUE(I.isIdentity());

  axom::numerics::Matrix<int> A(2, 3);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(0, 2) = 3;
  A(1, 0) = 4;
  A(1, 1) = 5;
  A(1, 2) = 6;
  EXPECT_FALSE(A.isIdentity());
}
