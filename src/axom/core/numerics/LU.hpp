// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_NUMERICS_LU_HPP_
#define AXOM_NUMERICS_LU_HPP_

#include "axom/core/utilities/Utilities.hpp"  // NearlyEqual(), swap() and abs()
#include "axom/core/memory_management.hpp"    // alloc() and free()
#include "axom/core/numerics/Matrix.hpp"      // for Matrix

// C/C++ includes
#include <cstring>  // for memcpy()

namespace axom
{
namespace numerics
{
enum ReturnCodes
{
  LU_SUCCESS,
  LU_SINGULAR_MATRIX,
  LU_NONSQUARE_MATRIX,
};

/// \name Matrix Operators
/// @{

/*!
 * \brief Perform LU-decomposition on a given square matrix, \f$ A \f$ using
 *  partial pivoting.
 *
 *  This method decomposes the square matrix \f$ A \f$ into a lower-triangular
 *  matrix \f$ L \f$ and an upper-triangular matrix \f$ U \f$, such that
 *  \f$ P A = L U \f$, where \f$ P \f$ is a permutation matrix that encodes
 *  any row interchanges. The resulting \f$ L \f$ and \f$ U \f$ matrices are
 *  stored in-place in the input matrix.
 *
 * \param [in,out] A the matrix to decompose. Matrix is decomposed in-place.
 * \param [in,out] pivots buffer to store pivoting, e.g., row-interchanges.
 * \return rc return code, LU_SUCCESS if operation is successful
 *
 * \note The method returns the following error codes if an error occurs:
 * <ul>
 *   <li> LU_NONSQUARE_MATRIX when input matrix is not square </li>
 *   <li> LU_SINGULAR_MATRIX when the input matrix is singular </li>
 * </ul>
 *
 * \note The matrix is decomposed in-place.
 *
 * \pre pivots != nullptr
 * \pre pivots must be able to hold A.getNumRows() entries
 * \pre A.isSquare()==true
 */
template <typename T>
int lu_decompose(Matrix<T>& A, int* pivots);

/*!
 * \brief Solve the system \f$ Ax=b \f$ by back-substitution, where, A is an LU
 *  decomposed matrix produced via a call to lu_decompose().
 *
 * \param [in] A an LU decomposed square matrix ( output from lu_decompose() )
 * \param [in] pivots stores row-interchanges ( output from lu_decompose() )
 * \param [in] b the vector on the right-hand side.
 * \param [out] x the solution vector ( computed ).
 * \return rc return code, LU_SUCCESS if operation is successful
 *
 * \pre A.isSquare()==true
 * \pre pivots != nullptr
 * \pre b != nullptr
 * \pre x != nullptr
 */
template <typename T>
int lu_solve(const Matrix<T>& A, const int* pivots, const T* b, T* x);

/// @}

} /* end namespace numerics */
} /* end namespace axom */

//------------------------------------------------------------------------------
// implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace numerics
{
template <typename T>
int lu_decompose(Matrix<T>& LU, int* pivots)
{
  // Sanity Checks
  assert("pre: pivots buffer is NULL" && (pivots != nullptr));

  if(!LU.isSquare())
  {
    return LU_NONSQUARE_MATRIX;
  }

  const int size = LU.getNumRows();

  for(IndexType i = 0; i < size; ++i)
  {
    // descend down the ith column and find pivot
    T max_element = utilities::abs(LU(i, i));  // stores max element
    pivots[i] = i;                             // row of max element
    for(IndexType j = i + 1; j < size; ++j)
    {
      T abs_value = utilities::abs(LU(j, i));
      if(max_element < abs_value)
      {
        max_element = abs_value;
        pivots[i] = j;
      }

    }  // END for all rows

    // swap rows to place the max element on the diagonal LU(i,i)
    if(pivots[i] != i)
    {
      LU.swapRows(i, pivots[i]);
    }

    if(utilities::isNearlyEqual(LU(i, i), static_cast<T>(0)))
    {
      return LU_SINGULAR_MATRIX;
    }  // END if

    // scale upper triangular entries by the diagonal
    const T scale_factor = static_cast<T>(1) / LU(i, i);

    for(IndexType j = i + 1; j < size; ++j)
    {
      LU(i, j) *= scale_factor;
    }

    // update sub-matrix by subtracting the upper triangular part
    for(IndexType irow = i + 1; irow < size; ++irow)
    {
      for(IndexType jcol = i + 1; jcol < size; ++jcol)
      {
        LU(irow, jcol) -= LU(irow, i) * LU(i, jcol);
      }  // END for all columns
    }  // END for all rows

  }  // END for all columns

  return LU_SUCCESS;
}

//------------------------------------------------------------------------------
template <typename T>
int lu_solve(const Matrix<T>& A, const int* pivots, const T* b, T* x)
{
  // Sanity checks
  assert("pre: pivots buffer is NULL!" && (pivots != nullptr));
  assert("pre: rhs vector is NULL!" && (b != nullptr));
  assert("pre: solution vector is NULL!" && (x != nullptr));

  if(!A.isSquare())
  {
    return LU_NONSQUARE_MATRIX;
  }

  const int size = A.getNumRows();
  T* rhs = axom::allocate<T>(size);
  memcpy(rhs, b, size * sizeof(T));

  // forward-solve L part (top-to-bottom)
  for(IndexType i = 0; i < size; ++i)
  {
    // account for row-interchanges
    if(pivots[i] != i)
    {
      utilities::swap(rhs[i], rhs[pivots[i]]);
    }

    x[i] = rhs[i];

    for(IndexType j = 0; j < i; ++j)
    {
      x[i] -= A(i, j) * x[j];
    }

    x[i] /= A(i, i);
  }  // END for

  // back-substitute U part (bottom-to-top)
  for(IndexType i = size - 1; i >= 0; --i)
  {
    for(IndexType j = i + 1; j < size; ++j)
    {
      x[i] -= A(i, j) * x[j];
    }  // END for j
  }  // END for i

  axom::deallocate(rhs);
  return LU_SUCCESS;
}

} /* end namespace numerics */
} /* end namespace axom */

#endif /* AXOM_NUMERICS_LU_HPP_ */
