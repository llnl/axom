// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MATRIX_NORMS_HPP_
#define AXOM_MATRIX_NORMS_HPP_

#include "axom/core/numerics/Matrix.hpp"      // for numerics::Matrix
#include "axom/core/utilities/Utilities.hpp"  // for utilities::abs()

namespace axom
{
namespace numerics
{
namespace internal
{
/*!
 * \brief Computes the p1-norm of a given \f$ M \times N \f$ matrix
 *
 * \param [in] A the matrix whose norm is computed
 * \return norm the computed p1-norm of the matrix
 *
 * \note The p1-norm is computed by taking the maximum absolute column sum,
 *  given by:
 *  \f[
 *     ||A||_1 = \max\limits_{1 \le j \le N}( \sum_{i=1}^M | a_{ij} | )
 *  \f]
 *
 * \pre A.getNumRows() >= 2
 * \pre A.getNumCols() >= 2
 */
template <typename T>
inline T matrix_p1_norm(const Matrix<T>& A)
{
  const int numRows = A.getNumRows();
  const int numCols = A.getNumColumns();

  if(numRows < 2 || numCols < 2)
  {
    /* short-circuit and return -1.0 to indicate error */
    return -1.0;
  }

  T norm = -1.0;
  for(IndexType j = 0; j < numCols; ++j)
  {
    const T* col = A.getColumn(j);
    T col_sum = 0.0;
    for(IndexType i = 0; i < numRows; ++i)
    {
      col_sum += utilities::abs(col[i]);
    }  // END for all rows

    norm = (col_sum > norm) ? col_sum : norm;
  }  // END for all columns

  return norm;
}

/*!
 * \brief Computes the infinity-norm of a given \f$ M \times N \f$ matrix
 *
 * \param [in] A the matrix whose norm is computed
 * \return norm the computed infinity-norm of the matrix.
 *
 * \note The infinity-norm is computed by taking the maximum absolute row sum,
 *  given by:
 *  \f[
 *    ||A||_\infty = \max\limits_{1 \le j \le M}( \sum_{i=1}^N | a_{ij} | )
 *  \f]
 *
 * \pre A.getNumRows() >= 2
 * \pre A.getNumCols() >= 2
 */
template <typename T>
inline T matrix_infty_norm(const Matrix<T>& A)
{
  const int numRows = A.getNumRows();
  const int numCols = A.getNumColumns();

  if(numRows < 2 || numCols < 2)
  {
    /* short-circuit and return -1.0 to indicate error */
    return -1.0;
  }

  T norm = -1.0;
  for(IndexType i = 0; i < numRows; ++i)
  {
    // The lines bracketed by rowsum_start and rowsum_end (prepended with
    // an underscore) are used in the Sphinx documentation of Matrix.
    // _rowsum_start
    IndexType p = 0;
    IndexType N = 0;
    const T* row = A.getRow(i, p, N);

    T row_sum = 0.0;
    for(IndexType j = 0; j < N; j += p)
    {
      row_sum += utilities::abs(row[j]);
    }  // END for all columns
    // _rowsum_end

    norm = (row_sum > norm) ? row_sum : norm;
  }  // END for all rows

  return norm;
}

/*!
 * \brief Computes the frobenius norm of a given \f$ M \times N \f$ matrix
 *
 * \param [in] A the matrix whose norm is computed
 * \return norm the computed frobenius norm of the matrix.
 *
 * \note The frobenius norm of an \f$ M \times N \f$ matrix, \f$ A \f$ is
 *  defined as follows:
 *  \f[
 *    ||A||_F = \sqrt{ \sum_{i=1}^M \sum_{j=1}^N ({a_{ij}})^2  }
 *  \f]
 *
 * \pre A.getNumRows() >= 2
 * \pre A.getNumCols() >= 2
 */
template <typename T>
inline T matrix_frobenious_norm(const Matrix<T>& A)
{
  AXOM_STATIC_ASSERT_MSG(std::is_floating_point<T>::value,
                         "T is required to be a floating type for computing the frobenious norm");

  const int numRows = A.getNumRows();
  const int numCols = A.getNumColumns();

  if(numRows < 2 || numCols < 2)
  {
    /* short-circuit and return -1.0 to indicate error */
    return -1.0;
  }

  T norm = 0.0;
  for(IndexType i = 0; i < numRows; ++i)
  {
    for(IndexType j = 0; j < numCols; ++j)
    {
      const T abs_a_ij = utilities::abs(A(i, j));
      norm += abs_a_ij * abs_a_ij;
    }  // END for all columns
  }  // END for all rows

  norm = sqrt(norm);
  return norm;
}

} /* end namespace internal */
} /* end namespace numerics */
} /* end namespace axom */

#endif /* AXOM_MATRIX_NORMS_HPP_ */
