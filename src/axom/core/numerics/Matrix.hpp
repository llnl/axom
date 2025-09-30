// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/core/memory_management.hpp"

#include "axom/fmt.hpp"

// C/C++ includes
#include <cassert>
#include <cstring>
#include <iostream>

#ifndef AXOM_MATRIX_HPP_
  #define AXOM_MATRIX_HPP_

namespace axom
{
namespace numerics
{
// Forward Declaration
template <typename T>
class Matrix;

/// \name Overloaded Matrix Operators
/// @{

/*!
 * \brief Overloaded output stream operator. Outputs the matrix coefficients
 *  in to the given output stream.
 *
 * \param [in,out] os output stream object.
 * \param [in] A user-supplied matrix instance.
 * \return os the updated output stream object.
 */
template <typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& A);

/*!
 * \brief Checks if two matrices are component-wise equal
 *
 * \param [in] lhs first matrix
 * \param [in] rhs second matrix
 * \return status true if lhs==rhs, otherwise, false.
 */
template <typename T>
bool operator==(const Matrix<T>& lhs, const Matrix<T>& rhs);

/*!
 * \brief Checks if two matrices are not component-wise equal
 *
 * \param [in] lhs first matrix
 * \param [in] rhs second matrix
 * \return status true if lhs!=rhs, otherwise, false.
 */
template <typename T>
bool operator!=(const Matrix<T>& lhs, const Matrix<T>& rhs);
/// @}

/// \name Matrix Operators
/// @{

/*!
 * \brief Extracts the lower triangular part of a square matrix.
 *
 * \param [in] A a square matrix
 * \param [in] unit_diagonal indicates if the diagonal entries are implicitly
 *  set to ones (1s) or copied from the original matrix, A. Default is false.
 *
 * \return L a new matrix instance of the same dimensions as A consisting of
 *  the lower triangular part of A and all the rest of its entries set to
 *  zero.
 *
 * \pre A.isSquare() == true
 * \post upper-triangular of matrix L is set to zeros (0s)
 */
template <typename T>
Matrix<T> lower_triangular(const Matrix<T>& A, bool unit_diagonal = false);

/*!
 * \brief Extract the upper triangular part of a square matrix.
 *
 * \param [in] A a square matrix
 * \param [in] unit_diagonal indicates if the diagonal entries are implicitly
 *  set to ones (1s) or copied from the original matrix, A. Default is true.
 *
 * \return U a new matrix instance of the same dimensions as A consisting of
 *  the upper triangulart part of A and all the rest of its entries set to
 *  zero.
 *
 * \pre A.isSquare() == true
 * \post lower-triangular of matrix U is set to zeros (0s)
 */
template <typename T>
Matrix<T> upper_triangular(const Matrix<T>& A, bool unit_diagonal = true);

/// @}

/*!
 * \accelerated
 * \class Matrix
 *
 * \brief The Matrix class is used to represent \f$ M \times N \f$ matrices. It
 *  provides common matrix operations and allows accessing matrix elements in
 *  a more natural way, using row and column indices, regardless of the
 *  underlying flat array storage layout.
 *
 * \note The underlying storage layout is column-major.
 *
 * Basic usage example:
 * \code
 * int main( int argc, char**argv )
 * {
 *   // Allocate a 5,5 matrix
 *   Matrix< double > A(5,5);
 *
 *   const int nrows = A.getNumRows(); // nrows=5
 *   const int ncols = A.getNumColumns(); // ncols=5
 *
 *   // loop over the elements of the matrix, row by row
 *   for ( IndexType i=0; i < nrows; ++i ) {
 *      for ( IndexType j=0; j < ncols; ++j ) {
 *
 *        A( i,j ) = getValue( );
 *
 *      } // END for all columns
 *   } // END for all rows
 *   return 0;
 * }
 * \endcode
 *
 * \tparam T the underlying matrix data type, e.g., float, double, etc.
 */
template <typename T>
class Matrix
{
public:
  /*!
   * \brief Default constructor
   */
  Matrix() : m_rows(0), m_cols(0), m_data(nullptr), m_usingExternal(false) { }

  /*!
   * \brief Constructor, creates a Matrix with the given rows and columns.
   *
   * \param [in] rows the number of rows in the matrix.
   * \param [in] cols the number of columns in the matrix.
   * \param [in] val optional argument to initialize the entries of the matrix.
   *  If not supplied, the default is zero.
   *
   * \pre rows >= 1
   * \pre cols >= 1
   */
  Matrix(int rows, int cols, T val = static_cast<T>(0));

  /*!
   * \brief Array constructor, creates a Matrix with the given rows and columns
   *  and initializes its entries with the data from the supplied array.
   *
   * \param [in] rows the number of rows in the matrix
   * \param [in] cols the number of columns in the matrix
   * \param [in] data pointer to user-supplied buffer to initialize the matrix.
   * \param [in] useExternal optional flag that indicates that this matrix
   *  instance should not make a deep copy of the data. Default is false.
   *
   * \pre rows >= 1
   * \pre cols >= 1
   * \pre data != nullptr
   */
  AXOM_HOST_DEVICE Matrix(int rows, int cols, T* data, bool useExternal = false);

  /*!
   * \brief Copy constructor.
   * \param [in] m the matrix instance that is being passed.
   */
  Matrix(const Matrix<T>& m);

  /*!
   * \brief Destructor.
   */
  AXOM_HOST_DEVICE ~Matrix();

  /*!
   * \brief Check to see if the matrix is square.
   * \return status true iff this instance is a square matrix, else, false.
   */
  bool isSquare() const { return (m_rows == m_cols); };

  /*!
   * \brief Checks if the matrix is empty
   * \return status true iff the matrix is empty, else, false.
   */
  bool empty() const { return (m_rows * m_cols == 0); };

  /*!
   * \brief Checks to see if the matrix has an external buffer.
   * \return status true iff the matrix has an external buffer, else, false.
   */
  bool usesExternalBuffer() const { return m_usingExternal; };

  /*!
   * \brief Returns the number of rows in the matrix..
   * \return numRows the number of rows in the matrix.
   */
  int getNumRows() const { return m_rows; };

  /*!
   * \brief Returns the number of columns in the matrix.
   * \return numCols the number of columns in the matrix.
   */
  AXOM_HOST_DEVICE int getNumColumns() const { return m_cols; };

  /*!
   * \brief Returns the size of the diagonal.
   * \return N the size of the diagonal.
   * \note For non-square matrices, \f$ N=min(num\_rows, num\_cols) \f$
   */
  int getDiagonalSize() const { return (m_rows < m_cols) ? m_rows : m_cols; };

  /*!
   * \brief Returns the diagonal entries of this Matrix instance.
   *
   * \param [in] diagonal user-supplied buffer to store the diagonal entries.
   *
   * \pre diagonal != nullptr
   * \pre diagonal must have sufficient storage for all the diagonal entries.
   *
   * \note For non-square matrices, this method will retrieve all the entries
   *  along the main diagonal, \f$ \alpha_{ii} \in \mathcal{A} \f$, where,
   *  \f$ i \in [0,N] \f$, \f$ N=min(num\_rows, num\_cols) \f$
   */
  void getDiagonal(T* diagonal) const;

  /*!
   * \brief Assigns <em>val</em> to all entries in the diagonal.
   * \param [in] val value to assign to all diagonal entries.
   *
   * \note For non-square matrices, this method will fill all the entries
   *  along the main diagonal, \f$ \alpha_{ii} \in \mathcal{A} \f$, where,
   *  \f$ i \in [0,N] \f$, \f$ N=min(num\_rows, num\_cols) \f$
   */
  void fillDiagonal(const T& val);

  /*!
   * \brief Assigns <em>val</em> to all elements in the given matrix row.
   * \param [in] i the row index.
   * \param [in] val value to assign to all elements in the given row.
   * \pre i >= 0 && i < m_rows
   */
  void fillRow(IndexType i, const T& val);

  /*!
   * \brief Assigns <em>val</em> to all elements in the given matrix column.
   * \param [in] j the column index.
   * \param [in] val value to assign to all elements in the given column.
   * \pre j >= 0 && j < m_cols
   */
  void fillColumn(IndexType j, const T& val);

  /*!
   * \brief Assigns <em>val</em> to all elements of the matrix.
   * \param [in] val value to assign to all elements of the matrix.
   */
  void fill(const T& val);

  /*!
   * \brief Swaps the rows of this matrix instance.
   * \param [in] i index of the first row to swap
   * \param [in] j index of the second row to swap
   * \pre i >= 0 && i < m_rows
   * \pre j >= 0 && j < m_rows
   */
  void swapRows(IndexType i, IndexType j);

  /*!
   * \brief Swaps the columns of this matrix instance.
   * \param [in] i index of the first column to swap
   * \param [in] j index of the second column to swap
   * \pre i >= 0 && i < m_cols
   * \pre j >= 0 && j < m_cols
   */
  void swapColumns(IndexType i, IndexType j);

  /// \name Random Access Operators
  /// @{

  /*!
   * \brief Given an \f$ M \times N \f$ matrix, \f$ \mathcal{A} \f$, return
   *  a const reference to matrix element \f$ \alpha_{ij} \f$
   *
   * \param [in] i the row index of the matrix element,
   * \param [in] j the column index of the matrix element.
   * \return \f$ \alpha_{ij} \f$ const reference to matrix element
   *
   * \pre i >= 0 && i < m_rows
   * \pre j >= 0 && j < m_cols
   */
  AXOM_HOST_DEVICE const T& operator()(IndexType i, IndexType j) const;

  /*!
   * \brief Given an \f$ M \times N \f$ matrix, \f$ \mathcal{A} \f$, return a
   *  reference to matrix element \f$ \alpha_{ij} \f$
   *
   * \param [in] i the row index of the matrix element.
   * \param [in] j the column index of the matrix element.
   * \return \f$ \alpha_{ij} \f$ const reference to matrix element
   *
   * \pre i >= 0 && i < m_rows
   * \pre j >= 0 && j < m_cols
   */
  AXOM_HOST_DEVICE T& operator()(IndexType i, IndexType j);

  /*!
   * \brief Returns a const pointer to the  \f$ jth \f$ column of an
   *  \f$ M \times N \f$ matrix, \f$ \mathcal{A} \f$
   *
   * \param [in] j the jth column of the matrix.
   * \return \f$ \mathcal{A}_{*j} \f$ pointer to the \f$ jth \f$ column.
   * \pre j >= 0 && j < m_cols
   *
   * \note Example Usage:
   * \code
   *  ...
   *  Matrix < double > A( MROWS,NCOLS );
   *  const double* column = A.getColumn( j );
   *  for ( IndexType i=0; i < MROWS; ++i ) {
   *     std::cout << column[ i ] << " ";
   *  }
   *  ...
   * \endcode
   */
  AXOM_HOST_DEVICE const T* getColumn(IndexType j) const;

  /*!
   * \brief Returns pointer to the \f$ jth \f$ column of an \f$ M \times N \f$
   *  matrix, \f$ \mathcal{A} \f$
   *
   * \param [in] j the jth column of the matrix.
   * \return \f$ \mathcal{A}_{*j} \f$ pointer to the \f$ jth \f$ column.
   * \pre j >= 0 && j < m_cols
   *
   * \note Example Usage:
   * \code
   *  ...
   *  Matrix < double > A( MROWS,NCOLS );
   *  double* column = A.getColumn( j );
   *  for ( IndexType i=0; i < MROWS; ++i ) {
   *     column[ i ] = newval;
   *  }
   *  ...
   * \endcode
   */
  AXOM_HOST_DEVICE T* getColumn(IndexType j);

  /*!
   * \brief Returns a const pointer for strided access along the main diagonal.
   *
   * \param [out] p stride used to access elements along the main diagonal.
   * \param [out] N upper-bound used to loop over the main diagonal entries.
   * \return diag pointer along the main diagonal.
   *
   * \post p = m_rows+1
   * \post N = this->getDiaonalSize()*m_rows
   *
   * \note For non-square matrices, the entries along the main diagonal, are
   *  given by \f$ \alpha_{ii} \in \mathcal{A} \f$, where,
   *  \f$ i \in [0,N] \f$, \f$ N=min(num\_rows, num\_cols) \f$
   *
   * \note Example Usage:
   * \code
   *  ...
   *  Matrix< double > A (MROWS,NCOLS);
   *
   *  IndexType  p = 0;
   *  IndexThype N = 0;
   *  const double* diag  = A.getDiagonal( p, N );
   *
   *  for ( IndexType i=0; i < N; i+=p ) {
   *     std::cout << diag[ i ];
   *  }
   *  ...
   * \endcode
   */
  const T* getDiagonal(IndexType& p, IndexType& N) const;

  /*!
   * \brief Returns a pointer for strided access along the main diagonal.
   *
   * \param [out] p stride used to access elements along the main diagonal.
   * \param [out] N upper-bound used to loop over the main diagonal entris.
   * \return diag pointer along the main diagonal.
   *
   * \post p = m_rows+1
   * \post N = this->getDiagonalSize()*m_rows
   *
   * \note For non-square matrices, the entries along the main diagonal, are
   *  given by \f$ \alpha_{ii} \in \mathcal{A} \f$, where,
   *  \f$ i \in [0,N] \f$, \f$ N=min(num\_rows, num\_cols) \f$
   *
   * \note Example Usage:
   * \code
   *  ...
   *  Matrix< double > A (MROWS,NCOLS);
   *
   *  IndexType p = 0;
   *  IndexType N = 0;
   *  const double* diag  = A.getDiagonal( p,N );
   *
   *  for ( IndexType i=0; i < N; i+=p ) {
   *     diag[ i ] = newval;
   *  }
   *  ...
   * \endcode
   */
  T* getDiagonal(IndexType& p, IndexType& N);

  /*!
   * \brief Returns a const pointer to the \f$ ith \f$ row of an
   *  \f$ M \times N \f$ matrix, \f$ \mathcal{A} \f$
   *
   * \param [in]  i index to the \f$ ith \f$ row of the matrix
   * \param [out] p stride used to access row elements
   * \param [out] N upper-bound to loop over
   * \return \f$ \mathcal{A}_{i*} \f$ pointer to the \f$ ith \f$ row.
   *
   * \pre i >= 0 && i < m_rows
   * \post p == m_rows
   *
   * \note Example Usage:
   * \code
   *   ...
   *   Matrix< double > A( MROWS,NCOLS );
   *
   *   IndexType p = 0;
   *   IndexType N = 0;
   *   const double* row = A.getRow( irow, p, N );
   *
   *   for ( IndexType j=0; j < N; j+=p ) {
   *      std::cout << row[ j ]
   *   }
   *   ...
   * \endcode
   */
  const T* getRow(IndexType i, IndexType& p, IndexType& N) const;

  /*!
   * \brief Returns a pointer to the \f$ ith \f$ row of an  \f$ M \times N \f$
   *  matrix, \f$ \mathcal{A} \f$
   *
   * \param [in]  i index to the \f$ ith \f$ row of the matrix
   * \param [out] p stride used to access row elements
   * \param [out] N upper-bound to loop over
   * \return \f$ \mathcal{A}_{i*} \f$ pointer to the \f$ ith \f$ row.
   *
   * \pre i >= 0 && i < m_rows
   * \post p == m_rows
   *
   * \note Example Usage:
   * \code
   *   ...
   *   Matrix< double > A( MROWS,NCOLS );
   *
   *   IndexType p = 0;
   *   IndexType N = 0;
   *   double* row = A.getRow( irow, p, N );
   *
   *   for ( IndexType j=0; j < N; j+=p ) {
   *      row[ j ] = newval;
   *   }
   *   ...
   * \endcode
   *
   */
  T* getRow(IndexType i, IndexType& p, IndexType& N);

  /*!
   * \brief Returns a const pointer to the raw data.
   *
   * \return ptr pointer to the raw data.
   * \note The raw data are stored in column-major layout.
   * \post ptr != nullptr
   */
  const T* data() const;

  /*!
   * \brief Returns pointer to the raw data.
   *
   * \return ptr pointer to the raw data.
   * \note The raw data are stored in column-major layout.
   * \post ptr != nullptr
   */
  T* data();

  /// @}

  /// \name Overloaded Operators
  /// @{

  /*!
   * \brief Overloaded assignment operator.
   *
   * \param [in] rhs matrix instance on the right-hand side.
   * \return M a copy of the matrix instance in rhs.
   */
  Matrix<T>& operator=(const Matrix<T>& rhs);

  /*!
   * \brief Determines whether matrix is identity matrix.
   *
   * \return True if matrix is identity, false otherwise.
   */
  bool isIdentity() const;

  /// @}

  /// \name Static Methods
  /// @{

  /*!
   * \brief Returns an <em> identity matrix </em> \f$ \mathcal{I}_n \f$
   *
   * \param [in] n the size of the identity matrix.
   * \return M the identity matrix of size \f$ \mathcal{I}_n \f$
   *
   * \pre n >= 1
   * \post M.isSquare()==true
   * \post M.getNumRows() == M.getNumCols() == n
   */
  static Matrix<T> identity(int n);

  /*!
   * \brief Returns a <em>zero</em> matrix, \f$ \mathcal{A} \f$
   *
   * \param [in] nrows the number of rows in the matrix.
   * \param [in] ncols the number of columns in the matrix.
   * \return \f$ \mathcal{A} \f$ a zero matrix.
   *
   * \pre nrows >= 1
   * \pre ncols >= 1
   * \post \f$ \alpha_{ij}=0 \forall \alpha_{ij} \in \mathcal{A} \f$
   */
  static Matrix<T> zeros(int nrows, int ncols);

  /*!
   * \brief Returns a <em>unity</em> matrix, \f$ \mathcal{A} \f$
   *
   * \param [in] nrows the number of rows in the matrix.
   * \param [in] ncols the number of columns in the matrix.
   * \return \f$ \mathcal{A} \f$ a zero matrix.
   *
   * \pre nrows >= 1
   * \pre ncols >= 1
   * \post \f$ \alpha_{ij}=1 \forall \alpha_{ij} \in \mathcal{A} \f$
   */
  static Matrix<T> ones(int nrows, int ncols);

  /// @}

private:
  /// \name Private Helper Methods
  /// @{

  /*!
   * \brief Copies the matrix into this matrix instance.
   * \param [in] rhs matrix on the right-hand side.
   */
  void copy(const Matrix<T>& rhs);

  /*!
   * \brief Deallocates all matrix data.
   */
  AXOM_HOST_DEVICE void clear();

  /// @}

  /// \name Private Data Members
  /// @{

  int m_rows;           /*!< the number of rows in the matrix */
  int m_cols;           /*!< the number of columns in the matrix */
  T* m_data;            /*!< raw storage buffer for the matrix data */
  bool m_usingExternal; /*!< indicates if an external buffer is used */

  /// @}
};

} /* end namespace numerics */
} /* end namespace axom */

//------------------------------------------------------------------------------
// Matrix class implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace numerics
{
template <typename T>
Matrix<T>::Matrix(int rows, int cols, T val) : m_rows(rows)
                                             , m_cols(cols)
                                             , m_usingExternal(false)
{
  // sanity checks
  assert(m_rows > 0);
  assert(m_cols > 0);

  m_data = allocate<T>(m_rows * m_cols);
  this->fill(val);
}

//-----------------------------------------------------------------------------
template <typename T>
AXOM_HOST_DEVICE Matrix<T>::Matrix(int rows, int cols, T* data, bool external)
  : m_rows(rows)
  , m_cols(cols)
  , m_usingExternal(external)
{
  assert(data != nullptr);

  if(m_usingExternal)
  {
    m_data = data;
  }
  else
  {
  #if defined(AXOM_DEVICE_CODE)
    assert(false);
  #else
    const int nitems = m_rows * m_cols;
    m_data = allocate<T>(nitems);
    memcpy(m_data, data, nitems * sizeof(T));
  #endif
  }
}

//-----------------------------------------------------------------------------
template <typename T>
Matrix<T>::Matrix(const Matrix<T>& rhs)
{
  m_usingExternal = false;
  m_rows = m_cols = 0;
  m_data = nullptr;
  this->copy(rhs);
}

//-----------------------------------------------------------------------------
template <typename T>
AXOM_HOST_DEVICE Matrix<T>::~Matrix()
{
  this->clear();
}

//-----------------------------------------------------------------------------
template <typename T>
void Matrix<T>::getDiagonal(T* diagonal) const
{
  assert(diagonal != nullptr);

  const int N = this->getDiagonalSize();
  const int p = m_rows + 1;
  for(IndexType i = 0, j = 0; i < N; ++i, j += p)
  {
    diagonal[i] = m_data[j];
  }
}

//-----------------------------------------------------------------------------
template <typename T>
void Matrix<T>::fillDiagonal(const T& val)
{
  const int N = this->getDiagonalSize() * m_rows;
  const int p = m_rows + 1;
  for(IndexType i = 0; i < N; i += p)
  {
    m_data[i] = val;
  }
}

//-----------------------------------------------------------------------------
template <typename T>
void Matrix<T>::fillRow(IndexType i, const T& val)
{
  assert((i >= 0) && (i < m_rows));

  const int N = (m_cols - 1) * m_rows + i + 1;
  const int p = m_rows;
  for(IndexType j = i; j < N; j += p)
  {
    m_data[j] = val;
  }
}

//-----------------------------------------------------------------------------
template <typename T>
void Matrix<T>::fillColumn(IndexType j, const T& val)
{
  assert((j >= 0) && (j < m_cols));

  const IndexType offset = j * m_rows;
  for(IndexType i = 0; i < m_rows; ++i)
  {
    m_data[offset + i] = val;
  }
}

//-----------------------------------------------------------------------------
template <typename T>
void Matrix<T>::fill(const T& val)
{
  const int nitems = m_rows * m_cols;
  for(IndexType i = 0; i < nitems; ++i)
  {
    m_data[i] = val;
  }
}

//-----------------------------------------------------------------------------
template <typename T>
void Matrix<T>::swapRows(IndexType irow, IndexType jrow)
{
  assert((irow >= 0) && (irow < m_rows));
  assert((jrow >= 0) && (jrow < m_rows));

  if(irow == jrow)
  {
    /* short-circuit */
    return;
  }

  // convenience reference to *this
  Matrix<T>& A = *this;

  for(IndexType k = 0; k < m_cols; ++k)
  {
    utilities::swap(A(irow, k), A(jrow, k));
  }
}

//-----------------------------------------------------------------------------
template <typename T>
void Matrix<T>::swapColumns(IndexType icol, IndexType jcol)
{
  assert((icol >= 0) && (icol < m_cols));
  assert((jcol >= 0) && (jcol < m_cols));

  if(icol == jcol)
  {
    /* short-circuit */
    return;
  }

  T* icol_data = this->getColumn(icol);
  T* jcol_data = this->getColumn(jcol);
  for(IndexType i = 0; i < m_rows; ++i)
  {
    const T temp = icol_data[i];
    icol_data[i] = jcol_data[i];
    jcol_data[i] = temp;
  }
}

//-----------------------------------------------------------------------------
template <typename T>
AXOM_HOST_DEVICE const T& Matrix<T>::operator()(IndexType i, IndexType j) const
{
  assert((i >= 0) && (i < m_rows));
  assert((j >= 0) && (j < m_cols));
  return m_data[j * m_rows + i];
}

//-----------------------------------------------------------------------------
template <typename T>
AXOM_HOST_DEVICE T& Matrix<T>::operator()(IndexType i, IndexType j)
{
  assert((i >= 0) && (i < m_rows));
  assert((j >= 0) && (j < m_cols));
  return m_data[j * m_rows + i];
}

//-----------------------------------------------------------------------------
template <typename T>
AXOM_HOST_DEVICE const T* Matrix<T>::getColumn(IndexType j) const
{
  assert((j >= 0) && (j < m_cols));
  return &m_data[j * m_rows];
}

//-----------------------------------------------------------------------------
template <typename T>
AXOM_HOST_DEVICE T* Matrix<T>::getColumn(IndexType j)
{
  assert((j >= 0) && (j < m_cols));
  return &m_data[j * m_rows];
}

//-----------------------------------------------------------------------------
template <typename T>
const T* Matrix<T>::getRow(IndexType i, IndexType& p, IndexType& N) const
{
  assert((i >= 0) && (i < m_rows));
  p = m_rows;
  N = (m_cols - 1) * m_rows + i + 1;
  return &m_data[i];
}

//-----------------------------------------------------------------------------
template <typename T>
T* Matrix<T>::getRow(IndexType i, IndexType& p, IndexType& N)
{
  assert((i >= 0) && (i < m_rows));
  p = m_rows;
  N = (m_cols - 1) * m_rows + i + 1;
  return &m_data[i];
}

//-----------------------------------------------------------------------------
template <typename T>
const T* Matrix<T>::getDiagonal(IndexType& p, IndexType& N) const
{
  p = m_rows + 1;
  N = this->getDiagonalSize() * m_rows;
  return &m_data[0];
}

//-----------------------------------------------------------------------------
template <typename T>
T* Matrix<T>::getDiagonal(IndexType& p, IndexType& N)
{
  p = m_rows + 1;
  N = this->getDiagonalSize() * m_rows;
  return &m_data[0];
}

//-----------------------------------------------------------------------------
template <typename T>
const T* Matrix<T>::data() const
{
  return m_data;
}

//-----------------------------------------------------------------------------
template <typename T>
T* Matrix<T>::data()
{
  return m_data;
}

//-----------------------------------------------------------------------------
template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& rhs)
{
  if(this != &rhs)
  {
    if((m_rows != rhs.m_rows) || (m_cols != rhs.m_cols))
    {
      this->clear();
    }

    this->copy(rhs);
  }

  return *this;
}

//-----------------------------------------------------------------------------
template <typename T>
bool Matrix<T>::isIdentity() const
{
  bool ok = (m_rows == m_cols);
  if(ok)
  {
    for(IndexType i = 0; i < m_rows && ok; ++i)
    {
      for(IndexType j = 0; j < m_cols && ok; ++j)
      {
        T expected = (i == j) ? 1.0 : 0.0;
        ok = (this->operator()(i, j) == expected);
      }  // END for all columns
    }  // END for all rows
  }
  return ok;
}

//-----------------------------------------------------------------------------
// STATIC METHODS
//-----------------------------------------------------------------------------
template <typename T>
Matrix<T> Matrix<T>::identity(int n)
{
  Matrix<T> In(n, n);

  for(IndexType i = 0; i < n; ++i)
  {
    for(IndexType j = 0; j < n; ++j)
    {
      In(i, j) = static_cast<T>((i == j) ? 1.0 : 0.0);
    }  // END for all columns
  }  // END for all rows

  return (In);
}

//-----------------------------------------------------------------------------
template <typename T>
Matrix<T> Matrix<T>::zeros(int nrows, int ncols)
{
  Matrix<T> M(nrows, ncols);
  M.fill(static_cast<T>(0));
  return (M);
}

//-----------------------------------------------------------------------------
template <typename T>
Matrix<T> Matrix<T>::ones(int nrows, int ncols)
{
  Matrix<T> M(nrows, ncols);
  M.fill(static_cast<T>(1));
  return (M);
}

//-----------------------------------------------------------------------------
// PRIVATE HELPER METHODS
//-----------------------------------------------------------------------------
template <typename T>
void Matrix<T>::copy(const Matrix<T>& rhs)
{
  bool do_allocate = m_usingExternal || (m_rows != rhs.m_rows) || (m_cols != rhs.m_cols);

  if(do_allocate)
  {
    assert(m_data == nullptr);

    m_rows = rhs.m_rows;
    m_cols = rhs.m_cols;
    m_data = allocate<T>(m_rows * m_cols);
  }

  assert(m_rows == rhs.m_rows);
  assert(m_cols == rhs.m_cols);
  assert(m_data != nullptr);

  const int nitems = m_rows * m_cols;
  const int bytesize = nitems * sizeof(T);
  memcpy(m_data, rhs.m_data, bytesize);
}

//-----------------------------------------------------------------------------
template <typename T>
AXOM_HOST_DEVICE void Matrix<T>::clear()
{
  #if defined(AXOM_DEVICE_CODE)
  assert(m_usingExternal);
  #else
  if(!m_usingExternal)
  {
    deallocate(m_data);
  }
  #endif

  m_rows = m_cols = 0;
}

//-----------------------------------------------------------------------------
// FREE METHODS
//-----------------------------------------------------------------------------

template <typename T>
bool operator==(const Matrix<T>& lhs, const Matrix<T>& rhs)
{
  const int R = lhs.getNumRows();
  const int C = lhs.getNumColumns();

  if(R != rhs.getNumRows() || C != rhs.getNumColumns())
  {
    return false;
  }

  using axom::utilities::isNearlyEqual;
  constexpr double EPS = 1e-12;

  for(int i = 0; i < R * C; ++i)
  {
    if(!isNearlyEqual(lhs.data()[i], rhs.data()[i], EPS))
    {
      return false;
    }
  }

  return true;
}

template <typename T>
bool operator!=(const Matrix<T>& lhs, const Matrix<T>& rhs)
{
  return !(lhs == rhs);
}

//-----------------------------------------------------------------------------
template <typename T>
Matrix<T> lower_triangular(const Matrix<T>& A, bool unit_diagonal)
{
  assert(A.isSquare());

  if(!A.isSquare())
  {
    return Matrix<T>::zeros(1, 1);
  }

  const int N = A.getNumRows();

  Matrix<T> L(N, N);

  for(IndexType i = 0; i < N; ++i)
  {
    for(IndexType j = 0; j < N; ++j)
    {
      if(i == j)
      {
        L(i, j) = (unit_diagonal) ? static_cast<T>(1) : A(i, j);
      }
      else
      {
        L(i, j) = (i > j) ? A(i, j) : static_cast<T>(0);
      }

    }  // END for all columns
  }  // END for all rows

  return (L);
}

//-----------------------------------------------------------------------------
template <typename T>
Matrix<T> upper_triangular(const Matrix<T>& A, bool unit_diagonal)
{
  assert(A.isSquare());

  if(!A.isSquare())
  {
    return Matrix<T>::zeros(1, 1);
  }

  const int N = A.getNumRows();

  Matrix<T> U(N, N);

  for(IndexType i = 0; i < N; ++i)
  {
    for(IndexType j = 0; j < N; ++j)
    {
      if(i == j)
      {
        U(i, j) = (unit_diagonal) ? static_cast<T>(1) : A(i, j);
      }
      else
      {
        U(i, j) = (i < j) ? A(i, j) : static_cast<T>(0);
      }

    }  // END for all columns
  }  // END for all rows

  return (U);
}

//-----------------------------------------------------------------------------
template <typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& M)
{
  const int nrows = M.getNumRows();
  const int ncols = M.getNumColumns();

  for(IndexType i = 0; i < nrows; ++i)
  {
    os << "[ ";
    for(int j = 0; j < ncols; ++j)
    {
      os << M(i, j) << " ";
    }  // END for all j

    os << "]\n";

  }  // END for all i

  return (os);
}

} /* end namespace numerics */
} /* end namespace axom */

/// Overload to format a numerics::Matrix using fmt
template <typename T>
struct axom::fmt::formatter<axom::numerics::Matrix<T>> : ostream_formatter
{ };

#endif /* AXOM_MATRIX_HPP_ */
