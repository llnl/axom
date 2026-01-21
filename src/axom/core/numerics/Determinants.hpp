// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_NUMERICS_DETERMINANTS_HPP_
#define AXOM_NUMERICS_DETERMINANTS_HPP_

#include "axom/core/numerics/LU.hpp"      // for lu_decompose()
#include "axom/core/numerics/Matrix.hpp"  // for Matrix

#include "axom/core/Macros.hpp"

namespace axom
{
namespace numerics
{
// clang-format off
/// \name Matrix Operators
/// @{

/*!
 * \brief Computes a 2X2 determinant of the given matrix.
 * \accelerated
 * \param [in] a00 matrix element
 * \param [in] a01 matrix element
 * \param [in] a10 matrix element
 * \param [in] a11 matrix element
 * \return det the determinant of the 2X2 matrix
 */
template <typename real>
inline AXOM_HOST_DEVICE 
real determinant(const real& a00, const real& a01,
                 const real& a10, const real& a11)
{
  return axom::utilities::fma(a00, a11, -a01*a10);
}

/*!
 * \brief Computes the 3x3 determinant of the given matrix.
 * \param [in] a00 matrix element
 * \param [in] a01 matrix element
 * \param [in] a02 matrix element
 * \param [in] a10 matrix element
 * \param [in] a11 matrix element
 * \param [in] a12 matrix element
 * \param [in] a20 matrix element
 * \param [in] a21 matrix element
 * \param [in] a22 matrix element
 * \return det the determinant of the 3x3 matrix.
 */
template <typename real>
inline AXOM_HOST_DEVICE
real determinant(const real& a00,  const real& a01,  const real& a02,
                 const real& a10,  const real& a11,  const real& a12,
                 const real& a20,  const real& a21,  const real& a22)
{
  const real m01 = axom::utilities::fma(a00, a11, - a10 * a01);
  const real m02 = axom::utilities::fma(a00, a21, - a20 * a01);
  const real m12 = axom::utilities::fma(a10, a21, - a20 * a11);

  // det = m01 * a22 - m02 * a12 + m12 * a02;
  const real t = axom::utilities::fma(-m02, a12, m12 * a02);
  return axom::utilities::fma(m01, a22, t);
}

/*!
 * \brief Computes the 4x4 determinant of the given matrix.
 *
 * \accelerated
 * \param [in] a00 matrix element
 * \param [in] a01 matrix element
 * \param [in] a02 matrix element
 * \param [in] a03 matrix element
 * \param [in] a10 matrix element
 * \param [in] a11 matrix element
 * \param [in] a12 matrix element
 * \param [in] a13 matrix element
 * \param [in] a20 matrix element
 * \param [in] a21 matrix element
 * \param [in] a22 matrix element
 * \param [in] a23 matrix element
 * \param [in] a30 matrix element
 * \param [in] a31 matrix element
 * \param [in] a32 matrix element
 * \param [in] a33 matrix element
 * \return det the determinant of the 4x4 matrix
 */
template <typename real>
AXOM_HOST_DEVICE
inline real determinant(
  const real& a00, const real& a01, const real& a02, const real& a03,
  const real& a10, const real& a11, const real& a12, const real& a13,
  const real& a20, const real& a21, const real& a22, const real& a23,
  const real& a30, const real& a31, const real& a32, const real& a33)
{
  /// 2x2 minors
  const real m01 = axom::utilities::fma(a10, a01, - a00 * a11);
  const real m02 = axom::utilities::fma(a20, a01, - a00 * a21);
  const real m03 = axom::utilities::fma(a30, a01, - a00 * a31);
  const real m12 = axom::utilities::fma(a20, a11, - a10 * a21);
  const real m13 = axom::utilities::fma(a30, a11, - a10 * a31);
  const real m23 = axom::utilities::fma(a30, a21, - a20 * a31);

  /// 3x3 minors
  // m012 = m12*a02 - m02*a12 + m01*a22
  const real m012 = axom::utilities::fma(m01, a22, axom::utilities::fma(m12, a02,  -m02 * a12));  
  // m013 = m13*a02 - m03*a12 + m01*a32
  const real m013 = axom::utilities::fma(m01, a32, axom::utilities::fma(m13, a02,  -m03 * a12));  
  // m023 = m23*a02 - m03*a22 + m02*a32
  const real m023 = axom::utilities::fma(m02, a32, axom::utilities::fma(m23, a02,  -m03 * a22));    
  // m123 = m23*a12 - m13*a22 + m12*a32
  const real m123 = axom::utilities::fma(m12, a32, axom::utilities::fma(m23, a12,  -m13 * a22));

  /// det = m123*a03 - m023*a13 + m013*a23 - m012*a33
  // Chain to reduce rounding in the alternating sum.
  const real t0  = axom::utilities::fma(-m023, a13,  m123 * a03); // m123*a03 - m023*a13
  const real t1  = axom::utilities::fma( m013, a23,  t0);         // + m013*a23
  const real det = axom::utilities::fma(-m012, a33,  t1);         // - m012*a33
  return det;
}

/*!
 * \brief Computes the determinant of the given square matrix.
 * \param [in] A an \f$ N \times N \f$ input matrix
 * \return det the computed determinant.
 * \note if \f$ A \f$ is not square or empty, this function will return 0.0
 */
template <typename real>
real determinant(const Matrix<real>& A)
{
  real det = 0.0;

  if(!A.isSquare() || A.empty())
  {
    /* short-circuit */
    return det;
  }

  const int N = A.getNumColumns();
  if(N == 1)
  {
    det = A(0, 0);
  }
  else if(N == 2)
  {

    det = determinant(A(0,0), A(0,1),
                      A(1,0), A(1,1));

  }
  else if(N == 3)
  {

    det = determinant(A(0,0), A(0,1), A(0,2),
                      A(1,0), A(1,1), A(1,2),
                      A(2,0), A(2,1), A(2,2));

  }
  else if(N == 4)
  {

    det = determinant(A(0,0), A(0,1), A(0,2), A(0,3),
                      A(1,0), A(1,1), A(1,2), A(1,3),
                      A(2,0), A(2,1), A(2,2), A(2,3),
                      A(3,0), A(3,1), A(3,2), A(3,3));

  }
  else
  {
    Matrix<real> lu = A;
    int* pivots = new int[N];

    int rc = lu_decompose(lu, pivots);
    if(rc == LU_SUCCESS)
    {
      // count number of row interchanges
      int row_interchanges = 0;
      for(int i = 0; i < N; ++i)
      {
        row_interchanges += (pivots[i] != i) ? 1 : 0;
      }  // END for all rows

      det = ((row_interchanges & 1) == 0) ? 1.0 : -1.0;
      for(int i = 0; i < N; ++i)
      {
        det *= lu(i, i);
      }
    }

    delete[] pivots;
  }

  return (det);
}
/// @}
// clang-format on

} /* namespace numerics */
} /* namespace axom */

#endif
