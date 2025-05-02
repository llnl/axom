// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_NUMERICS_TRANSFORMS_HPP_
#define AXOM_NUMERICS_TRANSFORMS_HPP_

#include "axom/config.hpp"
#include "axom/core/numerics/Matrix.hpp"

#include <cassert>
#include <cmath>

namespace axom
{
namespace numerics
{
namespace transforms
{

/*!
 * \brief Return rotation matrix about X axis.
 *
 * \param angleRad The angle to rotate in radians.
 * \param ndims The number of dimension to make for the Matrix. It needs to be 3 or 4.
 *
 * \return A Matrix containing the rotation transform.
 */
template <typename T = double>
Matrix<T> xRotation(double angleRad, int ndims = 3)
{
  assert(ndims >= 3 && ndims <= 4);
  Matrix<T> M = Matrix<T>::identity(ndims);
  const T cosA = static_cast<T>(cos(angleRad));
  const T sinA = static_cast<T>(sin(angleRad));
  M(1, 1) = cosA;
  M(2, 1) = sinA;

  M(1, 2) = -sinA;
  M(2, 2) = cosA;

  return M;
}

/*!
 * \brief Return rotation matrix about Y axis.
 *
 * \param angleRad The angle to rotate in radians.
 * \param ndims The number of dimension to make for the Matrix. It needs to be 3 or 4.
 *
 * \return A Matrix containing the rotation transform.
 */
template <typename T = double>
Matrix<T> yRotation(double angleRad, int ndims = 3)
{
  assert(ndims >= 3 && ndims <= 4);
  Matrix<T> M = Matrix<T>::identity(ndims);
  const T cosA = static_cast<T>(cos(angleRad));
  const T sinA = static_cast<T>(sin(angleRad));
  M(0, 0) = cosA;
  M(2, 0) = -sinA;

  M(0, 2) = sinA;
  M(2, 2) = cosA;

  return M;
}

/*!
 * \brief Return rotation matrix about Z axis.
 *
 * \param angleRad The angle to rotate in radians.
 * \param ndims The number of dimension to make for the Matrix. It needs to be 2, 3, or 4.
 *
 * \return A Matrix containing the rotation transform.
 */
template <typename T = double>
Matrix<T> zRotation(double angleRad, int ndims = 3)
{
  assert(ndims >= 2 && ndims <= 4);
  Matrix<T> M = Matrix<T>::identity(ndims);
  const T cosA = static_cast<T>(cos(angleRad));
  const T sinA = static_cast<T>(sin(angleRad));
  M(0, 0) = cosA;
  M(1, 0) = sinA;

  M(0, 1) = -sinA;
  M(1, 1) = cosA;

  return M;
}

/*!
 * \brief Return scaling matrix.
 *
 * \param s The scaling value.
 * \param ndims The number of dimension to make for the Matrix. It needs to be 2, 3, or 4.
 *
 * \return A Matrix containing the scaling transform.
 */
template <typename T = double>
Matrix<T> scale(T s, int ndims = 3)
{
  assert(ndims >= 2 && ndims <= 4);
  Matrix<T> M = Matrix<T>::identity(ndims);
  for(int i = 0; i < std::min(3, ndims); i++)
  {
    M(i, i) = s;
  }

  return M;
}

/*!
 * \brief Return scaling matrix.
 *
 * \param sx The scaling value in x.
 * \param sy The scaling value in y.
 * \param sz The scaling value in z.
 * \param ndims The number of dimension to make for the Matrix. It needs to be 2, 3, or 4.
 *
 * \return A Matrix containing the scaling transform.
 */
template <typename T = double>
Matrix<T> scale(T sx, T sy, T sz, int ndims = 3)
{
  assert(ndims >= 2 && ndims <= 4);
  Matrix<T> M = Matrix<T>::identity(ndims);
  M(0, 0) = sx;
  M(1, 1) = sy;
  if(ndims > 2)
  {
    M(2, 2) = sz;
  }
  return M;
}

/*!
 * \brief Return scaling matrix.
 *
 * \param sx The scaling value in x.
 * \param sy The scaling value in y.
 * \param ndims The number of dimension to make for the Matrix. It needs to be 2, 3, or 4.
 *
 * \return A Matrix containing the scaling transform.
 */
template <typename T = double>
Matrix<T> scale(T sx, T sy, int ndims = 3)
{
  return scale(sx, sy, 1., ndims);
}

/*!
 * \brief Return translation matrix.
 *
 * \param tx The translation in x.
 * \param ty The translation in y.
 *
 * \return A Matrix containing the translation transform.
 */
template <typename T = double>
Matrix<T> translate(T tx, T ty)
{
  Matrix<T> M = Matrix<T>::identity(3);
  M(0, 2) = tx;
  M(1, 2) = ty;
  return M;
}

template <typename T = double>
Matrix<T> translate(T tx, T ty, T tz)
{
  Matrix<T> M = Matrix<T>::identity(4);
  M(0, 3) = tx;
  M(1, 3) = ty;
  M(2, 3) = tz;
  return M;
}

}  // end namespace transforms
}  // end namespace numerics
}  // end namespace axom

#endif
