// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_COORDINATE_TRANSFORMER_HPP
#define AXOM_PRIMAL_COORDINATE_TRANSFORMER_HPP

#include "axom/core/numerics/Matrix.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/core/numerics/floating_point_limits.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"

namespace axom
{
namespace primal
{
namespace experimental
{
/*!
 * @brief 3D Coordinate transformation facilitating the placement of
 * geometries whose parameters can't easily describe it.
 *
 * The transformations may be described as translations, rotations
 * and arbitrary operators transforms on homogeneous coordinates.
 * These matrices should be 4x4 and have the last row values
 * [0, 0, 0, 1].
 *
 * To efficiently allow for large numbers of CoordinateTransformer
 * objects, this class is a POD.  Hence we don't use axom::Matrix
 * to store the matrix.
 *
 * Only supporting 3D coordinates presently.
 * This class is a new utility.  It is subject to change.
*/
template <typename T = double>
class CoordinateTransformer
{
public:
  /*!
   * @brief Default constructor sets an identity transformer.
   */
  CoordinateTransformer()
    : m_P {Vectr {1., 0., 0.}, Vectr {0., 1., 0.}, Vectr {0., 0., 1.}}
    , m_v {0., 0., 0.}
  { }

  /*!
   * @brief Copy constructor.
   */
  AXOM_HOST_DEVICE CoordinateTransformer(const CoordinateTransformer& other) { copyIn(other); }

  CoordinateTransformer& operator=(const CoordinateTransformer& other)
  {
    copyIn(other);
    return *this;
  }

  /*!
   * @brief Constructor sets the 4x4 transformation matrix.
   * @param matrix [in] The transformation matrix for homogeneous
   * coordinates.
   *
   * The last row of \c matrix must be be \c [0,0,0,1].
   */
  CoordinateTransformer(const numerics::Matrix<T>& matrix) { setMatrix(matrix); }

  /*!
   * @brief Construct transformer that moves 4 starting points to 4
   * destination points.
   * @see setByTerminusPts.
   */
  AXOM_HOST_DEVICE CoordinateTransformer(const primal::Point<T, 3>* startPts,
                                         const primal::Point<T, 3>* destPts)
  {
    setByTerminusPts(startPts, destPts);
  }

  /*!
   * @brief Set the matrix, discarding the current transformation.
   * @param matrix [in] The transformation matrix for homogeneous
   * coordinates.
   */
  void setMatrix(const numerics::Matrix<T>& matrix)
  {
    // Assert that matrix is a transformation in homogeneous coordinates.
    SLIC_ASSERT(matrix.getNumRows() == 4);
    SLIC_ASSERT(matrix.getNumColumns() == 4);
    SLIC_ASSERT(matrix(3, 0) == 0.0);
    SLIC_ASSERT(matrix(3, 1) == 0.0);
    SLIC_ASSERT(matrix(3, 2) == 0.0);
    SLIC_ASSERT(matrix(3, 3) == 1.0);
    for(int c = 0; c < 3; ++c)
    {
      for(int r = 0; r < 3; ++r)
      {
        m_P[r][c] = matrix(r, c);
      }
    }
    for(int r = 0; r < 3; ++r)
    {
      m_v[r] = matrix(r, 3);
    }
  }

  /*!
   * @brief Compute the 4x4 matrix to transform
   * 4 points to 4 specified destinations.
   * @param [in] startPts Four starting points.
   * @param [in] destPts Four destination points.
   *
   * The four starting points must define a non-degenerate volume.
   * Else the transformation is ill-defined and the transformer is set
   * to invalid.  If you think of each 4 points as the vertices of a
   * tetrahedron, the transformation moves the starting tet to the
   * destination tet.
   */
  AXOM_HOST_DEVICE void setByTerminusPts(const primal::Point<T, 3>* startPts,
                                         const primal::Point<T, 3>* destPts)
  {
    // Compute last 3 points relative to first points.
    Vectr sPts[3];
    Vectr dPts[3];
    for(int p = 0; p < 3; ++p)
    {
      sPts[p] = startPts[p + 1].array() - startPts[0].array();
      dPts[p] = destPts[p + 1].array() - destPts[0].array();
    }

    Matrx pMat;
    Matrx qMat;
    for(int r = 0; r < 3; ++r)
    {
      for(int c = 0; c < 3; ++c)
      {
        pMat[r][c] = sPts[c][r];
        qMat[r][c] = dPts[c][r];
      }
    }

    invertMatrx(pMat);
    if(std::isnan(pMat[0][0]))
    {
      setInvalid();
    }
    else
    {
      multMatrx(qMat, pMat, m_P);
      multMatrxVectr(m_P, startPts[0].array(), m_v);
      m_v = destPts[0].array() - m_v;
    }
  }

  /*!
   * @brief Set to invalid value.
   *
   * Validity can be checked with isValid().
   */
  AXOM_HOST_DEVICE void setInvalid() { m_P[0][0] = std::numeric_limits<T>::quiet_NaN(); }

  //! @brief Whether transformer is valid.
  AXOM_HOST_DEVICE bool isValid() { return !std::isnan(m_P[0][0]); }

  /*!
   * @brief Get the matrix for the transformation.
   */
  numerics::Matrix<T> getMatrix()
  {
    numerics::Matrix<T> rval(4, 4, 0.0);
    for(int c = 0; c < 3; ++c)
    {
      for(int r = 0; r < 3; ++r)
      {
        rval(r, c) = m_P[r][c];
      }
    }
    for(int r = 0; r < 3; ++r)
    {
      rval(r, 3) = m_v[r];
    }
    rval(3, 3) = 1.0;
    return rval;
  }

  /*!
   * @brief Apply a matrix transform to the current transformation.
   * @param matrix [in] The transformation matrix for homogeneous
   * coordinates.
   */
  void applyMatrix(const numerics::Matrix<T>& matrix)
  {
    numerics::Matrix<T> current = getMatrix();
    numerics::Matrix<T> updated(4, 4);
    axom::numerics::matrix_multiply(matrix, current, updated);
    setMatrix(updated);
  }

  //! @brief Apply a 3D translation to the current transformation.
  void applyTranslation(const axom::primal::Vector<T, 3>& d) { applyTranslation(d.array()); }

  //! @brief Add a 3D translation to the current transformation.
  void applyTranslation(const axom::NumericArray<T, 3>& d)
  {
    m_v[0] += d[0];
    m_v[1] += d[1];
    m_v[2] += d[2];
  }

  /*!
   * @brief Apply a 3D rotation to the current transformation.
   *
   * The rotation is defined by 2 vectors, start and end, and is not
   * unique.  The chosen rotation axis is the cross product of the
   * start and end vectors.
   *
   * @param start [in] Starting direction
   * @param end [in] Ending direction
   */
  void applyRotation(const axom::primal::Vector<T, 3>& start, const axom::primal::Vector<T, 3>& end)
  {
    Vector<T, 3> s = start.unitVector();
    Vector<T, 3> e = end.unitVector();
    Vector<T, 3> u;  // Rotation vector, the cross product of start and end.
    numerics::cross_product(s.data(), e.data(), u.data());
    const T sinT = u.norm();
    const T cosT = numerics::dot_product(s.data(), e.data(), 3);

    if(utilities::isNearlyEqual(sinT, 0.0))
    {
      // Degenerate: end is parallel to start.
      // angle near 0 (identity transform) or pi.
      if(cosT < 0)
      {
        setInvalid();  // Transformation is ill-defined
      }
      return;
    }

    u.array() /= sinT;  // Make u a unit vector.
    applyRotation(u, sinT, cosT);
  }

  /*!
   * @brief Add a 3D rotation to the current transformation.  The
   * rotation is given as a rotation axis (unit vector) and an angle.
   *
   * @param u [in] Rotation axis, a unit vector
   * @param angle [in] Rotation angle
   */
  void applyRotation(const axom::primal::Vector<T, 3>& u, T angle)
  {
    T sinT = sin(angle);
    T cosT = cos(angle);
    applyRotation(u, sinT, cosT);
  }

  /*!
   * @brief Add a 3D rotation to the current transformation.  The
   * rotation is given as a rotation axis (unit vector) and angle as
   * sine and cosine.
   *
   * @param u [in] Rotation axis, a unit vector
   * @param sinT [in] Sine of rotation angle
   * @param cosT [in] Cosine of rotation angle
   */
  void applyRotation(const axom::primal::Vector<T, 3>& u, T sinT, T cosT)
  {
    SLIC_ASSERT(
      axom::utilities::isNearlyEqual(u.squared_norm(),
                                     1.0,
                                     10 * axom::numerics::floating_point_limits<T>::epsilon()));
    const T ccosT = 1 - cosT;

    Matrx P;  // 3D rotation matrix.
    P[0][0] = u[0] * u[0] * ccosT + cosT;
    P[0][1] = u[0] * u[1] * ccosT - u[2] * sinT;
    P[0][2] = u[0] * u[2] * ccosT + u[1] * sinT;
    P[1][0] = u[1] * u[0] * ccosT + u[2] * sinT;
    P[1][1] = u[1] * u[1] * ccosT + cosT;
    P[1][2] = u[1] * u[2] * ccosT - u[0] * sinT;
    P[2][0] = u[2] * u[0] * ccosT - u[1] * sinT;
    P[2][1] = u[2] * u[1] * ccosT + u[0] * sinT;
    P[2][2] = u[2] * u[2] * ccosT + cosT;

    // Multiply P*m_P, saving in pNew, then copy back to m_P.
    Matrx pNew = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for(int r = 0; r < 3; ++r)
    {
      for(int c = 0; c < 3; ++c)
      {
        for(int k = 0; k < 3; ++k)
        {
          pNew[r][c] += P[r][k] * m_P[k][c];
        }
      }
    }
    for(int r = 0; r < 3; ++r)
    {
      for(int c = 0; c < 3; ++c)
      {
        m_P[r][c] = pNew[r][c];
      }
    }

    // Change to m_v.
    Vectr vOld = m_v;
    multMatrxVectr(P, vOld, m_v);
  }

  //! @brief Get a transformed 3D Point.
  AXOM_HOST_DEVICE axom::primal::Point<T, 3> getTransformed(const axom::primal::Point<T, 3>& pt) const
  {
    axom::primal::Point<T, 3> rval = pt;
    transform(rval[0], rval[1], rval[2]);
    return rval;
  }

  //! @brief Get a transformed 3D coordinate.
  AXOM_HOST_DEVICE axom::NumericArray<T, 3> getTransformed(const axom::NumericArray<T, 3>& pt) const
  {
    axom::NumericArray<T, 3> rval = pt;
    transform(rval[0], rval[1], rval[2]);
    return rval;
  }

  //! @brief Transform a 3D coordinate in place.
  AXOM_HOST_DEVICE void transform(axom::NumericArray<T, 3>& pt) const
  {
    transform(pt[0], pt[1], pt[2]);
  }

  //! @brief Transform a 3D coordinate in place.
  AXOM_HOST_DEVICE void transform(T& x, T& y, T& z) const
  {
    Vectr tmpPt = {x, y, z};
    const auto& P(m_P);  // shorthand
    x = P[0][0] * tmpPt[0] + P[0][1] * tmpPt[1] + P[0][2] * tmpPt[2] + m_v[0];
    y = P[1][0] * tmpPt[0] + P[1][1] * tmpPt[1] + P[1][2] * tmpPt[2] + m_v[1];
    z = P[2][0] * tmpPt[0] + P[2][1] * tmpPt[1] + P[2][2] * tmpPt[2] + m_v[2];
  }

  /*!
   * @brief Invert the transformation in place.
   *
   * Using a special inverse formula for 4x4 matrices with last row [0,0,0,1].
   * @verbatim
   * Minv = [ Pinv - Pinv*v ]
   *        [  0       1    ]
   * @endverbatim
   */
  AXOM_HOST_DEVICE void invert()
  {
    invertMatrx(m_P);
    auto negV = -m_v;
    multMatrxVectr(m_P, negV, m_v);
  }

  //! brief Get the inverse transformer.
  CoordinateTransformer getInverse() const
  {
    auto rval = *this;
    rval.invert();
    return rval;
  }

private:
  //!@brief 3-vector or 3D coordinates, depending on context.
  using Vectr = axom::NumericArray<T, 3>;

  //!@brief 3x3 transformation matrix with row-major storage.
  using Matrx = Vectr[3];

  /*
   * The 4x4 matrix is saved as a 3x3 matrix P and a 3x1 vector v.
   * Last row is not stored because it's always [0,0,0,1].
   * M = [ P v ]
   *     [ 0 1 ]
   *
   * Store m_P[0][0] = std::numeric_limits<T>::quiet_NaN() to set this
   * CoordinateTransformer as invalid.  See \a setInvalid()
   */
  Matrx m_P;
  Vectr m_v;

  // Invert a Matrx, or set first value to NaN if not invertible.
  AXOM_HOST_DEVICE static void invertMatrx(Matrx& m)
  {
    T a = m[0][0];
    T b = m[0][1];
    T c = m[0][2];
    T d = m[1][0];
    T e = m[1][1];
    T f = m[1][2];
    T g = m[2][0];
    T h = m[2][1];
    T i = m[2][2];
    T det = axom::numerics::determinant(a, b, c, d, e, f, g, h, i);
    if(det != 0.0)
    {
      T detInv = 1 / det;
      m[0][0] = detInv * (e * i - f * h);
      m[0][1] = -detInv * (b * i - c * h);
      m[0][2] = detInv * (b * f - c * e);
      m[1][0] = -detInv * (d * i - f * g);
      m[1][1] = detInv * (a * i - c * g);
      m[1][2] = -detInv * (a * f - c * d);
      m[2][0] = detInv * (d * h - e * g);
      m[2][1] = -detInv * (a * h - b * g);
      m[2][2] = detInv * (a * e - b * d);
    }
    else
    {
      m[0][0] = std::numeric_limits<T>::quiet_NaN();
    }
  }

  AXOM_HOST_DEVICE static void multMatrx(const Matrx& A, const Matrx& B, Matrx& prod)
  {
    for(int r = 0; r < 3; ++r)
    {
      for(int c = 0; c < 3; ++c)
      {
        prod[r][c] = 0.0;
        for(int k = 0; k < 3; ++k)
        {
          prod[r][c] += A[r][k] * B[k][c];
        }
      }
    }
  }

  AXOM_HOST_DEVICE static void multMatrxVectr(const Matrx& A, const Vectr& b, Vectr& prod)
  {
    for(int r = 0; r < 3; ++r)
    {
      prod[r] = dotProd(A[r], b);
    }
  }

  AXOM_HOST_DEVICE T static dotProd(const Vectr& u, const Vectr& v)
  {
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  }

  AXOM_HOST_DEVICE void copyIn(const CoordinateTransformer& other)
  {
    m_v = other.m_v;
    for(int r = 0; r < 3; ++r)
    {
      m_P[r] = other.m_P[r];
    }
  }
};

}  // namespace experimental
}  // namespace primal
}  // namespace axom

#endif
