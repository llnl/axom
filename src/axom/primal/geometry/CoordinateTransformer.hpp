// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_COORDINATE_TRANSFORMER_HPP
#define AXOM_PRIMAL_COORDINATE_TRANSFORMER_HPP

#include "axom/core/numerics/Matrix.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"

namespace axom
{
namespace primal
{
/*!
  @brief Coordinate transformation facilitating the placement of
  geometries whose parameters can't easily describe it.

  The transformations may be described as translations, rotations
  and arbitrary operators transforms on homogeneous coordinates.
  These matrices should be 4x4 and have the last row values
  [0, 0, 0, 1].

  Only supporting 3D coordinates presently.
  This class is a new utility.  It is subject to change.
*/
template<typename T = double>
class CoordinateTransformer
{
public:
  /*!
    @brief Default constructor sets an identity transformer.
  */
  AXOM_HOST_DEVICE CoordinateTransformer()
  : m_P{ {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.} }
  , m_v{0., 0., 0.}
  { }

  /*!
    @brief Copy constructor.
  */
  AXOM_HOST_DEVICE CoordinateTransformer(const CoordinateTransformer& other)
  {
    copyIn(other);
  }

  AXOM_HOST_DEVICE CoordinateTransformer& operator=(const CoordinateTransformer& other)
  {
    copyIn(other);
    return *this;
  }

  /*!
    @brief Constructor sets the 4x4 transformation matrix.
    @param matrix [in] The transformation matrix for homogeneous
    coordinates.
  */
  AXOM_HOST_DEVICE CoordinateTransformer(const numerics::Matrix<double>& matrix)
  {
    setMatrix(matrix);
  }

  /*!
    @brief Set the matrix, discarding the current transformation.
    @param matrix [in] The transformation matrix for homogeneous
    coordinates.
  */
  void setMatrix(const numerics::Matrix<double>& matrix)
  {
    // Assert that matrix is a transformation in homogeneous coordinates.
    SLIC_ASSERT(matrix.getNumRows() == 4);
    SLIC_ASSERT(matrix.getNumColumns() == 4);
    SLIC_ASSERT(matrix(3,0) == 0.0);
    SLIC_ASSERT(matrix(3,1) == 0.0);
    SLIC_ASSERT(matrix(3,2) == 0.0);
    SLIC_ASSERT(matrix(3,3) == 1.0);
    for (int c = 0; c < 3; ++c )
    {
      for (int r = 0; r < 3; ++r )
      {
        m_P[r][c] = matrix(r, c);
      }
    }
    for (int r = 0; r < 3; ++r )
    {
      m_v[r] = matrix(r,3);
    }
  }

  /*!
    @brief Get the matrix for the transformation.
  */
  numerics::Matrix<double> getMatrix()
  {
    numerics::Matrix<double> rval(4, 4, 0.0);
    for (int c = 0; c < 3; ++c )
    {
      for (int r = 0; r < 3; ++r )
      {
        rval(r, c) = m_P[r][c];
      }
    }
    for (int r = 0; r < 3; ++r )
    {
      rval(r,3) = m_v[r];
    }
  }

  /*!
    @brief Add a matrix transform to the current transformation.
    @param matrix [in] The transformation matrix for homogeneous
    coordinates.
  */
  void addMatrix(const numerics::Matrix<double>& matrix)
  {
    numerics::Matrix<double> current = getMatrix();
    numerics::Matrix<double> updated(4, 4);
    axom::numerics::matrix_multiply(matrix, current, updated);
    setMatrix(updated);
  }

  //! @brief Add a 3D translation to the current transformation.
  void addTranslation(const axom::primal::Vector<T, 3>& d)
  {
    addTranslation(d.array());
  }

  //! @brief Add a 3D translation to the current transformation.
  void addTranslation(const axom::NumericArray<T, 3>& d)
  {
    m_v[0] += d[0];
    m_v[1] += d[1];
    m_v[2] += d[2];
  }

  /*!
    @brief Add a 3D rotation to the current transformation.
    The rotation is equivalent to rotating a vector
    from one direction to another.

    @param start [in] Starting direction
    @param end [in] Ending direction
  */
  void addRotation(const axom::primal::Vector<T, 3>& start,
                   const axom::primal::Vector<T, 3>& end)
  {
    // Note that the rotation matrix is not unique.
    Vector<T, 3> s = start.unitVector();
    Vector<T, 3> e = end.unitVector();
    Vector<T, 3> u;  // Rotation vector, the cross product of start and end.
    numerics::cross_product(s.data(), e.data(), u.data());
    T sinT = u.norm();
    T cosT = numerics::dot_product(s.data(), e.data(), 3);

    // Degenerate: end parallel to start, angle near 0 or pi.
    if(utilities::isNearlyEqual(sinT, 0.0))
    {
      if(cosT < 0)
      {
        // Negative identity transform.  Change signs.
        for(int r=0; r<3; ++r)
        {
          m_v[r] = -m_v[r];
          for(int c=0; c<3; ++c)
          {
            m_P[r][c] = -m_P[r][c];
          }
        }
      }
      return;
    }

    u = u.unitVector();
    privateAddRotation(u, sinT, cosT);
  }

  /*!
    @brief Add a 3D rotation to the current transformation.
    The rotation is given as a rotation axis and an angle.

    @param axisDir [in]
    @param angle [in]
  */
  void addRotation(const axom::primal::Vector<T, 3>& axisDir,
                   T angle)
  {
    SLIC_ASSERT(axisDir.squared_norm() > 1e-20);
    auto unitDir = axisDir.unitVector();
    T sinT = sin(angle);
    T cosT = cos(angle);
    privateAddRotation(unitDir, sinT, cosT);
  }

  //!@brief Add rotation, given unit vector and angle as sine and cosine.
  void privateAddRotation(const axom::primal::Vector<T, 3>& u, T sinT, T cosT)
  {
  double ccosT = 1 - cosT;

  double P[3][3]; // 3D rotation matrix.
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
  double pNew[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  for(int r=0; r<3; ++r)
  {
    for(int c=0; c<3; ++c)
    {
      for(int k=0; k<3; ++k)
      {
        pNew[r][c] += P[r][k] * m_P[k][c];
      }
    }
  }
  for(int r=0; r<3; ++r)
  {
    for(int c=0; c<3; ++c)
    {
      m_P[r][c] = pNew[r][c];
    }
  }

  // Change to m_v.
  T vOld[3] = {m_v[0], m_v[1], m_v[2]};
  m_v[0] = P[0][0]*vOld[0] + P[0][1]*vOld[1] + P[0][2]*vOld[2];
  m_v[1] = P[1][0]*vOld[0] + P[1][1]*vOld[1] + P[1][2]*vOld[2];
  m_v[2] = P[2][0]*vOld[0] + P[2][1]*vOld[1] + P[2][2]*vOld[2];
  }

  //! @brief Get a trransformed 3D Point.
  AXOM_HOST_DEVICE axom::primal::Point<T, 3> getTransformed(const axom::primal::Point<T, 3>& pt) const
  {
    axom::primal::Point<T, 3> rval = pt;
    transform(rval[0], rval[1], rval[2]);
    return rval;
  }

  //! @brief Get a trransformed 3D Vector.
  AXOM_HOST_DEVICE axom::primal::Vector<T, 3> getTransformed(const axom::primal::Vector<T, 3>& pt) const
  {
    axom::primal::Vector<T, 3> rval = pt;
    transform(rval[0], rval[1], rval[2]);
    return rval;
  }

  //! @brief Get a trransformed 3D coordinate.
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
    double tmpPt[3] = {x, y, z};
    const auto& P(m_P); // shorthand
    x = P[0][0]*tmpPt[0] + P[0][1]*tmpPt[1] + P[0][2]*tmpPt[2] + m_v[0];
    y = P[1][0]*tmpPt[0] + P[1][1]*tmpPt[1] + P[1][2]*tmpPt[2] + m_v[1];
    z = P[2][0]*tmpPt[0] + P[2][1]*tmpPt[1] + P[2][2]*tmpPt[2] + m_v[2];
  }

  /*!
    @brief Invert the transformation in place.

    We use a special inverse formula for 4x4 matrices with last row [0,0,0,1].
    @verbatim
    Minv = [ Pinv - Pinv*v ]
           [  0       1    ]
    @endverbatim
  */
  AXOM_HOST_DEVICE void invert()
  {
    double a = m_P[0][0];
    double b = m_P[0][1];
    double c = m_P[0][2];
    double d = m_P[1][0];
    double e = m_P[1][1];
    double f = m_P[1][2];
    double g = m_P[2][0];
    double h = m_P[2][1];
    double i = m_P[2][2];
    double detP = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
    double detPInv = 1/detP;
    m_P[0][0] =  detPInv * (e*i - f*h);
    m_P[0][1] = -detPInv * (b*i - c*h);
    m_P[0][2] =  detPInv * (b*f - c*e);
    m_P[1][0] = -detPInv * (d*i - f*g);
    m_P[1][1] =  detPInv * (a*i - c*g);
    m_P[1][2] = -detPInv * (a*f - c*d);
    m_P[2][0] =  detPInv * (d*h - e*g);
    m_P[2][1] = -detPInv * (a*h - b*g);
    m_P[2][2] =  detPInv * (a*e - b*d);
    double v0 = m_v[0];
    double v1 = m_v[1];
    double v2 = m_v[2];
    m_v[0] = -(m_P[0][0] * v0 + m_P[0][1] * v1 + m_P[0][2] * v2);
    m_v[1] = -(m_P[1][0] * v0 + m_P[1][1] * v1 + m_P[1][2] * v2);
    m_v[2] = -(m_P[2][0] * v0 + m_P[2][1] * v1 + m_P[2][2] * v2);
  }

  //! brief Get the inverse transformer.
  CoordinateTransformer getInverse() const
  {
    auto rval = *this;
    rval.invert();
    return rval;
  }

private:
  /*
    The 4x4 matrix is saved as a 3x3 matrix P and a 3x1 vector V.
    Last row is not stored because it's always [0,0,0,1].
    M = [ P v ]
        [ 0 1 ]
  */
  double m_P[3][3];
  double m_v[3];

  AXOM_HOST_DEVICE void copyIn(const CoordinateTransformer& other)
  {
    for(int r = 0; r < 3; ++r)
    {
      m_v[r] = other.m_v[r];
      for(int c = 0; c < 3; ++c)
      {
        m_P[r][c] = other.m_P[r][c];
      }
    }
  }
};

}  // namespace primal
}  // namespace axom

#endif
