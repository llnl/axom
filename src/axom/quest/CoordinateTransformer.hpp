// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_COORDINATE_TRANSFORMER_HPP
#define AXOM_QUEST_COORDINATE_TRANSFORMER_HPP

#include "axom/core/numerics/Matrix.hpp"

#if defined(AXOM_USE_MPI)
  #include "mpi.h"
#endif

namespace axom
{
namespace quest
{
/*!
  @brief For applying coordinate transformation using a
  4x4 transformation matrix operating on homogeneous coordinates.

  Transformation matrix should have the last row values
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
    @brief Constructor sets the 4x4 transformation matrix.
    @param matrix [in] The transformation matrix for homogeneous
    coordinates.
  */
  AXOM_HOST_DEVICE CoordinateTransformer(const numerics::Matrix<double>& matrix)
  {
    setMatrix(matrix);
  }

  /*!
    @brief Set the matrix, discarding the current one.
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

  //! @brief Transform a 3D coordinate.
  AXOM_HOST_DEVICE axom::primal::Point<T, 3> getTransform(const axom::primal::Point<T, 3>& pt) const
  {
    axom::primal::Point<T, 3> rval = pt;
    transform(rval[0], rval[1], rval[2]);
    return rval;
  }

  //! @brief Transform a 3D coordinate in place.
  AXOM_HOST_DEVICE void transform(axom::primal::Point<T, 3>& pt) const
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
    Minv = [ Pinv - Pinv*v ]
           [  0       1    ]
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
    m_v[0] = -m_P[0][0] * v0 + m_P[0][1] * v1 + m_P[0][2] * v2;
    m_v[1] = -m_P[1][0] * v0 + m_P[1][1] * v1 + m_P[1][2] * v2;
    m_v[2] = -m_P[2][0] * v0 + m_P[2][1] * v1 + m_P[2][2] * v2;
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
};

}  // namespace quest
}  // namespace axom

#endif
