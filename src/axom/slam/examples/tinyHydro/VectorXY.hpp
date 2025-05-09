// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// 2D vector for XY geometry
// Fri Nov 21 10:50:53 PST 2014
#include <math.h>

#ifndef	_VECTORXY_H
#define	_VECTORXY_H 1

namespace tinyHydro {

  class VectorXY
  {
  public:
    explicit VectorXY(double xx = 0., double yy = 0.) : x(xx), y(yy) {}

    double x;
    double y;

  public:

    //----------------------------------------------
    /// Addition operations
    VectorXY& accum(const VectorXY & b)
    {
      x += b.x;
      y += b.y;
      return *this;
    }

    VectorXY& operator+=(const VectorXY & b)
    {
      return accum(b);
    }

    VectorXY add(const VectorXY& a) const
    {
      return VectorXY(x + a.x, y + a.y);
    }

    //----------------------------------------------
    /// Subtraction operations
    VectorXY& elim(const VectorXY& b)
    {
      x -= b.x;
      y -= b.y;
      return *this;
    }

    VectorXY& operator-=(const VectorXY & b)
    {
      return elim(b);
    }


    VectorXY sub(const VectorXY& b) const
    {
      return VectorXY(x - b.x, y - b.y);
    }

    VectorXY perp() const
    {
      return VectorXY(-y, x);
    }

    //----------------------------------------------
    // Scalar multiplication operations

    VectorXY& scale(const double s)
    {
      x *= s;
      y *= s;
      return *this;
    }

    VectorXY& operator*=(const double s)
    {
      return scale(s);
    }

    VectorXY mult(double s) const
    {
      return VectorXY(s * x, s * y);
    }


    //----------------------------------------------
    // Vector product
    double dot(const VectorXY & v) const
    {
      return x * v.x + y * v.y;
    }

    // Perp-dot product of a 2D vector -- equivalent to 'return perp().dot(v)'
    double cross(const VectorXY & v) const
    {
      return x * v.y - y * v.x;
    }

    //----------------------------------------------
    double mag() const
    {
      return sqrt(mag2());
    }

    // Magnitude of the vector -- Equivalent to dot(*this)
    double mag2() const
    {
      return x * x + y * y;
    }

//----------------------------------------------
  };


// Free functions
  inline VectorXY operator+(const VectorXY& a, const VectorXY& b)
  {
    VectorXY v(a);

    v += b;
    return v;
  }

  inline VectorXY operator-(const VectorXY& a, const VectorXY& b)
  {
    VectorXY v(a);

    v -= b;
    return v;
  }

  inline VectorXY operator*(const VectorXY& a, const double s)
  {
    VectorXY v(a);

    v *= s;
    return v;
  }
  inline VectorXY operator*(const double s, const VectorXY& a)
  {
    VectorXY v(a);

    v *= s;
    return v;
  }




} // end namespace tinyHydro

#endif
