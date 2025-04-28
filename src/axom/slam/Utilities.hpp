// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file
 * \brief A few utility functions used by the SLAM component.
 */
#ifndef SLAM_UTILITIES_H_
#define SLAM_UTILITIES_H_

#include "axom/core.hpp"
#include "axom/fmt.hpp"

#include <string>
#include <iostream>

namespace axom
{
namespace slam
{
using DefaultPositionType = axom::IndexType;
using DefaultElementType = axom::IndexType;

class NotImplementedException
{ };

namespace util
{
/** \brief A helper class to print the name of a few types */
template <typename T>
struct TypeToString
{
  static std::string to_string() { return "<unspecialized>"; }
};

/** \brief A helper class to print the name of integers as 'int' */
template <>
struct TypeToString<int>
{
  static std::string to_string() { return "int"; }
};

/** \brief A helper class to print the name of doubles as 'double' */
template <>
struct TypeToString<double>
{
  static std::string to_string() { return "double"; }
};

/**
 * \brief A simple 3D point class similar to primal's point class,
 * with some basic Point/Vector functionalities
 * 
 * \note This is needed for internal testing in slam (which does not depend on primal)
 */
template <typename DataType = double>
struct Point3
{
  Point3(const DataType& x, const DataType& y, const DataType& z) : m_x(x), m_y(y), m_z(z) { }

  Point3(const DataType* d) : m_x(d[0]), m_y(d[1]), m_z(d[2]) { }
  Point3() : m_x(DataType()), m_y(DataType()), m_z(DataType()) { }

  /// Distance from origin
  DataType radius() const { return std::sqrt(m_x * m_x + m_y * m_y + m_z * m_z); }

  Point3& operator+=(const Point3& pt)
  {
    m_x += pt.m_x;
    m_y += pt.m_y;
    m_z += pt.m_z;
    return *this;
  }

  Point3& operator-=(const Point3& pt)
  {
    m_x -= pt.m_x;
    m_y -= pt.m_y;
    m_z -= pt.m_z;
    return *this;
  }

  Point3& operator*=(const DataType& sc)
  {
    m_x *= sc;
    m_y *= sc;
    m_z *= sc;
    return *this;
  }

  template <typename T>
  Point3& operator/=(const T& sc)
  {
    return operator*=(1. / sc);
  }

  /// access the xyz values
  DataType operator[](unsigned int i) const
  {
    DataType vals[] = {m_x, m_y, m_z};
    return vals[i];
  }

  bool operator==(const Point3& pt) const
  {
    return pt.m_x == m_x && pt.m_y == m_y && pt.m_z == m_z;
  }

  friend std::ostream& operator<<(std::ostream& os, const Point3& pt)
  {
    return os << "(" << pt.m_x << "," << pt.m_y << "," << pt.m_z << ")";
  }

  DataType m_x, m_y, m_z;
};

template <typename T>
Point3<T> operator+(const Point3<T>& pt1, const Point3<T>& pt2)
{
  Point3<T> pt(pt1);
  pt += pt2;
  return pt;
}
template <typename T>
Point3<T> operator-(const Point3<T>& pt1, const Point3<T>& pt2)
{
  Point3<T> pt(pt1);
  pt += pt2;
  return pt;
}
template <typename T>
Point3<T> normalize(const Point3<T>& pt)
{
  Point3<T> pr(pt);
  T r = pt.radius();
  if(r > 0)
  {
    pr /= r;
  }
  return pr;
}

template <typename T>
Point3<T> cross(const Point3<T>& p1, const Point3<T>& p2)
{
  Point3<T> r;
  r.m_x = p1.m_y * p2.m_z - p1.m_z * p2.m_y;
  r.m_y = p1.m_z * p2.m_x - p1.m_x * p2.m_z;
  r.m_z = p1.m_x * p2.m_y - p1.m_y * p2.m_x;
  r = normalize(r);
  return r;
}
template <typename T>
T distance(const Point3<T>& pt1, const Point3<T>& pt2)
{
  return Point3<T>(pt1 - pt2).radius();
}

}  // end namespace util
}  // end namespace slam
}  // end namespace axom

/// Overload to format an axom::slam::util::Point3 using fmt
template <typename DataType>
struct axom::fmt::formatter<axom::slam::util::Point3<DataType>> : ostream_formatter
{ };

#endif  //  SLAM_UTILITIES_H_
