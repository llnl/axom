// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file in_sphere.hpp
 *
 * \brief Consists of methods that test whether a given query point is
 * inside the unique sphere circumscribing a 2D triangle or a 3D tetrahedron.
 *
 * This is a well known computational geometry primitive.  For reference,
 * see Section 3.1.6.4 in "Real-time collision detection" by C. Ericson.
 */

#ifndef AXOM_PRIMAL_IN_SPHERE_H_
#define AXOM_PRIMAL_IN_SPHERE_H_

#include "axom/core.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"

namespace axom
{
namespace primal
{
/*!
 * \brief Tests whether a query point lies inside a 2D triangle's circumcircle
 *
 * A triangle's circumcircle is the unique circle (i.e. a 2-sphere) that
 * passes through each of its three vertices.
 *
 * \param [in] q the query point
 * \param [in] p0 the first vertex of the triangle
 * \param [in] p1 the second vertex of the triangle
 * \param [in] p2 the third vertex of the triangle
 * \param [in] EPS tolerance for determining if \a q is on the boundary. Default: 1e-8.
 * \return true if the point is inside the circumcircle, false if it is on
 * the circle's boundary or outside the circle
 */
template <typename T>
inline bool in_sphere(const Point<T, 2>& q,
                      const Point<T, 2>& p0,
                      const Point<T, 2>& p1,
                      const Point<T, 2>& p2,
                      double EPS = 1e-8)
{
  const auto ba = p1 - p0;
  const auto ca = p2 - p0;
  const auto qa = q - p0;

  // clang-format off
  const double det = axom::numerics::determinant(
    ba[0], ba[1], ba.squared_norm(),
    ca[0], ca[1], ca.squared_norm(),
    qa[0], qa[1], qa.squared_norm());
  // clang-format on

  return axom::utilities::isNearlyEqual(det, 0., EPS) ? false : (det < 0);
}

/*!
 * \brief Tests whether a query point lies inside a 2D triangle's circumcircle
 *
 * \param [in] q the query point
 * \param [in] tri the triangle
 * \param [in] EPS tolerance for determining if \a q is on the boundary. Default: 1e-8.
 * \see in_sphere
 */
template <typename T>
inline bool in_sphere(const Point<T, 2>& q, const Triangle<T, 2>& tri, double EPS = 1e-8)
{
  return in_sphere(q, tri[0], tri[1], tri[2], EPS);
}

/*!
 * \brief Tests whether a query point lies inside a 3D tetrahedron's
 * circumsphere
 *
 * A tetrahedron's circumsphere is the unique sphere that passes through each
 * of its four vertices.
 *
 * \param [in] q the query point
 * \param [in] p0 the first vertex of the tetrahedron
 * \param [in] p1 the second vertex of the tetrahedron
 * \param [in] p2 the third vertex of the tetrahedron
 * \param [in] p3 the fourth vertex of the tetrahedron
 * \param [in] EPS tolerance for determining if \a q is on the boundary. Default: 1e-8.
 * \return true if the point is inside the circumsphere, false if it is on
 * the sphere's boundary or outside the sphere
 */
template <typename T>
inline bool in_sphere(const Point<T, 3>& q,
                      const Point<T, 3>& p0,
                      const Point<T, 3>& p1,
                      const Point<T, 3>& p2,
                      const Point<T, 3>& p3,
                      double EPS = 1e-8)
{
  const auto ba = p1 - p0;
  const auto ca = p2 - p0;
  const auto da = p3 - p0;
  const auto qa = q - p0;

  // clang-format off
  const double det = axom::numerics::determinant(
    ba[0], ba[1], ba[2], ba.squared_norm(),
    ca[0], ca[1], ca[2], ca.squared_norm(),
    da[0], da[1], da[2], da.squared_norm(),
    qa[0], qa[1], qa[2], qa.squared_norm());
  // clang-format on

  return axom::utilities::isNearlyEqual(det, 0., EPS) ? false : (det < 0);
}

/*!
 * \brief Tests whether a query point lies inside a 3D tetrahedron's circumsphere
 *
 * \param [in] q the query point
 * \param [in] tet the tetrahedron
 * \param [in] EPS tolerance for determining if \a q is on the boundary. Default: 1e-8.
 * \see in_sphere
 */
template <typename T>
inline bool in_sphere(const Point<T, 3>& q, const Tetrahedron<T, 3>& tet, double EPS = 1e-8)
{
  return in_sphere(q, tet[0], tet[1], tet[2], tet[3], EPS);
}

/*!
 * \brief Tests whether a bounding box lies inside a 2D sphere
 * 
 * \param [in] bb the bounding box
 * \param [in] circle the sphere
 */
template <typename T>
inline bool in_sphere(const BoundingBox<T, 2>& bb, const Sphere<T, 2>& circle)
{
  // Check if any corner of the bounding box is outside the sphere.
  //  This version requires a lot of multiplications, but no square roots
  //  and a higher likelihood of early returns.
  auto the_max = bb.getMax();
  auto the_min = bb.getMin();
  auto radius = circle.getRadius();
  auto center = circle.getCenter();

  if((center[0] - the_min[0]) * (center[0] - the_min[0]) +
       (center[1] - the_min[1]) * (center[1] - the_min[1]) >
     radius * radius)
  {
    return false;
  }

  if((center[0] - the_max[0]) * (center[0] - the_max[0]) +
       (center[1] - the_min[1]) * (center[1] - the_min[1]) >
     radius * radius)
  {
    return false;
  }

  if((center[0] - the_min[0]) * (center[0] - the_min[0]) +
       (center[1] - the_max[1]) * (center[1] - the_max[1]) >
     radius * radius)
  {
    return false;
  }

  if((center[0] - the_max[0]) * (center[0] - the_max[0]) +
       (center[1] - the_max[1]) * (center[1] - the_max[1]) >
     radius * radius)
  {
    return false;
  }

  return true;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_IN_SPHERE_H_
