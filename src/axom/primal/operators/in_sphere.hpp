// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
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
#include "axom/primal/geometry/OrientationResult.hpp"

#include "axom/primal/operators/detail/predicate_determinants.hpp"

namespace axom
{
namespace primal
{
namespace robust
{
/*!
 * \brief Classifies a query point against a 2D triangle's circumcircle.
 *
 * \return ON_NEGATIVE_SIDE if inside, ON_POSITIVE_SIDE if outside, ON_BOUNDARY otherwise.
 */
template <typename T>
inline int in_sphere(const Point<T, 2>& q,
                     const Point<T, 2>& p0,
                     const Point<T, 2>& p1,
                     const Point<T, 2>& p2,
                     double EPS = 1e-8)
{
  const double det = detail::in_sphere_determinant(q, p0, p1, p2);
  if(axom::utilities::isNearlyEqual(det, 0., EPS))
  {
    return primal::ON_BOUNDARY;
  }

  return det < 0. ? primal::ON_NEGATIVE_SIDE : primal::ON_POSITIVE_SIDE;
}

template <typename T>
inline int in_sphere(const Point<T, 2>& q, const Triangle<T, 2>& tri, double EPS = 1e-8)
{
  return robust::in_sphere(q, tri[0], tri[1], tri[2], EPS);
}

/*!
 * \brief Classifies a query point against a 3D tetrahedron's circumsphere.
 *
 * \return ON_NEGATIVE_SIDE if inside, ON_POSITIVE_SIDE if outside, ON_BOUNDARY otherwise.
 */
template <typename T>
inline int in_sphere(const Point<T, 3>& q,
                     const Point<T, 3>& p0,
                     const Point<T, 3>& p1,
                     const Point<T, 3>& p2,
                     const Point<T, 3>& p3,
                     double EPS = 1e-8)
{
  const double det = detail::in_sphere_determinant(q, p0, p1, p2, p3);
  if(axom::utilities::isNearlyEqual(det, 0., EPS))
  {
    return primal::ON_BOUNDARY;
  }

  return det < 0. ? primal::ON_NEGATIVE_SIDE : primal::ON_POSITIVE_SIDE;
}

template <typename T>
inline int in_sphere(const Point<T, 3>& q, const Tetrahedron<T, 3>& tet, double EPS = 1e-8)
{
  return robust::in_sphere(q, tet[0], tet[1], tet[2], tet[3], EPS);
}

}  // namespace robust

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
 * \param [in] includeBoundary if true, points on the circumcircle are treated
 *  as inside. Default: false.
 * \return true if the point is inside the circumcircle, false if it is on
 * the circle's boundary or outside the circle
 */
template <typename T>
inline bool in_sphere(const Point<T, 2>& q,
                      const Point<T, 2>& p0,
                      const Point<T, 2>& p1,
                      const Point<T, 2>& p2,
                      double EPS = 1e-8,
                      bool includeBoundary = false)
{
  const int res = robust::in_sphere(q, p0, p1, p2, EPS);
  return includeBoundary ? (res != primal::ON_POSITIVE_SIDE) : (res == primal::ON_NEGATIVE_SIDE);
}

/*!
 * \brief Tests whether a query point lies inside a 2D triangle's circumcircle
 *
 * \param [in] q the query point
 * \param [in] tri the triangle
 * \param [in] EPS tolerance for determining if \a q is on the boundary. Default: 1e-8.
 * \param [in] includeBoundary if true, points on the circumcircle are treated
 *  as inside. Default: false.
 * \see in_sphere
 */
template <typename T>
inline bool in_sphere(const Point<T, 2>& q,
                      const Triangle<T, 2>& tri,
                      double EPS = 1e-8,
                      bool includeBoundary = false)
{
  return in_sphere(q, tri[0], tri[1], tri[2], EPS, includeBoundary);
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
 * \param [in] includeBoundary if true, points on the circumsphere are treated
 *  as inside. Default: false.
 * \return true if the point is inside the circumsphere, false if it is on
 * the sphere's boundary or outside the sphere
 */
template <typename T>
inline bool in_sphere(const Point<T, 3>& q,
                      const Point<T, 3>& p0,
                      const Point<T, 3>& p1,
                      const Point<T, 3>& p2,
                      const Point<T, 3>& p3,
                      double EPS = 1e-8,
                      bool includeBoundary = false)
{
  const int res = robust::in_sphere(q, p0, p1, p2, p3, EPS);
  return includeBoundary ? (res != primal::ON_POSITIVE_SIDE) : (res == primal::ON_NEGATIVE_SIDE);
}

/*!
 * \brief Tests whether a query point lies inside a 3D tetrahedron's circumsphere
 *
 * \param [in] q the query point
 * \param [in] tet the tetrahedron
 * \param [in] EPS tolerance for determining if \a q is on the boundary. Default: 1e-8.
 * \param [in] includeBoundary if true, points on the circumsphere are treated
 *  as inside. Default: false.
 * \see in_sphere
 */
template <typename T>
inline bool in_sphere(const Point<T, 3>& q,
                      const Tetrahedron<T, 3>& tet,
                      double EPS = 1e-8,
                      bool includeBoundary = false)
{
  return in_sphere(q, tet[0], tet[1], tet[2], tet[3], EPS, includeBoundary);
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
