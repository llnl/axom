// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file intersect_bezier_impl.hpp
 *
 * This file provides helper functions for testing the intersection
 * of Bezier curves with other Bezier curves and other geometric objects
 */

#ifndef AXOM_PRIMAL_INTERSECT_BEZIER_IMPL_HPP_
#define AXOM_PRIMAL_INTERSECT_BEZIER_IMPL_HPP_

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"

#include "axom/primal/operators/in_sphere.hpp"
#include "axom/primal/operators/detail/intersect_impl.hpp"
#include "axom/primal/operators/detail/intersect_ray_impl.hpp"

#include <vector>
#include <math.h>

namespace axom
{
namespace primal
{
namespace detail
{
//---------------------------- FUNCTION DECLARATIONS ---------------------------

/*!
 * \brief Recursive function to find the intersections between two Bezier curves
 *
 * \param [in] c1 The first input curve
 * \param [in] c2 The second input curves
 * \param [out] sp Parametric coordinates of intersections in \a c1
 * \param [out] tp Parametric coordinates of intersections in \a c2
 * \param [in] sq_tol The squared tolerance parameter for determining if a
 * Bezier curve is linear
 * \param [in] order1 The order of \a c1
 * \param [in] order2 The order of \a c2
 * \param [in] s_offset The offset in parameter space for \a c1
 * \param [in] s_scale The scale in parameter space for \a c1
 * \param [in] t_offset The offset in parameter space for \a c2
 * \param [in] t_scale The scale in parameter space for \a c2
 *
 * Bezier curves can only intersect when their bounding boxes intersect.
 * The base case of the recursion is when we can approximate the curves as
 * line segments, where we directly find their intersections. Otherwise,
 * check for intersections recursively after bisecting one of the curves.
 *
 * \note A BezierCurve is parametrized in [0,1). The scale and offset parameters
 * are used to track the local curve parameters during subdivisions
 *
 * \return True if the two curves intersect, False otherwise
 * \sa intersect_bezier
 */
template <typename T>
bool intersect_bezier_curves(const BezierCurve<T, 2> &c1,
                             const BezierCurve<T, 2> &c2,
                             axom::Array<T> &sp,
                             axom::Array<T> &tp,
                             double sq_tol,
                             int order1,
                             int order2,
                             double s_offset,
                             double s_scale,
                             double t_offset,
                             double t_scale);

/*!
 * \brief Tests intersection of two line segments defined by
 * their end points (a,b) and (c,d)
 *
 * \param [in] a,d,c,b the endpoints of the segments
 * \param [out] s,t The parametrized s and t values at which intersection occurs
 * Range of output values for \a s and \a t is [0,1).
 *
 * \return True, if the two line segments intersect, false otherwise.
 * When the two line segments intersect, their intersection parameters
 * are stored in output parameters \a s and \a t, respectively, such that,
 * if \a p is the intersection point,
 * \f$    p = a + s * (b-a) = c + t * (d-c).  \f$
 *
 * \note We consider both line segments to have half-open intervals.
 * As such, the we do not consider the lines to intersect if they do so
 * at the endpoints \a b or \d, respectively.
 *
 * \note This function assumes the all intersections have multiplicity
 * one, i.e. there are no points at which the curves and their derivatives
 * both intersect. Thus, the function does not find tangencies.
 *
 * \note This function assumes two dimensional curves in a plane.
 *
 * \note This function does not properly handle collinear lines
 */
template <typename T>
bool intersect_2d_linear(const Point<T, 2> &a,
                         const Point<T, 2> &b,
                         const Point<T, 2> &c,
                         const Point<T, 2> &d,
                         T &s,
                         T &t);

/*!
 * \brief Recursive function to find intersections between a ray and a Bezier curve
 *
 * \param [in] r The input ray
 * \param [in] c The input curve
 * \param [out] rp Parametric coordinates of intersections in \a r [0, inf)
 * \param [out] cp Parametric coordinates of intersections in \a c [0, 1)
 * \param [in] sq_tol The squared tolerance parameter for distances in physical space
 * \param [in] EPS The tolerance parameter for distances in parameter space
 * \param [in] order The order of \a c
 * \param [in] c_offset The offset in parameter space for \a c
 * \param [in] c_scale The scale in parameter space for \a c
 * \param [in] isHalfOpen If true, ignore intersections at t=1 in the parameter space of the curve
 *
 * A ray can only intersect a Bezier curve if it intersects its bounding box
 * The base case of the recursion is when we can approximate the curves parametrically with
 * line segments, where we directly find their intersection with the ray. Otherwise,
 * check for intersections recursively after bisecting the curve.
 *
 * \note A BezierCurve is parametrized in [0,1] if isHalfOpen is false, and [0, 1) if true.
 * The scale and offset parameters are used to track the local curve parameters during subdivisions
 *
 * \note This function can't be used to identify tangents at local a min/max
 *   of Bezier curves.
 * 
 * \return True if the ray and curve intersect, False otherwise
 * \sa intersect_bezier
 */
template <typename T>
bool intersect_ray_bezier(const Ray<T, 2> &r,
                          const BezierCurve<T, 2> &c,
                          axom::Array<T> &rp,
                          axom::Array<T> &cp,
                          double sq_tol,
                          double EPS,
                          int order,
                          double c_offset,
                          double c_scale,
                          bool isHalfOpen = true);

/*!
 * \brief Recursive function to find intersections between a 2D sphere (circle)
 *  and a Bezier curve
 *
 * \param [in] circle The input circle
 * \param [in] curve The input curve
 * \param [out] circle_params Parametric coordinates of intersections in \a circle [0, 2pi)
 * \param [out] curve_params Parametric coordinates of intersections in \a curve [0, 1)
 * \param [in] sq_tol The squared tolerance parameter for distances in physical space
 * \param [in] EPS The tolerance parameter for distances in parameter space
 * \param [in] order The order of \a c
 * \param [in] c_offset The offset in parameter space for \a c
 * \param [in] c_scale The scale in parameter space for \a c
 *
 * A circle can only intersect a Bezier curve if it intersects its bounding box
 * The base case of the recursion is when we can approximate the curve parametrically with
 * a line segment, where we directly find its intersection with the circle. Otherwise,
 * check for intersections recursively after bisecting the curve.
 *
 * \note A BezierCurve is parametrized in [0,1). The scale and offset parameters
 * are used to track the local curve parameters during subdivisions
 *
 * \note This function can't be used to identify tangents at local a min/max
 *   of Bezier curves.
 * 
 * \return True if the two curves intersect, False otherwise
 * \sa intersect_bezier
 */
template <typename T>
bool intersect_circle_bezier(const Sphere<T, 2> &circle,
                             const BezierCurve<T, 2> &curve,
                             axom::Array<T> &circle_params,
                             axom::Array<T> &curve_params,
                             double sq_tol,
                             double EPS,
                             int order,
                             double c_offset,
                             double c_scale);

/*!
 * \brief Tests intersection of a line and a cirlce
 *
 * \param [in] a, b the endpoints of a segment which defines the line
 * \param [out] c1, c2, t1, t2 The parametrized curve values (c) and 
 *   line values (t) at which intersection occurs.
 * Range of output values for \a c is [0, 2pi) and \a t is [-inf, inf).
 *
 * \return True, if the line segment intersects, false otherwise.
 *
 * \note This function can't be used to identify tangents at local a min/max
 *   of Bezier curves.
 */
template <typename T>
bool intersect_2d_circle_line(const Sphere<T, 2> &circ,
                              const Point<T, 2> &a,
                              const Point<T, 2> &b,
                              T &c1,
                              T &c2,
                              T &t1,
                              T &t2,
                              double EPS);
//------------------------------ IMPLEMENTATIONS ------------------------------

template <typename T>
bool intersect_bezier_curves(const BezierCurve<T, 2> &c1,
                             const BezierCurve<T, 2> &c2,
                             axom::Array<T> &sp,
                             axom::Array<T> &tp,
                             double sq_tol,
                             int order1,
                             int order2,
                             double s_offset,
                             double s_scale,
                             double t_offset,
                             double t_scale)
{
  using BCurve = BezierCurve<T, 2>;

  // Check bounding boxes to short-circuit the intersection
  if(!intersect(c1.boundingBox(), c2.boundingBox()))
  {
    return false;
  }

  bool foundIntersection = false;

  if(c1.isLinear(sq_tol, true) && c2.isLinear(sq_tol, true))
  {
    T s, t;
    if(intersect_2d_linear(c1[0], c1[order1], c2[0], c2[order2], s, t))
    {
      sp.push_back(s_offset + s_scale * s);
      tp.push_back(t_offset + t_scale * t);
      foundIntersection = true;
    }
  }
  else
  {
    constexpr double splitVal = 0.5;
    constexpr double scaleFac = 0.5;

    BCurve c3(order1);
    BCurve c4(order1);
    c1.split(splitVal, c3, c4);

    s_scale *= scaleFac;

    // Note: we want to find all intersections, so don't short-circuit
    if(intersect_bezier_curves(c2, c3, tp, sp, sq_tol, order2, order1, t_offset, t_scale, s_offset, s_scale))
    {
      foundIntersection = true;
    }
    if(intersect_bezier_curves(c2,
                               c4,
                               tp,
                               sp,
                               sq_tol,
                               order2,
                               order1,
                               t_offset,
                               t_scale,
                               s_offset + s_scale,
                               s_scale))
    {
      foundIntersection = true;
    }
  }

  return foundIntersection;
}

template <typename T>
bool intersect_2d_linear(const Point<T, 2> &a,
                         const Point<T, 2> &b,
                         const Point<T, 2> &c,
                         const Point<T, 2> &d,
                         T &s,
                         T &t)
{
  // Implementation inspired by Section 5.1.9.1 of
  // C. Ericson's Real-Time Collision Detection book

  // Note: Uses exact floating point comparisons since the subdivision algorithm
  // provides both sides of the line segments for interior curve points.

  // compute signed areas of endpoints of segment (c,d) w.r.t. segment (a,b)
  const auto area1 = twoDcross(a, b, c);
  const auto area2 = twoDcross(a, b, d);

  // early return if both have same orientation, or if d is collinear w/ (a,b)
  if(area2 == 0. || (area1 * area2) > 0.)
  {
    return false;
  }

  // compute signed areas of endpoints of segment (a,b) w.r.t. segment (c,d)
  const auto area3 = twoDcross(c, d, a);
  const auto area4 = area3 + area1 - area2;  // equivalent to twoDcross(c,d,b)

  // early return if both have same orientation, or if b is collinear w/ (c,d)
  if(area4 == 0. || (area3 * area4) > 0.)
  {
    return false;
  }

  // Compute intersection parameters using linear interpolation
  // Divisions are safe due to early return conditions
  s = area3 / (area3 - area4);
  t = area1 / (area1 - area2);

  return true;
}

template <typename T>
bool intersect_ray_bezier(const Ray<T, 2> &r,
                          const BezierCurve<T, 2> &c,
                          axom::Array<T> &rp,
                          axom::Array<T> &cp,
                          double sq_tol,
                          double EPS,
                          int order,
                          double c_offset,
                          double c_scale,
                          bool isHalfOpen)
{
  using BCurve = BezierCurve<T, 2>;

  // Check bounding box to short-circuit the intersection
  T r0, s0;

  // Need to expand the bounding box, since this ray-bb intersection routine
  //  only parameterizes the ray on (0, inf)
  T tmin = axom::numerics::floating_point_limits<T>::min();
  T tmax = axom::numerics::floating_point_limits<T>::max();

  if(!detail::intersect_ray(r, c.boundingBox().expand(10 * EPS), tmin, tmax, EPS))
  {
    return false;
  }

  bool foundIntersection = false;

  // For the base case, represent the Bezier curve as a line segment
  if(c.isLinear(sq_tol, true))
  {
    Segment<T, 2> seg(c[0], c[order]);

    // Need to check intersection with zero tolerance
    //  to handle cases where `intersect` treats the ray as collinear
    bool foundIntersection = detail::intersect_ray(r, seg, r0, s0, EPS);
    if(foundIntersection && (!isHalfOpen || s0 < 1.0 - EPS))
    {
      rp.push_back(r0);
      cp.push_back(c_offset + c_scale * s0);

      return true;
    }
  }
  else
  {
    constexpr double splitVal = 0.5;
    constexpr double scaleFac = 0.5;

    BCurve c1(order);
    BCurve c2(order);
    c.split(splitVal, c1, c2);
    c_scale *= scaleFac;

    // Note: we want to find all intersections, so don't short-circuit
    if(intersect_ray_bezier(r, c1, rp, cp, sq_tol, EPS, order, c_offset, c_scale, isHalfOpen))
    {
      foundIntersection = true;
    }
    if(intersect_ray_bezier(r, c2, rp, cp, sq_tol, EPS, order, c_offset + c_scale, c_scale, isHalfOpen))
    {
      foundIntersection = true;
    }
  }

  return foundIntersection;
}

template <typename T>
bool intersect_circle_bezier(const Sphere<T, 2> &circle,
                             const BezierCurve<T, 2> &curve,
                             axom::Array<T> &circle_p,
                             axom::Array<T> &curve_p,
                             double sq_tol,
                             double EPS,
                             int order,
                             double c_offset,
                             double c_scale)
{
  using BCurve = BezierCurve<T, 2>;

  // Other function put here to avoid circular dependency
  primal::BoundingBox<T, 2> bb = curve.boundingBox().scale(1.1);
  if(!detail::intersect_circle_bbox(circle, bb) || in_sphere(bb, circle))
  {
    return false;
  }

  bool foundIntersection = false;

  if(curve.isLinear(sq_tol, true))
  {
    T c1, c2, t1, t2;
    if(intersect_2d_circle_line(circle, curve[0], curve[order], c1, c2, t1, t2, EPS))
    {
      if(t1 >= -EPS && t1 < 1.0 - EPS)
      {
        circle_p.push_back(axom::utilities::isNearlyEqual(c1, 2.0 * M_PI, EPS) ? 0.0 : c1);
        curve_p.push_back(c_offset + c_scale * t1);
        foundIntersection = true;
      }

      if(t2 >= -EPS && t2 < 1.0 - EPS)
      {
        circle_p.push_back(axom::utilities::isNearlyEqual(c2, 2.0 * M_PI, EPS) ? 0.0 : c2);
        curve_p.push_back(c_offset + c_scale * t2);
        foundIntersection = true;
      }
    }
  }
  else
  {
    constexpr double splitVal = 0.5;
    constexpr double scaleFac = 0.5;

    BCurve c1(order);
    BCurve c2(order);
    curve.split(splitVal, c1, c2);
    c_scale *= scaleFac;

    // Note: we want to find all intersections, so don't short-circuit
    if(intersect_circle_bezier(circle, c1, circle_p, curve_p, sq_tol, EPS, order, c_offset, c_scale))
    {
      foundIntersection = true;
    }
    if(intersect_circle_bezier(circle, c2, circle_p, curve_p, sq_tol, EPS, order, c_offset + c_scale, c_scale))
    {
      foundIntersection = true;
    }
  }

  return foundIntersection;
}

template <typename T>
bool intersect_2d_circle_line(const Sphere<T, 2> &circ,
                              const Point<T, 2> &a,
                              const Point<T, 2> &b,
                              T &c1,
                              T &c2,
                              T &t1,
                              T &t2,
                              double EPS)
{
  T dx = b[0] - a[0];
  T dy = b[1] - a[1];
  T dr = std::sqrt(dx * dx + dy * dy);

  T D = (a[0] - circ.getCenter()[0]) * (b[1] - circ.getCenter()[1]) -
    (b[0] - circ.getCenter()[0]) * (a[1] - circ.getCenter()[1]);

  T disc = circ.getRadius() * circ.getRadius() * dr * dr - D * D;

  // Identify near-tangent cases by computing the distance from the line to the circle center,
  //  which should be equal to the radius of the circle if the line is tangent.
  if(axom::utilities::isNearlyEqual(std::abs(D) / dr, circ.getRadius(), EPS))
  {
    // Find the closest point on the line to the circle center
    Point<T, 2> cp(a);
    double num = (circ.getCenter()[0] - a[0]) * dx + (circ.getCenter()[1] - a[1]) * dy;

    cp[0] += num * dx / (dr * dr);
    cp[1] += num * dy / (dr * dr);

    c1 = std::atan2(cp[1] - circ.getCenter()[1], cp[0] - circ.getCenter()[0]);
    c1 = (c1 < 0.0) ? c1 + 2.0 * M_PI : c1;

    // Find the distance along the line to the closest point
    if(std::abs(dx) > std::abs(dy))
    {
      t1 = (cp[0] - a[0]) / dx;
    }
    else
    {
      t1 = (cp[1] - a[1]) / dy;
    }

    c2 = 0.0;
    t2 = -1.0;

    return true;
  }
  else if(disc < 0.0)
  {
    // No intersection
    return false;
  }

  disc = std::sqrt(disc);

  T x1 = (D * dy + (dy < 0 ? -1 : 1) * dx * disc) / (dr * dr);
  T x2 = (D * dy - (dy < 0 ? -1 : 1) * dx * disc) / (dr * dr);

  T y1 = (-D * dx + std::abs(dy) * disc) / (dr * dr);
  T y2 = (-D * dx - std::abs(dy) * disc) / (dr * dr);

  // Find the parameter values for the circle
  c1 = std::atan2(y1, x1);
  c1 = (c1 < 0.0) ? c1 + 2.0 * M_PI : c1;

  c2 = std::atan2(y2, x2);
  c2 = (c2 < 0.0) ? c2 + 2.0 * M_PI : c2;

  // Find the parameter values for the line
  if(std::abs(dx) > std::abs(dy))
  {
    t1 = (x1 - a[0] + circ.getCenter()[0]) / dx;
    t2 = (x2 - a[0] + circ.getCenter()[0]) / dx;
  }
  else
  {
    t1 = (y1 - a[1] + circ.getCenter()[1]) / dy;
    t2 = (y2 - a[1] + circ.getCenter()[1]) / dy;
  }

  return true;
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif  // AXOM_PRIMAL_INTERSECT_BEZIER_IMPL_HPP_
