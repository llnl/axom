// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file intersect.hpp
 *
 * \brief Consists of functions to test intersection among geometric primitives.
 */

#ifndef AXOM_PRIMAL_INTERSECT_HPP_
#define AXOM_PRIMAL_INTERSECT_HPP_

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/utilities/Utilities.hpp"

#include "axom/primal/constants.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"
#include "axom/primal/geometry/Plane.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/Ray.hpp"
#include "axom/primal/geometry/Line.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Sphere.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/geometry/NURBSCurve.hpp"
#include "axom/primal/geometry/NURBSPatch.hpp"

#include "axom/primal/operators/detail/intersect_impl.hpp"
#include "axom/primal/operators/detail/intersect_ray_impl.hpp"
#include "axom/primal/operators/detail/intersect_bounding_box_impl.hpp"
#include "axom/primal/operators/detail/intersect_bezier_impl.hpp"
#include "axom/primal/operators/detail/intersect_patch_impl.hpp"

namespace axom
{
namespace primal
{
/// \name Triangle Intersection Routines
/// \accelerated
/// @{

/*!
 * \brief Tests if 3D Triangles t1 and t2 intersect.
 * 
 * \param [in] t1 The first triangle
 * \param [in] t2 The second triangle
 * \param [in] includeBoundary Indicates if boundaries should be considered
 * when detecting intersections (default: false)
 * \param [in] EPS Tolerance for determining intersections (default: 1E-8)
 * \return status true iff t1 intersects with t2, otherwise, false.
 *
 * If parameter \a includeBoundary is false (default), this function will
 * return true if the interior of t1 intersects the interior of t2.  To include
 * triangle boundaries in intersections, specify \a includeBoundary as true.
 */
template <typename T>
AXOM_HOST_DEVICE bool intersect(const Triangle<T, 3>& t1,
                                const Triangle<T, 3>& t2,
                                bool includeBoundary = false,
                                double EPS = 1E-08)
{
  return detail::intersect_tri3D_tri3D<T>(t1, t2, includeBoundary, EPS);
}

/*!
 * \brief Tests if 2D Triangles t1 and t2 intersect.
 * \param [in] t1 The first triangle
 * \param [in] t2 The second triangle
 * \param [in] includeBoundary Indicates if boundaries should be considered
 * when detecting intersections (default: false)
 * \param [in] EPS Tolerance for determining intersections (default: 1E-8)
 * \return status true iff t1 intersects with t2, otherwise, false.
 *
 * If parameter \a includeBoundary is false (default), this function will
 * return true if the interior of t1 intersects the interior of t2.  To include
 * triangle boundaries in intersections, specify \a includeBoundary as true.
 */
template <typename T>
bool intersect(const Triangle<T, 2>& t1,
               const Triangle<T, 2>& t2,
               bool includeBoundary = false,
               double EPS = 1E-08)
{
  return detail::intersect_tri2D_tri2D<T>(t1, t2, includeBoundary, EPS);
}

/*!
 * \brief Determines if a triangle and a bounding box intersect
 * \param [in] tri user-supplied triangle (with three vertices).
 * \param [in] bb user-supplied axis aligned bounding box.
 * \return true iff tri intersects with bb, otherwise, false.
 */
template <typename T>
AXOM_HOST_DEVICE bool intersect(const Triangle<T, 3>& tri, const BoundingBox<T, 3>& bb)
{
  return detail::intersect_tri_bbox(tri, bb);
}

/*!
 * \brief Determines if a 3D triangle intersects a 3D ray.
 * \param [in] tri A 3D triangle
 * \param [in] ray A 3D ray
 * \return true iff tri intersects with ray, otherwise, false.
 */
template <typename T>
bool intersect(const Triangle<T, 3>& tri, const Ray<T, 3>& ray)
{
  T t = T();
  Point<double, 3> p;
  return detail::intersect_tri_ray(tri, ray, t, p);
}

/*!
 * \brief Determines if a 3D triangle intersects a 3D ray.
 * \param [in] tri A 3D triangle
 * \param [in] ray A 3D ray
 * \param [out] t Intersection point of tri and R, w.r.t. parametrization of R
 * \note If there is an intersection, the intersection point is:  R.at(t)
 * \return true iff tri intersects with ray, otherwise, false.
 */
template <typename T>
bool intersect(const Triangle<T, 3>& tri, const Ray<T, 3>& ray, T& t)
{
  Point<double, 3> p;
  return detail::intersect_tri_ray(tri, ray, t, p);
}

/*!
 * \brief Determines if a 3D triangle intersects a 3D ray.
 * \param [in] tri A 3D triangle
 * \param [in] ray A 3D ray
 * \param [out] t Intersection point of tri and R, w.r.t. parametrization of R
 * \param [out] p Intersection point of tri and R, in barycentric coordinates
 *   relative to tri.
 * \note If there is an intersection, the intersection point is:  R.at(t)
 * \return true iff tri intersects with ray, otherwise, false.
 * \note \a t and \a p only valid when function returns true
 */
template <typename T>
AXOM_HOST_DEVICE bool intersect(const Triangle<T, 3>& tri,
                                const Ray<T, 3>& ray,
                                T& t,
                                Point<double, 3>& p)
{
  bool retval = detail::intersect_tri_ray(tri, ray, t, p);

  if(retval)
  {
    // Add a small EPS to avoid dividing by zero
    double normalizer = p[0] + p[1] + p[2] + primal::PRIMAL_TINY;
    p.array() *= 1. / normalizer;
  }

  return retval;
}

/*!
 * \brief Determines if a 3D triangle intersects a 3D segment.
 * \param [in] tri A 3D triangle
 * \param [in] seg A 3D line segment
 * \return true iff tri intersects with seg, otherwise, false.
 */
template <typename T>
AXOM_HOST_DEVICE bool intersect(const Triangle<T, 3>& tri, const Segment<T, 3>& seg)
{
  T t = T();
  Point<double, 3> p;
  return detail::intersect_tri_segment(tri, seg, t, p);
}

/*!
 * \brief Determines if a 3D triangle intersects a 3D segment.
 * \param [in] tri A 3D triangle
 * \param [in] seg A 3D line segment
 * \param [out] t Intersection point of tri and seg, w.r.t. seg's
 *  parametrization
 * \return true iff tri intersects with seg, otherwise, false.
 */
template <typename T>
bool intersect(const Triangle<T, 3>& tri, const Segment<T, 3>& seg, T& t)
{
  Point<double, 3> p;
  return detail::intersect_tri_segment(tri, seg, t, p);
}

/*!
 * \brief Determines if a 3D triangle intersects a 3D segment.
 * \param [in] tri A 3D triangle
 * \param [in] seg A 3D line segment
 * \param [out] t Intersection point of tri and seg, w.r.t. seg's
 *  parametrization
 * \param [out] p Intersection point of tri and R, in barycentric coordinates
 *   relative to tri.
 * \note If there is an intersection, the intersection point pt is:
 *                     pt = seg.source() + t * ( seg.dest() - seg.target() )
 * \return true iff tri intersects with seg, otherwise, false.
 * \note \a t and \a p only valid when function returns true
 */
template <typename T>
bool intersect(const Triangle<T, 3>& tri, const Segment<T, 3>& seg, T& t, Point<double, 3>& p)
{
  bool retval = detail::intersect_tri_segment(tri, seg, t, p);

  if(retval)
  {
    // Add a small EPS to avoid dividing by zero
    double normalizer = p[0] + p[1] + p[2] + primal::PRIMAL_TINY;
    p.array() *= 1. / normalizer;
  }

  return retval;
}

/// @}

/// \name Ray Intersection Routines
/// @{

/*!
 * \brief Computes the intersection of the given ray, R, with the segment, S.
 *
 * \param [in] R the specified ray
 * \param [in] S the segment to check
 * \param [out] ray_param parametric coordinate of intersection along R, valid only if status=true.
 * \param [out] seg_param parametric coordinate of intersection along S, valid only if status=true.
 * \param [in] EPS tolerance for intersection tests
 *
 * \return status true iff R intersects with S, otherwise, false.
 *
 * \see primal::Ray
 * \see primal::Segment
 */
template <typename T>
AXOM_HOST_DEVICE bool intersect(const Ray<T, 2>& R,
                                const Segment<T, 2>& S,
                                T& ray_param,
                                T& seg_param,
                                const T EPS = 1e-8)
{
  return detail::intersect_ray(R, S, ray_param, seg_param, EPS);
}

/*!
 * \brief Computes the intersection of the given ray, R, with the segment, S.
 *
 * \param [in] R the specified ray
 * \param [in] S the segment to check
 * \param [out] ray_param parametric coordinate of intersection along R, valid only if status=true
 *
 * \note If you need to specify a tolerance for the intersection tests, please use the overload
 * of this function with two [OUT] parameters (\a ray_param and \a seg_param)
 *
 * \return status true iff R intersects with S, otherwise, false.
 *
 * \see primal::Ray
 * \see primal::Segment
 */
template <typename T>
bool intersect(const Ray<T, 2>& R, const Segment<T, 2>& S, T& ray_param)
{
  T seg_param;
  return intersect(R, S, ray_param, seg_param);
}

/*!
 * \brief Computes the intersection of the given ray, R, with the segment, S.
 *
 * \param [in] R the specified ray
 * \param [in] S the segment to check
 *
 * \note If you need to specify a tolerance for the intersection tests, please use the overload
 * of this function with two [OUT] parameters (\a ray_param and \a seg_param)
 *
 * \return status true iff R intersects with S, otherwise, false.
 *
 * \see primal::Ray
 * \see primal::Segment
 */
template <typename T>
AXOM_HOST_DEVICE bool intersect(const Ray<T, 2>& R, const Segment<T, 2>& S)
{
  T ray_param;
  T seg_param;
  return intersect(R, S, ray_param, seg_param);
}

/*!
 * \brief Computes the intersection of the given ray, R, with the segment, S.
 *
 * \param [in] R the specified ray
 * \param [in] S the segment to check
 * \param [out] ip the intersection point on S, valid only if status=true.
 * \param [in] EPS tolerance for intersection tests
 *
 * \return status true iff R intersects with S, otherwise, false.
 *
 * \see primal::Ray
 * \see primal::Segment
 * \see primal::Point
 */
template <typename T>
bool intersect(const Ray<T, 2>& R, const Segment<T, 2>& S, Point<T, 2>& ip, const T EPS = 1e-8)
{
  T ray_param;
  T seg_param;
  if(detail::intersect_ray(R, S, ray_param, seg_param, EPS))
  {
    ip = R.at(ray_param);
    return true;
  }
  return false;
}

/*!
 * \brief Computes the intersection of the given ray, R, with the Box, bb.
 *
 * \param [in] R the specified ray
 * \param [in] bb the user-supplied axis-aligned bounding box
 *
 * \param [out] ip the intersection point with minimum parameter value where R intersects bb.
 *
 * \return status true iff bb intersects with R, otherwise, false.
 *
 * \see primal::Ray
 * \see primal::Segment
 * \see primal::BoundingBox
 *
 * \note Computes Ray Box intersection using the slab method from pg 180 of
 *  Real Time Collision Detection by Christer Ericson.
 */
template <typename T, int DIM>
AXOM_HOST_DEVICE bool intersect(const Ray<T, DIM>& R, const BoundingBox<T, DIM>& bb, Point<T, DIM>& ip)
{
  return detail::intersect_ray(R, bb, ip);
}

/*!
 * \brief Computes the intersection of the given ray, R, with the Box, bb.
 *
 * \param [in] R the specified ray
 * \param [in] bb the user-supplied axis-aligned bounding box
 *
 * \return status true iff bb intersects with R, otherwise, false.
 *
 * \see primal::Ray
 * \see primal::Segment
 * \see primal::BoundingBox
 *
 * \note Computes Ray Box intersection using the slab method from pg 180 of
 *  Real Time Collision Detection by Christer Ericson.
 */
template <typename T, int DIM>
AXOM_HOST_DEVICE bool intersect(const Ray<T, DIM>& R, const BoundingBox<T, DIM>& bb)
{
  AXOM_STATIC_ASSERT(std::is_floating_point<T>::value);

  T tmin = axom::numerics::floating_point_limits<T>::min();
  T tmax = axom::numerics::floating_point_limits<T>::max();

  const T EPS = numerics::floating_point_limits<T>::epsilon();

  return detail::intersect_ray(R, bb, tmin, tmax, EPS);
}

/*!
 * \brief Computes the intersection of the given line, L, with the Box, bb.
 *
 * \param [in] L the specified line (two-sided ray)
 * \param [in] bb the user-supplied axis-aligned bounding box
 *
 * \param [out] ip the intersection point with minimum parameter value where L intersects bb.
 *
 * \return status true iff bb intersects with R, otherwise, false.
 *
 * \see primal::Line
 * \see primal::Segment
 * \see primal::BoundingBox
 *
 * \note Computes Ray Box intersection using the slab method from pg 180 of
 *  Real Time Collision Detection by Christer Ericson.
 */
template <typename T, int DIM>
AXOM_HOST_DEVICE bool intersect(const Line<T, DIM>& L, const BoundingBox<T, DIM>& bb, Point<T, DIM>& ip)
{
  return detail::intersect_line(L, bb, ip);
}

/*!
 * \brief Computes the intersection of the given line, L, with the Box, bb.
 *
 * \param [in] L the specified line (two-sided ray)
 * \param [in] bb the user-supplied axis-aligned bounding box
 *
 * \return status true iff bb intersects with R, otherwise, false.
 *
 * \see primal::Line
 * \see primal::BoundingBox
 *
 * \note Computes Ray Box intersection using the slab method from pg 180 of
 *  Real Time Collision Detection by Christer Ericson.
 */
template <typename T, int DIM>
AXOM_HOST_DEVICE bool intersect(const Line<T, DIM>& L, const BoundingBox<T, DIM>& bb)
{
  AXOM_STATIC_ASSERT(std::is_floating_point<T>::value);

  T tmin = -axom::numerics::floating_point_limits<T>::max();
  T tmax = axom::numerics::floating_point_limits<T>::max();

  const T EPS = numerics::floating_point_limits<T>::epsilon();

  return detail::intersect_ray(Ray<T, DIM>(L.origin(), L.direction()), bb, tmin, tmax, EPS);
}
/// @}

/// \name Segment-BoundingBox Intersection Routines
/// @{

/*!
 * \brief Computes the intersection of the given segment, \a S, with the Box, \a bb.
 *     If an intersection is found, output parameter \a ip contains an intersection point
 * \return status true iff \a bb intersects with \a S, otherwise, false.
 *
 * \note The intersection between segment \a S and box \a bb intersect, will, in general,
 * be a along a (1D) subset of segment \a S. One variant of this function returns the two
 * parametric coordinates of the intersections along \a S found while determining 
 * if there is a valid intersection. Another variant returns an intersection point along \a S
 * Specifically, it is the point of smallest parametric coordinate that is contained in \a bb
 * (i.e. with parameter \a tmin). These are only valid when the function returns true
 * 
 * Computes Segment-Box intersection using the slab method from pg 180 of
 * Real Time Collision Detection by Christer Ericson.
 */

/// This variant returns the two parametric coordinates of the intersection segment as OUT parameters
template <typename T, int DIM>
AXOM_HOST_DEVICE bool intersect(const Segment<T, DIM>& S,
                                const BoundingBox<T, DIM>& bb,
                                T& tmin,
                                T& tmax,
                                const double& EPS = 1e-8)
{
  const T segLength = S.length();
  tmin = static_cast<T>(0);
  tmax = static_cast<T>(segLength);

  bool intersects = segLength > 0. && detail::intersect_ray(Ray<T, DIM>(S), bb, tmin, tmax, EPS);

  // Scale parametric coordinates with respect to the segment
  if(intersects)
  {
    tmin /= segLength;
    tmax /= segLength;
  }

  return intersects;
}

/// This variant returns a point within the intersection as an OUT parameters
template <typename T, int DIM>
bool intersect(const Segment<T, DIM>& S,
               const BoundingBox<T, DIM>& bb,
               Point<T, DIM>& ip,
               const double& EPS = 1e-8)
{
  T tmin, tmax;
  if(intersect(S, bb, tmin, tmax, EPS))
  {
    ip = S.at(tmin);
    return true;
  }
  return false;
}

template <typename T, int DIM>
AXOM_HOST_DEVICE bool intersect(const Segment<T, DIM>& S, const BoundingBox<T, DIM>& bb)
{
  T tmin, tmax;
  return intersect(S, bb, tmin, tmax);
}

/// @}

/// \name Axis-Aligned Bounding Box Intersection Routines
/// @{

/*!
 * \brief Determines if two axis aligned bounding boxes intersect
 * \param [in] bb1 user-supplied axis aligned bounding box.
 * \param [in] bb2 user-supplied axis aligned bounding box.
 * \return true iff bb1 intersects with bb2, otherwise, false.
 */
template <typename T, int DIM>
AXOM_HOST_DEVICE bool intersect(const BoundingBox<T, DIM>& bb1, const BoundingBox<T, DIM>& bb2)
{
  return bb1.intersectsWith(bb2);
}

/// @}

/// \name Sphere Intersection Routines
/// @{

/*!
 * \brief Determines if two spheres intersect.
 *
 * \param [in] s1 user-supplied sphere object to check for intersection.
 * \param [in] s2 user-supplied sphere object to check for intersection.
 * \param [in] TOL tolerance used for intersection check (optional)
 *
 * \note If TOL is not supplied, the default is 1.e-9.
 *
 * \return status true iff s1 intersects with s2, otherwise, false.
 */
template <typename T, int DIM>
bool intersect(const Sphere<T, DIM>& s1, const Sphere<T, DIM>& s2, double TOL = 1.e-9)
{
  return s1.intersectsWith(s2, TOL);
}

/*!
 * \brief Determines if a 2D sphere intersects (overlaps with) a bounding box 
 *
 * \param [in] circle The sphere to check for intersection
 * \param [in] bb The bounding box to check for intersection
 */
template <typename T>
bool intersect(const Sphere<T, 2>& circle, const BoundingBox<T, 2>& bb)
{
  return detail::intersect_circle_bbox(circle, bb);
}

/*!
 * \brief Determines if a 2D sphere intersects a NURBS Curve
 * 
 * \param [in] circle The sphere to check for intersection
 * \param [in] curve The NURBS curve to check for intersection
 * \param [out] circle_params The parameter values of the sphere at the intersection points
 * \param [out] curve_params The parameter values of the curve at the intersection points
 * \param [in] tol Tolerance parameter for physical distances
 * \param [in] EPS Tolerance parameter for parameter-space distances
 * 
 * \return True if the sphere intersects the curve, false otherwise
 */
template <typename T>
bool intersect(const Sphere<T, 2>& circle,
               const NURBSCurve<T, 2>& curve,
               axom::Array<T>& circle_params,
               axom::Array<T>& curve_params,
               double tol = 1e-8,
               double EPS = 1e-8)
{
  const double sq_tol = tol * tol;
  bool foundIntersection = false;

  // Extract the Bezier curves of the NURBS curve
  auto beziers = curve.extractBezier();
  axom::Array<T> knot_vals = curve.getKnots().getUniqueKnots();

  // Check each Bezier segment for intersection
  for(int i = 0; i < beziers.size(); ++i)
  {
    axom::Array<T> temp_curve_p;
    axom::Array<T> temp_circle_p;
    detail::intersect_circle_bezier(circle,
                                    beziers[i],
                                    temp_circle_p,
                                    temp_curve_p,
                                    sq_tol,
                                    EPS,
                                    beziers[i].getOrder(),
                                    0.,
                                    1.);

    foundIntersection |= temp_curve_p.size() > 0;

    // Scale the intersection parameters back into the span of the NURBS curve
    for(int j = 0; j < temp_curve_p.size(); ++j)
    {
      curve_params.push_back(knot_vals[i] + temp_curve_p[j] * (knot_vals[i + 1] - knot_vals[i]));
      circle_params.push_back(temp_circle_p[j]);
    }
  }

  return foundIntersection;
}
/// @}

/// \name Oriented Bounding Box Intersection Routines
/// @{

template <typename T>
bool intersect(const OrientedBoundingBox<T, 1>& b1, const OrientedBoundingBox<T, 1>& b2)
{
  return detail::intersect_obb1D_obb1D(b1, b2);
}

/*!
 * \brief Determines if a 2D OBB intersects a 2D OBB.
 * \param [in] b1 A 2D OrientedBoundingBox
 * \param [in] b2 A 2D OrientedBoundingBox
 * \return true iff b1 intersects with b2, otherwise, false.
 */
template <typename T>
bool intersect(const OrientedBoundingBox<T, 2>& b1, const OrientedBoundingBox<T, 2>& b2)
{
  return detail::intersect_obb2D_obb2D(b1, b2);
}

/*!
 * \brief Determines if a 3D OBB intersects a 3D OBB.
 * \param [in] b1 A 3D OrientedBoundingBox
 * \param [in] b2 A 3D OrientedBoundingBox
 * \param [in] EPS error tolerance for intersection
 * \return true iff b1 intersects with b2, otherwise, false.
 */
template <typename T>
bool intersect(const OrientedBoundingBox<T, 3>& b1,
               const OrientedBoundingBox<T, 3>& b2,
               double EPS = 1E-4)
{
  return detail::intersect_obb3D_obb3D(b1, b2, EPS);
}

/// @}

/// \name Bezier Curve Intersection Routines
/// @{

/*!
 * \brief Tests if two Bezier Curves \a c1 and \a c2 intersect.
 * \return status true iff \a c1 intersects \a c2, otherwise false.
 *
 * \param [in] c1 the first BezierCurve, parametrized in [0,1)
 * \param [in] c2 the second BezierCurve, parametrized in [0,1)
 * \param [out] sp vector of parameter space intersection points for \a c1
 * \param [out] tp vector of parameter space intersection points for \a c2
 * \param [in] tol tolerance parameter for determining if a curve can
 * be approximated by a line segment.
 * \return True if the curves intersect, false otherwise. Intersection
 * parameters are stored in \a sp and \a tp
 *
 * Finds all intersection points between the two curves.
 *
 * \note This function assumes two dimensional curves in a plane.
 *
 * \note This function assumes that the curves are in general position.
 * Specifically, we assume that all intersections are at points and that
 * the curves don't overlap.
 *
 * \note This function assumes the all intersections have multiplicity
 * one, i.e. there are no points at which the curves and their derivatives
 * both intersect. Thus, the function does not find tangencies.
 *
 * \note This function assumes that the curves are half-open, i.e. they
 * contain their first endpoint, but not their last endpoint. Thus, the
 * curves do not intersect at \f$ s==1 \f$ or at \f$ t==1 \f$.
 */
template <typename T>
bool intersect(const BezierCurve<T, 2>& c1,
               const BezierCurve<T, 2>& c2,
               axom::Array<T>& sp,
               axom::Array<T>& tp,
               double tol = 1E-8)
{
  const double offset = 0.;
  const double scale = 1.;

  // for efficiency, linearity check actually uses a squared tolerance
  const double sq_tol = tol * tol;

  return detail::intersect_bezier_curves(c1,
                                         c2,
                                         sp,
                                         tp,
                                         sq_tol,
                                         c1.getOrder(),
                                         c2.getOrder(),
                                         offset,
                                         scale,
                                         offset,
                                         scale);
}

/*!
 * \brief Function to find intersections between a ray and a Bezier curve
 *
 * \param [in] r The input ray
 * \param [in] c The input curve
 * \param [out] rp Parametric coordinates of intersections in \a r [0, inf)
 * \param [out] cp Parametric coordinates of intersections in \a c [0, 1)
 * \param [in] tol Tolerance parameter for physical distances
 * \param [in] EPS Tolerance parameter for parameter-space distances
 * 
 * \note A BezierCurve is parametrized in [0,1). This function assumes the all
 *  intersections have multiplicity one, i.e. the function does not find tangencies.
 * 
 * \return True if the ray intersects the Bezier curve, False otherwise
 */
template <typename T>
bool intersect(const Ray<T, 2>& r,
               const BezierCurve<T, 2>& c,
               axom::Array<T>& rp,
               axom::Array<T>& cp,
               double tol = 1E-8,
               double EPS = 1E-8)
{
  const double offset = 0.;
  const double scale = 1.;

  // for efficiency, linearity check actually uses a squared tolerance
  const double sq_tol = tol * tol;
  const bool isHalfOpen = true;

  return detail::intersect_ray_bezier(r, c, rp, cp, sq_tol, EPS, c.getOrder(), offset, scale, isHalfOpen);
}

/*!
 * \brief Function to find intersections between a ray and a NURBS curve
 *
 * \param [in] r The input ray
 * \param [in] n The input curve
 * \param [out] rp Parametric coordinates of intersections in \a r [0, inf)
 * \param [out] cp Parametric coordinates of intersections in the knot span of \a n
 * \param [in] tol Tolerance parameter for physical distances
 * \param [in] EPS Tolerance parameter for parameter-space distances
 *
 * \note Assumes the NURBS curve is parameterized on a half-open interval [a, b),
 *  and assumes the all intersections have multiplicity one, i.e. the function does not find tangencies.
 * 
 * \return True if the ray intersects the NURBS curve, False otherwise
 */
template <typename T>
bool intersect(const Ray<T, 2>& r,
               const NURBSCurve<T, 2>& n,
               axom::Array<T>& rp,
               axom::Array<T>& np,
               double tol = 1E-8,
               double EPS = 1E-8)
{
  // Check a bounding box of the entire NURBS first
  Point<T, 2> ip;
  if(!intersect(r, n.boundingBox(), ip))
  {
    return false;
  }

  bool foundIntersection = false;

  // Decompose the NURBS curve into Bezier segments
  auto beziers = n.extractBezier();
  axom::Array<T> knot_vals = n.getKnots().getUniqueKnots();

  // Check each Bezier segment, and scale the intersection parameters
  //  back into the span of the original NURBS curve
  for(int i = 0; i < beziers.size(); ++i)
  {
    axom::Array<T> rc, nc;
    intersect(r, beziers[i], rc, nc, tol, EPS);

    foundIntersection |= !rc.empty();

    // Scale the intersection parameters back into the span of the NURBS curve
    for(int j = 0; j < rc.size(); ++j)
    {
      rp.push_back(rc[j]);
      np.push_back(axom::utilities::lerp(knot_vals[i], knot_vals[i + 1], nc[j]));
    }
  }

  return foundIntersection;
}
/// @}

/// \name Plane Intersection Routines
/// @{

/*!
 * \brief Determines if a 3D plane intersects a 3D bounding box.
 *        By default (checkOverlaps is false), checks if |s| <= r, 
 *        where "s" is the distance of the bounding box center to the plane,
 *        and "r" is the projected radius of the bounding box along the line
 *        parallel to the plane normal and going through the box center.
 *        If checkOverlaps is true, checks if |s| < r,
 *        where the bounding box overlaps both half spaces of the plane.
 * \param [in] p A 3D plane
 * \param [in] bb A 3D bounding box
 * \param [in] checkOverlaps If true, checks if bounding box overlaps both 
 *             halfspaces of the plane.
 *             Otherwise, overlap of both halfspaces is not guaranteed.
 *             Default is false.
 * \param [in] EPS tolerance parameter for determining if "s"
 *             is just within min/max of "r".
 * \return true iff plane intersects with bounding box, otherwise, false.
 *
 * \note Uses method from pg 164 of 
 *       Real Time Collision Detection by Christer Ericson.
 */
template <typename T>
AXOM_HOST_DEVICE bool intersect(const Plane<T, 3>& p,
                                const BoundingBox<T, 3>& bb,
                                bool checkOverlaps = false,
                                double EPS = 1E-08)
{
  return detail::intersect_plane_bbox(p, bb, checkOverlaps, EPS);
}

/*!
 * \brief Determines if a plane intersects a segment.
 * \param [in] plane A plane
 * \param [in] seg A line segment
 * \param [out] t Intersection point of plane and seg, w.r.t. seg's
 *  parametrization
 * \param [in] EPS tolerance parameter for determining if 0.0 <= t <= 1.0
 * \note If there is an intersection, the intersection point pt is:
 *                     pt = seg.at(t)
 * \return true iff plane intersects with seg, otherwise, false.
 * \note \a t is only valid when function returns true
 *
 * \note Uses method from pg 176 of 
 *       Real Time Collision Detection by Christer Ericson.
 */
template <typename T, int DIM>
AXOM_HOST_DEVICE bool intersect(const Plane<T, DIM>& plane,
                                const Segment<T, DIM>& seg,
                                T& t,
                                const double& EPS = 1e-12)
{
  return detail::intersect_plane_seg(plane, seg, t, EPS);
}

/*!
 * \brief Determines if a 3D plane intersects a tetrahedron.
 *
 * \param [in] p A 3D plane
 * \param [in] tet A 3D tetrahedron
 * \param [out] intersection A polygon containing the intersection.
 *
 * \return true if plane intersects with tetrahedron, otherwise, false.
 *
 * \note If no intersection is found, the output polygon will be empty.
 *       If the plane intersects at a tetrahedron vertex, the polygon
 *       will contain duplicated points.
 */
template <typename T>
AXOM_HOST_DEVICE bool intersect(const Plane<T, 3>& p,
                                const Tetrahedron<T, 3>& tet,
                                Polygon<T, 3>& intersection)
{
  return detail::intersect_plane_tet3d(p, tet, intersection);
}

/// @}

/*!
 * \brief Determines if a ray intersects a Bezier patch.
 * \param [in] ray The ray to intersect with the patch.
 * \param [in] patch The Bezier patch to intersect with the ray.
 * \param [out] t The t parameter(s) of intersection point(s).
 * \param [out] u The u parameter(s) of intersection point(s).
 * \param [out] v The v parameter(s) of intersection point(s).
 * \param [in] tol The tolerance for intersection (for physical distances).
 * \param [in] EPS The tolerance for intersection (for parameter distances).
 * \param [in] isHalfOpen True if the patch is parameterized in [0,1)^2.
 * \param [out] success False if an early return occured.
 * 
 * For bilinear patches, implements GARP algorithm from Chapter 8 of Ray Tracing Gems (2019)
 * For higher order patches, intersections are found through recursive subdivison
 *  until the subpatch is approximated by a bilinear patch.
 * Assumes that the ray is not tangent to the patch, and that the intersection
 *  is not at a point of degeneracy for which there are *infinitely* many intersections.
 * 
 * \note This method clears the input vectors t, u, v
 *  
 * \warning This function returns early if we record excessive intersections.
 *    This implies the patch is degenerate at the point of intersection.
 * 
 * \return True if an intersection is found, false otherwise
 */
template <typename T>
bool intersect(const Ray<T, 3>& ray,
               const BezierPatch<T, 3>& patch,
               axom::Array<T>& t,
               axom::Array<T>& u,
               axom::Array<T>& v,
               double tol,
               double EPS,
               bool isHalfOpen,
               bool& success)
{
  // Clear the input arrays
  t.clear();
  u.clear();
  v.clear();

  bool isRay = true;

  const int order_u = patch.getOrder_u();
  const int order_v = patch.getOrder_v();
  success = true;

  if(order_u < 1 || order_v < 1)
  {
    // Patch has no surface area, ergo no intersections
    return false;
  }
  else if(order_u == 1 && order_v == 1)
  {
    // Store the candidate intersections
    StaticArray<T, 2> tc, uc, vc;
    primal::Line<T, 3> line(ray.origin(), ray.direction());
    detail::intersect_line_bilinear_patch(line,
                                          patch(0, 0),
                                          patch(order_u, 0),
                                          patch(order_u, order_v),
                                          patch(0, order_v),
                                          tc,
                                          uc,
                                          vc,
                                          EPS,
                                          isRay);

    detail::select_candidates(tc, uc, vc, t, u, v, EPS, isHalfOpen);
  }
  else
  {
    // Store the candidate intersections
    axom::Array<T> tc, uc, vc;
    primal::Line<T, 3> line(ray.origin(), ray.direction());

    double u_offset = 0., v_offset = 0.;
    double u_scale = 1., v_scale = 1.;
    // For efficiency, linearity check actually uses a squared tolerance
    const double sq_tol = tol * tol;

    detail::intersect_line_patch(line,
                                 patch,
                                 tc,
                                 uc,
                                 vc,
                                 order_u,
                                 order_v,
                                 u_offset,
                                 u_scale,
                                 v_offset,
                                 v_scale,
                                 sq_tol,
                                 EPS,
                                 isRay,
                                 success);

    // Don't select candidates if we're in a failure state
    if(success)
    {
      detail::select_candidates(tc, uc, vc, t, u, v, EPS, isHalfOpen);
    }
    else
    {
      t = tc;
      u = uc;
      v = vc;
    }
  }

  return !t.empty();
}

///! \brief Determines if a ray intersects a Bezier patch with fewer exposed parameters.
template <typename T>
bool intersect(const Ray<T, 3>& ray,
               const BezierPatch<T, 3>& patch,
               axom::Array<T>& t,
               axom::Array<T>& u,
               axom::Array<T>& v,
               double tol = 1e-8,
               double EPS = 1e-8,
               bool isHalfOpen = false)
{
  bool success = true;
  return intersect(ray, patch, t, u, v, tol, EPS, isHalfOpen, success);
}

/*!
 * \brief Determines if a line (two-sided ray) intersects a Bezier patch.
 * \param [in] line The line to intersect with the patch.
 * \param [in] patch The Bezier patch to intersect with the ray.
 * \param [out] t The t parameter(s) of intersection point(s).
 * \param [out] u The u parameter(s) of intersection point(s).
 * \param [out] v The v parameter(s) of intersection point(s).
 * \param [in] tol The tolerance for intersection (for physical distances).
 * \param [in] EPS The tolerance for intersection (for parameter distances).
 * \param [in] isHalfOpen True if the patch is parameterized in [0,1)^2.
 * \param [out] success False if an early return occured.
 * 
 * For bilinear patches, implements GARP algorithm from Chapter 8 of Ray Tracing Gems (2019)
 * For higher order patches, intersections are found through recursive subdivison
 *  until the subpatch is approximated by a bilinear patch.
 * Assumes that the line is not tangent to the patch, and that the intersection
 *  is not at a point of degeneracy for which there are *infinitely* many intersections.
 * 
 * \note This method clears the input vectors t, u, v
 *  
 * \warning This function returns early if we record excessive intersections.
 *    This implies the patch is degenerate at the point of intersection.
 * 
 * \return True if an intersection is found, false otherwise
 */
template <typename T>
bool intersect(const Line<T, 3>& line,
               const BezierPatch<T, 3>& patch,
               axom::Array<T>& t,
               axom::Array<T>& u,
               axom::Array<T>& v,
               double tol,
               double EPS,
               bool isHalfOpen,
               bool& success)
{
  // Clear the input arrays
  t.clear();
  u.clear();
  v.clear();

  bool isRay = false;

  const int order_u = patch.getOrder_u();
  const int order_v = patch.getOrder_v();
  success = true;

  if(order_u < 1 || order_v < 1)
  {
    // Patch has no surface area, ergo no intersections
    return false;
  }
  else if(order_u == 1 && order_v == 1)
  {
    // Store the candidate intersections
    StaticArray<T, 2> tc, uc, vc;
    detail::intersect_line_bilinear_patch(line,
                                          patch(0, 0),
                                          patch(order_u, 0),
                                          patch(order_u, order_v),
                                          patch(0, order_v),
                                          tc,
                                          uc,
                                          vc,
                                          EPS,
                                          isRay);

    detail::select_candidates(tc, uc, vc, t, u, v, EPS, isHalfOpen);
  }
  else
  {
    // Store the candidate intersections
    axom::Array<T> tc, uc, vc;

    double u_offset = 0., v_offset = 0.;
    double u_scale = 1., v_scale = 1.;
    // For efficiency, linearity check actually uses a squared tolerance
    const double sq_tol = tol * tol;

    detail::intersect_line_patch(line,
                                 patch,
                                 tc,
                                 uc,
                                 vc,
                                 order_u,
                                 order_v,
                                 u_offset,
                                 u_scale,
                                 v_offset,
                                 v_scale,
                                 sq_tol,
                                 EPS,
                                 isRay,
                                 success);

    // Don't select candidates if we're in a failure state
    if(success)
    {
      detail::select_candidates(tc, uc, vc, t, u, v, EPS, isHalfOpen);
    }
    else
    {
      t = tc;
      u = uc;
      v = vc;
    }
  }

  return !t.empty();
}

///! \brief Determines if a line intersects a Bezier patch with fewer exposed parameters.
template <typename T>
bool intersect(const Line<T, 3>& line,
               const BezierPatch<T, 3>& patch,
               axom::Array<T>& t,
               axom::Array<T>& u,
               axom::Array<T>& v,
               double tol = 1e-8,
               double EPS = 1e-8,
               bool isHalfOpen = false)
{
  bool success = true;
  return intersect(line, patch, t, u, v, tol, EPS, isHalfOpen, success);
}

/*! 
 * \brief Determines if a ray intersects a NURBS patch.
 * \param [in] ray The ray to intersect with the patch.
 * \param [in] patch The NURBS patch to intersect with the ray.
 * \param [out] t The t parameter(s) of intersection point(s).
 * \param [out] u The u parameter(s) of intersection point(s).
 * \param [out] v The v parameter(s) of intersection point(s).
 * \param [in] tol The tolerance for intersection (for physical distances).
 * \param [in] EPS The tolerance for intersection (for parameter distances).
 * \param [in] countUntrimmed True if intersections with the untrimmed patch should also be recorded.
 * \param [in] isHalfOpen True if the patch is parameterized in [0,1)^2.
 * \param [out] success False if an early return occured.
 * 
 * Perform Bezier extraction and record intersections with each patch.
 * After intersections are recorded, parameter points located outside the trimming
 *  curves are pruned from the list (unless specified by `countUntrimmed`).
 * 
 * \note This method clears the input vectors t, u, v
 *  
 * \warning This function returns early if we record excessive intersections.
 *    This implies the patch is degenerate at the point of intersection.
 * 
 * \return True if an intersection is found, false otherwise
 */
template <typename T>
bool intersect(const Ray<T, 3>& ray,
               const NURBSPatch<T, 3>& patch,
               axom::Array<T>& t,
               axom::Array<T>& u,
               axom::Array<T>& v,
               double tol,
               double EPS,
               bool countUntrimmed,
               bool isHalfOpen,
               bool& success)
{
  // Clear the input arrays
  t.clear();
  u.clear();
  v.clear();

  // Check a bounding box of the entire NURBS first
  Point<T, 3> ip;
  if(!intersect(ray, patch.boundingBox(), ip))
  {
    return false;
  }

  // Decompose the NURBS patch into Bezier patches
  auto beziers = patch.extractBezier();

  axom::Array<T> knot_vals_u = patch.getKnots_u().getUniqueKnots();
  axom::Array<T> knot_vals_v = patch.getKnots_v().getUniqueKnots();

  const auto num_knot_span_u = knot_vals_u.size() - 1;
  const auto num_knot_span_v = knot_vals_v.size() - 1;

  // Store candidate intersections
  axom::Array<T> tc, uc, vc;

  bool isBezierHalfOpen = false;
  success = true;

  // Check each Bezier patch, and scale the intersection parameters
  //   back into the span of the original NURBS patch
  // If we're already in a failure state, skip the remaining patches
  for(int i = 0; i < num_knot_span_u; ++i)
  {
    for(int j = 0; success && j < num_knot_span_v; ++j)
    {
      auto& bezier = beziers[i * num_knot_span_v + j];

      // Store candidate intersections from each Bezier patch
      axom::Array<T> tcc, ucc, vcc;

      intersect(ray, bezier, tcc, ucc, vcc, tol, EPS, isBezierHalfOpen, success);

      // Scale the intersection parameters back into the span of the NURBS patch
      for(int k = 0; k < tcc.size(); ++k)
      {
        tc.push_back(tcc[k]);
        uc.push_back(axom::utilities::lerp(knot_vals_u[i], knot_vals_u[i + 1], ucc[k]));
        vc.push_back(axom::utilities::lerp(knot_vals_v[j], knot_vals_v[j + 1], vcc[k]));
      }
    }
  }

  // Skip second pass if in a failure state
  if(!success)
  {
    t = tc;
    u = uc;
    v = vc;

    return !t.empty();
  }

  // Do a second pass to remove duplicates from uc, vc
  const double sq_EPS = EPS * EPS;

  // The number of reported intersection points will be small,
  //  so we don't need to fully sort the list

  double max_u_knot = patch.getKnots_u()[patch.getKnots_u().getNumKnots() - 1];
  double max_v_knot = patch.getKnots_v()[patch.getKnots_v().getNumKnots() - 1];

  // Don't de-duplicate if we're in a failure state
  for(int i = 0; i < tc.size(); ++i)
  {
    // Also remove any intersections on the half-interval boundaries
    if(isHalfOpen && (uc[i] >= max_u_knot - EPS || vc[i] >= max_v_knot - EPS))
    {
      continue;
    }

    // Also remove any intersections that are trimmed out
    if(!countUntrimmed && !patch.isVisible(uc[i], vc[i]))
    {
      continue;
    }

    Point<T, 2> uv({uc[i], vc[i]});

    bool foundDuplicate = false;
    for(int j = i + 1; !foundDuplicate && j < tc.size(); ++j)
    {
      if(squared_distance(uv, Point<T, 2>({uc[j], vc[j]})) < sq_EPS)
      {
        foundDuplicate = true;
      }
    }

    if(!foundDuplicate)
    {
      t.push_back(tc[i]);
      u.push_back(uc[i]);
      v.push_back(vc[i]);
    }
  }

  return !t.empty();
}

///! \brief Determines if a ray intersects a NURBS patch with fewer exposed parameters.
template <typename T>
bool intersect(const Ray<T, 3>& ray,
               const NURBSPatch<T, 3>& patch,
               axom::Array<T>& t,
               axom::Array<T>& u,
               axom::Array<T>& v,
               double tol = 1e-8,
               double EPS = 1e-8,
               bool countUntrimmed = true,
               bool isHalfOpen = false)
{
  bool success = true;
  return intersect(ray, patch, t, u, v, tol, EPS, countUntrimmed, isHalfOpen, success);
}

/*! 
 * \brief Determines if a line (two-sided ray) intersects a NURBS patch.
 * \param [in] line The line to intersect with the patch.
 * \param [in] patch The NURBS patch to intersect with the line.
 * \param [out] t The t parameter(s) of intersection point(s).
 * \param [out] u The u parameter(s) of intersection point(s).
 * \param [out] v The v parameter(s) of intersection point(s).
 * \param [in] tol The tolerance for intersection (for physical distances).
 * \param [in] EPS The tolerance for intersection (for parameter distances).
 * \param [in] countUntrimmed True if intersections with the untrimmed patch should also be recorded.
 * \param [in] isHalfOpen True if the patch is parameterized in [0,1)^2.
 * \param [out] success False if an early return occured.
 * 
 * Perform Bezier extraction and record intersections with each patch.
 * After intersections are recorded, parameter points located outside the trimming
 *  curves are pruned from the list (unless specified by `countUntrimmed`).
 *  
 * \note This method clears the input vectors t, u, v
 *  
 * \warning This function returns early if we record excessive intersections.
 *    This implies the patch is degenerate at the point of intersection.
 * 
 * \return True if an intersection is found, false otherwise
 */
template <typename T>
bool intersect(const Line<T, 3>& line,
               const NURBSPatch<T, 3>& patch,
               axom::Array<T>& t,
               axom::Array<T>& u,
               axom::Array<T>& v,
               double tol,
               double EPS,
               bool countUntrimmed,
               bool isHalfOpen,
               bool& success)
{
  // Clear the input arrays
  t.clear();
  u.clear();
  v.clear();

  // Check a bounding box of the entire NURBS first
  Point<T, 3> ip;
  if(!intersect(line, patch.boundingBox().expand(10 * tol), ip))
  {
    return false;
  }

  // Decompose the NURBS patch into Bezier patches
  auto beziers = patch.extractBezier();

  axom::Array<T> knot_vals_u = patch.getKnots_u().getUniqueKnots();
  axom::Array<T> knot_vals_v = patch.getKnots_v().getUniqueKnots();

  const auto num_knot_span_u = knot_vals_u.size() - 1;
  const auto num_knot_span_v = knot_vals_v.size() - 1;

  // Store candidate intersections
  axom::Array<T> tc, uc, vc;

  bool isBezierHalfOpen = false;
  success = true;

  // Check each Bezier patch, and scale the intersection parameters
  //  back into the span of the original NURBS patch
  // If we're already in a failure state, skip the remaining patches
  for(int i = 0; i < num_knot_span_u; ++i)
  {
    for(int j = 0; success && j < num_knot_span_v; ++j)
    {
      auto& bezier = beziers[i * num_knot_span_v + j];

      // Store candidate intersections from each Bezier patch
      axom::Array<T> tcc, ucc, vcc;

      // If we're already in a failure state, skip the remaining patches
      intersect(line, bezier, tcc, ucc, vcc, tol, EPS, isBezierHalfOpen, success);

      // Scale the intersection parameters back into the span of the NURBS patch
      for(int k = 0; k < tcc.size(); ++k)
      {
        tc.push_back(tcc[k]);
        uc.push_back(axom::utilities::lerp(knot_vals_u[i], knot_vals_u[i + 1], ucc[k]));
        vc.push_back(axom::utilities::lerp(knot_vals_v[j], knot_vals_v[j + 1], vcc[k]));
      }
    }
  }

  // Skip second pass if in a failure state
  if(!success)
  {
    t = tc;
    u = uc;
    v = vc;

    return !t.empty();
  }

  // Do a second pass to remove duplicates from uc, vc
  const double sq_EPS = EPS * EPS;

  // The number of reported intersection points will be small,
  //  so we don't need to fully sort the list

  double max_u_knot = patch.getKnots_u()[patch.getKnots_u().getNumKnots() - 1];
  double max_v_knot = patch.getKnots_v()[patch.getKnots_v().getNumKnots() - 1];

  // Don't de-duplicate if we're in a failure state
  for(int i = 0; i < tc.size(); ++i)
  {
    // Also remove any intersections on the half-interval boundaries
    if(isHalfOpen && (uc[i] >= max_u_knot - EPS || vc[i] >= max_v_knot - EPS))
    {
      continue;
    }

    // Also remove any intersections that are trimmed out
    if(!countUntrimmed && !patch.isVisible(uc[i], vc[i]))
    {
      continue;
    }

    Point<T, 2> uv({uc[i], vc[i]});

    bool foundDuplicate = false;
    for(int j = i + 1; !foundDuplicate && j < tc.size(); ++j)
    {
      if(squared_distance(uv, Point<T, 2>({uc[j], vc[j]})) < sq_EPS)
      {
        foundDuplicate = true;
      }
    }

    if(!foundDuplicate)
    {
      t.push_back(tc[i]);
      u.push_back(uc[i]);
      v.push_back(vc[i]);
    }
  }

  return !t.empty();
}

///! \brief Determines if a line intersects a NURBS patch with fewer exposed parameters.
template <typename T>
bool intersect(const Line<T, 3>& line,
               const NURBSPatch<T, 3>& patch,
               axom::Array<T>& t,
               axom::Array<T>& u,
               axom::Array<T>& v,
               double tol = 1e-8,
               double EPS = 1e-8,
               bool countUntrimmed = true,
               bool isHalfOpen = false)
{
  bool success = true;
  return intersect(line, patch, t, u, v, tol, EPS, countUntrimmed, isHalfOpen, success);
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_INTERSECT_HPP_
