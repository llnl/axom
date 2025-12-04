// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file winding_number.hpp
 *
 * \brief Consists of methods to compute the generalized winding number (GWN) 
 *        for points with respect to various geometric objects.
 */

#ifndef AXOM_PRIMAL_WINDING_NUMBER_HPP_
#define AXOM_PRIMAL_WINDING_NUMBER_HPP_

// Axom includes
#include "axom/core.hpp"
#include "axom/config.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/Polyhedron.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/NURBSCurve.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/geometry/NURBSPatch.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"

#include "axom/primal/operators/detail/winding_number_2d_impl.hpp"
#include "axom/primal/operators/detail/winding_number_2d_memoization.hpp"

#include "axom/primal/operators/detail/winding_number_3d_impl.hpp"
#include "axom/primal/operators/detail/winding_number_3d_memoization.hpp"

// C++ includes
#include <cmath>

namespace axom
{
namespace primal
{
/*!
 * \brief Compute the GWN for a 2D point wrt a 2D line segment
 *
 * \param [in] q The query point to test
 * \param [in] s The line segment
 * \param [in] edge_tol The tolerance at which a point is on the line
 *
 * \return The GWN
 */
template <typename T>
double winding_number(const Point<T, 2>& q, const Segment<T, 2>& s, double edge_tol = 1e-8)
{
  bool dummy_isOnEdge = false;
  return detail::linear_winding_number(q, s[0], s[1], dummy_isOnEdge, edge_tol);
}

/*!
 * \brief Compute the winding number for a 2D point wrt a 2D triangle
 *
 * \param [in] q The query point to test
 * \param [in] tri The triangle
 * \param [in] includeBoundary If true, points on the boundary are considered interior.
 * \param [in] edge_tol The tolerance at which a point is on the line
 *
 * The triangle is assumed to be closed, so the winding number is an integer
 * 
 * \return The integer winding number
 */
template <typename T>
int winding_number(const Point<T, 2>& q,
                   const Triangle<T, 2>& tri,
                   bool includeBoundary = false,
                   double edge_tol = 1e-8)
{
  return winding_number(q,
                        Polygon<T, 2>(axom::Array<Point<T, 2>>({tri[0], tri[1], tri[2]})),
                        includeBoundary,
                        edge_tol);
}

/*!
 * \brief Computes the winding number for a 2D point wrt a 2D polygon
 *
 * \param [in] R The query point to test
 * \param [in] P The Polygon object to test for containment
 * \param [out] isOnEdge An optional return parameter if the point is on the boundary
 * \param [in] includeBoundary If true, points on the boundary are considered interior
 * \param [in] edge_tol The distance at which a point is considered on the boundary
 * 
 * \return The integer winding number
 */
template <typename T>
int winding_number(const Point<T, 2>& R,
                   const Polygon<T, 2>& P,
                   bool& isOnEdge,
                   bool includeBoundary = false,
                   double edge_tol = 1e-8)
{
  return detail::polygon_winding_number(R, P, isOnEdge, includeBoundary, edge_tol);
}

/*!
 * \brief Computes the winding number for a 2D point wrt a 2D polygon
 *
 * \param [in] R The query point to test
 * \param [in] P The Polygon object to test for containment
 * \param [in] includeBoundary If true, points on the boundary are considered interior
 * \param [in] edge_tol The distance at which a point is considered on the boundary
 * 
 * Computes the integer winding number for a polygon without an additional
 *  return parameter for whether the point is on the boundary.
 * 
 * \return The integer winding number
 */
template <typename T>
int winding_number(const Point<T, 2>& R,
                   const Polygon<T, 2>& P,
                   bool includeBoundary = false,
                   double edge_tol = 1e-8)
{
  bool isOnEdge = false;
  return detail::polygon_winding_number(R, P, isOnEdge, includeBoundary, edge_tol);
}

/*!
 * \brief Computes the GWN for a 2D point wrt a 2D NURBS curve
 *
 * \param [in] q The query point to test
 * \param [in] n The NURBS curve object 
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const NURBSCurve<T, 2>& n,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  bool dummy_isOnCurve = false;
  return detail::nurbs_winding_number(q, n, dummy_isOnCurve, edge_tol, EPS);
}

/*!
 * \brief Computes the GWN for a 2D point wrt a 2D NURBS curve
 *
 * \param [in] q The query point to test
 * \param [in] bezier The Bezier curve object 
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const BezierCurve<T, 2>& bezier,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  bool dummy_isOnCurve = false;
  return detail::bezier_winding_number(q, bezier, dummy_isOnCurve, edge_tol, EPS);
}

/*!
 * \brief Computes the GWN for a 2D point wrt to a 2D curved polygon
 *
 * \tparam CurveType The BezierCurve, NURBSCurve, or NURBSCurveGWNCache which represents the curve
 * \param [in] q The query point to test
 * \param [in] cpoly The CurvedPolygon object
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the GWN for the curved polygon by summing the GWN for each curved edge
 * 
 * \return The GWN.
 */
template <typename T, typename CurveType>
double winding_number(const Point<T, 2>& q,
                      const CurvedPolygon<CurveType>& cpoly,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  double ret_val = 0.0;
  for(int i = 0; i < cpoly.numEdges(); i++)
  {
    ret_val += winding_number(q, cpoly[i], edge_tol, EPS);
  }

  return ret_val;
}

/*!
 * \brief Computes the GWN for a 2D point wrt to a collection of 2D Bezier curves
 *
 * \param [in] q The query point to test
 * \param [in] carray The array of Bezier curves
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Sums the GWN at `q` for each curved edge
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const axom::Array<BezierCurve<T, 2>>& carray,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  bool dummy_isOnCurve = false;
  double ret_val = 0.0;
  for(int i = 0; i < carray.size(); i++)
  {
    ret_val += detail::bezier_winding_number(q, carray[i], dummy_isOnCurve, edge_tol, EPS);
  }

  return ret_val;
}

/*!
 * \brief Computes the GWN for a 2D point wrt to a collection of 2D NURBS curves
 *
 * \param [in] q The query point to test
 * \param [in] narray The array of NURBS curves
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Sums the GWN at `q` for each curved edge
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const axom::Array<NURBSCurve<T, 2>>& narray,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  bool dummy_isOnCurve = false;
  double ret_val = 0.0;
  for(int i = 0; i < narray.size(); i++)
  {
    ret_val += detail::nurbs_winding_number(q, narray[i], dummy_isOnCurve, edge_tol, EPS);
  }

  return ret_val;
}

/*!
 * \brief Computes the GWN for a 2D point wrt memoized data for a 2D NURBS curve 
 *
 * \param [in] query The query point to test
 * \param [in] nurbs_cache The NURBS curve cache data object containing memoized values
 * \param [out] isOnCurve Set to true is the query point is on the curve
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 2>& query,
                      const detail::NURBSCurveGWNCache<T>& nurbs_cache,
                      bool& isOnCurve,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  // Early return is possible for most points + curves
  if(!nurbs_cache.boundingBox().contains(query))
  {
    return detail::linear_winding_number(query,
                                         nurbs_cache.getInitPoint(),
                                         nurbs_cache.getEndPoint(),
                                         isOnCurve,
                                         edge_tol);
  }

  double gwn = 0.0;

  bool isOnThisCurve = false;
  isOnCurve = false;

  for(int n = 0; n < nurbs_cache.getNumKnotSpans(); ++n)
  {
    gwn +=
      detail::bezier_winding_number_memoized(query, nurbs_cache, n, 0, 0, isOnThisCurve, edge_tol, EPS);
    isOnCurve = isOnCurve || isOnThisCurve;
  }

  return gwn;
}

//! \brief Overload without optional return parameter
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const detail::NURBSCurveGWNCache<T>& nurbs_cache,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  bool dummy_isOnCurve;
  return winding_number(q, nurbs_cache, dummy_isOnCurve, edge_tol, EPS);
}

/*!
 * \brief Computes the GWN for a 2D point wrt an array of memoized data for 2D NURBS curves
 *
 * \param [in] query The query point to test
 * \param [in] nurbs_curve_arr The array of memoized curve objects
 * \param [out] isOnCurve Set to true is the query point is on the curve
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 2>& query,
                      const axom::Array<detail::NURBSCurveGWNCache<T>>& nurbs_curve_arr,
                      bool& isOnCurve,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  double gwn = 0;
  isOnCurve = false;
  for(int i = 0; i < nurbs_curve_arr.size(); ++i)
  {
    bool isOnThisCurve = false;
    gwn += winding_number(query, nurbs_curve_arr[i], isOnThisCurve, edge_tol, EPS);
    isOnCurve = isOnCurve || isOnThisCurve;
  }

  return gwn;
}

//! \brief Overload without optional return parameter
template <typename T>
double winding_number(const Point<T, 2>& query,
                      const axom::Array<detail::NURBSCurveGWNCache<T>>& nurbs_curve_arr,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  bool dummy_isOnThisCurve = false;
  return winding_number(query, nurbs_curve_arr, dummy_isOnThisCurve, edge_tol, EPS);
}

/*!
 * \brief Computes the GWN for an array of 2D points wrt an array of cached data for 2D NURBS curves
 *
 * \param [in] query_arr The array of query points to test
 * \param [in] nurbs_curve_arr The array of memoized curve objects
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * \note This method is accelerated via memoization, i.e. dynamically caching and reusing intermediate
 *   values for each curve across query points
 * 
 * \return The array of GWN values.
 */
template <typename T>
axom::Array<double> winding_number(const axom::Array<Point<T, 2>>& query_arr,
                                   const axom::Array<detail::NURBSCurveGWNCache<T>>& nurbs_curve_arr,
                                   double edge_tol = 1e-8,
                                   double EPS = 1e-8)
{
  bool dummy_isOnCurve;
  axom::Array<double> ret_val(query_arr.size());
  for(int n = 0; n < query_arr.size(); ++n)
  {
    ret_val[n] = 0.0;

    for(int i = 0; i < nurbs_curve_arr.size(); ++i)
    {
      ret_val[n] += detail::bezier_winding_number_memoized(query_arr[n],
                                                           nurbs_curve_arr[i],
                                                           dummy_isOnCurve,
                                                           edge_tol,
                                                           EPS);
    }
  }

  return ret_val;
}

/*!
 * \brief Computes the GWN for an array of 2D points wrt an array of generic 2D curves
 *
 * \tparam CurveType The BezierCurve, NURBSCurve, or NURBSCurveGWNCache which represents the curve
 * \param [in] query_arr The array of query point to test
 * \param [in] curve_arr The array of curve objects
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * \note This method is accelerated via memoization, i.e. dynamically caching and reusing intermediate
 *   values for each curve across query points
 *
 * \return The array of GWN values.
 */
template <typename T, typename CurveType>
axom::Array<double> winding_number(const axom::Array<Point<T, 2>>& query_arr,
                                   const axom::Array<CurveType>& curve_arr,
                                   double edge_tol = 1e-8,
                                   double EPS = 1e-8)
{
  axom::Array<detail::NURBSCurveGWNCache<T>> cache_arr(0, curve_arr.size());

  for(int i = 0; i < curve_arr.size(); ++i)
  {
    cache_arr.emplace_back(detail::NURBSCurveGWNCache<T>(curve_arr[i], edge_tol));
  }

  return winding_number(query_arr, cache_arr, edge_tol, EPS);
}

/*!
 * \brief Computes the GWN for an array of 2D points wrt to a 2D curved polygon
 *
 * \tparam CurveType The BezierCurve, NURBSCurve, or NURBSCurveGWNCache which represents the curve
 * \param [in] q_arr The array of query points to test
 * \param [in] cpoly The CurvedPolygon object of generic curves
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the GWN for the curved polygon by summing the GWN for each curved edge
 * 
 * \note This method is accelerated via memoization, i.e. dynamically caching and reusing intermediate
 *   values for each curve across query points
 * 
 * \return The GWN.
 */
template <typename T, typename CurveType>
axom::Array<double> winding_number(const axom::Array<Point<T, 2>>& q_arr,
                                   const CurvedPolygon<CurveType>& cpoly,
                                   double edge_tol = 1e-8,
                                   double EPS = 1e-8)
{
  axom::Array<detail::NURBSCurveGWNCache<T>> cache_arr(0, cpoly.numEdges());

  for(int i = 0; i < cpoly.numEdges(); ++i)
  {
    cache_arr.emplace_back(detail::NURBSCurveGWNCache<T>(cpoly[i], edge_tol));
  }

  axom::Array<double> ret_val(q_arr.size());
  for(int n = 0; n < q_arr.size(); ++n)
  {
    ret_val[n] = winding_number(q_arr[n], cache_arr, edge_tol, EPS);
  }

  return ret_val;
}

/*!
 * \brief Computes the GWN for an array of 2D points wrt to a 2D curved polygon
 *
 * \param [in] q_arr The array of query points to test
 * \param [in] cpoly The CurvedPolygon object of NURBS curves with cached GWN data
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the GWN for the curved polygon by summing the GWN for each curved edge
 * 
 * \warning Because the cache isKdiscarded immediately after computation,
 *  this method is not accelerated by memoization
 * 
 * \return The GWN.
 */
template <typename T>
axom::Array<double> winding_number(const axom::Array<Point<T, 2>>& q_arr,
                                   const CurvedPolygon<detail::NURBSCurveGWNCache<T>>& cpoly,
                                   double edge_tol = 1e-8,
                                   double EPS = 1e-8)
{
  axom::Array<double> ret_val(q_arr.size());
  for(int n = 0; n < q_arr.size(); ++n)
  {
    ret_val[n] = winding_number(q_arr[n], cpoly, edge_tol, EPS);
  }

  return ret_val;
}

///@{
//! @name Winding number operations between 3D points and primitives

/*!
 * \brief Computes the GWN for a 3D point wrt a 3D triangle
 *
 * \param [in] q The query point to test
 * \param [in] tri The 3D Triangle object
 * \param [out] isOnFace An optional return parameter if the point is on the triangle
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the GWN as the solid angle modulo 4pi using the formula from 
 *  Oosterom, Strackee, "The Solid Angle of a Plane Triangle" 
 *  IEEE Transactions on Biomedical Engineering, Vol BME-30, No. 2, February 1983
 * with extra adjustments if the triangle takes up a full octant
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 3>& q,
                      const Triangle<T, 3>& tri,
                      bool& isOnFace,
                      const double edge_tol = 1e-8,
                      const double EPS = 1e-8)
{
  constexpr double gwn_modulo = 0.5 * M_1_PI;

  using Vec3 = Vector<T, 3>;

  if(tri.area() == 0)
  {
    return 0;
  }

  const Vec3 a = tri[0] - q;
  const Vec3 b = tri[1] - q;
  const Vec3 c = tri[2] - q;

  // Compute norms. Possibly return early
  const double a_norm = a.norm();
  const double b_norm = b.norm();
  const double c_norm = c.norm();

  if(a_norm < edge_tol || b_norm < edge_tol || c_norm < edge_tol)
  {
    return 0;
  }

  const double num = Vec3::scalar_triple_product(a, b, c);
  if(axom::utilities::isNearlyEqual(num, 0.0, EPS))
  {
    isOnFace = true;
    return 0;
  }

  const double denom =
    a_norm * b_norm * c_norm + a_norm * b.dot(c) + b_norm * a.dot(c) + c_norm * a.dot(b);

  // Handle direct cases where argument to atan is undefined
  if(axom::utilities::isNearlyEqual(denom, 0.0, EPS))
  {
    return (num > 0) ? 0.25 : -0.25;
  }

  // Note: denom==0 and num==0 handled above
  if(denom > 0)
  {
    return gwn_modulo * atan(num / denom);
  }
  else
  {
    return (num > 0) ? gwn_modulo * atan(num / denom) + 0.5 : gwn_modulo * atan(num / denom) - 0.5;
  }
}

/*!
 * \brief Computes the GWN for a 3D point wrt a 3D triangle
 *
 * \param [in] q The query point to test
 * \param [in] tri The 3D Triangle object
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the GWN for the triangle without an additional return parameter
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 3>& q,
                      const Triangle<T, 3>& tri,
                      const double edge_tol = 1e-8,
                      const double EPS = 1e-8)
{
  bool isOnFace = false;
  return winding_number(q, tri, isOnFace, edge_tol, EPS);
}

/*!
 * \brief Computes the GWN for a 3D point wrt a 3D planar polygon
 *
 * \param [in] q The query point to test
 * \param [in] poly The Polygon object
 * \param [out] isOnFace Return variable to show if the point is on the polygon
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * \pre Assumes the polygon is planar. Otherwise, a meaningless value is returned.
 * 
 * Triangulates the polygon and computes the triangular GWN for each component
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 3>& q,
                      const Polygon<T, 3>& poly,
                      bool& isOnFace,
                      const double edge_tol = 1e-8,
                      const double EPS = 1e-8)
{
  const int num_verts = poly.numVertices();
  if(num_verts < 3)
  {
    return 0;
  }

  double wn = 0.0;
  for(int i = 0; i < num_verts - 2; ++i)
  {
    wn +=
      winding_number(q, Triangle<T, 3>(poly[0], poly[i + 1], poly[i + 2]), isOnFace, edge_tol, EPS);
  }

  return wn;
}

/*!
 * \brief Computes the GWN for a 3D point wrt a 3D planar polygon
 *
 * \param [in] q The query point to test
 * \param [in] poly The Polygon object
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * \pre Assumes the polygon is planar. Otherwise, a meaningless value is returned.
 * 
 * Computes the GWN for the polygon without an additional return parameter
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 3>& q,
                      const Polygon<T, 3>& poly,
                      const double edge_tol = 1e-8,
                      const double EPS = 1e-8)
{
  bool isOnFace = false;
  return winding_number(q, poly, isOnFace, edge_tol, EPS);
}

/*!
 * \brief Computes the winding number for a 3D point wrt a 3D convex polyhedron
 *
 * \param [in] q The query point to test
 * \param [in] poly The Polyhedron object
 * \param [in] includeBoundary If true, points on the boundary are considered interior.
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * \pre Expects the polyhedron to be convex and closed so that the returned value is an integer.
 * 
 * Computes the faces of the polyhedron and computes the GWN for each.
 * The sum is then rounded to the nearest integer, as the shape is assumed to be closed.
 * 
 * \return The integer winding number.
 */
template <typename T>
int winding_number(const Point<T, 3>& q,
                   const Polyhedron<T, 3>& poly,
                   bool includeBoundary = false,
                   double edge_tol = 1e-8,
                   double EPS = 1e-8)
{
  SLIC_ASSERT(poly.hasNeighbors());
  const int num_verts = poly.numVertices();

  axom::Array<int> faces(num_verts * num_verts), face_size(2 * num_verts), face_offset(2 * num_verts);
  axom::IndexType face_count;

  poly.getFaces(faces.data(), face_size.data(), face_offset.data(), face_count);

  bool isOnFace = false;
  double wn = 0;
  for(axom::IndexType i = 0; i < face_count; ++i)
  {
    const int N = face_size[i];
    const int i_offset = face_offset[i];
    Polygon<T, 3> the_face(N);
    for(int j = 0; j < N; ++j)
    {
      the_face.addVertex(poly[faces[i_offset + j]]);
    }

    wn += winding_number(q, the_face, isOnFace, edge_tol, EPS);

    if(isOnFace)
    {
      return includeBoundary;
    }
  }

  return std::lround(wn);
}

/*!
 * \brief Computes the GWN for a 3D point wrt a 3D NURBS patch with precomputed data
 *
 * \param [in] query The query point to test
 * \param [in] nurbs The NURBS patch object with data
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] ls_tol The tolerance for the line-surface intersection routine
 * \param [in] quad_tol The maximum relative error allowed in the quadrature
 * \param [in] disk_size The size of extracted disks as a percent of parameter bbox diagonal
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * Computes the generalized winding number for a NURBS patch using Stokes theorem.
 *
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 3>& query,
                      const detail::NURBSPatchGWNCache<T>& nurbs,
                      const double edge_tol = 1e-8,
                      const double ls_tol = 1e-8,
                      const double quad_tol = 1e-8,
                      const double disk_size = 0.01,
                      const double EPS = 1e-8)
{
  // Select the cast direction as an average normal of the untrimmed surface
  auto cast_direction = nurbs.getAverageNormal();
  if(cast_direction.norm() < 1e-10)
  {
    // ...unless the average direction is zero
    double theta = axom::utilities::random_real(0.0, 2 * M_PI);
    double u = axom::utilities::random_real(-1.0, 1.0);
    cast_direction = Vector<T, 3> {sin(theta) * sqrt(1 - u * u), cos(theta) * sqrt(1 - u * u), u};
  }
  else
  {
    cast_direction = cast_direction.unitVector();
  }

  return detail::nurbs_winding_number(query,
                                      nurbs,
                                      cast_direction,
                                      edge_tol,
                                      ls_tol,
                                      quad_tol,
                                      disk_size,
                                      EPS);
}

/*!
 * \brief Computes the GWN for a 3D point wrt a generic 3D surface object
 *
 * \tparam SurfaceType The BezierPatch or NURBSPatch which represents the surface
 * \param [in] query The query point to test
 * \param [in] surf The BezierPatch or NURBSPatch object
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] ls_tol The tolerance for the line-surface intersection routine
 * \param [in] quad_tol The maximum relative error allowed in the quadrature
 * \param [in] disk_size The size of extracted disks as a percent of parameter bbox diagonal
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * \warning Because the constructed cache object is discarded after gwn evaluation,
 *  this method is not accelerated via memoization
 * 
 * \return The GWN.
 */
template <typename T, typename SurfaceType>
double winding_number(const Point<T, 3>& query,
                      const SurfaceType& surf,
                      const double edge_tol = 1e-8,
                      const double ls_tol = 1e-8,
                      const double quad_tol = 1e-8,
                      const double disk_size = 0.01,
                      const double EPS = 1e-8)
{
  return winding_number(query,
                        detail::NURBSPatchGWNCache<T>(surf),
                        edge_tol,
                        ls_tol,
                        quad_tol,
                        disk_size,
                        EPS);
}

/*!
 * \brief Computes the GWN for an array of 3D point wrt an array of NURBS patch cache data
 *
 * \param [in] query_arr The query point to test
 * \param [in] nurbs_arr Array of NURBSPatchGWNCache object containing intermediate values
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] ls_tol The tolerance for the line-surface intersection routine
 * \param [in] quad_tol The maximum relative error allowed in the quadrature
 * \param [in] disk_size The size of extracted disks as a percent of parameter bbox diagonal
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * Computes the generalized winding number for a NURBS patch using Stokes theorem.
 * 
 * \note This method is accelerated via memoization, i.e. dynamically caching and reusing intermediate
 *   values for each curve across query points
 * 
 * \return The array of GWN values.
 */
template <typename T>
axom::Array<double> winding_number(const axom::Array<Point<T, 3>>& query_arr,
                                   const axom::Array<detail::NURBSPatchGWNCache<T>>& nurbs_arr,
                                   const double edge_tol = 1e-8,
                                   const double ls_tol = 1e-8,
                                   const double quad_tol = 1e-8,
                                   const double disk_size = 0.01,
                                   const double EPS = 1e-8)
{
  // Pull precomputed cast directions for each patch
  axom::Array<Vector<T, 3>> cast_direction_arr(0, nurbs_arr.size());
  for(int i = 0; i < nurbs_arr.size(); ++i)
  {
    // Select the cast direction as an average normal of the untrimmed surface
    cast_direction_arr.emplace_back(nurbs_arr[i].getAverageNormal());
    if(cast_direction_arr[i].norm() < 1e-10)
    {
      // ...unless the average direction is zero
      double theta = axom::utilities::random_real(0.0, 2 * M_PI);
      double u = axom::utilities::random_real(-1.0, 1.0);
      cast_direction_arr[i] =
        Vector<T, 3> {sin(theta) * sqrt(1 - u * u), cos(theta) * sqrt(1 - u * u), u};
    }
    else
    {
      cast_direction_arr[i] = cast_direction_arr[i].unitVector();
    }
  }

  axom::Array<double> ret_val(query_arr.size());
  for(int n = 0; n < query_arr.size(); ++n)
  {
    ret_val[n] = 0.0;

    for(int i = 0; i < nurbs_arr.size(); ++i)
    {
      ret_val[n] += detail::nurbs_winding_number(query_arr[n],
                                                 nurbs_arr[i],
                                                 cast_direction_arr[i],
                                                 edge_tol,
                                                 ls_tol,
                                                 quad_tol,
                                                 disk_size,
                                                 EPS);
    }
  }

  return ret_val;
}

/*!
 * \brief Computes the GWN for an array an array of 3D points wrt an array of generic surfaces
 *
 * \tparam SurfaceType The BezierPatch or NURBSPatch which represents the surface
 * \param [in] query_arr The query point to test
 * \param [in] surf_arr Array of NURBSPatch or BezierPatch objects
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] ls_tol The tolerance for the line-surface intersection routine
 * \param [in] quad_tol The maximum relative error allowed in the quadrature
 * \param [in] disk_size The size of extracted disks as a percent of parameter bbox diagonal
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * Computes the generalized winding number for a NURBS patch using Stokes theorem.
 * 
 * \note This method is accelerated via memoization, i.e. dynamically caching and reusing intermediate
 *   values for each curve across query points
 * 
 * \return The array of GWN values.
 */
template <typename T, typename SurfaceType>
axom::Array<double> winding_number(const axom::Array<Point<T, 3>>& query_arr,
                                   const axom::Array<SurfaceType>& surf_arr,
                                   const double edge_tol = 1e-8,
                                   const double ls_tol = 1e-8,
                                   const double quad_tol = 1e-8,
                                   const double disk_size = 0.01,
                                   const double EPS = 1e-8)
{
  // Precompute the expansions and cast directions for each patch
  axom::Array<detail::NURBSPatchGWNCache<T>> nurbs_cache_arr(0, surf_arr.size());
  for(int i = 0; i < surf_arr.size(); ++i)
  {
    nurbs_cache_arr.emplace_back(detail::NURBSPatchGWNCache<T>(surf_arr[i]));
  }

  return winding_number(query_arr, nurbs_cache_arr, edge_tol, ls_tol, quad_tol, disk_size, EPS);
}
///@}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_WINDING_NUMBER_H_
