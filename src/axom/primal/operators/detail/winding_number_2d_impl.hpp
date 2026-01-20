// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef PRIMAL_WINDING_NUMBER_2D_IMPL_HPP_
#define PRIMAL_WINDING_NUMBER_2D_IMPL_HPP_

// Axom includes
#include "axom/config.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/NURBSCurve.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/operators/is_convex.hpp"
#include "axom/primal/operators/squared_distance.hpp"

#include "axom/primal/operators/detail/winding_number_2d_memoization.hpp"

// C++ includes
#include <math.h>

namespace axom
{
namespace primal
{
namespace detail
{
/*!
 * \brief Computes the winding number for a 2D point wrt a 2D polygon
 *
 * \param [in] R The query point to test
 * \param [in] P The Polygon object to test for containment
 * \param [out] isOnEdge An optional return parameter if the point is on the boundary
 * \param [in] includeBoundary If true, points on the boundary are considered interior
 * \param [in] edge_tol The distance at which a point is considered on the boundary
 * 
 * Uses an adapted ray-casting approach that counts quarter-rotation
 * of vertices around the query point. Current policy is to return 1 on edges
 * without strict inclusion, 0 on edges with strict inclusion.
 *
 * The polygon is assumed to be closed, so the winding number is an integer
 * 
 * Directly uses algorithm in 
 * Kai Hormann, Alexander Agathos, "The point in polygon problem for arbitrary polygons"
 * Computational Geometry, Volume 20, Issue 3, 2001,
 * 
 * \return The integer winding number
 */
template <typename T>
int polygon_winding_number(const Point<T, 2>& R,
                           const Polygon<T, 2>& P,
                           bool& isOnEdge,
                           bool includeBoundary,
                           double edge_tol)
{
  const int nverts = P.numVertices();
  const double edge_tol_2 = edge_tol * edge_tol;
  isOnEdge = false;

  int winding_num = 0;
  for(int i = 0; i < nverts; i++)
  {
    int j = (i == nverts - 1) ? 0 : i + 1;

    // Check if the point is on the edge up to some tolerance
    if(squared_distance(R, Segment<T, 2>(P[i], P[j])) <= edge_tol_2)
    {
      isOnEdge = true;
      return includeBoundary ? 1 : 0;
    }

    // Check if edge crosses horizontal line
    if((P[i][1] < R[1]) != (P[j][1] < R[1]))
    {
      if(P[i][0] >= R[0])
      {
        if(P[j][0] > R[0])
        {
          winding_num += 2 * (P[j][1] > P[i][1]) - 1;
        }
        else
        {
          // clang-format off
          double det = axom::numerics::determinant(P[i][0] - R[0], P[j][0] - R[0],
                                                   P[i][1] - R[1], P[j][1] - R[1]);
          // clang-format on

          // Check if edge intersects horitonal ray to the right of R
          if((det > 0) == (P[j][1] > P[i][1]))
          {
            winding_num += 2 * (P[j][1] > P[i][1]) - 1;
          }
        }
      }
      else
      {
        if(P[j][0] > R[0])
        {
          // clang-format off
          double det = axom::numerics::determinant(P[i][0] - R[0], P[j][0] - R[0],
                                                   P[i][1] - R[1], P[j][1] - R[1]);
          // clang-format on

          // Check if edge intersects horitonal ray to the right of R
          if((det > 0) == (P[j][1] > P[i][1]))
          {
            winding_num += 2 * (P[j][1] > P[i][1]) - 1;
          }
        }
      }
    }
  }

  return winding_num;
}

/*!
 * \brief Compute the GWN at a 2D point wrt a 2D line segment
 *
 * \param [in] q The query point to test
 * \param [in] c0 The initial point of the line segment
 * \param [in] c1 The terminal point of the line segment
 * \param [out] isOnEdge Return flag if the point is on the boundary
 * \param [in] edge_tol The tolerance at which a point is on the line
 *
 * The GWN for a 2D point with respect to a 2D straight line
 * is the signed angle subtended by the query point to each endpoint.
 * Colinear points return 0 for their GWN.
 *
 * \return The GWN
 */
template <typename T>
double linear_winding_number(const Point<T, 2>& q,
                             const Point<T, 2>& c0,
                             const Point<T, 2>& c1,
                             bool& isOnEdge,
                             double edge_tol)
{
  const auto V0 = c0 - q;
  const auto V1 = c1 - q;
  const auto seg_vec = c0 - c1;

  // Measures (twice) the signed area of the triangle with vertices q, c0, c1
  const double det_01 = axom::numerics::determinant(V0[0], V1[0], V0[1], V1[1]);

  // Compute distance from line connecting endpoints to query
  isOnEdge = false;
  if(const double tol_sq = edge_tol * edge_tol, seg_sq = seg_vec.squared_norm();
     det_01 * det_01 <= tol_sq * seg_sq)
  {
    // Check for degenerate line segment; treat as a point
    if(axom::utilities::isNearlyEqual(seg_sq, 0.0))
    {
      isOnEdge = (V0.squared_norm() <= tol_sq);
    }
    else if(const double proj_val = V0.dot(seg_vec) / seg_sq;
            (proj_val >= 0 - edge_tol) && (proj_val <= 1 + edge_tol))
    {
      isOnEdge = true;
    }

    return 0;
  }

  constexpr double gwn_modulo = 0.5 * M_1_PI;
  return gwn_modulo * atan2(det_01, V0.dot(V1));
}

/*!
 * \brief Compute the GWN at either endpoint of a 
 *        2D Bezier curve with a convex control polygon
 *
 * \param [in] q The query point
 * \param [in] c The BezierCurve object to compute the winding number along
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance for isNearlyZero
 * \pre Control polygon for c must be convex
 * \pre The query point must be on one of the endpoints
 *
 * The GWN for a Bezier curve with a convex control polygon is
 * given by the signed angle between the tangent vector at that endpoint and
 * the vector in the direction of the other endpoint. 
 * 
 * See Algorithm 2 in
 *  Jacob Spainhour, David Gunderman, and Kenneth Weiss. 2024. 
 *  Robust Containment Queries over Collections of Rational Parametric Curves via Generalized Winding Numbers. 
 *  ACM Trans. Graph. 43, 4, Article 38 (July 2024)
 * 
 * The query can be located on both endpoints if it is closed, in which case
 * the angle is that between the tangent lines at both endpoints
 * 
 * \return The GWN
 */
template <typename T>
double convex_endpoint_winding_number(const Point<T, 2>& q,
                                      const BezierCurve<T, 2>& c,
                                      double edge_tol,
                                      double EPS)
{
  const int ord = c.getOrder();
  if(ord == 1)
  {
    return 0;
  }

  const double edge_tol_sq = edge_tol * edge_tol;

  // Verify that the shape is convex, and that the query point is at an endpoint
  SLIC_ASSERT(is_convex(Polygon<T, 2>(c.getControlPoints()), EPS));
  SLIC_ASSERT((squared_distance(q, c[0]) <= edge_tol_sq) ||
              (squared_distance(q, c[ord]) <= edge_tol_sq));

  // Need to find vectors that subtend the entire curve.
  //   We must ignore duplicate nodes
  int idx;
  for(idx = 0; idx <= ord; ++idx)
  {
    if(squared_distance(q, c[idx]) > edge_tol_sq)
    {
      break;
    }
  }
  const auto V1 = c[idx] - q;

  for(idx = ord; idx >= 0; --idx)
  {
    if(squared_distance(q, c[idx]) > edge_tol_sq)
    {
      break;
    }
  }
  auto V2 = c[idx] - q;

  // Measures the signed area of the triangle spanned by V1 and V2
  double tri_area = axom::numerics::determinant(V1[0], V2[0], V1[1], V2[1]);

  // This means the bounding vectors are anti-parallel.
  //  Parallel tangents can't happen with nontrivial convex control polygons
  if((ord > 3) && axom::utilities::isNearlyEqual(tri_area, 0.0, EPS))
  {
    for(int i = 1; i < ord; ++i)
    {
      // Need to find the first non-parallel control node
      V2 = c[i] - q;
      tri_area = axom::numerics::determinant(V1[0], V2[0], V1[1], V2[1]);

      // Because we are convex, a single non-collinear vertex tells us the orientation
      if(!axom::utilities::isNearlyEqual(tri_area, 0.0, EPS))
      {
        return (tri_area > 0) ? 0.5 : -0.5;
      }
    }

    // If all vectors are parallel, the curve is linear and return 0
    return 0;
  }

  // Compute signed angle between vectors
  constexpr double gwn_modulo = 0.5 * M_1_PI;
  return gwn_modulo * atan2(tri_area, V1.dot(V2));
}

/*!
 * \brief Computes the GWN for a 2D point wrt cached data for a 2D NURBS curve
 *
 * \param [in] q The query point to test
 * \param [in] nurbs_cache The NURBS curve cached data object
 * \param [in] bezier_idx The index of the bezier curve in the extracted NURBS curve
 * \param [in] refinement_level The current subdivision level for the Bezier curve
 * \param [in] refinement_index The index of the specific subdivision at the specified level
 * \param [out] isOnCurve An returned flag if the point is on the curve
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the GWN by decomposing into rational Bezier curves
 *  and summing the resulting GWNs. Far-away curves can be evaluated
 *  without decomposition using direct formula.
 * 
 * \return The GWN.
 */
template <typename T>
double bezier_winding_number_memoized(const Point<T, 2>& q,
                                      const NURBSCurveGWNCache<T>& nurbs_cache,
                                      int bezier_idx,
                                      int refinement_level,
                                      int refinement_index,
                                      bool& isOnCurve,
                                      double edge_tol = 1e-8,
                                      double EPS = 1e-8)
{
  // Early exit for degenerate curves
  const int deg = nurbs_cache.getDegree();
  if(deg <= 0)
  {
    return 0.0;
  }

  auto& bezier_data =
    nurbs_cache.getSubdivisionData(bezier_idx, refinement_level, refinement_index, edge_tol);
  auto& bezier_curve = bezier_data.getCurve();

  // If outside a bounding box, the curve can be treated as linear between its endpoints
  if(!bezier_data.getBoundingBox().contains(q) || bezier_curve.isLinear(EPS))
  {
    return detail::linear_winding_number(q, bezier_curve[0], bezier_curve[deg], isOnCurve, edge_tol);
  }

  // If the control polygon is convex, we can handle coincidence with a direct formula.
  //  In this case, accessing subdivision bounding boxes is fast enough to skip checking a polygon
  if(bezier_data.isConvexControlPolygon() &&
     (squared_distance(q, bezier_curve[0]) <= edge_tol * edge_tol ||
      squared_distance(q, bezier_curve[deg]) <= edge_tol * edge_tol))
  {
    isOnCurve = true;
    return convex_endpoint_winding_number(q, bezier_curve, edge_tol, EPS);
  }

  const auto gwn_half_1 = detail::bezier_winding_number_memoized(q,
                                                                 nurbs_cache,
                                                                 bezier_idx,
                                                                 refinement_level + 1,
                                                                 2 * refinement_index,
                                                                 isOnCurve,
                                                                 edge_tol,
                                                                 EPS);
  const auto gwn_half_2 = detail::bezier_winding_number_memoized(q,
                                                                 nurbs_cache,
                                                                 bezier_idx,
                                                                 refinement_level + 1,
                                                                 2 * refinement_index + 1,
                                                                 isOnCurve,
                                                                 edge_tol,
                                                                 EPS);
  return gwn_half_1 + gwn_half_2;
}

/*!
 * \brief Recursively construct a polygon with the same *integer* winding number 
 *        as the closed original Bezier curve at a given query point.
 *
 * \param [in] q The query point at which to compute winding number
 * \param [in] c A BezierCurve subcurve of the curve along which to compute the winding number
 * \param [in] isConvexControlPolygon Boolean flag if the input Bezier subcurve 
                                      is already convex
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance for nonphysical distances, used in
 *   isLinear, isNearlyZero, is_convex
 * \param [out] approximating_polygon The Polygon that, by termination of recursion,
 *   has the same integer winding number as the original closed curve
 * \param [out] endpoint_gwn A running sum for the exact GWN if the point is at the 
 *   endpoint of a subcurve
 * \param [out] isOnPolygonEndpoint A flag to indicate if the query is on a vertex of the polygon
 *
 * By the termination of the recursive algorithm, `approximating_polygon` contains
 *  a polygon that has the same *integer* winding number as the original curve.
 * 
 * Upon entering this algorithm, the closing line of `c` is already an
 *  edge of the approximating polygon.
 * If q is outside a convex shape that contains the entire curve, the 
 *  integer winding number for the *closed* curve `c` is zero, 
 *  and the algorithm terminates.
 * If the shape is not convex or we're inside it, instead add the midpoint 
 *  as a vertex and repeat the algorithm. 
 */
template <typename T>
void construct_approximating_polygon(const Point<T, 2>& q,
                                     const BezierCurve<T, 2>& c,
                                     bool isConvexControlPolygon,
                                     double edge_tol,
                                     double EPS,
                                     Polygon<T, 2>& approximating_polygon,
                                     double& endpoint_gwn,
                                     bool& isOnPolygonEndpoint)
{
  const int ord = c.getOrder();

  // Simplest convex shape containing c is its bounding box
  if(!c.boundingBox().expand(edge_tol).contains(q))
  {
    return;
  }

  // Use linearity as base case for recursion
  if(c.isLinear(EPS))
  {
    return;
  }

  // Check if our control polygon is convex.
  //  If so, all subsequent control polygons will be convex as well
  Polygon<T, 2> controlPolygon(c.getControlPoints());
  const bool includeBoundary = true;
  bool isOnEdge = false;

  if(!isConvexControlPolygon)
  {
    isConvexControlPolygon = is_convex(controlPolygon, EPS);
  }

  // Formulas for winding number only work if shape is convex
  if(isConvexControlPolygon)
  {
    // Bezier curves are always contained in their convex control polygon
    if(polygon_winding_number(q, controlPolygon, isOnEdge, includeBoundary, edge_tol) == 0)
    {
      return;
    }

    // If the query point is at either endpoint...
    if(const double edge_tol_sq = edge_tol * edge_tol;
       squared_distance(q, c[0]) <= edge_tol_sq || squared_distance(q, c[ord]) <= edge_tol_sq)
    {
      // ...we can use a direct formula for the GWN at the endpoint
      endpoint_gwn += convex_endpoint_winding_number(q, c, edge_tol, EPS);
      isOnPolygonEndpoint = true;

      return;
    }
  }

  // Recursively split curve until query is outside some known convex region
  BezierCurve<T, 2> c1, c2;
  c.split(0.5, c1, c2);

  construct_approximating_polygon(q,
                                  c1,
                                  isConvexControlPolygon,
                                  edge_tol,
                                  EPS,
                                  approximating_polygon,
                                  endpoint_gwn,
                                  isOnPolygonEndpoint);
  approximating_polygon.addVertex(c2[0]);
  construct_approximating_polygon(q,
                                  c2,
                                  isConvexControlPolygon,
                                  edge_tol,
                                  EPS,
                                  approximating_polygon,
                                  endpoint_gwn,
                                  isOnPolygonEndpoint);

  return;
}

/*!
 * \brief Computes the GWN for a 2D point wrt a 2D Bezier curve
 *
 * \param [in] q The query point to test
 * \param [in] c The Bezier curve object 
 * \param [out] isOnCurve An returned flag if the point is on the curve
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the GWN using a recursive, bisection algorithm
 * that constructs a polygon with the same *integer* WN as 
 * the curve closed with a linear segment. The *generalized* WN 
 * of the closing line is then subtracted from the integer WN to
 * return the GWN of the original curve.
 *  
 * Nearly-linear Bezier curves are the base case for recursion.
 * 
 * See Algorithm 2 in
 *  Jacob Spainhour, David Gunderman, and Kenneth Weiss. 2024. 
 *  Robust Containment Queries over Collections of Rational Parametric Curves via Generalized Winding Numbers. 
 *  ACM Trans. Graph. 43, 4, Article 38 (July 2024)
 * 
 * \return The GWN.
 */
template <typename T>
double bezier_winding_number(const Point<T, 2>& q,
                             const BezierCurve<T, 2>& c,
                             bool& isOnCurve,
                             double edge_tol = 1e-8,
                             double EPS = 1e-8)
{
  const int ord = c.getOrder();
  if(ord <= 0)
  {
    return 0.0;
  }

  // Early return is possible for most points + curves
  if(!c.boundingBox().expand(edge_tol).contains(q) || c.isLinear(EPS))
  {
    return detail::linear_winding_number(q, c[0], c[ord], isOnCurve, edge_tol);
  }

  // The first vertex of the polygon is the t=0 point of the curve
  Polygon<T, 2> approximating_polygon(1);
  approximating_polygon.addVertex(c[0]);

  // Need to keep a running total of the GWN to account for
  //  the winding number of coincident points
  double gwn = 0.0;
  bool isOnPolygonEndpoint = false;  // Indicates coincidence with curve
  detail::construct_approximating_polygon(q,
                                          c,
                                          false,
                                          edge_tol,
                                          EPS,
                                          approximating_polygon,
                                          gwn,
                                          isOnPolygonEndpoint);

  // The last vertex of the polygon is the t=1 point of the curve
  approximating_polygon.addVertex(c[ord]);

  // Compute the integer winding number of the closed curve
  bool isOnPolygonEdge = false;
  double closed_curve_wn =
    detail::polygon_winding_number(q, approximating_polygon, isOnPolygonEdge, false, edge_tol);

  // Compute the fractional value of the closed curve
  bool isOnClosure = false;
  const int n = approximating_polygon.numVertices();
  const double closure_wn = detail::linear_winding_number(q,
                                                          approximating_polygon[n - 1],
                                                          approximating_polygon[0],
                                                          isOnClosure,
                                                          edge_tol);

  // If the point is on the boundary of the approximating polygon,
  //  or coincident with the curve (isOnPolygonEndpoint), then winding_number<polygon>
  //  doesn't return the right half-integer. Have to go edge-by-edge.
  if(isOnPolygonEndpoint || isOnPolygonEdge)
  {
    closed_curve_wn = closure_wn;
    for(int i = 1; i < n; ++i)
    {
      bool isOnThisEdge = false;
      closed_curve_wn += detail::linear_winding_number(q,
                                                       approximating_polygon[i - 1],
                                                       approximating_polygon[i],
                                                       isOnThisEdge,
                                                       edge_tol);
    }
  }

  isOnCurve = isOnPolygonEndpoint;
  return gwn + closed_curve_wn - closure_wn;
}

/*!
 * \brief Computes the GWN for a 2D point wrt a 2D NURBS curve
 *
 * \param [in] q The query point to test
 * \param [in] n The NURBS curve object 
 * \param [out] isOnCurve An returned flag if the point is on the curve
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the GWN by decomposing into rational Bezier curves
 *  and summing the resulting GWNs. Far-away curves can be evaluated
 *  without decomposition using direct formula.
 * 
 * \return The GWN.
 */
template <typename T>
double nurbs_winding_number(const Point<T, 2>& q,
                            const NURBSCurve<T, 2>& nurbs,
                            bool& isOnCurve,
                            double edge_tol = 1e-8,
                            double EPS = 1e-8)
{
  const int deg = nurbs.getDegree();
  if(deg <= 0)
  {
    return 0.0;
  }

  // Early return is possible for most points + curves
  if(!nurbs.boundingBox().expand(edge_tol).contains(q))
  {
    return detail::linear_winding_number(q,
                                         nurbs[0],
                                         nurbs[nurbs.getNumControlPoints() - 1],
                                         isOnCurve,
                                         edge_tol);
  }

  // Decompose the NURBS curve into Bezier segments
  auto beziers = nurbs.extractBezier();

  // Compute the GWN for each Bezier segment
  double gwn = 0.0;
  for(int i = 0; i < beziers.size(); i++)
  {
    bool isOnThisCurve = false;
    gwn += detail::bezier_winding_number(q, beziers[i], isOnThisCurve, edge_tol, EPS);
    isOnCurve = isOnCurve || isOnThisCurve;
  }

  return gwn;
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif
