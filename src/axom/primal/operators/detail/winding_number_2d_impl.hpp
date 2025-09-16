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

// C++ includes
#include <math.h>

// MFEM includes
#ifdef AXOM_USE_MFEM
  #include "mfem.hpp"
#endif

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
 * \param [in] isOnEdge An optional return parameter if the point is on the boundary
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

/*
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
  Vector<T, 2> V0(q, c0);
  Vector<T, 2> V1(q, c1);

  // clang-format off
  // Measures the signed area of the triangle with vertices q, c0, c1
  double tri_area = axom::numerics::determinant(V0[0] - V1[0], V1[0], 
                                                V0[1] - V1[1], V1[1]);
  // clang-format on

  double segment_sq_len = (V0 - V1).squared_norm();

  // Compute distance from line connecting endpoints to query
  isOnEdge = false;
  if(tri_area * tri_area <= edge_tol * edge_tol * segment_sq_len)
  {
    // Project the query point onto the line segment to see if they're coincident
    if(axom::utilities::isNearlyEqual(segment_sq_len, 0.0))
    {
      isOnEdge = (V0.squared_norm() <= edge_tol * edge_tol);
    }
    else
    {
      double proj_val = Vector<T, 2>::dot_product(V0, V0 - V1) / segment_sq_len;
      if((proj_val >= 0 - edge_tol) && (proj_val <= 1 + edge_tol))
      {
        isOnEdge = true;
      }
    }

    return 0;
  }

  // Compute signed angle between vectors
  double dotprod =
    axom::utilities::clampVal(Vector<T, 2>::dot_product(V0.unitVector(), V1.unitVector()), -1.0, 1.0);

  return 0.5 * M_1_PI * acos(dotprod) * ((tri_area > 0) ? 1 : -1);
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

  double edge_tol_sq = edge_tol * edge_tol;

  // Verify that the shape is convex, and that the query point is at an endpoint
  SLIC_ASSERT(is_convex(Polygon<T, 2>(c.getControlPoints()), EPS));
  SLIC_ASSERT((squared_distance(q, c[0]) <= edge_tol_sq) ||
              (squared_distance(q, c[ord]) <= edge_tol_sq));

  int idx;

  // Need to find vectors that subtend the entire curve.
  //   We must ignore duplicate nodes
  for(idx = 0; idx <= ord; ++idx)
  {
    if(squared_distance(q, c[idx]) > edge_tol_sq)
    {
      break;
    }
  }
  Vector<T, 2> V1(q, c[idx]);

  for(idx = ord; idx >= 0; --idx)
  {
    if(squared_distance(q, c[idx]) > edge_tol_sq)
    {
      break;
    }
  }
  Vector<T, 2> V2(q, c[idx]);

  // clang-format off
  // Measures the signed area of the triangle spanned by V1 and V2
  double tri_area = axom::numerics::determinant(V1[0] - V2[0], V2[0], 
                                                V1[1] - V2[1], V2[1]);
  // clang-format on

  // This means the bounding vectors are anti-parallel.
  //  Parallel tangents can't happen with nontrivial convex control polygons
  if((ord > 3) && axom::utilities::isNearlyEqual(tri_area, 0.0, EPS))
  {
    for(int i = 1; i < ord; ++i)
    {
      // Need to find the first non-parallel control node
      V2 = Vector<T, 2>(q, c[i]);

      // clang-format off
      tri_area = axom::numerics::determinant(V1[0] - V2[0], V2[0], 
                                             V1[1] - V2[1], V2[1]);
      // clang-format on

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
  double dotprod =
    axom::utilities::clampVal(Vector<T, 2>::dot_product(V1.unitVector(), V2.unitVector()), -1.0, 1.0);
  return 0.5 * M_1_PI * acos(dotprod) * ((tri_area > 0) ? 1 : -1);
}

/*!
 * \brief Computes the GWN for a 2D point wrt cached data for a 2D NURBS curve
 *
 * \param [in] q The query point to test
 * \param [in] nurbs_cache The NURBS curve cached data object 
 * \param [in] isOnEdge An returned flag if the point is on the curve
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

  auto& bezier_data = nurbs_cache.getSubdivisionData(bezier_idx, refinement_level, refinement_index);
  auto& bezier_curve = bezier_data.curve;

  // If outside a bounding box, the curve can be treated as linear between its endpoints
  if(!bezier_curve.boundingBox().expand(edge_tol).contains(q) || bezier_curve.isLinear(EPS))
  {
    return detail::linear_winding_number(q, bezier_curve[0], bezier_curve[deg], isOnCurve, edge_tol);
  }

  // If the control polygon is convex, we can handle coincidence with a direct formula
  if(bezier_data.isConvexControlPolygon &&
     (squared_distance(q, bezier_curve[0]) <= edge_tol * edge_tol ||
      squared_distance(q, bezier_curve[deg]) <= edge_tol * edge_tol))
  {
    isOnCurve = true;
    return convex_endpoint_winding_number(q, bezier_curve, edge_tol, EPS);
  }

  return detail::bezier_winding_number_memoized(q,
                                                nurbs_cache,
                                                bezier_idx,
                                                refinement_level + 1,
                                                2 * refinement_index,
                                                isOnCurve,
                                                edge_tol,
                                                EPS) +
    detail::bezier_winding_number_memoized(q,
                                           nurbs_cache,
                                           bezier_idx,
                                           refinement_level + 1,
                                           2 * refinement_index + 1,
                                           isOnCurve,
                                           edge_tol,
                                           EPS);
}

/*!
 * \brief Computes the GWN for a 2D point wrt a 2D Bezier curve
 *
 * \param [in] q The query point to test
 * \param [in] bezier The Bezier curve object 
 * \param [in] isOnEdge An returned flag if the point is on the curve
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Evaluates the GWN through recursive subdivision of the curve until each component is "far-away"
 *  i.e. the query it outside its bounding box or convex control polygon, at which point a direct formula for the GWN 
 *  of that component can be used.
 *  
 * \return The GWN.
 */
template <typename T>
double bezier_winding_number(const Point<T, 2>& q,
                             const BezierCurve<T, 2>& bezier,
                             bool isConvexControlPolygon,
                             bool& isOnCurve,
                             double edge_tol = 1e-8,
                             double EPS = 1e-8)
{
  // Early exit for degenerate curves
  const int ord = bezier.getOrder();
  if(ord <= 0)
  {
    return 0.0;
  }

  // If outside a bounding box, the curve can be treated as linear between its endpoints
  if(!bezier.boundingBox().expand(edge_tol).contains(q) || bezier.isLinear(EPS))
  {
    return detail::linear_winding_number(q, bezier[0], bezier[ord], isOnCurve, edge_tol);
  }

  // If the control polygon isn't convex, we need to subdivide regardless
  auto controlPolygon = Polygon<T, 2>(bezier.getControlPoints());
  if(!isConvexControlPolygon)
  {
    isConvexControlPolygon = is_convex(controlPolygon, EPS);
  }
  // Otherwise, we can check that we're outside the convex control polygon
  //  and handle coincidence with a direct formula
  else
  {
    constexpr bool includeBoundary = true;
    bool isOnEdge = false;
    if(polygon_winding_number(q, controlPolygon, isOnEdge, includeBoundary, edge_tol) == 0)
    {
      return detail::linear_winding_number(q, bezier[0], bezier[ord], isOnCurve, edge_tol);
    }

    if(squared_distance(q, bezier[0]) <= edge_tol * edge_tol ||
       squared_distance(q, bezier[ord]) <= edge_tol * edge_tol)
    {
      isOnCurve = true;
      return convex_endpoint_winding_number(q, bezier, edge_tol, EPS);
    }
  }

  BezierCurve<T, 2> b1, b2;
  bezier.split(0.5, b1, b2);

  return detail::bezier_winding_number(q, b1, isConvexControlPolygon, isOnCurve, edge_tol, EPS) +
    detail::bezier_winding_number(q, b2, isConvexControlPolygon, isOnCurve, edge_tol, EPS);
}

/*!
 * \brief Computes the GWN for a 2D point wrt a 2D NURBS curve
 *
 * \param [in] q The query point to test
 * \param [in] n The NURBS curve object 
 * \param [in] isOnEdge An returned flag if the point is on the curve
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
    gwn += detail::bezier_winding_number(q, beziers[i], false, isOnThisCurve, edge_tol, EPS);
    isOnCurve = isOnCurve || isOnThisCurve;
  }

  return gwn;
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif
