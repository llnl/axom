// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_COMPUTE_MOMENTS_HPP_
#define AXOM_PRIMAL_COMPUTE_MOMENTS_HPP_

/*!
 * \file compute_moments.hpp
 *
 * \brief Consists of a set of methods to compute areas/volumes and centroids 
 * for Polygon and CurvedPolygon objects composed of nonrational BezierCurve objects
 */

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"
#include "axom/primal/operators/detail/compute_moments_impl.hpp"

#include <vector>
#include <map>

namespace axom
{
namespace primal
{
/*!
   * \brief Calculates the sector area of a planar, nonrational Bezier Curve
   *
   * The sector area is the area between the curve and the origin.
   * The equation and derivation is described in:
   *  Ueda, K. "Signed area of sectors between spline curves and the origin"
   *  IEEE International Conference on Information Visualization, 1999.
   */
template <typename T>
T sector_area(const primal::BezierCurve<T, 2>& curve)
{
  // Algorithm works only on nonrational Bezier curves
  SLIC_ASSERT(!curve.isRational());

  // Weights for each polynomial order are precomputed and memoized
  static detail::MemoizedSectorAreaWeights<T> s_weights;

  T A = 0;
  const int ord = curve.getOrder();
  const auto& weights = s_weights.getWeights(ord);

  for(int p = 0; p <= ord; ++p)
  {
    for(int q = 0; q <= ord; ++q)
    {
      A += weights(p, q) * curve[p][1] * curve[q][0];
    }
  }
  return A;
}

/*!
   * \brief Calculates the sector centroid of a planar, nonrational Bezier Curve
   *
   * This is the centroid of the region between the curve and the origin.
   * The equation and derivation are generalizations of:
   *  Ueda, K. "Signed area of sectors between spline curves and the origin"
   *  IEEE International Conference on Information Visualization, 1999.
   */
template <typename T>
primal::Point<T, 2> sector_centroid(const primal::BezierCurve<T, 2>& curve)
{
  // Algorithm works only on nonrational Bezier curves
  SLIC_ASSERT(!curve.isRational());

  // Weights for each polynomial order's centroid are precomputed and memoized
  static detail::MemoizedSectorCentroidWeights<T> s_weights;

  T Mx = 0;
  T My = 0;
  const int ord = curve.getOrder();
  for(int r = 0; r <= ord; ++r)
  {
    const auto& weights_r = s_weights.getWeights(ord, r);
    for(int p = 0; p <= ord; ++p)
    {
      for(int q = 0; q <= ord; ++q)
      {
        Mx += weights_r(p, q) * curve[p][1] * curve[q][0] * curve[r][0];
        My += weights_r(p, q) * curve[p][1] * curve[q][0] * curve[r][1];
      }
    }
  }
  return primal::Point<T, 2> {Mx, My};
}

/// \brief Returns the area enclosed by the CurvedPolygon
template <typename T>
T area(const primal::CurvedPolygon<BezierCurve<T, 2>>& poly, double tol = 1e-8)
{
  const int ngon = poly.numEdges();
  T A = 0.0;
  if(!poly.isClosed(1e3 * tol))
  {
    SLIC_DEBUG("Warning! The area is 0 because the curved polygon is not closed.");
    return A;
  }
  else
  {
    for(int ed = 0; ed < ngon; ++ed)
    {
      // Algorithm works only on nonrational Bezier curves
      SLIC_ASSERT(!poly[ed].isRational());

      A += primal::sector_area(poly[ed]);
    }
    return A;
  }
}

/// \brief Returns the centroid of the CurvedPolygon
template <typename T>
primal::Point<T, 2> centroid(const primal::CurvedPolygon<BezierCurve<T, 2>>& poly, double tol = 1e-8)
{
  using PointType = primal::Point<T, 2>;

  const int ngon = poly.numEdges();
  PointType M;

  if(!poly.isClosed(1e3 * tol))
  {
    SLIC_DEBUG("Warning! The moments are 0 because the curved polygon is not closed.");
    return M;
  }
  else
  {
    const T A = area(poly, tol);
    if(A != 0.)
    {
      for(int ed = 0; ed < ngon; ++ed)
      {
        // Algorithm works only on nonrational Bezier curves
        SLIC_ASSERT(!poly[ed].isRational());

        M.array() += primal::sector_centroid(poly[ed]).array();
      }
      M.array() /= A;
    }
    return M;
  }
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_COMPUTE_MOMENTS_HPP_
