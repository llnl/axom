// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file in_curved_polygon.hpp
 *
 * \brief Consists of methods that test whether a given query point is
 * inside a given curved polygon.
 *
 * Uses an adaptive winding number calculation
 */

#ifndef AXOM_PRIMAL_IN_CURVED_POLYGON_HPP_
#define AXOM_PRIMAL_IN_CURVED_POLYGON_HPP_

// Axom includes
#include "axom/config.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"

#include "axom/primal/operators/winding_number.hpp"

namespace axom
{
namespace primal
{
/*!
 * \brief Robustly determine if query point is interior to a curved polygon
 *
 * \tparam CurveType The BezierCurve, NURBSCurve, or NURBSCurveGWNCache which represents the curve
 * \param [in] query The query point to test
 * \param [in] cpoly The CurvedPolygon object to test for containment
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * Determines containment using the (rounded) winding number with respect
 * to the given curved polygon. This algorithm is robust, as the winding number is rounded 
 * in the final in/out determination. 
 * Different protocols determine containment from the winding number differently.
 *   Nonzero Rule:  If the winding number is nonzero, the point is interior.
 *   Even/Odd rule: If the winding number is odd, it is interior. Exterior otherwise.
 *
 * \return A boolean value indicating containment.
 */
template <typename T, typename CurveType>
bool in_curved_polygon(const Point<T, 2>& query,
                       const CurvedPolygon<CurveType>& cpoly,
                       bool useNonzeroRule = true,
                       double edge_tol = 1e-8,
                       double EPS = 1e-8)
{
  return useNonzeroRule ? (std::lround(winding_number(query, cpoly, edge_tol, EPS)) != 0)
                        : (std::lround(winding_number(query, cpoly, edge_tol, EPS)) % 2) != 0;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_IN_CURVED_POLYGON_H_
