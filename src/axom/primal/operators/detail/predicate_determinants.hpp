// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file predicate_determinants.hpp
 *
 * \brief Low-level determinant helpers for common computational geometry predicates.
 *
 * These routines centralize the core determinant computations used by primal
 * predicates (e.g. orientation and in-sphere).  They return the raw determinant
 * values without interpreting tolerances or mapping to OrientationResult.
 */

#ifndef AXOM_PRIMAL_PREDICATE_DETERMINANTS_HPP_
#define AXOM_PRIMAL_PREDICATE_DETERMINANTS_HPP_

#include "axom/core/numerics/Determinants.hpp"

#include "axom/primal/geometry/Point.hpp"

namespace axom
{
namespace primal
{
namespace detail
{

/*!
 * \brief Returns the raw 2D orientation determinant for three points.
 *
 * This determinant is twice the signed area of the triangle (a,b,c).
 */
template <typename T>
inline double orientation_determinant(const Point<T, 2>& a, const Point<T, 2>& b, const Point<T, 2>& c)
{
  const auto ba = b - a;
  const auto ca = c - a;

  return axom::numerics::determinant(ba[0], ba[1], ca[0], ca[1]);
}

/*!
 * \brief Returns the raw 3D orientation determinant for four points.
 *
 * This determinant is six times the signed volume of the tetrahedron (a,b,c,d).
 */
template <typename T>
inline double orientation_determinant(const Point<T, 3>& a,
                                      const Point<T, 3>& b,
                                      const Point<T, 3>& c,
                                      const Point<T, 3>& d)
{
  const auto ba = b - a;
  const auto ca = c - a;
  const auto da = d - a;

  // clang-format off
  return axom::numerics::determinant(
    ba[0], ba[1], ba[2],
    ca[0], ca[1], ca[2],
    da[0], da[1], da[2]);
  // clang-format on
}

/*!
 * \brief Returns the raw in-sphere determinant for a 2D triangle circumcircle test.
 *
 * The sign convention matches primal::in_sphere(): a negative determinant means
 * the query point is inside the circumcircle for a consistently oriented input.
 */
template <typename T>
inline double in_sphere_determinant(const Point<T, 2>& q,
                                    const Point<T, 2>& p0,
                                    const Point<T, 2>& p1,
                                    const Point<T, 2>& p2)
{
  const auto ba = p1 - p0;
  const auto ca = p2 - p0;
  const auto qa = q - p0;

  // clang-format off
  return axom::numerics::determinant(
    ba[0], ba[1], ba.squared_norm(),
    ca[0], ca[1], ca.squared_norm(),
    qa[0], qa[1], qa.squared_norm());
  // clang-format on
}

/*!
 * \brief Returns the raw in-sphere determinant for a 3D tetrahedron circumsphere test.
 *
 * The sign convention matches primal::in_sphere(): a negative determinant means
 * the query point is inside the circumsphere for a consistently oriented input.
 */
template <typename T>
inline double in_sphere_determinant(const Point<T, 3>& q,
                                    const Point<T, 3>& p0,
                                    const Point<T, 3>& p1,
                                    const Point<T, 3>& p2,
                                    const Point<T, 3>& p3)
{
  const auto ba = p1 - p0;
  const auto ca = p2 - p0;
  const auto da = p3 - p0;
  const auto qa = q - p0;

  // clang-format off
  return axom::numerics::determinant(
    ba[0], ba[1], ba[2], ba.squared_norm(),
    ca[0], ca[1], ca[2], ca.squared_norm(),
    da[0], da[1], da[2], da.squared_norm(),
    qa[0], qa[1], qa[2], qa.squared_norm());
  // clang-format on
}

}  // namespace detail
}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_PREDICATE_DETERMINANTS_HPP_
