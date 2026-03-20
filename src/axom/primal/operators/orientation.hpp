// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file orientation.hpp
 *
 * \brief Consists of a set of templated (overloaded) routines used to calculate
 *  the orientation of a given point to another geometric entity.
 *
 */

#ifndef AXOM_PRIMAL_ORIENTATION_HPP_
#define AXOM_PRIMAL_ORIENTATION_HPP_

#include "axom/core/numerics/Determinants.hpp"
#include "axom/core/utilities/Utilities.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/OrientationResult.hpp"

#include "axom/primal/operators/detail/predicate_determinants.hpp"

#include "axom/slic/interface/slic.hpp"

namespace axom
{
namespace primal
{
namespace robust
{
/*!
 * \brief Computes the orientation of a point \a p with respect to an oriented triangle \a tri.
 *
 * \return ON_BOUNDARY if within tolerance, ON_POSITIVE_SIDE / ON_NEGATIVE_SIDE otherwise.
 */
template <typename T>
inline int orientation(const Point<T, 3>& p, const Triangle<T, 3>& tri, double EPS = 1e-9)
{
  const double det = detail::orientation_determinant(p, tri[0], tri[1], tri[2]);

  if(axom::utilities::isNearlyEqual(det, 0., EPS))
  {
    return primal::ON_BOUNDARY;
  }

  // Preserve existing convention: det < 0 implies ON_POSITIVE_SIDE.
  return det < 0. ? primal::ON_POSITIVE_SIDE : primal::ON_NEGATIVE_SIDE;
}

/*!
 * \brief Computes the orientation of a point \a p with respect to an oriented segment \a seg.
 *
 * \return ON_BOUNDARY if within tolerance, ON_POSITIVE_SIDE / ON_NEGATIVE_SIDE otherwise.
 */
template <typename T>
inline int orientation(const Point<T, 2>& p, const Segment<T, 2>& seg, double EPS = 1e-9)
{
  const double det = detail::orientation_determinant(p, seg[0], seg[1]);

  if(axom::utilities::isNearlyEqual(det, 0., EPS))
  {
    return primal::ON_BOUNDARY;
  }

  // Preserve existing convention: det < 0 implies ON_POSITIVE_SIDE.
  return det < 0. ? primal::ON_POSITIVE_SIDE : primal::ON_NEGATIVE_SIDE;
}

}  // namespace robust

/*!
 * \brief Computes the orientation of a point \a p with respect to an
 *  oriented triangle \a tri
 * \param [in] p the query point
 * \param [in] tri an oriented triangle
 * \param [in] EPS a tolerance for determining if \a p and \a tri are coplanar
 * \return The orientation of \a p with respect to \a tri
 * \note The triangle lies in a plane that divides space into the positive
 * half-space and the negative half-space. The triangle's normal vector
 * points in the positive half-space.
 * The return value of this routine can be one of the following:
 * <ul>
 *  <li> ON_BOUNDARY, when \a p is coplanar with \a tri (within tolerance \a EPS)</li>
 *  <li> ON_POSITIVE_SIDE, when \a p lies in the positive half-space </li>
 *  <li> ON_NEGATIVE_SIDE, when \a p lies in the negative half-space </li>
 * </ul>
 * \sa OrientationResult
 */
template <typename T>
inline int orientation(const Point<T, 3>& p, const Triangle<T, 3>& tri, double EPS = 1e-9)
{
  return robust::orientation(p, tri, EPS);
}

/*!
 * \brief Computes the orientation of a point \a p with respect to an
 *  oriented segment
 * \param [in] p the query point
 * \param [in] seg an oriented segment
 * \param [in] EPS a tolerance for determining if \a p and \a seg are coplanar
 * \return The orientation of \a p with respect to \a seg
 * \note The segment lies in a plane that divides space into the positive
 * half-space and the negative half-space. The segment's normal vector
 * points in the positive half-space.
 * The return value can be one of the following:
 * <ul>
 *  <li> ON_BOUNDARY, when \a p is coplanar with \a seg (within tolerance \a EPS)</li>
 *  <li> ON_POSITIVE_SIDE, when \a p lies in the positive half-space </li>
 *  <li> ON_NEGATIVE_SIDE, when \a p lies in the negative half-space </li>
 * </ul>
 * \sa OrientationResult
 */
template <typename T>
inline int orientation(const Point<T, 2>& p, const Segment<T, 2>& seg, double EPS = 1e-9)
{
  return robust::orientation(p, seg, EPS);
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_ORIENTATION_HPP_
