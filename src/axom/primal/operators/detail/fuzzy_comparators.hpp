// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file fuzzy_comparators.hpp
 *
 * This file provides helper functions for fuzzy comparisons
 */

#ifndef AXOM_PRIMAL_FUZZY_COMPARATORS_HPP_
#define AXOM_PRIMAL_FUZZY_COMPARATORS_HPP_

#include "axom/core/Macros.hpp"
#include "axom/core/utilities/Utilities.hpp"

namespace axom
{
namespace primal
{
namespace detail
{
/// \brief Checks if x > y, within a specified tolerance.
AXOM_HOST_DEVICE
inline bool isGt(double x, double y, double EPS = 1e-12)
{
  return ((x > y) && !(axom::utilities::isNearlyEqual(x, y, EPS)));
}

/// \brief Checks if x < y, within a specified tolerance.
AXOM_HOST_DEVICE
inline bool isLt(double x, double y, double EPS = 1e-12)
{
  return ((x < y) && !(axom::utilities::isNearlyEqual(x, y, EPS)));
}

/// \brief Checks if x <= y, within a specified tolerance.
AXOM_HOST_DEVICE
inline bool isLeq(double x, double y, double EPS = 1e-12) { return !(isGt(x, y, EPS)); }

/// \brief Checks if x >= y, within a specified tolerance.
AXOM_HOST_DEVICE
inline bool isGeq(double x, double y, double EPS = 1e-12) { return !(isLt(x, y, EPS)); }

/*!
 * \brief Checks if x < y, or possibly x == y, within a specified tolerance.
 *
 * The check for equality is controlled by parameter includeEqual.  This
 * lets users specify at compile time whether triangles intersecting only on
 * border points are reported as intersecting or not.
 *
 * Supports checkEdge and checkVertex
 */
AXOM_HOST_DEVICE
inline bool isLpeq(double x, double y, bool includeEqual = false, double EPS = 1e-12)
{
  if(includeEqual && axom::utilities::isNearlyEqual(x, y, EPS))
  {
    return true;
  }

  return isLt(x, y, EPS);
}

/*!
 * \brief Checks if x > y, or possibly x == y, within a specified tolerance.
 *
 * The check for equality is controlled by parameter includeEqual.  This
 * lets users specify at compile time whether triangles intersecting only on
 * border points are reported as intersecting or not.
 *
 * Supports checkEdge and checkVertex
 */
AXOM_HOST_DEVICE
inline bool isGpeq(double x, double y, bool includeEqual = false, double EPS = 1e-12)
{
  if(includeEqual && axom::utilities::isNearlyEqual(x, y, EPS))
  {
    return true;
  }

  return isGt(x, y, EPS);
}

}  // namespace detail
}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_FUZZY_COMPARATORS_HPP_