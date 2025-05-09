// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_NUMERICS_FLOATING_POINT_LIMITS_HPP_
#define AXOM_NUMERICS_FLOATING_POINT_LIMITS_HPP_

#include "axom/core/Macros.hpp"  // for Axom macros

#include <type_traits>  // for std::is_floating_point
#include <cfloat>       // for floating point limits

namespace axom
{
namespace numerics
{
/*!
 * \accelerated
 * \brief Traits class for floating point types to query limit information.
 *  The floating_point_limits class is traits class providing a standardized
 *  and portable way to query the following information on either host/device:
 *  <ul>
 *   <li> lowest()  -- returns the lowest finite value </li>
 *   <li> min()     -- returns the smallest finite value </li>
 *   <li> max()     -- returns the largest finite value </li>
 *   <li> epsilon() -- returns difference between 1.0 and the next value </li>
 *  </ul>
 *
 *  \tparam T the floating point type, e.g., float, double, or long double.
 */
template <typename T>
struct floating_point_limits
{
  AXOM_STATIC_ASSERT_MSG(std::is_floating_point<T>::value,
                         "floating_point_limits< T > must be used with a floating type!");
};

//------------------------------------------------------------------------------
template <>
struct floating_point_limits<float>
{
  AXOM_HOST_DEVICE
  static constexpr float lowest() { return -FLT_MAX; };

  AXOM_HOST_DEVICE
  static constexpr float min() { return FLT_MIN; };

  AXOM_HOST_DEVICE
  static constexpr float max() { return FLT_MAX; };

  AXOM_HOST_DEVICE
  static constexpr float epsilon() { return FLT_EPSILON; };
};

//------------------------------------------------------------------------------
template <>
struct floating_point_limits<double>
{
  AXOM_HOST_DEVICE
  static constexpr double lowest() { return -DBL_MAX; };

  AXOM_HOST_DEVICE
  static constexpr double min() { return DBL_MIN; };

  AXOM_HOST_DEVICE
  static constexpr double max() { return DBL_MAX; };

  AXOM_HOST_DEVICE
  static constexpr double epsilon() { return DBL_EPSILON; };
};

//------------------------------------------------------------------------------
template <>
struct floating_point_limits<long double>
{
  AXOM_HOST_DEVICE
  static constexpr long double lowest() { return -LDBL_MAX; };

  AXOM_HOST_DEVICE
  static constexpr long double min() { return LDBL_MIN; };

  AXOM_HOST_DEVICE
  static constexpr long double max() { return LDBL_MAX; };

  AXOM_HOST_DEVICE
  static constexpr long double epsilon() { return LDBL_EPSILON; };
};

} /* namespace numerics */
} /* namespace axom */

#endif /* AXOM_NUMERICS_FLOATING_POINT_LIMITS_HPP_ */
