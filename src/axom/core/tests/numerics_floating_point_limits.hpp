// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// axom includes
#include "axom/core/Macros.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/core/numerics/floating_point_limits.hpp"

// gtest includes
#include "gtest/gtest.h"

// C/C++ includes
#include <type_traits>  // for std::is_floating_point

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
template <typename T>
void check_type_limits(const std::string& typeName)
{
  AXOM_STATIC_ASSERT(std::is_floating_point<T>::value);

  SCOPED_TRACE("Testing [" + typeName + "]");

  const T EPS = axom::numerics::floating_point_limits<T>::epsilon();
  const T LOWEST = axom::numerics::floating_point_limits<T>::lowest();
  const T MIN = axom::numerics::floating_point_limits<T>::min();
  const T MAX = axom::numerics::floating_point_limits<T>::max();

  const T STD_EPS = std::numeric_limits<T>::epsilon();
  const T STD_LOWEST = std::numeric_limits<T>::lowest();
  const T STD_MIN = std::numeric_limits<T>::min();
  const T STD_MAX = std::numeric_limits<T>::max();

  EXPECT_TRUE(axom::utilities::isNearlyEqual(LOWEST, STD_LOWEST, EPS));
  EXPECT_TRUE(axom::utilities::isNearlyEqual(MIN, STD_MIN, EPS));
  EXPECT_TRUE(axom::utilities::isNearlyEqual(MAX, STD_MAX, EPS));
  EXPECT_TRUE(axom::utilities::isNearlyEqual(EPS, STD_EPS, EPS));
}

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(numerics_floating_point_limits, consistency_with_standard_numeric_limits)
{
  check_type_limits<float>("float");
  check_type_limits<double>("double");
#if !defined(AXOM_DEVICE_CODE)
  check_type_limits<long double>("long double");
#endif
}
