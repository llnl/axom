// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"
#include "axom/slic.hpp"

#include "axom/primal/operators/detail/intersect_impl.hpp"

namespace primal = axom::primal;

TEST(primal_intersection_impl, fuzzy_comparisons)
{
  constexpr double eps = 0.1;

  SLIC_INFO("This test demonstrates the fuzzy comparison"
            << " operators used in quest's intersection tests"
            << " using a large epsilon (" << eps << ")");

  // Testing 'greater than' function
  // Note: Excludes the eps window
  EXPECT_TRUE(primal::detail::isGt(1.15, 1.0, eps));
  // --
  EXPECT_FALSE(primal::detail::isGt(1.05, 1.0, eps));
  EXPECT_FALSE(primal::detail::isGt(1.0, 1.0, eps));
  EXPECT_FALSE(primal::detail::isGt(0.95, 1.0, eps));
  EXPECT_FALSE(primal::detail::isGt(0.85, 1.0, eps));

  // Testing 'less than or equal' function
  // Note: Includes eps window
  // Note: isLeq() has opposite value of isGt()
  EXPECT_FALSE(primal::detail::isLeq(1.15, 1.0, eps));
  // --
  EXPECT_TRUE(primal::detail::isLeq(1.05, 1.0, eps));
  EXPECT_TRUE(primal::detail::isLeq(1.0, 1.0, eps));
  EXPECT_TRUE(primal::detail::isLeq(0.95, 1.0, eps));
  EXPECT_TRUE(primal::detail::isLeq(0.85, 1.0, eps));

  // Testing 'less than' function
  // Note: Excludes eps window
  EXPECT_FALSE(primal::detail::isLt(1.15, 1.0, eps));
  EXPECT_FALSE(primal::detail::isLt(1.05, 1.0, eps));
  EXPECT_FALSE(primal::detail::isLt(1.0, 1.0, eps));
  EXPECT_FALSE(primal::detail::isLt(0.95, 1.0, eps));
  // --
  EXPECT_TRUE(primal::detail::isLt(0.85, 1.0, eps));

  // Testing 'greater than or equal' function
  // Note: Includes eps window
  // Note: isGeq() has opposite value of isLt()
  EXPECT_TRUE(primal::detail::isGeq(1.15, 1.0, eps));
  EXPECT_TRUE(primal::detail::isGeq(1.05, 1.0, eps));
  EXPECT_TRUE(primal::detail::isGeq(1.0, 1.0, eps));
  EXPECT_TRUE(primal::detail::isGeq(0.95, 1.0, eps));
  // --
  EXPECT_FALSE(primal::detail::isGeq(0.85, 1.0, eps));

  EXPECT_FALSE(primal::detail::isLpeq(1.15, 1.0, true, eps));
  // --
  EXPECT_TRUE(primal::detail::isLpeq(1.05, 1.0, true, eps));
  EXPECT_TRUE(primal::detail::isLpeq(1.0, 1.0, true, eps));
  EXPECT_TRUE(primal::detail::isLpeq(0.95, 1.0, true, eps));
  EXPECT_TRUE(primal::detail::isLpeq(0.85, 1.0, true, eps));

  EXPECT_TRUE(primal::detail::isGpeq(1.15, 1.0, true, eps));
  EXPECT_TRUE(primal::detail::isGpeq(1.05, 1.0, true, eps));
  EXPECT_TRUE(primal::detail::isGpeq(1.0, 1.0, true, eps));
  EXPECT_TRUE(primal::detail::isGpeq(0.95, 1.0, true, eps));
  // --
  EXPECT_FALSE(primal::detail::isGpeq(0.85, 1.0, true, eps));

  EXPECT_FALSE(primal::detail::isLpeq(1.15, 1.0, false, eps));
  EXPECT_FALSE(primal::detail::isLpeq(1.05, 1.0, false, eps));
  EXPECT_FALSE(primal::detail::isLpeq(1.0, 1.0, false, eps));
  EXPECT_FALSE(primal::detail::isLpeq(0.95, 1.0, false, eps));
  // --
  EXPECT_TRUE(primal::detail::isLpeq(0.85, 1.0, false, eps));

  EXPECT_TRUE(primal::detail::isGpeq(1.15, 1.0, false, eps));
  // --
  EXPECT_FALSE(primal::detail::isGpeq(1.05, 1.0, false, eps));
  EXPECT_FALSE(primal::detail::isGpeq(1.0, 1.0, false, eps));
  EXPECT_FALSE(primal::detail::isGpeq(0.95, 1.0, false, eps));
  EXPECT_FALSE(primal::detail::isGpeq(0.85, 1.0, false, eps));
}

TEST(primal_intersection_impl, zero_count)
{
  constexpr double eps = 0.1;

  int expectedCount = 0;

  EXPECT_EQ(primal::detail::countZeros(1, 1, 1, eps), expectedCount);
  EXPECT_EQ(primal::detail::countZeros(-1, 1, -5, eps), expectedCount);
  EXPECT_EQ(primal::detail::countZeros(0.15, -0.15, 1.3, eps), expectedCount);

  expectedCount = 1;

  EXPECT_EQ(primal::detail::countZeros(-1, 1, 0, eps), expectedCount);
  EXPECT_EQ(primal::detail::countZeros(-0.2, -0.05, 3, eps), expectedCount);
  EXPECT_EQ(primal::detail::countZeros(0.07, -0.5, 3, eps), expectedCount);

  expectedCount = 2;

  EXPECT_EQ(primal::detail::countZeros(1, 0.05, 0, eps), expectedCount);
  EXPECT_EQ(primal::detail::countZeros(0.05, -0.02, -0.12, eps), expectedCount);
  EXPECT_EQ(primal::detail::countZeros(0.07, -0.5, 0, eps), expectedCount);

  expectedCount = 3;

  EXPECT_EQ(primal::detail::countZeros(0.01, 0.05, 0, eps), expectedCount);
  EXPECT_EQ(primal::detail::countZeros(-0.01, -0.02, 0.05, eps), expectedCount);

  EXPECT_TRUE(primal::detail::oneZeroOthersMatch(0, 1, 0.2, eps));
  EXPECT_TRUE(primal::detail::oneZeroOthersMatch(-1, 0.02, -9, eps));
  EXPECT_TRUE(primal::detail::oneZeroOthersMatch(-1, -0.9, -0.05, eps));
  EXPECT_FALSE(primal::detail::oneZeroOthersMatch(0, -0.2, 0.2, eps));
  EXPECT_FALSE(primal::detail::oneZeroOthersMatch(-1, 0.02, 5, eps));
  EXPECT_FALSE(primal::detail::oneZeroOthersMatch(-1, 0.02, 0.05, eps));
  EXPECT_FALSE(primal::detail::oneZeroOthersMatch(-1, 2, 5, eps));
  EXPECT_FALSE(primal::detail::oneZeroOthersMatch(-0.01, 0.05, 0, eps));

  EXPECT_TRUE(primal::detail::twoZeros(0, 0, 1, eps));
  EXPECT_TRUE(primal::detail::twoZeros(0.05, -0.05, 1, eps));
  EXPECT_TRUE(primal::detail::twoZeros(-0.05, -1.5, 0.05, eps));
  EXPECT_TRUE(primal::detail::twoZeros(1, 0.05, 0.01, eps));
  EXPECT_FALSE(primal::detail::twoZeros(0, 0.05, 0.01, eps));
  EXPECT_FALSE(primal::detail::twoZeros(0, 0.05, 0.01, eps));
  EXPECT_FALSE(primal::detail::twoZeros(1, 0.05, 0.15, eps));
  EXPECT_FALSE(primal::detail::twoZeros(1, 5, 0.15, eps));
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
