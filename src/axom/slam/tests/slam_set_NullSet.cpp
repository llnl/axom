// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_set_NullSet.cpp
 *
 * Unit tests for the NullSet class in Slam.
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/slam/Set.hpp"
#include "axom/slam/NullSet.hpp"

TEST(slam_set_nullset, construct)
{
  axom::slam::NullSet<> ns;

  // Test function: isValid()
  EXPECT_TRUE(ns.isValid());

  // Test functions: empty() and size()
  EXPECT_TRUE(ns.empty());
  EXPECT_EQ(0, ns.size());

  // Test that the parent of a nullset is itself
  EXPECT_EQ(ns, *ns.parentSet());

  // A Slam NullSet is not a subset of another set
  EXPECT_FALSE(ns.isSubset());
}

TEST(slam_set_nullset, subscript_fails)
{
  SLIC_INFO("Testing subscript access on NullSet" << " -- code is expected to assert and die.");

  using SetPosition = axom::slam::DefaultPositionType;
  axom::slam::NullSet<SetPosition> n;

  EXPECT_EQ(n.size(), SetPosition()) << "size of null set is defined to be zero";

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode,
  // so this test will only fail in debug mode

  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED(n[0], "") << "subscript operator on null set asserts";
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();
  return result;
}
