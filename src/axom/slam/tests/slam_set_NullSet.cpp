// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
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

#if AXOM_SLIC_DEBUG_MACROS_ENABLED
  // This check is only active when Slic debug macros are enabled.

  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED(n[0], "") << "subscript operator on null set asserts";
#else
  SLIC_INFO("Skipped assertion failure check because Slic debug macros are disabled.");
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
