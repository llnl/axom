// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/bump/tests/blueprint_testing_helpers.hpp"
#include "axom/mir/tests/mir_equiz2d_impl.hpp"

axom::blueprint::testing::TestApplication TestApp;

TEST(mir_equiz2d_hip, equiz_uniform_unibuffer_hip)
{
  AXOM_ANNOTATE_SCOPE("equiz_uniform_unibuffer_hip");
  test_equiz_uniform_unibuffer<hip_exec>();
}

TEST(mir_equiz2d_hip, equiz_polygonal_unibuffer_hip)
{
  AXOM_ANNOTATE_SCOPE("equiz_polygonal_unibuffer_hip");
  test_Polygonal_MIR<hip_exec>::test("equiz_polygonal_unibuffer");
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
