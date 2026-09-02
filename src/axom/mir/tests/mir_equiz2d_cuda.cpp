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

TEST(mir_equiz2d_cuda, equiz_uniform_unibuffer_cuda)
{
  AXOM_ANNOTATE_SCOPE("equiz_uniform_unibuffer_cuda");
  test_equiz_uniform_unibuffer<cuda_exec>();
}

TEST(mir_equiz2d_cuda, equiz_polygonal_unibuffer_cuda)
{
  AXOM_ANNOTATE_SCOPE("equiz_polygonal_unibuffer_cuda");
  test_Polygonal_MIR<cuda_exec>::test("equiz_polygonal_unibuffer");
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
