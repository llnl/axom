// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/bump/tests/blueprint_testing_helpers.hpp"
#include "axom/mir/tests/mir_equiz3d_impl.hpp"

axom::blueprint::testing::TestApplication TestApp;

TEST(mir_equiz3d_omp, equiz_hex_unibuffer_omp)
{
  AXOM_ANNOTATE_SCOPE("equiz_hex_unibuffer_omp");
  braid3d_mat_test<omp_exec>("hexs", "unibuffer", "equiz_hex_unibuffer");
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
