// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/mir/tests/mir_coupled_impl.hpp"
#include "axom/bump/tests/blueprint_testing_helpers.hpp"

axom::blueprint::testing::TestApplication TestApp;

TEST(mir_coupled_seq, coupling_2D_sz0_ss0_seq)
{
  AXOM_ANNOTATE_SCOPE("coupling_2D_sz0_ss0_seq");
  test_coupling<seq_exec>::test2D("coupling_2D_sz0_ss0", false, false);
}
TEST(mir_coupled_seq, coupling_2D_sz0_ss1_seq)
{
  AXOM_ANNOTATE_SCOPE("coupling_2D_sz0_ss1_seq");
  test_coupling<seq_exec>::test2D("coupling_2D_sz0_ss1", false, true);
}
TEST(mir_coupled_seq, coupling_2D_sz1_ss0_seq)
{
  AXOM_ANNOTATE_SCOPE("coupling_2D_sz1_ss0_seq");
  test_coupling<seq_exec>::test2D("coupling_2D_sz1_ss0", true, false);
}
TEST(mir_coupled_seq, coupling_2D_sz1_ss1_seq)
{
  AXOM_ANNOTATE_SCOPE("coupling_2D_sz1_ss1_seq");
  test_coupling<seq_exec>::test2D("coupling_2D_sz1_ss1", true, true);
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
