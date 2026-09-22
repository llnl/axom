// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/bump/tests/blueprint_testing_helpers.hpp"
#include "axom/mir/tests/mir_coupled3d_impl.hpp"

axom::blueprint::testing::TestApplication TestApp;

TEST(mir_coupled3d_hip, coupling_3d_hip)
{
  AXOM_ANNOTATE_SCOPE("coupling_3d_hip");
  const bool cleanMesh = true;
  test_coupling<hip_exec>::test("coupling_3d", cleanMesh);
}

TEST(mir_coupled3d_hip, coupling_3d_noclean_hip)
{
  AXOM_ANNOTATE_SCOPE("coupling_3d_noclean_hip");
  const bool cleanMesh = false;
  test_coupling<hip_exec>::test("coupling_3d", cleanMesh);
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
