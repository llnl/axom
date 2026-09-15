// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/bump/tests/blueprint_testing_helpers.hpp"
#include "axom/mir/tests/mir_elvira3d_impl.hpp"

axom::blueprint::testing::TestApplication TestApp;

TEST(mir_elvira3d_omp, elvira3d_unibuffer_omp)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_omp");
  const bool selectZones = false;
  const bool pointMesh = false;
  test_Elvira3D<omp_exec>::test("elvira3d_unibuffer", selectZones, pointMesh);
}
TEST(mir_elvira3d_omp, elvira3d_unibuffer_sel_omp)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_sel_omp");
  const bool selectZones = true;
  const bool pointMesh = false;
  test_Elvira3D<omp_exec>::test("elvira3d_unibuffer_sel", selectZones, pointMesh);
}

TEST(mir_elvira3d_omp, elvira3d_unibuffer_pm_omp)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_pm_omp");
  const bool selectZones = false;
  const bool pointMesh = true;
  test_Elvira3D<omp_exec>::test("elvira3d_unibuffer_pm", selectZones, pointMesh);
}

TEST(mir_elvira3d_omp, elvira3d_unibuffer_sel_pm_omp)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_sel_pm_omp");
  const bool selectZones = true;
  const bool pointMesh = true;
  test_Elvira3D<omp_exec>::test("elvira3d_unibuffer_sel_pm", selectZones, pointMesh);
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
