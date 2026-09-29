// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/bump/tests/blueprint_testing_helpers.hpp"
#include "axom/mir/tests/mir_elvira2d_impl.hpp"

axom::blueprint::testing::TestApplication TestApp;

TEST(mir_elvira2d_seq, options)
{
  conduit::Node n_options;

  axom::mir::ELVIRAOptions opts(n_options);
  EXPECT_FALSE(opts.pointmesh());
  EXPECT_FALSE(opts.plane());

  n_options["plane"] = 1;
  n_options["pointmesh"] = 1;

  EXPECT_TRUE(opts.pointmesh());
  EXPECT_TRUE(opts.plane());
}
//------------------------------------------------------------------------------
TEST(mir_elvira2d_seq, elvira_uniform_unibuffer_seq)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_seq");
  const bool selectZones = false;
  const bool pointMesh = false;
  braid2d_mat_test<seq_exec>::test("uniform",
                                   "unibuffer",
                                   "elvira_uniform_unibuffer",
                                   selectZones,
                                   pointMesh);
  // Run 2 domain example
  {
    const int nDomains = 2;
    const bool cleanMats = false;
    braid2d_mat_test<seq_exec>::test("uniform",
                                     "unibuffer",
                                     "elvira_uniform_unibuffer",
                                     selectZones,
                                     pointMesh,
                                     cleanMats,
                                     nDomains);
  }

  // Run clean mats example.
  {
    const bool cleanMats = true;
    braid2d_mat_test<seq_exec>::test("uniform",
                                     "unibuffer",
                                     "elvira_uniform_unibuffer_clean",
                                     selectZones,
                                     pointMesh,
                                     cleanMats);
  }
}

TEST(mir_elvira2d_seq, elvira_uniform_unibuffer_sel_seq)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_sel_seq");
  const bool selectZones = true;
  const bool pointMesh = false;
  braid2d_mat_test<seq_exec>::test("uniform",
                                   "unibuffer",
                                   "elvira_uniform_unibuffer_sel",
                                   selectZones,
                                   pointMesh);

  // Run clean mats example with selected zones.
  {
    const bool cleanMats = true;
    braid2d_mat_test<seq_exec>::test("uniform",
                                     "unibuffer",
                                     "elvira_uniform_unibuffer_sel_clean",
                                     selectZones,
                                     pointMesh,
                                     cleanMats);
  }
}

TEST(mir_elvira2d_seq, elvira_uniform_unibuffer_seq_pm)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_pm_seq");
  const bool selectZones = false;
  const bool pointMesh = true;
  braid2d_mat_test<seq_exec>::test("uniform",
                                   "unibuffer",
                                   "elvira_uniform_unibuffer_pm",
                                   selectZones,
                                   pointMesh);
}

TEST(mir_elvira2d_seq, elvira_uniform_unibuffer_sel_pm_seq)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_sel_pm_seq");
  const bool selectZones = true;
  const bool pointMesh = true;
  braid2d_mat_test<seq_exec>::test("uniform",
                                   "unibuffer",
                                   "elvira_uniform_unibuffer_sel_pm",
                                   selectZones,
                                   pointMesh);
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
