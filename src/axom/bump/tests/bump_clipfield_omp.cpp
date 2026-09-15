// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/bump/tests/blueprint_testing_helpers.hpp"
#include "axom/bump/tests/bump_clipfield_impl.hpp"

axom::blueprint::testing::TestApplication TestApp;

TEST(bump_clipfield_omp, unique_omp) { test_unique<omp_exec>::test(); }
TEST(bump_clipfield_omp, onetet_omp)
{
  conduit::Node hostMesh;
  axom::blueprint::testing::data::make_one_tet(hostMesh);
  test_one_shape<omp_exec, axom::bump::views::TetShape<int>>(hostMesh, "one_tet");
}
TEST(bump_clipfield_omp, onepyr_omp)
{
  conduit::Node hostMesh;
  axom::blueprint::testing::data::make_one_pyr(hostMesh);
  test_one_shape<omp_exec, axom::bump::views::PyramidShape<int>>(hostMesh, "one_pyr");
}
TEST(bump_clipfield_omp, onewdg_omp)
{
  conduit::Node hostMesh;
  axom::blueprint::testing::data::make_one_wdg(hostMesh);
  test_one_shape<omp_exec, axom::bump::views::WedgeShape<int>>(hostMesh, "one_wdg");
}
TEST(bump_clipfield_omp, onehex_omp)
{
  conduit::Node hostMesh;
  axom::blueprint::testing::data::make_one_hex(hostMesh);
  test_one_shape<omp_exec, axom::bump::views::HexShape<int>>(hostMesh, "one_hex");
}
TEST(bump_clipfield_omp, uniform2d_omp) { braid2d_clip_test<omp_exec>("uniform", "uniform2d"); }
TEST(bump_clipfield_omp, rectilinear2d_omp)
{
  braid_rectilinear_clip_test<omp_exec, 2>("rectilinear2d");
}
TEST(bump_clipfield_omp, rectilinear3d_omp)
{
  braid_rectilinear_clip_test<omp_exec, 3>("rectilinear3d");
}
TEST(bump_clipfield_omp, strided_structured_2d_omp)
{
  conduit::Node options;
  options["field"] = "vert_vals";
  options["value"] = 6.5;
  options["inside"] = 1;
  options["outside"] = 1;
  strided_structured_clip_test<omp_exec, 2>("strided_structured_2d", options);

  options["selectedZones"].set(std::vector<axom::IndexType> {{0, 2, 3, 5}});
  strided_structured_clip_test<omp_exec, 2>("strided_structured_2d_sel", options);
}
TEST(bump_clipfield_omp, tet_omp)
{
  braid3d_clip_test<omp_exec, axom::bump::views::TetShape<int>>("tets", "tet");
}
TEST(bump_clipfield_omp, pyramid_omp)
{
  braid3d_clip_test<omp_exec, axom::bump::views::PyramidShape<int>>("pyramids", "pyr");
}
TEST(bump_clipfield_omp, wedge_omp)
{
  braid3d_clip_test<omp_exec, axom::bump::views::WedgeShape<int>>("wedges", "wdg");
}
TEST(bump_clipfield_omp, hex_omp)
{
  braid3d_clip_test<omp_exec, axom::bump::views::HexShape<int>>("hexs", "hex");
}
TEST(bump_clipfield_omp, mixed_omp) { braid3d_mixed_clip_test<omp_exec>("mixed"); }
TEST(bump_clipfield_omp, pointmerging_omp) { point_merge_test<omp_exec>::test(); }
TEST(bump_clipfield_omp, selectedzones_omp) { test_selectedzones<omp_exec>::test(); }

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
