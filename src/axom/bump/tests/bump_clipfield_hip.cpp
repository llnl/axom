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

TEST(bump_clipfield_hip, unique_hip) { test_unique<hip_exec>::test(); }
TEST(bump_clipfield_hip, onetet_hip)
{
  conduit::Node hostMesh;
  axom::blueprint::testing::data::make_one_tet(hostMesh);
  test_one_shape<hip_exec, axom::bump::views::TetShape<int>>(hostMesh, "one_tet");
}
TEST(bump_clipfield_hip, onepyr_hip)
{
  conduit::Node hostMesh;
  axom::blueprint::testing::data::make_one_pyr(hostMesh);
  test_one_shape<hip_exec, axom::bump::views::PyramidShape<int>>(hostMesh, "one_pyr");
}
TEST(bump_clipfield_hip, onewdg_hip)
{
  conduit::Node hostMesh;
  axom::blueprint::testing::data::make_one_wdg(hostMesh);
  test_one_shape<hip_exec, axom::bump::views::WedgeShape<int>>(hostMesh, "one_wdg");
}
TEST(bump_clipfield_hip, onehex_hip)
{
  conduit::Node hostMesh;
  axom::blueprint::testing::data::make_one_hex(hostMesh);
  test_one_shape<hip_exec, axom::bump::views::HexShape<int>>(hostMesh, "one_hex");
}
TEST(bump_clipfield_hip, uniform2d_hip) { braid2d_clip_test<hip_exec>("uniform", "uniform2d"); }
TEST(bump_clipfield_hip, rectilinear2d_hip)
{
  braid_rectilinear_clip_test<hip_exec, 2>("rectilinear2d");
}
TEST(bump_clipfield_hip, rectilinear3d_hip)
{
  braid_rectilinear_clip_test<hip_exec, 3>("rectilinear3d");
}
TEST(bump_clipfield_hip, strided_structured_2d_hip)
{
  conduit::Node options;
  options["field"] = "vert_vals";
  options["value"] = 6.5;
  options["inside"] = 1;
  options["outside"] = 1;
  strided_structured_clip_test<hip_exec, 2>("strided_structured_2d", options);

  options["selectedZones"].set(std::vector<axom::IndexType> {{0, 2, 3, 5}});
  strided_structured_clip_test<hip_exec, 2>("strided_structured_2d_sel", options);
}
TEST(bump_clipfield_hip, tet_hip)
{
  braid3d_clip_test<hip_exec, axom::bump::views::TetShape<int>>("tets", "tet");
}
TEST(bump_clipfield_hip, pyramid_hip)
{
  braid3d_clip_test<hip_exec, axom::bump::views::PyramidShape<int>>("pyramids", "pyr");
}
TEST(bump_clipfield_hip, wedge_hip)
{
  braid3d_clip_test<hip_exec, axom::bump::views::WedgeShape<int>>("wedges", "wdg");
}
TEST(bump_clipfield_hip, hex_hip)
{
  braid3d_clip_test<hip_exec, axom::bump::views::HexShape<int>>("hexs", "hex");
}
TEST(bump_clipfield_hip, mixed_hip) { braid3d_mixed_clip_test<hip_exec>("mixed"); }
TEST(bump_clipfield_hip, pointmerging_hip) { point_merge_test<hip_exec>::test(); }
TEST(bump_clipfield_hip, selectedzones_hip) { test_selectedzones<hip_exec>::test(); }

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
