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

TEST(bump_clipfield_cuda, unique_cuda) { test_unique<cuda_exec>::test(); }
TEST(bump_clipfield_cuda, onetet_cuda)
{
  conduit::Node hostMesh;
  axom::blueprint::testing::data::make_one_tet(hostMesh);
  test_one_shape<cuda_exec, axom::bump::views::TetShape<int>>(hostMesh, "one_tet");
}
TEST(bump_clipfield_cuda, onepyr_cuda)
{
  conduit::Node hostMesh;
  axom::blueprint::testing::data::make_one_pyr(hostMesh);
  test_one_shape<cuda_exec, axom::bump::views::PyramidShape<int>>(hostMesh, "one_pyr");
}
TEST(bump_clipfield_cuda, onewdg_cuda)
{
  conduit::Node hostMesh;
  axom::blueprint::testing::data::make_one_wdg(hostMesh);
  test_one_shape<cuda_exec, axom::bump::views::WedgeShape<int>>(hostMesh, "one_wdg");
}
TEST(bump_clipfield_cuda, onehex_cuda)
{
  conduit::Node hostMesh;
  axom::blueprint::testing::data::make_one_hex(hostMesh);
  test_one_shape<cuda_exec, axom::bump::views::HexShape<int>>(hostMesh, "one_hex");
}
TEST(bump_clipfield_cuda, uniform2d_cuda) { braid2d_clip_test<cuda_exec>("uniform", "uniform2d"); }
TEST(bump_clipfield_cuda, rectilinear2d_cuda)
{
  braid_rectilinear_clip_test<cuda_exec, 2>("rectilinear2d");
}
TEST(bump_clipfield_cuda, rectilinear3d_cuda)
{
  braid_rectilinear_clip_test<cuda_exec, 3>("rectilinear3d");
}
TEST(bump_clipfield_cuda, strided_structured_2d_cuda)
{
  conduit::Node options;
  options["field"] = "vert_vals";
  options["value"] = 6.5;
  options["inside"] = 1;
  options["outside"] = 1;
  strided_structured_clip_test<cuda_exec, 2>("strided_structured_2d", options);

  options["selectedZones"].set(std::vector<axom::IndexType> {{0, 2, 3, 5}});
  strided_structured_clip_test<cuda_exec, 2>("strided_structured_2d_sel", options);
}
TEST(bump_clipfield_cuda, tet_cuda)
{
  braid3d_clip_test<cuda_exec, axom::bump::views::TetShape<int>>("tets", "tet");
}
TEST(bump_clipfield_cuda, pyramid_cuda)
{
  braid3d_clip_test<cuda_exec, axom::bump::views::PyramidShape<int>>("pyramids", "pyr");
}
TEST(bump_clipfield_cuda, wedge_cuda)
{
  braid3d_clip_test<cuda_exec, axom::bump::views::WedgeShape<int>>("wedges", "wdg");
}
TEST(bump_clipfield_cuda, hex_cuda)
{
  braid3d_clip_test<cuda_exec, axom::bump::views::HexShape<int>>("hexs", "hex");
}
TEST(bump_clipfield_cuda, mixed_cuda) { braid3d_mixed_clip_test<cuda_exec>("mixed"); }
TEST(bump_clipfield_cuda, pointmerging_cuda) { point_merge_test<cuda_exec>::test(); }
TEST(bump_clipfield_cuda, selectedzones_cuda) { test_selectedzones<cuda_exec>::test(); }

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
