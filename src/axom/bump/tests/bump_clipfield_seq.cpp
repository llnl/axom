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

TEST(bump_clipfield_seq, options)
{
  int nzones = 6;

  conduit::Node options;
  axom::bump::extraction::FieldOptions opts(options);

  options["field"] = "distance";
  EXPECT_EQ(opts.field(), options["field"].as_string());

  EXPECT_EQ(opts.value(), 0.);
  options["value"] = 2.5f;
  EXPECT_EQ(opts.value(), 2.5f);

  EXPECT_EQ(opts.topologyName("default"), "default");
  options["topologyName"] = "topo";
  EXPECT_EQ(opts.topologyName("default"), "topo");

  EXPECT_EQ(opts.coordsetName("default"), "default");
  options["coordsetName"] = "coords";
  EXPECT_EQ(opts.coordsetName("default"), "coords");

  EXPECT_EQ(opts.colorField(), "color");
  options["colorField"] = "custom_color";
  EXPECT_EQ(opts.colorField(), "custom_color");

  EXPECT_TRUE(opts.inside());
  options["inside"] = 1;
  EXPECT_TRUE(opts.inside());
  options["inside"] = 0;
  EXPECT_FALSE(opts.inside());

  EXPECT_FALSE(opts.outside());
  options["outside"] = 1;
  EXPECT_TRUE(opts.outside());
  options["outside"] = 0;
  EXPECT_FALSE(opts.outside());

  // The clip field has to be present
  conduit::Node n_fields;
  n_fields["distance/topology"] = "topo";
  n_fields["distance/association"] = "vertex";
  n_fields["distance/values"].set(std::vector<float> {0., 1., 2., 3.});

  // There are currently no fields in the options. fields should just return the clip field.
  std::map<std::string, std::string> fields;
  auto have_fields = opts.fields(fields);
  EXPECT_FALSE(have_fields);
  EXPECT_EQ(fields.size(), 0);

  // Add an empty fields node so we select NO fields.
  (void)options["fields"];
  have_fields = opts.fields(fields);
  EXPECT_TRUE(have_fields);
  EXPECT_EQ(fields.size(), 0);

  // Add some fields
  options["fields/distance"] = "distance";
  options["fields/source"] = "destination";
  options["fields/same"] = 1;
  have_fields = opts.fields(fields);
  EXPECT_TRUE(have_fields);
  EXPECT_EQ(fields.size(), 3);
  int i = 0;
  for(auto it = fields.begin(); it != fields.end(); it++, i++)
  {
    if(i == 0)
    {
      EXPECT_EQ(it->first, "distance");
      EXPECT_EQ(it->second, "distance");
    }
    else if(i == 1)
    {
      EXPECT_EQ(it->first, "same");
      EXPECT_EQ(it->second, "same");
    }
    else if(i == 2)
    {
      EXPECT_EQ(it->first, "source");
      EXPECT_EQ(it->second, "destination");
    }
  }

  // There are no "selectedZones" in the options. We should get nzones values from 0 onward.
  bump::SelectedZones<seq_exec> selectedZones(nzones, options);
  auto selectedZonesView = selectedZones.view();
  EXPECT_EQ(selectedZonesView.size(), 6);
  EXPECT_EQ(selectedZonesView[0], 0);
  EXPECT_EQ(selectedZonesView[1], 1);
  EXPECT_EQ(selectedZonesView[2], 2);
  EXPECT_EQ(selectedZonesView[3], 3);
  EXPECT_EQ(selectedZonesView[4], 4);
  EXPECT_EQ(selectedZonesView[5], 5);

  // Put some "selectedZones" in the options.
  options["selectedZones"].set(std::vector<axom::IndexType> {5, 4, 3});
  bump::SelectedZones<seq_exec> selectedZones2(nzones, options);
  selectedZonesView = selectedZones2.view();
  EXPECT_EQ(selectedZonesView.size(), 3);
  EXPECT_EQ(selectedZonesView[0], 3);
  EXPECT_EQ(selectedZonesView[1], 4);
  EXPECT_EQ(selectedZonesView[2], 5);
}
TEST(bump_clipfield_seq, blend_group_builder)
{
  using IndexType = axom::IndexType;
  using KeyType = std::uint64_t;

  /*

  We'll make 2 quads

  3      4      5
  *--8---*------*
  |      |      |
  |  6   9      |
  |      |      |
  *--7---*------*
  0      1      2

  */
  axom::Array<IndexType> blendGroups {{8, 5}};
  axom::Array<IndexType> blendGroupsLen {
    {/*zone 0*/ 4 + 1 + 1 + 1 + 1 + 2 + 2 + 2, /*zone 1*/ 1 + 1 + 1 + 1 + 2}};
  axom::Array<IndexType> blendGroupOffsets {{0, 8}};
  axom::Array<IndexType> blendOffsets {{0, blendGroupsLen[0]}};

  axom::Array<KeyType> blendNames {{/*zone 0*/ 6, 0, 1, 3, 4, 7, 8, 9, /*zone 1*/ 1, 2, 4, 5, 9}};
  axom::Array<IndexType> blendGroupSizes {
    {/*zone 0*/ 4, 1, 1, 1, 1, 2, 2, 2, /*zone 1*/ 1, 1, 1, 1, 2}};
  axom::Array<IndexType> blendGroupStart {
    {/*zone 0*/ 0, 4, 5, 6, 7, 8, 10, 12, /*zone 1*/ 13, 14, 15, 16, 18}};
  axom::Array<IndexType> blendIds {{
    /*zone 0*/
    0,
    1,
    2,
    3,  // 6 (bgname) // 0 (bgindex)
    0,  // 0          // 1
    1,  // 1          // 2
    3,  // 3          // 3
    4,  // 4          // 4
    0,
    1,  // 7          // 5
    3,
    4,  // 8          // 6
    1,
    4,  // 9          // 7
    /*zone 1*/
    1,  // 1          // 8
    2,  // 2          // 9
    4,  // 4          // 10
    5,  // 5          // 11
    1,
    4  // 9          // 12
  }};
  axom::Array<float> blendCoeff {{/*zone 0*/
                                  0.25,
                                  0.25,
                                  0.25,
                                  0.25,
                                  1.,
                                  1.,
                                  1.,
                                  1.,
                                  0.5,
                                  0.5,
                                  0.5,
                                  0.5,
                                  0.5,
                                  0.5,
                                  /*zone 1*/
                                  1.,
                                  1.,
                                  1.,
                                  1.,
                                  0.5,
                                  0.5}};
  axom::Array<KeyType> blendUniqueNames {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}};
  axom::Array<KeyType> blendUniqueIndices {{1, 2, 9, 3, 4, 11, 0, 5, 6, 7}};

  using NamingPolicyView = typename axom::bump::HashNaming<axom::IndexType>::View;

  axom::bump::extraction::BlendGroupBuilder<seq_exec, NamingPolicyView> builder;
  builder.setBlendGroupSizes(blendGroups.view(), blendGroupsLen.view());
  builder.setBlendGroupOffsets(blendOffsets.view(), blendGroupOffsets.view());
  builder.setBlendViews(blendNames.view(),
                        blendGroupSizes.view(),
                        blendGroupStart.view(),
                        blendIds.view(),
                        blendCoeff.view());

  //std::cout << "-------- zone 0 --------" << std::endl;
  auto z0 = builder.blendGroupsForZone(0);
  EXPECT_EQ(z0.numGroups(), 8);
  IndexType index = 0;
  for(IndexType i = 0; i < z0.numGroups(); i++, index++)
  {
    //z0.print(std::cout);
    EXPECT_EQ(z0.ids().size(), blendGroupSizes[index]);

    z0++;
  }

  //std::cout << "-------- zone 1 --------" << std::endl;
  auto z1 = builder.blendGroupsForZone(1);
  EXPECT_EQ(z1.numGroups(), 5);
  for(IndexType i = 0; i < z1.numGroups(); i++, index++)
  {
    //z1.print(std::cout);
    EXPECT_EQ(z1.ids().size(), blendGroupSizes[index]);
    z1++;
  }
}
TEST(bump_clipfield_seq, sort_values)
{
  constexpr int MaxSize = 15;
  for(int n = 1; n < MaxSize; n++)
  {
    for(int trial = 1; trial <= n; trial++)
    {
      auto values = makeUnsortedArray(n);
      axom::utilities::Sorting<int, MaxSize>::sort(values.data(), values.size());
      EXPECT_TRUE(increasing(values));
    }
  }
}
TEST(bump_clipfield_seq, make_name)
{
  axom::bump::HashNaming<int> naming;

  for(int n = 1; n < 14; n++)
  {
    // Make a set of scrambled ids.
    auto values = makeRandomArray(n);
    // Compute the name for that list of ids.
    auto name = naming.makeName(values.data(), n);

    for(int trial = 0; trial < 1000; trial++)
    {
      // Scramble the id list.
      auto values2 = permute(values);
      // Compute the name for that list of ids.
      auto name2 = naming.makeName(values2.data(), n);

      // The names for the 2 scrambled lists of numbers should be the same.
      EXPECT_EQ(name, name2);
    }
  }
}
TEST(bump_clipfield_seq, unique_seq) { test_unique<seq_exec>::test(); }
TEST(bump_clipfield_seq, onetet_seq)
{
  conduit::Node hostMesh;
  axom::blueprint::testing::data::make_one_tet(hostMesh);
  test_one_shape<seq_exec, axom::bump::views::TetShape<int>>(hostMesh, "one_tet");
}
TEST(bump_clipfield_seq, onepyr_seq)
{
  conduit::Node hostMesh;
  axom::blueprint::testing::data::make_one_pyr(hostMesh);
  test_one_shape<seq_exec, axom::bump::views::PyramidShape<int>>(hostMesh, "one_pyr");
}
TEST(bump_clipfield_seq, onewdg_seq)
{
  conduit::Node hostMesh;
  axom::blueprint::testing::data::make_one_wdg(hostMesh);
  test_one_shape<seq_exec, axom::bump::views::WedgeShape<int>>(hostMesh, "one_wdg");
}
TEST(bump_clipfield_seq, onehex_seq)
{
  conduit::Node hostMesh;
  axom::blueprint::testing::data::make_one_hex(hostMesh);
  test_one_shape<seq_exec, axom::bump::views::HexShape<int>>(hostMesh, "one_hex");
}
TEST(bump_clipfield_seq, uniform2d_seq) { braid2d_clip_test<seq_exec>("uniform", "uniform2d"); }
TEST(bump_clipfield_seq, rectilinear2d_seq)
{
  braid_rectilinear_clip_test<seq_exec, 2>("rectilinear2d");
}
TEST(bump_clipfield_seq, rectilinear3d_seq)
{
  braid_rectilinear_clip_test<seq_exec, 3>("rectilinear3d");
}
TEST(bump_clipfield_seq, strided_structured_2d_seq)
{
  conduit::Node options;
  options["field"] = "vert_vals";
  options["value"] = 6.5;
  options["inside"] = 1;
  options["outside"] = 1;
  strided_structured_clip_test<seq_exec, 2>("strided_structured_2d", options);

  options["selectedZones"].set(std::vector<axom::IndexType> {{0, 2, 3, 5}});
  strided_structured_clip_test<seq_exec, 2>("strided_structured_2d_sel", options);
}
TEST(bump_clipfield_seq, tet_seq)
{
  braid3d_clip_test<seq_exec, axom::bump::views::TetShape<int>>("tets", "tet");
}
TEST(bump_clipfield_seq, pyramid_seq)
{
  braid3d_clip_test<seq_exec, axom::bump::views::PyramidShape<int>>("pyramids", "pyr");
}
TEST(bump_clipfield_seq, wedge_seq)
{
  braid3d_clip_test<seq_exec, axom::bump::views::WedgeShape<int>>("wedges", "wdg");
}
TEST(bump_clipfield_seq, hex_seq)
{
  braid3d_clip_test<seq_exec, axom::bump::views::HexShape<int>>("hexs", "hex");
}
TEST(bump_clipfield_seq, mixed_seq) { braid3d_mixed_clip_test<seq_exec>("mixed"); }
TEST(bump_clipfield_seq, pointmerging_seq) { point_merge_test<seq_exec>::test(); }
TEST(bump_clipfield_seq, selectedzones_seq) { test_selectedzones<seq_exec>::test(); }

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
