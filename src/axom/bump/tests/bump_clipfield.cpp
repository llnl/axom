// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/bump.hpp"
#include "axom/bump/tests/blueprint_testing_data_helpers.hpp"
#include "axom/bump/tests/blueprint_testing_helpers.hpp"

#include <conduit/conduit_relay_io_blueprint.hpp>
#include <cmath>
#include <cstdlib>

namespace bump = axom::bump;
namespace views = axom::bump::views;
namespace utils = axom::bump::utilities;

std::string baselineDirectory()
{
  return pjoin(dataDirectory(), "bump", "regression", "bump_clipfield");
}

//------------------------------------------------------------------------------
// Global test application object.
axom::blueprint::testing::TestApplication TestApp;

//------------------------------------------------------------------------------
TEST(bump_clipfield, options)
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

TEST(bump_clipfield, blend_group_builder)
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

//------------------------------------------------------------------------------
template <typename ArrayType>
bool increasing(const ArrayType &arr)
{
  bool retval = true;
  for(size_t i = 1; i < arr.size(); i++) retval &= (arr[i] >= arr[i - 1]);
  return retval;
}

template <typename ArrayType>
bool decreasing(const ArrayType &arr)
{
  bool retval = true;
  for(size_t i = 1; i < arr.size(); i++) retval &= (arr[i] <= arr[i - 1]);
  return retval;
}

std::vector<int> permute(const std::vector<int> &input)
{
  std::vector<int> values, indices;
  std::vector<double> order;

  values.resize(input.size());
  indices.resize(input.size());
  order.resize(input.size());

  std::iota(indices.begin(), indices.end(), 0);
  for(size_t i = 0; i < input.size(); i++)
  {
    order[i] = axom::utilities::random_real(0., 1.);
  }
  std::sort(indices.begin(), indices.end(), [&](int a, int b) { return order[a] < order[b]; });
  for(size_t i = 0; i < input.size(); i++) values[i] = input[indices[i]];
  return values;
}

std::vector<int> makeUnsortedArray(int n)
{
  std::vector<int> values;
  values.resize(n);
  std::iota(values.begin(), values.end(), 0);
  return permute(values);
}

std::vector<int> makeRandomArray(int n)
{
  constexpr double largestId = static_cast<double>(1 << 28);
  std::vector<int> values;
  values.resize(n);
  for(int i = 0; i < n; i++)
  {
    values[i] = static_cast<int>(axom::utilities::random_real(0., largestId));
  }
  return values;
}
std::vector<double> makeRandomDoubleArray(int n)
{
  std::vector<double> values;
  values.resize(n);
  for(int i = 0; i < n; i++)
  {
    values[i] = axom::utilities::random_real(0., 1000.);
  }
  return values;
}

//------------------------------------------------------------------------------
TEST(bump_clipfield, sort_values)
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

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_unique
{
  static void test()
  {
    /*
    8---9---10--11
    |   |   |   |
    4---5---6---7
    |   |   |   |
    0---1---2---3
    */
    // _bump_utilities_unique_begin
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    axom::Array<int> ids {
      {0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 4, 5, 9, 8, 5, 6, 10, 9, 6, 7, 11, 10}};
    EXPECT_EQ(ids.size(), 24);
    EXPECT_EQ(ids.view().size(), 24);

    // host->device
    axom::Array<int> deviceIds(ids.size(), ids.size(), allocatorID);
    axom::copy(deviceIds.data(), ids.data(), sizeof(int) * ids.size());
    EXPECT_EQ(deviceIds.size(), 24);

    // Make unique ids.
    axom::Array<int> uIds;
    axom::Array<axom::IndexType> uIndices;
    bump::Unique<ExecSpace, int>::execute(deviceIds.view(), uIds, uIndices);
    // _bump_utilities_unique_end

    // device->host
    axom::Array<int> hostuIds(uIds.size());
    axom::Array<axom::IndexType> hostuIndices(uIndices.size());
    axom::copy(hostuIds.data(), uIds.data(), sizeof(int) * uIds.size());
    axom::copy(hostuIndices.data(), uIndices.data(), sizeof(axom::IndexType) * uIndices.size());

    // compare results
    EXPECT_EQ(hostuIds.size(), 12);
    EXPECT_EQ(hostuIndices.size(), 12);
    for(axom::IndexType i = 0; i < hostuIds.size(); i++)
    {
      EXPECT_EQ(hostuIds[i], i);
      EXPECT_EQ(hostuIds[i], ids[hostuIndices[i]]);
    }
  }
};

TEST(bump_clipfield, unique_seq) { test_unique<seq_exec>::test(); }
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_clipfield, unique_omp) { test_unique<omp_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_clipfield, unique_cuda) { test_unique<cuda_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_clipfield, unique_hip) { test_unique<hip_exec>::test(); }
#endif

//------------------------------------------------------------------------------
TEST(bump_clipfield, make_name)
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

//------------------------------------------------------------------------------
template <typename ExecSpace, typename ShapeType>
void test_one_shape(const conduit::Node &hostMesh, const std::string &name)
{
  using TopoView = axom::bump::views::UnstructuredTopologySingleShapeView<ShapeType>;
  using CoordsetView = axom::bump::views::ExplicitCoordsetView<float, 3>;

  // Copy mesh to device
  conduit::Node deviceMesh;
  utils::copy<ExecSpace>(deviceMesh, hostMesh);

  // _bump_utilities_clipfield_begin
  // Make views for the device mesh.
  conduit::Node &n_x = deviceMesh.fetch_existing("coordsets/coords/values/x");
  conduit::Node &n_y = deviceMesh.fetch_existing("coordsets/coords/values/y");
  conduit::Node &n_z = deviceMesh.fetch_existing("coordsets/coords/values/z");
  axom::ArrayView<float> xView(static_cast<float *>(n_x.data_ptr()),
                               n_x.dtype().number_of_elements());
  axom::ArrayView<float> yView(static_cast<float *>(n_y.data_ptr()),
                               n_y.dtype().number_of_elements());
  axom::ArrayView<float> zView(static_cast<float *>(n_z.data_ptr()),
                               n_z.dtype().number_of_elements());
  CoordsetView coordsetView(xView, yView, zView);

  conduit::Node &n_conn = deviceMesh.fetch_existing("topologies/topo/elements/connectivity");
  axom::ArrayView<int> connView(static_cast<int *>(n_conn.data_ptr()),
                                n_conn.dtype().number_of_elements());
  TopoView topoView(connView);

  // Clip the data
  conduit::Node deviceClipMesh, options;
  axom::bump::extraction::ClipField<ExecSpace, TopoView, CoordsetView> clipper(topoView,
                                                                               coordsetView);
  options["field"] = "distance";
  options["value"] = 0.;
  options["inside"] = 1;
  options["outside"] = 1;
  clipper.execute(deviceMesh, options, deviceClipMesh);
  // _bump_utilities_clipfield_end

  // Copy device->host
  conduit::Node hostClipMesh;
  utils::copy<seq_exec>(hostClipMesh, deviceClipMesh);

  TestApp.saveVisualization(name, hostClipMesh);

  // Handle baseline comparison.
  EXPECT_TRUE(TestApp.test<ExecSpace>(name, hostClipMesh));
}

template <typename ShapeType>
void test_one_shape_exec(const conduit::Node &hostMesh, const std::string &name)
{
  test_one_shape<seq_exec, ShapeType>(hostMesh, name);

#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  test_one_shape<omp_exec, ShapeType>(hostMesh, name);
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  test_one_shape<cuda_exec, ShapeType>(hostMesh, name);
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  test_one_shape<hip_exec, ShapeType>(hostMesh, name);
#endif
}

TEST(bump_clipfield, onetet)
{
  conduit::Node hostMesh;
  axom::blueprint::testing::data::make_one_tet(hostMesh);
  test_one_shape_exec<axom::bump::views::TetShape<int>>(hostMesh, "one_tet");
}

TEST(bump_clipfield, onepyr)
{
  conduit::Node hostMesh;
  axom::blueprint::testing::data::make_one_pyr(hostMesh);
  test_one_shape_exec<axom::bump::views::PyramidShape<int>>(hostMesh, "one_pyr");
}

TEST(bump_clipfield, onewdg)
{
  conduit::Node hostMesh;
  axom::blueprint::testing::data::make_one_wdg(hostMesh);
  test_one_shape_exec<axom::bump::views::WedgeShape<int>>(hostMesh, "one_wdg");
}

TEST(bump_clipfield, onehex)
{
  conduit::Node hostMesh;
  axom::blueprint::testing::data::make_one_hex(hostMesh);
  test_one_shape_exec<axom::bump::views::HexShape<int>>(hostMesh, "one_hex");
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
void braid2d_clip_test(const std::string &type, const std::string &name)
{
  using Indexing = axom::bump::views::StructuredIndexing<axom::IndexType, 2>;
  using TopoView = axom::bump::views::StructuredTopologyView<Indexing>;
  using CoordsetView = axom::bump::views::UniformCoordsetView<double, 2>;

  axom::StackArray<axom::IndexType, 2> dims {10, 10};
  axom::StackArray<axom::IndexType, 2> zoneDims {dims[0] - 1, dims[1] - 1};

  // Create the data
  conduit::Node hostMesh, deviceMesh;
  axom::blueprint::testing::data::braid(type, dims, hostMesh);
  utils::copy<ExecSpace>(deviceMesh, hostMesh);
  TestApp.saveVisualization(name + "_orig", hostMesh);

  // Create views
  axom::StackArray<double, 2> origin {0., 0.}, spacing {1., 1.};
  CoordsetView coordsetView(dims, origin, spacing);
  TopoView topoView(Indexing {zoneDims});

  // Create options to control the clipping.
  const std::string clipTopoName("cliptopo");
  conduit::Node options;
  options["field"] = "distance";
  options["inside"] = 1;
  options["outside"] = 1;
  options["topologyName"] = clipTopoName;
  options["coordsetName"] = "clipcoords";
  options["fields/braid"] = "new_braid";
  options["fields/radial"] = "new_radial";

  // Clip the data
  conduit::Node deviceClipMesh;
  axom::bump::extraction::ClipField<ExecSpace, TopoView, CoordsetView> clipper(topoView,
                                                                               coordsetView);
  clipper.execute(deviceMesh, options, deviceClipMesh);

  // Copy device->host
  conduit::Node hostClipMesh;
  utils::copy<seq_exec>(hostClipMesh, deviceClipMesh);

  TestApp.saveVisualization(name, hostClipMesh);

  // Handle baseline comparison.
  EXPECT_TRUE(TestApp.test<ExecSpace>(name, hostClipMesh));

  // Now, take the clipped mesh and clip it again using a mixed topology view.
  using ExpCoordsetView = axom::bump::views::ExplicitCoordsetView<double, 2>;
  const auto xView =
    utils::make_array_view<double>(deviceClipMesh.fetch_existing("coordsets/clipcoords/values/x"));
  const auto yView =
    utils::make_array_view<double>(deviceClipMesh.fetch_existing("coordsets/clipcoords/values/y"));
  ExpCoordsetView expCoordsetView(xView, yView);

  conduit::Node &n_device_topo = deviceClipMesh.fetch_existing("topologies/" + clipTopoName);
  const auto connView =
    utils::make_array_view<axom::IndexType>(n_device_topo.fetch_existing("elements/connectivity"));

  options["field"] = "new_braid";
  options["value"] = 1.;
  options["fields"].reset();
  options["fields/new_braid"] = "new_braid2";
  options["fields/color"] = "new_color";
  options["fields/new_radial"] = "new_radial2";

  const auto shape = n_device_topo.fetch_existing("elements/shape").as_string();

  conduit::Node deviceClipMesh2;
  if(shape == "polygonal")
  {
    const auto sizesView =
      utils::make_array_view<axom::IndexType>(n_device_topo.fetch_existing("elements/sizes"));
    const auto offsetsView =
      utils::make_array_view<axom::IndexType>(n_device_topo.fetch_existing("elements/offsets"));

    using PolyTopoView =
      views::UnstructuredTopologySingleShapeView<views::PolygonShape<axom::IndexType>>;
    PolyTopoView polyTopoView(connView, sizesView, offsetsView);

    // Clip the data
    axom::bump::extraction::ClipField<ExecSpace, PolyTopoView, ExpCoordsetView> polyClipper(
      polyTopoView,
      expCoordsetView);
    polyClipper.execute(deviceClipMesh, options, deviceClipMesh2);
  }
  else if(shape == "quad")
  {
    // Depending on optimizations, we might get a mesh with just quads.
    using QuadTopoView =
      axom::bump::views::UnstructuredTopologySingleShapeView<axom::bump::views::QuadShape<axom::IndexType>>;
    QuadTopoView quadTopoView(connView);

    // Clip the data
    axom::bump::extraction::ClipField<ExecSpace, QuadTopoView, ExpCoordsetView> quadClipper(
      quadTopoView,
      expCoordsetView);
    quadClipper.execute(deviceClipMesh, options, deviceClipMesh2);
  }
  else
  {
    std::cout << "The test got an unexpected shape " << shape << std::endl;
    FAIL();
  }

  // Copy device->host
  conduit::Node hostClipMesh2;
  utils::copy<seq_exec>(hostClipMesh2, deviceClipMesh2);

  TestApp.saveVisualization(name + "_clip2", hostClipMesh2);

  // Handle baseline comparison.
  EXPECT_TRUE(TestApp.test<ExecSpace>(name + "_clip2", hostClipMesh2));
}

TEST(bump_clipfield, uniform2d)
{
  braid2d_clip_test<seq_exec>("uniform", "uniform2d");

#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  braid2d_clip_test<omp_exec>("uniform", "uniform2d");
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  braid2d_clip_test<cuda_exec>("uniform", "uniform2d");
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  braid2d_clip_test<hip_exec>("uniform", "uniform2d");
#endif
}

//------------------------------------------------------------------------------
template <typename ExecSpace, int NDIMS>
void braid_rectilinear_clip_test(const std::string &name)
{
  using Indexing = axom::bump::views::StructuredIndexing<axom::IndexType, NDIMS>;
  using TopoView = axom::bump::views::StructuredTopologyView<Indexing>;

  axom::StackArray<axom::IndexType, NDIMS> dims, zoneDims;
  for(int i = 0; i < NDIMS; i++)
  {
    dims[i] = 10;
    zoneDims[i] = dims[i] - 1;
  }

  // Create the data
  conduit::Node hostMesh, deviceMesh;
  axom::blueprint::testing::data::braid("rectilinear", dims, hostMesh);
  TestApp.saveVisualization(name + "_orig", hostMesh);

  // host->device
  utils::copy<ExecSpace>(deviceMesh, hostMesh);

  // Create views
  auto coordsetView = axom::bump::views::make_rectilinear_coordset<double, NDIMS>::view(
    deviceMesh["coordsets/coords"]);
  using CoordsetView = decltype(coordsetView);
  TopoView topoView(Indexing {zoneDims});

  // Create options to control the clipping.
  conduit::Node options;
  options["field"] = "distance";
  options["inside"] = 1;
  options["outside"] = 1;

  // Clip the data
  conduit::Node deviceClipMesh;
  axom::bump::extraction::ClipField<ExecSpace, TopoView, CoordsetView> clipper(topoView,
                                                                               coordsetView);
  clipper.execute(deviceMesh, options, deviceClipMesh);

  // Copy device->host
  conduit::Node hostClipMesh;
  utils::copy<seq_exec>(hostClipMesh, deviceClipMesh);

  TestApp.saveVisualization(name, hostClipMesh);

  // Handle baseline comparison.
  EXPECT_TRUE(TestApp.test<ExecSpace>(name, hostClipMesh));
}

TEST(bump_clipfield, rectilinear2d)
{
  braid_rectilinear_clip_test<seq_exec, 2>("rectilinear2d");

#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  braid_rectilinear_clip_test<omp_exec, 2>("rectilinear2d");
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  braid_rectilinear_clip_test<cuda_exec, 2>("rectilinear2d");
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  braid_rectilinear_clip_test<hip_exec, 2>("rectilinear2d");
#endif
}

TEST(bump_clipfield, rectilinear3d)
{
  braid_rectilinear_clip_test<seq_exec, 3>("rectilinear3d");

#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  braid_rectilinear_clip_test<omp_exec, 3>("rectilinear3d");
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  braid_rectilinear_clip_test<cuda_exec, 3>("rectilinear3d");
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  braid_rectilinear_clip_test<hip_exec, 3>("rectilinear3d");
#endif
}

//------------------------------------------------------------------------------
template <typename ExecSpace, int NDIMS>
void strided_structured_clip_test(const std::string &name, const conduit::Node &options)
{
  // Create the data
  conduit::Node hostMesh, deviceMesh;
  axom::blueprint::testing::data::strided_structured<NDIMS>(hostMesh);
  TestApp.saveVisualization(name + "_orig", hostMesh);

  conduit::Node deviceOptions, deviceClipMesh, hostClipMesh;

  // host->device
  utils::copy<ExecSpace>(deviceMesh, hostMesh);
  utils::copy<ExecSpace>(deviceOptions, options);

  // Create views
  const conduit::Node &n_coordset = deviceMesh["coordsets/coords"];
  const conduit::Node &n_topo = deviceMesh["topologies/mesh"];
  auto coordsetView = axom::bump::views::make_explicit_coordset<double, 2>::view(n_coordset);
  auto topoView = axom::bump::views::make_strided_structured_topology<2>::view(n_topo);

  using CoordsetView = decltype(coordsetView);
  using TopoView = decltype(topoView);

  // Clip the data
  axom::bump::extraction::ClipField<ExecSpace, TopoView, CoordsetView> clipper(topoView,
                                                                               coordsetView);
  clipper.execute(deviceMesh, deviceOptions, deviceClipMesh);

  // device->host
  utils::copy<seq_exec>(hostClipMesh, deviceClipMesh);

  TestApp.saveVisualization(name, hostClipMesh);

  // Handle baseline comparison.
  EXPECT_TRUE(TestApp.test<ExecSpace>(name, hostClipMesh));
}

void strided_structured_clip_test_exec(const std::string &name, const conduit::Node &options)
{
  strided_structured_clip_test<seq_exec, 2>(name, options);

#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  strided_structured_clip_test<omp_exec, 2>(name, options);
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  strided_structured_clip_test<cuda_exec, 2>(name, options);
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  strided_structured_clip_test<hip_exec, 2>(name, options);
#endif
}

TEST(bump_clipfield, strided_structured_2d)
{
  // Create options to control the clipping.
  conduit::Node options;
  options["field"] = "vert_vals";
  options["value"] = 6.5;
  options["inside"] = 1;
  options["outside"] = 1;
  strided_structured_clip_test_exec("strided_structured_2d", options);

  // Clip strided structure on some selected zones.
  options["selectedZones"].set(std::vector<axom::IndexType> {{0, 2, 3, 5}});
  strided_structured_clip_test_exec("strided_structured_2d_sel", options);
}

//------------------------------------------------------------------------------
template <typename ExecSpace, typename ShapeType>
void braid3d_clip_test(const std::string &type, const std::string &name)
{
  using TopoView = axom::bump::views::UnstructuredTopologySingleShapeView<ShapeType>;
  using CoordsetView = axom::bump::views::ExplicitCoordsetView<double, 3>;

  // Create the data
  const axom::StackArray<axom::IndexType, 3> dims {10, 10, 10};
  conduit::Node hostMesh, deviceMesh;
  axom::blueprint::testing::data::braid(type, dims, hostMesh);
  utils::copy<ExecSpace>(deviceMesh, hostMesh);
  TestApp.saveVisualization(name + "_orig", hostMesh);

  // Create views
  conduit::Node &n_x = deviceMesh.fetch_existing("coordsets/coords/values/x");
  conduit::Node &n_y = deviceMesh.fetch_existing("coordsets/coords/values/y");
  conduit::Node &n_z = deviceMesh.fetch_existing("coordsets/coords/values/z");
  const axom::ArrayView<double> x(static_cast<double *>(n_x.data_ptr()),
                                  n_x.dtype().number_of_elements());
  const axom::ArrayView<double> y(static_cast<double *>(n_y.data_ptr()),
                                  n_y.dtype().number_of_elements());
  const axom::ArrayView<double> z(static_cast<double *>(n_z.data_ptr()),
                                  n_z.dtype().number_of_elements());
  CoordsetView coordsetView(x, y, z);

  conduit::Node &n_conn = deviceMesh.fetch_existing("topologies/mesh/elements/connectivity");
  const axom::ArrayView<int> conn(static_cast<int *>(n_conn.data_ptr()),
                                  n_conn.dtype().number_of_elements());
  TopoView topoView(conn);

  // Create options to control the clipping.
  conduit::Node options;
  options["field"] = "distance";
  options["inside"] = 1;
  options["outside"] = 0;

  // Clip the data
  conduit::Node deviceClipMesh;
  axom::bump::extraction::ClipField<ExecSpace, TopoView, CoordsetView> clipper(topoView,
                                                                               coordsetView);
  clipper.execute(deviceMesh, options, deviceClipMesh);

  // Copy device->host
  conduit::Node hostClipMesh;
  utils::copy<seq_exec>(hostClipMesh, deviceClipMesh);

  TestApp.saveVisualization(name, hostClipMesh);

  // Handle baseline comparison.
  EXPECT_TRUE(TestApp.test<ExecSpace>(name, hostClipMesh));
}

/// Execute the braid3d test for a single shape on multiple ExecSpaces
template <typename ShapeType>
void braid3d_clip_test_exec(const std::string &type, const std::string &name)
{
  braid3d_clip_test<seq_exec, ShapeType>(type, name);

#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  braid3d_clip_test<omp_exec, ShapeType>(type, name);
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  braid3d_clip_test<cuda_exec, ShapeType>(type, name);
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  braid3d_clip_test<hip_exec, ShapeType>(type, name);
#endif
}

TEST(bump_clipfield, tet)
{
  braid3d_clip_test_exec<axom::bump::views::TetShape<int>>("tets", "tet");
}

TEST(bump_clipfield, pyramid)
{
  braid3d_clip_test_exec<axom::bump::views::PyramidShape<int>>("pyramids", "pyr");
}

TEST(bump_clipfield, wedge)
{
  braid3d_clip_test_exec<axom::bump::views::WedgeShape<int>>("wedges", "wdg");
}

TEST(bump_clipfield, hex)
{
  braid3d_clip_test_exec<axom::bump::views::HexShape<int>>("hexs", "hex");
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
void braid3d_mixed_clip_test(const std::string &name)
{
  using CoordType = float;
  using ConnType = int;
  using TopoView = axom::bump::views::UnstructuredTopologyMixedShapeView<ConnType>;
  using CoordsetView = axom::bump::views::ExplicitCoordsetView<CoordType, 3>;

  // Create the data
  conduit::Node hostMesh, deviceMesh;
  axom::blueprint::testing::data::mixed3d(hostMesh);
  utils::copy<ExecSpace>(deviceMesh, hostMesh);
  TestApp.saveVisualization(name + "_orig", hostMesh);

  // Create views
  conduit::Node &n_x = deviceMesh.fetch_existing("coordsets/coords/values/x");
  conduit::Node &n_y = deviceMesh.fetch_existing("coordsets/coords/values/y");
  conduit::Node &n_z = deviceMesh.fetch_existing("coordsets/coords/values/z");
  const axom::ArrayView<CoordType> x(static_cast<CoordType *>(n_x.data_ptr()),
                                     n_x.dtype().number_of_elements());
  const axom::ArrayView<CoordType> y(static_cast<CoordType *>(n_y.data_ptr()),
                                     n_y.dtype().number_of_elements());
  const axom::ArrayView<CoordType> z(static_cast<CoordType *>(n_z.data_ptr()),
                                     n_z.dtype().number_of_elements());
  CoordsetView coordsetView(x, y, z);

  conduit::Node &n_device_topo = deviceMesh.fetch_existing("topologies/mesh");
  conduit::Node &n_conn = n_device_topo.fetch_existing("elements/connectivity");
  conduit::Node &n_shapes = n_device_topo.fetch_existing("elements/shapes");
  conduit::Node &n_sizes = n_device_topo.fetch_existing("elements/sizes");
  conduit::Node &n_offsets = n_device_topo.fetch_existing("elements/offsets");
  axom::ArrayView<ConnType> connView(static_cast<ConnType *>(n_conn.data_ptr()),
                                     n_conn.dtype().number_of_elements());
  axom::ArrayView<ConnType> shapesView(static_cast<ConnType *>(n_shapes.data_ptr()),
                                       n_shapes.dtype().number_of_elements());
  axom::ArrayView<ConnType> sizesView(static_cast<ConnType *>(n_sizes.data_ptr()),
                                      n_sizes.dtype().number_of_elements());
  axom::ArrayView<ConnType> offsetsView(static_cast<ConnType *>(n_offsets.data_ptr()),
                                        n_offsets.dtype().number_of_elements());

  // Make the shape map.
  axom::Array<axom::IndexType> values, ids;
  auto shapeMap = axom::bump::views::buildShapeMap(n_device_topo,
                                                   values,
                                                   ids,
                                                   axom::execution_space<ExecSpace>::allocatorID());

  TopoView topoView(connView, shapesView, sizesView, offsetsView, shapeMap);

  // Create options to control the clipping.
  conduit::Node options;
  options["field"] = "distance";
  options["value"] = 12.f;
  options["inside"] = 1;
  options["outside"] = 0;

  // Clip the data
  conduit::Node deviceClipMesh;
  axom::bump::extraction::ClipField<ExecSpace, TopoView, CoordsetView> clipper(topoView,
                                                                               coordsetView);
  clipper.execute(deviceMesh, options, deviceClipMesh);

  // Copy device->host
  conduit::Node hostClipMesh;
  utils::copy<seq_exec>(hostClipMesh, deviceClipMesh);

  TestApp.saveVisualization(name, hostClipMesh);

  // Handle baseline comparison.
  EXPECT_TRUE(TestApp.test<ExecSpace>(name, hostClipMesh));
}

TEST(bump_clipfield, mixed_seq) { braid3d_mixed_clip_test<seq_exec>("mixed"); }
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_clipfield, mixed_omp) { braid3d_mixed_clip_test<omp_exec>("mixed"); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_clipfield, mixed_cuda) { braid3d_mixed_clip_test<cuda_exec>("mixed"); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_clipfield, mixed_hip) { braid3d_mixed_clip_test<hip_exec>("mixed"); }
#endif

//------------------------------------------------------------------------------
template <typename Container1, typename Container2>
void compare_values(const Container1 &c1, const Container2 &c2)
{
  EXPECT_EQ(c1.size(), c2.number_of_elements());
  for(size_t i = 0; i < c1.size(); i++)
  {
    EXPECT_EQ(c1[i], c2[i]);
  }
}

template <typename ExecSpace>
struct point_merge_test
{
  static void create(conduit::Node &hostMesh)
  {
    hostMesh["coordsets/coords/type"] = "explicit";
    hostMesh["coordsets/coords/values/x"].set(
      std::vector<float> {{0., 1., 2., 0., 1., 2., 0., 1., 2.}});
    hostMesh["coordsets/coords/values/y"].set(
      std::vector<float> {{0., 0., 0., 1., 1., 1., 2., 2., 2.}});
    hostMesh["topologies/mesh/type"] = "unstructured";
    hostMesh["topologies/mesh/coordset"] = "coords";
    hostMesh["topologies/mesh/elements/shape"] = "quad";
    hostMesh["topologies/mesh/elements/connectivity"].set(
      std::vector<int> {{0, 1, 4, 3, 1, 2, 5, 4, 3, 4, 7, 6, 4, 5, 8, 7}});
    hostMesh["topologies/mesh/elements/sizes"].set(std::vector<int> {{4, 4, 4, 4}});
    hostMesh["topologies/mesh/elements/offsets"].set(std::vector<int> {{0, 4, 8, 12}});
    hostMesh["fields/clip/topology"] = "mesh";
    hostMesh["fields/clip/association"] = "vertex";
    hostMesh["fields/clip/values"].set(std::vector<float> {{1., 1., 0.5, 1., 1., 0., 0.5, 0., 0.5}});
  }

  static void test()
  {
    conduit::Node hostMesh;
    create(hostMesh);
    TestApp.saveVisualization("pointmerge_orig", hostMesh);

    // host->device
    conduit::Node deviceMesh;
    utils::copy<ExecSpace>(deviceMesh, hostMesh);

    // Set up views for the mesh.
    using CoordsetView = axom::bump::views::ExplicitCoordsetView<float, 2>;
    CoordsetView coordsetView(
      utils::make_array_view<float>(deviceMesh.fetch_existing("coordsets/coords/values/x")),
      utils::make_array_view<float>(deviceMesh.fetch_existing("coordsets/coords/values/y")));

    using TopologyView =
      axom::bump::views::UnstructuredTopologySingleShapeView<axom::bump::views::QuadShape<int>>;
    TopologyView topologyView(utils::make_array_view<int>(
      deviceMesh.fetch_existing("topologies/mesh/elements/connectivity")));

    // Clip
    conduit::Node options, deviceClipMesh;
    options["field"] = "clip";
    options["value"] = 0.5;
    using Clip = axom::bump::extraction::ClipField<ExecSpace, TopologyView, CoordsetView>;
    Clip clip(topologyView, coordsetView);
    clip.execute(deviceMesh, options, deviceClipMesh);

    // device->host
    conduit::Node hostClipMesh;
    utils::copy<axom::SEQ_EXEC>(hostClipMesh, deviceClipMesh);
    //printNode(hostClipMesh);

    TestApp.saveVisualization("pointmerge", hostClipMesh);

    // Check that the points were merged when making the new mesh.
    std::vector<float> x {{2.0, 2.0, 0.0, 1.0, 2.0, 1.5, 1.0}};
    std::vector<float> y {{0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 1.5}};
    EXPECT_EQ(x.size(), 7);
    EXPECT_EQ(y.size(), 7);
    for(size_t i = 0; i < x.size(); i++)
    {
      EXPECT_FLOAT_EQ(hostClipMesh["coordsets/coords/values/x"].as_float_accessor()[i], x[i]);
      EXPECT_FLOAT_EQ(hostClipMesh["coordsets/coords/values/y"].as_float_accessor()[i], y[i]);
    }

    // Check that we git tris and a pentagon.
    EXPECT_EQ(hostClipMesh["topologies/mesh/elements/shape"].as_string(), "polygonal");
    std::vector<int> sizes {{3, 3, 5}};
    std::vector<int> offsets {{0, 4, 8}};
    compare_values(sizes, hostClipMesh["topologies/mesh/elements/sizes"].as_int_accessor());
    compare_values(offsets, hostClipMesh["topologies/mesh/elements/offsets"].as_int_accessor());
  }
};

TEST(bump_clipfield, pointmerging_seq) { point_merge_test<seq_exec>::test(); }
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_clipfield, pointmerging_omp) { point_merge_test<omp_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_clipfield, pointmerging_cuda) { point_merge_test<cuda_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_clipfield, pointmerging_hip) { point_merge_test<hip_exec>::test(); }
#endif

//------------------------------------------------------------------------------

template <typename ExecSpace>
struct test_selectedzones
{
  static void test()
  {
    conduit::Node hostMesh;
    create(hostMesh);

    // host->device
    conduit::Node deviceMesh;
    utils::copy<ExecSpace>(deviceMesh, hostMesh);

    // Wrap the data in views.
    auto coordsetView = axom::bump::views::make_rectilinear_coordset<conduit::float64, 2>::view(
      deviceMesh["coordsets/coords"]);
    using CoordsetView = decltype(coordsetView);

    auto topologyView =
      axom::bump::views::make_rectilinear_topology<2>::view(deviceMesh["topologies/mesh"]);
    using TopologyView = decltype(topologyView);

    conduit::Node hostOptions;
    hostOptions["selectedZones"].set(std::vector<axom::IndexType> {{1, 3, 4, 5, 7}});
    hostOptions["inside"] = 1;
    hostOptions["outside"] = 1;
    hostOptions["field"] = "zero";

    conduit::Node deviceOptions, deviceResult;
    utils::copy<ExecSpace>(deviceOptions, hostOptions);

    axom::bump::extraction::ClipField<ExecSpace, TopologyView, CoordsetView> clip(topologyView,
                                                                                  coordsetView);
    clip.execute(deviceMesh, deviceOptions, deviceResult);

    // device->host
    conduit::Node hostResult;
    utils::copy<seq_exec>(hostResult, deviceResult);

    TestApp.saveVisualization("selectedzones1", hostResult);

    // Handle baseline comparison.
    EXPECT_TRUE(TestApp.test<ExecSpace>("selectedzones1", hostResult));

    //---------------------
    // Try a different clip
    hostOptions["outside"] = 0;
    hostOptions["field"] = "radial";
    hostOptions["value"] = 3.2;
    utils::copy<ExecSpace>(deviceOptions, hostOptions);
    deviceResult.reset();
    clip.execute(deviceMesh, deviceOptions, deviceResult);

    // device->host
    utils::copy<seq_exec>(hostResult, deviceResult);

    TestApp.saveVisualization("selectedzones2", hostResult);

    EXPECT_TRUE(TestApp.test<ExecSpace>("selectedzones2", hostResult));
  }

  static void create(conduit::Node &mesh)
  {
    /*
      12--13--14--15
      |   | x |   |
      8---9--10---11
      | x | x | x |   x=selected zones
      4---5---6---7
      |   | x |   |
      0---1---2---3
      */
    const char *yaml = R"xx(
coordsets:
  coords:
    type: rectilinear
    values:
      x: [0., 1., 2., 3.]
      y: [0., 1., 2., 3.]
topologies:
  mesh:
    type: rectilinear
    coordset: coords
fields:
  zero:
    topology: mesh
    association: vertex
    values: [0., 0., 0., 0.,  0., 0., 0., 0.,  0., 0., 0., 0.,  0., 0., 0., 0.]
  radial:
    topology: mesh
    association: vertex
    values: [0.0, 1.0, 2.0, 3.0, 1.0, 1.414, 2.236, 3.162, 2.0, 2.236, 2.828, 3.605, 3.0, 3.162, 3.605, 4.242]
)xx";

    mesh.parse(yaml);
  }
};

TEST(bump_clipfield, selectedzones_seq) { test_selectedzones<seq_exec>::test(); }
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_clipfield, selectedzones_omp) { test_selectedzones<omp_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_clipfield, selectedzones_cuda) { test_selectedzones<cuda_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_clipfield, selectedzones_hip) { test_selectedzones<hip_exec>::test(); }
#endif

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
