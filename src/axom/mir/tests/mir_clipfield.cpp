// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/mir.hpp"
#include "axom/mir/tests/mir_testing_data_helpers.hpp"

#include <conduit/conduit_relay_io_blueprint.hpp>
#include <cmath>
#include <cstdlib>

namespace bputils = axom::mir::utilities::blueprint;

//------------------------------------------------------------------------------

// Uncomment to make Conduit errors hang.
//#define DEBUGGING_TEST_CASES

// Uncomment to generate baselines
//#define AXOM_TESTING_GENERATE_BASELINES

// Uncomment to save visualization files for debugging (when making baselines)
//#define AXOM_TESTING_SAVE_VISUALIZATION

#include "axom/mir/tests/mir_testing_helpers.hpp"

std::string baselineDirectory()
{
  return pjoin(dataDirectory(), "mir", "regression", "mir_clipfield");
}
//------------------------------------------------------------------------------
TEST(mir_clipfield, options)
{
  int nzones = 6;

  conduit::Node options;
  axom::mir::clipping::ClipOptions opts(options);

  options["clipField"] = "distance";
  EXPECT_EQ(opts.clipField(), options["clipField"].as_string());

  EXPECT_EQ(opts.clipValue(), 0.);
  options["clipValue"] = 2.5f;
  EXPECT_EQ(opts.clipValue(), 2.5f);

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
  bputils::SelectedZones<seq_exec> selectedZones(nzones, options);
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
  bputils::SelectedZones<seq_exec> selectedZones2(nzones, options);
  selectedZonesView = selectedZones2.view();
  EXPECT_EQ(selectedZonesView.size(), 3);
  EXPECT_EQ(selectedZonesView[0], 3);
  EXPECT_EQ(selectedZonesView[1], 4);
  EXPECT_EQ(selectedZonesView[2], 5);
}

TEST(mir_clipfield, blend_group_builder)
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

  using NamingPolicyView = typename axom::mir::utilities::HashNaming<axom::IndexType>::View;

  axom::mir::clipping::BlendGroupBuilder<seq_exec, NamingPolicyView> builder;
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
TEST(mir_clipfield, sort_values)
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
    // _mir_utilities_unique_begin
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
    axom::mir::utilities::Unique<ExecSpace, int>::execute(deviceIds.view(), uIds, uIndices);
    // _mir_utilities_unique_end

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

TEST(mir_clipfield, unique_seq) { test_unique<seq_exec>::test(); }
#if defined(AXOM_USE_OPENMP)
TEST(mir_clipfield, unique_omp) { test_unique<omp_exec>::test(); }
#endif
#if defined(AXOM_USE_CUDA)
TEST(mir_clipfield, unique_cuda) { test_unique<cuda_exec>::test(); }
#endif
#if defined(AXOM_USE_HIP)
TEST(mir_clipfield, unique_hip) { test_unique<hip_exec>::test(); }
#endif

//------------------------------------------------------------------------------
TEST(mir_clipfield, make_name)
{
  axom::mir::utilities::HashNaming<int> naming;

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
  using TopoView = axom::mir::views::UnstructuredTopologySingleShapeView<ShapeType>;
  using CoordsetView = axom::mir::views::ExplicitCoordsetView<float, 3>;

  // Copy mesh to device
  conduit::Node deviceMesh;
  bputils::copy<ExecSpace>(deviceMesh, hostMesh);

  // _mir_utilities_clipfield_begin
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
  axom::mir::clipping::ClipField<ExecSpace, TopoView, CoordsetView> clipper(topoView, coordsetView);
  options["clipField"] = "distance";
  options["clipValue"] = 0.;
  options["inside"] = 1;
  options["outside"] = 1;
  clipper.execute(deviceMesh, options, deviceClipMesh);
  // _mir_utilities_clipfield_end

  // Copy device->host
  conduit::Node hostClipMesh;
  bputils::copy<seq_exec>(hostClipMesh, deviceClipMesh);

  // Handle baseline comparison.
  std::string baselineName(yamlRoot(name));
  const auto paths = baselinePaths<ExecSpace>();
#if defined(AXOM_TESTING_GENERATE_BASELINES)
  saveBaseline(paths, baselineName, hostClipMesh);
#else
  EXPECT_TRUE(compareBaseline(paths, baselineName, hostClipMesh));
#endif
}

template <typename ShapeType>
void test_one_shape_exec(const conduit::Node &hostMesh, const std::string &name)
{
  test_one_shape<seq_exec, ShapeType>(hostMesh, name);

#if defined(AXOM_USE_OPENMP)
  test_one_shape<omp_exec, ShapeType>(hostMesh, name);
#endif

#if defined(AXOM_USE_CUDA)
  test_one_shape<cuda_exec, ShapeType>(hostMesh, name);
#endif

#if defined(AXOM_USE_HIP)
  test_one_shape<hip_exec, ShapeType>(hostMesh, name);
#endif
}

TEST(mir_clipfield, onetet)
{
  conduit::Node hostMesh;
  axom::mir::testing::data::make_one_tet(hostMesh);
  test_one_shape_exec<axom::mir::views::TetShape<int>>(hostMesh, "one_tet");
}

TEST(mir_clipfield, onepyr)
{
  conduit::Node hostMesh;
  axom::mir::testing::data::make_one_pyr(hostMesh);
  test_one_shape_exec<axom::mir::views::PyramidShape<int>>(hostMesh, "one_pyr");
}

TEST(mir_clipfield, onewdg)
{
  conduit::Node hostMesh;
  axom::mir::testing::data::make_one_wdg(hostMesh);
  test_one_shape_exec<axom::mir::views::WedgeShape<int>>(hostMesh, "one_wdg");
}

TEST(mir_clipfield, onehex)
{
  conduit::Node hostMesh;
  axom::mir::testing::data::make_one_hex(hostMesh);
  test_one_shape_exec<axom::mir::views::HexShape<int>>(hostMesh, "one_hex");
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
void braid2d_clip_test(const std::string &type, const std::string &name)
{
  using Indexing = axom::mir::views::StructuredIndexing<axom::IndexType, 2>;
  using TopoView = axom::mir::views::StructuredTopologyView<Indexing>;
  using CoordsetView = axom::mir::views::UniformCoordsetView<double, 2>;

  axom::StackArray<axom::IndexType, 2> dims {10, 10};
  axom::StackArray<axom::IndexType, 2> zoneDims {dims[0] - 1, dims[1] - 1};

  // Create the data
  conduit::Node hostMesh, deviceMesh;
  axom::mir::testing::data::braid(type, dims, hostMesh);
  bputils::copy<ExecSpace>(deviceMesh, hostMesh);
#if defined(AXOM_TESTING_SAVE_VISUALIZATION) && defined(AXOM_USE_HDF5)
  conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig", "hdf5");
#endif
#if defined(DEBUGGING_TEST_CASES)
  std::cout << "---------------------------------- input mesh "
               "-------------------------------------"
            << std::endl;
  printNode(hostMesh);
  std::cout << "---------------------------------------------------------------"
               "--------------------"
            << std::endl;
#endif

  // Create views
  axom::StackArray<double, 2> origin {0., 0.}, spacing {1., 1.};
  CoordsetView coordsetView(dims, origin, spacing);
  TopoView topoView(Indexing {zoneDims});

  // Create options to control the clipping.
  const std::string clipTopoName("cliptopo");
  conduit::Node options;
  options["clipField"] = "distance";
  options["inside"] = 1;
  options["outside"] = 1;
  options["topologyName"] = clipTopoName;
  options["coordsetName"] = "clipcoords";
  options["fields/braid"] = "new_braid";
  options["fields/radial"] = "new_radial";

  // Clip the data
  conduit::Node deviceClipMesh;
  axom::mir::clipping::ClipField<ExecSpace, TopoView, CoordsetView> clipper(topoView, coordsetView);
  clipper.execute(deviceMesh, options, deviceClipMesh);

  // Copy device->host
  conduit::Node hostClipMesh;
  bputils::copy<seq_exec>(hostClipMesh, deviceClipMesh);

#if defined(DEBUGGING_TEST_CASES)
  std::cout << "---------------------------------- clipped uniform "
               "----------------------------------"
            << std::endl;
  printNode(hostClipMesh);
  std::cout << "---------------------------------------------------------------"
               "----------------------"
            << std::endl;
#endif

  // Handle baseline comparison.
  {
    std::string baselineName(yamlRoot(name));
    const auto paths = baselinePaths<ExecSpace>();
#if defined(AXOM_TESTING_GENERATE_BASELINES)
    saveBaseline(paths, baselineName, hostClipMesh);
#else
    EXPECT_TRUE(compareBaseline(paths, baselineName, hostClipMesh));
#endif
  }

  // Now, take the clipped mesh and clip it again using a mixed topology view.
  using ExpCoordsetView = axom::mir::views::ExplicitCoordsetView<double, 2>;
  const auto xView = bputils::make_array_view<double>(
    deviceClipMesh.fetch_existing("coordsets/clipcoords/values/x"));
  const auto yView = bputils::make_array_view<double>(
    deviceClipMesh.fetch_existing("coordsets/clipcoords/values/y"));
  ExpCoordsetView expCoordsetView(xView, yView);

  conduit::Node &n_device_topo = deviceClipMesh.fetch_existing("topologies/" + clipTopoName);
  const auto connView = bputils::make_array_view<axom::IndexType>(
    n_device_topo.fetch_existing("elements/connectivity"));

  options["clipField"] = "new_braid";
  options["clipValue"] = 1.;
  options["fields"].reset();
  options["fields/new_braid"] = "new_braid2";
  options["fields/color"] = "new_color";
  options["fields/new_radial"] = "new_radial2";

  conduit::Node deviceClipMixedMesh;
  if(n_device_topo.has_path("elements/shape") &&
     n_device_topo.fetch_existing("elements/shape").as_string() == "mixed")
  {
    auto shapesView =
      bputils::make_array_view<axom::IndexType>(n_device_topo.fetch_existing("elements/shapes"));
    const auto sizesView =
      bputils::make_array_view<axom::IndexType>(n_device_topo.fetch_existing("elements/sizes"));
    const auto offsetsView =
      bputils::make_array_view<axom::IndexType>(n_device_topo.fetch_existing("elements/offsets"));

    // Make the shape map.
    volatile int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    axom::Array<axom::IndexType> values, ids;
    auto shapeMap = axom::mir::views::buildShapeMap(n_device_topo, values, ids, allocatorID);

    using MixedTopoView = axom::mir::views::UnstructuredTopologyMixedShapeView<axom::IndexType>;
    MixedTopoView mixedTopoView(connView, shapesView, sizesView, offsetsView, shapeMap);

    // Clip the data
    axom::mir::clipping::ClipField<ExecSpace, MixedTopoView, ExpCoordsetView> mixedClipper(
      mixedTopoView,
      expCoordsetView);
    mixedClipper.execute(deviceClipMesh, options, deviceClipMixedMesh);
  }
  else
  {
    // Depending on optimizations, we might get a mesh with just quads.
    using QuadTopoView =
      axom::mir::views::UnstructuredTopologySingleShapeView<axom::mir::views::QuadShape<axom::IndexType>>;
    QuadTopoView quadTopoView(connView);

    // Clip the data
    axom::mir::clipping::ClipField<ExecSpace, QuadTopoView, ExpCoordsetView> quadClipper(
      quadTopoView,
      expCoordsetView);
    quadClipper.execute(deviceClipMesh, options, deviceClipMixedMesh);
  }

  // Copy device->host
  conduit::Node hostClipMixedMesh;
  bputils::copy<seq_exec>(hostClipMixedMesh, deviceClipMixedMesh);
#if defined(DEBUGGING_TEST_CASES)
  std::cout << "---------------------------------- clipped mixed "
               "----------------------------------"
            << std::endl;
  printNode(hostClipMixedMesh);
  std::cout << "---------------------------------------------------------------"
               "--------------------"
            << std::endl;
#endif

  // Handle baseline comparison.
  {
    std::string baselineName(yamlRoot(name + "_mixed"));
    const auto paths = baselinePaths<ExecSpace>();
#if defined(AXOM_TESTING_GENERATE_BASELINES)
    saveBaseline(paths, baselineName, hostClipMixedMesh);
#else
    EXPECT_TRUE(compareBaseline(paths, baselineName, hostClipMixedMesh));
#endif
  }
}

TEST(mir_clipfield, uniform2d)
{
  braid2d_clip_test<seq_exec>("uniform", "uniform2d");

#if defined(AXOM_USE_OPENMP)
  braid2d_clip_test<omp_exec>("uniform", "uniform2d");
#endif

#if defined(AXOM_USE_CUDA)
  braid2d_clip_test<cuda_exec>("uniform", "uniform2d");
#endif

#if defined(AXOM_USE_HIP)
  braid2d_clip_test<hip_exec>("uniform", "uniform2d");
#endif
}

//------------------------------------------------------------------------------
template <typename ExecSpace, int NDIMS>
void braid_rectilinear_clip_test(const std::string &name)
{
  using Indexing = axom::mir::views::StructuredIndexing<axom::IndexType, NDIMS>;
  using TopoView = axom::mir::views::StructuredTopologyView<Indexing>;

  axom::StackArray<axom::IndexType, NDIMS> dims, zoneDims;
  for(int i = 0; i < NDIMS; i++)
  {
    dims[i] = 10;
    zoneDims[i] = dims[i] - 1;
  }

  // Create the data
  conduit::Node hostMesh, deviceMesh;
  axom::mir::testing::data::braid("rectilinear", dims, hostMesh);
#if defined(AXOM_TESTING_SAVE_VISUALIZATION) && defined(AXOM_USE_HDF5)
  conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig", "hdf5");
#endif

  // host->device
  bputils::copy<ExecSpace>(deviceMesh, hostMesh);

  // Create views
  auto coordsetView = axom::mir::views::make_rectilinear_coordset<double, NDIMS>::view(
    deviceMesh["coordsets/coords"]);
  using CoordsetView = decltype(coordsetView);
  TopoView topoView(Indexing {zoneDims});

  // Create options to control the clipping.
  conduit::Node options;
  options["clipField"] = "distance";
  options["inside"] = 1;
  options["outside"] = 1;

  // Clip the data
  conduit::Node deviceClipMesh;
  axom::mir::clipping::ClipField<ExecSpace, TopoView, CoordsetView> clipper(topoView, coordsetView);
  clipper.execute(deviceMesh, options, deviceClipMesh);

  // Copy device->host
  conduit::Node hostClipMesh;
  bputils::copy<seq_exec>(hostClipMesh, deviceClipMesh);

  // Handle baseline comparison.
  {
    std::string baselineName(yamlRoot(name));
    const auto paths = baselinePaths<ExecSpace>();
#if defined(AXOM_TESTING_GENERATE_BASELINES)
    saveBaseline(paths, baselineName, hostClipMesh);
#else
    EXPECT_TRUE(compareBaseline(paths, baselineName, hostClipMesh));
#endif
  }
}

TEST(mir_clipfield, rectilinear2d)
{
  braid_rectilinear_clip_test<seq_exec, 2>("rectilinear2d");

#if defined(AXOM_USE_OPENMP)
  braid_rectilinear_clip_test<omp_exec, 2>("rectilinear2d");
#endif

#if defined(AXOM_USE_CUDA)
  braid_rectilinear_clip_test<cuda_exec, 2>("rectilinear2d");
#endif

#if defined(AXOM_USE_HIP)
  braid_rectilinear_clip_test<hip_exec, 2>("rectilinear2d");
#endif
}

TEST(mir_clipfield, rectilinear3d)
{
  braid_rectilinear_clip_test<seq_exec, 3>("rectilinear3d");

#if defined(AXOM_USE_OPENMP)
  braid_rectilinear_clip_test<omp_exec, 3>("rectilinear3d");
#endif

#if defined(AXOM_USE_CUDA)
  braid_rectilinear_clip_test<cuda_exec, 3>("rectilinear3d");
#endif

#if defined(AXOM_USE_HIP)
  braid_rectilinear_clip_test<hip_exec, 3>("rectilinear3d");
#endif
}

//------------------------------------------------------------------------------
template <typename ExecSpace, int NDIMS>
void strided_structured_clip_test(const std::string &name, const conduit::Node &options)
{
  // Create the data
  conduit::Node hostMesh, deviceMesh;
  axom::mir::testing::data::strided_structured<NDIMS>(hostMesh);
  //hostMesh.print();
#if defined(AXOM_TESTING_SAVE_VISUALIZATION) && defined(AXOM_USE_HDF5)
  conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig", "hdf5");
  conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig_yaml", "yaml");
#endif

  conduit::Node deviceOptions, deviceClipMesh, hostClipMesh;

  // host->device
  bputils::copy<ExecSpace>(deviceMesh, hostMesh);
  bputils::copy<ExecSpace>(deviceOptions, options);

  // Create views
  auto coordsetView =
    axom::mir::views::make_explicit_coordset<double, 2>::view(deviceMesh["coordsets/coords"]);
  auto topoView = axom::mir::views::make_strided_structured<2>::view(deviceMesh["topologies/mesh"]);

  using CoordsetView = decltype(coordsetView);
  using TopoView = decltype(topoView);

  // Clip the data
  axom::mir::clipping::ClipField<ExecSpace, TopoView, CoordsetView> clipper(topoView, coordsetView);
  clipper.execute(deviceMesh, deviceOptions, deviceClipMesh);

  // device->host
  bputils::copy<seq_exec>(hostClipMesh, deviceClipMesh);

  // Handle baseline comparison.
  {
    std::string baselineName(yamlRoot(name));
    const auto paths = baselinePaths<ExecSpace>();
#if defined(AXOM_TESTING_GENERATE_BASELINES)
    saveBaseline(paths, baselineName, hostClipMesh);
#else
    EXPECT_TRUE(compareBaseline(paths, baselineName, hostClipMesh));
#endif
  }
}

void strided_structured_clip_test_exec(const std::string &name, const conduit::Node &options)
{
  strided_structured_clip_test<seq_exec, 2>(name, options);

#if defined(AXOM_USE_OPENMP)
  strided_structured_clip_test<omp_exec, 2>(name, options);
#endif

#if defined(AXOM_USE_CUDA)
  strided_structured_clip_test<cuda_exec, 2>(name, options);
#endif

#if defined(AXOM_USE_HIP)
  strided_structured_clip_test<hip_exec, 2>(name, options);
#endif
}

TEST(mir_clipfield, strided_structured_2d)
{
  // Create options to control the clipping.
  conduit::Node options;
  options["clipField"] = "vert_vals";
  options["clipValue"] = 6.5;
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
  using TopoView = axom::mir::views::UnstructuredTopologySingleShapeView<ShapeType>;
  using CoordsetView = axom::mir::views::ExplicitCoordsetView<double, 3>;

  // Create the data
  const axom::StackArray<axom::IndexType, 3> dims {10, 10, 10};
  conduit::Node hostMesh, deviceMesh;
  axom::mir::testing::data::braid(type, dims, hostMesh);
  bputils::copy<ExecSpace>(deviceMesh, hostMesh);
#if defined(AXOM_TESTING_SAVE_VISUALIZATION) && defined(AXOM_USE_HDF5)
  conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig", "hdf5");
#endif

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
  options["clipField"] = "distance";
  options["inside"] = 1;
  options["outside"] = 0;

  // Clip the data
  conduit::Node deviceClipMesh;
  axom::mir::clipping::ClipField<ExecSpace, TopoView, CoordsetView> clipper(topoView, coordsetView);
  clipper.execute(deviceMesh, options, deviceClipMesh);

  // Copy device->host
  conduit::Node hostClipMesh;
  bputils::copy<seq_exec>(hostClipMesh, deviceClipMesh);

  // Handle baseline comparison.
  std::string baselineName(yamlRoot(name));
  const auto paths = baselinePaths<ExecSpace>();
#if defined(AXOM_TESTING_GENERATE_BASELINES)
  saveBaseline(paths, baselineName, hostClipMesh);
#else
  EXPECT_TRUE(compareBaseline(paths, baselineName, hostClipMesh));
#endif
}

/// Execute the braid3d test for a single shape on multiple ExecSpaces
template <typename ShapeType>
void braid3d_clip_test_exec(const std::string &type, const std::string &name)
{
  braid3d_clip_test<seq_exec, ShapeType>(type, name);

#if defined(AXOM_USE_OPENMP)
  braid3d_clip_test<omp_exec, ShapeType>(type, name);
#endif

#if defined(AXOM_USE_CUDA)
  braid3d_clip_test<cuda_exec, ShapeType>(type, name);
#endif

#if defined(AXOM_USE_HIP)
  braid3d_clip_test<hip_exec, ShapeType>(type, name);
#endif
}

TEST(mir_clipfield, tet) { braid3d_clip_test_exec<axom::mir::views::TetShape<int>>("tets", "tet"); }

TEST(mir_clipfield, pyramid)
{
  braid3d_clip_test_exec<axom::mir::views::PyramidShape<int>>("pyramids", "pyr");
}

TEST(mir_clipfield, wedge)
{
  braid3d_clip_test_exec<axom::mir::views::WedgeShape<int>>("wedges", "wdg");
}

TEST(mir_clipfield, hex) { braid3d_clip_test_exec<axom::mir::views::HexShape<int>>("hexs", "hex"); }

//------------------------------------------------------------------------------
template <typename ExecSpace>
void braid3d_mixed_clip_test(const std::string &name)
{
  using CoordType = float;
  using ConnType = int;
  using TopoView = axom::mir::views::UnstructuredTopologyMixedShapeView<ConnType>;
  using CoordsetView = axom::mir::views::ExplicitCoordsetView<CoordType, 3>;

  // Create the data
  conduit::Node hostMesh, deviceMesh;
  axom::mir::testing::data::mixed3d(hostMesh);
  bputils::copy<ExecSpace>(deviceMesh, hostMesh);
#if defined(AXOM_TESTING_SAVE_VISUALIZATION) && defined(AXOM_USE_HDF5)
  conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig", "hdf5");
#endif

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
  auto shapeMap = axom::mir::views::buildShapeMap(n_device_topo,
                                                  values,
                                                  ids,
                                                  axom::execution_space<ExecSpace>::allocatorID());

  TopoView topoView(connView, shapesView, sizesView, offsetsView, shapeMap);

  // Create options to control the clipping.
  conduit::Node options;
  options["clipField"] = "distance";
  options["clipValue"] = 12.f;
  options["inside"] = 1;
  options["outside"] = 0;

  // Clip the data
  conduit::Node deviceClipMesh;
  axom::mir::clipping::ClipField<ExecSpace, TopoView, CoordsetView> clipper(topoView, coordsetView);
  clipper.execute(deviceMesh, options, deviceClipMesh);

  // Copy device->host
  conduit::Node hostClipMesh;
  bputils::copy<seq_exec>(hostClipMesh, deviceClipMesh);

  // Handle baseline comparison.
  std::string baselineName(yamlRoot(name));
  const auto paths = baselinePaths<ExecSpace>();
#if defined(AXOM_TESTING_GENERATE_BASELINES)
  saveBaseline(paths, baselineName, hostClipMesh);
#else
  EXPECT_TRUE(compareBaseline(paths, baselineName, hostClipMesh));
#endif
}

TEST(mir_clipfield, mixed_seq) { braid3d_mixed_clip_test<seq_exec>("mixed"); }
#if defined(AXOM_USE_OPENMP)
TEST(mir_clipfield, mixed_omp) { braid3d_mixed_clip_test<omp_exec>("mixed"); }
#endif
#if defined(AXOM_USE_CUDA)
TEST(mir_clipfield, mixed_cuda) { braid3d_mixed_clip_test<cuda_exec>("mixed"); }
#endif
#if defined(AXOM_USE_HIP)
TEST(mir_clipfield, mixed_hip) { braid3d_mixed_clip_test<hip_exec>("mixed"); }
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

    // host->device
    conduit::Node deviceMesh;
    bputils::copy<ExecSpace>(deviceMesh, hostMesh);

    // Set up views for the mesh.
    using CoordsetView = axom::mir::views::ExplicitCoordsetView<float, 2>;
    CoordsetView coordsetView(
      bputils::make_array_view<float>(deviceMesh.fetch_existing("coordsets/coords/values/x")),
      bputils::make_array_view<float>(deviceMesh.fetch_existing("coordsets/coords/values/y")));

    using TopologyView =
      axom::mir::views::UnstructuredTopologySingleShapeView<axom::mir::views::QuadShape<int>>;
    TopologyView topologyView(bputils::make_array_view<int>(
      deviceMesh.fetch_existing("topologies/mesh/elements/connectivity")));

    // Clip
    conduit::Node options, deviceClipMesh;
    options["clipField"] = "clip";
    options["clipValue"] = 0.5;
    using Clip = axom::mir::clipping::ClipField<ExecSpace, TopologyView, CoordsetView>;
    Clip clip(topologyView, coordsetView);
    clip.execute(deviceMesh, options, deviceClipMesh);

    // device->host
    conduit::Node hostClipMesh;
    bputils::copy<axom::SEQ_EXEC>(hostClipMesh, deviceClipMesh);
    //printNode(hostClipMesh);

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

    // Check that the degenerate quads were turned into triangles.
    std::vector<int> shapes {{2, 2, 3, 2}};
    std::vector<int> sizes {{3, 3, 4, 3}};
    std::vector<int> offsets {{0, 4, 8, 12}};
    compare_values(shapes, hostClipMesh["topologies/mesh/elements/shapes"].as_int_accessor());
    compare_values(sizes, hostClipMesh["topologies/mesh/elements/sizes"].as_int_accessor());
    compare_values(offsets, hostClipMesh["topologies/mesh/elements/offsets"].as_int_accessor());
  }
};

TEST(mir_clipfield, pointmerging_seq) { point_merge_test<seq_exec>::test(); }
#if defined(AXOM_USE_OPENMP)
TEST(mir_clipfield, pointmerging_omp) { point_merge_test<omp_exec>::test(); }
#endif
#if defined(AXOM_USE_CUDA)
TEST(mir_clipfield, pointmerging_cuda) { point_merge_test<cuda_exec>::test(); }
#endif
#if defined(AXOM_USE_HIP)
TEST(mir_clipfield, pointmerging_hip) { point_merge_test<hip_exec>::test(); }
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
    bputils::copy<ExecSpace>(deviceMesh, hostMesh);

    // Wrap the data in views.
    auto coordsetView = axom::mir::views::make_rectilinear_coordset<conduit::float64, 2>::view(
      deviceMesh["coordsets/coords"]);
    using CoordsetView = decltype(coordsetView);

    auto topologyView = axom::mir::views::make_rectilinear<2>::view(deviceMesh["topologies/mesh"]);
    using TopologyView = decltype(topologyView);

    conduit::Node hostOptions;
    hostOptions["selectedZones"].set(std::vector<axom::IndexType> {{1, 3, 4, 5, 7}});
    hostOptions["inside"] = 1;
    hostOptions["outside"] = 1;
    hostOptions["clipField"] = "zero";

    conduit::Node deviceOptions, deviceResult;
    bputils::copy<ExecSpace>(deviceOptions, hostOptions);

    axom::mir::clipping::ClipField<ExecSpace, TopologyView, CoordsetView> clip(topologyView,
                                                                               coordsetView);
    clip.execute(deviceMesh, deviceOptions, deviceResult);

    // device->host
    conduit::Node hostResult;
    bputils::copy<seq_exec>(hostResult, deviceResult);
    //printNode(hostResult);

    // Handle baseline comparison.
    const auto paths = baselinePaths<ExecSpace>();
    std::string baselineName(yamlRoot("selectedzones1"));
#if defined(AXOM_TESTING_GENERATE_BASELINES)
    saveBaseline(paths, baselineName, hostResult);
#else
    EXPECT_TRUE(compareBaseline(paths, baselineName, hostResult));
#endif

    //---------------------
    // Try a different clip
    hostOptions["outside"] = 0;
    hostOptions["clipField"] = "radial";
    hostOptions["clipValue"] = 3.2;
    bputils::copy<ExecSpace>(deviceOptions, hostOptions);
    deviceResult.reset();
    clip.execute(deviceMesh, deviceOptions, deviceResult);

    // device->host
    bputils::copy<seq_exec>(hostResult, deviceResult);
    //printNode(hostResult);

    baselineName = yamlRoot("selectedzones2");
#if defined(AXOM_TESTING_GENERATE_BASELINES)
    saveBaseline(paths, baselineName, hostResult);
#else
    EXPECT_TRUE(compareBaseline(paths, baselineName, hostResult));
#endif
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

TEST(mir_clipfield, selectedzones_seq) { test_selectedzones<seq_exec>::test(); }
#if defined(AXOM_USE_OPENMP)
TEST(mir_clipfield, selectedzones_omp) { test_selectedzones<omp_exec>::test(); }
#endif
#if defined(AXOM_USE_CUDA)
TEST(mir_clipfield, selectedzones_cuda) { test_selectedzones<cuda_exec>::test(); }
#endif
#if defined(AXOM_USE_HIP)
TEST(mir_clipfield, selectedzones_hip) { test_selectedzones<hip_exec>::test(); }
#endif

//------------------------------------------------------------------------------
void conduit_debug_err_handler(const std::string &s1, const std::string &s2, int i1)
{
  std::cout << "s1=" << s1 << ", s2=" << s2 << ", i1=" << i1 << std::endl;
  // This is on purpose.
  while(1)
    ;
}

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::CLI::App app;
#if defined(AXOM_USE_CALIPER)
  std::string annotationMode("none");
  app.add_option("--caliper", annotationMode)
    ->description(
      "caliper annotation mode. Valid options include 'none' and 'report'. "
      "Use 'help' to see full list.")
    ->capture_default_str()
    ->check(axom::utilities::ValidCaliperMode);
#endif
  bool handlerEnabled = false;
  app.add_flag("--handler", handlerEnabled, "Enable Conduit handler.");

  // Parse command line options.
  try
  {
    app.parse(argc, argv);

#if defined(AXOM_USE_CALIPER)
    axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(annotationMode);
#endif

    axom::slic::SimpleLogger logger;  // create & initialize test logger,
    if(handlerEnabled)
    {
      conduit::utils::set_error_handler(conduit_debug_err_handler);
    }

    result = RUN_ALL_TESTS();
  }
  catch(axom::CLI::CallForHelp &e)
  {
    std::cout << app.help() << std::endl;
    result = 0;
  }
  catch(axom::CLI::ParseError &e)
  {
    // Handle other parsing errors
    std::cerr << e.what() << std::endl;
    result = app.exit(e);
  }

  return result;
}
