// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/mir.hpp"
#include "axom/primal.hpp"
#include "axom/mir/tests/mir_testing_data_helpers.hpp"

#include <cmath>

namespace bputils = axom::mir::utilities::blueprint;

//------------------------------------------------------------------------------

//#define DEBUGGING_TEST_CASES

// Uncomment to generate baselines
//#define AXOM_TESTING_GENERATE_BASELINES

// Uncomment to save visualization files for debugging (when making baselines)
//#define AXOM_TESTING_SAVE_VISUALIZATION

#include "axom/mir/tests/mir_testing_helpers.hpp"

std::string baselineDirectory() { return pjoin(dataDirectory(), "mir", "regression", "mir_views"); }
//------------------------------------------------------------------------------

TEST(mir_views, shape2conduitName)
{
  EXPECT_STREQ(axom::mir::views::LineShape<int>::name(), "line");
  EXPECT_STREQ(axom::mir::views::LineShape<long>::name(), "line");

  EXPECT_STREQ(axom::mir::views::TriShape<int>::name(), "tri");
  EXPECT_STREQ(axom::mir::views::TriShape<long>::name(), "tri");

  EXPECT_STREQ(axom::mir::views::QuadShape<int>::name(), "quad");
  EXPECT_STREQ(axom::mir::views::QuadShape<long>::name(), "quad");

  EXPECT_STREQ(axom::mir::views::TetShape<int>::name(), "tet");
  EXPECT_STREQ(axom::mir::views::TetShape<long>::name(), "tet");

  EXPECT_STREQ(axom::mir::views::PyramidShape<int>::name(), "pyramid");
  EXPECT_STREQ(axom::mir::views::PyramidShape<long>::name(), "pyramid");

  EXPECT_STREQ(axom::mir::views::WedgeShape<int>::name(), "wedge");
  EXPECT_STREQ(axom::mir::views::WedgeShape<long>::name(), "wedge");

  EXPECT_STREQ(axom::mir::views::HexShape<int>::name(), "hex");
  EXPECT_STREQ(axom::mir::views::HexShape<long>::name(), "hex");
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_node_to_arrayview
{
  static int constexpr sum(int n)
  {
    int s = 0;
    for(int i = 0; i < n; i++) s += i;
    return s;
  }

  static void test()
  {
    std::vector<int> dtypes {conduit::DataType::INT8_ID,
                             conduit::DataType::INT16_ID,
                             conduit::DataType::INT32_ID,
                             conduit::DataType::INT64_ID,
                             conduit::DataType::UINT8_ID,
                             conduit::DataType::UINT16_ID,
                             conduit::DataType::UINT32_ID,
                             conduit::DataType::UINT64_ID,
                             conduit::DataType::FLOAT32_ID,
                             conduit::DataType::FLOAT64_ID};
    constexpr int n = 16;
    axom::mir::utilities::blueprint::ConduitAllocateThroughAxom<ExecSpace> c2a;
    for(int dtype : dtypes)
    {
      // Make a node and fill it with data.
      conduit::Node n_data;
      n_data.set_allocator(c2a.getConduitAllocatorID());
      n_data.set(conduit::DataType(dtype, n));

      int sumValues = 0;
      axom::mir::views::Node_to_ArrayView(n_data,
                                          [&](auto dataView) { sumValues = testBody(dataView, n); });

      EXPECT_EQ(sumValues, sum(n));
    }
  }

  template <typename DataView>
  static int testBody(DataView dataView, int n)
  {
    using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;
    using value_type = typename DataView::value_type;

    std::cout << axom::mir::views::array_view_traits<DataView>::name() << std::endl;

    // Make sure we can store values in dataView
    axom::for_all<ExecSpace>(
      n,
      AXOM_LAMBDA(axom::IndexType index) { dataView[index] = static_cast<value_type>(index); });

    // Read the values and sum them.
    RAJA::ReduceSum<reduce_policy, value_type> sumValues_reduce(0);
    axom::for_all<ExecSpace>(
      n,
      AXOM_LAMBDA(axom::IndexType index) { sumValues_reduce += dataView[index]; });
    return static_cast<int>(sumValues_reduce.get());
  }
};

TEST(mir_views, node_to_arrayview_seq) { test_node_to_arrayview<seq_exec>::test(); }
#if defined(AXOM_USE_OPENMP)
TEST(mir_views, node_to_arrayview_omp) { test_node_to_arrayview<omp_exec>::test(); }
#endif
#if defined(AXOM_USE_CUDA)
TEST(mir_views, node_to_arrayview_cuda) { test_node_to_arrayview<cuda_exec>::test(); }
#endif
#if defined(AXOM_USE_HIP)
TEST(mir_views, node_to_arrayview_hip) { test_node_to_arrayview<hip_exec>::test(); }
#endif

//------------------------------------------------------------------------------
TEST(mir_views, explicit_coordsetview)
{
  axom::Array<float> x {{0., 1., 2., 3., 4., 5.}};
  axom::Array<float> y {{10., 11., 12., 13., 14., 15.}};
  axom::Array<float> z {{20., 21., 22., 23., 24., 25.}};

  axom::mir::views::ExplicitCoordsetView<float, 2> view2d(x.view(), y.view());
  EXPECT_EQ(view2d.size(), 6);
  for(axom::IndexType i = 0; i < view2d.size(); i++)
  {
    axom::primal::Point<float, 2> P({x[i], y[i]});
    EXPECT_EQ(view2d.getPoint(i), P);
    EXPECT_EQ(view2d[i], P);
  }

  axom::mir::views::ExplicitCoordsetView<float, 3> view3d(x.view(), y.view(), z.view());
  EXPECT_EQ(view3d.size(), 6);
  for(axom::IndexType i = 0; i < view3d.size(); i++)
  {
    axom::primal::Point<float, 3> P({x[i], y[i], z[i]});
    EXPECT_EQ(view3d.getPoint(i), P);
    EXPECT_EQ(view3d[i], P);
  }
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_structured_topology_view_rectilinear
{
  static void test()
  {
    conduit::Node hostMesh;
    create(hostMesh);

    // host->device
    conduit::Node deviceMesh;
    bputils::copy<ExecSpace>(deviceMesh, hostMesh);

    // Make results view on device.
    constexpr int nzones = 9;
    axom::Array<axom::IndexType> results(nzones,
                                         nzones,
                                         axom::execution_space<ExecSpace>::allocatorID());
    auto resultsView = results.view();

    // Execute the kernel for each zone (find max node number in zone).
    auto topoView = axom::mir::views::make_rectilinear<2>::view(deviceMesh["topologies/mesh"]);
    axom::for_all<ExecSpace>(
      topoView.numberOfZones(),
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
        const auto zone = topoView.zone(zoneIndex);
        axom::IndexType m = -1;
        for(const auto &id : zone.getIds())
        {
          m = axom::utilities::max(static_cast<axom::IndexType>(id), m);
        }
        resultsView[zoneIndex] = m;
      });

    // device->host
    axom::Array<axom::IndexType> hostResults(nzones,
                                             nzones,
                                             axom::execution_space<axom::SEQ_EXEC>::allocatorID());
    axom::copy(hostResults.data(), results.data(), nzones * sizeof(axom::IndexType));

    // Compare.
    const axom::IndexType expected[] = {5, 6, 7, 9, 10, 11, 13, 14, 15};
    for(int i = 0; i < nzones; i++)
    {
      EXPECT_EQ(hostResults[i], expected[i]);
    }
  }

  static void create(conduit::Node &mesh)
  {
    std::vector<int> dims {4, 4};
    axom::mir::testing::data::braid("rectilinear", dims, mesh);
  }
};

TEST(mir_views, stopo_rectilinear_2d_seq)
{
  test_structured_topology_view_rectilinear<seq_exec>::test();
}
#if defined(AXOM_USE_OPENMP)
TEST(mir_views, stopo_rectilinear_2d_omp)
{
  test_structured_topology_view_rectilinear<omp_exec>::test();
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(mir_views, stopo_rectilinear_2d_cuda)
{
  test_structured_topology_view_rectilinear<cuda_exec>::test();
}
#endif
#if defined(AXOM_USE_HIP)
TEST(mir_views, stopo_rectilinear_2d_hip)
{
  test_structured_topology_view_rectilinear<hip_exec>::test();
}
#endif

//------------------------------------------------------------------------------
struct test_strided_structured
{
  static void test()
  {
    conduit::Node hostMesh;
    axom::mir::testing::data::strided_structured<2>(hostMesh);
    //  hostMesh.print();

    axom::mir::views::dispatch_explicit_coordset(hostMesh["coordsets/coords"], [&](auto coordsetView) {
      axom::mir::views::dispatch_structured_topology<axom::mir::views::select_dimensions(2)>(
        hostMesh["topologies/mesh"],
        [&](const std::string &AXOM_UNUSED_PARAM(shape), auto topoView) {
          execute(coordsetView, topoView);
        });
    });
  }

  template <typename CoordsetView, typename TopologyView>
  static void execute(CoordsetView coordsetView, TopologyView topoView)
  {
    using ExecSpace = seq_exec;

    // These are the expected node ids for this strided structured mesh.
    // clang-format off
    const axom::Array<int> expectedNodes {{16, 17, 24, 23,
                                           17, 18, 25, 24,
                                           18, 19, 26, 25,
                                           23, 24, 31, 30,
                                           24, 25, 32, 31,
                                           25, 26, 33, 32}};
    // clang-format on
    auto expectedNodesView = expectedNodes.view();
    axom::IndexType n4 = expectedNodesView.size();

    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    axom::Array<int> actualNodes(n4, n4, allocatorID);
    axom::Array<int> logicalNodes(n4 * 2, n4 * 2, allocatorID);
    auto actualNodesView = actualNodes.view();
    auto logicalNodesView = logicalNodes.view();

    // Traverse the zones in the mesh and gather node ids
    axom::for_all<ExecSpace>(
      topoView.numberOfZones(),
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
        const auto zone = topoView.zone(zoneIndex);
        const auto nodeIndexing = topoView.indexing().expand();

        // Get node ids for zone.
        const auto ids = zone.getIds();
        for(axom::IndexType i = 0; i < ids.size(); i++)
        {
          actualNodesView[zoneIndex * 4 + i] = ids[i];

          // Get the logical local id for the id.
          const auto index = nodeIndexing.GlobalToLocal(ids[i]);
          const auto logical = nodeIndexing.IndexToLogicalIndex(index);
          logicalNodesView[(zoneIndex * 4 + i) * 2 + 0] = logical[0];
          logicalNodesView[(zoneIndex * 4 + i) * 2 + 1] = logical[1];
        }
      });

    for(axom::IndexType i = 0; i < n4; i++)
    {
      EXPECT_EQ(expectedNodesView[i], actualNodesView[i]);
    }

    // Check coordinates
    for(axom::IndexType i = 0; i < n4; i++)
    {
      const auto id = actualNodesView[i];

      // Get coordinate from coordsetView.
      const auto pt = coordsetView[id];

      // Get the logical local id for the id.
      const auto logicalI = logicalNodesView[i * 2 + 0];
      const auto logicalJ = logicalNodesView[i * 2 + 1];

      // Expected coordinate
      double x = (3. + 1. / 3.) * static_cast<double>(logicalI - 1);
      const double yvals[] = {-2, 2, 6};
      double y = yvals[logicalJ];

      const double dx = pt[0] - x;
      const double dy = pt[1] - y;
      double d = sqrt(dx * dx + dy * dy);

      EXPECT_TRUE(d < 1.e-10);
    }
  }
};

TEST(mir_views, strided_structured_seq) { test_strided_structured::test(); }

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_braid2d_mat
{
  static void test(const std::string &type,
                   const std::string &mattype,
#if defined(AXOM_TESTING_SAVE_VISUALIZATION)
                   const std::string &name
#else
                   const std::string &AXOM_UNUSED_PARAM(name)
#endif
  )
  {
    namespace bputils = axom::mir::utilities::blueprint;
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    axom::StackArray<axom::IndexType, 2> dims {10, 10};
    axom::StackArray<axom::IndexType, 2> zoneDims {dims[0] - 1, dims[1] - 1};
    const axom::IndexType nzones = zoneDims[0] * zoneDims[1];

    // Create the data
    conduit::Node hostMesh, deviceMesh;
    axom::mir::testing::data::braid(type, dims, hostMesh);
    axom::mir::testing::data::make_matset(mattype, "mesh", zoneDims, hostMesh);
    axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);
#if defined(AXOM_TESTING_SAVE_VISUALIZATION) && defined(AXOM_USE_HDF5)
    conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig", "hdf5");
#endif

    if(mattype == "unibuffer")
    {
      // clang-format off
      // _mir_views_matsetview_begin
      using MatsetView = axom::mir::views::UnibufferMaterialView<int, float, 3>;
      MatsetView matsetView;
      matsetView.set(bputils::make_array_view<int>(deviceMesh["matsets/mat/material_ids"]),
                     bputils::make_array_view<float>(deviceMesh["matsets/mat/volume_fractions"]),
                     bputils::make_array_view<int>(deviceMesh["matsets/mat/sizes"]),
                     bputils::make_array_view<int>(deviceMesh["matsets/mat/offsets"]),
                     bputils::make_array_view<int>(deviceMesh["matsets/mat/indices"]));
      // _mir_views_matsetview_end
      // clang-format on
      test_matsetview(nzones, matsetView, allocatorID);
    }
    else if(mattype == "multibuffer")
    {
      axom::mir::views::dispatch_material_multibuffer(
        deviceMesh["matsets/mat"],
        [&](auto matsetView) { test_matsetview(nzones, matsetView, allocatorID); });
    }
    else if(mattype == "element_dominant")
    {
      axom::mir::views::dispatch_material_element_dominant(
        deviceMesh["matsets/mat"],
        [&](auto matsetView) { test_matsetview(nzones, matsetView, allocatorID); });
    }
    else if(mattype == "material_dominant")
    {
      axom::mir::views::dispatch_material_material_dominant(
        deviceMesh["matsets/mat"],
        [&](auto matsetView) { test_matsetview(nzones, matsetView, allocatorID); });
    }
  }

  template <typename MatsetView>
  static void test_matsetview(axom::IndexType nzones, MatsetView matsetView, int allocatorID)
  {
    // These values are used in that material_map.
    constexpr int MATA = 22;
    constexpr int MATB = 66;
    constexpr int MATC = 33;
    // The zone ids that are being queried.
    const int zoneids[] = {0, 36, 40};

    // clang-format off
    int results[] = {/*nzones*/ static_cast<int>(nzones),
                     /*contains mat*/ 0, 1, 0, /*nmats in zone*/ 1, /*ids.size*/ 1, /*mats in zone*/ MATB, -1, -1,
                     /*contains mat*/ 1, 1, 0, /*nmats in zone*/ 2, /*ids.size*/ 2, /*mats in zone*/ MATA, MATB, -1,
                     /*contains mat*/ 1, 1, 1, /*nmats in zone*/ 3, /*ids.size*/ 3, /*mats in zone*/ MATA, MATB, MATC};
    // clang-format on
    constexpr int nTestZones = sizeof(zoneids) / sizeof(int);

    // Get zoneids into zoneidsView for device.
    axom::Array<int> zoneidsArray(nTestZones, nTestZones, allocatorID);
    axom::copy(zoneidsArray.data(), zoneids, sizeof(int) * nTestZones);
    auto zoneidsView = zoneidsArray.view();

    // Allocate results array on device.
    constexpr int nResults = sizeof(results) / sizeof(int);
    axom::Array<int> resultsArrayDevice(nResults, nResults, allocatorID);
    auto resultsView = resultsArrayDevice.view();

    // Fill in resultsView on the device.
    constexpr int nResultsPerZone = 8;
    axom::for_all<ExecSpace>(
      nTestZones,
      AXOM_LAMBDA(axom::IndexType index) {
        if(index == 0)
        {
          // Compute number of zones here since some views need to look inside
          // data to determine the number of zones.
          resultsView[0] = matsetView.numberOfZones();
        }
        const int offset = 1 + nResultsPerZone * index;
        // contains mat
        resultsView[offset + 0] = matsetView.zoneContainsMaterial(zoneidsView[index], MATA) ? 1 : 0;
        resultsView[offset + 1] = matsetView.zoneContainsMaterial(zoneidsView[index], MATB) ? 1 : 0;
        resultsView[offset + 2] = matsetView.zoneContainsMaterial(zoneidsView[index], MATC) ? 1 : 0;
        // nmats in zone
        resultsView[offset + 3] = matsetView.numberOfMaterials(zoneidsView[index]);

        typename MatsetView::IDList ids {};
        typename MatsetView::VFList vfs {};
        // ids.size
        matsetView.zoneMaterials(zoneidsView[index], ids, vfs);
        resultsView[offset + 4] = ids.size();
        // mats in zone
        for(axom::IndexType i = 0; i < 3; i++)
        {
          resultsView[offset + 5 + i] = (i < ids.size()) ? static_cast<int>(ids[i]) : -1;
        }
      });
    // Get containsView data to the host and compare results
    std::vector<int> resultsHost(nResults);
    axom::copy(resultsHost.data(), resultsView.data(), sizeof(int) * nResults);
    for(int i = 0; i < nResults; i++)
    {
      EXPECT_EQ(results[i], resultsHost[i]);
    }
  }
};

// Unibuffer
TEST(mir_views, matset_unibuffer_seq)
{
  test_braid2d_mat<seq_exec>::test("uniform", "unibuffer", "uniform2d_unibuffer");
}
#if defined(AXOM_USE_OPENMP)
TEST(mir_views, matset_unibuffer_omp)
{
  test_braid2d_mat<omp_exec>::test("uniform", "unibuffer", "uniform2d_unibuffer");
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(mir_views, matset_unibuffer_cuda)
{
  test_braid2d_mat<cuda_exec>::test("uniform", "unibuffer", "uniform2d_unibuffer");
}
#endif
#if defined(AXOM_USE_HIP)
TEST(mir_views, matset_unibuffer_hip)
{
  test_braid2d_mat<hip_exec>::test("uniform", "unibuffer", "uniform2d_unibuffer");
}
#endif

// Multibuffer
TEST(mir_views, matset_multibuffer_seq)
{
  test_braid2d_mat<seq_exec>::test("uniform", "multibuffer", "uniform2d_multibuffer");
}
#if defined(AXOM_USE_OPENMP)
TEST(mir_views, matset_multibuffer_omp)
{
  test_braid2d_mat<omp_exec>::test("uniform", "multibuffer", "uniform2d_multibuffer");
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(mir_views, matset_multibuffer_cuda)
{
  test_braid2d_mat<cuda_exec>::test("uniform", "multibuffer", "uniform2d_multibuffer");
}
#endif
#if defined(AXOM_USE_HIP)
TEST(mir_views, matset_multibuffer_hip)
{
  test_braid2d_mat<hip_exec>::test("uniform", "multibuffer", "uniform2d_multibuffer");
}
#endif

// Element-dominant
TEST(mir_views, matset_element_dominant_seq)
{
  test_braid2d_mat<seq_exec>::test("uniform", "element_dominant", "uniform2d_element_dominant");
}
#if defined(AXOM_USE_OPENMP)
TEST(mir_views, matset_element_dominant_omp)
{
  test_braid2d_mat<omp_exec>::test("uniform", "element_dominant", "uniform2d_element_dominant");
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(mir_views, matset_element_dominant_cuda)
{
  test_braid2d_mat<cuda_exec>::test("uniform", "element_dominant", "uniform2d_element_dominant");
}
#endif
#if defined(AXOM_USE_HIP)
TEST(mir_views, matset_element_dominant_hip)
{
  test_braid2d_mat<hip_exec>::test("uniform", "element_dominant", "uniform2d_element_dominant");
}
#endif

// Material-dominant
TEST(mir_views, matset_material_dominant_seq)
{
  test_braid2d_mat<seq_exec>::test("uniform", "material_dominant", "uniform2d_material_dominant");
}
#if defined(AXOM_USE_OPENMP)
TEST(mir_views, matset_material_dominant_omp)
{
  test_braid2d_mat<omp_exec>::test("uniform", "material_dominant", "uniform2d_material_dominant");
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(mir_views, matset_material_dominant_cuda)
{
  test_braid2d_mat<cuda_exec>::test("uniform", "material_dominant", "uniform2d_material_dominant");
}
#endif
#if defined(AXOM_USE_HIP)
TEST(mir_views, matset_material_dominant_hip)
{
  test_braid2d_mat<hip_exec>::test("uniform", "material_dominant", "uniform2d_material_dominant");
}
#endif

TEST(mir_views, matset_multibuffer)
{
  const char *yaml = R"(
matsets:
  matset:
    topology: topology
    volume_fractions:
      a:
        values: [0, 0, 0, 0.33, 0, 0.3]  # [0, 0, 0, a1, 0, a0]
        indices: [5, 3]
      b:
        values: [0, 0.7, 1., 0.67, 0]    # [0, b0, b2, b1, 0]
        indices: [1, 3, 2]
    material_map: # (optional)
      a: 0
      b: 1
)";
  conduit::Node matsets;
  matsets.parse(yaml);
  const conduit::Node &n_matset = matsets["matsets/matset"];
  axom::mir::views::dispatch_material_multibuffer(n_matset, [&](auto matsetView) {
    using IDList = typename decltype(matsetView)::IDList;
    using VFList = typename decltype(matsetView)::VFList;
    using VFType = typename VFList::value_type;

    EXPECT_EQ(matsetView.numberOfZones(), 3);
    EXPECT_EQ(matsetView.numberOfMaterials(0), 2);
    EXPECT_EQ(matsetView.numberOfMaterials(1), 2);
    EXPECT_EQ(matsetView.numberOfMaterials(2), 1);

    IDList m0, m1, m2;
    VFList vf0, vf1, vf2;
    matsetView.zoneMaterials(0, m0, vf0);
    EXPECT_EQ(m0.size(), 2);
    EXPECT_EQ(vf0.size(), 2);
    EXPECT_EQ(m0[0], 0);
    EXPECT_EQ(m0[1], 1);
    EXPECT_EQ(vf0[0], 0.3);
    EXPECT_EQ(vf0[1], 0.7);

    VFType vf;
    EXPECT_TRUE(matsetView.zoneContainsMaterial(0, 0, vf));
    EXPECT_EQ(vf, 0.3);
    EXPECT_TRUE(matsetView.zoneContainsMaterial(0, 1, vf));
    EXPECT_EQ(vf, 0.7);

    matsetView.zoneMaterials(1, m1, vf1);
    EXPECT_EQ(m1.size(), 2);
    EXPECT_EQ(vf1.size(), 2);
    EXPECT_EQ(m1[0], 0);
    EXPECT_EQ(m1[1], 1);
    EXPECT_EQ(vf1[0], 0.33);
    EXPECT_EQ(vf1[1], 0.67);

    EXPECT_TRUE(matsetView.zoneContainsMaterial(1, 0, vf));
    EXPECT_EQ(vf, 0.33);
    EXPECT_TRUE(matsetView.zoneContainsMaterial(1, 1, vf));
    EXPECT_EQ(vf, 0.67);

    matsetView.zoneMaterials(2, m2, vf2);
    EXPECT_EQ(m2.size(), 1);
    EXPECT_EQ(m2[0], 1);
    EXPECT_EQ(vf2[0], 1);

    EXPECT_FALSE(matsetView.zoneContainsMaterial(2, 0, vf));
    EXPECT_EQ(vf, 0.);
    EXPECT_TRUE(matsetView.zoneContainsMaterial(2, 1, vf));
    EXPECT_EQ(vf, 1.);
  });
}

//------------------------------------------------------------------------------
#if defined(DEBUGGING_TEST_CASES)
void conduit_debug_err_handler(const std::string &s1, const std::string &s2, int i1)
{
  std::cout << "s1=" << s1 << ", s2=" << s2 << ", i1=" << i1 << std::endl;
  // This is on purpose.
  while(1)
    ;
}
#endif
//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger,
#if defined(DEBUGGING_TEST_CASES)
  conduit::utils::set_error_handler(conduit_debug_err_handler);
#endif
  result = RUN_ALL_TESTS();
  return result;
}
