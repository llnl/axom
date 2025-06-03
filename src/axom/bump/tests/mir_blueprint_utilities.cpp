// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/mir.hpp"

#include "axom/mir/tests/mir_testing_helpers.hpp"
#include "axom/mir/tests/mir_testing_data_helpers.hpp"

#include <iostream>
#include <algorithm>

namespace mir = axom::mir;
namespace views = axom::mir::views;
namespace bputils = axom::mir::utilities::blueprint;

std::string baselineDirectory()
{
  return pjoin(dataDirectory(), "mir", "regression", "mir_blueprint_utilities");
}

//------------------------------------------------------------------------------
// Global test application object.
MIRTestApplication TestApp;

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_conduit_allocate
{
  static void test()
  {
    axom::mir::utilities::blueprint::ConduitAllocateThroughAxom<ExecSpace> c2a;
    EXPECT_TRUE(c2a.getConduitAllocatorID() > 0);

    constexpr int nValues = 100;
    conduit::Node n;
    n.set_allocator(c2a.getConduitAllocatorID());
    n.set(conduit::DataType::int32(nValues));

    // Make sure we can store some values into the data that were allocated.
    auto nview = bputils::make_array_view<int>(n);
    axom::for_all<ExecSpace>(
      nValues,
      AXOM_LAMBDA(axom::IndexType index) { nview[index] = index; });

    EXPECT_EQ(n.dtype().number_of_elements(), nValues);

    // Get the values back to the host.
    std::vector<int> hostValues(nValues);
    axom::copy(hostValues.data(), n.data_ptr(), sizeof(int) * nValues);

    // Check that the values were set.
    for(int i = 0; i < nValues; i++)
    {
      EXPECT_EQ(hostValues[i], i);
    }
  }
};

TEST(mir_blueprint_utilities, allocate_seq) { test_conduit_allocate<seq_exec>::test(); }
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(mir_blueprint_utilities, allocate_omp) { test_conduit_allocate<omp_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(mir_blueprint_utilities, allocate_cuda) { test_conduit_allocate<cuda_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(mir_blueprint_utilities, allocate_hip) { test_conduit_allocate<hip_exec>::test(); }
#endif

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_copy_braid
{
  static void test()
  {
    conduit::Node hostMesh;
    create(hostMesh);

    // host->device
    conduit::Node deviceMesh;
    axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);

    // Run some minmax operations on device (proves that the data was in the right place) and check the results.

    constexpr double eps = 1.e-7;

    auto x = axom::mir::utilities::blueprint::MinMax<ExecSpace, double>::execute(
      deviceMesh["coordsets/coords/values/x"]);
    //std::cout << std::setw(16) << "x={" << x.first << ", " << x.second << "}\n";
    EXPECT_NEAR(x.first, -10., eps);
    EXPECT_NEAR(x.second, 10., eps);

    auto y = axom::mir::utilities::blueprint::MinMax<ExecSpace, double>::execute(
      deviceMesh["coordsets/coords/values/y"]);
    //std::cout << std::setw(16) << "y={" << y.first << ", " << y.second << "}\n";
    EXPECT_NEAR(y.first, -10., eps);
    EXPECT_NEAR(y.second, 10., eps);

    auto c = axom::mir::utilities::blueprint::MinMax<ExecSpace, double>::execute(
      deviceMesh["topologies/mesh/elements/connectivity"]);
    //std::cout << std::setw(16) << "conn={" << c.first << ", " << c.second << "}\n";
    EXPECT_NEAR(c.first, 0., eps);
    EXPECT_NEAR(c.second, 999., eps);

    auto r = axom::mir::utilities::blueprint::MinMax<ExecSpace, double>::execute(
      deviceMesh["fields/radial/values"]);
    //std::cout << std::setw(16) << "radial={" << r.first << ", " << r.second << "}\n";
    EXPECT_NEAR(r.first, 19.2450089729875, eps);
    EXPECT_NEAR(r.second, 173.205080756888, eps);
  }

  static void create(conduit::Node &mesh)
  {
    const int d[3] = {10, 10, 10};
    conduit::blueprint::mesh::examples::braid("hexs", d[0], d[1], d[2], mesh);
  }
};

TEST(mir_blueprint_utilities, copy_seq) { test_copy_braid<seq_exec>::test(); }

#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(mir_blueprint_utilities, copy_omp) { test_copy_braid<omp_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(mir_blueprint_utilities, copy_cuda) { test_copy_braid<cuda_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(mir_blueprint_utilities, copy_hip) { test_copy_braid<hip_exec>::test(); }
#endif

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_make_unstructured
{
  static void test()
  {
    conduit::Node hostMesh;
    create(hostMesh);

    // host->device
    conduit::Node deviceMesh;
    bputils::copy<ExecSpace>(deviceMesh, hostMesh);

    // _mir_utilities_makeunstructured_begin
    conduit::Node deviceResult;
    bputils::MakeUnstructured<ExecSpace> uns;
    uns.execute(deviceMesh["topologies/mesh"], deviceMesh["coordsets/coords"], "mesh", deviceResult);
    // _mir_utilities_makeunstructured_end

    // device->host
    conduit::Node hostResult;
    bputils::copy<axom::SEQ_EXEC>(hostResult, deviceResult);

    TestApp.saveVisualization("unstructured", hostResult);

    // Handle baseline comparison.
    EXPECT_TRUE(TestApp.test<ExecSpace>("unstructured", hostResult));
  }

  static void create(conduit::Node &mesh)
  {
    std::vector<int> dims {4, 4};
    axom::mir::testing::data::braid("uniform", dims, mesh);
  }
};

TEST(mir_blueprint_utilities, make_unstructured_seq) { test_make_unstructured<seq_exec>::test(); }
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(mir_blueprint_utilities, make_unstructured_omp) { test_make_unstructured<omp_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(mir_blueprint_utilities, make_unstructured_cuda) { test_make_unstructured<cuda_exec>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(mir_blueprint_utilities, make_unstructured_hip) { test_make_unstructured<hip_exec>::test(); }
#endif

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_recenter_field
{
  static void test()
  {
    conduit::Node hostMesh;
    create(hostMesh);

    // host -> device
    conduit::Node deviceMesh;
    axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);
    // _mir_utilities_recenterfield_begin
    const conduit::Node &deviceTopo = deviceMesh["topologies/mesh"];
    const conduit::Node &deviceCoordset = deviceMesh["coordsets/coords"];

    // Make a node to zone relation on the device.
    conduit::Node deviceRelation;
    axom::mir::utilities::blueprint::NodeToZoneRelationBuilder<ExecSpace> n2z;
    n2z.execute(deviceTopo, deviceCoordset, deviceRelation);

    // Recenter a field zonal->nodal on the device
    axom::mir::utilities::blueprint::RecenterField<ExecSpace> r;
    r.execute(deviceMesh["fields/easy_zonal"], deviceRelation, deviceMesh["fields/z2n"]);

    // Recenter a field nodal->zonal on the device. (The elements are an o2m relation)
    r.execute(deviceMesh["fields/z2n"],
              deviceMesh["topologies/mesh/elements"],
              deviceMesh["fields/n2z"]);
    // _mir_utilities_recenterfield_end

    // device -> host
    conduit::Node hostResultMesh;
    axom::mir::utilities::blueprint::copy<seq_exec>(hostResultMesh, deviceMesh);

    // Print the results.
    //printNode(hostResultMesh);

    const float n2z_result[] = {1., 2., 4., 5., 4., 5., 7., 8., 7., 8., 10., 11.};
    for(size_t i = 0; i < (sizeof(n2z_result) / sizeof(float)); i++)
    {
      EXPECT_EQ(n2z_result[i], hostResultMesh["fields/z2n/values"].as_float_accessor()[i]);
    }
    const float z2n_result[] = {3., 4.5, 6., 6., 7.5, 9.};
    for(size_t i = 0; i < (sizeof(z2n_result) / sizeof(float)); i++)
    {
      EXPECT_EQ(z2n_result[i], hostResultMesh["fields/n2z/values"].as_float_accessor()[i]);
    }
  }

  static void create(conduit::Node &mesh)
  {
    /*
      8---9--10--11
      |   |   |   |
      4---5---6---7
      |   |   |   |
      0---1---2---3
      */
    axom::StackArray<int, 2> dims {{4, 3}};
    axom::mir::testing::data::braid("quads", dims, mesh);
    mesh["topologies/mesh/elements/sizes"].set(std::vector<int> {{4, 4, 4, 4, 4, 4}});
    mesh["topologies/mesh/elements/offsets"].set(std::vector<int> {{0, 4, 8, 12, 16, 20}});
    mesh["fields/easy_zonal/topology"] = "mesh";
    mesh["fields/easy_zonal/association"] = "element";
    mesh["fields/easy_zonal/values"].set(std::vector<float> {{1, 3, 5, 7, 9, 11}});
  }
};

TEST(mir_blueprint_utilities, recenterfield_seq) { test_recenter_field<seq_exec>::test(); }

#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(mir_blueprint_utilities, recenterfield_omp) { test_recenter_field<omp_exec>::test(); }
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(mir_blueprint_utilities, recenterfield_cuda) { test_recenter_field<cuda_exec>::test(); }
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(mir_blueprint_utilities, recenterfield_hip) { test_recenter_field<hip_exec>::test(); }
#endif

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_extractzones
{
  static void test()
  {
    conduit::Node hostMesh;
    create(hostMesh);

    // host->device
    conduit::Node deviceMesh;
    bputils::copy<ExecSpace>(deviceMesh, hostMesh);

    axom::Array<axom::IndexType> ids {{1, 3, 5}};
    const auto nzones = ids.size();
    axom::Array<axom::IndexType> selectedZones(nzones,
                                               nzones,
                                               axom::execution_space<ExecSpace>::allocatorID());
    axom::copy(selectedZones.data(), ids.data(), nzones * sizeof(axom::IndexType));

    // Wrap the data in views.
    auto coordsetView =
      views::make_explicit_coordset<conduit::float64, 2>::view(deviceMesh["coordsets/coords"]);
    using CoordsetView = decltype(coordsetView);

    using TopologyView = views::UnstructuredTopologySingleShapeView<views::QuadShape<conduit::int64>>;
    TopologyView topoView(
      bputils::make_array_view<conduit::int64>(deviceMesh["topologies/mesh/elements/connectivity"]),
      bputils::make_array_view<conduit::int64>(deviceMesh["topologies/mesh/elements/sizes"]),
      bputils::make_array_view<conduit::int64>(deviceMesh["topologies/mesh/elements/offsets"]));

    // Pull out selected zones
    bputils::ExtractZones<ExecSpace, TopologyView, CoordsetView> extract(topoView, coordsetView);
    conduit::Node options, newDeviceMesh;
    options["topology"] = "mesh";
    extract.execute(selectedZones.view(), deviceMesh, options, newDeviceMesh);

    // device->host
    conduit::Node newHostMesh;
    bputils::copy<axom::SEQ_EXEC>(newHostMesh, newDeviceMesh);

    //printNode(newHostMesh);

    // Check some of the key arrays
    const axom::Array<conduit::int64> connectivity {{0, 1, 4, 3, 2, 3, 7, 6, 4, 5, 9, 8}};
    const axom::Array<conduit::int64> sizes {{4, 4, 4}};
    const axom::Array<conduit::int64> offsets {{0, 4, 8}};
    const axom::Array<conduit::float64> x {{1.0, 2.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0}};
    const axom::Array<conduit::float64> y {{0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0}};
    const axom::Array<conduit::float64> zonal {{1.0, 3.0, 5.0}};
    const axom::Array<conduit::float64> nodal {{1.0, 2.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0}};

    EXPECT_TRUE(compare_views(connectivity.view(),
                              bputils::make_array_view<conduit::int64>(
                                newHostMesh["topologies/mesh/elements/connectivity"])));
    EXPECT_TRUE(compare_views(
      sizes.view(),
      bputils::make_array_view<conduit::int64>(newHostMesh["topologies/mesh/elements/sizes"])));
    EXPECT_TRUE(compare_views(
      offsets.view(),
      bputils::make_array_view<conduit::int64>(newHostMesh["topologies/mesh/elements/offsets"])));
    EXPECT_TRUE(compare_views(
      x.view(),
      bputils::make_array_view<conduit::float64>(newHostMesh["coordsets/coords/values/x"])));
    EXPECT_TRUE(compare_views(
      y.view(),
      bputils::make_array_view<conduit::float64>(newHostMesh["coordsets/coords/values/y"])));
    EXPECT_TRUE(compare_views(
      zonal.view(),
      bputils::make_array_view<conduit::float64>(newHostMesh["fields/zonal/values"])));
    EXPECT_TRUE(compare_views(
      nodal.view(),
      bputils::make_array_view<conduit::float64>(newHostMesh["fields/nodal/values"])));

    // Do the material too.
    using MatsetView = views::UnibufferMaterialView<conduit::int64, conduit::float64, 3>;
    MatsetView matsetView;
    matsetView.set(
      bputils::make_array_view<conduit::int64>(deviceMesh["matsets/mat1/material_ids"]),
      bputils::make_array_view<conduit::float64>(deviceMesh["matsets/mat1/volume_fractions"]),
      bputils::make_array_view<conduit::int64>(deviceMesh["matsets/mat1/sizes"]),
      bputils::make_array_view<conduit::int64>(deviceMesh["matsets/mat1/offsets"]),
      bputils::make_array_view<conduit::int64>(deviceMesh["matsets/mat1/indices"]));

    // Pull out selected zones
    bputils::ExtractZonesAndMatset<ExecSpace, TopologyView, CoordsetView, MatsetView> extractM(
      topoView,
      coordsetView,
      matsetView);
    newDeviceMesh.reset();
    extractM.execute(selectedZones.view(), deviceMesh, options, newDeviceMesh);

    // device->host
    newHostMesh.reset();
    bputils::copy<axom::SEQ_EXEC>(newHostMesh, newDeviceMesh);

    // Check some of the key arrays in the sliced material
    const axom::Array<conduit::int64> mat_sizes {{2, 1, 2}};
    const axom::Array<conduit::int64> mat_offsets {{0, 2, 3}};
    const axom::Array<conduit::int64> mat_indices {{0, 1, 2, 3, 4}};
    const axom::Array<conduit::int64> mat_material_ids {{1, 2, 2, 2, 3}};
    const axom::Array<conduit::float64> mat_volume_fractions {{0.5, 0.5, 1.0, 0.8, 0.2}};

    EXPECT_TRUE(
      compare_views(mat_sizes.view(),
                    bputils::make_array_view<conduit::int64>(newHostMesh["matsets/mat1/sizes"])));
    EXPECT_TRUE(
      compare_views(mat_offsets.view(),
                    bputils::make_array_view<conduit::int64>(newHostMesh["matsets/mat1/offsets"])));
    EXPECT_TRUE(
      compare_views(mat_indices.view(),
                    bputils::make_array_view<conduit::int64>(newHostMesh["matsets/mat1/indices"])));
    EXPECT_TRUE(compare_views(
      mat_material_ids.view(),
      bputils::make_array_view<conduit::int64>(newHostMesh["matsets/mat1/material_ids"])));
    EXPECT_TRUE(compare_views(
      mat_volume_fractions.view(),
      bputils::make_array_view<conduit::float64>(newHostMesh["matsets/mat1/volume_fractions"])));
  }

  static void create(conduit::Node &hostMesh)
  {
    /*
      8-------9------10------11
      |  2/1  | 1/0.1 | 2/0.8 |
      |       | 2/0.5 | 3/0.2 |
      |       | 3/0.4 |       |
      4-------5-------6-------7
      |       | 1/0.5 | 1/0.2 |
      |  1/1  | 2/0.5 | 2/0.8 |
      |       |       |       |
      0-------1-------2-------3
      */
    const char *yaml = R"xx(
coordsets:
  coords:
    type: explicit
    values:
      x: [0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3.]
      y: [0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2.]
topologies:
  mesh:
    coordset: coords
    type: unstructured
    elements:
      shape: quad
      connectivity: [0,1,5,4, 1,2,6,5, 2,3,7,6, 4,5,9,8, 5,6,10,9, 6,7,11,10]
      sizes: [4,4,4,4,4,4]
      offsets: [0,4,8,12,16,20]
fields:
  zonal:
    topology: mesh
    association: element
    values: [0.,1.,2.,3.,4.,5.]
  nodal:
    topology: mesh
    association: vertex
    values: [0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.]
matsets:
  mat1:
    topology: mesh
    material_map:
      a: 1
      b: 2
      c: 3
    material_ids: [1, 1,2, 1,2, 2, 1,2,3, 2,3]
    volume_fractions: [1., 0.5,0.5, 0.2,0.8, 1., 0.1,0.5,0.4, 0.8,0.2]
    indices: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    sizes: [1, 2, 2, 1, 3, 2]
    offsets: [0, 1, 3, 5, 6, 9]
)xx";

    hostMesh.parse(yaml);
  }
};

TEST(mir_blueprint_utilities, extractzones_seq) { test_extractzones<seq_exec>::test(); }

#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(mir_blueprint_utilities, extractzones_omp) { test_extractzones<omp_exec>::test(); }
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(mir_blueprint_utilities, extractzones_cuda) { test_extractzones<cuda_exec>::test(); }
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(mir_blueprint_utilities, extractzones_hip) { test_extractzones<hip_exec>::test(); }
#endif

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_extractzones_polyhedral
{
  static void test(const std::string &name, bool selectZones)
  {
    constexpr int MAXMATERIALS = 5;
    const int gridSize = 7;
    const int numCircles = 2;
    SLIC_ASSERT(numCircles + 1 <= MAXMATERIALS);

    conduit::Node hostMesh;
    create(gridSize, numCircles, hostMesh);

    // Save visualization, if enabled.
    TestApp.saveVisualization(name + "_orig", hostMesh);

    // host->device
    conduit::Node deviceMesh;
    bputils::copy<ExecSpace>(deviceMesh, hostMesh);

    // Make selected zones.
    axom::Array<axom::IndexType> ids;
    for(int k = 0; k < gridSize; k++)
    {
      for(int j = 0; j < gridSize; j++)
      {
        for(int i = 0; i < gridSize; i++)
        {
          // Halfway up on Y, checkerboard the zones
          if(selectZones && (j > gridSize / 2) && ((i + j + k) % 2 == 0))
          {
            continue;
          }

          int zoneIndex = k * gridSize * gridSize + j * gridSize + i;
          ids.push_back(zoneIndex);
        }
      }
    }

    // Put selected zones on device.
    const auto nzones = ids.size();
    axom::Array<axom::IndexType> selectedZones(nzones,
                                               nzones,
                                               axom::execution_space<ExecSpace>::allocatorID());
    axom::copy(selectedZones.data(), ids.data(), nzones * sizeof(axom::IndexType));

    const conduit::Node &n_coordset = deviceMesh["coordsets/coords"];
    const conduit::Node &n_topology = deviceMesh["topologies/mesh"];
    const conduit::Node &n_matset = deviceMesh["matsets/mat"];

    // Wrap the data in views.
    auto coordsetView = views::make_explicit_coordset<float, 3>::view(n_coordset);
    using CoordsetView = decltype(coordsetView);

    auto topologyView = views::make_structured_topology<3>::view(n_topology);
    using TopologyView = decltype(topologyView);
    using IndexingPolicy = typename TopologyView::IndexingPolicy;

    auto matsetView = views::make_unibuffer_matset<int, float, MAXMATERIALS>::view(n_matset);
    using MatsetView = decltype(matsetView);

    // Pull out selected zones as polyhedral zones
    bputils::ExtractZonesAndMatsetPolyhedral<ExecSpace, IndexingPolicy, CoordsetView, MatsetView>
      extract(topologyView, coordsetView, matsetView);
    conduit::Node newDeviceMesh, options;
    extract.execute(selectedZones.view(), deviceMesh, options, newDeviceMesh);

    // device->host
    conduit::Node newHostMesh;
    bputils::copy<axom::SEQ_EXEC>(newHostMesh, newDeviceMesh);

    // Save visualization, if enabled.
    TestApp.saveVisualization(name, newHostMesh);

    // Test against baseline.
    EXPECT_TRUE(TestApp.test<ExecSpace>(name, newHostMesh));
  }

  static void create(int gridSize, int numCircles, conduit::Node &hostMesh)
  {
    AXOM_ANNOTATE_SCOPE("generate");
    mir::MeshTester tester;
    tester.setStructured(true);
    tester.initTestCaseSix(gridSize, numCircles, hostMesh);
  }
};

TEST(mir_blueprint_utilities, extractzones_polyhedral_seq)
{
  test_extractzones_polyhedral<seq_exec>::test("extractzones_polyhedral", false);
}
TEST(mir_blueprint_utilities, extractzones_polyhedral_sel_seq)
{
  test_extractzones_polyhedral<seq_exec>::test("extractzones_polyhedral_sel", true);
}

#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(mir_blueprint_utilities, extractzones_polyhedral_omp)
{
  test_extractzones_polyhedral<omp_exec>::test("extractzones_polyhedral", false);
}
TEST(mir_blueprint_utilities, extractzones_polyhedral_sel_omp)
{
  test_extractzones_polyhedral<omp_exec>::test("extractzones_polyhedral_sel", true);
}
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(mir_blueprint_utilities, extractzones_polyhedral_cuda)
{
  test_extractzones_polyhedral<cuda_exec>::test("extractzones_polyhedral", false);
}
TEST(mir_blueprint_utilities, extractzones_polyhedral_sel_cuda)
{
  test_extractzones_polyhedral<cuda_exec>::test("extractzones_polyhedral_sel", true);
}
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(mir_blueprint_utilities, extractzones_polyhedral_hip)
{
  test_extractzones_polyhedral<hip_exec>::test("extractzones_polyhedral", false);
}
TEST(mir_blueprint_utilities, extractzones_polyhedral_sel_hip)
{
  test_extractzones_polyhedral<hip_exec>::test("extractzones_polyhedral_sel", true);
}
#endif

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_zonelistbuilder
{
  static void test()
  {
    conduit::Node hostMesh;
    create(hostMesh);

    // host->device
    conduit::Node deviceMesh;
    bputils::copy<ExecSpace>(deviceMesh, hostMesh);

    // Wrap the data in views.
    auto coordsetView =
      views::make_rectilinear_coordset<conduit::float64, 2>::view(deviceMesh["coordsets/coords"]);

    auto topologyView = views::make_rectilinear_topology<2>::view(deviceMesh["topologies/mesh"]);
    using TopologyView = decltype(topologyView);

    // Do the material too.
    using MatsetView = views::UnibufferMaterialView<conduit::int64, conduit::float64, 2>;
    MatsetView matsetView;
    matsetView.set(
      bputils::make_array_view<conduit::int64>(deviceMesh["matsets/mat1/material_ids"]),
      bputils::make_array_view<conduit::float64>(deviceMesh["matsets/mat1/volume_fractions"]),
      bputils::make_array_view<conduit::int64>(deviceMesh["matsets/mat1/sizes"]),
      bputils::make_array_view<conduit::int64>(deviceMesh["matsets/mat1/offsets"]),
      bputils::make_array_view<conduit::int64>(deviceMesh["matsets/mat1/indices"]));

    // Determine the list of clean and mixed zones (taking into account #mats at the nodes)
    bputils::ZoneListBuilder<ExecSpace, TopologyView, MatsetView> zlb(topologyView, matsetView);
    axom::Array<axom::IndexType> clean, mixed;
    zlb.execute(coordsetView.numberOfNodes(), clean, mixed);

    conduit::Node deviceData;
    deviceData["clean"].set_external(clean.data(), clean.size());
    deviceData["mixed"].set_external(mixed.data(), mixed.size());

    // device->host
    conduit::Node hostData;
    bputils::copy<axom::SEQ_EXEC>(hostData, deviceData);

    // Compare expected
    const axom::Array<axom::IndexType> cleanResult {{0, 1, 2, 3, 4, 8, 12}};
    const axom::Array<axom::IndexType> mixedResult {{5, 6, 7, 9, 10, 11, 13, 14, 15}};
    EXPECT_TRUE(compare_views(cleanResult.view(),
                              bputils::make_array_view<axom::IndexType>(hostData["clean"])));
    EXPECT_TRUE(compare_views(mixedResult.view(),
                              bputils::make_array_view<axom::IndexType>(hostData["mixed"])));

    // Try selecting a subset of the zones.
    axom::Array<axom::IndexType> ids {{2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 14}};
    axom::Array<axom::IndexType> selectedZones(ids.size(),
                                               ids.size(),
                                               axom::execution_space<ExecSpace>::allocatorID());
    axom::copy(selectedZones.data(), ids.data(), ids.size() * sizeof(axom::IndexType));
    zlb.execute(coordsetView.numberOfNodes(), selectedZones.view(), clean, mixed);

    deviceData["clean"].set_external(clean.data(), clean.size());
    deviceData["mixed"].set_external(mixed.data(), mixed.size());

    // device->host
    bputils::copy<axom::SEQ_EXEC>(hostData, deviceData);

    // Compare expected
    const axom::Array<axom::IndexType> cleanResult2 {{2, 3, 8, 12}};
    const axom::Array<axom::IndexType> mixedResult2 {{6, 7, 9, 10, 11, 13, 14}};
    EXPECT_TRUE(compare_views(cleanResult2.view(),
                              bputils::make_array_view<axom::IndexType>(hostData["clean"])));
    EXPECT_TRUE(compare_views(mixedResult2.view(),
                              bputils::make_array_view<axom::IndexType>(hostData["mixed"])));
  }

  static void create(conduit::Node &hostMesh)
  {
    /*
    20------21-------22-------23-------24
    |  1/1  |  1/1   |  1/.5  |  2/1.  |    1/1 = mat#1, vf=1.0
    |       |        |  2/.5  |        |
    |z12    |z13     |z14     |z15     |
    15------16-------17-------18-------19
    |  1/1  |  1/1   |  1/0.7 |  1/.5  |
    |       |        |  2/0.3 |  2/.5  |
    |z8     |z9      |z10     |z11     |
    10------11-------12-------13-------14
    |  1/1  |  1/1   |  1/1   |  1/1   |
    |       |        |        |        |
    |z4     |z5      |z6      |z7      |
    5-------6--------7--------8--------9
    |  1/1  |  1/1   |  1/1   |  1/1   |
    |       |        |        |        |
    |z0     |z1      |z2      |z3      |
    0-------1--------2--------3--------4
    */
    const char *yaml = R"xx(
coordsets:
  coords:
    type: rectilinear
    values:
      x: [0., 1., 2., 3., 4.]
      y: [0., 1., 2., 3., 4.]
topologies:
  mesh:
    type: rectilinear
    coordset: coords
matsets:
  mat1:
    topology: mesh
    material_map:
      a: 1
      b: 2
    material_ids: [1, 1, 1, 1,    1, 1, 1, 1,   1, 1, 1, 2, 1, 2,   1, 1, 1, 2, 2]
    volume_fractions: [1., 1., 1., 1.,    1., 1., 1., 1.,   1., 1., 0.7, 0.3, .5, 0.5,   1., 1., 0.5, 0.5, 1.]
    indices: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
    sizes: [1, 1, 1, 1,   1, 1, 1, 1,   1, 1, 2, 2,   1, 1, 2, 1]
    offsets: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15, 16, 18]
)xx";

    hostMesh.parse(yaml);
  }
};

TEST(mir_blueprint_utilities, zonelistbuilder_seq)
{
  AXOM_ANNOTATE_SCOPE("zonelistbuilder_seq");
  test_zonelistbuilder<seq_exec>::test();
}
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(mir_blueprint_utilities, zonelistbuilder_omp)
{
  AXOM_ANNOTATE_SCOPE("zonelistbuilder_omp");
  test_zonelistbuilder<omp_exec>::test();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(mir_blueprint_utilities, zonelistbuilder_cuda)
{
  AXOM_ANNOTATE_SCOPE("zonelistbuilder_cuda");
  test_zonelistbuilder<cuda_exec>::test();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(mir_blueprint_utilities, zonelistbuilder_hip)
{
  AXOM_ANNOTATE_SCOPE("zonelistbuilder_hip");
  test_zonelistbuilder<hip_exec>::test();
}
#endif

//------------------------------------------------------------------------------

template <typename ExecSpace>
struct test_makezonecenters
{
  static void test()
  {
    conduit::Node hostMesh;
    create(hostMesh);

    // host->device
    conduit::Node deviceMesh;
    bputils::copy<ExecSpace>(deviceMesh, hostMesh);

    const conduit::Node &n_rmesh = deviceMesh["topologies/rmesh"];
    auto rmeshView = views::make_rectilinear_topology<2>::view(n_rmesh);
    testTopo(deviceMesh, rmeshView, n_rmesh);

    const conduit::Node &n_umesh = deviceMesh["topologies/umesh"];
    views::UnstructuredTopologySingleShapeView<views::QuadShape<conduit::index_t>> umeshView(
      bputils::make_array_view<conduit::index_t>(n_umesh["elements/connectivity"]),
      bputils::make_array_view<conduit::index_t>(n_umesh["elements/sizes"]),
      bputils::make_array_view<conduit::index_t>(n_umesh["elements/offsets"]));
    testTopo(deviceMesh, umeshView, n_umesh);
  }

  template <typename TopologyView>
  static void testTopo(const conduit::Node &deviceMesh,
                       const TopologyView &topoView,
                       const conduit::Node &n_topo)
  {
    const conduit::Node &n_coordset = deviceMesh["coordsets/coords"];
    auto coordsetView = views::make_rectilinear_coordset<double, 2>::view(n_coordset);
    using CoordsetView = decltype(coordsetView);

    bputils::MakeZoneCenters<ExecSpace, TopologyView, CoordsetView> zc(topoView, coordsetView);
    conduit::Node n_field;
    zc.execute(n_topo, n_coordset, n_field);

    // device->host
    conduit::Node n_hostField;
    bputils::copy<axom::SEQ_EXEC>(n_hostField, n_field);

    const double eps = 1.e-9;
    const double res_x[] = {0.5, 1.5, 2.5, 0.5, 1.5, 2.5, 0.5, 1.5, 2.5};
    const double res_y[] = {0.5, 0.5, 0.5, 1.5, 1.5, 1.5, 2.5, 2.5, 2.5};
    const auto x = n_hostField["values/x"].as_double_accessor();
    const auto y = n_hostField["values/y"].as_double_accessor();
    EXPECT_EQ(x.number_of_elements(), 9);
    EXPECT_EQ(y.number_of_elements(), 9);
    for(int i = 0; i < 9; i++)
    {
      EXPECT_NEAR(x[i], res_x[i], eps);
      EXPECT_NEAR(y[i], res_y[i], eps);
    }
  }

  static void create(conduit::Node &hostMesh)
  {
    /*
    12------13-------14-------15
    |       |        |        |
    |       |        |        |
    |z6     |z7      |z8      |
    8-------9--------10-------11
    |       |        |        |
    |       |        |        |
    |z3     |z4      |z5      |
    4-------5--------6--------7
    |       |        |        |
    |       |        |        |
    |z0     |z1      |z2      |
    0-------1--------2--------3
    */
    const char *yaml = R"xx(
coordsets:
  coords:
    type: rectilinear
    values:
      x: [0., 1., 2., 3.]
      y: [0., 1., 2., 3.]
topologies:
  rmesh:
    type: rectilinear
    coordset: coords
  umesh:
    type: unstructured
    coordset: coords
    elements:
      shape: "quad"
      connectivity: [0,1,5,4, 1,2,6,5, 2,3,7,6, 4,5,9,8, 5,6,10,9, 6,7,11,10, 8,9,13,12, 9,10,14,13, 10,11,15,14]
      sizes: [4,4,4,4,4,4,4,4,4]
      offsets: [0,4,8,12,16,20,24,28,32]
)xx";

    hostMesh.parse(yaml);
  }
};

TEST(mir_blueprint_utilities, makezonecenters_seq)
{
  AXOM_ANNOTATE_SCOPE("makezonecenters_seq");
  test_makezonecenters<seq_exec>::test();
}
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(mir_blueprint_utilities, makezonecenters_omp)
{
  AXOM_ANNOTATE_SCOPE("makezonecenters_omp");
  test_makezonecenters<omp_exec>::test();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(mir_blueprint_utilities, makezonecenters_cuda)
{
  AXOM_ANNOTATE_SCOPE("makezonecenters_cuda");
  test_makezonecenters<cuda_exec>::test();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(mir_blueprint_utilities, makezonecenters_hip)
{
  AXOM_ANNOTATE_SCOPE("makezonecenters_hip");
  test_makezonecenters<hip_exec>::test();
}
#endif

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_mergecoordsetpoints
{
  static void test()
  {
    conduit::Node hostMesh;
    create(hostMesh);

    // host->device
    conduit::Node deviceMesh;
    bputils::copy<ExecSpace>(deviceMesh, hostMesh);

    conduit::Node &n_coordset = deviceMesh["coordsets/coords"];
    auto coordsetView = views::make_explicit_coordset<double, 2>::view(n_coordset);
    using CoordsetView = decltype(coordsetView);

    bputils::MergeCoordsetPoints<ExecSpace, CoordsetView> mcp(coordsetView);
    conduit::Node n_newCoordset;
    axom::Array<axom::IndexType> selectedIds, old2new;
    // Make anything closer than 0.005 match
    const double eps = 0.005;
    conduit::Node n_options;
    n_options["tolerance"] = eps;
    mcp.execute(n_coordset, n_options, selectedIds, old2new);

    // Stash the selectedIds and old2new in the coordset node so we can bring it to the host easier.
    n_coordset["selectedIds"].set_external(selectedIds.data(), selectedIds.size());
    n_coordset["old2new"].set_external(old2new.data(), old2new.size());

    // device->host
    conduit::Node n_hostCoordset;
    bputils::copy<axom::SEQ_EXEC>(n_hostCoordset, n_coordset);
    //printNode(n_hostCoordset);

    // Compare results.
    const double res_x[] = {2.0, -0.0001, 1.0001, 0.0, 2.0, 1.0, 2.0, 1.0, 0.0};
    const double res_y[] = {2.0, 1.0, 1.0, 0.0001, -0.0001, 0.0, 1.0001, 2.0001, 2.0};
    const int res_old2new[] = {3, 5, 2, 1, 5, 4, 6, 2, 1, 2, 7, 8, 2, 6, 0, 7};

    const auto x = n_hostCoordset["values/x"].as_double_accessor();
    const auto y = n_hostCoordset["values/y"].as_double_accessor();
    const auto o2n = n_hostCoordset["old2new"].as_int_accessor();
    EXPECT_EQ(x.number_of_elements(), 9);
    EXPECT_EQ(y.number_of_elements(), 9);
    EXPECT_EQ(o2n.number_of_elements(), 16);
    for(int i = 0; i < 9; i++)
    {
      EXPECT_NEAR(x[i], res_x[i], eps);
      EXPECT_NEAR(y[i], res_y[i], eps);
    }
    for(int i = 0; i < 16; i++)
    {
      EXPECT_EQ(o2n[i], res_old2new[i]);
    }
  }

  static void create(conduit::Node &hostMesh)
  {
    /*
    We have nodes that are given such that each zone corner is repeated and may have some
    small tolerance variations. We want to make sure the nodes get joined.

    8-------9--------10
    |       |        |
    |       |        |
    |z3     |z4      |
    4-------5--------6
    |       |        |
    |       |        |
    |z0     |z1      |
    0-------1--------2
    */
    const char *yaml = R"xx(
coordsets:
  coords:
    type: explicit
    values:
      x: [0.,     1.0001, 1., 0.,       1.,  2.,     2., 1.0001,       -0.0001, 1., 1.,     0.,       1.0001, 2.,     2., 1.]
      y: [0.0001, 0.,     1., 1.,       0., -0.0001, 1., 1.,            1.,     1., 2.0001, 2.,       1.,     1.0001, 2., 2.0001]
)xx";

    hostMesh.parse(yaml);
  }
};

TEST(mir_blueprint_utilities, mergecoordsetpoints_seq)
{
  AXOM_ANNOTATE_SCOPE("mergecoordsetpoints_seq");
  test_mergecoordsetpoints<seq_exec>::test();
}
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(mir_blueprint_utilities, mergecoordsetpoints_omp)
{
  AXOM_ANNOTATE_SCOPE("mergecoordsetpoints_omp");
  test_mergecoordsetpoints<omp_exec>::test();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(mir_blueprint_utilities, mergecoordsetpoints_cuda)
{
  AXOM_ANNOTATE_SCOPE("mergecoordsetpoints_cuda");
  test_mergecoordsetpoints<cuda_exec>::test();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(mir_blueprint_utilities, mergecoordsetpoints_hip)
{
  AXOM_ANNOTATE_SCOPE("mergecoordsetpoints_hip");
  test_mergecoordsetpoints<hip_exec>::test();
}
#endif

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_makepointmesh
{
  static void test()
  {
    conduit::Node hostMesh;
    create(hostMesh);

    // host->device
    conduit::Node deviceMesh;
    bputils::copy<ExecSpace>(deviceMesh, hostMesh);

    const conduit::Node &n_coordset = deviceMesh["coordsets/coords"];
    const conduit::Node &n_topology = deviceMesh["topologies/mesh"];

    // Wrap the data in views.
    auto coordsetView = views::make_explicit_coordset<conduit::float64, 2>::view(n_coordset);
    using CoordsetView = decltype(coordsetView);
    auto topoView =
      views::make_unstructured_single_shape_topology<views::QuadShape<conduit::int64>>::view(
        n_topology);
    using TopologyView = decltype(topoView);

    bputils::MakePointMesh<ExecSpace, TopologyView, CoordsetView> pm(topoView, coordsetView);
    conduit::Node options, newDeviceMesh;
    options["topologyName"] = "pointmesh";
    options["coordsetName"] = "points";
    for(int i = 0; i < 2; i++)
    {
      if(i == 1)
      {
        axom::Array<axom::IndexType> ids {{1, 3, 5}};
        const auto nzones = ids.size();
        axom::Array<axom::IndexType> selectedZones(nzones,
                                                   nzones,
                                                   axom::execution_space<ExecSpace>::allocatorID());
        axom::copy(selectedZones.data(), ids.data(), nzones * sizeof(axom::IndexType));

        pm.execute(selectedZones.view(), n_topology, n_coordset, options, newDeviceMesh);

        // device->host
        conduit::Node newHostMesh;
        bputils::copy<axom::SEQ_EXEC>(newHostMesh, newDeviceMesh);

        EXPECT_TRUE(newHostMesh.has_path("coordsets/points"));
        EXPECT_TRUE(newHostMesh.has_path("topologies/pointmesh"));
        EXPECT_EQ(newHostMesh["topologies/pointmesh/elements/shape"].as_string(), "point");

        // Compare answer
        const axom::Array<conduit::float64> x {{1.5, 0.5, 2.5}};
        const axom::Array<conduit::float64> y {{0.5, 1.5, 1.5}};
        const axom::Array<conduit::int64> connectivity {{0, 1, 2}};
        const axom::Array<conduit::int64> sizes {{1, 1, 1}};
        const axom::Array<conduit::int64> offsets {{0, 1, 2}};
        compare(newHostMesh, x, y, connectivity, sizes, offsets);
      }
      else
      {
        pm.execute(n_topology, n_coordset, options, newDeviceMesh);

        // device->host
        conduit::Node newHostMesh;
        bputils::copy<axom::SEQ_EXEC>(newHostMesh, newDeviceMesh);

        EXPECT_TRUE(newHostMesh.has_path("coordsets/points"));
        EXPECT_TRUE(newHostMesh.has_path("topologies/pointmesh"));
        EXPECT_EQ(newHostMesh["topologies/pointmesh/elements/shape"].as_string(), "point");

        // Compare answer
        const axom::Array<conduit::float64> x {{0.5, 1.5, 2.5, 0.5, 1.5, 2.5}};
        const axom::Array<conduit::float64> y {{0.5, 0.5, 0.5, 1.5, 1.5, 1.5}};
        const axom::Array<conduit::int64> connectivity {{0, 1, 2, 3, 4, 5}};
        const axom::Array<conduit::int64> sizes {{1, 1, 1, 1, 1, 1}};
        const axom::Array<conduit::int64> offsets {{0, 1, 2, 3, 4, 5}};
        compare(newHostMesh, x, y, connectivity, sizes, offsets);
      }
    }
  }

  static void compare(const conduit::Node &n_mesh,
                      const axom::Array<conduit::float64> &x,
                      const axom::Array<conduit::float64> &y,
                      const axom::Array<conduit::int64> &connectivity,
                      const axom::Array<conduit::int64> &sizes,
                      const axom::Array<conduit::int64> &offsets)
  {
    EXPECT_TRUE(compare_views(
      x.view(),
      bputils::make_array_view<conduit::float64>(n_mesh["coordsets/points/values/x"])));
    EXPECT_TRUE(compare_views(
      y.view(),
      bputils::make_array_view<conduit::float64>(n_mesh["coordsets/points/values/y"])));
    EXPECT_TRUE(compare_views(connectivity.view(),
                              bputils::make_array_view<conduit::int64>(
                                n_mesh["topologies/pointmesh/elements/connectivity"])));
    EXPECT_TRUE(compare_views(
      sizes.view(),
      bputils::make_array_view<conduit::int64>(n_mesh["topologies/pointmesh/elements/sizes"])));
    EXPECT_TRUE(compare_views(
      offsets.view(),
      bputils::make_array_view<conduit::int64>(n_mesh["topologies/pointmesh/elements/offsets"])));
  }

  static void create(conduit::Node &hostMesh)
  {
    /*
      8-------9------10------11
      |  2/1  | 1/0.1 | 2/0.8 |
      |       | 2/0.5 | 3/0.2 |
      |       | 3/0.4 |       |
      4-------5-------6-------7
      |       | 1/0.5 | 1/0.2 |
      |  1/1  | 2/0.5 | 2/0.8 |
      |       |       |       |
      0-------1-------2-------3
      */
    const char *yaml = R"xx(
coordsets:
  coords:
    type: explicit
    values:
      x: [0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3.]
      y: [0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2.]
topologies:
  mesh:
    coordset: coords
    type: unstructured
    elements:
      shape: quad
      connectivity: [0,1,5,4, 1,2,6,5, 2,3,7,6, 4,5,9,8, 5,6,10,9, 6,7,11,10]
      sizes: [4,4,4,4,4,4]
      offsets: [0,4,8,12,16,20]
)xx";

    hostMesh.parse(yaml);
  }
};

TEST(mir_blueprint_utilities, makepointmesh_seq)
{
  AXOM_ANNOTATE_SCOPE("makepointmesh_seq");
  test_makepointmesh<seq_exec>::test();
}
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(mir_blueprint_utilities, makepointmesh_omp)
{
  AXOM_ANNOTATE_SCOPE("makepointmesh_omp");
  test_makepointmesh<omp_exec>::test();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(mir_blueprint_utilities, makepointmesh_cuda)
{
  AXOM_ANNOTATE_SCOPE("makepointmesh_cuda");
  test_makepointmesh<cuda_exec>::test();
}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(mir_blueprint_utilities, makepointmesh_hip)
{
  AXOM_ANNOTATE_SCOPE("makepointmesh_hip");
  test_makepointmesh<hip_exec>::test();
}
#endif

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
