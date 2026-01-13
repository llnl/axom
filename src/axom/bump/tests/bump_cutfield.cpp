// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/bump.hpp"
#include "axom/bump/tests/blueprint_testing_data_helpers.hpp"
#include "axom/bump/tests/blueprint_testing_helpers.hpp"

#include <conduit/conduit_blueprint_mesh_examples_tiled.hpp>
#include <conduit/conduit_relay_io_blueprint.hpp>
#include <cmath>
#include <cstdlib>

namespace bump = axom::bump;
namespace views = axom::bump::views;
namespace utils = axom::bump::utilities;

std::string baselineDirectory()
{
  return pjoin(dataDirectory(), "bump", "regression", "bump_cutfield");
}

//------------------------------------------------------------------------------
// Global test application object.
axom::blueprint::testing::TestApplication TestApp;

//------------------------------------------------------------------------------

template <typename ExecSpace, int NDIMS>
struct test_cutfield
{
  static void test()
  {
    const std::string name(axom::fmt::format("cutfield_{}D", NDIMS));

    conduit::Node hostMesh;
    initialize(hostMesh);

    TestApp.saveVisualization(name + "_orig", hostMesh);

    // host->device
    conduit::Node deviceMesh;
    utils::copy<ExecSpace>(deviceMesh, hostMesh);

    // _bump_utilities_cutfield_begin
    // Wrap the data in views.
    auto coordsetView = axom::bump::views::make_rectilinear_coordset<conduit::float64, NDIMS>::view(
      deviceMesh["coordsets/coords"]);
    using CoordsetView = decltype(coordsetView);

    auto topologyView =
      axom::bump::views::make_rectilinear_topology<NDIMS>::view(deviceMesh["topologies/mesh"]);
    using TopologyView = decltype(topologyView);

    conduit::Node hostOptions;
    hostOptions["field"] = "braid";
    hostOptions["value"] = 1.;

    conduit::Node deviceOptions, deviceResult;
    utils::copy<ExecSpace>(deviceOptions, hostOptions);
    axom::bump::extraction::CutField<ExecSpace, TopologyView, CoordsetView> iso(topologyView,
                                                                                coordsetView);
    iso.execute(deviceMesh, deviceOptions, deviceResult);
    // _bump_utilities_cutfield_end

    // device->host
    conduit::Node hostResult;
    utils::copy<seq_exec>(hostResult, deviceResult);

    TestApp.saveVisualization(name + "_braid", hostResult);

    // Remove some fields to make the baseline smaller
    hostResult.remove("fields/vel");
    hostResult.remove("fields/distance");
    hostResult.remove("fields/braid");

    // Handle baseline comparison.
    EXPECT_TRUE(TestApp.test<ExecSpace>(name + "_braid", hostResult));

    //---------------------
    // Try a different clip
    hostOptions["field"] = "gyroid";
    hostOptions["value"] = 0;
    utils::copy<ExecSpace>(deviceOptions, hostOptions);
    hostResult.reset();
    deviceResult.reset();
    iso.execute(deviceMesh, deviceOptions, deviceResult);

    // device->host
    utils::copy<seq_exec>(hostResult, deviceResult);

    // Save a field for the number of sides in each shape.
    hostResult["fields/sides/topology"] = "mesh";
    hostResult["fields/sides/association"] = "element";
    hostResult["fields/sides/values"].set_external(hostResult["topologies/mesh/elements/sizes"]);

    TestApp.saveVisualization(name + "_gyroid", hostResult);

    // Remove some fields to make the baseline smaller
    hostResult.remove("fields/vel");
    hostResult.remove("fields/distance");
    hostResult.remove("fields/gyroid");
    hostResult.remove("fields/sides");

    EXPECT_TRUE(TestApp.test<ExecSpace>(name + "_gyroid", hostResult));
  }

  static void initialize(conduit::Node &mesh)
  {
    const axom::IndexType N = 20;
    const axom::StackArray<axom::IndexType, 3> dims {N, N, (NDIMS > 2) ? N : 0};
    const auto maxZ = axom::utilities::max(dims[2], axom::IndexType {1});
    const auto nnodes = dims[0] * dims[1] * maxZ;

    // Create the data
    axom::blueprint::testing::data::braid("rectilinear", dims, mesh);

    // Add a gyroid field
    mesh["fields/gyroid/topology"] = "mesh";
    mesh["fields/gyroid/association"] = "vertex";
    mesh["fields/gyroid/values"].set(conduit::DataType::float64(nnodes));
    auto gyroid = mesh["fields/gyroid/values"].as_float64_ptr();

    auto coordsetView = axom::bump::views::make_rectilinear_coordset<conduit::float64, NDIMS>::view(
      mesh["coordsets/coords"]);
    for(axom::IndexType k = 0; k < maxZ; k++)
    {
      const auto kOffset = k * dims[0] * dims[1];
      for(axom::IndexType j = 0; j < dims[1]; j++)
      {
        for(axom::IndexType i = 0; i < dims[0]; i++)
        {
          const auto index = kOffset + j * dims[0] + i;
          const auto pt = coordsetView[index];
          const double scale = 0.5;
          const auto x = scale * pt[0];
          const auto y = scale * pt[1];
          const auto z = scale * ((NDIMS == 2) ? 0. : pt[2]);
          gyroid[index] = sin(x) * cos(y) + sin(y) * cos(z) + sin(z) * cos(x);
        }
      }
    }
  }

  static void test_polygonal()
  {
    const std::string name("cutfield_2D_polygonal");

    conduit::Node hostMesh;
    initialize_polygonal(hostMesh);

    TestApp.saveVisualization(name + "_orig", hostMesh);

    // host->device
    conduit::Node deviceMesh;
    utils::copy<ExecSpace>(deviceMesh, hostMesh);

    // Wrap the data in views.
    auto coordsetView = axom::bump::views::make_explicit_coordset<conduit::float64, NDIMS>::view(
      deviceMesh["coordsets/coords"]);
    using CoordsetView = decltype(coordsetView);

    auto topologyView = axom::bump::views::make_unstructured_single_shape_topology<
      views::PolygonShape<conduit::index_t>>::view(deviceMesh["topologies/mesh"]);
    using TopologyView = decltype(topologyView);

    conduit::Node hostOptions, deviceOptions, deviceResult;
    hostOptions["field"] = "gyroid";
    hostOptions["value"] = 0;
    utils::copy<ExecSpace>(deviceOptions, hostOptions);
    axom::bump::extraction::CutField<ExecSpace, TopologyView, CoordsetView> iso(topologyView,
                                                                                coordsetView);
    iso.execute(deviceMesh, deviceOptions, deviceResult);

    // device->host
    conduit::Node hostResult;
    utils::copy<seq_exec>(hostResult, deviceResult);

    TestApp.saveVisualization(name + "_gyroid", hostResult);

    EXPECT_TRUE(TestApp.test<ExecSpace>(name + "_gyroid", hostResult));
  }

  static void initialize_polygonal(conduit::Node &n_mesh)
  {
    // This is a tile definition for a tile that contains polygons with 3-8 sides.
    static const char *tile = R"(
coordsets:
  coords:
    type: explicit
    values:
      x: [0., 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1.,      # row 0: points 0-9
          0.25, 0.35, 0.65, 0.75,                              # row 1: points 10-13
          0., 1.,                                              # row 2: points 14,15
          0.25, 0.35, 0.65, 0.75,                              # row 3: points 16-19
          0.1, 0.2, 0.4, 0.6, 0.8, 0.9,                        # row 4: points 20-25 (middle)
          0.25, 0.35, 0.65, 0.75,                              # row 5: points 26-29
          0., 1.,                                              # row 6: points 30-31
          0.25, 0.35, 0.65, 0.75,                              # row 7: points 32-35
          0., 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1.       # row 8: points 36-45
         ]
      y: [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,              # row 0: points 0-9
          0.125, 0.125, 0.125, 0.125,                          # row 1: points 10-13
          0.25, 0.25,                                          # row 2: points 14,15
          0.375, 0.375, 0.375, 0.375,                          # row 3: points 16-19
          0.5, 0.5, 0.5, 0.5, 0.5, 0.5,                        # row 4: points 20-25 (middle)
          0.625, 0.625, 0.625, 0.625,                          # row 5: points 26-29
          0.75, 0.75,                                          # row 6: points 30-31
          0.875, 0.875, 0.875, 0.875,                          # row 7: points 32-35
          1., 1., 1., 1., 1., 1., 1., 1., 1., 1.               # row 8: points 36-45
         ]
topologies:
  mesh:
    type: unstructured
    coordset: coords
    elements:
      shape: polygonal
      connectivity: [0, 1, 14,
                     1,2,10,16,21,20,14,
                     2, 3, 10,
                     3,11,17,16,10,
                     3,4,11,
                     4,5,12,18,23,22,17,11,
                     5,6,12,
                     6,13,19,18,12,
                     6,7,13,
                     7,8,15,25,24,19,13,
                     8,9,15,
                     14,20,30,
                     16,17,22,27,26,21,
                     18,19,24,29,28,23,
                     15,31,25,
                     30,37,36,
                     20,21,26,30,
                     30,26,32,37,
                     32,38,37,
                     32,39,38,
                     26,27,33,39,32,
                     33,40,39,
                     22,23,28,34,41,40,33,27,
                     34,42,41,
                     28,29,35,42,34,
                     35,43,42,
                     24,25,31,44,43,35,29,
                     31,45,44
                    ]
      sizes: [3, 7, 3, 5, 3, 8, 3, 5, 3, 7, 3, 3, 6, 6, 3, 3, 4, 4, 3, 3, 5, 3, 8, 3, 5, 3, 7, 3]
      offsets: [0, 3, 10, 13, 18, 21, 29, 32, 37, 40, 47, 50, 53, 59, 65, 68, 71, 75, 79, 82, 85, 90, 93, 101, 104, 109, 112, 119]
left: [0,14,30,36]
right: [9,15,31,45]
bottom: [0,1,2,3,4,5,6,7,8,9]
top: [36,37,38,39,40,41,42,43,44,45]
translate:
  x: 1.
  y: 1.
    )";

    conduit::Node n_options;
    n_options["tile"].parse(tile);
    const int N = 10;
    conduit::blueprint::mesh::examples::tiled(N, N, 0, n_mesh, n_options);
    conduit::blueprint::mesh::topology::unstructured::generate_offsets_inline(
      n_mesh["topologies/mesh"]);

    const auto xc = n_mesh["coordsets/coords/values/x"].as_double_accessor();
    const auto yc = n_mesh["coordsets/coords/values/y"].as_double_accessor();
    const auto nnodes = xc.number_of_elements();

    // Add a gyroid field
    n_mesh["fields/gyroid/topology"] = "mesh";
    n_mesh["fields/gyroid/association"] = "vertex";
    n_mesh["fields/gyroid/values"].set(conduit::DataType::float64(nnodes));
    auto gyroid = n_mesh["fields/gyroid/values"].as_float64_ptr();
    for(axom::IndexType i = 0; i < nnodes; i++)
    {
      const double scale = 2.;
      const double x = scale * xc[i];
      const double y = scale * yc[i];
      const double z = 0.;
      gyroid[i] = sin(x) * cos(y) + sin(y) * cos(z) + sin(z) * cos(x);
    }
  }
};

TEST(bump_cutfield, cutfield_2D_seq) { test_cutfield<seq_exec, 2>::test(); }
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_cutfield, cutfield_2D_omp) { test_cutfield<omp_exec, 2>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_cutfield, cutfield_2D_cuda) { test_cutfield<cuda_exec, 2>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_cutfield, cutfield_2D_hip) { test_cutfield<hip_exec, 2>::test(); }
#endif

TEST(bump_cutfield, cutfield_2D_polygonal_seq) { test_cutfield<seq_exec, 2>::test_polygonal(); }
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_cutfield, cutfield_2D_polygonal_omp) { test_cutfield<omp_exec, 2>::test_polygonal(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_cutfield, cutfield_2D_polygonal_cuda) { test_cutfield<cuda_exec, 2>::test_polygonal(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_cutfield, cutfield_2D_polygonal_hip) { test_cutfield<hip_exec, 2>::test_polygonal(); }
#endif

TEST(bump_cutfield, cutfield_3D_seq) { test_cutfield<seq_exec, 3>::test(); }
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_cutfield, cutfield_3D_omp) { test_cutfield<omp_exec, 3>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_cutfield, cutfield_3D_cuda) { test_cutfield<cuda_exec, 3>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_cutfield, cutfield_3D_hip) { test_cutfield<hip_exec, 3>::test(); }
#endif

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
