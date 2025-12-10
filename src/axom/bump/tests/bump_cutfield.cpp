// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/bump.hpp"
#include "axom/bump/tests/blueprint_testing_data_helpers.hpp"
#include "axom/bump/tests/blueprint_testing_helpers.hpp"

#include "axom/bump/MapBasedNaming.hpp"

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
#if 0
    // For Debugging
    axom::bump::extraction::CutField<ExecSpace,
                                     TopologyView,
                                     CoordsetView,
                                     axom::bump::extraction::FieldIntersector<ExecSpace, TopologyView, CoordsetView>,
                                     axom::bump::MapBasedNaming<axom::IndexType>> iso(topologyView, coordsetView);
#else
    axom::bump::extraction::CutField<ExecSpace, TopologyView, CoordsetView> iso(topologyView,
                                                                                coordsetView);
#endif
    iso.execute(deviceMesh, deviceOptions, deviceResult);

    // device->host
    conduit::Node hostResult;
    utils::copy<seq_exec>(hostResult, deviceResult);

    TestApp.saveVisualization(name + "_braid", hostResult);

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

    TestApp.saveVisualization(name + "_gyroid", hostResult);

    EXPECT_TRUE(TestApp.test<ExecSpace>(name + "_gyroid", hostResult));
  }

  static void initialize(conduit::Node &mesh)
  {
    const axom::IndexType N = 20;
    const axom::StackArray<axom::IndexType, 3> dims {N, N, (NDIMS > 2) ? N : 0};
    const auto maxZ = axom::utilities::max(dims[2], axom::IndexType{1});
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
          gyroid[index] = sin(x) * cos(y) + 
                          sin(y) * cos(z) + 
                          sin(z) * cos(x);
        }
      }
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
