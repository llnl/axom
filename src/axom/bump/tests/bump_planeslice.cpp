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
  return pjoin(dataDirectory(), "bump", "regression", "bump_planeslice");
}

//------------------------------------------------------------------------------
// Global test application object.
axom::blueprint::testing::TestApplication TestApp;

template <int NDIMS>
struct make_planes;

template <>
struct make_planes<2>
{
  using PlaneType = axom::primal::Plane<double, 2>;
  using PointType = typename PlaneType::PointType;
  using VectorType = typename PlaneType::VectorType;

  static std::map<std::string, PlaneType> planes()
  {
    std::map<std::string, PlaneType> p;
    p["x=1"] = PlaneType(VectorType {{1., 0.}}, PointType {{1., 0.}});
    p["y=1"] = PlaneType(VectorType {{0., 1.}}, PointType {{0., 1.}});
    p["free"] = PlaneType(VectorType {{1., 1.}}, PointType {{1.5, 2.5}});
    return p;
  }
};

template <>
struct make_planes<3>
{
  using PlaneType = axom::primal::Plane<double, 3>;
  using PointType = typename PlaneType::PointType;
  using VectorType = typename PlaneType::VectorType;

  static std::map<std::string, PlaneType> planes()
  {
    std::map<std::string, PlaneType> p;
    p["x=0"] = PlaneType(VectorType {{1., 0., 0.}}, PointType {{0., 0., 0.}});
    p["y=0"] = PlaneType(VectorType {{0., 1., 0.}}, PointType {{0., 0., 0.}});
    p["z=0"] = PlaneType(VectorType {{0., 0., 1.}}, PointType {{0., 0., 0.}});
    p["free"] = PlaneType(VectorType {{1., 1., 1.}}, PointType {{-1.1, -1.2, -1.3}});
    return p;
  }
};

//------------------------------------------------------------------------------
template <typename ExecSpace, int NDIMS>
struct test_planeslice
{
  static void test()
  {
    const std::string name(axom::fmt::format("planeslice_{}D", NDIMS));

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

    const auto p = make_planes<NDIMS>::planes();
    for(auto it = p.begin(); it != p.end(); it++)
    {
      //_bump_utilities_planeslice_begin
      // Encode the plane in the options.
      conduit::Node hostOptions;
      hostOptions["topology"] = "mesh";
      hostOptions["normal"].set(it->second.getNormal().data(), NDIMS);
      auto origin = it->second.getNormal() * it->second.getOffset();
      hostOptions["origin"].set(origin.data(), NDIMS);

      conduit::Node deviceOptions, deviceResult;
      utils::copy<ExecSpace>(deviceOptions, hostOptions);

      axom::bump::extraction::PlaneSlice<ExecSpace, TopologyView, CoordsetView> slice(topologyView,
                                                                                      coordsetView);
      slice.execute(deviceMesh, deviceOptions, deviceResult);
      //_bump_utilities_planeslice_end

      // device->host
      conduit::Node hostResult;
      utils::copy<seq_exec>(hostResult, deviceResult);

      TestApp.saveVisualization(name + "_" + it->first, hostResult);

      // Remove some fields to make the baseline smaller
      hostResult.remove("fields/distance");

      // Handle baseline comparison.
      EXPECT_TRUE(TestApp.test<ExecSpace>(name + "_" + it->first, hostResult));
    }
  }

  static void initialize(conduit::Node &mesh)
  {
    const axom::IndexType N = 10;
    const axom::StackArray<axom::IndexType, 3> dims {N, N, (NDIMS > 2) ? N : 0};

    // Create the data
    axom::blueprint::testing::data::braid("rectilinear", dims, mesh);
  }
};

//------------------------------------------------------------------------------

TEST(bump_planeslice, planeslice_2D_seq) { test_planeslice<seq_exec, 2>::test(); }
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_planeslice, planeslice_2D_omp) { test_planeslice<omp_exec, 2>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_planeslice, planeslice_2D_cuda) { test_planeslice<cuda_exec, 2>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_planeslice, planeslice_2D_hip) { test_planeslice<hip_exec, 2>::test(); }
#endif

TEST(bump_planeslice, planeslice_3D_seq) { test_planeslice<seq_exec, 3>::test(); }
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(bump_planeslice, planeslice_3D_omp) { test_planeslice<omp_exec, 3>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(bump_planeslice, planeslice_3D_cuda) { test_planeslice<cuda_exec, 3>::test(); }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(bump_planeslice, planeslice_3D_hip) { test_planeslice<hip_exec, 3>::test(); }
#endif

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
