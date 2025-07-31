// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/bump.hpp"
#include "axom/primal.hpp"
#include "axom/bump/tests/blueprint_testing_data_helpers.hpp"
#include "axom/bump/tests/blueprint_testing_helpers.hpp"

namespace bump = axom::bump;
namespace utils = axom::bump::utilities;
namespace views = axom::bump::views;

std::string baselineDirectory()
{
  return pjoin(dataDirectory(), "bump", "regression", "bump_make_polyhedral_topology");
}

//------------------------------------------------------------------------------
// Global test application object.
axom::blueprint::testing::TestApplication TestApp;

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct make_polyhedral
{
  static void initialize(const std::string &type, conduit::Node &n_mesh)
  {
    axom::StackArray<axom::IndexType, 3> dims {4, 4, 4};
    axom::blueprint::testing::data::braid(type, dims, n_mesh);
  }

  static void test(const std::string &type, const std::string &name)
  {
    // Create the data
    conduit::Node hostMesh, deviceMesh;
    initialize(type, hostMesh);
    utils::copy<ExecSpace>(deviceMesh, hostMesh);
    TestApp.saveVisualization(name + "_orig", hostMesh);

    //_bump_utilities_makepolyhedraltopology_begin
    // Run the algorithm
    const conduit::Node &n_input = deviceMesh["topologies/mesh"];
    conduit::Node &n_output = deviceMesh["topologies/polymesh"];
    if(type == "uniform")
    {
      auto topologyView = views::make_uniform_topology<3>::view(n_input);
      using TopologyView = decltype(topologyView);
      using ConnectivityType = typename TopologyView::ConnectivityType;

      bump::MakePolyhedralTopology<ExecSpace, TopologyView> mp(topologyView);
      mp.execute(n_input, n_output);
      bump::MergePolyhedralFaces<ExecSpace, ConnectivityType>::execute(n_output);
    }
    //_bump_utilities_makepolyhedraltopology_end
    else if(type == "tets")
    {
      using TetShape = views::TetShape<int>;
      auto topologyView = views::make_unstructured_single_shape_topology<TetShape>::view(n_input);
      using TopologyView = decltype(topologyView);
      using ConnectivityType = typename TopologyView::ConnectivityType;

      bump::MakePolyhedralTopology<ExecSpace, TopologyView> mp(topologyView);
      mp.execute(n_input, n_output);
      bump::MergePolyhedralFaces<ExecSpace, ConnectivityType>::execute(n_output);
    }
    else if(type == "pyramids")
    {
      using PyramidShape = views::PyramidShape<int>;
      auto topologyView = views::make_unstructured_single_shape_topology<PyramidShape>::view(n_input);
      using TopologyView = decltype(topologyView);
      using ConnectivityType = typename TopologyView::ConnectivityType;

      bump::MakePolyhedralTopology<ExecSpace, TopologyView> mp(topologyView);
      mp.execute(n_input, n_output);
      bump::MergePolyhedralFaces<ExecSpace, ConnectivityType>::execute(n_output);
    }
    else if(type == "wedges")
    {
      using WedgeShape = views::WedgeShape<int>;
      auto topologyView = views::make_unstructured_single_shape_topology<WedgeShape>::view(n_input);
      using TopologyView = decltype(topologyView);
      using ConnectivityType = typename TopologyView::ConnectivityType;

      bump::MakePolyhedralTopology<ExecSpace, TopologyView> mp(topologyView);
      mp.execute(n_input, n_output);
      bump::MergePolyhedralFaces<ExecSpace, ConnectivityType>::execute(n_output);
    }
    else if(type == "hexs")
    {
      using HexShape = views::HexShape<int>;
      auto topologyView = views::make_unstructured_single_shape_topology<HexShape>::view(n_input);
      using TopologyView = decltype(topologyView);
      using ConnectivityType = typename TopologyView::ConnectivityType;

      bump::MakePolyhedralTopology<ExecSpace, TopologyView> mp(topologyView);
      mp.execute(n_input, n_output);
      bump::MergePolyhedralFaces<ExecSpace, ConnectivityType>::execute(n_output);
    }
    else
    {
      SLIC_ERROR(axom::fmt::format("Unsupported shape in the test: {}", type));
    }

    // device->host
    conduit::Node hostOutputMesh;
    utils::copy<seq_exec>(hostOutputMesh, deviceMesh);

    TestApp.saveVisualization(name, hostOutputMesh);

    // Handle baseline comparison.
    constexpr double tolerance = 2.6e-06;
    EXPECT_TRUE(TestApp.test<ExecSpace>(name, hostOutputMesh, tolerance));
  }
};

//------------------------------------------------------------------------------
TEST(bump_make_polyhedral_topology, uniform_seq)
{
  AXOM_ANNOTATE_SCOPE("uniform_seq");
  make_polyhedral<seq_exec>::test("uniform", "make_polyhedral_uniform");
}

TEST(bump_make_polyhedral_topology, tets_seq)
{
  AXOM_ANNOTATE_SCOPE("tets_seq");
  make_polyhedral<seq_exec>::test("tets", "make_polyhedral_tets");
}

TEST(bump_make_polyhedral_topology, pyramids_seq)
{
  AXOM_ANNOTATE_SCOPE("pyramids_seq");
  make_polyhedral<seq_exec>::test("pyramids", "make_polyhedral_pyramids");
}

TEST(bump_make_polyhedral_topology, wedges_seq)
{
  AXOM_ANNOTATE_SCOPE("wedges_seq");
  make_polyhedral<seq_exec>::test("wedges", "make_polyhedral_wedges");
}

TEST(bump_make_polyhedral_topology, hexs_seq)
{
  AXOM_ANNOTATE_SCOPE("hexs_seq");
  make_polyhedral<seq_exec>::test("hexs", "make_polyhedral_hexs");
}

#if defined(AXOM_USE_OPENMP)
TEST(bump_make_polyhedral_topology, uniform_omp)
{
  AXOM_ANNOTATE_SCOPE("uniform_omp");
  make_polyhedral<omp_exec>::test("uniform", "make_polyhedral_uniform");
}

TEST(bump_make_polyhedral_topology, tets_omp)
{
  AXOM_ANNOTATE_SCOPE("tets_omp");
  make_polyhedral<omp_exec>::test("tets", "make_polyhedral_tets");
}

TEST(bump_make_polyhedral_topology, pyramids_omp)
{
  AXOM_ANNOTATE_SCOPE("pyramids_omp");
  make_polyhedral<omp_exec>::test("pyramids", "make_polyhedral_pyramids");
}

TEST(bump_make_polyhedral_topology, wedges_omp)
{
  AXOM_ANNOTATE_SCOPE("wedges_omp");
  make_polyhedral<omp_exec>::test("wedges", "make_polyhedral_wedges");
}

TEST(bump_make_polyhedral_topology, hexs_omp)
{
  AXOM_ANNOTATE_SCOPE("hexs_omp");
  make_polyhedral<omp_exec>::test("hexs", "make_polyhedral_hexs");
}
#endif

#if defined(AXOM_USE_CUDA)
TEST(bump_make_polyhedral_topology, uniform_cuda)
{
  AXOM_ANNOTATE_SCOPE("uniform_cuda");
  make_polyhedral<cuda_exec>::test("uniform", "make_polyhedral_uniform");
}

TEST(bump_make_polyhedral_topology, tets_cuda)
{
  AXOM_ANNOTATE_SCOPE("tets_cuda");
  make_polyhedral<cuda_exec>::test("tets", "make_polyhedral_tets");
}

TEST(bump_make_polyhedral_topology, pyramids_cuda)
{
  AXOM_ANNOTATE_SCOPE("pyramids_cuda");
  make_polyhedral<cuda_exec>::test("pyramids", "make_polyhedral_pyramids");
}

TEST(bump_make_polyhedral_topology, wedges_cuda)
{
  AXOM_ANNOTATE_SCOPE("wedges_cuda");
  make_polyhedral<cuda_exec>::test("wedges", "make_polyhedral_wedges");
}

TEST(bump_make_polyhedral_topology, hexs_cuda)
{
  AXOM_ANNOTATE_SCOPE("hexs_cuda");
  make_polyhedral<cuda_exec>::test("hexs", "make_polyhedral_hexs");
}
#endif

#if defined(AXOM_USE_HIP)
TEST(bump_make_polyhedral_topology, uniform_hip)
{
  AXOM_ANNOTATE_SCOPE("uniform_hip");
  make_polyhedral<hip_exec>::test("uniform", "make_polyhedral_uniform");
}

TEST(bump_make_polyhedral_topology, tets_hip)
{
  AXOM_ANNOTATE_SCOPE("tets_hip");
  make_polyhedral<hip_exec>::test("tets", "make_polyhedral_tets");
}

TEST(bump_make_polyhedral_topology, pyramids_hip)
{
  AXOM_ANNOTATE_SCOPE("pyramids_hip");
  make_polyhedral<hip_exec>::test("pyramids", "make_polyhedral_pyramids");
}

TEST(bump_make_polyhedral_topology, wedges_hip)
{
  AXOM_ANNOTATE_SCOPE("wedges_hip");
  make_polyhedral<hip_exec>::test("wedges", "make_polyhedral_wedges");
}

TEST(bump_make_polyhedral_topology, hexs_hip)
{
  AXOM_ANNOTATE_SCOPE("hexs_hip");
  make_polyhedral<hip_exec>::test("hexs", "make_polyhedral_hexs");
}
#endif

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
