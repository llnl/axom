// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/mir.hpp"
#include "axom/primal.hpp"
#include "axom/mir/tests/mir_testing_data_helpers.hpp"
#include "axom/mir/tests/mir_testing_helpers.hpp"

std::string baselineDirectory()
{
  return pjoin(dataDirectory(), "mir", "regression", "mir_elvira");
}

//------------------------------------------------------------------------------
// Global test application object.
MIRTestApplication TestApp;

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct braid2d_mat_test
{
  static void initialize(const std::string &type, const std::string &mattype, conduit::Node &n_mesh)
  {
    axom::StackArray<axom::IndexType, 2> dims {10, 10};
    axom::StackArray<axom::IndexType, 2> zoneDims {dims[0] - 1, dims[1] - 1};
    axom::mir::testing::data::braid(type, dims, n_mesh);
    axom::mir::testing::data::make_matset(mattype, "mesh", zoneDims, n_mesh);
  }

  // Select a chunk of clean and mixed zones.
  static void selectZones(conduit::Node &n_options)
  {
    n_options["selectedZones"].set(std::vector<axom::IndexType> {30, 31, 32, 39, 40, 41, 48, 49, 50});
  }

  static void test(const std::string &type,
                   const std::string &mattype,
                   const std::string &name,
                   bool selectedZones = false)
  {
    namespace bputils = axom::mir::utilities::blueprint;

    // Create the data
    conduit::Node hostMesh, deviceMesh;
    initialize(type, mattype, hostMesh);
    axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);
    TestApp.saveVisualization(name + "_orig", hostMesh);

    // _elvira_mir_start
    // Make views.
    auto coordsetView =
      axom::mir::views::make_uniform_coordset<2>::view(deviceMesh["coordsets/coords"]);
    auto topologyView =
      axom::mir::views::make_uniform_topology<2>::view(deviceMesh["topologies/mesh"]);
    using CoordsetView = decltype(coordsetView);
    using TopologyView = decltype(topologyView);
    using IndexingPolicy = typename TopologyView::IndexingPolicy;

    conduit::Node deviceMIRMesh;
    if(mattype == "unibuffer")
    {
      auto matsetView =
        axom::mir::views::make_unibuffer_matset<int, float, 3>::view(deviceMesh["matsets/mat"]);
      using MatsetView = decltype(matsetView);

      using MIR = axom::mir::ElviraAlgorithm<ExecSpace, IndexingPolicy, CoordsetView, MatsetView>;
      MIR m(topologyView, coordsetView, matsetView);
      conduit::Node options;
      options["matset"] = "mat";
      options["normal"] = 1;
      if(selectedZones)
      {
        selectZones(options);
      }
      m.execute(deviceMesh, options, deviceMIRMesh);
    }
    // _elvira_mir_end

    // device->host
    conduit::Node hostMIRMesh;
    axom::mir::utilities::blueprint::copy<seq_exec>(hostMIRMesh, deviceMIRMesh);

    TestApp.saveVisualization(name, hostMIRMesh);

    // Handle baseline comparison.
    constexpr double tolerance = 2.6e-06;
    EXPECT_TRUE(TestApp.test<ExecSpace>(name, hostMIRMesh, tolerance));
  }
};

//------------------------------------------------------------------------------
TEST(mir_elvira, elvira_uniform_unibuffer_seq)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_seq");
  braid2d_mat_test<seq_exec>::test("uniform", "unibuffer", "elvira_uniform_unibuffer");
}

TEST(mir_elvira, elvira_uniform_unibuffer_sel_seq)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_sel_seq");
  constexpr bool selectZones = true;
  braid2d_mat_test<seq_exec>::test("uniform", "unibuffer", "elvira_uniform_unibuffer_sel", selectZones);
}

#if defined(AXOM_USE_OPENMP)
TEST(mir_elvira, elvira_uniform_unibuffer_omp)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_omp");
  braid2d_mat_test<omp_exec>::test("uniform", "unibuffer", "elvira_uniform_unibuffer");
}

TEST(mir_elvira, elvira_uniform_unibuffer_sel_omp)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_sel_omp");
  constexpr bool selectZones = true;
  braid2d_mat_test<omp_exec>::test("uniform", "unibuffer", "elvira_uniform_unibuffer_sel", selectZones);
}
#endif

#if defined(AXOM_USE_CUDA)
TEST(mir_elvira, elvira_uniform_unibuffer_cuda)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_cuda");
  braid2d_mat_test<cuda_exec>::test("uniform", "unibuffer", "elvira_uniform_unibuffer");
}

TEST(mir_elvira, elvira_uniform_unibuffer_sel_cuda)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_sel_cuda");
  constexpr bool selectZones = true;
  braid2d_mat_test<cuda_exec>::test("uniform", "unibuffer", "elvira_uniform_unibuffer_sel", selectZones);
}
#endif

#if defined(AXOM_USE_HIP)
TEST(mir_elvira, elvira_uniform_unibuffer_hip)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_hip");
  braid2d_mat_test<hip_exec>::test("uniform", "unibuffer", "elvira_uniform_unibuffer");
}

TEST(mir_elvira, elvira_uniform_unibuffer_sel_hip)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_sel_hip");
  constexpr bool selectZones = true;
  braid2d_mat_test<hip_exec>::test("uniform", "unibuffer", "elvira_uniform_unibuffer_sel", selectZones);
}
#endif

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
