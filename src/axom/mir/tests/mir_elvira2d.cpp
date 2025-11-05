// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/bump.hpp"
#include "axom/mir.hpp"
#include "axom/primal.hpp"
#include "axom/bump/tests/blueprint_testing_data_helpers.hpp"
#include "axom/bump/tests/blueprint_testing_helpers.hpp"

#include <cstdlib>

std::string baselineDirectory()
{
  return pjoin(dataDirectory(), "mir", "regression", "mir_elvira");
}

namespace utils = axom::bump::utilities;

//------------------------------------------------------------------------------
// Global test application object.
axom::blueprint::testing::TestApplication TestApp;

//------------------------------------------------------------------------------
template <typename ExecSpace>
constexpr int maxAttempts()
{
  int n = 1;
#if defined(AXOM_USE_HIP)
  if constexpr(axom::execution_space<ExecSpace>::onDevice())
  {
    n = 2;
  }
#endif
  return n;
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct braid2d_mat_test
{
  static void initialize(const std::string &type,
                         const std::string &mattype,
                         bool cleanMats,
                         conduit::Node &n_mesh)
  {
    axom::StackArray<axom::IndexType, 2> dims {10, 10};
    axom::StackArray<axom::IndexType, 2> zoneDims {dims[0] - 1, dims[1] - 1};
    axom::blueprint::testing::data::braid(type, dims, n_mesh);
    axom::blueprint::testing::data::make_matset(mattype, "mesh", zoneDims, cleanMats, n_mesh);
  }

  // Select a chunk of clean and mixed zones.
  static void selectZones(conduit::Node &n_options)
  {
    n_options["selectedZones"].set(std::vector<axom::IndexType> {30, 31, 32, 39, 40, 41, 48, 49, 50});
  }

  static void test(const std::string &type,
                   const std::string &mattype,
                   const std::string &name,
                   bool selectedZones = false,
                   bool pointMesh = false,
                   bool cleanMats = false,
                   int nDomains = 1)
  {
    // Create the data (1+ domains)
    conduit::Node hostMesh, deviceMesh;
    for(int dom = 0; dom < nDomains; dom++)
    {
      const std::string domainName = axom::fmt::format("domain_{:07}", dom);
      conduit::Node &hostDomain = (nDomains > 1) ? hostMesh[domainName] : hostMesh;

      initialize(type, mattype, cleanMats, hostDomain);
      TestApp.saveVisualization(name + "_orig", hostDomain);
    }
    utils::copy<ExecSpace>(deviceMesh, hostMesh);

    // NOTE: As a workaround on HIP, we can attempt the test more than once.
    bool pass = false;
    for(int attempt = 0; attempt < maxAttempts<ExecSpace>(); attempt++)
    {
      // Do MIR on all domains.
      bool thisAttempt = true;
      for(int dom = 0; dom < nDomains; dom++)
      {
        const std::string domainName = axom::fmt::format("domain_{:07}", dom);
        conduit::Node &deviceDomain = (nDomains > 1) ? deviceMesh[domainName] : deviceMesh;

        // _elvira_mir_start
        namespace views = axom::bump::views;
        // Make views.
        auto coordsetView = views::make_uniform_coordset<2>::view(deviceDomain["coordsets/coords"]);
        auto topologyView = views::make_uniform_topology<2>::view(deviceDomain["topologies/mesh"]);
        using CoordsetView = decltype(coordsetView);
        using TopologyView = decltype(topologyView);
        using IndexingPolicy = typename TopologyView::IndexingPolicy;

        conduit::Node deviceMIRMesh;
        if(mattype == "unibuffer")
        {
          auto matsetView =
            views::make_unibuffer_matset<int, float, 3>::view(deviceDomain["matsets/mat"]);
          using MatsetView = decltype(matsetView);

          using MIR = axom::mir::ElviraAlgorithm<ExecSpace, IndexingPolicy, CoordsetView, MatsetView>;
          MIR m(topologyView, coordsetView, matsetView);
          conduit::Node options;
          options["matset"] = "mat";
          options["plane"] = 1;
          options["pointmesh"] = pointMesh ? 1 : 0;
          if(cleanMats)
          {
            // Set the output names
            options["topologyName"] = "postmir_topology";
            options["coordsetName"] = "postmir_coords";
            options["matsetName"] = "postmir_matset";
          }
          if(selectedZones)
          {
            selectZones(options);
          }
          m.execute(deviceDomain, options, deviceMIRMesh);
        }
        // _elvira_mir_end

        // device->host
        conduit::Node hostMIRMesh;
        utils::copy<seq_exec>(hostMIRMesh, deviceMIRMesh);

        // Verify the hostMIRMesh to look for errors.
        conduit::Node info;
        bool verifyOK = conduit::blueprint::mesh::verify(hostMIRMesh, info);
        if(!verifyOK)
        {
          printNode(hostMIRMesh);
          info.print();
        }
        EXPECT_TRUE(verifyOK);

        TestApp.saveVisualization(name, hostMIRMesh);

        // Handle baseline comparison.
        constexpr double tolerance = 2.6e-06;
        TestApp.setVerbose(attempt >= maxAttempts<ExecSpace>() - 1);
        thisAttempt &= TestApp.test<ExecSpace>(name, hostMIRMesh, tolerance);
      }

      // See whether this attempt worked.
      pass = thisAttempt;
      if(thisAttempt)
      {
        break;
      }
      else if(attempt + 1 < maxAttempts<ExecSpace>())
      {
        SLIC_WARNING(axom::fmt::format("Attempting test {} again!", name));
      }
    }
    EXPECT_TRUE(pass);
    TestApp.setVerbose(true);

    reset();
  }

  /*
   * \brief This function runs a simple kernel and synchronizes, both of which seem
   *        to be needed to work around intermittent failures on HIP platforms.
   *        Users can set the "NO_RESET" environment variable to skip this workaround.
   */
  static void reset()
  {
    if(getenv("NO_RESET") == nullptr)
    {
      const axom::IndexType N = 10000;
      axom::Array<double> arr(N, N, axom::execution_space<ExecSpace>::allocatorID());
      auto arrView = arr.view();
      for(int i = 0; i < 2; i++)
      {
        axom::for_all<ExecSpace>(
          N,
          AXOM_LAMBDA(axom::IndexType index) { arrView[index] = index * index; });
      }
      axom::synchronize<ExecSpace>();
    }
  }
};

//------------------------------------------------------------------------------
TEST(mir_elvira, options)
{
  conduit::Node n_options;

  axom::mir::ELVIRAOptions opts(n_options);
  EXPECT_FALSE(opts.pointmesh());
  EXPECT_FALSE(opts.plane());

  n_options["plane"] = 1;
  n_options["pointmesh"] = 1;

  EXPECT_TRUE(opts.pointmesh());
  EXPECT_TRUE(opts.plane());
}

//------------------------------------------------------------------------------
TEST(mir_elvira, elvira_uniform_unibuffer_seq)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_seq");
  const bool selectZones = false;
  const bool pointMesh = false;
  braid2d_mat_test<seq_exec>::test("uniform",
                                   "unibuffer",
                                   "elvira_uniform_unibuffer",
                                   selectZones,
                                   pointMesh);
  // Run 2 domain example
  {
    const int nDomains = 2;
    const bool cleanMats = false;
    braid2d_mat_test<seq_exec>::test("uniform",
                                     "unibuffer",
                                     "elvira_uniform_unibuffer",
                                     selectZones,
                                     pointMesh,
                                     cleanMats,
                                     nDomains);
  }

  // Run clean mats example.
  {
    const bool cleanMats = true;
    braid2d_mat_test<seq_exec>::test("uniform",
                                     "unibuffer",
                                     "elvira_uniform_unibuffer_clean",
                                     selectZones,
                                     pointMesh,
                                     cleanMats);
  }
}

TEST(mir_elvira, elvira_uniform_unibuffer_sel_seq)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_sel_seq");
  const bool selectZones = true;
  const bool pointMesh = false;
  braid2d_mat_test<seq_exec>::test("uniform",
                                   "unibuffer",
                                   "elvira_uniform_unibuffer_sel",
                                   selectZones,
                                   pointMesh);

  // Run clean mats example with selected zones.
  {
    const bool cleanMats = true;
    braid2d_mat_test<seq_exec>::test("uniform",
                                     "unibuffer",
                                     "elvira_uniform_unibuffer_sel_clean",
                                     selectZones,
                                     pointMesh,
                                     cleanMats);
  }
}

TEST(mir_elvira, elvira_uniform_unibuffer_seq_pm)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_pm_seq");
  const bool selectZones = false;
  const bool pointMesh = true;
  braid2d_mat_test<seq_exec>::test("uniform",
                                   "unibuffer",
                                   "elvira_uniform_unibuffer_pm",
                                   selectZones,
                                   pointMesh);
}

TEST(mir_elvira, elvira_uniform_unibuffer_sel_pm_seq)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_sel_pm_seq");
  const bool selectZones = true;
  const bool pointMesh = true;
  braid2d_mat_test<seq_exec>::test("uniform",
                                   "unibuffer",
                                   "elvira_uniform_unibuffer_sel_pm",
                                   selectZones,
                                   pointMesh);
}

#if defined(AXOM_USE_OPENMP)
TEST(mir_elvira, elvira_uniform_unibuffer_omp)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_omp");
  const bool selectZones = false;
  const bool pointMesh = false;
  braid2d_mat_test<omp_exec>::test("uniform",
                                   "unibuffer",
                                   "elvira_uniform_unibuffer",
                                   selectZones,
                                   pointMesh);
  // Run 2 domain example
  {
    const int nDomains = 2;
    const bool cleanMats = false;
    braid2d_mat_test<omp_exec>::test("uniform",
                                     "unibuffer",
                                     "elvira_uniform_unibuffer",
                                     selectZones,
                                     pointMesh,
                                     cleanMats,
                                     nDomains);
  }
  // Run clean mats example.
  {
    const bool cleanMats = true;
    braid2d_mat_test<seq_exec>::test("uniform",
                                     "unibuffer",
                                     "elvira_uniform_unibuffer_clean",
                                     selectZones,
                                     pointMesh,
                                     cleanMats);
  }
}

TEST(mir_elvira, elvira_uniform_unibuffer_sel_omp)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_sel_omp");
  const bool selectZones = true;
  const bool pointMesh = false;
  braid2d_mat_test<omp_exec>::test("uniform",
                                   "unibuffer",
                                   "elvira_uniform_unibuffer_sel",
                                   selectZones,
                                   pointMesh);
  // Run clean mats example with selected zones.
  {
    const bool cleanMats = true;
    braid2d_mat_test<omp_exec>::test("uniform",
                                     "unibuffer",
                                     "elvira_uniform_unibuffer_sel_clean",
                                     selectZones,
                                     pointMesh,
                                     cleanMats);
  }
}

TEST(mir_elvira, elvira_uniform_unibuffer_pm_omp)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_pm_omp");
  const bool selectZones = false;
  const bool pointMesh = true;
  braid2d_mat_test<omp_exec>::test("uniform",
                                   "unibuffer",
                                   "elvira_uniform_unibuffer_pm",
                                   selectZones,
                                   pointMesh);
}

TEST(mir_elvira, elvira_uniform_unibuffer_sel_pm_omp)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_sel_pm_omp");
  const bool selectZones = true;
  const bool pointMesh = true;
  braid2d_mat_test<omp_exec>::test("uniform",
                                   "unibuffer",
                                   "elvira_uniform_unibuffer_sel_pm",
                                   selectZones,
                                   pointMesh);
}
#endif

#if defined(AXOM_USE_CUDA)
TEST(mir_elvira, elvira_uniform_unibuffer_cuda)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_cuda");
  const bool selectZones = false;
  const bool pointMesh = false;
  braid2d_mat_test<cuda_exec>::test("uniform",
                                    "unibuffer",
                                    "elvira_uniform_unibuffer",
                                    selectZones,
                                    pointMesh);
  // Run 2 domain example
  {
    const int nDomains = 2;
    const bool cleanMats = false;
    braid2d_mat_test<cuda_exec>::test("uniform",
                                      "unibuffer",
                                      "elvira_uniform_unibuffer",
                                      selectZones,
                                      pointMesh,
                                      cleanMats,
                                      nDomains);
  }
  // Run clean mats example.
  {
    const bool cleanMats = true;
    braid2d_mat_test<seq_exec>::test("uniform",
                                     "unibuffer",
                                     "elvira_uniform_unibuffer_clean",
                                     selectZones,
                                     pointMesh,
                                     cleanMats);
  }
}

TEST(mir_elvira, elvira_uniform_unibuffer_sel_cuda)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_sel_cuda");
  const bool selectZones = true;
  const bool pointMesh = false;
  braid2d_mat_test<cuda_exec>::test("uniform",
                                    "unibuffer",
                                    "elvira_uniform_unibuffer_sel",
                                    selectZones,
                                    pointMesh);

  // Run clean mats example with selected zones.
  {
    const bool cleanMats = true;
    braid2d_mat_test<cuda_exec>::test("uniform",
                                      "unibuffer",
                                      "elvira_uniform_unibuffer_sel_clean",
                                      selectZones,
                                      pointMesh,
                                      cleanMats);
  }
}

TEST(mir_elvira, elvira_uniform_unibuffer_pm_cuda)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_pm_cuda");
  const bool selectZones = false;
  const bool pointMesh = true;
  braid2d_mat_test<cuda_exec>::test("uniform",
                                    "unibuffer",
                                    "elvira_uniform_unibuffer_pm",
                                    selectZones,
                                    pointMesh);
}

TEST(mir_elvira, elvira_uniform_unibuffer_sel_pm_cuda)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_sel_pm_cuda");
  const bool selectZones = true;
  const bool pointMesh = true;
  braid2d_mat_test<cuda_exec>::test("uniform",
                                    "unibuffer",
                                    "elvira_uniform_unibuffer_sel_pm",
                                    selectZones,
                                    pointMesh);
}
#endif

#if defined(AXOM_USE_HIP)
TEST(mir_elvira, elvira_uniform_unibuffer_hip)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_hip");
  const bool selectZones = false;
  const bool pointMesh = false;
  braid2d_mat_test<hip_exec>::test("uniform",
                                   "unibuffer",
                                   "elvira_uniform_unibuffer",
                                   selectZones,
                                   pointMesh);
  // Run 2 domain example
  {
    const int nDomains = 2;
    const bool cleanMats = false;
    braid2d_mat_test<hip_exec>::test("uniform",
                                     "unibuffer",
                                     "elvira_uniform_unibuffer",
                                     selectZones,
                                     pointMesh,
                                     cleanMats,
                                     nDomains);
  }
  // Run clean mats example.
  {
    const bool cleanMats = true;
    braid2d_mat_test<seq_exec>::test("uniform",
                                     "unibuffer",
                                     "elvira_uniform_unibuffer_clean",
                                     selectZones,
                                     pointMesh,
                                     cleanMats);
  }
}

TEST(mir_elvira, elvira_uniform_unibuffer_sel_hip)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_sel_hip");
  const bool selectZones = true;
  const bool pointMesh = false;
  braid2d_mat_test<hip_exec>::test("uniform",
                                   "unibuffer",
                                   "elvira_uniform_unibuffer_sel",
                                   selectZones,
                                   pointMesh);
  // Run clean mats example with selected zones.
  {
    const bool cleanMats = true;
    braid2d_mat_test<hip_exec>::test("uniform",
                                     "unibuffer",
                                     "elvira_uniform_unibuffer_sel_clean",
                                     selectZones,
                                     pointMesh,
                                     cleanMats);
  }
}

TEST(mir_elvira, elvira_uniform_unibuffer_pm_hip)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_pm_hip");
  const bool selectZones = false;
  const bool pointMesh = true;
  braid2d_mat_test<hip_exec>::test("uniform",
                                   "unibuffer",
                                   "elvira_uniform_unibuffer_pm",
                                   selectZones,
                                   pointMesh);
}

TEST(mir_elvira, elvira_uniform_unibuffer_sel_pm_hip)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_sel_pm_hip");
  const bool selectZones = true;
  const bool pointMesh = true;
  braid2d_mat_test<hip_exec>::test("uniform",
                                   "unibuffer",
                                   "elvira_uniform_unibuffer_sel_pm",
                                   selectZones,
                                   pointMesh);
}
#endif

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
