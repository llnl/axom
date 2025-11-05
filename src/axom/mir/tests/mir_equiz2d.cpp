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

namespace utils = axom::bump::utilities;
namespace views = axom::bump::views;

std::string baselineDirectory() { return pjoin(dataDirectory(), "mir", "regression", "mir_equiz"); }

//------------------------------------------------------------------------------
// Global test application object.
axom::blueprint::testing::TestApplication TestApp;

//------------------------------------------------------------------------------
TEST(mir_equiz, miralgorithm)
{
  axom::mir::MIRAlgorithm *m = nullptr;
  EXPECT_EQ(m, nullptr);
}

//------------------------------------------------------------------------------
TEST(mir_equiz, materialinformation)
{
  conduit::Node matset;
  matset["material_map/a"] = 1;
  matset["material_map/b"] = 2;
  matset["material_map/c"] = 0;

  auto mi = axom::bump::views::materials(matset);
  EXPECT_EQ(mi.size(), 3);
  EXPECT_EQ(mi[0].number, 1);
  EXPECT_EQ(mi[0].name, "a");

  EXPECT_EQ(mi[1].number, 2);
  EXPECT_EQ(mi[1].name, "b");

  EXPECT_EQ(mi[2].number, 0);
  EXPECT_EQ(mi[2].name, "c");
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
void braid2d_mat_test(const std::string &type,
                      const std::string &mattype,
                      const std::string &name,
                      int nDomains = 1)
{
  axom::StackArray<axom::IndexType, 2> dims {10, 10};
  axom::StackArray<axom::IndexType, 2> zoneDims {dims[0] - 1, dims[1] - 1};

  // Create the data (make 1+ domains of the same thing)
  conduit::Node hostMesh, deviceMesh;
  for(int dom = 0; dom < nDomains; dom++)
  {
    const std::string domainName = axom::fmt::format("domain_{:07}", dom);
    conduit::Node &hostDomain = (nDomains > 1) ? hostMesh[domainName] : hostMesh;
    const bool cleanMats = false;
    axom::blueprint::testing::data::braid(type, dims, hostDomain);
    axom::blueprint::testing::data::make_matset(mattype, "mesh", zoneDims, cleanMats, hostDomain);
    TestApp.saveVisualization(name + "_orig", hostDomain);
  }

  // host->device
  utils::copy<ExecSpace>(deviceMesh, hostMesh);

  for(int dom = 0; dom < nDomains; dom++)
  {
    const std::string domainName = axom::fmt::format("domain_{:07}", dom);
    conduit::Node &deviceDomain = (nDomains > 1) ? deviceMesh[domainName] : deviceMesh;

    // Make views.
    auto coordsetView = views::make_uniform_coordset<2>::view(deviceDomain["coordsets/coords"]);
    auto topologyView = views::make_uniform_topology<2>::view(deviceDomain["topologies/mesh"]);
    using CoordsetView = decltype(coordsetView);
    using TopologyView = decltype(topologyView);

    conduit::Node deviceMIRDomain;
    if(mattype == "unibuffer")
    {
      // clang-format off
      using MatsetView = views::UnibufferMaterialView<int, float, 3>;
      MatsetView matsetView;
      matsetView.set(utils::make_array_view<int>(deviceDomain["matsets/mat/material_ids"]),
                     utils::make_array_view<float>(deviceDomain["matsets/mat/volume_fractions"]),
                     utils::make_array_view<int>(deviceDomain["matsets/mat/sizes"]),
                     utils::make_array_view<int>(deviceDomain["matsets/mat/offsets"]),
                     utils::make_array_view<int>(deviceDomain["matsets/mat/indices"]));
      // clang-format on

      using MIR = axom::mir::EquiZAlgorithm<ExecSpace, TopologyView, CoordsetView, MatsetView>;
      MIR m(topologyView, coordsetView, matsetView);
      conduit::Node options;
      options["matset"] = "mat";
      m.execute(deviceDomain, options, deviceMIRDomain);
    }

    // device->host for the current domain
    conduit::Node hostMIRDomain;
    utils::copy<seq_exec>(hostMIRDomain, deviceMIRDomain);

    TestApp.saveVisualization(name, hostMIRDomain);

    // Handle baseline comparison.
    constexpr double tolerance = 2.6e-06;
    EXPECT_TRUE(TestApp.test<ExecSpace>(name, hostMIRDomain, tolerance));
  }
}

//------------------------------------------------------------------------------
TEST(mir_equiz, equiz_uniform_unibuffer_seq)
{
  AXOM_ANNOTATE_SCOPE("equiz_uniform_unibuffer_seq");
  braid2d_mat_test<seq_exec>("uniform", "unibuffer", "equiz_uniform_unibuffer");
  braid2d_mat_test<seq_exec>("uniform", "unibuffer", "equiz_uniform_unibuffer", 2);
}

#if defined(AXOM_USE_OPENMP)
TEST(mir_equiz, equiz_uniform_unibuffer_omp)
{
  AXOM_ANNOTATE_SCOPE("equiz_uniform_unibuffer_omp");
  braid2d_mat_test<omp_exec>("uniform", "unibuffer", "equiz_uniform_unibuffer");
  braid2d_mat_test<omp_exec>("uniform", "unibuffer", "equiz_uniform_unibuffer", 2);
}
#endif

#if defined(AXOM_USE_CUDA)
TEST(mir_equiz, equiz_uniform_unibuffer_cuda)
{
  AXOM_ANNOTATE_SCOPE("equiz_uniform_unibuffer_cuda");
  braid2d_mat_test<cuda_exec>("uniform", "unibuffer", "equiz_uniform_unibuffer");
  braid2d_mat_test<cuda_exec>("uniform", "unibuffer", "equiz_uniform_unibuffer", 2);
}
#endif

#if defined(AXOM_USE_HIP)
TEST(mir_equiz, equiz_uniform_unibuffer_hip)
{
  AXOM_ANNOTATE_SCOPE("equiz_uniform_unibuffer_hip");
  braid2d_mat_test<hip_exec>("uniform", "unibuffer", "equiz_uniform_unibuffer");
  braid2d_mat_test<hip_exec>("uniform", "unibuffer", "equiz_uniform_unibuffer", 2);
}
#endif

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
