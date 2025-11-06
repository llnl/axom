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
template <typename ExecSpace>
void braid3d_mat_test(const std::string &type, const std::string &mattype, const std::string &name)
{
  axom::StackArray<axom::IndexType, 3> dims {11, 11, 11};
  axom::StackArray<axom::IndexType, 3> zoneDims {dims[0] - 1, dims[1] - 1, dims[2] - 1};

  // Create the data
  const bool cleanMats = false;
  conduit::Node hostMesh, deviceMesh;
  axom::blueprint::testing::data::braid(type, dims, hostMesh);
  axom::blueprint::testing::data::make_matset(mattype, "mesh", zoneDims, cleanMats, hostMesh);
  utils::copy<ExecSpace>(deviceMesh, hostMesh);
  TestApp.saveVisualization(name + "_orig", hostMesh);

  // Make views.
  auto coordsetView =
    views::make_explicit_coordset<double, 3>::view(deviceMesh["coordsets/coords"]);
  using CoordsetView = decltype(coordsetView);

  using ShapeType = views::HexShape<int>;
  using TopologyView = views::UnstructuredTopologySingleShapeView<ShapeType>;
  auto connView = utils::make_array_view<int>(deviceMesh["topologies/mesh/elements/connectivity"]);
  TopologyView topologyView(connView);

  conduit::Node deviceMIRMesh;
  if(mattype == "unibuffer")
  {
    // clang-format off
    using MatsetView = views::UnibufferMaterialView<int, float, 3>;
    MatsetView matsetView;
    matsetView.set(utils::make_array_view<int>(deviceMesh["matsets/mat/material_ids"]),
                   utils::make_array_view<float>(deviceMesh["matsets/mat/volume_fractions"]),
                   utils::make_array_view<int>(deviceMesh["matsets/mat/sizes"]),
                   utils::make_array_view<int>(deviceMesh["matsets/mat/offsets"]),
                   utils::make_array_view<int>(deviceMesh["matsets/mat/indices"]));
    // clang-format on

    using MIR = axom::mir::EquiZAlgorithm<ExecSpace, TopologyView, CoordsetView, MatsetView>;
    MIR m(topologyView, coordsetView, matsetView);
    conduit::Node options;
    options["matset"] = "mat";
    m.execute(deviceMesh, options, deviceMIRMesh);
  }

  // device->host
  conduit::Node hostMIRMesh;
  utils::copy<seq_exec>(hostMIRMesh, deviceMIRMesh);

  TestApp.saveVisualization(name, hostMIRMesh);

  // Handle baseline comparison.
  constexpr double tolerance = 1.7e-6;
  EXPECT_TRUE(TestApp.test<ExecSpace>(name, hostMIRMesh, tolerance));
}

//------------------------------------------------------------------------------
TEST(mir_equiz, equiz_hex_unibuffer_seq)
{
  AXOM_ANNOTATE_SCOPE("equiz_explicit_hex_seq");
  braid3d_mat_test<seq_exec>("hexs", "unibuffer", "equiz_hex_unibuffer");
}

#if defined(AXOM_USE_OPENMP)
TEST(mir_equiz, equiz_hex_unibuffer_omp)
{
  AXOM_ANNOTATE_SCOPE("equiz_hex_unibuffer_omp");
  braid3d_mat_test<omp_exec>("hexs", "unibuffer", "equiz_hex_unibuffer");
}
#endif

#if defined(AXOM_USE_CUDA)
TEST(mir_equiz, equiz_hex_unibuffer_cuda)
{
  AXOM_ANNOTATE_SCOPE("equiz_hex_unibuffer_cuda");
  braid3d_mat_test<cuda_exec>("hexs", "unibuffer", "equiz_hex_unibuffer");
}
#endif

#if defined(AXOM_USE_HIP)
TEST(mir_equiz, equiz_hex_unibuffer_hip)
{
  AXOM_ANNOTATE_SCOPE("equiz_hex_unibuffer_hip");
  braid3d_mat_test<hip_exec>("hexs", "unibuffer", "equiz_hex_unibuffer");
}
#endif

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
