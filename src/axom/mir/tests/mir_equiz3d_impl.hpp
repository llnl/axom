// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/*! \file mir_equiz3d_impl.hpp
 *  \brief Implementation shared by the EquiZ 3D execution-policy tests.
 */

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/bump.hpp"
#include "axom/mir.hpp"
#include "axom/primal.hpp"
#include "axom/bump/tests/blueprint_testing_data_helpers.hpp"
#include "axom/bump/tests/blueprint_testing_helpers.hpp"

namespace utils = axom::bump::utilities;
namespace views = axom::bump::views;

inline std::string baselineDirectory()
{
  return pjoin(dataDirectory(), "mir", "regression", "mir_equiz");
}

//------------------------------------------------------------------------------
// Global test application object.
extern axom::blueprint::testing::TestApplication TestApp;

//------------------------------------------------------------------------------
template <typename ExecSpace>
void braid3d_mat_test(const std::string& type, const std::string& mattype, const std::string& name)
{
  axom::StackArray<axom::IndexType, 3> dims {11, 11, 11};
  axom::StackArray<axom::IndexType, 3> zoneDims {dims[0] - 1, dims[1] - 1, dims[2] - 1};

  // Create the data
  const bool cleanMats = false;
  conduit::Node hostMesh, deviceMesh;
  axom::blueprint::testing::data::braid(type, dims, hostMesh);
  const bool makeMixedField = false;  // for now
  axom::blueprint::testing::data::make_matset(mattype,
                                              "mesh",
                                              zoneDims,
                                              cleanMats,
                                              makeMixedField,
                                              hostMesh);
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
    options["verbose"] = 1;
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
