// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/*! \file mir_elvira2d_impl.hpp
 *  \brief Shared implementation and registrations for the Elvira 2D
 *         execution-policy tests.
 */

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/bump.hpp"
#include "axom/mir.hpp"
#include "axom/primal.hpp"
#include "axom/bump/tests/blueprint_testing_data_helpers.hpp"
#include "axom/bump/tests/blueprint_testing_helpers.hpp"

#include <cstdlib>

inline std::string baselineDirectory()
{
  return pjoin(dataDirectory(), "mir", "regression", "mir_elvira");
}

namespace utils = axom::bump::utilities;
namespace views = axom::bump::views;

//------------------------------------------------------------------------------
// Global test application object.
extern axom::blueprint::testing::TestApplication TestApp;

//------------------------------------------------------------------------------
template <typename ExecSpace>
double sumValues(axom::ArrayView<double> values)
{
  axom::ReduceSum<ExecSpace, double> reduce_values(0.);
  axom::for_all<ExecSpace>(values.size(),
                           [=] AXOM_HOST_DEVICE(axom::IndexType i) { reduce_values += values[i]; });
  return reduce_values.get();
}

//------------------------------------------------------------------------------
inline void makeSingleZoneCoordset2D(conduit::Node& n_mesh,
                                     const std::string& coordsetName,
                                     const axom::StackArray<conduit::index_t, 2>& dims)
{
  SLIC_ASSERT(dims[0] >= 2 && dims[1] >= 2);

  const conduit::index_t nnodes = dims[0] * dims[1];
  conduit::Node& n_coordset = n_mesh["coordsets/" + coordsetName];
  n_coordset["type"] = "explicit";
  n_coordset["values/x"].set(conduit::DataType::float32(nnodes));
  n_coordset["values/y"].set(conduit::DataType::float32(nnodes));

  float* x = n_coordset["values/x"].as_float_ptr();
  float* y = n_coordset["values/y"].as_float_ptr();

  conduit::index_t index = 0;
  for(conduit::index_t j = 0; j < dims[1]; j++)
  {
    for(conduit::index_t i = 0; i < dims[0]; i++, index++)
    {
      x[index] = static_cast<float>(i);
      y[index] = static_cast<float>(j);
    }
  }
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct single_zone_elvira2d_test
{
  static void initialize(bool stridedStructured, bool cleanMats, conduit::Node& n_mesh)
  {
    constexpr conduit::index_t NLEFT = 2;
    constexpr conduit::index_t NRIGHT = 1;

    const axom::StackArray<axom::IndexType, 2> zoneDims {1, 1};
    const axom::StackArray<conduit::index_t, 2> nodeDims =
      stridedStructured ? axom::StackArray<conduit::index_t, 2> {zoneDims[0] + 1 + NLEFT + NRIGHT,
                                                                 zoneDims[1] + 1 + NLEFT + NRIGHT}
                        : axom::StackArray<conduit::index_t, 2> {zoneDims[0] + 1, zoneDims[1] + 1};

    makeSingleZoneCoordset2D(n_mesh, "coords", nodeDims);

    conduit::Node& n_topology = n_mesh["topologies/mesh"];
    n_topology["type"] = "structured";
    n_topology["coordset"] = "coords";
    n_topology["elements/dims/i"] = zoneDims[0];
    n_topology["elements/dims/j"] = zoneDims[1];
    if(stridedStructured)
    {
      const axom::StackArray<conduit::index_t, 2> offsets {NLEFT, NLEFT};
      const axom::StackArray<conduit::index_t, 2> strides {1, nodeDims[0]};
      n_topology["elements/dims/offsets"].set(offsets.data(), offsets.size());
      n_topology["elements/dims/strides"].set(strides.data(), strides.size());
    }

    axom::blueprint::testing::data::make_matset("unibuffer",
                                                "mesh",
                                                zoneDims,
                                                cleanMats,
                                                false,
                                                n_mesh);
  }

  static void test(bool stridedStructured, bool cleanMats)
  {
    conduit::Node hostMesh, deviceMesh;
    initialize(stridedStructured, cleanMats, hostMesh);
    utils::copy<ExecSpace>(deviceMesh, hostMesh);

    const conduit::Node& n_coordset = deviceMesh["coordsets/coords"];
    const conduit::Node& n_topology = deviceMesh["topologies/mesh"];
    const conduit::Node& n_matset = deviceMesh["matsets/mat"];

    auto coordsetView = views::make_explicit_coordset<float, 2>::view(n_coordset);
    auto matsetView = views::make_unibuffer_matset<int, float, 3>::view(n_matset);
    const auto expectedFragments = matsetView.numberOfMaterials(0);
    if(cleanMats)
    {
      EXPECT_EQ(expectedFragments, 1);
    }
    else
    {
      EXPECT_GT(expectedFragments, 1);
    }

    conduit::Node deviceMIRMesh;
    conduit::Node options;
    options["matset"] = "mat";
    options["plane"] = 1;
    options["pointmesh"] = 0;
    options["verbose"] = 1;

    if(stridedStructured)
    {
      auto topologyView = views::make_strided_structured_topology<2>::view(n_topology);
      execute(topologyView, coordsetView, matsetView, deviceMesh, options, deviceMIRMesh);
    }
    else
    {
      auto topologyView = views::make_structured_topology<2>::view(n_topology);
      execute(topologyView, coordsetView, matsetView, deviceMesh, options, deviceMIRMesh);
    }

    conduit::Node hostMIRMesh;
    utils::copy<seq_exec>(hostMIRMesh, deviceMIRMesh);

    conduit::Node info;
    const bool verifyOK = conduit::blueprint::mesh::verify(hostMIRMesh, info);
    if(!verifyOK)
    {
      printNode(hostMIRMesh);
      info.print();
    }
    EXPECT_TRUE(verifyOK);

    const conduit::Node& n_mir_topology = deviceMIRMesh["topologies/mesh"];
    const conduit::Node& n_mir_coordset = deviceMIRMesh["coordsets/coords"];
    auto mirCoordsetView = views::make_explicit_coordset<float, 2>::view(n_mir_coordset);

    if(cleanMats)
    {
      if(stridedStructured)
      {
        auto mirTopoView = views::make_strided_structured_topology<2>::view(n_mir_topology);
        EXPECT_EQ(mirTopoView.numberOfZones(), expectedFragments);
        axom::bump::MakeZoneVolumes<ExecSpace, decltype(mirTopoView), decltype(mirCoordsetView)> makeArea(
          mirTopoView,
          mirCoordsetView);
        makeArea.execute(n_mir_topology, n_mir_coordset, deviceMIRMesh["fields/area"]);
      }
      else
      {
        auto mirTopoView = views::make_structured_topology<2>::view(n_mir_topology);
        EXPECT_EQ(mirTopoView.numberOfZones(), expectedFragments);
        axom::bump::MakeZoneVolumes<ExecSpace, decltype(mirTopoView), decltype(mirCoordsetView)> makeArea(
          mirTopoView,
          mirCoordsetView);
        makeArea.execute(n_mir_topology, n_mir_coordset, deviceMIRMesh["fields/area"]);
      }
    }
    else
    {
      auto mirTopoView =
        views::make_unstructured_single_shape_topology<views::PolygonShape<axom::IndexType>>::view(
          n_mir_topology);
      EXPECT_EQ(mirTopoView.numberOfZones(), expectedFragments);
      axom::bump::MakeZoneVolumes<ExecSpace, decltype(mirTopoView), decltype(mirCoordsetView)> makeArea(
        mirTopoView,
        mirCoordsetView);
      makeArea.execute(n_mir_topology, n_mir_coordset, deviceMIRMesh["fields/area"]);
    }

    const auto areaView = utils::make_array_view<double>(deviceMIRMesh["fields/area/values"]);
    EXPECT_NEAR(sumValues<ExecSpace>(areaView), 1., 1.e-6);
  }

private:
  template <typename TopologyView, typename CoordsetView, typename MatsetView>
  static void execute(const TopologyView& topologyView,
                      const CoordsetView& coordsetView,
                      const MatsetView& matsetView,
                      conduit::Node& deviceMesh,
                      const conduit::Node& options,
                      conduit::Node& deviceMIRMesh)
  {
    using IndexingPolicy = typename TopologyView::IndexingPolicy;
    using MIR = axom::mir::ElviraAlgorithm<ExecSpace, IndexingPolicy, CoordsetView, MatsetView>;
    MIR m(topologyView, coordsetView, matsetView);
    m.execute(deviceMesh, options, deviceMIRMesh);
  }
};

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
  static void initialize(const std::string& type,
                         const std::string& mattype,
                         bool cleanMats,
                         conduit::Node& n_mesh)
  {
    axom::StackArray<axom::IndexType, 2> dims {10, 10};
    axom::StackArray<axom::IndexType, 2> zoneDims {dims[0] - 1, dims[1] - 1};
    axom::blueprint::testing::data::braid(type, dims, n_mesh);
    const bool makeMixedField = false;  // for now
    axom::blueprint::testing::data::make_matset(mattype,
                                                "mesh",
                                                zoneDims,
                                                cleanMats,
                                                makeMixedField,
                                                n_mesh);
  }

  // Select a chunk of clean and mixed zones.
  static void selectZones(conduit::Node& n_options)
  {
    n_options["selectedZones"].set(std::vector<axom::IndexType> {30, 31, 32, 39, 40, 41, 48, 49, 50});
  }

  static void test(const std::string& type,
                   const std::string& mattype,
                   const std::string& name,
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
      conduit::Node& hostDomain = (nDomains > 1) ? hostMesh[domainName] : hostMesh;

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
        conduit::Node& deviceDomain = (nDomains > 1) ? deviceMesh[domainName] : deviceMesh;

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
          options["verbose"] = 1;
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
        axom::for_all<ExecSpace>(N, [=] AXOM_HOST_DEVICE(axom::IndexType index) {
          arrView[index] = index * index;
        });
      }
      axom::synchronize<ExecSpace>();
    }
  }
};

//------------------------------------------------------------------------------
