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

namespace bump = axom::bump;
namespace utils = axom::bump::utilities;
namespace views = axom::bump::views;

std::string baselineDirectory()
{
  return pjoin(dataDirectory(), "mir", "regression", "mir_elvira3d");
}

//------------------------------------------------------------------------------
// Global test application object.
axom::blueprint::testing::TestApplication TestApp;

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_Elvira3D
{
  static const int gridSize = 10;
  static const int numSpheres = 2;

  static void initialize(conduit::Node &n_mesh)
  {
    AXOM_ANNOTATE_SCOPE("initialize");
    axom::bump::data::MeshTester M;
    M.setStructured(true);
    M.initTestCaseSix(gridSize, numSpheres, n_mesh);
  }

  // Select a chunk of zones.
  static int selectZones(conduit::Node &n_options)
  {
    std::vector<axom::IndexType> selected;
    for(int k = 0; k < gridSize; k++)
    {
      for(int j = 0; j < gridSize; j++)
      {
        for(int i = 0; i < gridSize; i++)
        {
          // Save all but an octant so we can see inside.
          if(i > gridSize / 3 && j > gridSize / 3 && k > gridSize / 3)
          {
            continue;
          }
          auto idx = k * gridSize * gridSize + j * gridSize + i;
          selected.push_back(idx);
        }
      }
    }
    n_options["selectedZones"].set(selected);
    return static_cast<int>(selected.size());
  }

  static void test(const std::string &name, bool selectedZones = false, bool pointMesh = false)
  {
    const double expectedVolume = gridSize * gridSize * gridSize;
    double mirExpectedVolume = expectedVolume;

    // Create the data
    conduit::Node hostMesh, deviceMesh;
    initialize(hostMesh);
    {
      AXOM_ANNOTATE_SCOPE("host_to_device");
      utils::copy<ExecSpace>(deviceMesh, hostMesh);
    }
    // Save visualization, if enabled.
    TestApp.saveVisualization(name + "_orig", hostMesh);

    //--------------------------------------------------------------------------
    const conduit::Node &n_coordset = deviceMesh.fetch_existing("coordsets/coords");
    const conduit::Node &n_topology = deviceMesh.fetch_existing("topologies/mesh");
    const conduit::Node &n_matset = deviceMesh.fetch_existing("matsets/mat");

    // Make views.
    auto coordsetView = views::make_explicit_coordset<float, 3>::view(n_coordset);
    using CoordsetView = decltype(coordsetView);

    auto topologyView = views::make_structured_topology<3>::view(n_topology);
    using TopologyView = decltype(topologyView);
    using IndexingPolicy = typename TopologyView::IndexingPolicy;

    auto matsetView = views::make_unibuffer_matset<int, float, 3>::view(n_matset);
    using MatsetView = decltype(matsetView);

    // Do MIR
    using MIR = axom::mir::ElviraAlgorithm<ExecSpace, IndexingPolicy, CoordsetView, MatsetView>;
    MIR m(topologyView, coordsetView, matsetView);
    conduit::Node deviceMIRMesh;
    conduit::Node options;
    options["matset"] = "mat";
    options["plane"] = pointMesh ? 1 : 0;
    options["pointmesh"] = pointMesh ? 1 : 0;
    // Be more lenient in how far away points are in order to combine them.
    options["point_tolerance"] = 1.e-4;
    if(selectedZones)
    {
      // The MIR volume changes when the set of selected zones changes.
      mirExpectedVolume = selectZones(options);
    }
    m.execute(deviceMesh, options, deviceMIRMesh);

    if(pointMesh)
    {
      comparePointMesh(name, deviceMIRMesh);
    }
    else
    {
      compare(name,
              selectedZones,
              deviceMesh,
              topologyView,
              coordsetView,
              n_topology,
              n_coordset,
              deviceMIRMesh,
              expectedVolume,
              mirExpectedVolume);
    }
  }

  static void comparePointMesh(const std::string &name, const conduit::Node &deviceMIRMesh)
  {
    // device->host
    conduit::Node hostMIRMesh;
    {
      AXOM_ANNOTATE_SCOPE("device_to_host");
      utils::copy<seq_exec>(hostMIRMesh, deviceMIRMesh);
    }

    // Save visualization of MIR mesh, if enabled.
    TestApp.saveVisualization(name, hostMIRMesh);

    // Handle baseline comparison.
    constexpr double tolerance = 8.e-06;
    EXPECT_TRUE(TestApp.test<ExecSpace>(name, hostMIRMesh, tolerance));
  }

  template <typename TopologyView, typename CoordsetView>
  static void compare(const std::string &name,
                      bool selectedZones,
                      conduit::Node &deviceMesh,
                      const TopologyView &topologyView,
                      const CoordsetView &coordsetView,
                      const conduit::Node &n_topology,
                      const conduit::Node &n_coordset,
                      conduit::Node &deviceMIRMesh,
                      double expectedVolume,
                      double mirExpectedVolume)
  {
    //--------------------------------------------------------------------------
    // Compute volumes for original mesh as a field.
    AXOM_ANNOTATE_BEGIN("volume");
    bump::MakeZoneVolumes<ExecSpace, TopologyView, CoordsetView> origZV(topologyView, coordsetView);
    origZV.execute(n_topology, n_coordset, deviceMesh["fields/volume"]);

    //--------------------------------------------------------------------------
    // Compute volumes for MIR mesh as a field.
    conduit::Node &n_mir_coordset = deviceMIRMesh["coordsets/coords"];
    auto mirCoordsetView = views::make_explicit_coordset<float, 3>::view(n_mir_coordset);
    using MirCoordsetView = decltype(mirCoordsetView);

    // Make polyhedral topology view.
    const conduit::Node &n_mir_topology = deviceMIRMesh["topologies/mesh"];
    auto mirTopoView =
      views::make_unstructured_polyhedral_topology<axom::IndexType>::view(n_mir_topology);
    using MirTopologyView = decltype(mirTopoView);

    bump::MakeZoneVolumes<ExecSpace, MirTopologyView, MirCoordsetView> mirZV(mirTopoView,
                                                                             mirCoordsetView);
    mirZV.execute(n_mir_topology, n_mir_coordset, deviceMIRMesh["fields/volume"]);
    AXOM_ANNOTATE_END("volume");

    //--------------------------------------------------------------------------
    // device->host
    conduit::Node hostMIRMesh;
    {
      AXOM_ANNOTATE_SCOPE("device_to_host");
      utils::copy<seq_exec>(hostMIRMesh, deviceMIRMesh);
    }

    // Save visualization of MIR mesh, if enabled.
    TestApp.saveVisualization(name, hostMIRMesh);

    //--------------------------------------------------------------------------
    // Handle baseline comparison.
#if 0
    // NOTE: Comparing against this baseline is turned off for now. Rather than
    //       compare the Conduit nodes, we skip compare material volumes
    //       before/after MIR.
    constexpr double tolerance = 2.6e-06;
    EXPECT_TRUE(TestApp.test<ExecSpace>(name, hostMIRMesh, tolerance));
#endif
    const conduit::Node &n_matset = deviceMesh["matsets/mat"];
    auto matsetView = views::make_unibuffer_matset<int, float, 3>::view(n_matset);

    const conduit::Node &n_mir_matset = deviceMIRMesh["matsets/mat"];
    auto mirMatsetView = views::make_unibuffer_matset<int, float, 3>::view(n_mir_matset);

    //--------------------------------------------------------------------------
    // Compute the total volumes on the original and MIR meshes.
    constexpr double tolerance = 3.e-5;

    const auto orig_volume = utils::make_array_view<double>(deviceMesh["fields/volume/values"]);
    EXPECT_NEAR(expectedVolume, variableSum(orig_volume), tolerance);
    const auto origMatInfo = views::materials(n_matset);
    const auto origTotalVolumes = sumMaterialVolumes(matsetView, orig_volume, origMatInfo);

    const auto mir_volume = utils::make_array_view<double>(deviceMIRMesh["fields/volume/values"]);
    EXPECT_NEAR(mirExpectedVolume, variableSum(mir_volume), tolerance);
    const auto mirMatInfo = views::materials(n_mir_matset);
    const auto mirTotalVolumes = sumMaterialVolumes(mirMatsetView, mir_volume, mirMatInfo);

    //--------------------------------------------------------------------------
    // comparisons
    EXPECT_EQ(origTotalVolumes.size(), mirTotalVolumes.size());
    // Expected values for the total volumes when we use selected zones.
    const double selectedZonesTotalVolume[] = {35.73998360180606, 292.443992377414, 455.8160238981947};
    double volumeSums[2] = {0., 0.};
    for(size_t i = 0; i < origTotalVolumes.size(); i++)
    {
      const double origTotalVol = selectedZones ? selectedZonesTotalVolume[i] : origTotalVolumes[i];

      SLIC_INFO(
        axom::fmt::format("Material {}: origVF = {}, mirVF = {}", i, origTotalVol, mirTotalVolumes[i]));

      volumeSums[0] += origTotalVol;
      volumeSums[1] += mirTotalVolumes[i];
    }
    EXPECT_NEAR(volumeSums[0], selectedZones ? mirExpectedVolume : expectedVolume, tolerance);
    EXPECT_NEAR(volumeSums[1], mirExpectedVolume, tolerance);
    for(size_t i = 0; i < origTotalVolumes.size(); i++)
    {
      const double origTotalVol = selectedZones ? selectedZonesTotalVolume[i] : origTotalVolumes[i];

      EXPECT_NEAR(origTotalVol, mirTotalVolumes[i], tolerance);
    }
  }

  /*!
   * \brief Sums the input array view.
   *
   * \param var An array view to sum.
   *
   * \return The sum of the input array view.
   */
  static double variableSum(axom::ArrayView<double> var)
  {
    axom::ReduceSum<ExecSpace, double> reduceVar(0.);
    axom::for_all<ExecSpace>(var.size(), AXOM_LAMBDA(axom::IndexType i) { reduceVar += var[i]; });
    return reduceVar.get();
  }

  /*!
   * \brief Compute the total volumes for the materials in the mesh.
   *
   * \param matsetView The matset view that contains the material data (on device)
   * \param zoneVolumes An array view that contains zone volumes (on device)
   * \param matInfo The material information for the mesh.
   *
   * \return A vector of volumes on the host, 1 element per material.
   */
  template <typename MatsetView>
  static std::vector<double> sumMaterialVolumes(MatsetView matsetView,
                                                axom::ArrayView<double> zoneVolumes,
                                                const views::MaterialInformation &matInfo)
  {
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    AXOM_ANNOTATE_SCOPE("sumMaterialVolumes");

    // Make a sorted list of material numbers on device.
    const int nmats = static_cast<int>(matInfo.size());
    axom::Array<int> sortedIdsHost(matInfo.size());
    int mi = 0;
    for(const auto &mat : matInfo)
    {
      sortedIdsHost[mi++] = mat.number;
    }
    axom::utilities::Sorting<int>::sort(sortedIdsHost.data(), sortedIdsHost.size());
    axom::Array<int> sortedIds(nmats, nmats, allocatorID);
    axom::copy(sortedIds.data(), sortedIdsHost.data(), nmats * sizeof(int));
    auto sortedIdsView = sortedIds.view();

    // Compute the total volumes for each material.
    axom::Array<double> totalVolume(nmats, nmats, allocatorID);
    auto totalVolumeView = totalVolume.view();
    axom::for_all<ExecSpace>(nmats, AXOM_LAMBDA(axom::IndexType i) { totalVolumeView[i] = 0.; });
    axom::for_all<ExecSpace>(
      matsetView.numberOfZones(),
      AXOM_LAMBDA(axom::IndexType zi) {
        // Add the material volumes to the total volumes.
        const auto end = matsetView.endZone(zi);
        for(auto zoneMat = matsetView.beginZone(zi); zoneMat != end; zoneMat++)
        {
          auto index = axom::utilities::binary_search(sortedIdsView, zoneMat.material_id());
          // RelWithDebInfo workaround - "sortedIdsView.size()" substitutes lambda capture device failure for "nmats"
          SLIC_ASSERT(index >= 0 && index < sortedIdsView.size());

          // Use an atomic to sum the value.
          axom::atomicAdd<ExecSpace>(totalVolumeView.data() + index,
                                     zoneVolumes[zi] * zoneMat.volume_fraction());
        }
      });

    // Move results back to host.
    std::vector<double> hostVolumes(nmats, 0.);
    axom::copy(hostVolumes.data(), totalVolume.data(), nmats * sizeof(double));

    return hostVolumes;
  }
};

//------------------------------------------------------------------------------
TEST(mir_elvira3d, elvira3d_unibuffer_seq)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_seq");
  const bool selectZones = false;
  const bool pointMesh = false;
  test_Elvira3D<seq_exec>::test("elvira3d_unibuffer", selectZones, pointMesh);
}

TEST(mir_elvira3d, elvira3d_unibuffer_sel_seq)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_sel_seq");
  const bool selectZones = true;
  const bool pointMesh = false;
  test_Elvira3D<seq_exec>::test("elvira3d_unibuffer_sel", selectZones, pointMesh);
}

TEST(mir_elvira3d, elvira3d_unibuffer_pm_seq)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_pm_seq");
  const bool selectZones = false;
  const bool pointMesh = true;
  test_Elvira3D<seq_exec>::test("elvira3d_unibuffer_pm", selectZones, pointMesh);
}

TEST(mir_elvira3d, elvira3d_unibuffer_sel_pm_seq)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_sel_pm_seq");
  const bool selectZones = true;
  const bool pointMesh = true;
  test_Elvira3D<seq_exec>::test("elvira3d_unibuffer_sel_pm", selectZones, pointMesh);
}

#if defined(AXOM_USE_OPENMP)
TEST(mir_elvira3d, elvira3d_unibuffer_omp)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_omp");
  const bool selectZones = false;
  const bool pointMesh = false;
  test_Elvira3D<omp_exec>::test("elvira3d_unibuffer", selectZones, pointMesh);
}

TEST(mir_elvira3d, elvira3d_unibuffer_sel_omp)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_sel_omp");
  const bool selectZones = true;
  const bool pointMesh = false;
  test_Elvira3D<omp_exec>::test("elvira3d_unibuffer_sel", selectZones, pointMesh);
}

TEST(mir_elvira3d, elvira3d_unibuffer_pm_omp)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_pm_omp");
  const bool selectZones = false;
  const bool pointMesh = true;
  test_Elvira3D<omp_exec>::test("elvira3d_unibuffer_pm", selectZones, pointMesh);
}

TEST(mir_elvira3d, elvira3d_unibuffer_sel_pm_omp)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_sel_pm_omp");
  const bool selectZones = true;
  const bool pointMesh = true;
  test_Elvira3D<omp_exec>::test("elvira3d_unibuffer_sel_pm", selectZones, pointMesh);
}
#endif

#if defined(AXOM_USE_CUDA)
TEST(mir_elvira3d, elvira3d_unibuffer_cuda)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_cuda");
  const bool selectZones = false;
  const bool pointMesh = false;
  test_Elvira3D<cuda_exec>::test("elvira3d_unibuffer", selectZones, pointMesh);
}

TEST(mir_elvira3d, elvira3d_unibuffer_sel_cuda)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_sel_cuda");
  const bool selectZones = true;
  const bool pointMesh = false;
  test_Elvira3D<cuda_exec>::test("elvira3d_unibuffer_sel", selectZones, pointMesh);
}

TEST(mir_elvira3d, elvira3d_unibuffer_pm_cuda)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_pm_cuda");
  const bool selectZones = false;
  const bool pointMesh = true;
  test_Elvira3D<cuda_exec>::test("elvira3d_unibuffer_pm", selectZones, pointMesh);
}

TEST(mir_elvira3d, elvira3d_unibuffer_sel_pm_cuda)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_sel_pm_cuda");
  const bool selectZones = true;
  const bool pointMesh = true;
  test_Elvira3D<cuda_exec>::test("elvira3d_unibuffer_sel_pm", selectZones, pointMesh);
}
#endif

#if defined(AXOM_USE_HIP)
TEST(mir_elvira3d, elvira3d_unibuffer_hip)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_hip");
  const bool selectZones = false;
  const bool pointMesh = false;
  test_Elvira3D<hip_exec>::test("elvira3d_unibuffer", selectZones, pointMesh);
}

TEST(mir_elvira3d, elvira3d_unibuffer_sel_hip)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_sel_hip");
  const bool selectZones = true;
  const bool pointMesh = false;
  test_Elvira3D<hip_exec>::test("elvira3d_unibuffer_sel", selectZones, pointMesh);
}

TEST(mir_elvira3d, elvira3d_unibuffer_pm_hip)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_pm_hip");
  const bool selectZones = false;
  const bool pointMesh = true;
  test_Elvira3D<hip_exec>::test("elvira3d_unibuffer_pm", selectZones, pointMesh);
}

TEST(mir_elvira3d, elvira3d_unibuffer_sel_pm_hip)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_sel_pm_hip");
  const bool selectZones = true;
  const bool pointMesh = true;
  test_Elvira3D<hip_exec>::test("elvira3d_unibuffer_sel_pm", selectZones, pointMesh);
}
#endif

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
