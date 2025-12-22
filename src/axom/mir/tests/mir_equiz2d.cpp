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
namespace bump = axom::bump;

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
                      int nDomains,
                      bool selectedZones,
                      bool cleanMats)
{
  axom::StackArray<axom::IndexType, 2> dims {10, 10};
  axom::StackArray<axom::IndexType, 2> zoneDims {dims[0] - 1, dims[1] - 1};

  // Create the data (make 1+ domains of the same thing)
  conduit::Node hostMesh, deviceMesh;
  for(int dom = 0; dom < nDomains; dom++)
  {
    const std::string domainName = axom::fmt::format("domain_{:07}", dom);
    conduit::Node &hostDomain = (nDomains > 1) ? hostMesh[domainName] : hostMesh;
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
      if(cleanMats)
      {
        // Set the output names
        options["topologyName"] = "postmir_topology";
        options["coordsetName"] = "postmir_coords";
        options["matsetName"] = "postmir_matset";
      }
      if(selectedZones)
      {
        options["selectedZones"].set(
          std::vector<axom::IndexType> {30, 31, 32, 39, 40, 41, 48, 49, 50});
      }
      m.execute(deviceDomain, options, deviceMIRDomain);
    }

    // device->host for the current domain
    conduit::Node hostMIRDomain;
    utils::copy<seq_exec>(hostMIRDomain, deviceMIRDomain);

    // Verify the hostMIRMesh to look for errors.
    conduit::Node info;
    bool verifyOK = conduit::blueprint::mesh::verify(hostMIRDomain, info);
    if(!verifyOK)
    {
      printNode(hostMIRDomain);
      info.print();
    }
    EXPECT_TRUE(verifyOK);

    TestApp.saveVisualization(name, hostMIRDomain);

    // Handle baseline comparison.
    constexpr double tolerance = 2.6e-06;
    EXPECT_TRUE(TestApp.test<ExecSpace>(name, hostMIRDomain, tolerance));
  }
}
//------------------------------------------------------------------------------
/*!
 * \brief Tests the MIR+TopologyMapper on polygonal geometry.
 *
 *        1. Make polygonal mesh 1 with clean matset
 *        2. Make mesh 2, a rotated version of polygonal mesh 1
 *        3. Map material from mesh 1 onto mesh 2
 *        4. Run MIR on mesh 2
 */
template <typename ExecSpace>
class test_Polygonal_MIR
{
public:
  static constexpr conduit::index_t NLEVELS = 4;
  static constexpr int MAX_MATERIALS = NLEVELS + 1;

  static void test(const std::string &name)
  {
    // Make the 2D input mesh.
    conduit::Node n_mesh;
    initialize(n_mesh);

    // host->device
    conduit::Node n_dev;
    utils::copy<ExecSpace>(n_dev, n_mesh);

    mapping_target2(n_dev);
    mir_target2(n_dev);

    // device->host
    conduit::Node hostResult;
    utils::copy<seq_exec>(hostResult, n_dev);

    EXPECT_EQ(countBadMaterialZones(hostResult["matsets/target2_matset"]), 0);

    TestApp.saveVisualization(name, hostResult);

    // Handle baseline comparison.
    EXPECT_TRUE(TestApp.test<ExecSpace>(name, hostResult));
  }

  static void initialize(conduit::Node &n_mesh)
  {
    // Make polygonal geometry
    const conduit::index_t nz = 1;
    conduit::blueprint::mesh::examples::polytess(NLEVELS, nz, n_mesh);

    // Make a matset from the level field.
    conduit::Node &n_matset = n_mesh["matsets/mat"];
    n_matset["topology"] = "topo";
    for(int mat = 1; mat <= NLEVELS; mat++)
    {
      n_matset[axom::fmt::format("material_map/mat{}", mat)] = mat;
    }
    const auto values = n_mesh["fields/level/values"].as_int_accessor();
    const int nzones = values.number_of_elements();
    n_matset["material_ids"].set(conduit::DataType::int32(nzones));
    n_matset["indices"].set(conduit::DataType::int32(nzones));
    n_matset["sizes"].set(conduit::DataType::int32(nzones));
    n_matset["offsets"].set(conduit::DataType::int32(nzones));
    n_matset["volume_fractions"].set(conduit::DataType::float32(nzones));

    auto material_ids = n_matset["material_ids"].as_int32_ptr();
    auto indices = n_matset["indices"].as_int32_ptr();
    auto sizes = n_matset["sizes"].as_int32_ptr();
    auto offsets = n_matset["offsets"].as_int32_ptr();
    auto volume_fractions = n_matset["volume_fractions"].as_float32_ptr();
    for(int i = 0; i < nzones; i++)
    {
      material_ids[i] = values[i];
      indices[i] = i;
      sizes[i] = 1;
      offsets[i] = i;
      volume_fractions[i] = 1.f;
    }

    make_target2(n_mesh);
  }

  static void make_target2(conduit::Node &n_mesh)
  {
    const auto x = n_mesh["coordsets/coords/values/x"].as_float64_accessor();
    const auto y = n_mesh["coordsets/coords/values/y"].as_float64_accessor();

    // Make a rotated copy of the input topo mesh.
    conduit::Node &target2_coords = n_mesh["coordsets/target2_coords"];
    target2_coords["type"] = "explicit";
    target2_coords["values/x"].set(conduit::DataType::float64(x.number_of_elements()));
    target2_coords["values/y"].set(conduit::DataType::float64(y.number_of_elements()));
    auto xp = target2_coords["values/x"].as_float64_ptr();
    auto yp = target2_coords["values/y"].as_float64_ptr();

    const double A = M_PI / 16.;
    const double sinA = sin(A);
    const double cosA = cos(A);
    const double M[2][2] = {{cosA, -sinA}, {sinA, cosA}};
    for(conduit::index_t i = 0; i < x.number_of_elements(); i++)
    {
      xp[i] = M[0][0] * x[i] + M[0][1] * y[i];
      yp[i] = M[1][0] * x[i] + M[1][1] * y[i];
    }

    n_mesh["topologies/target2"].set(n_mesh["topologies/topo"]);
    n_mesh["topologies/target2/coordset"] = "target2_coords";
  }

  static void mapping_target2(conduit::Node &n_dev)
  {
    // Wrap polygonal mesh in views.
    auto srcCoordset = views::make_explicit_coordset<double, 2>::view(n_dev["coordsets/coords"]);
    using SrcCoordsetView = decltype(srcCoordset);

    const conduit::Node &n_srcTopo = n_dev["topologies/topo"];
    auto srcTopo =
      views::make_unstructured_single_shape_topology<views::PolygonShape<std::uint64_t>>::view(
        n_srcTopo);
    using SrcTopologyView = decltype(srcTopo);

    const conduit::Node &n_srcMatset = n_dev["matsets/mat"];
    auto srcMatset = views::make_unibuffer_matset<int, float, MAX_MATERIALS>::view(n_srcMatset);
    using SrcMatsetView = decltype(srcMatset);

    // Wrap target2 mesh in views.
    auto targetCoordset =
      views::make_explicit_coordset<double, 2>::view(n_dev["coordsets/target2_coords"]);
    using TargetCoordsetView = decltype(targetCoordset);

    const conduit::Node &n_targetTopo = n_dev["topologies/target2"];
    auto targetTopo =
      views::make_unstructured_single_shape_topology<views::PolygonShape<std::uint64_t>>::view(
        n_targetTopo);
    using TargetTopologyView = decltype(targetTopo);

    // Make new VFs via mapper.
    constexpr int MAX_VERTS = 16;  // Use a larger MAX_VERTS to handle oct-oct clipping.
    using Mapper = bump::TopologyMapper<ExecSpace,
                                        SrcTopologyView,
                                        SrcCoordsetView,
                                        SrcMatsetView,
                                        TargetTopologyView,
                                        TargetCoordsetView,
                                        MAX_VERTS>;
    Mapper mapper(srcTopo, srcCoordset, srcMatset, targetTopo, targetCoordset);
    conduit::Node n_opts;
    n_opts["source/matsetName"] = "mat";
    n_opts["target/topologyName"] = "target2";
    n_opts["target/matsetName"] = "target2_matset";
    mapper.execute(n_dev, n_opts, n_dev);
  }

  static void mir_target2(conduit::Node &n_dev)
  {
    // Wrap target2 mesh in views.
    auto coordsetView =
      views::make_explicit_coordset<double, 2>::view(n_dev["coordsets/target2_coords"]);
    using CoordsetView = decltype(coordsetView);

    const conduit::Node &n_targetTopo = n_dev["topologies/target2"];
    auto topologyView =
      views::make_unstructured_single_shape_topology<views::PolygonShape<std::uint64_t>>::view(
        n_targetTopo);
    using TopologyView = decltype(topologyView);

    const conduit::Node &n_targetMatset = n_dev["matsets/target2_matset"];
    auto matsetView = views::make_unibuffer_matset<int, float, MAX_MATERIALS>::view(n_targetMatset);
    using MatsetView = decltype(matsetView);

    // Do MIR (into new node)
    conduit::Node n_mir;
    using MIR = axom::mir::EquiZAlgorithm<ExecSpace, TopologyView, CoordsetView, MatsetView>;
    MIR m(topologyView, coordsetView, matsetView);
    conduit::Node options;
    options["matset"] = "target2_matset";
    options["matsetName"] = "mir_matset";
    m.execute(n_dev, options, n_mir);

    // Move the MIR output to the n_dev node.
    n_dev["coordsets/mir_coords"].move(n_mir["coordsets/target2_coords"]);
    n_dev["topologies/mir"].move(n_mir["topologies/target2"]);
    n_dev["topologies/mir/coordset"] = "mir_coords";
    n_dev["matsets/mir_matset"].move(n_mir["matsets/mir_matset"]);
    n_dev["matsets/mir_matset/topology"] = "mir";
    n_dev["fields/originalElements"].move(n_mir["fields/originalElements"]);
    n_dev["fields/originalElements/topology"] = "mir";
  }

  static int countBadMaterialZones(const conduit::Node &matset, double eps = 1.e-4)
  {
    const auto volume_fractions = utils::make_array_view<float>(matset["volume_fractions"]);
    //const auto material_ids = utils::make_array_view<int>(matset["material_ids"]);
    const auto indices = utils::make_array_view<int>(matset["indices"]);
    const auto sizes = utils::make_array_view<int>(matset["sizes"]);
    const auto offsets = utils::make_array_view<int>(matset["offsets"]);

    const int nzones = sizes.size();
    int badZones = 0;
    for(int zi = 0; zi < nzones; zi++)
    {
      int matsThisZone = sizes[zi];
      int offset = offsets[zi];

      // What is the total VF for the zone?
      double vfSum = 0.;
      for(int m = 0; m < matsThisZone; m++)
      {
        const int index = indices[offset + m];
        vfSum += volume_fractions[index];
      }

      if(fabs(1. - vfSum) > eps)
      {
        badZones++;
      }
    }
    return badZones;
  }
};

//------------------------------------------------------------------------------
template <typename ExecSpace>
void test_equiz_uniform_unibuffer()
{
  {
    const bool selectedZones = false;
    const bool cleanMats = false;
    braid2d_mat_test<ExecSpace>("uniform",
                                "unibuffer",
                                "equiz_uniform_unibuffer",
                                1,
                                selectedZones,
                                cleanMats);
    braid2d_mat_test<ExecSpace>("uniform",
                                "unibuffer",
                                "equiz_uniform_unibuffer",
                                2,
                                selectedZones,
                                cleanMats);
  }
  {
    const bool selectedZones = true;
    const bool cleanMats = false;
    braid2d_mat_test<ExecSpace>("uniform",
                                "unibuffer",
                                "equiz_uniform_unibuffer_sel",
                                1,
                                selectedZones,
                                cleanMats);
  }
  {
    const bool selectedZones = false;
    const bool cleanMats = true;
    braid2d_mat_test<ExecSpace>("uniform",
                                "unibuffer",
                                "equiz_uniform_unibuffer_clean",
                                1,
                                selectedZones,
                                cleanMats);
  }
  {
    const bool selectedZones = true;
    const bool cleanMats = true;
    braid2d_mat_test<ExecSpace>("uniform",
                                "unibuffer",
                                "equiz_uniform_unibuffer_sel_clean",
                                1,
                                selectedZones,
                                cleanMats);
  }
}

//------------------------------------------------------------------------------
TEST(mir_equiz, equiz_uniform_unibuffer_seq)
{
  AXOM_ANNOTATE_SCOPE("equiz_uniform_unibuffer_seq");
  test_equiz_uniform_unibuffer<seq_exec>();
}

#if defined(AXOM_USE_OPENMP)
TEST(mir_equiz, equiz_uniform_unibuffer_omp)
{
  AXOM_ANNOTATE_SCOPE("equiz_uniform_unibuffer_omp");
  test_equiz_uniform_unibuffer<omp_exec>();
}
#endif

#if defined(AXOM_USE_CUDA)
TEST(mir_equiz, equiz_uniform_unibuffer_cuda)
{
  AXOM_ANNOTATE_SCOPE("equiz_uniform_unibuffer_cuda");
  test_equiz_uniform_unibuffer<cuda_exec>();
}
#endif

#if defined(AXOM_USE_HIP)
TEST(mir_equiz, equiz_uniform_unibuffer_hip)
{
  AXOM_ANNOTATE_SCOPE("equiz_uniform_unibuffer_hip");
  test_equiz_uniform_unibuffer<hip_exec>();
}
#endif

//------------------------------------------------------------------------------
TEST(mir_equiz, equiz_polygonal_unibuffer_seq)
{
  AXOM_ANNOTATE_SCOPE("equiz_polygonal_unibuffer_seq");
  test_Polygonal_MIR<seq_exec>::test("equiz_polygonal_unibuffer");
}

#if defined(AXOM_USE_OPENMP)
TEST(mir_equiz, equiz_polygonal_unibuffer_omp)
{
  AXOM_ANNOTATE_SCOPE("equiz_polygonal_unibuffer_omp");
  test_Polygonal_MIR<omp_exec>::test("equiz_polygonal_unibuffer");
}
#endif

#if defined(AXOM_USE_CUDA)
TEST(mir_equiz, equiz_polygonal_unibuffer_cuda)
{
  AXOM_ANNOTATE_SCOPE("equiz_polygonal_unibuffer_cuda");
  test_Polygonal_MIR<cuda_exec>::test("equiz_polygonal_unibuffer");
}
#endif

#if defined(AXOM_USE_HIP)
TEST(mir_equiz, equiz_polygonal_unibuffer_hip)
{
  AXOM_ANNOTATE_SCOPE("equiz_polygonal_unibuffer_hip");
  test_Polygonal_MIR<hip_exec>::test("equiz_polygonal_unibuffer");
}
#endif

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
