// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/mir.hpp"
#include "axom/mir/tests/mir_testing_data_helpers.hpp"

#include <conduit/conduit_relay_io_blueprint.hpp>
#include <cmath>
#include <cstdlib>

namespace bputils = axom::mir::utilities::blueprint;

//------------------------------------------------------------------------------

// Uncomment to generate baselines
//#define AXOM_TESTING_GENERATE_BASELINES

// Uncomment to save visualization files for debugging (when making baselines)
//#define AXOM_TESTING_SAVE_VISUALIZATION

#include "axom/mir/tests/mir_testing_helpers.hpp"

std::string baselineDirectory()
{
  return pjoin(dataDirectory(), "mir", "regression", "mir_topology_mapper");
}

//------------------------------------------------------------------------------
/*!

coarse

%12-------%13-------%14-------%15
|         |         |         |
|         |         |         |
|         |         |         |
|         |         |         |
|         |         |         |
%8--------%9--------%10-------%11
|         |         |         |
|         |         |         |
|         |         |         |
|         |         |         |
|         |         |         |
%4--------%5--------%6--------%7
|         |         |         |
|         |         |         |
|         |         |         |
|         |         |         |
|         |         |         |
%---------%---------%---------%
0         1         2         3


postmir

%18-------%19-------%20-------%21
|\   (3)  |         |         |
|  \      |         |         |
|    \    |         |         |
|      \  |         |         |
|    (2) \|     (3) |    (3)  |
%14-------%15-------%16-------%17
|\   (2)  | \   (3) |    (3)  |
|  \      |   \     |         |
|    \    | (2)/%- -%12- - - -%13
| (1)  \  |  /  11  |         |
|        \|/   (2)  |    (2)  |
%7--------%8--------%9--------%10
|\   (1)  | \  (2)  |    (2)  |
|  \      |   \     |         |
|    \    | (1)/%- -%5 - - - -%6
|  (0) \  |  /  4   |         |
|        \|/   (1)  |    (1)  |
%---------%---------%---------%
0         1         2         3


fine - refines coarse with equal sized quads.

If fine has 2x2 refinement, it looks like this:

%----%----%----%----%----%----%
|2:.5|    |    |    |    |    |
|3:.5| 3  | 3  | 3  | 3  | 3  |
%----%----%----%----%----%----%
|    |2:.5|    |    |    |    |
| 2  |3:.5| 3  | 3  | 3  | 3  |
%----%----%----%----%----%----%
|1:.5|    |2:.5|    |    |    |
|2:.5| 2  |3:.5| 3  | 3  | 3  |
%----%----%----%----%----%----%
|    |1:.5|    |    |    |    |
| 1  |2:.5| 2  | 2  | 2  | 2  |
%----%----%----%----%----%----%
|0:.5|    |1:.5|    |    |    |
|1:.5| 1  |2:.5| 2  | 2  | 2  |
%----%----%----%----%----%----%
|    |0:.5|    |    |    |    |
| 0  |1:.5| 1  | 1  | 1  | 1  |
%----%----%----%----%----%----%

*/
const char *yaml = R"(
coordsets:
  coarse_coords:
    type: explicit
    values:
      x: [0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3.]
      y: [0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2., 3., 3., 3., 3.]
  postmir_coords:
    type: explicit
    values:
      x: [0., 1., 2., 3., 1.5, 2., 3., 0., 1., 2., 3., 1.5, 2., 3., 0., 1., 2., 3., 0., 1., 2., 3]
      y: [0., 0., 0., 0., 0.5, 0.5, 0.5, 1., 1., 1., 1., 1.5, 1.5, 1.5, 2., 2., 2., 2., 3., 3., 3., 3]
topologies:
  coarse:
    type: unstructured
    coordset: coarse_coords
    elements:
      shape: quad
      connectivity: [0,1,5,4, 1,2,6,5, 2,3,7,6, 4,5,9,8, 5,6,10,9, 6,7,11,10, 8,9,13,12, 9,10,14,13, 10,11,15,14]
      sizes: [4,4,4, 4,4,4, 4,4,4]
      offsets: [0,4,8, 12,16,20, 24,28,32]
  postmir:
    type: unstructured
    coordset: postmir_coords
    elements:
      shape: mixed
      shapes: [3,3,3,4,4,4,4, 3,3,3,4,4,4,4, 3,3,4,4]
      connectivity: [0,1,7, 1,8,7, 1,4,8, 1,2,5,4, 4,5,9,8, 2,3,6,5, 5,6,10,9, 7,8,14, 8,15,14, 8,11,15, 8,9,12,11, 11,12,16,15, 9,10,13,12, 12,13,17,16, 14,15,18, 15,19,18, 15,16,20,19, 16,17,21,20]
      sizes: [3,3,3,4,4,4,4, 3,3,3,4,4,4,4, 3,3,4,4]
      offsets: [0, 3, 6, 9, 13, 17, 21, 25, 28, 31, 34, 38, 42, 46, 50, 53, 56, 60]
      shape_map:
        quad: 4
        tri: 3
matsets:
  coarse_matset:
   topology: coarse
   material_map:
     mat0: 0
     mat1: 1
     mat2: 2
     mat3: 3
   material_ids: [0,1, 1,2, 1,2, 1,2, 2,3, 2,3, 2,3, 3, 3]
   volume_fractions: [0.5,0.5, 0.625,0.375, 0.5,0.5, 0.5,0.5, 0.625,0.375, 0.5,0.5, 0.5,0.5, 1., 1.]
   indices: [0,1, 2,3, 4,5, 6,7, 8,9, 10,11, 12,13, 14, 15]
   sizes: [2, 2, 2, 2, 2, 2, 2, 1, 1]
   offsets: [0, 2, 4, 6, 8, 10, 12, 13, 14]
  postmir_matset:
   topology: postmir
   material_map:
     mat0: 0
     mat1: 1
     mat2: 2
     mat3: 3
   material_ids: [0, 1, 1, 1, 2, 1, 2, 1, 2, 2, 2, 3, 2, 3, 2, 3, 3, 3]
   volume_fractions: [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.]
   indices: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
   sizes: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
   offsets: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
)";

//------------------------------------------------------------------------------
/*!
 * \brief Make a fine mesh coordset and topology that represents a refined mesh.
 *
 * \param n_mesh The Conduit node that will contain the mesh.
 * \param coordsetName The name of the new fine mesh's coordset.
 * \param topoName The name of the new fine mesh topology.
 * \param nx The number of nodes in the X direction.
 * \param ny The number of nodes in the Y direction.
 * \param refinement The number of refinements to make from the coarse to fine levels.
 */
void make_fine(conduit::Node &n_mesh,
               const std::string &coordsetName,
               const std::string &topoName,
               int nx,
               int ny,
               int refinement)
{
  int nxr = (nx - 1) * refinement + 1;
  int nyr = (ny - 1) * refinement + 1;

  const int nnodes = nxr * nyr;
  std::vector<double> xc, yc;
  xc.reserve(nnodes);
  yc.reserve(nnodes);

  constexpr double y1 = 3.;
  constexpr double x1 = 3.;

  for(int j = 0; j < nyr; j++)
  {
    double tj = double(j) / double(nyr - 1);
    double y = tj * y1;
    for(int i = 0; i < nxr; i++)
    {
      double ti = double(i) / double(nxr - 1);
      double x = ti * x1;
      xc.push_back(x);
      yc.push_back(y);
    }
  }

  std::vector<int> conn, sizes, offsets;
  int cxr = nxr - 1;
  int cyr = nyr - 1;
  int offset = 0;
  for(int j = 0; j < cyr; j++)
  {
    for(int i = 0; i < cxr; i++)
    {
      conn.push_back(j * nxr + i);
      conn.push_back(j * nxr + i + 1);
      conn.push_back((j + 1) * nxr + i + 1);
      conn.push_back((j + 1) * nxr + i);
      sizes.push_back(4);
      offsets.push_back(offset);
      offset += 4;
    }
  }

  conduit::Node &n_coordset = n_mesh["coordsets/" + coordsetName];
  n_coordset["type"] = "explicit";
  n_coordset["values/x"].set(xc);
  n_coordset["values/y"].set(yc);

  conduit::Node &n_topo = n_mesh["topologies/" + topoName];
  n_topo["type"] = "unstructured";
  n_topo["coordset"] = coordsetName;
  n_topo["elements/shape"] = "quad";
  n_topo["elements/connectivity"].set(conn);
  n_topo["elements/sizes"].set(sizes);
  n_topo["elements/offsets"].set(offsets);
}

//------------------------------------------------------------------------------
/*!
 * \brief Tests the TopologyMapper.
 *
 * \note TODO: Test selected zone lists on source and target.
 */
template <typename ExecSpace>
class test_TopologyMapper
{
public:
  static void test2D()
  {
    // Make the 2D input mesh.
    conduit::Node n_mesh;
    initialize(n_mesh);

    // host->device
    conduit::Node n_dev;
    axom::mir::utilities::blueprint::copy<ExecSpace>(n_dev, n_mesh);

    mapping2D(n_dev);

    // device->host
    conduit::Node hostResult;
    bputils::copy<seq_exec>(hostResult, n_dev);

#if defined(AXOM_TESTING_SAVE_VISUALIZATION)
    conduit::relay::io::blueprint::save_mesh(hostResult, "test2D", "hdf5");
    conduit::relay::io::save(hostResult, "test2D.yaml", "yaml");
#endif

    // Handle baseline comparison.
    const auto paths = baselinePaths<ExecSpace>();
    std::string baselineName(yamlRoot("test2D"));
#if defined(AXOM_TESTING_GENERATE_BASELINES)
    saveBaseline(paths, baselineName, hostResult);
#else
    EXPECT_TRUE(compareBaseline(paths, baselineName, hostResult));
#endif
  }

  static void test3D()
  {
    // Make the 2D input mesh.
    conduit::Node n_mesh;
    initialize(n_mesh);

    // host->device
    conduit::Node n_dev;
    axom::mir::utilities::blueprint::copy<ExecSpace>(n_dev, n_mesh);

    // Extrude relevant meshes into 3D.
    extrude(n_dev);

    mapping3D(n_dev);

    // device->host
    conduit::Node hostResult;
    bputils::copy<seq_exec>(hostResult, n_dev);

#if defined(AXOM_TESTING_SAVE_VISUALIZATION)
    conduit::relay::io::blueprint::save_mesh(hostResult, "test3D", "hdf5");
    conduit::relay::io::save(hostResult, "test3D.yaml", "yaml");
#endif

    // Handle baseline comparison.
    const auto paths = baselinePaths<ExecSpace>();
    std::string baselineName(yamlRoot("test3D"));
#if defined(AXOM_TESTING_GENERATE_BASELINES)
    saveBaseline(paths, baselineName, hostResult);
#else
    EXPECT_TRUE(compareBaseline(paths, baselineName, hostResult));
#endif
  }

private:
  static constexpr int refinement = 4;

  static void initialize(conduit::Node &n_mesh)
  {
    // Make the 2D input mesh.
    n_mesh.parse(yaml);

    // Make a
    int nx = 4;
    int ny = 4;
    make_fine(n_mesh, "fine_coords", "fine", nx, ny, refinement);

    // workaround (make shape map)
    n_mesh["topologies/postmir/elements/shape_map/quad"] = 4;
    n_mesh["topologies/postmir/elements/shape_map/tri"] = 3;
  }

  /*!
   * \brief Extrude postmir and fine meshes.
   *
   * \param n_dev The Conduit node that contains the input meshes and will contain
   *              the output meshes. The data needs to be in the right memory for
   *              the ExecutionSpace.
   */
  static void extrude(conduit::Node &n_dev)
  {
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    // Wrap coarse/post_mir mesh in views.
    using SrcCoordsetView = axom::mir::views::ExplicitCoordsetView<double, 2>;
    using SrcTopologyView = axom::mir::views::UnstructuredTopologyMixedShapeView<conduit::index_t>;
    SrcCoordsetView srcCoordset(
      bputils::make_array_view<double>(n_dev["coordsets/postmir_coords/values/x"]),
      bputils::make_array_view<double>(n_dev["coordsets/postmir_coords/values/y"]));

    axom::Array<axom::IndexType> shapeValues, shapeIds;
    const conduit::Node &n_srcTopo = n_dev["topologies/postmir"];
    auto shapeMap = axom::mir::views::buildShapeMap(n_srcTopo, shapeValues, shapeIds, allocatorID);
    SrcTopologyView srcTopo(
      bputils::make_array_view<conduit::index_t>(n_srcTopo["elements/connectivity"]),
      bputils::make_array_view<conduit::index_t>(n_srcTopo["elements/shapes"]),
      bputils::make_array_view<conduit::index_t>(n_srcTopo["elements/sizes"]),
      bputils::make_array_view<conduit::index_t>(n_srcTopo["elements/offsets"]),
      shapeMap);

    // Wrap fine mesh in views.
    using TargetCoordsetView = axom::mir::views::ExplicitCoordsetView<double, 2>;
    using TargetShapeType = axom::mir::views::QuadShape<int>;
    using TargetTopologyView = axom::mir::views::UnstructuredTopologySingleShapeView<TargetShapeType>;

    TargetCoordsetView targetCoordset(
      bputils::make_array_view<double>(n_dev["coordsets/fine_coords/values/x"]),
      bputils::make_array_view<double>(n_dev["coordsets/fine_coords/values/y"]));
    const conduit::Node &n_targetTopo = n_dev["topologies/fine"];
    TargetTopologyView targetTopo(
      bputils::make_array_view<int>(n_targetTopo["elements/connectivity"]),
      bputils::make_array_view<int>(n_targetTopo["elements/sizes"]),
      bputils::make_array_view<int>(n_targetTopo["elements/offsets"]));

    //_mir_utilities_extrudemesh_begin
    // Make new VFs via mapper.
    const int coarseNodesInZ = 4;
    using SrcExtruder = bputils::ExtrudeMesh<ExecSpace, SrcTopologyView, SrcCoordsetView>;
    SrcExtruder srcExt(srcTopo, srcCoordset);
    conduit::Node n_opts;
    n_opts["nz"] = coarseNodesInZ;
    n_opts["z0"] = 0.;
    n_opts["z1"] = 3.;
    n_opts["topologyName"] = "postmir";
    n_opts["outputTopologyName"] = "epm";  // epm = "Extruded Post MIR"
    n_opts["outputCoordsetName"] = "epm_coords";
    n_opts["outputMatsetName"] = "epm_matset";
    srcExt.execute(n_dev, n_opts, n_dev);
    //_mir_utilities_extrudemesh_end

    using TargetExtruder = bputils::ExtrudeMesh<ExecSpace, TargetTopologyView, TargetCoordsetView>;
    TargetExtruder targetExt(targetTopo, targetCoordset);
    int fineNodesInZ = (coarseNodesInZ - 1) * refinement + 1;
    conduit::Node n_opts2;
    n_opts2["nz"] = fineNodesInZ;
    n_opts2["z0"] = 0.;
    n_opts2["z1"] = 3.;
    n_opts2["topologyName"] = "fine";
    n_opts2["outputTopologyName"] = "efm";  // epm = "Extruded Fine Mesh"
    n_opts2["outputCoordsetName"] = "efm_coords";
    targetExt.execute(n_dev, n_opts2, n_dev);
  }

  static void mapping2D(conduit::Node &n_dev)
  {
    namespace bputils = axom::mir::utilities::blueprint;

    // Wrap coarse/post_mir mesh in views.
    using SrcCoordsetView = axom::mir::views::ExplicitCoordsetView<double, 2>;
    using SrcTopologyView = axom::mir::views::UnstructuredTopologyMixedShapeView<conduit::index_t>;
    SrcCoordsetView srcCoordset(
      bputils::make_array_view<double>(n_dev["coordsets/postmir_coords/values/x"]),
      bputils::make_array_view<double>(n_dev["coordsets/postmir_coords/values/y"]));

    axom::Array<axom::IndexType> shapeValues, shapeIds;
    const conduit::Node &n_srcTopo = n_dev["topologies/postmir"];
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    auto shapeMap = axom::mir::views::buildShapeMap(n_srcTopo, shapeValues, shapeIds, allocatorID);
    SrcTopologyView srcTopo(
      bputils::make_array_view<conduit::index_t>(n_srcTopo["elements/connectivity"]),
      bputils::make_array_view<conduit::index_t>(n_srcTopo["elements/shapes"]),
      bputils::make_array_view<conduit::index_t>(n_srcTopo["elements/sizes"]),
      bputils::make_array_view<conduit::index_t>(n_srcTopo["elements/offsets"]),
      shapeMap);

    const conduit::Node &n_srcMatset = n_dev["matsets/postmir_matset"];
    auto srcMatset =
      axom::mir::views::make_unibuffer_matset<std::int64_t, double, 4>::view(n_srcMatset);
    using SrcMatsetView = decltype(srcMatset);

    // Wrap fine mesh in views.
    using TargetCoordsetView = axom::mir::views::ExplicitCoordsetView<double, 2>;
    using TargetShapeType = axom::mir::views::QuadShape<int>;
    using TargetTopologyView = axom::mir::views::UnstructuredTopologySingleShapeView<TargetShapeType>;

    TargetCoordsetView targetCoordset(
      bputils::make_array_view<double>(n_dev["coordsets/fine_coords/values/x"]),
      bputils::make_array_view<double>(n_dev["coordsets/fine_coords/values/y"]));
    const conduit::Node &n_targetTopo = n_dev["topologies/fine"];
    TargetTopologyView targetTopo(
      bputils::make_array_view<int>(n_targetTopo["elements/connectivity"]),
      bputils::make_array_view<int>(n_targetTopo["elements/sizes"]),
      bputils::make_array_view<int>(n_targetTopo["elements/offsets"]));

    // _mir_utilities_topologymapper_begin
    // Make new VFs via mapper.
    using Mapper =
      bputils::TopologyMapper<ExecSpace, SrcTopologyView, SrcCoordsetView, SrcMatsetView, TargetTopologyView, TargetCoordsetView>;
    Mapper mapper(srcTopo, srcCoordset, srcMatset, targetTopo, targetCoordset);
    conduit::Node n_opts;
    n_opts["source/matsetName"] = "postmir_matset";
    n_opts["target/topologyName"] = "fine";
    n_opts["target/matsetName"] = "fine_matset";
    mapper.execute(n_dev, n_opts, n_dev);
    // _mir_utilities_topologymapper_end
  }

  static void mapping3D(conduit::Node &n_dev)
  {
    // Wrap coarse/post_mir mesh in views.
    using SrcCoordsetView = axom::mir::views::ExplicitCoordsetView<double, 3>;
    using SrcTopologyView = axom::mir::views::UnstructuredTopologyMixedShapeView<conduit::index_t>;
    SrcCoordsetView srcCoordset(
      bputils::make_array_view<double>(n_dev["coordsets/epm_coords/values/x"]),
      bputils::make_array_view<double>(n_dev["coordsets/epm_coords/values/y"]),
      bputils::make_array_view<double>(n_dev["coordsets/epm_coords/values/z"]));

    axom::Array<axom::IndexType> shapeValues, shapeIds;
    const conduit::Node &n_srcTopo = n_dev["topologies/epm"];
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    auto shapeMap = axom::mir::views::buildShapeMap(n_srcTopo, shapeValues, shapeIds, allocatorID);
    SrcTopologyView srcTopo(
      bputils::make_array_view<conduit::index_t>(n_srcTopo["elements/connectivity"]),
      bputils::make_array_view<conduit::index_t>(n_srcTopo["elements/shapes"]),
      bputils::make_array_view<conduit::index_t>(n_srcTopo["elements/sizes"]),
      bputils::make_array_view<conduit::index_t>(n_srcTopo["elements/offsets"]),
      shapeMap);

    const conduit::Node &n_srcMatset = n_dev["matsets/epm_matset"];
    auto srcMatset =
      axom::mir::views::make_unibuffer_matset<std::int64_t, double, 4>::view(n_srcMatset);
    using SrcMatsetView = decltype(srcMatset);

    // Wrap fine mesh in views.
    using TargetCoordsetView = axom::mir::views::ExplicitCoordsetView<double, 3>;
    using TargetShapeType = axom::mir::views::HexShape<int>;
    using TargetTopologyView = axom::mir::views::UnstructuredTopologySingleShapeView<TargetShapeType>;

    TargetCoordsetView targetCoordset(
      bputils::make_array_view<double>(n_dev["coordsets/efm_coords/values/x"]),
      bputils::make_array_view<double>(n_dev["coordsets/efm_coords/values/y"]),
      bputils::make_array_view<double>(n_dev["coordsets/efm_coords/values/z"]));
    const conduit::Node &n_targetTopo = n_dev["topologies/efm"];
    TargetTopologyView targetTopo(
      bputils::make_array_view<int>(n_targetTopo["elements/connectivity"]),
      bputils::make_array_view<int>(n_targetTopo["elements/sizes"]),
      bputils::make_array_view<int>(n_targetTopo["elements/offsets"]));

    // Make new VFs via mapper.
    using Mapper =
      bputils::TopologyMapper<ExecSpace, SrcTopologyView, SrcCoordsetView, SrcMatsetView, TargetTopologyView, TargetCoordsetView>;
    Mapper mapper(srcTopo, srcCoordset, srcMatset, targetTopo, targetCoordset);
    conduit::Node n_opts;
    n_opts["source/matsetName"] = "epm_matset";
    n_opts["target/topologyName"] = "efm";
    n_opts["target/matsetName"] = "efm_matset";
    mapper.execute(n_dev, n_opts, n_dev);
  }
};

//------------------------------------------------------------------------------
TEST(mir_topology_mapper, TopologyMapper_2D_seq)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_2D_seq");
  test_TopologyMapper<seq_exec>::test2D();
}
#if defined(AXOM_USE_OPENMP)
TEST(mir_topology_mapper, TopologyMapper_2D_omp)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_2D_omp");
  test_TopologyMapper<omp_exec>::test2D();
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(mir_topology_mapper, TopologyMapper_2D_cuda)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_2D_cuda");
  test_TopologyMapper<cuda_exec>::test2D();
}
#endif
#if defined(AXOM_USE_HIP)
TEST(mir_topology_mapper, TopologyMapper_2D_hip)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_2D_hip");
  test_TopologyMapper<hip_exec>::test2D();
}
#endif

//------------------------------------------------------------------------------
TEST(mir_topology_mapper, TopologyMapper_3D_seq)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_3D_seq");
  test_TopologyMapper<seq_exec>::test3D();
}
#if defined(AXOM_USE_OPENMP)
TEST(mir_topology_mapper, TopologyMapper_3D_omp)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_3D_omp");
  test_TopologyMapper<omp_exec>::test3D();
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(mir_topology_mapper, TopologyMapper_3D_cuda)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_3D_cuda");
  test_TopologyMapper<cuda_exec>::test3D();
}
#endif
#if defined(AXOM_USE_HIP)
TEST(mir_topology_mapper, TopologyMapper_3D_hip)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_3D_hip");
  test_TopologyMapper<hip_exec>::test3D();
}
#endif

//------------------------------------------------------------------------------
void conduit_debug_err_handler(const std::string &s1, const std::string &s2, int i1)
{
  std::cout << "s1=" << s1 << ", s2=" << s2 << ", i1=" << i1 << std::endl;
  // This is on purpose.
  while(1)
    ;
}

//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  // Define command line options.
  axom::CLI::App app;
#if defined(AXOM_USE_CALIPER)
  std::string annotationMode("none");
  app.add_option("--caliper", annotationMode)
    ->description(
      "caliper annotation mode. Valid options include 'none' and 'report'. "
      "Use 'help' to see full list.")
    ->capture_default_str()
    ->check(axom::utilities::ValidCaliperMode);
#endif
  bool handlerEnabled = false;
  app.add_flag("--handler", handlerEnabled, "Enable Conduit handler.");

  // Parse command line options.
  try
  {
    app.parse(argc, argv);

#if defined(AXOM_USE_CALIPER)
    axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(annotationMode);
#endif

    axom::slic::SimpleLogger logger;  // create & initialize test logger,
    if(handlerEnabled)
    {
      conduit::utils::set_error_handler(conduit_debug_err_handler);
    }

    result = RUN_ALL_TESTS();
  }
  catch(axom::CLI::CallForHelp &e)
  {
    std::cout << app.help() << std::endl;
    result = 0;
  }
  catch(axom::CLI::ParseError &e)
  {
    // Handle other parsing errors
    std::cerr << e.what() << std::endl;
    result = app.exit(e);
  }

  return result;
}
