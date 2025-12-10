// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/bump.hpp"
#include "axom/bump/tests/blueprint_testing_data_helpers.hpp"
#include "axom/bump/tests/blueprint_testing_helpers.hpp"

#include <conduit/conduit_relay_io_blueprint.hpp>
#include <cmath>
#include <cstdlib>

namespace bump = axom::bump;
namespace utils = axom::bump::utilities;
namespace views = axom::bump::views;

std::string baselineDirectory()
{
  return pjoin(dataDirectory(), "bump", "regression", "bump_topology_mapper");
}

//------------------------------------------------------------------------------
// Global test application object.
axom::blueprint::testing::TestApplication TestApp;

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
 * \param extents The mesh extents {x0, x1, y0, y1}.
 * \param nx The number of nodes in the X direction.
 * \param ny The number of nodes in the Y direction.
 * \param refinement The number of refinements to make from the coarse to fine levels.
 */
void make_fine(conduit::Node &n_mesh,
               const std::string &coordsetName,
               const std::string &topoName,
               const double *extents,
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

  for(int j = 0; j < nyr; j++)
  {
    double tj = double(j) / double(nyr - 1);
    double y = extents[2] + tj * (extents[3] - extents[2]);
    for(int i = 0; i < nxr; i++)
    {
      double ti = double(i) / double(nxr - 1);
      double x = extents[0] + ti * (extents[1] - extents[0]);
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
    utils::copy<ExecSpace>(n_dev, n_mesh);

    mapping2D(n_dev);

    // device->host
    conduit::Node hostResult;
    utils::copy<seq_exec>(hostResult, n_dev);

    TestApp.saveVisualization("test2D", hostResult);

    // Handle baseline comparison.
    EXPECT_TRUE(TestApp.test<ExecSpace>("test2D", hostResult));
  }

  static void test3D()
  {
    // Make the 2D input mesh.
    conduit::Node n_mesh;
    initialize(n_mesh);

    // host->device
    conduit::Node n_dev;
    utils::copy<ExecSpace>(n_dev, n_mesh);

    // Extrude relevant meshes into 3D.
    extrude(n_dev);

    mapping3D(n_dev);

    // device->host
    conduit::Node hostResult;
    utils::copy<seq_exec>(hostResult, n_dev);

    TestApp.saveVisualization("test3D", hostResult);

    // Handle baseline comparison.
    EXPECT_TRUE(TestApp.test<ExecSpace>("test3D", hostResult));
  }

  static void testPolyhedral()
  {
    // Make the 2D input mesh.
    conduit::Node n_mesh;
    initialize(n_mesh);

    // host->device
    conduit::Node n_dev;
    utils::copy<ExecSpace>(n_dev, n_mesh);

    // Extrude relevant meshes into 3D.
    extrude(n_dev);

    // Make the source mesh polyhedral
    makePolyhedral(n_dev);

    mappingPolyhedral(n_dev);

    // device->host
    conduit::Node hostResult;
    utils::copy<seq_exec>(hostResult, n_dev);

    TestApp.saveVisualization("testPH", hostResult);

    // Handle baseline comparison.
    EXPECT_TRUE(TestApp.test<ExecSpace>("testPH", hostResult));
  }

private:
  static constexpr int refinement = 4;

  static void initialize(conduit::Node &n_mesh)
  {
    // Make the 2D input mesh.
    n_mesh.parse(yaml);

    // Make a fine mesh
    const double extents[] = {0., 3., 0., 3.};
    const int nx = 4;
    const int ny = 4;
    make_fine(n_mesh, "fine_coords", "fine", extents, nx, ny, refinement);

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
    auto srcCoordset =
      views::make_explicit_coordset<double, 2>::view(n_dev["coordsets/postmir_coords"]);
    using SrcCoordsetView = decltype(srcCoordset);

    using SrcTopologyView = views::UnstructuredTopologyMixedShapeView<conduit::index_t>;
    axom::Array<axom::IndexType> shapeValues, shapeIds;
    const conduit::Node &n_srcTopo = n_dev["topologies/postmir"];
    auto shapeMap = views::buildShapeMap(n_srcTopo, shapeValues, shapeIds, allocatorID);
    SrcTopologyView srcTopo(
      utils::make_array_view<conduit::index_t>(n_srcTopo["elements/connectivity"]),
      utils::make_array_view<conduit::index_t>(n_srcTopo["elements/shapes"]),
      utils::make_array_view<conduit::index_t>(n_srcTopo["elements/sizes"]),
      utils::make_array_view<conduit::index_t>(n_srcTopo["elements/offsets"]),
      shapeMap);

    // Wrap fine mesh in views.
    auto targetCoordset =
      views::make_explicit_coordset<double, 2>::view(n_dev["coordsets/fine_coords"]);
    using TargetCoordsetView = decltype(targetCoordset);

    const conduit::Node &n_targetTopo = n_dev["topologies/fine"];
    using TargetShapeType = views::QuadShape<int>;
    auto targetTopo =
      views::make_unstructured_single_shape_topology<TargetShapeType>::view(n_targetTopo);
    using TargetTopologyView = decltype(targetTopo);

    //_bump_utilities_extrudemesh_begin
    // Make new VFs via mapper.
    const int coarseNodesInZ = 4;
    using SrcExtruder = bump::ExtrudeMesh<ExecSpace, SrcTopologyView, SrcCoordsetView>;
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
    //_bump_utilities_extrudemesh_end

    using TargetExtruder = bump::ExtrudeMesh<ExecSpace, TargetTopologyView, TargetCoordsetView>;
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
    // Wrap coarse/post_mir mesh in views.
    auto srcCoordset =
      views::make_explicit_coordset<double, 2>::view(n_dev["coordsets/postmir_coords"]);
    using SrcCoordsetView = decltype(srcCoordset);

    using SrcTopologyView = views::UnstructuredTopologyMixedShapeView<conduit::index_t>;
    axom::Array<axom::IndexType> shapeValues, shapeIds;
    const conduit::Node &n_srcTopo = n_dev["topologies/postmir"];
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    auto shapeMap = views::buildShapeMap(n_srcTopo, shapeValues, shapeIds, allocatorID);
    SrcTopologyView srcTopo(
      utils::make_array_view<conduit::index_t>(n_srcTopo["elements/connectivity"]),
      utils::make_array_view<conduit::index_t>(n_srcTopo["elements/shapes"]),
      utils::make_array_view<conduit::index_t>(n_srcTopo["elements/sizes"]),
      utils::make_array_view<conduit::index_t>(n_srcTopo["elements/offsets"]),
      shapeMap);

    const conduit::Node &n_srcMatset = n_dev["matsets/postmir_matset"];
    auto srcMatset = views::make_unibuffer_matset<std::int64_t, double, 4>::view(n_srcMatset);
    using SrcMatsetView = decltype(srcMatset);

    // Wrap fine mesh in views.
    auto targetCoordset =
      views::make_explicit_coordset<double, 2>::view(n_dev["coordsets/fine_coords"]);
    using TargetCoordsetView = decltype(targetCoordset);

    const conduit::Node &n_targetTopo = n_dev["topologies/fine"];
    using TargetShapeType = views::QuadShape<int>;
    auto targetTopo =
      views::make_unstructured_single_shape_topology<TargetShapeType>::view(n_targetTopo);
    using TargetTopologyView = decltype(targetTopo);

    // _bump_utilities_topologymapper_begin
    // Make new VFs via mapper.
    using Mapper =
      bump::TopologyMapper<ExecSpace, SrcTopologyView, SrcCoordsetView, SrcMatsetView, TargetTopologyView, TargetCoordsetView>;
    Mapper mapper(srcTopo, srcCoordset, srcMatset, targetTopo, targetCoordset);
    conduit::Node n_opts;
    n_opts["source/matsetName"] = "postmir_matset";
    n_opts["target/topologyName"] = "fine";
    n_opts["target/matsetName"] = "fine_matset";
    mapper.execute(n_dev, n_opts, n_dev);
    // _bump_utilities_topologymapper_end
  }

  static void mapping3D(conduit::Node &n_dev)
  {
    // Wrap coarse/post_mir mesh in views.
    auto srcCoordset =
      views::make_explicit_coordset<double, 3>::view(n_dev["coordsets/epm_coords"]);
    using SrcCoordsetView = decltype(srcCoordset);
    using SrcTopologyView = views::UnstructuredTopologyMixedShapeView<conduit::index_t>;
    axom::Array<axom::IndexType> shapeValues, shapeIds;
    const conduit::Node &n_srcTopo = n_dev["topologies/epm"];
    EXPECT_EQ(n_srcTopo["type"].as_string(), "unstructured");
    EXPECT_EQ(n_srcTopo["elements/shape"].as_string(), "mixed");

    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    auto shapeMap = views::buildShapeMap(n_srcTopo, shapeValues, shapeIds, allocatorID);
    const auto srcConnView =
      utils::make_array_view<conduit::index_t>(n_srcTopo["elements/connectivity"]);
    const auto srcShapesView =
      utils::make_array_view<conduit::index_t>(n_srcTopo["elements/shapes"]);
    const auto srcSizesView = utils::make_array_view<conduit::index_t>(n_srcTopo["elements/sizes"]);
    const auto srcOffsetsView =
      utils::make_array_view<conduit::index_t>(n_srcTopo["elements/offsets"]);
    // Check sizes with the current version of this test.
    EXPECT_TRUE(srcSizesView.size() == 54);
    EXPECT_TRUE(srcSizesView.size() == srcOffsetsView.size());
    EXPECT_TRUE(srcSizesView.size() == srcShapesView.size());

    SrcTopologyView srcTopo(srcConnView, srcShapesView, srcSizesView, srcOffsetsView, shapeMap);

    const conduit::Node &n_srcMatset = n_dev["matsets/epm_matset"];
    auto srcMatset = views::make_unibuffer_matset<std::int64_t, double, 4>::view(n_srcMatset);
    using SrcMatsetView = decltype(srcMatset);

    // Wrap fine mesh in views.
    auto targetCoordset =
      views::make_explicit_coordset<double, 3>::view(n_dev["coordsets/efm_coords"]);
    using TargetCoordsetView = decltype(targetCoordset);

    const conduit::Node &n_targetTopo = n_dev["topologies/efm"];
    using TargetShapeType = views::HexShape<int>;
    auto targetTopo =
      views::make_unstructured_single_shape_topology<TargetShapeType>::view(n_targetTopo);
    using TargetTopologyView = decltype(targetTopo);

    // Make new VFs via mapper.
    using Mapper =
      bump::TopologyMapper<ExecSpace, SrcTopologyView, SrcCoordsetView, SrcMatsetView, TargetTopologyView, TargetCoordsetView>;
    Mapper mapper(srcTopo, srcCoordset, srcMatset, targetTopo, targetCoordset);
    conduit::Node n_opts;
    n_opts["source/matsetName"] = "epm_matset";
    n_opts["target/topologyName"] = "efm";
    n_opts["target/matsetName"] = "efm_matset";
    mapper.execute(n_dev, n_opts, n_dev);
  }

  static void makePolyhedral(conduit::Node &n_dev)
  {
    // Wrap coarse/epm mesh in a view.
    using SrcTopologyView = views::UnstructuredTopologyMixedShapeView<conduit::index_t>;
    axom::Array<axom::IndexType> shapeValues, shapeIds;
    const conduit::Node &n_srcTopo = n_dev["topologies/epm"];
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    auto shapeMap = views::buildShapeMap(n_srcTopo, shapeValues, shapeIds, allocatorID);
    SrcTopologyView srcTopo(
      utils::make_array_view<conduit::index_t>(n_srcTopo["elements/connectivity"]),
      utils::make_array_view<conduit::index_t>(n_srcTopo["elements/shapes"]),
      utils::make_array_view<conduit::index_t>(n_srcTopo["elements/sizes"]),
      utils::make_array_view<conduit::index_t>(n_srcTopo["elements/offsets"]),
      shapeMap);

    // Turn the source mesh "epm" polyhedral and store in phmesh.
    conduit::Node &n_phTopo = n_dev["topologies/phmesh"];
    bump::MakePolyhedralTopology<ExecSpace, SrcTopologyView> mph(srcTopo);
    mph.execute(n_srcTopo, n_phTopo);
    bump::MergePolyhedralFaces<ExecSpace, conduit::index_t>::execute(n_phTopo);

    // Copy epm_matset to phmatset.
    utils::copy<ExecSpace>(n_dev["matsets/ph_matset"], n_dev["matsets/epm_matset"]);
    n_dev["matsets/ph_matset/topology"] = "phmesh";
  }

  static void mappingPolyhedral(conduit::Node &n_dev)
  {
    // Wrap coarse/post_mir mesh in views.
    auto srcCoordset =
      views::make_explicit_coordset<double, 3>::view(n_dev["coordsets/epm_coords"]);
    using SrcCoordsetView = decltype(srcCoordset);

    // Make polyhedral topology view.
    const conduit::Node &n_srcTopo = n_dev["topologies/phmesh"];
    auto srcTopo = views::make_unstructured_polyhedral_topology<conduit::index_t>::view(n_srcTopo);
    using SrcTopologyView = decltype(srcTopo);

    const conduit::Node &n_srcMatset = n_dev["matsets/ph_matset"];
    auto srcMatset = views::make_unibuffer_matset<std::int64_t, double, 4>::view(n_srcMatset);
    using SrcMatsetView = decltype(srcMatset);

    // Wrap fine mesh in views.
    auto targetCoordset =
      views::make_explicit_coordset<double, 3>::view(n_dev["coordsets/efm_coords"]);
    using TargetCoordsetView = decltype(targetCoordset);

    const conduit::Node &n_targetTopo = n_dev["topologies/efm"];
    using TargetShapeType = views::HexShape<int>;
    auto targetTopo =
      views::make_unstructured_single_shape_topology<TargetShapeType>::view(n_targetTopo);
    using TargetTopologyView = decltype(targetTopo);

    // Make new VFs via mapper.
    using Mapper =
      bump::TopologyMapper<ExecSpace, SrcTopologyView, SrcCoordsetView, SrcMatsetView, TargetTopologyView, TargetCoordsetView>;
    Mapper mapper(srcTopo, srcCoordset, srcMatset, targetTopo, targetCoordset);
    conduit::Node n_opts;
    n_opts["source/matsetName"] = "ph_matset";
    n_opts["target/topologyName"] = "efm";
    n_opts["target/matsetName"] = "efm_matset";
    mapper.execute(n_dev, n_opts, n_dev);
  }
};

//------------------------------------------------------------------------------
/*!
 * \brief Tests the TopologyMapper on polygonal geometry.
 */
template <typename ExecSpace>
class test_TopologyMapper_Polygonal
{
public:
  static void test()
  {
    // Make the 2D input mesh.
    conduit::Node n_mesh;
    initialize(n_mesh);

    // host->device
    conduit::Node n_dev;
    utils::copy<ExecSpace>(n_dev, n_mesh);

    mapping_target1(n_dev);
    mapping_target2(n_dev);

    // device->host
    conduit::Node hostResult;
    utils::copy<seq_exec>(hostResult, n_dev);

    EXPECT_EQ(countBadMaterialZones(hostResult["matsets/target1_matset"]), 0);
    EXPECT_EQ(countBadMaterialZones(hostResult["matsets/target2_matset"]), 0);

    TestApp.saveVisualization("test_poly", hostResult);

    // Handle baseline comparison.
    EXPECT_TRUE(TestApp.test<ExecSpace>("test_poly", hostResult));
  }

  static void initialize(conduit::Node &n_mesh)
  {
    // Make polygonal geometry
    const conduit::index_t nlevels = 4;
    const conduit::index_t nz = 1;
    conduit::blueprint::mesh::examples::polytess(nlevels, nz, n_mesh);

    // Make a matset from the level field.
    conduit::Node &n_matset = n_mesh["matsets/mat"];
    n_matset["topology"] = "topo";
    for(int mat = 1; mat <= nlevels; mat++)
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

    make_target1(n_mesh);
    make_target2(n_mesh);
  }

  static void make_target1(conduit::Node &n_mesh)
  {
    // Make a quad mesh
    double extents[] = {-6.32843, 6.32843, -6.32843, 6.32843};
    make_fine(n_mesh, "target1_coords", "target1", extents, 100, 100, 1);
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

  static void mapping_target1(conduit::Node &n_dev)
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
    auto srcMatset = views::make_unibuffer_matset<int, float, 4>::view(n_srcMatset);
    using SrcMatsetView = decltype(srcMatset);

    // Wrap target1 mesh in views.
    auto targetCoordset =
      views::make_explicit_coordset<double, 2>::view(n_dev["coordsets/target1_coords"]);
    using TargetCoordsetView = decltype(targetCoordset);

    const conduit::Node &n_targetTopo = n_dev["topologies/target1"];
    using TargetShapeType = views::QuadShape<int>;
    auto targetTopo =
      views::make_unstructured_single_shape_topology<TargetShapeType>::view(n_targetTopo);
    using TargetTopologyView = decltype(targetTopo);

    // Make new VFs via mapper.
    using Mapper =
      bump::TopologyMapper<ExecSpace, SrcTopologyView, SrcCoordsetView, SrcMatsetView, TargetTopologyView, TargetCoordsetView>;
    Mapper mapper(srcTopo, srcCoordset, srcMatset, targetTopo, targetCoordset);
    conduit::Node n_opts;
    n_opts["source/matsetName"] = "mat";
    n_opts["target/topologyName"] = "target1";
    n_opts["target/matsetName"] = "target1_matset";
    mapper.execute(n_dev, n_opts, n_dev);
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
    auto srcMatset = views::make_unibuffer_matset<int, float, 4>::view(n_srcMatset);
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
TEST(bump_topology_mapper, TopologyMapper_2D_seq)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_2D_seq");
  test_TopologyMapper<seq_exec>::test2D();
}
#if defined(AXOM_USE_OPENMP)
TEST(bump_topology_mapper, TopologyMapper_2D_omp)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_2D_omp");
  test_TopologyMapper<omp_exec>::test2D();
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(bump_topology_mapper, TopologyMapper_2D_cuda)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_2D_cuda");
  test_TopologyMapper<cuda_exec>::test2D();
}
#endif
#if defined(AXOM_USE_HIP)
TEST(bump_topology_mapper, TopologyMapper_2D_hip)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_2D_hip");
  test_TopologyMapper<hip_exec>::test2D();
}
#endif

//------------------------------------------------------------------------------
TEST(bump_topology_mapper, TopologyMapper_3D_seq)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_3D_seq");
  test_TopologyMapper<seq_exec>::test3D();
}
#if defined(AXOM_USE_OPENMP)
TEST(bump_topology_mapper, TopologyMapper_3D_omp)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_3D_omp");
  test_TopologyMapper<omp_exec>::test3D();
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(bump_topology_mapper, TopologyMapper_3D_cuda)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_3D_cuda");
  test_TopologyMapper<cuda_exec>::test3D();
}
#endif
#if defined(AXOM_USE_HIP)
TEST(bump_topology_mapper, TopologyMapper_3D_hip)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_3D_hip");
  test_TopologyMapper<hip_exec>::test3D();
}
#endif

//------------------------------------------------------------------------------
TEST(bump_topology_mapper, TopologyMapper_Polyhedral_seq)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_Polyhedral_seq");
  test_TopologyMapper<seq_exec>::testPolyhedral();
}
#if defined(AXOM_USE_OPENMP)
TEST(bump_topology_mapper, TopologyMapper_Polyhedral_omp)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_Polyhedral_omp");
  test_TopologyMapper<omp_exec>::testPolyhedral();
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(bump_topology_mapper, TopologyMapper_Polyhedral_cuda)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_Polyhedral_cuda");
  test_TopologyMapper<cuda_exec>::testPolyhedral();
}
#endif
#if defined(AXOM_USE_HIP)
TEST(bump_topology_mapper, TopologyMapper_Polyhedral_hip)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_Polyhedral_hip");
  test_TopologyMapper<hip_exec>::testPolyhedral();
}
#endif

//------------------------------------------------------------------------------
TEST(bump_topology_mapper, TopologyMapper_Polygonal_seq)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_Polygonal_seq");
  test_TopologyMapper_Polygonal<seq_exec>::test();
}
#if defined(AXOM_USE_OPENMP)
TEST(bump_topology_mapper, TopologyMapper_Polygonal_omp)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_Polygonal_omp");
  test_TopologyMapper_Polygonal<omp_exec>::test();
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(bump_topology_mapper, TopologyMapper_Polygonal_cuda)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_Polygonal_cuda");
  test_TopologyMapper_Polygonal<cuda_exec>::test();
}
#endif
#if defined(AXOM_USE_HIP)
TEST(bump_topology_mapper, TopologyMapper_Polygonal_hip)
{
  AXOM_ANNOTATE_SCOPE("TopologyMapper_Polygonal_hip");
  test_TopologyMapper_Polygonal<hip_exec>::test();
}
#endif

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
