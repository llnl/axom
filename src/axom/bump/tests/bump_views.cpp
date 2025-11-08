// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/bump.hpp"
#include "axom/primal.hpp"
#include "axom/bump/tests/blueprint_testing_data_helpers.hpp"
#include "axom/bump/tests/blueprint_testing_helpers.hpp"

#include <cmath>

namespace utils = axom::bump::utilities;
namespace views = axom::bump::views;

std::string baselineDirectory()
{
  return pjoin(dataDirectory(), "bump", "regression", "bump_views");
}

//------------------------------------------------------------------------------
// Global test application object.
axom::blueprint::testing::TestApplication TestApp;

//------------------------------------------------------------------------------
TEST(bump_views, shape2conduitName)
{
  EXPECT_STREQ(axom::bump::views::LineShape<int>::name(), "line");
  EXPECT_STREQ(axom::bump::views::LineShape<long>::name(), "line");

  EXPECT_STREQ(axom::bump::views::TriShape<int>::name(), "tri");
  EXPECT_STREQ(axom::bump::views::TriShape<long>::name(), "tri");

  EXPECT_STREQ(axom::bump::views::QuadShape<int>::name(), "quad");
  EXPECT_STREQ(axom::bump::views::QuadShape<long>::name(), "quad");

  EXPECT_STREQ(axom::bump::views::PolygonShape<int>::name(), "polygonal");
  EXPECT_STREQ(axom::bump::views::PolygonShape<long>::name(), "polygonal");

  EXPECT_STREQ(axom::bump::views::TetShape<int>::name(), "tet");
  EXPECT_STREQ(axom::bump::views::TetShape<long>::name(), "tet");

  EXPECT_STREQ(axom::bump::views::PyramidShape<int>::name(), "pyramid");
  EXPECT_STREQ(axom::bump::views::PyramidShape<long>::name(), "pyramid");

  EXPECT_STREQ(axom::bump::views::WedgeShape<int>::name(), "wedge");
  EXPECT_STREQ(axom::bump::views::WedgeShape<long>::name(), "wedge");

  EXPECT_STREQ(axom::bump::views::HexShape<int>::name(), "hex");
  EXPECT_STREQ(axom::bump::views::HexShape<long>::name(), "hex");
}

//------------------------------------------------------------------------------
template <typename ShapeType, typename VariableShapeType>
void compareShapes(const ShapeType &shape1, const VariableShapeType &shape2)
{
  using ConnType = typename ShapeType::ConnectivityType;

  EXPECT_EQ(shape1.dimension(), shape2.dimension());

  EXPECT_EQ(shape1.numberOfNodes(), shape2.numberOfNodes());
  for(axom::IndexType i = 0; i < shape1.numberOfNodes(); i++)
  {
    EXPECT_EQ(shape1.getId(i), shape2.getId(i));
  }

  EXPECT_EQ(shape1.numberOfEdges(), shape2.numberOfEdges());
  for(axom::IndexType i = 0; i < shape1.numberOfEdges(); i++)
  {
    const auto edge1 = shape1.getEdge(i);
    const auto edge2 = shape1.getEdge(i);
    EXPECT_EQ(edge1, edge2);
  }

  ConnType face1[8], face2[8];
  EXPECT_EQ(shape1.numberOfFaces(), shape2.numberOfFaces());
  for(axom::IndexType i = 0; i < shape1.numberOfFaces(); i++)
  {
    EXPECT_EQ(shape1.numberOfNodesInFace(i), shape2.numberOfNodesInFace(i));

    axom::IndexType numIds1 = 0, numIds2 = 0;
    shape1.getFace(i, face1, numIds1);
    shape2.getFace(i, face2, numIds2);
    EXPECT_EQ(numIds1, numIds2);
    for(axom::IndexType j = 0; j < numIds1; j++)
    {
      EXPECT_EQ(face1[j], face2[j]);
    }
  }
}

TEST(bump_views, shape_faces)
{
  using ConnType = int;

  ConnType face[5];
  axom::IndexType numIds;

  // Point
  ConnType point_ids[] = {10};
  views::PointShape<ConnType> pointShape(axom::ArrayView<ConnType>(point_ids, 1));
  EXPECT_EQ(pointShape.numberOfFaces(), 0);
  EXPECT_EQ(pointShape.getId(0), point_ids[0]);
  views::VariableShape<ConnType> pointVarShape(views::Point_ShapeID,
                                               axom::ArrayView<ConnType>(point_ids, 1));
  compareShapes(pointShape, pointVarShape);

  // Line
  ConnType line_ids[] = {10, 20};
  views::LineShape<ConnType> lineShape(axom::ArrayView<ConnType>(line_ids, 2));
  EXPECT_EQ(lineShape.numberOfFaces(), 0);
  lineShape.getFace(0, face, numIds);
  EXPECT_EQ(numIds, 0);
  EXPECT_EQ(lineShape.getId(0), line_ids[0]);
  EXPECT_EQ(lineShape.getId(1), line_ids[1]);
  views::VariableShape<ConnType> lineVarShape(views::Line_ShapeID,
                                              axom::ArrayView<ConnType>(line_ids, 2));
  compareShapes(lineShape, lineVarShape);

  // Tri
  ConnType tri_ids[] = {10, 20, 30};
  views::TriShape<ConnType> triShape(axom::ArrayView<ConnType>(tri_ids, 3));
  EXPECT_EQ(triShape.numberOfFaces(), 1);
  triShape.getFace(0, face, numIds);
  EXPECT_EQ(numIds, 3);
  EXPECT_EQ(face[0], tri_ids[0]);
  EXPECT_EQ(face[1], tri_ids[1]);
  EXPECT_EQ(face[2], tri_ids[2]);
  EXPECT_EQ(triShape.getId(0), tri_ids[0]);
  EXPECT_EQ(triShape.getId(1), tri_ids[1]);
  EXPECT_EQ(triShape.getId(2), tri_ids[2]);
  views::VariableShape<ConnType> triVarShape(views::Tri_ShapeID,
                                             axom::ArrayView<ConnType>(tri_ids, 3));
  compareShapes(triShape, triVarShape);

  // Quad
  ConnType quad_ids[] = {10, 20, 30, 40};
  views::QuadShape<ConnType> quadShape(axom::ArrayView<ConnType>(quad_ids, 4));
  EXPECT_EQ(quadShape.numberOfFaces(), 1);
  quadShape.getFace(0, face, numIds);
  EXPECT_EQ(numIds, 4);
  EXPECT_EQ(face[0], quad_ids[0]);
  EXPECT_EQ(face[1], quad_ids[1]);
  EXPECT_EQ(face[2], quad_ids[2]);
  EXPECT_EQ(face[3], quad_ids[3]);
  EXPECT_EQ(quadShape.getId(0), quad_ids[0]);
  EXPECT_EQ(quadShape.getId(1), quad_ids[1]);
  EXPECT_EQ(quadShape.getId(2), quad_ids[2]);
  EXPECT_EQ(quadShape.getId(3), quad_ids[3]);
  views::VariableShape<ConnType> quadVarShape(views::Quad_ShapeID,
                                              axom::ArrayView<ConnType>(quad_ids, 4));
  compareShapes(quadShape, quadVarShape);

  // Polygon
  ConnType polygon_ids[] = {10, 20, 30, 40, 50};
  views::PolygonShape<ConnType> polyShape(axom::ArrayView<ConnType>(polygon_ids, 5));
  EXPECT_EQ(polyShape.numberOfFaces(), 1);
  polyShape.getFace(0, face, numIds);
  EXPECT_EQ(numIds, 5);
  EXPECT_EQ(face[0], polygon_ids[0]);
  EXPECT_EQ(face[1], polygon_ids[1]);
  EXPECT_EQ(face[2], polygon_ids[2]);
  EXPECT_EQ(face[3], polygon_ids[3]);
  EXPECT_EQ(face[4], polygon_ids[4]);
  EXPECT_EQ(polyShape.getId(0), polygon_ids[0]);
  EXPECT_EQ(polyShape.getId(1), polygon_ids[1]);
  EXPECT_EQ(polyShape.getId(2), polygon_ids[2]);
  EXPECT_EQ(polyShape.getId(3), polygon_ids[3]);
  EXPECT_EQ(polyShape.getId(4), polygon_ids[4]);
  views::VariableShape<ConnType> polyVarShape(views::Polygon_ShapeID,
                                              axom::ArrayView<ConnType>(polygon_ids, 5));
  compareShapes(polyShape, polyVarShape);

  // Tet
  ConnType tet_ids[] = {10, 20, 30, 40};
  views::TetShape<ConnType> tetShape(axom::ArrayView<ConnType>(tet_ids, 4));
  EXPECT_EQ(tetShape.numberOfFaces(), 4);
  axom::IndexType tet_nids[4] = {3, 3, 3, 3};
  ConnType tet_faces[4][3] = {{10, 20, 40}, {20, 30, 40}, {30, 10, 40}, {10, 30, 20}};
  for(int f = 0; f < 4; f++)
  {
    tetShape.getFace(f, face, numIds);
    EXPECT_EQ(numIds, tet_nids[f]);
    for(axom::IndexType i = 0; i < tet_nids[f]; i++)
    {
      EXPECT_EQ(face[i], tet_faces[f][i]);
    }
  }
  views::VariableShape<ConnType> tetVarShape(views::Tet_ShapeID,
                                             axom::ArrayView<ConnType>(tet_ids, 4));
  compareShapes(tetShape, tetVarShape);

  // Pyramid
  ConnType pyr_ids[] = {10, 20, 30, 40, 50};
  views::PyramidShape<ConnType> pyrShape(axom::ArrayView<ConnType>(pyr_ids, 5));
  EXPECT_EQ(pyrShape.numberOfFaces(), 5);
  axom::IndexType pyr_nids[5] = {4, 3, 3, 3, 3};
  ConnType pyr_faces[5][4] = {{40, 30, 20, 10},
                              {10, 20, 50, -1},
                              {20, 30, 50, -1},
                              {30, 40, 50, -1},
                              {40, 10, 50, -1}};
  for(int f = 0; f < 4; f++)
  {
    pyrShape.getFace(f, face, numIds);
    EXPECT_EQ(numIds, pyr_nids[f]);
    for(axom::IndexType i = 0; i < pyr_nids[f]; i++)
    {
      EXPECT_EQ(face[i], pyr_faces[f][i]);
    }
  }
  views::VariableShape<ConnType> pyrVarShape(views::Pyramid_ShapeID,
                                             axom::ArrayView<ConnType>(pyr_ids, 5));
  compareShapes(pyrShape, pyrVarShape);

  // Wedge
  ConnType wed_ids[] = {10, 20, 30, 40, 50, 60};
  views::WedgeShape<ConnType> wedShape(axom::ArrayView<ConnType>(wed_ids, 6));
  EXPECT_EQ(wedShape.numberOfFaces(), 5);
  axom::IndexType wed_nids[5] = {3, 3, 4, 4, 4};
  ConnType wed_faces[5][4] = {{10, 30, 20, -1},
                              {40, 50, 60, -1},
                              {10, 20, 50, 40},
                              {20, 30, 60, 50},
                              {30, 10, 40, 60}};
  for(int f = 0; f < 4; f++)
  {
    wedShape.getFace(f, face, numIds);
    EXPECT_EQ(numIds, wed_nids[f]);
    for(axom::IndexType i = 0; i < wed_nids[f]; i++)
    {
      EXPECT_EQ(face[i], wed_faces[f][i]);
    }
  }
  views::VariableShape<ConnType> wedVarShape(views::Wedge_ShapeID,
                                             axom::ArrayView<ConnType>(wed_ids, 6));
  compareShapes(wedShape, wedVarShape);

  // Hex
  ConnType hex_ids[] = {10, 20, 30, 40, 50, 60, 70, 80};
  views::HexShape<ConnType> hexShape(axom::ArrayView<ConnType>(hex_ids, 8));
  EXPECT_EQ(hexShape.numberOfFaces(), 6);
  axom::IndexType hex_nids[6] = {4, 4, 4, 4, 4, 4};
  ConnType hex_faces[6][4] = {{40, 10, 50, 80},
                              {20, 30, 70, 60},
                              {10, 20, 60, 50},
                              {40, 80, 70, 30},
                              {10, 40, 30, 20},
                              {50, 60, 70, 80}};
  for(int f = 0; f < 6; f++)
  {
    hexShape.getFace(f, face, numIds);
    EXPECT_EQ(numIds, hex_nids[f]);
    for(axom::IndexType i = 0; i < hex_nids[f]; i++)
    {
      EXPECT_EQ(face[i], hex_faces[f][i]);
    }
  }
  views::VariableShape<ConnType> hexVarShape(views::Hex_ShapeID,
                                             axom::ArrayView<ConnType>(hex_ids, 8));
  compareShapes(hexShape, hexVarShape);
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_node_to_arrayview
{
  static int constexpr sum(int n)
  {
    int s = 0;
    for(int i = 0; i < n; i++) s += i;
    return s;
  }

  static void test()
  {
    std::vector<int> dtypes {conduit::DataType::INT8_ID,
                             conduit::DataType::INT16_ID,
                             conduit::DataType::INT32_ID,
                             conduit::DataType::INT64_ID,
                             conduit::DataType::UINT8_ID,
                             conduit::DataType::UINT16_ID,
                             conduit::DataType::UINT32_ID,
                             conduit::DataType::UINT64_ID,
                             conduit::DataType::FLOAT32_ID,
                             conduit::DataType::FLOAT64_ID};
    constexpr int n = 16;
    utils::ConduitAllocateThroughAxom<ExecSpace> c2a;
    for(int dtype : dtypes)
    {
      // Make a node and fill it with data.
      conduit::Node n_data;
      n_data.set_allocator(c2a.getConduitAllocatorID());
      n_data.set(conduit::DataType(dtype, n));

      int sumValues = 0;
      axom::bump::views::Node_to_ArrayView(n_data, [&](auto dataView) {
        sumValues = testBody(dataView, n);
      });

      EXPECT_EQ(sumValues, sum(n));
    }
  }

  template <typename DataView>
  static int testBody(DataView dataView, int n)
  {
    using value_type = typename DataView::value_type;

    std::cout << axom::bump::views::array_view_traits<DataView>::name() << std::endl;

    // Make sure we can store values in dataView
    axom::for_all<ExecSpace>(
      n,
      AXOM_LAMBDA(axom::IndexType index) { dataView[index] = static_cast<value_type>(index); });

    // Read the values and sum them.
    axom::ReduceSum<ExecSpace, value_type> sumValues_reduce(0);
    axom::for_all<ExecSpace>(
      n,
      AXOM_LAMBDA(axom::IndexType index) { sumValues_reduce += dataView[index]; });
    return static_cast<int>(sumValues_reduce.get());
  }
};

TEST(bump_views, node_to_arrayview_seq) { test_node_to_arrayview<seq_exec>::test(); }
#if defined(AXOM_USE_OPENMP)
TEST(bump_views, node_to_arrayview_omp) { test_node_to_arrayview<omp_exec>::test(); }
#endif
#if defined(AXOM_USE_CUDA)
TEST(bump_views, node_to_arrayview_cuda) { test_node_to_arrayview<cuda_exec>::test(); }
#endif
#if defined(AXOM_USE_HIP)
TEST(bump_views, node_to_arrayview_hip) { test_node_to_arrayview<hip_exec>::test(); }
#endif

//------------------------------------------------------------------------------
TEST(bump_views, explicit_coordsetview)
{
  axom::Array<float> x {{0., 1., 2., 3., 4., 5.}};
  axom::Array<float> y {{10., 11., 12., 13., 14., 15.}};
  axom::Array<float> z {{20., 21., 22., 23., 24., 25.}};

  axom::bump::views::ExplicitCoordsetView<float, 2> view2d(x.view(), y.view());
  EXPECT_EQ(view2d.size(), 6);
  for(axom::IndexType i = 0; i < view2d.size(); i++)
  {
    axom::primal::Point<float, 2> P({x[i], y[i]});
    EXPECT_EQ(view2d.getPoint(i), P);
    EXPECT_EQ(view2d[i], P);
  }

  axom::bump::views::ExplicitCoordsetView<float, 3> view3d(x.view(), y.view(), z.view());
  EXPECT_EQ(view3d.size(), 6);
  for(axom::IndexType i = 0; i < view3d.size(); i++)
  {
    axom::primal::Point<float, 3> P({x[i], y[i], z[i]});
    EXPECT_EQ(view3d.getPoint(i), P);
    EXPECT_EQ(view3d[i], P);
  }
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_structured_topology_view_rectilinear
{
  static void test()
  {
    conduit::Node hostMesh;
    create(hostMesh);

    // host->device
    conduit::Node deviceMesh;
    utils::copy<ExecSpace>(deviceMesh, hostMesh);

    // Make results view on device.
    constexpr int nzones = 9;
    axom::Array<axom::IndexType> results(nzones,
                                         nzones,
                                         axom::execution_space<ExecSpace>::allocatorID());
    auto resultsView = results.view();

    // Execute the kernel for each zone (find max node number in zone).
    auto topoView =
      axom::bump::views::make_rectilinear_topology<2>::view(deviceMesh["topologies/mesh"]);
    axom::for_all<ExecSpace>(
      topoView.numberOfZones(),
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
        const auto zone = topoView.zone(zoneIndex);
        axom::IndexType m = -1;
        for(const auto &id : zone.getIds())
        {
          m = axom::utilities::max(static_cast<axom::IndexType>(id), m);
        }
        resultsView[zoneIndex] = m;
      });

    // device->host
    axom::Array<axom::IndexType> hostResults(nzones,
                                             nzones,
                                             axom::execution_space<axom::SEQ_EXEC>::allocatorID());
    axom::copy(hostResults.data(), results.data(), nzones * sizeof(axom::IndexType));

    // Compare.
    const axom::IndexType expected[] = {5, 6, 7, 9, 10, 11, 13, 14, 15};
    for(int i = 0; i < nzones; i++)
    {
      EXPECT_EQ(hostResults[i], expected[i]);
    }
  }

  static void create(conduit::Node &mesh)
  {
    std::vector<int> dims {4, 4};
    axom::blueprint::testing::data::braid("rectilinear", dims, mesh);
  }
};

TEST(bump_views, stopo_rectilinear_2d_seq)
{
  test_structured_topology_view_rectilinear<seq_exec>::test();
}
#if defined(AXOM_USE_OPENMP)
TEST(bump_views, stopo_rectilinear_2d_omp)
{
  test_structured_topology_view_rectilinear<omp_exec>::test();
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(bump_views, stopo_rectilinear_2d_cuda)
{
  test_structured_topology_view_rectilinear<cuda_exec>::test();
}
#endif
#if defined(AXOM_USE_HIP)
TEST(bump_views, stopo_rectilinear_2d_hip)
{
  test_structured_topology_view_rectilinear<hip_exec>::test();
}
#endif

//------------------------------------------------------------------------------
struct test_strided_structured
{
  static void test()
  {
    conduit::Node hostMesh;
    axom::blueprint::testing::data::strided_structured<2>(hostMesh);
    //  hostMesh.print();

    axom::bump::views::dispatch_explicit_coordset(hostMesh["coordsets/coords"], [&](auto coordsetView) {
      axom::bump::views::dispatch_structured_topology<axom::bump::views::select_dimensions(2)>(
        hostMesh["topologies/mesh"],
        [&](const std::string &AXOM_UNUSED_PARAM(shape), auto topoView) {
          execute(coordsetView, topoView);
        });
    });
  }

  template <typename CoordsetView, typename TopologyView>
  static void execute(CoordsetView coordsetView, TopologyView topoView)
  {
    using ExecSpace = seq_exec;

    // These are the expected node ids for this strided structured mesh.
    // clang-format off
    const axom::Array<int> expectedNodes {{16, 17, 24, 23,
                                           17, 18, 25, 24,
                                           18, 19, 26, 25,
                                           23, 24, 31, 30,
                                           24, 25, 32, 31,
                                           25, 26, 33, 32}};
    // clang-format on
    auto expectedNodesView = expectedNodes.view();
    axom::IndexType n4 = expectedNodesView.size();

    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    axom::Array<int> actualNodes(n4, n4, allocatorID);
    axom::Array<int> logicalNodes(n4 * 2, n4 * 2, allocatorID);
    auto actualNodesView = actualNodes.view();
    auto logicalNodesView = logicalNodes.view();

    // Traverse the zones in the mesh and gather node ids
    axom::for_all<ExecSpace>(
      topoView.numberOfZones(),
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
        const auto zone = topoView.zone(zoneIndex);
        const auto nodeIndexing = topoView.indexing().expand();

        // Get node ids for zone.
        const auto ids = zone.getIds();
        for(axom::IndexType i = 0; i < ids.size(); i++)
        {
          actualNodesView[zoneIndex * 4 + i] = ids[i];

          // Get the logical local id for the id.
          const auto index = nodeIndexing.GlobalToLocal(ids[i]);
          const auto logical = nodeIndexing.IndexToLogicalIndex(index);
          logicalNodesView[(zoneIndex * 4 + i) * 2 + 0] = logical[0];
          logicalNodesView[(zoneIndex * 4 + i) * 2 + 1] = logical[1];
        }
      });

    for(axom::IndexType i = 0; i < n4; i++)
    {
      EXPECT_EQ(expectedNodesView[i], actualNodesView[i]);
    }

    // Check coordinates
    for(axom::IndexType i = 0; i < n4; i++)
    {
      const auto id = actualNodesView[i];

      // Get coordinate from coordsetView.
      const auto pt = coordsetView[id];

      // Get the logical local id for the id.
      const auto logicalI = logicalNodesView[i * 2 + 0];
      const auto logicalJ = logicalNodesView[i * 2 + 1];

      // Expected coordinate
      double x = (3. + 1. / 3.) * static_cast<double>(logicalI - 1);
      const double yvals[] = {-2, 2, 6};
      double y = yvals[logicalJ];

      const double dx = pt[0] - x;
      const double dy = pt[1] - y;
      double d = sqrt(dx * dx + dy * dy);

      EXPECT_TRUE(d < 1.e-10);
    }
  }
};

TEST(bump_views, strided_structured_seq) { test_strided_structured::test(); }

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_braid2d_mat
{
  static void test(const std::string &type, const std::string &mattype, const std::string &name)
  {
    namespace utils = axom::bump::utilities;
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    axom::StackArray<axom::IndexType, 2> dims {10, 10};
    axom::StackArray<axom::IndexType, 2> zoneDims {dims[0] - 1, dims[1] - 1};
    const axom::IndexType nzones = zoneDims[0] * zoneDims[1];

    // Create the data
    const bool cleanMats = false;
    conduit::Node hostMesh, deviceMesh;
    axom::blueprint::testing::data::braid(type, dims, hostMesh);
    axom::blueprint::testing::data::make_matset(mattype, "mesh", zoneDims, cleanMats, hostMesh);
    utils::copy<ExecSpace>(deviceMesh, hostMesh);
    TestApp.saveVisualization(name + "_orig", hostMesh);

    if(mattype == "unibuffer")
    {
      // clang-format off
      // _bump_views_matsetview_begin
      using MatsetView = axom::bump::views::UnibufferMaterialView<int, float, 3>;
      MatsetView matsetView;
      matsetView.set(utils::make_array_view<int>(deviceMesh["matsets/mat/material_ids"]),
                     utils::make_array_view<float>(deviceMesh["matsets/mat/volume_fractions"]),
                     utils::make_array_view<int>(deviceMesh["matsets/mat/sizes"]),
                     utils::make_array_view<int>(deviceMesh["matsets/mat/offsets"]),
                     utils::make_array_view<int>(deviceMesh["matsets/mat/indices"]));
      // _bump_views_matsetview_end
      // clang-format on
      test_matsetview(nzones, matsetView, allocatorID);
    }
    else if(mattype == "multibuffer")
    {
      axom::bump::views::dispatch_material_multibuffer(
        deviceMesh["matsets/mat"],
        [&](auto matsetView) { test_matsetview(nzones, matsetView, allocatorID); });
    }
    else if(mattype == "element_dominant")
    {
      axom::bump::views::dispatch_material_element_dominant(
        deviceMesh["matsets/mat"],
        [&](auto matsetView) { test_matsetview(nzones, matsetView, allocatorID); });
    }
    else if(mattype == "material_dominant")
    {
      axom::bump::views::dispatch_material_material_dominant(
        deviceMesh["matsets/mat"],
        [&](auto matsetView) { test_matsetview(nzones, matsetView, allocatorID); });
    }
  }

  template <typename MatsetView>
  static void test_matsetview(axom::IndexType nzones, MatsetView matsetView, int allocatorID)
  {
    // These values are used in that material_map.
    constexpr int MATA = 22;
    constexpr int MATB = 66;
    constexpr int MATC = 33;
    // The zone ids that are being queried.
    const int zoneids[] = {0, 36, 40};

    // clang-format off
    int results[] = {/*nzones*/ static_cast<int>(nzones),
                     /*contains mat*/ 0, 1, 0, /*nmats in zone*/ 1, /*ids.size*/ 1, /*mats in zone*/ MATB, -1, -1,
                     /*contains mat*/ 1, 1, 0, /*nmats in zone*/ 2, /*ids.size*/ 2, /*mats in zone*/ MATA, MATB, -1,
                     /*contains mat*/ 1, 1, 1, /*nmats in zone*/ 3, /*ids.size*/ 3, /*mats in zone*/ MATA, MATB, MATC};
    // clang-format on
    constexpr int nTestZones = sizeof(zoneids) / sizeof(int);

    // Get zoneids into zoneidsView for device.
    axom::Array<int> zoneidsArray(nTestZones, nTestZones, allocatorID);
    axom::copy(zoneidsArray.data(), zoneids, sizeof(int) * nTestZones);
    auto zoneidsView = zoneidsArray.view();

    // Allocate results array on device.
    constexpr int nResults = sizeof(results) / sizeof(int);
    axom::Array<int> resultsArrayDevice(nResults, nResults, allocatorID);
    auto resultsView = resultsArrayDevice.view();

    // Fill in resultsView on the device.
    constexpr int nResultsPerZone = 8;
    axom::for_all<ExecSpace>(
      nTestZones,
      AXOM_LAMBDA(axom::IndexType index) {
        if(index == 0)
        {
          // Compute number of zones here since some views need to look inside
          // data to determine the number of zones.
          resultsView[0] = matsetView.numberOfZones();
        }
        const int offset = 1 + nResultsPerZone * index;
        // contains mat
        resultsView[offset + 0] = matsetView.zoneContainsMaterial(zoneidsView[index], MATA) ? 1 : 0;
        resultsView[offset + 1] = matsetView.zoneContainsMaterial(zoneidsView[index], MATB) ? 1 : 0;
        resultsView[offset + 2] = matsetView.zoneContainsMaterial(zoneidsView[index], MATC) ? 1 : 0;
        // nmats in zone
        resultsView[offset + 3] = matsetView.numberOfMaterials(zoneidsView[index]);

        typename MatsetView::IDList ids {};
        typename MatsetView::VFList vfs {};
        // ids.size
        matsetView.zoneMaterials(zoneidsView[index], ids, vfs);
        resultsView[offset + 4] = ids.size();
        // mats in zone
        for(axom::IndexType i = 0; i < 3; i++)
        {
          resultsView[offset + 5 + i] = (i < ids.size()) ? static_cast<int>(ids[i]) : -1;
        }
      });
    // Get containsView data to the host and compare results
    std::vector<int> resultsHost(nResults);
    axom::copy(resultsHost.data(), resultsView.data(), sizeof(int) * nResults);
    for(int i = 0; i < nResults; i++)
    {
      EXPECT_EQ(results[i], resultsHost[i]);
    }

    // Test iterators.
    test_matsetview_iterators(matsetView, allocatorID);
  }

  template <typename MatsetView>
  static void test_matsetview_iterators(MatsetView matsetView, int allocatorID)
  {
    using ZoneIndex = typename MatsetView::ZoneIndex;
    // Allocate results array on device.
    const int nResults = matsetView.numberOfZones();
    axom::Array<int> resultsArrayDevice(nResults, nResults, allocatorID);
    auto resultsView = resultsArrayDevice.view();

    axom::for_all<ExecSpace>(
      matsetView.numberOfZones(),
      AXOM_LAMBDA(axom::IndexType index) {
        typename MatsetView::IDList ids {};
        typename MatsetView::VFList vfs {};
        matsetView.zoneMaterials(index, ids, vfs);

        // Get the end iterator for the zone.
        const auto end = matsetView.endZone(index);

        int eq_count = 0;
        int count = 0;

        // Make sure the iterator is for the right zone.
        eq_count += (end.zoneIndex() == static_cast<ZoneIndex>(index)) ? 1 : 0;
        count++;

        // Make sure incrementing the last iterator has no effect.
        auto end2 = end;
        end2++;
        eq_count += (end == end2) ? 1 : 0;
        count++;

        // Make sure the iterator order is the same as for the values we got from zoneMaterials().
        int i = 0;
        for(auto it = matsetView.beginZone(index); it != end; it++, i++)
        {
          eq_count += (vfs[i] == it.volume_fraction() && ids[i] == it.material_id()) ? 1 : 0;
          count++;
        }

        // Test ArrayView version of zoneMaterials().
        using IndexType = typename MatsetView::IndexType;
        using FloatType = typename MatsetView::FloatType;
        constexpr int ARRAY_SIZE = 10;
        IndexType idStorage[ARRAY_SIZE];
        FloatType vfStorage[ARRAY_SIZE];
        axom::ArrayView<IndexType> idView(idStorage, ARRAY_SIZE);
        axom::ArrayView<FloatType> vfView(vfStorage, ARRAY_SIZE);
        const auto nmats = matsetView.zoneMaterials(index, idView, vfView);
        eq_count += (nmats == ids.size()) ? 1 : 0;
        count++;
        for(axom::IndexType j = 0; j < nmats; j++)
        {
          eq_count += (vfs[j] == vfView[j] && ids[j] == idView[j]) ? 1 : 0;
          count++;
        }

        resultsView[index] = (eq_count == count) ? 1 : 0;
      });

    // Get containsView data to the host and compare results
    std::vector<int> resultsHost(nResults);
    axom::copy(resultsHost.data(), resultsView.data(), sizeof(int) * nResults);
    for(int i = 0; i < nResults; i++)
    {
      EXPECT_EQ(resultsHost[i], 1);
    }
  }
};

// Unibuffer
TEST(bump_views, matset_unibuffer_seq)
{
  test_braid2d_mat<seq_exec>::test("uniform", "unibuffer", "uniform2d_unibuffer");
}
#if defined(AXOM_USE_OPENMP)
TEST(bump_views, matset_unibuffer_omp)
{
  test_braid2d_mat<omp_exec>::test("uniform", "unibuffer", "uniform2d_unibuffer");
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(bump_views, matset_unibuffer_cuda)
{
  test_braid2d_mat<cuda_exec>::test("uniform", "unibuffer", "uniform2d_unibuffer");
}
#endif
#if defined(AXOM_USE_HIP)
TEST(bump_views, matset_unibuffer_hip)
{
  test_braid2d_mat<hip_exec>::test("uniform", "unibuffer", "uniform2d_unibuffer");
}
#endif

// Multibuffer
TEST(bump_views, matset_multibuffer_seq)
{
  test_braid2d_mat<seq_exec>::test("uniform", "multibuffer", "uniform2d_multibuffer");
}
#if defined(AXOM_USE_OPENMP)
TEST(bump_views, matset_multibuffer_omp)
{
  test_braid2d_mat<omp_exec>::test("uniform", "multibuffer", "uniform2d_multibuffer");
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(bump_views, matset_multibuffer_cuda)
{
  test_braid2d_mat<cuda_exec>::test("uniform", "multibuffer", "uniform2d_multibuffer");
}
#endif
#if defined(AXOM_USE_HIP)
TEST(bump_views, matset_multibuffer_hip)
{
  test_braid2d_mat<hip_exec>::test("uniform", "multibuffer", "uniform2d_multibuffer");
}
#endif

// Element-dominant
TEST(bump_views, matset_element_dominant_seq)
{
  test_braid2d_mat<seq_exec>::test("uniform", "element_dominant", "uniform2d_element_dominant");
}
#if defined(AXOM_USE_OPENMP)
TEST(bump_views, matset_element_dominant_omp)
{
  test_braid2d_mat<omp_exec>::test("uniform", "element_dominant", "uniform2d_element_dominant");
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(bump_views, matset_element_dominant_cuda)
{
  test_braid2d_mat<cuda_exec>::test("uniform", "element_dominant", "uniform2d_element_dominant");
}
#endif
#if defined(AXOM_USE_HIP)
TEST(bump_views, matset_element_dominant_hip)
{
  test_braid2d_mat<hip_exec>::test("uniform", "element_dominant", "uniform2d_element_dominant");
}
#endif

// Material-dominant
TEST(bump_views, matset_material_dominant_seq)
{
  test_braid2d_mat<seq_exec>::test("uniform", "material_dominant", "uniform2d_material_dominant");
}
#if defined(AXOM_USE_OPENMP)
TEST(bump_views, matset_material_dominant_omp)
{
  test_braid2d_mat<omp_exec>::test("uniform", "material_dominant", "uniform2d_material_dominant");
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(bump_views, matset_material_dominant_cuda)
{
  test_braid2d_mat<cuda_exec>::test("uniform", "material_dominant", "uniform2d_material_dominant");
}
#endif
#if defined(AXOM_USE_HIP)
TEST(bump_views, matset_material_dominant_hip)
{
  test_braid2d_mat<hip_exec>::test("uniform", "material_dominant", "uniform2d_material_dominant");
}
#endif

TEST(bump_views, matset_multibuffer)
{
  const char *yaml = R"(
matsets:
  matset:
    topology: topology
    volume_fractions:
      a:
        values: [0, 0, 0, 0.33, 0, 0.3]  # [0, 0, 0, a1, 0, a0]
        indices: [5, 3]
      b:
        values: [0, 0.7, 1., 0.67, 0]    # [0, b0, b2, b1, 0]
        indices: [1, 3, 2]
    material_map: # (optional)
      a: 0
      b: 1
)";
  conduit::Node matsets;
  matsets.parse(yaml);
  const conduit::Node &n_matset = matsets["matsets/matset"];
  axom::bump::views::dispatch_material_multibuffer(n_matset, [&](auto matsetView) {
    using IDList = typename decltype(matsetView)::IDList;
    using VFList = typename decltype(matsetView)::VFList;
    using VFType = typename VFList::value_type;

    EXPECT_EQ(matsetView.numberOfZones(), 3);
    EXPECT_EQ(matsetView.numberOfMaterials(0), 2);
    EXPECT_EQ(matsetView.numberOfMaterials(1), 2);
    EXPECT_EQ(matsetView.numberOfMaterials(2), 1);

    IDList m0, m1, m2;
    VFList vf0, vf1, vf2;
    matsetView.zoneMaterials(0, m0, vf0);
    EXPECT_EQ(m0.size(), 2);
    EXPECT_EQ(vf0.size(), 2);
    EXPECT_EQ(m0[0], 0);
    EXPECT_EQ(m0[1], 1);
    EXPECT_EQ(vf0[0], 0.3);
    EXPECT_EQ(vf0[1], 0.7);

    VFType vf;
    EXPECT_TRUE(matsetView.zoneContainsMaterial(0, 0, vf));
    EXPECT_EQ(vf, 0.3);
    EXPECT_TRUE(matsetView.zoneContainsMaterial(0, 1, vf));
    EXPECT_EQ(vf, 0.7);

    matsetView.zoneMaterials(1, m1, vf1);
    EXPECT_EQ(m1.size(), 2);
    EXPECT_EQ(vf1.size(), 2);
    EXPECT_EQ(m1[0], 0);
    EXPECT_EQ(m1[1], 1);
    EXPECT_EQ(vf1[0], 0.33);
    EXPECT_EQ(vf1[1], 0.67);

    EXPECT_TRUE(matsetView.zoneContainsMaterial(1, 0, vf));
    EXPECT_EQ(vf, 0.33);
    EXPECT_TRUE(matsetView.zoneContainsMaterial(1, 1, vf));
    EXPECT_EQ(vf, 0.67);

    matsetView.zoneMaterials(2, m2, vf2);
    EXPECT_EQ(m2.size(), 1);
    EXPECT_EQ(m2[0], 1);
    EXPECT_EQ(vf2[0], 1);

    EXPECT_FALSE(matsetView.zoneContainsMaterial(2, 0, vf));
    EXPECT_EQ(vf, 0.);
    EXPECT_TRUE(matsetView.zoneContainsMaterial(2, 1, vf));
    EXPECT_EQ(vf, 1.);
  });
}

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
