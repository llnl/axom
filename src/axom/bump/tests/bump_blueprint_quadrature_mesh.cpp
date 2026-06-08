// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/bump/GenerateQuadratureMesh.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"
#include "axom/bump/views/dispatch_coordset.hpp"
#include "axom/bump/views/dispatch_topology.hpp"

#include "conduit.hpp"
#include "conduit_blueprint.hpp"

#include <type_traits>

namespace
{

template <typename T, typename U>
bool compareArrayView(axom::ArrayView<T> lhs, axom::ArrayView<U> rhs)
{
  if(lhs.size() != rhs.size())
  {
    return false;
  }

  for(axom::IndexType i = 0; i < lhs.size(); ++i)
  {
    if(lhs[i] != rhs[i])
    {
      return false;
    }
  }
  return true;
}

void setNodeValues(conduit::Node& node, axom::ArrayView<const double> values)
{
  node.set(conduit::DataType::float64(values.size()));
  auto* data = node.as_float64_ptr();
  for(axom::IndexType i = 0; i < values.size(); ++i)
  {
    data[i] = values[i];
  }
}

void setNodeValues(conduit::Node& node, axom::ArrayView<const conduit::index_t> values)
{
  node.set(conduit::DataType::index_t(values.size()));
  auto* data = node.as_index_t_ptr();
  for(axom::IndexType i = 0; i < values.size(); ++i)
  {
    data[i] = values[i];
  }
}

conduit::Node makeQuadMesh(const std::string& topoName = "mesh")
{
  conduit::Node mesh;

  mesh["coordsets/coords/type"] = "explicit";
  const axom::Array<double> x {{0., 1., 0., 1.}};
  const axom::Array<double> y {{0., 0., 1., 1.}};
  setNodeValues(mesh["coordsets/coords/values/x"], x.view());
  setNodeValues(mesh["coordsets/coords/values/y"], y.view());

  mesh["topologies"][topoName]["type"] = "unstructured";
  mesh["topologies"][topoName]["coordset"] = "coords";
  mesh["topologies"][topoName]["elements/shape"] = "quad";
  const axom::Array<conduit::index_t> connectivity {{0, 1, 3, 2}};
  setNodeValues(mesh["topologies"][topoName]["elements/connectivity"], connectivity.view());

  return mesh;
}

conduit::Node makeHexMesh(const std::string& topoName = "mesh")
{
  conduit::Node mesh;

  mesh["coordsets/coords/type"] = "explicit";
  const axom::Array<double> x {{0., 1., 0., 1., 0., 1., 0., 1.}};
  const axom::Array<double> y {{0., 0., 1., 1., 0., 0., 1., 1.}};
  const axom::Array<double> z {{0., 0., 0., 0., 1., 1., 1., 1.}};
  setNodeValues(mesh["coordsets/coords/values/x"], x.view());
  setNodeValues(mesh["coordsets/coords/values/y"], y.view());
  setNodeValues(mesh["coordsets/coords/values/z"], z.view());

  mesh["topologies"][topoName]["type"] = "unstructured";
  mesh["topologies"][topoName]["coordset"] = "coords";
  mesh["topologies"][topoName]["elements/shape"] = "hex";
  const axom::Array<conduit::index_t> connectivity {{0, 1, 3, 2, 4, 5, 7, 6}};
  setNodeValues(mesh["topologies"][topoName]["elements/connectivity"], connectivity.view());

  return mesh;
}

conduit::Node makeDistortedQuadMesh(const std::string& topoName = "mesh")
{
  conduit::Node mesh;

  mesh["coordsets/coords/type"] = "explicit";
  const axom::Array<double> x {{0., 2., 0., 1.}};
  const axom::Array<double> y {{0., 0., 1., 1.}};
  setNodeValues(mesh["coordsets/coords/values/x"], x.view());
  setNodeValues(mesh["coordsets/coords/values/y"], y.view());

  mesh["topologies"][topoName]["type"] = "unstructured";
  mesh["topologies"][topoName]["coordset"] = "coords";
  mesh["topologies"][topoName]["elements/shape"] = "quad";
  const axom::Array<conduit::index_t> connectivity {{0, 1, 3, 2}};
  setNodeValues(mesh["topologies"][topoName]["elements/connectivity"], connectivity.view());

  return mesh;
}

conduit::Node makeStructuredQuadMesh(const std::string& topoName = "mesh")
{
  conduit::Node mesh;

  mesh["coordsets/coords/type"] = "explicit";
  const axom::Array<double> x {{0., 1., 0., 1.}};
  const axom::Array<double> y {{0., 0., 1., 1.}};
  setNodeValues(mesh["coordsets/coords/values/x"], x.view());
  setNodeValues(mesh["coordsets/coords/values/y"], y.view());

  mesh["topologies"][topoName]["type"] = "structured";
  mesh["topologies"][topoName]["coordset"] = "coords";
  mesh["topologies"][topoName]["elements/shape"] = "quad";
  mesh["topologies"][topoName]["elements/dims/i"] = 1;
  mesh["topologies"][topoName]["elements/dims/j"] = 1;

  return mesh;
}

void generateQuadratureMesh(conduit::Node& mesh,
                            const std::string& topologyName,
                            axom::ArrayView<int> sampleResolution,
                            axom::numerics::QuadratureType quadratureType)
{
  namespace views = axom::bump::views;

  const int allocatorId = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  auto ruleX =
    axom::numerics::get_quadrature_rule(quadratureType, sampleResolution[0], allocatorId);
  auto ruleY =
    axom::numerics::get_quadrature_rule(quadratureType, sampleResolution[1], allocatorId);
  const int nz = sampleResolution.size() > 2 ? sampleResolution[2] : 1;
  auto ruleZ = axom::numerics::get_quadrature_rule(quadratureType, nz, allocatorId);

  const conduit::Node& topoNode =
    mesh.fetch_existing("topologies").fetch_existing(topologyName);
  const conduit::Node& coordsetNode =
    mesh.fetch_existing("coordsets").fetch_existing(topoNode.fetch_existing("coordset").as_string());

  views::dispatch_explicit_coordset(coordsetNode, [&](auto coordsetView) {
    using CoordsetView = typename std::decay<decltype(coordsetView)>::type;
    constexpr int supportedShapes = (CoordsetView::dimension() == 2)
      ? views::select_shapes(views::Quad_ShapeID)
      : views::select_shapes(views::Hex_ShapeID);

    views::dispatch_topology<views::select_dimensions(CoordsetView::dimension()), supportedShapes>(
      topoNode,
      [&](const auto&, auto topoView) {
        axom::bump::GenerateQuadratureMesh<axom::SEQ_EXEC,
                                           decltype(topoView),
                                           CoordsetView>
          generator(topoView, coordsetView);
        generator.setAllocatorID(allocatorId);
        generator.execute(topoNode,
                          coordsetNode,
                          "quadrature_points",
                          "quadrature_points",
                          "originalElements",
                          "quadratureWeights",
                          "quadraturePhysicalWeights",
                          ruleX,
                          ruleY,
                          ruleZ,
                          mesh);
      });
  });
}

TEST(bump_blueprint_quadrature_mesh, generate_closed_uniform_quad_mesh)
{
  conduit::Node mesh = makeQuadMesh();

  int sampleResolution[] = {2, 3};
  generateQuadratureMesh(mesh,
                         "mesh",
                         axom::ArrayView<int> {sampleResolution, 2},
                         axom::numerics::QuadratureType::ClosedUniform);

  conduit::Node info;
  EXPECT_TRUE(conduit::blueprint::mesh::verify(mesh, info)) << info.to_yaml();

  const conduit::Node& quadTopo = mesh["topologies/quadrature_points"];
  EXPECT_EQ(quadTopo["type"].as_string(), "unstructured");
  EXPECT_EQ(quadTopo["coordset"].as_string(), "quadrature_points");
  EXPECT_EQ(quadTopo["elements/shape"].as_string(), "point");

  namespace utils = axom::bump::utilities;
  const auto connView = utils::make_array_view<conduit::index_t>(
    mesh["topologies/quadrature_points/elements/connectivity"]);
  const auto sizesView =
    utils::make_array_view<conduit::index_t>(mesh["topologies/quadrature_points/elements/sizes"]);
  const auto offsetsView =
    utils::make_array_view<conduit::index_t>(mesh["topologies/quadrature_points/elements/offsets"]);
  const auto originalElementsView =
    utils::make_array_view<conduit::index_t>(mesh["fields/originalElements/values"]);
  const auto quadratureWeightsView =
    utils::make_array_view<double>(mesh["fields/quadratureWeights/values"]);
  const auto physicalQuadratureWeightsView =
    utils::make_array_view<double>(mesh["fields/quadraturePhysicalWeights/values"]);

  const axom::Array<double> expectedX {{0., 1., 0., 1., 0., 1.}};
  const axom::Array<double> expectedY {{0., 0., 0.5, 0.5, 1., 1.}};
  const axom::Array<conduit::index_t> expectedConn {{0, 1, 2, 3, 4, 5}};
  const axom::Array<conduit::index_t> expectedSizes {{1, 1, 1, 1, 1, 1}};
  const axom::Array<conduit::index_t> expectedOffsets {{0, 1, 2, 3, 4, 5}};
  const axom::Array<conduit::index_t> expectedOriginalElements {{0, 0, 0, 0, 0, 0}};
  const axom::Array<double> expectedWeights {
    {1. / 12., 1. / 12., 1. / 3., 1. / 3., 1. / 12., 1. / 12.}};

  axom::bump::views::dispatch_explicit_coordset(
    mesh["coordsets/quadrature_points"],
    [&](auto coordsetView) {
      for(axom::IndexType i = 0; i < expectedX.size(); ++i)
      {
        EXPECT_NEAR(coordsetView[i][0], expectedX[i], 1e-12);
        EXPECT_NEAR(coordsetView[i][1], expectedY[i], 1e-12);
      }
    });
  EXPECT_TRUE(compareArrayView(expectedConn.view(), connView));
  EXPECT_TRUE(compareArrayView(expectedSizes.view(), sizesView));
  EXPECT_TRUE(compareArrayView(expectedOffsets.view(), offsetsView));
  EXPECT_TRUE(compareArrayView(expectedOriginalElements.view(), originalElementsView));
  for(axom::IndexType i = 0; i < expectedWeights.size(); ++i)
  {
    EXPECT_NEAR(expectedWeights[i], quadratureWeightsView[i], 1e-12);
    EXPECT_NEAR(expectedWeights[i], physicalQuadratureWeightsView[i], 1e-12);
  }
}

TEST(bump_blueprint_quadrature_mesh, generate_open_uniform_hex_mesh)
{
  conduit::Node mesh = makeHexMesh();

  int sampleResolution[3] = {2, 1, 2};
  generateQuadratureMesh(mesh,
                         "mesh",
                         axom::ArrayView<int> {sampleResolution, 3},
                         axom::numerics::QuadratureType::OpenUniform);

  conduit::Node info;
  EXPECT_TRUE(conduit::blueprint::mesh::verify(mesh, info)) << info.to_yaml();

  EXPECT_TRUE(mesh.has_path("coordsets/quadrature_points/type")) << mesh.to_yaml();
  EXPECT_TRUE(mesh.has_path("coordsets/quadrature_points/values/x")) << mesh.to_yaml();
  EXPECT_TRUE(mesh.has_path("coordsets/quadrature_points/values/y")) << mesh.to_yaml();
  EXPECT_TRUE(mesh.has_path("coordsets/quadrature_points/values/z")) << mesh.to_yaml();

  namespace utils = axom::bump::utilities;
  const auto originalElementsView =
    utils::make_array_view<conduit::index_t>(mesh["fields/originalElements/values"]);
  const auto quadratureWeightsView =
    utils::make_array_view<double>(mesh["fields/quadratureWeights/values"]);
  const auto physicalQuadratureWeightsView =
    utils::make_array_view<double>(mesh["fields/quadraturePhysicalWeights/values"]);

  const axom::Array<double> expectedX {{1. / 3., 2. / 3., 1. / 3., 2. / 3.}};
  const axom::Array<double> expectedY {{0.5, 0.5, 0.5, 0.5}};
  const axom::Array<double> expectedZ {{1. / 3., 1. / 3., 2. / 3., 2. / 3.}};
  const axom::Array<conduit::index_t> expectedOriginalElements {{0, 0, 0, 0}};
  const axom::Array<double> expectedWeights {{0.25, 0.25, 0.25, 0.25}};

  axom::bump::views::dispatch_explicit_coordset(
    mesh["coordsets/quadrature_points"],
    [&](auto coordsetView) {
      for(axom::IndexType i = 0; i < expectedX.size(); ++i)
      {
        EXPECT_NEAR(coordsetView[i][0], expectedX[i], 1e-6);
        EXPECT_NEAR(coordsetView[i][1], expectedY[i], 1e-6);
        EXPECT_NEAR(coordsetView[i][2], expectedZ[i], 1e-6);
      }
    });
  EXPECT_TRUE(compareArrayView(expectedOriginalElements.view(), originalElementsView));
  for(axom::IndexType i = 0; i < expectedWeights.size(); ++i)
  {
    EXPECT_NEAR(expectedWeights[i], quadratureWeightsView[i], 1e-12);
    EXPECT_NEAR(expectedWeights[i], physicalQuadratureWeightsView[i], 1e-12);
  }
}

TEST(bump_blueprint_quadrature_mesh, generate_closed_uniform_structured_quad_mesh)
{
  conduit::Node mesh = makeStructuredQuadMesh();

  int sampleResolution[] = {2, 2};
  generateQuadratureMesh(mesh,
                         "mesh",
                         axom::ArrayView<int> {sampleResolution, 2},
                         axom::numerics::QuadratureType::ClosedUniform);

  conduit::Node info;
  EXPECT_TRUE(conduit::blueprint::mesh::verify(mesh, info)) << info.to_yaml();

  ASSERT_TRUE(mesh.has_path("coordsets/quadrature_points/values/x"));
  ASSERT_TRUE(mesh.has_path("fields/originalElements/values"));

  namespace utils = axom::bump::utilities;
  const auto originalElementsView =
    utils::make_array_view<conduit::index_t>(mesh["fields/originalElements/values"]);

  const axom::Array<double> expectedX {{0., 1., 0., 1.}};
  const axom::Array<double> expectedY {{0., 0., 1., 1.}};
  const axom::Array<conduit::index_t> expectedOriginalElements {{0, 0, 0, 0}};

  axom::bump::views::dispatch_explicit_coordset(
    mesh["coordsets/quadrature_points"],
    [&](auto coordsetView) {
      for(axom::IndexType i = 0; i < expectedX.size(); ++i)
      {
        EXPECT_NEAR(coordsetView[i][0], expectedX[i], 1e-12);
        EXPECT_NEAR(coordsetView[i][1], expectedY[i], 1e-12);
      }
    });
  EXPECT_TRUE(compareArrayView(expectedOriginalElements.view(), originalElementsView));
}

TEST(bump_blueprint_quadrature_mesh, mapped_zone_helper_computes_distorted_quad_measure_factor)
{
  conduit::Node mesh = makeDistortedQuadMesh();

  double lowerFactor = 0.;
  double upperFactor = 0.;

  axom::bump::views::dispatch_explicit_coordset(mesh["coordsets/coords"], [&](auto coordsetView) {
    axom::bump::views::dispatch_topology(mesh["topologies/mesh"], [&](const auto&, auto topoView) {
      const auto zone = topoView.zone(0);
      lowerFactor =
        axom::bump::detail::computePhysicalMeasureFactor(zone, coordsetView, 0.5, 0.0);
      upperFactor =
        axom::bump::detail::computePhysicalMeasureFactor(zone, coordsetView, 0.5, 1.0);
    });
  });

  EXPECT_NEAR(lowerFactor, 2.0, 1e-12);
  EXPECT_NEAR(upperFactor, 1.0, 1e-12);
}

}  // namespace
