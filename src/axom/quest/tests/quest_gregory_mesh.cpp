// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/GregoryMesh.hpp"

#include "axom/core/ArrayView.hpp"
#include "axom/slic.hpp"

#include "conduit_blueprint.hpp"
#include "conduit_node.hpp"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <variant>
#include <vector>

namespace
{
namespace quest = axom::quest;
using PointType = quest::GregoryMesh::PointType;
using VectorType = quest::GregoryMesh::VectorType;

// Create the explicit 3D coordset shared by the test topologies.
void setCoords(conduit::Node& mesh,
               const std::vector<double>& x,
               const std::vector<double>& y,
               const std::vector<double>& z)
{
  mesh["coordsets/coords/type"] = "explicit";
  mesh["coordsets/coords/values/x"].set(x);
  mesh["coordsets/coords/values/y"].set(y);
  mesh["coordsets/coords/values/z"].set(z);
}

// Add a fixed-shape topology whose connectivity indexes the shared coordset.
void setTopology(conduit::Node& mesh,
                 const std::string& name,
                 const std::string& shape,
                 const std::vector<int>& connectivity)
{
  conduit::Node& topology = mesh["topologies/" + name];
  topology["type"] = "unstructured";
  topology["coordset"] = "coords";
  topology["elements/shape"] = shape;
  topology["elements/connectivity"].set(connectivity);
}

// Compare every constructed control point, not just the input corner positions.
void expectTriangleGeometry(const quest::GregoryFacet& facet,
                            const std::array<PointType, 3>& points,
                            const std::array<VectorType, 3>& normals,
                            const char* context)
{
  SCOPED_TRACE(context);
  const axom::ArrayView<const PointType> coords(points.data(), 3);
  const axom::ArrayView<const VectorType> normal_view(normals.data(), 3);
  const axom::primal::GregoryTriangle<double> expected(coords, normal_view);
  const auto& actual = std::get<axom::primal::GregoryTriangle<double>>(facet.geometry);
  for(int i = 0; i < actual.NPTS; ++i)
  {
    for(int dimension = 0; dimension < 3; ++dimension)
    {
      EXPECT_NEAR(actual.getControlPoints()[i][dimension], expected.getControlPoints()[i][dimension], 1.e-12)
        << "control point " << i << ", dimension " << dimension;
    }
  }
}

TEST(quest_gregory_mesh, preserves_named_mixed_topologies_and_node_ids)
{
  conduit::Node mesh;
  setCoords(mesh,
            {0., 1., 0., 2., 2., 3., 4., 3., 4.},
            {0., 0., 1., 0., 1., 0., 0., 1., 1.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0.});

  // Set the mixed topology directly.
  conduit::Node& panel = mesh["topologies/panel"];
  panel["type"] = "unstructured";
  panel["coordset"] = "coords";
  panel["elements/shape"] = "mixed";
  panel["elements/shape_map/tri"] = 5;
  panel["elements/shape_map/quad"] = 9;
  panel["elements/shapes"].set(std::vector<int> {5, 9});
  panel["elements/sizes"].set(std::vector<int> {3, 4});
  panel["elements/offsets"].set(std::vector<int> {0, 3});
  panel["elements/connectivity"].set(std::vector<int> {0, 1, 2, 1, 3, 4, 2});

  // Set the quad topology with the helper
  setTopology(mesh, "island", "quad", {5, 6, 8, 7});

  conduit::Node info;
  ASSERT_TRUE(conduit::blueprint::mesh::verify(mesh, info)) << info.to_yaml();
  const quest::GregoryMesh gregory_mesh(mesh);

  // Surface and facet order, shapes, and node IDs should follow Blueprint connectivity.
  const auto& surfaces = gregory_mesh.getSurfaces();
  ASSERT_EQ(surfaces.size(), 2);
  EXPECT_EQ(surfaces[0].name, "panel");
  EXPECT_EQ(surfaces[1].name, "island");
  ASSERT_EQ(surfaces[0].facets.size(), 2);
  ASSERT_EQ(surfaces[1].facets.size(), 1);

  // Check the facet geometry
  const auto& tri = surfaces[0].facets[0];
  const auto& quad = surfaces[0].facets[1];
  EXPECT_EQ(tri.shape, quest::GregoryFacetShape::Tri);
  EXPECT_EQ(quad.shape, quest::GregoryFacetShape::Quad);
  EXPECT_EQ(tri.node_ids[0], 0);
  EXPECT_EQ(tri.node_ids[1], 1);
  EXPECT_EQ(tri.node_ids[2], 2);
  EXPECT_EQ(quad.node_ids[0], 1);
  EXPECT_EQ(quad.node_ids[1], 3);
  EXPECT_EQ(quad.node_ids[2], 4);
  EXPECT_EQ(quad.node_ids[3], 2);

  // Check the facet ordering
  EXPECT_EQ(surfaces[1].facets[0].shape, quest::GregoryFacetShape::Quad);
  EXPECT_EQ(surfaces[1].facets[0].node_ids[0], 5);
  EXPECT_EQ(surfaces[1].facets[0].node_ids[1], 6);
  EXPECT_EQ(surfaces[1].facets[0].node_ids[2], 8);
  EXPECT_EQ(surfaces[1].facets[0].node_ids[3], 7);

  // Check the variant status
  EXPECT_TRUE(std::holds_alternative<axom::primal::GregoryTriangle<double>>(tri.geometry));
  EXPECT_TRUE(std::holds_alternative<axom::primal::GregoryPatch<double>>(quad.geometry));
}

TEST(quest_gregory_mesh, facet_evaluation_and_bounds_dispatch)
{
  // Exercise the shape-dispatching facet API on one triangle and one quad.
  conduit::Node mesh;
  setCoords(mesh,
            {0., 1., 0., 2., 3., 3., 2.},
            {0., 0., 1., 0., 0., 1., 1.},
            {0., 0., 0., 0., 0., 0., 0.});
  setTopology(mesh, "triangle", "tri", {0, 1, 2});
  setTopology(mesh, "quadrilateral", "quad", {3, 4, 5, 6});

  // Check surface size and order
  const quest::GregoryMesh gregory_mesh(mesh);
  const auto& surfaces = gregory_mesh.getSurfaces();
  ASSERT_EQ(surfaces.size(), 2);
  ASSERT_EQ(surfaces[0].facets.size(), 1);
  ASSERT_EQ(surfaces[1].facets.size(), 1);

  // Check surface and facet geometry
  const auto& triangle = surfaces[0].facets[0];
  const auto& quadrilateral = surfaces[1].facets[0];
  EXPECT_EQ(triangle.numVertices(), 3);
  EXPECT_EQ(quadrilateral.numVertices(), 4);

  // Verify `evaluate` and `boundingBox` methods for the facet variant
  EXPECT_TRUE(triangle.boundingBox().contains(PointType {0., 0., 0.}));
  EXPECT_TRUE(quadrilateral.boundingBox().contains(PointType {3., 1., 0.}));
  EXPECT_TRUE(triangle.evaluate(0., 0.).isNearlyEqual(PointType {0., 0., 0.}, 1.e-12));
  EXPECT_TRUE(quadrilateral.evaluate(1., 1.).isNearlyEqual(PointType {3., 1., 0.}, 1.e-12));
}

TEST(quest_gregory_mesh, averages_normals_at_shared_nodes)
{
  // Two triangles meet along an edge with perpendicular face normals.
  conduit::Node mesh;
  setCoords(mesh, {0., 1., 0., 0.}, {0., 0., 1., 0.}, {0., 0., 0., 1.});
  setTopology(mesh, "bend", "tri", {0, 1, 2, 0, 3, 1});

  const quest::GregoryMesh gregory_mesh(mesh);
  ASSERT_EQ(gregory_mesh.getSurfaces().size(), 1);
  const auto& facets = gregory_mesh.getSurfaces()[0].facets;
  ASSERT_EQ(facets.size(), 2);
  const double component = 1. / std::sqrt(2.);
  const PointType p0 {0., 0., 0.};
  const PointType p1 {1., 0., 0.};
  const PointType p2 {0., 1., 0.};
  const PointType p3 {0., 0., 1.};
  const VectorType average {0., component, component};
  const VectorType up {0., 0., 1.};
  const VectorType forward {0., 1., 0.};

  // Shared corners get the normalized sum; unshared corners keep the face normal.
  expectTriangleGeometry(facets[0], {p0, p1, p2}, {average, average, up}, "bend first facet");
  expectTriangleGeometry(facets[1], {p0, p3, p1}, {average, forward, average}, "bend second facet");
}

TEST(quest_gregory_mesh, keeps_normals_separate_across_surfaces)
{
  // Reuse node IDs across topologies; their normals must be averaged separately.
  conduit::Node mesh;
  setCoords(mesh, {0., 1., 0., 0.}, {0., 0., 1., 0.}, {0., 0., 0., 1.});
  setTopology(mesh, "first", "tri", {0, 1, 2});
  setTopology(mesh, "second", "tri", {0, 3, 1});

  const quest::GregoryMesh separate_mesh(mesh);
  const auto& surfaces = separate_mesh.getSurfaces();
  ASSERT_EQ(surfaces.size(), 2);
  ASSERT_EQ(surfaces[0].facets.size(), 1);
  ASSERT_EQ(surfaces[1].facets.size(), 1);
  const PointType p0 {0., 0., 0.};
  const PointType p1 {1., 0., 0.};
  const PointType p2 {0., 1., 0.};
  const PointType p3 {0., 0., 1.};
  const VectorType up {0., 0., 1.};
  const VectorType forward {0., 1., 0.};

  // Different surfaces don't share normals in the average
  expectTriangleGeometry(surfaces[0].facets[0],
                         {p0, p1, p2},
                         {up, up, up},
                         "separate first surface");
  expectTriangleGeometry(surfaces[1].facets[0],
                         {p0, p3, p1},
                         {forward, forward, forward},
                         "separate second surface");
}

TEST(quest_gregory_mesh, disconnected_topology_keeps_one_name)
{
  // Disconnected triangles in one topology still form one named surface.
  conduit::Node mesh;
  setCoords(mesh, {0., 1., 0., 3., 4., 3.}, {0., 0., 1., 0., 0., 1.}, {0., 0., 0., 0., 0., 0.});
  setTopology(mesh, "parts", "tri", {0, 1, 2, 3, 4, 5});

  quest::GregoryMesh gregory_mesh(mesh);
  ASSERT_EQ(gregory_mesh.getSurfaces().size(), 1);
  EXPECT_EQ(gregory_mesh.getSurfaces()[0].name, "parts");
  EXPECT_EQ(gregory_mesh.getSurfaces()[0].facets.size(), 2);
}

TEST(quest_gregory_mesh, skips_unsupported_topology_and_keeps_later_surfaces)
{
  conduit::Node mesh;
  // Put the unsupported line first to exercise selection of a later surface coordset.
  setCoords(mesh, {0., 1., 0., 1.}, {0., 0., 1., 1.}, {0., 0., 0., 0.});
  setTopology(mesh, "edge", "line", {0, 1});
  setTopology(mesh, "first", "tri", {0, 1, 2});
  setTopology(mesh, "second", "tri", {1, 3, 2});

  conduit::Node info;
  ASSERT_TRUE(conduit::blueprint::mesh::verify(mesh, info)) << info.to_yaml();
  quest::GregoryMesh gregory_mesh(mesh);
  const auto& surfaces = gregory_mesh.getSurfaces();
  ASSERT_EQ(surfaces.size(), 2);
  EXPECT_EQ(surfaces[0].name, "first");
  EXPECT_EQ(surfaces[1].name, "second");
}

TEST(quest_gregory_mesh, retains_degenerate_triangle_and_later_surface)
{
  // For now, a collinear triangle should remain present without hiding a later valid surface.
  conduit::Node mesh;
  setCoords(mesh, {0., 1., 2., 0., 1., 0.}, {0., 0., 0., 0., 0., 1.}, {0., 0., 0., 0., 0., 0.});
  setTopology(mesh, "flat", "tri", {0, 1, 2});
  setTopology(mesh, "good", "tri", {3, 4, 5});

  const quest::GregoryMesh gregory_mesh(mesh);
  const auto& surfaces = gregory_mesh.getSurfaces();
  ASSERT_EQ(surfaces.size(), 2);
  EXPECT_EQ(surfaces[0].name, "flat");
  ASSERT_EQ(surfaces[0].facets.size(), 1);
  EXPECT_EQ(surfaces[0].facets[0].shape, quest::GregoryFacetShape::Tri);
  EXPECT_EQ(surfaces[0].facets[0].node_ids[0], 0);
  EXPECT_EQ(surfaces[0].facets[0].node_ids[1], 1);
  EXPECT_EQ(surfaces[0].facets[0].node_ids[2], 2);
  EXPECT_EQ(surfaces[1].name, "good");
  EXPECT_EQ(surfaces[1].facets.size(), 1);
}

TEST(quest_gregory_mesh, retains_degenerate_triangle_with_shared_vertices)
{
  // Keep the collinear facet in source order when neighboring facets share its nodes.
  conduit::Node mesh;
  setCoords(mesh, {0., 1., 2., 0.}, {0., 0., 0., 1.}, {0., 0., 0., 0.});
  setTopology(mesh, "parts", "tri", {0, 1, 2, 0, 1, 3, 1, 2, 3});

  quest::GregoryMesh gregory_mesh(mesh);
  const auto& surfaces = gregory_mesh.getSurfaces();
  ASSERT_EQ(surfaces.size(), 1);
  ASSERT_EQ(surfaces[0].facets.size(), 3);
  const auto& flat = surfaces[0].facets[0];
  EXPECT_EQ(flat.shape, quest::GregoryFacetShape::Tri);
  EXPECT_EQ(flat.node_ids[0], 0);
  EXPECT_EQ(flat.node_ids[1], 1);
  EXPECT_EQ(flat.node_ids[2], 2);
}

}  // namespace
