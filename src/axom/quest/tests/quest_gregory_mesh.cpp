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
  setTopology(mesh, "island", "quad", {5, 6, 8, 7});

  conduit::Node info;
  ASSERT_TRUE(conduit::blueprint::mesh::verify(mesh, info)) << info.to_yaml();
  const quest::GregoryMesh gregory_mesh(mesh);

  const auto& surfaces = gregory_mesh.getSurfaces();
  ASSERT_EQ(surfaces.size(), 2);
  EXPECT_EQ(surfaces[0].name, "panel");
  EXPECT_EQ(surfaces[1].name, "island");
  ASSERT_EQ(surfaces[0].facets.size(), 2);
  ASSERT_EQ(surfaces[1].facets.size(), 1);
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
  EXPECT_EQ(surfaces[1].facets[0].shape, quest::GregoryFacetShape::Quad);
  EXPECT_EQ(surfaces[1].facets[0].node_ids[0], 5);
  EXPECT_EQ(surfaces[1].facets[0].node_ids[1], 6);
  EXPECT_EQ(surfaces[1].facets[0].node_ids[2], 8);
  EXPECT_EQ(surfaces[1].facets[0].node_ids[3], 7);
  EXPECT_TRUE(std::holds_alternative<axom::primal::GregoryTriangle<double>>(tri.geometry));
  EXPECT_TRUE(std::holds_alternative<axom::primal::GregoryPatch<double>>(quad.geometry));
}

TEST(quest_gregory_mesh, owns_geometry_after_blueprint_is_destroyed)
{
  std::unique_ptr<quest::GregoryMesh> gregory_mesh;
  {
    conduit::Node mesh;
    setCoords(mesh,
              {0., 1., 0., 2., 3., 3., 2.},
              {0., 0., 1., 0., 0., 1., 1.},
              {0., 0., 0., 0., 0., 0., 0.});
    setTopology(mesh, "triangle", "tri", {0, 1, 2});
    setTopology(mesh, "quadrilateral", "quad", {3, 4, 5, 6});
    gregory_mesh = std::make_unique<quest::GregoryMesh>(mesh);
  }

  const auto& surfaces = gregory_mesh->getSurfaces();
  ASSERT_EQ(surfaces.size(), 2);
  ASSERT_EQ(surfaces[0].facets.size(), 1);
  ASSERT_EQ(surfaces[1].facets.size(), 1);
  const auto& triangle =
    std::get<axom::primal::GregoryTriangle<double>>(surfaces[0].facets[0].geometry);
  const auto& quadrilateral =
    std::get<axom::primal::GregoryPatch<double>>(surfaces[1].facets[0].geometry);
  EXPECT_EQ(triangle.getCorner(0), (PointType {0., 0., 0.}));
  EXPECT_EQ(triangle.getCorner(1), (PointType {1., 0., 0.}));
  EXPECT_EQ(triangle.getCorner(2), (PointType {0., 1., 0.}));
  EXPECT_EQ(quadrilateral.getCorner(0), (PointType {2., 0., 0.}));
  EXPECT_EQ(quadrilateral.getCorner(2), (PointType {3., 1., 0.}));
}

TEST(quest_gregory_mesh, facet_evaluation_and_bounds_dispatch)
{
  conduit::Node mesh;
  setCoords(mesh,
            {0., 1., 0., 2., 3., 3., 2.},
            {0., 0., 1., 0., 0., 1., 1.},
            {0., 0., 0., 0., 0., 0., 0.});
  setTopology(mesh, "triangle", "tri", {0, 1, 2});
  setTopology(mesh, "quadrilateral", "quad", {3, 4, 5, 6});

  const quest::GregoryMesh gregory_mesh(mesh);
  const auto& surfaces = gregory_mesh.getSurfaces();
  ASSERT_EQ(surfaces.size(), 2);
  ASSERT_EQ(surfaces[0].facets.size(), 1);
  ASSERT_EQ(surfaces[1].facets.size(), 1);
  const auto& triangle = surfaces[0].facets[0];
  const auto& quadrilateral = surfaces[1].facets[0];
  EXPECT_EQ(triangle.numVertices(), 3);
  EXPECT_EQ(quadrilateral.numVertices(), 4);
  EXPECT_TRUE(triangle.boundingBox().contains(PointType {0., 0., 0.}));
  EXPECT_TRUE(quadrilateral.boundingBox().contains(PointType {3., 1., 0.}));
  EXPECT_TRUE(triangle.evaluate(0., 0.).isNearlyEqual(PointType {0., 0., 0.}, 1.e-12));
  EXPECT_TRUE(quadrilateral.evaluate(1., 1.).isNearlyEqual(PointType {3., 1., 0.}, 1.e-12));
}

TEST(quest_gregory_mesh, averages_normals_at_shared_nodes)
{
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
  expectTriangleGeometry(facets[0], {p0, p1, p2}, {average, average, up}, "bend first facet");
  expectTriangleGeometry(facets[1], {p0, p3, p1}, {average, forward, average}, "bend second facet");
}

TEST(quest_gregory_mesh, keeps_normals_separate_across_surfaces)
{
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
  axom::slic::SimpleLogger logger;
  conduit::Node mesh;
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

TEST(quest_gregory_mesh, vtk_output_triangulates_facets_and_identifies_sources)
{
  conduit::Node mesh;
  setCoords(mesh,
            {0., 1., 0., 1., 2., 3., 3., 2.},
            {0., 0., 1., 1., 0., 0., 1., 1.},
            {0., 0., 0., 0., 0., 0., 0., 0.});
  setTopology(mesh, "triangles", "tri", {0, 1, 2, 1, 3, 2});
  setTopology(mesh, "quadrilateral", "quad", {4, 5, 6, 7});

  const quest::GregoryMesh gregory_mesh(mesh);
  const std::string file_path = "quest_gregory_mesh_output.vtk";
  gregory_mesh.writeVtkFile(file_path, 2);

  std::ifstream input(file_path);
  ASSERT_TRUE(input.good());
  const std::string vtk_contents(std::istreambuf_iterator<char> {input},
                                 std::istreambuf_iterator<char> {});
  input.close();
  EXPECT_EQ(std::remove(file_path.c_str()), 0);

  EXPECT_NE(vtk_contents.find("POINTS 21 double\n"), std::string::npos);
  EXPECT_NE(vtk_contents.find("CELLS 16 64\n"), std::string::npos);
  EXPECT_NE(vtk_contents.find("CELL_DATA 16\n"), std::string::npos);
  EXPECT_NE(vtk_contents.find("SCALARS surface_id int 1\n"), std::string::npos);
  EXPECT_NE(vtk_contents.find("SCALARS facet_id int 1\n"), std::string::npos);

  const std::string cell_types_header = "CELL_TYPES 16\n";
  const auto cell_types_start = vtk_contents.find(cell_types_header);
  ASSERT_NE(cell_types_start, std::string::npos);
  std::istringstream cell_types(vtk_contents.substr(cell_types_start + cell_types_header.size()));
  for(int cell_id = 0; cell_id < 16; ++cell_id)
  {
    int cell_type = -1;
    ASSERT_TRUE(static_cast<bool>(cell_types >> cell_type));
    EXPECT_EQ(cell_type, 5);  // VTK_TRIANGLE
  }

  const std::string surface_field_header = "SCALARS surface_id int 1\nLOOKUP_TABLE default\n";
  const auto surface_field_start = vtk_contents.find(surface_field_header);
  ASSERT_NE(surface_field_start, std::string::npos);
  std::istringstream surface_values(
    vtk_contents.substr(surface_field_start + surface_field_header.size()));
  for(int cell_id = 0; cell_id < 16; ++cell_id)
  {
    int surface_id = -1;
    ASSERT_TRUE(static_cast<bool>(surface_values >> surface_id));
    EXPECT_EQ(surface_id, cell_id < 8 ? 0 : 1);
  }

  const std::string facet_field_header = "SCALARS facet_id int 1\nLOOKUP_TABLE default\n";
  const auto facet_field_start = vtk_contents.find(facet_field_header);
  ASSERT_NE(facet_field_start, std::string::npos);
  std::istringstream facet_values(vtk_contents.substr(facet_field_start + facet_field_header.size()));
  for(int cell_id = 0; cell_id < 16; ++cell_id)
  {
    int facet_id = -1;
    ASSERT_TRUE(static_cast<bool>(facet_values >> facet_id));
    EXPECT_EQ(facet_id, cell_id < 8 ? cell_id / 4 : 2);
  }
}

TEST(quest_gregory_mesh, vtk_output_rejects_nonpositive_subdivisions)
{
  axom::slic::SimpleLogger logger;
  conduit::Node mesh;
  setCoords(mesh, {0., 1., 0.}, {0., 0., 1.}, {0., 0., 0.});
  setTopology(mesh, "triangle", "tri", {0, 1, 2});

  const quest::GregoryMesh gregory_mesh(mesh);
  axom::slic::ScopedAbortToThrow abort_guard;
  EXPECT_THROW(gregory_mesh.writeVtkFile("quest_gregory_mesh_invalid.vtk", 0),
               axom::slic::SlicAbortException);
}

TEST(quest_gregory_mesh, vtk_output_reports_file_open_failure)
{
  axom::slic::SimpleLogger logger;
  conduit::Node mesh;
  setCoords(mesh, {0., 1., 0.}, {0., 0., 1.}, {0., 0., 0.});
  setTopology(mesh, "triangle", "tri", {0, 1, 2});

  const quest::GregoryMesh gregory_mesh(mesh);
  axom::slic::ScopedAbortToThrow abort_guard;
  EXPECT_THROW(gregory_mesh.writeVtkFile("", 1), axom::slic::SlicAbortException);
}

TEST(quest_gregory_mesh, invalid_blueprint_returns_empty_mesh_when_abort_disabled)
{
  axom::slic::SimpleLogger logger;
  conduit::Node mesh;
  setCoords(mesh, {0., 1., 0.}, {0., 0., 1.}, {0., 0., 0.});
  setTopology(mesh, "good", "tri", {0, 1, 2});

  mesh["coordsets/coords/values/z"].reset();
  const bool abort_was_enabled = axom::slic::isAbortOnErrorsEnabled();
  axom::slic::disableAbortOnError();
  quest::GregoryMesh invalid_without_abort(mesh);
  axom::slic::setAbortOnError(abort_was_enabled);
  EXPECT_TRUE(invalid_without_abort.getSurfaces().empty());
}

TEST(quest_gregory_mesh, invalid_blueprint_reports_slic_error)
{
  axom::slic::SimpleLogger logger;
  conduit::Node mesh;
  setCoords(mesh, {0., 1., 0.}, {0., 0., 1.}, {0., 0., 0.});
  setTopology(mesh, "good", "tri", {0, 1, 2});
  mesh["coordsets/coords/values/z"].reset();

  axom::slic::ScopedAbortToThrow abort_guard;
  EXPECT_THROW(quest::GregoryMesh invalid(mesh), axom::slic::SlicAbortException);
}

TEST(quest_gregory_mesh, empty_mesh_reports_slic_error)
{
  axom::slic::SimpleLogger logger;
  conduit::Node mesh;
  axom::slic::ScopedAbortToThrow abort_guard;
  EXPECT_THROW(quest::GregoryMesh gregory_mesh(mesh), axom::slic::SlicAbortException);
}

TEST(quest_gregory_mesh, unsupported_only_topology_reports_slic_error)
{
  axom::slic::SimpleLogger logger;
  conduit::Node mesh;
  setCoords(mesh, {0., 1.}, {0., 0.}, {0., 0.});
  setTopology(mesh, "edge", "line", {0, 1});

  conduit::Node info;
  ASSERT_TRUE(conduit::blueprint::mesh::verify(mesh, info)) << info.to_yaml();
  axom::slic::ScopedAbortToThrow abort_guard;
  EXPECT_THROW(quest::GregoryMesh gregory_mesh(mesh), axom::slic::SlicAbortException);
}

TEST(quest_gregory_mesh, nonexplicit_coordinates_report_slic_error)
{
  axom::slic::SimpleLogger logger;
  conduit::Node mesh;
  mesh["coordsets/coords/type"] = "uniform";
  mesh["coordsets/coords/dims/i"] = 2;
  mesh["coordsets/coords/dims/j"] = 2;
  mesh["coordsets/coords/dims/k"] = 2;
  setTopology(mesh, "surface", "tri", {0, 1, 2});

  conduit::Node info;
  ASSERT_TRUE(conduit::blueprint::mesh::verify(mesh, info)) << info.to_yaml();
  axom::slic::ScopedAbortToThrow abort_guard;
  EXPECT_THROW(quest::GregoryMesh gregory_mesh(mesh), axom::slic::SlicAbortException);
}

}  // namespace
