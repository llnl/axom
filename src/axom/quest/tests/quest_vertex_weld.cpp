// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"
#include "axom/slic.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/quest/MeshTester.hpp"
#include "quest_test_utilities.hpp"

namespace
{
constexpr int DIM = 3;
constexpr double EPS = 1e-6;

using UMesh = axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>;
using Point3 = axom::primal::Point<double, 3>;

/*! Insert a vertex with coordinates (x,y,z) into \a mesh  */
void insertVertex(UMesh* mesh, double x, double y, double z) { mesh->appendNode(x, y, z); }

/*! Insert a triangle with vertex indices (v1,v2,v3) into \a mesh  */
void insertTriangle(UMesh* mesh, int v1, int v2, int v3)
{
  axom::IndexType indices[3] = {v1, v2, v3};
  mesh->appendCell(indices);
}

}  // namespace

//------------------------------------------------------------------------------
TEST(quest_vertex_weld, emptyMesh)
{
  SLIC_INFO("*** Tests welding function on an empty triangle mesh");

  UMesh* mesh = new UMesh(DIM, axom::mint::TRIANGLE);

  EXPECT_EQ(0, mesh->getNumberOfNodes());
  EXPECT_EQ(0, mesh->getNumberOfCells());

  axom::quest::weldTriMeshVertices(&mesh, EPS);

  EXPECT_EQ(0, mesh->getNumberOfNodes());
  EXPECT_EQ(0, mesh->getNumberOfCells());

  delete mesh;
  mesh = nullptr;
}

//------------------------------------------------------------------------------
TEST(quest_vertex_weld, onlyVertices)
{
  SLIC_INFO("*** Tests welding function on a triangle mesh" << " with vertices but no triangles.");

  UMesh* mesh = new UMesh(DIM, axom::mint::TRIANGLE);
  insertVertex(mesh, 1, 0, 0);
  insertVertex(mesh, 1, 1, 0);
  insertVertex(mesh, 1, 1, 1);

  EXPECT_EQ(3, mesh->getNumberOfNodes());
  EXPECT_EQ(0, mesh->getNumberOfCells());

  axom::quest::weldTriMeshVertices(&mesh, EPS);

  EXPECT_EQ(3, mesh->getNumberOfNodes());
  EXPECT_EQ(0, mesh->getNumberOfCells());

  delete mesh;
  mesh = nullptr;
}

//------------------------------------------------------------------------------
TEST(quest_vertex_weld, oneTriangle)
{
  SLIC_INFO("*** Tests welding function on a triangle mesh with one triangle.");

  UMesh* mesh = new UMesh(DIM, axom::mint::TRIANGLE);
  insertVertex(mesh, 1, 0, 0);
  insertVertex(mesh, 1, 1, 0);
  insertVertex(mesh, 1, 1, 1);

  EXPECT_EQ(3, mesh->getNumberOfNodes());

  insertTriangle(mesh, 0, 1, 2);
  EXPECT_EQ(1, mesh->getNumberOfCells());

  axom::quest::weldTriMeshVertices(&mesh, EPS);

  EXPECT_EQ(3, mesh->getNumberOfNodes());
  EXPECT_EQ(1, mesh->getNumberOfCells());

  delete mesh;
  mesh = nullptr;
}

//------------------------------------------------------------------------------
TEST(quest_vertex_weld, degenerateTriangle)
{
  SLIC_INFO("*** Tests welding function on a triangle mesh degenerate triangles.");

  UMesh* mesh = new UMesh(DIM, axom::mint::TRIANGLE);
  insertVertex(mesh, 1, 0, 0);
  insertVertex(mesh, 1, 1, 0);
  insertVertex(mesh, 1, 1, 1);
  insertVertex(mesh, 1.01, 1.01, 1.01);  // should be welded

  EXPECT_EQ(4, mesh->getNumberOfNodes());

  insertTriangle(mesh, 0, 1, 2);
  insertTriangle(mesh, 0, 1, 1);  // degenerate: only two distinct vertices
  insertTriangle(mesh, 1, 1, 1);  // degenerate: one one distinct vertex
  insertTriangle(mesh, 1, 2, 3);  // degenerate: verts 2 and 3 should be welded
  EXPECT_EQ(4, mesh->getNumberOfCells());

  const double eps = .1;
  axom::quest::weldTriMeshVertices(&mesh, eps);

  EXPECT_EQ(3, mesh->getNumberOfNodes());
  EXPECT_EQ(1, mesh->getNumberOfCells());

  delete mesh;
  mesh = nullptr;
}

//------------------------------------------------------------------------------
TEST(quest_vertex_weld, vertexAdjacentTrianglePair)
{
  SLIC_INFO("*** Tests welding function on a triangle mesh"
            << " with two triangles adjacent along a vertex.");

  UMesh* mesh = new UMesh(DIM, axom::mint::TRIANGLE);
  insertVertex(mesh, 1, 0, 0);
  insertVertex(mesh, 1, 1, 0);
  insertVertex(mesh, 1, 1, 1);

  insertVertex(mesh, 0, 1, 0);
  insertVertex(mesh, 0, 1, 1);
  insertVertex(mesh, 1, 1, 1);  // Duplicate

  EXPECT_EQ(6, mesh->getNumberOfNodes());

  insertTriangle(mesh, 0, 1, 2);
  insertTriangle(mesh, 3, 4, 5);
  EXPECT_EQ(2, mesh->getNumberOfCells());

  axom::quest::weldTriMeshVertices(&mesh, EPS);

  // Duplicated vertices should have been removed
  EXPECT_EQ(5, mesh->getNumberOfNodes());
  EXPECT_EQ(2, mesh->getNumberOfCells());

  delete mesh;
  mesh = nullptr;
}

//------------------------------------------------------------------------------
TEST(quest_vertex_weld, edgeAdjacentTrianglePair)
{
  SLIC_INFO("*** Tests welding function on a triangle mesh"
            << " with two triangles adjacent along an edge.");

  UMesh* mesh = new UMesh(DIM, axom::mint::TRIANGLE);
  insertVertex(mesh, 1, 0, 0);
  insertVertex(mesh, 1, 1, 0);
  insertVertex(mesh, 1, 1, 1);

  insertVertex(mesh, 0, 1, 0);
  insertVertex(mesh, 1, 1, 0);  // Duplicate
  insertVertex(mesh, 1, 1, 1);  // Duplicate

  EXPECT_EQ(6, mesh->getNumberOfNodes());

  insertTriangle(mesh, 0, 1, 2);
  insertTriangle(mesh, 3, 4, 5);
  EXPECT_EQ(2, mesh->getNumberOfCells());

  axom::quest::weldTriMeshVertices(&mesh, EPS);

  // Duplicated vertices should have been removed
  EXPECT_EQ(4, mesh->getNumberOfNodes());
  EXPECT_EQ(2, mesh->getNumberOfCells());

  delete mesh;
  mesh = nullptr;
}

//------------------------------------------------------------------------------
TEST(quest_vertex_weld, fuzzWeld)
{
  SLIC_INFO("*** Tests welding function on a triangle mesh"
            << " with a pair of triangles adjacent along an edge"
            << " whose vertices are welded");

  UMesh* mesh = new UMesh(DIM, axom::mint::TRIANGLE);
  insertVertex(mesh, 1, 0, 0);
  insertVertex(mesh, 1, 1, 0);
  insertVertex(mesh, 1, 1, 1);
  insertVertex(mesh, 0.99, 0.99, -0.01);  // should weld with 2nd vert
  insertVertex(mesh, 1.01, 1.01, 1.01);   // should weld with 3rd vert
  insertVertex(mesh, 0, 1, 0);

  EXPECT_EQ(6, mesh->getNumberOfNodes());

  insertTriangle(mesh, 0, 1, 2);
  insertTriangle(mesh, 3, 4, 5);
  EXPECT_EQ(2, mesh->getNumberOfCells());

  const double eps = .1;
  axom::quest::weldTriMeshVertices(&mesh, eps);

  EXPECT_EQ(4, mesh->getNumberOfNodes());
  EXPECT_EQ(2, mesh->getNumberOfCells());

  delete mesh;
  mesh = nullptr;
}

//------------------------------------------------------------------------------
TEST(quest_vertex_weld, disconnectedTrianglePair)
{
  SLIC_INFO("*** Tests welding function on a triangle mesh"
            << " with a pair of disconnected triangles.");

  UMesh* mesh = new UMesh(DIM, axom::mint::TRIANGLE);
  insertVertex(mesh, 1, 0, 0);  // verts for first triangle
  insertVertex(mesh, 1, 1, 0);
  insertVertex(mesh, 1, 1, 1);
  insertVertex(mesh, 0, -1, 0);  // verts for second triangle
  insertVertex(mesh, -1, -1, 0);
  insertVertex(mesh, -1, -1, -1);

  EXPECT_EQ(6, mesh->getNumberOfNodes());

  insertTriangle(mesh, 0, 1, 2);
  insertTriangle(mesh, 3, 4, 5);
  EXPECT_EQ(2, mesh->getNumberOfCells());

  axom::quest::weldTriMeshVertices(&mesh, EPS);

  // Duplicated vertices should have been removed
  EXPECT_EQ(6, mesh->getNumberOfNodes());
  EXPECT_EQ(2, mesh->getNumberOfCells());

  delete mesh;
  mesh = nullptr;
}

//------------------------------------------------------------------------------
TEST(quest_vertex_weld, indexedTetrahedron)
{
  SLIC_INFO("*** Tests welding function on an indexed tetrahedron");

  // Get the tetrahedron mesh.  It is already an indexed mesh.
  UMesh* mesh = static_cast<UMesh*>(axom::quest::utilities::make_tetrahedron_mesh());

  const int NV = 4;
  const int NT = 4;

  EXPECT_EQ(NV, mesh->getNumberOfNodes());
  EXPECT_EQ(NT, mesh->getNumberOfCells());

  axom::quest::weldTriMeshVertices(&mesh, EPS);

  // welding shouldn't change anything on this mesh
  EXPECT_EQ(NV, mesh->getNumberOfNodes());
  EXPECT_EQ(NT, mesh->getNumberOfCells());

  delete mesh;
  mesh = nullptr;
}

//------------------------------------------------------------------------------
TEST(quest_vertex_weld, indexedOctahedron)
{
  SLIC_INFO("*** Tests welding function on an indexed octahedron");

  // Get the octahedron mesh.  It is already an indexed mesh.
  UMesh* mesh = static_cast<UMesh*>(axom::quest::utilities::make_octahedron_mesh());

  const int NV = 6;
  const int NT = 8;

  EXPECT_EQ(NV, mesh->getNumberOfNodes());
  EXPECT_EQ(NT, mesh->getNumberOfCells());

  axom::quest::weldTriMeshVertices(&mesh, EPS);

  // welding shouldn't change anything on this mesh
  EXPECT_EQ(NV, mesh->getNumberOfNodes());
  EXPECT_EQ(NT, mesh->getNumberOfCells());

  delete mesh;
  mesh = nullptr;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
