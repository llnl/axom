// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Mint includes
#include "axom/mint/mesh/UniformMesh.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/mint/utils/su2_utils.hpp"

// Axom includes
#include "axom/core.hpp"
#include "axom/slic.hpp"

// gtest includes
#include "gtest/gtest.h"

namespace mint = axom::mint;

//------------------------------------------------------------------------------
// HELPER FUNCTIONS
//------------------------------------------------------------------------------
namespace
{
void check_cell_connectivity(mint::Mesh* m1, mint::Mesh* m2)
{
  EXPECT_TRUE(m1 != nullptr);
  EXPECT_TRUE(m2 != nullptr);

  const axom::IndexType ncells = m1->getNumberOfCells();
  EXPECT_EQ(ncells, m2->getNumberOfCells());

  axom::IndexType c1[mint::MAX_CELL_NODES];
  axom::IndexType c2[mint::MAX_CELL_NODES];

  for(axom::IndexType icell = 0; icell < ncells; ++icell)
  {
    EXPECT_TRUE(m1->getCellType(icell) == m2->getCellType(icell));

    axom::IndexType nnodes = m1->getNumberOfCellNodes(icell);
    EXPECT_EQ(nnodes, m2->getNumberOfCellNodes(icell));

    m1->getCellNodeIDs(icell, c1);
    m2->getCellNodeIDs(icell, c2);

    for(axom::IndexType inode = 0; inode < nnodes; ++inode)
    {
      EXPECT_EQ(c1[inode], c2[inode]);
    }  // END for all cell nodes

  }  // END for all cells
}

//------------------------------------------------------------------------------
void check_nodes(mint::Mesh* m1, mint::Mesh* m2)
{
  EXPECT_TRUE(m1 != nullptr);
  EXPECT_TRUE(m2 != nullptr);

  constexpr int NDIMS = 2;
  EXPECT_EQ(m1->getDimension(), NDIMS);
  EXPECT_EQ(m2->getDimension(), NDIMS);

  const axom::IndexType nnodes = m1->getNumberOfNodes();
  EXPECT_EQ(nnodes, m2->getNumberOfNodes());

  double n1[NDIMS];
  double n2[NDIMS];

  for(axom::IndexType inode = 0; inode < nnodes; ++inode)
  {
    m1->getNode(inode, n1);
    m2->getNode(inode, n2);

    for(int idim = 0; idim < NDIMS; ++idim)
    {
      EXPECT_DOUBLE_EQ(n1[idim], n2[idim]);
    }  // END for each dimension

  }  // END for all nodes
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(mint_su2_io, invalid_write_with_structured_mesh)
{
  const double lo[] = {0.0, 0.0, 0.0};
  const double hi[] = {2.0, 2.0, 2.0};
  constexpr axom::IndexType N = 5;

  mint::UniformMesh mesh(lo, hi, N, N, N);
  int rc = mint::write_su2(&mesh, "test.su2");
  EXPECT_TRUE(rc != 0);
}

//------------------------------------------------------------------------------
TEST(mint_su2_io, read_invalid_file)
{
  mint::Mesh* mesh = nullptr;
  int rc = mint::read_su2("some/invalid/file.su2", mesh);
  EXPECT_TRUE(rc != 0);
  EXPECT_TRUE(mesh == nullptr);
}

//------------------------------------------------------------------------------
TEST(mint_su2_io, write_read_mixed_cell_topology_mesh)
{
  constexpr int DIMENSION = 2;

  // Construct the mesh object
  mint::UnstructuredMesh<mint::MIXED_SHAPE> mesh(DIMENSION);

  // Append the mesh nodes
  const axom::IndexType n0 = mesh.appendNode(0.0, 0.0);
  const axom::IndexType n1 = mesh.appendNode(2.0, 0.0);
  const axom::IndexType n2 = mesh.appendNode(1.0, 1.0);
  const axom::IndexType n3 = mesh.appendNode(3.5, 1.0);
  const axom::IndexType n4 = mesh.appendNode(2.5, 2.0);
  const axom::IndexType n5 = mesh.appendNode(5.0, 0.0);

  // Append mesh cells
  const axom::IndexType c0[] = {n0, n1, n2};
  const axom::IndexType c1[] = {n1, n5, n3, n2};
  const axom::IndexType c2[] = {n3, n4, n2};

  mesh.appendCell(c0, mint::TRIANGLE);
  mesh.appendCell(c1, mint::QUAD);
  mesh.appendCell(c2, mint::TRIANGLE);

  // write an SU2 file
  axom::utilities::filesystem::TempFile tmpfile("mixed_cell_mesh", ".su2");
  const std::string& su2File = tmpfile.getPath();
  int rc = mint::write_su2(&mesh, su2File);
  EXPECT_EQ(rc, 0);

  // read it on a new mesh instance
  mint::Mesh* test_mesh = nullptr;
  rc = mint::read_su2(su2File, test_mesh);
  EXPECT_EQ(rc, 0);
  EXPECT_TRUE(test_mesh != nullptr);
  EXPECT_TRUE(test_mesh->isUnstructured());

  // test the mesh
  EXPECT_EQ(test_mesh->getDimension(), DIMENSION);
  EXPECT_EQ(test_mesh->getNumberOfNodes(), mesh.getNumberOfNodes());
  EXPECT_EQ(test_mesh->getNumberOfCells(), mesh.getNumberOfCells());
  EXPECT_EQ(test_mesh->hasMixedCellTypes(), mesh.hasMixedCellTypes());

  // test mesh nodes
  check_nodes(&mesh, test_mesh);

  // test mesh connectivity
  check_cell_connectivity(&mesh, test_mesh);

  // cleanup
  delete test_mesh;
}

//------------------------------------------------------------------------------
TEST(mint_su2_io, write_read_single_cell_topology_mesh)
{
  constexpr int DIMENSION = 2;
  constexpr mint::CellType CELL_TYPE = mint::TRIANGLE;

  // Construct the mesh object
  mint::UnstructuredMesh<mint::SINGLE_SHAPE> mesh(DIMENSION, CELL_TYPE);

  // Append the mesh nodes
  const axom::IndexType n0 = mesh.appendNode(0.0, 0.0);
  const axom::IndexType n1 = mesh.appendNode(2.0, 0.0);
  const axom::IndexType n2 = mesh.appendNode(1.0, 1.0);
  const axom::IndexType n3 = mesh.appendNode(3.5, 1.0);
  const axom::IndexType n4 = mesh.appendNode(2.5, 2.0);
  const axom::IndexType n5 = mesh.appendNode(5.0, 0.0);

  // Append mesh cells
  const axom::IndexType c0[] = {n1, n3, n2};
  const axom::IndexType c1[] = {n2, n0, n1};
  const axom::IndexType c2[] = {n3, n4, n2};
  const axom::IndexType c3[] = {n1, n5, n3};

  mesh.appendCell(c0);
  mesh.appendCell(c1);
  mesh.appendCell(c2);
  mesh.appendCell(c3);

  // write an SU2 file
  axom::utilities::filesystem::TempFile tempfile("simple_mesh", ".su2");
  const std::string& su2File = tempfile.getPath();
  int rc = mint::write_su2(&mesh, su2File);
  EXPECT_EQ(rc, 0);

  // read it on a new mesh instance
  mint::Mesh* test_mesh = nullptr;
  rc = mint::read_su2(su2File, test_mesh);
  EXPECT_EQ(rc, 0);
  EXPECT_TRUE(test_mesh != nullptr);
  EXPECT_TRUE(test_mesh->isUnstructured());

  // test the mesh
  EXPECT_EQ(test_mesh->getDimension(), DIMENSION);
  EXPECT_EQ(test_mesh->getNumberOfNodes(), mesh.getNumberOfNodes());
  EXPECT_EQ(test_mesh->getNumberOfCells(), mesh.getNumberOfCells());
  EXPECT_EQ(test_mesh->hasMixedCellTypes(), mesh.hasMixedCellTypes());

  // test mesh nodes
  check_nodes(&mesh, test_mesh);

  // test mesh connectivity
  check_cell_connectivity(&mesh, test_mesh);

  // cleanup
  delete test_mesh;
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();
  return result;
}
