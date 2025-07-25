// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/core/NumericLimits.hpp"
#include "axom/mint/utils/vtk_utils.hpp"
#include "axom/quest/io/ProEReader.hpp"
#include "axom/slic.hpp"

// gtest includes
#include "gtest/gtest.h"

#include <fstream>
#include <iostream>

namespace fs = axom::utilities::filesystem;

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
/*!
 * \brief Generates a Pro/E file consisting of a single tetrahedron
 * \param [in] file the temp file to generate.
 * \pre file.empty() == false
 */
void generate_pro_e_file(fs::TempFile& file)
{
  EXPECT_FALSE(file.getPath().empty());

  file.open();
  EXPECT_TRUE(file.is_open());

  file << "# Comment header to ignore" << std::endl;
  file << "# Another comment" << std::endl;

  // Number of nodes followed by number of tetrahedra
  file << "4 1" << std::endl;

  // Node ID followed by xyz coordinates
  file << "1 -1.0 0.0 0.0" << std::endl;
  file << "2 1.0 0.0 0.0" << std::endl;
  file << "3 0.0 1.0 0.0" << std::endl;
  file << "4 0.0 0.0 1.0" << std::endl;

  // Tetrahedron ID followed by corresponding Node IDs
  file << "1 1 2 3 4" << std::endl;

  file.close();
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
TEST(quest_pro_e_reader, read_missing_file)
{
  const std::string INVALID_FILE = "nonexistent_file.proe";
  axom::quest::ProEReader reader;
  reader.setFileName(INVALID_FILE);
  int status = reader.read();
  EXPECT_TRUE(status != 0);
}

//------------------------------------------------------------------------------
TEST(quest_pro_e_reader, read_to_invalid_mesh)
{
  const char* IGNORE_OUTPUT = ".*";

  // STEP 0: generate a temporary Pro/E file for testing
  fs::TempFile testFile("tet", ".proe");
  generate_pro_e_file(testFile);

  // STEP 1: constructs mesh object to read in the mesh to
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> trimesh(2, axom::mint::TRIANGLE);
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> hexmesh(3, axom::mint::HEX);

  // STEP 2: read in the Pro/E mesh data
  axom::quest::ProEReader reader;
  reader.setFileName(testFile.getPath());
  int status = reader.read();
  EXPECT_EQ(status, 0);

  // STEP 3: death tests

  // read the Pro/E mesh data to a 2D axom::mint::Mesh should fail
  EXPECT_DEATH_IF_SUPPORTED(reader.getMesh(&trimesh), IGNORE_OUTPUT);

  // read the Pro/E mesh data to a axom::mint::Mesh that has a different cell type
  EXPECT_DEATH_IF_SUPPORTED(reader.getMesh(&hexmesh), IGNORE_OUTPUT);
}

//------------------------------------------------------------------------------
TEST(quest_pro_e_reader, read_pro_e)
{
  const double x_expected[] = {-1.0, 1.0, 0.0, 0.0};
  const double y_expected[] = {0.0, 0.0, 1.0, 0.0};
  const double z_expected[] = {0.0, 0.0, 0.0, 1.0};

  // STEP 0: generate a temporary Pro/E file for testing
  fs::TempFile testFile("tet", ".proe");
  generate_pro_e_file(testFile);

  // STEP 1: create an Pro/E reader and read-in the mesh data
  axom::quest::ProEReader reader;
  reader.setFileName(testFile.getPath());
  int status = reader.read();
  EXPECT_EQ(status, 0);

  // STEP 2: reading the Pro/E mesh data into a axom::mint::Mesh
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> mesh(3, axom::mint::TET);
  reader.getMesh(&mesh);

  // STEP 3: ensure the mesh is what is expected
  EXPECT_EQ(mesh.getNumberOfCells(), 1);
  EXPECT_EQ(mesh.getNumberOfNodes(), 4);

  const double* x = mesh.getCoordinateArray(axom::mint::X_COORDINATE);
  const double* y = mesh.getCoordinateArray(axom::mint::Y_COORDINATE);
  const double* z = mesh.getCoordinateArray(axom::mint::Z_COORDINATE);
  EXPECT_TRUE(x != nullptr);
  EXPECT_TRUE(y != nullptr);
  EXPECT_TRUE(z != nullptr);

  axom::IndexType numNodes = mesh.getNumberOfNodes();
  for(axom::IndexType inode = 0; inode < numNodes; ++inode)
  {
    EXPECT_NEAR(x[inode], x_expected[inode], axom::numeric_limits<double>::epsilon());
    EXPECT_NEAR(y[inode], y_expected[inode], axom::numeric_limits<double>::epsilon());
    EXPECT_NEAR(z[inode], z_expected[inode], axom::numeric_limits<double>::epsilon());
  }
}

//------------------------------------------------------------------------------
TEST(quest_pro_e_reader, read_pro_e_invbbox)
{
  const double x_expected[] = {-1.0, 1.0, 0.0, 0.0};
  const double y_expected[] = {0.0, 0.0, 1.0, 0.0};
  const double z_expected[] = {0.0, 0.0, 0.0, 1.0};

  // STEP 0: generate a temporary Pro/E file for testing
  fs::TempFile testFile("tet", ".proe");
  generate_pro_e_file(testFile);

  // STEP 1: create an Pro/E reader and read-in the mesh data
  axom::quest::ProEReader reader;
  // invalid bounding box is the same as no bounding box: keep everything
  axom::quest::ProEReader::BBox3D invbbox;
  reader.setTetPredFromBoundingBox(invbbox, false);
  reader.setFileName(testFile.getPath());
  int status = reader.read();
  EXPECT_EQ(status, 0);

  // STEP 2: reading the Pro/E mesh data into a axom::mint::Mesh
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> mesh(3, axom::mint::TET);
  reader.getMesh(&mesh);

  // STEP 3: ensure the mesh is what is expected
  EXPECT_EQ(mesh.getNumberOfCells(), 1);
  EXPECT_EQ(mesh.getNumberOfNodes(), 4);

  const double* x = mesh.getCoordinateArray(axom::mint::X_COORDINATE);
  const double* y = mesh.getCoordinateArray(axom::mint::Y_COORDINATE);
  const double* z = mesh.getCoordinateArray(axom::mint::Z_COORDINATE);
  EXPECT_TRUE(x != nullptr);
  EXPECT_TRUE(y != nullptr);
  EXPECT_TRUE(z != nullptr);

  axom::IndexType numNodes = mesh.getNumberOfNodes();
  for(axom::IndexType inode = 0; inode < numNodes; ++inode)
  {
    EXPECT_NEAR(x[inode], x_expected[inode], axom::numeric_limits<double>::epsilon());
    EXPECT_NEAR(y[inode], y_expected[inode], axom::numeric_limits<double>::epsilon());
    EXPECT_NEAR(z[inode], z_expected[inode], axom::numeric_limits<double>::epsilon());
  }
}

//------------------------------------------------------------------------------
TEST(quest_pro_e_reader, read_pro_e_bbox_all)
{
  const double x_expected[] = {-1.0, 1.0, 0.0, 0.0};
  const double y_expected[] = {0.0, 0.0, 1.0, 0.0};
  const double z_expected[] = {0.0, 0.0, 0.0, 1.0};

  // STEP 0: generate a temporary Pro/E file for testing
  fs::TempFile testFile("tet", ".proe");
  generate_pro_e_file(testFile);

  // STEP 1: create an Pro/E reader and read-in the mesh data
  axom::quest::ProEReader reader;
  // A bounding box that catches all the points
  axom::quest::ProEReader::BBox3D bbox;
  bbox.addPoint(axom::quest::ProEReader::Point3D {-1.5, -0.5, -0.5});
  bbox.addPoint(axom::quest::ProEReader::Point3D {1.5, 1.5, 1.5});
  reader.setTetPredFromBoundingBox(bbox, false);
  reader.setFileName(testFile.getPath());
  int status = reader.read();
  EXPECT_EQ(status, 0);

  // STEP 2: reading the Pro/E mesh data into a axom::mint::Mesh
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> mesh(3, axom::mint::TET);
  reader.getMesh(&mesh);

  // STEP 3: ensure the mesh is what is expected
  EXPECT_EQ(mesh.getNumberOfCells(), 1);
  EXPECT_EQ(mesh.getNumberOfNodes(), 4);

  const double* x = mesh.getCoordinateArray(axom::mint::X_COORDINATE);
  const double* y = mesh.getCoordinateArray(axom::mint::Y_COORDINATE);
  const double* z = mesh.getCoordinateArray(axom::mint::Z_COORDINATE);
  EXPECT_TRUE(x != nullptr);
  EXPECT_TRUE(y != nullptr);
  EXPECT_TRUE(z != nullptr);

  axom::IndexType numNodes = mesh.getNumberOfNodes();
  for(axom::IndexType inode = 0; inode < numNodes; ++inode)
  {
    EXPECT_NEAR(x[inode], x_expected[inode], axom::numeric_limits<double>::epsilon());
    EXPECT_NEAR(y[inode], y_expected[inode], axom::numeric_limits<double>::epsilon());
    EXPECT_NEAR(z[inode], z_expected[inode], axom::numeric_limits<double>::epsilon());
  }
}

//------------------------------------------------------------------------------
TEST(quest_pro_e_reader, read_pro_e_bbox_some)
{
  const double x_expected[] = {-1.0, 1.0, 0.0, 0.0};
  const double y_expected[] = {0.0, 0.0, 1.0, 0.0};
  const double z_expected[] = {0.0, 0.0, 0.0, 1.0};

  // STEP 0: generate a temporary Pro/E file for testing
  fs::TempFile testFile("tet", ".proe");
  generate_pro_e_file(testFile);

  // STEP 1: create an Pro/E reader and read-in the mesh data
  axom::quest::ProEReader reader;
  // A bounding box that catches some of the points, so the single
  // tet does not get added.
  axom::quest::ProEReader::BBox3D bbox;
  bbox.addPoint(axom::quest::ProEReader::Point3D {-1.5, -0.5, -0.5});
  bbox.addPoint(axom::quest::ProEReader::Point3D {0, 1.5, 1.5});
  reader.setTetPredFromBoundingBox(bbox, false);
  reader.setFileName(testFile.getPath());
  int status = reader.read();
  EXPECT_EQ(status, 0);

  // STEP 2: reading the Pro/E mesh data into a axom::mint::Mesh
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> mesh(3, axom::mint::TET);
  reader.getMesh(&mesh);

  // STEP 3: ensure the mesh is what is expected
  EXPECT_EQ(mesh.getNumberOfCells(), 0);
  EXPECT_EQ(mesh.getNumberOfNodes(), 4);

  const double* x = mesh.getCoordinateArray(axom::mint::X_COORDINATE);
  const double* y = mesh.getCoordinateArray(axom::mint::Y_COORDINATE);
  const double* z = mesh.getCoordinateArray(axom::mint::Z_COORDINATE);
  EXPECT_TRUE(x != nullptr);
  EXPECT_TRUE(y != nullptr);
  EXPECT_TRUE(z != nullptr);

  axom::IndexType numNodes = mesh.getNumberOfNodes();
  for(axom::IndexType inode = 0; inode < numNodes; ++inode)
  {
    EXPECT_NEAR(x[inode], x_expected[inode], axom::numeric_limits<double>::epsilon());
    EXPECT_NEAR(y[inode], y_expected[inode], axom::numeric_limits<double>::epsilon());
    EXPECT_NEAR(z[inode], z_expected[inode], axom::numeric_limits<double>::epsilon());
  }
}

//------------------------------------------------------------------------------
TEST(quest_pro_e_reader, read_pro_e_bbox_some_incl)
{
  const double x_expected[] = {-1.0, 1.0, 0.0, 0.0};
  const double y_expected[] = {0.0, 0.0, 1.0, 0.0};
  const double z_expected[] = {0.0, 0.0, 0.0, 1.0};

  // STEP 0: generate a temporary Pro/E file for testing
  fs::TempFile testFile("tet", ".proe");
  generate_pro_e_file(testFile);

  // STEP 1: create an Pro/E reader and read-in the mesh data
  axom::quest::ProEReader reader;
  // A bounding box that catches some of the points.  We'll catch
  // all tets that have at least one corner in the bbox.
  axom::quest::ProEReader::BBox3D bbox;
  bbox.addPoint(axom::quest::ProEReader::Point3D {-1.5, -0.5, -0.5});
  bbox.addPoint(axom::quest::ProEReader::Point3D {0, 1.5, 1.5});
  reader.setTetPredFromBoundingBox(bbox);
  reader.setFileName(testFile.getPath());
  int status = reader.read();
  EXPECT_EQ(status, 0);

  // STEP 2: reading the Pro/E mesh data into a axom::mint::Mesh
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> mesh(3, axom::mint::TET);
  reader.getMesh(&mesh);

  // STEP 3: ensure the mesh is what is expected
  EXPECT_EQ(mesh.getNumberOfCells(), 1);
  EXPECT_EQ(mesh.getNumberOfNodes(), 4);

  const double* x = mesh.getCoordinateArray(axom::mint::X_COORDINATE);
  const double* y = mesh.getCoordinateArray(axom::mint::Y_COORDINATE);
  const double* z = mesh.getCoordinateArray(axom::mint::Z_COORDINATE);
  EXPECT_TRUE(x != nullptr);
  EXPECT_TRUE(y != nullptr);
  EXPECT_TRUE(z != nullptr);

  axom::IndexType numNodes = mesh.getNumberOfNodes();
  for(axom::IndexType inode = 0; inode < numNodes; ++inode)
  {
    EXPECT_NEAR(x[inode], x_expected[inode], axom::numeric_limits<double>::epsilon());
    EXPECT_NEAR(y[inode], y_expected[inode], axom::numeric_limits<double>::epsilon());
    EXPECT_NEAR(z[inode], z_expected[inode], axom::numeric_limits<double>::epsilon());
  }
}

//------------------------------------------------------------------------------
TEST(quest_pro_e_reader, read_pro_e_external)
{
  constexpr axom::IndexType N_NODES = 4;
  constexpr axom::IndexType N_TETS = 1;
  const double x_expected[] = {-1.0, 1.0, 0.0, 0.0};
  const double y_expected[] = {0.0, 0.0, 1.0, 0.0};
  const double z_expected[] = {0.0, 0.0, 0.0, 1.0};

  double xin[] = {-1.0, -1.0, -1.0, -1.0};
  double yin[] = {-1.0, -1.0, -1.0, -1.0};
  double zin[] = {-1.0, -1.0, -1.0, -1.0};

  axom::IndexType conn[] = {-1, -1, -1, -1};

  // STEP 0: generate a temporary Pro/E file for testing
  fs::TempFile testFile("tet", ".proe");
  generate_pro_e_file(testFile);

  // STEP 1: create a Pro/E reader and read-in the mesh data
  axom::quest::ProEReader reader;
  reader.setFileName(testFile.getPath());
  int status = reader.read();
  EXPECT_EQ(status, 0);

  // STEP 2: reading the Pro/E mesh data into a mint::Mesh
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>
    mesh(axom::mint::TET, N_TETS, conn, N_NODES, xin, yin, zin);
  EXPECT_EQ(mesh.getNumberOfCells(), N_TETS);
  EXPECT_EQ(mesh.getNumberOfNodes(), N_NODES);

  reader.getMesh(&mesh);

  const double* x = mesh.getCoordinateArray(axom::mint::X_COORDINATE);
  const double* y = mesh.getCoordinateArray(axom::mint::Y_COORDINATE);
  const double* z = mesh.getCoordinateArray(axom::mint::Z_COORDINATE);
  EXPECT_TRUE(x != nullptr);
  EXPECT_TRUE(y != nullptr);
  EXPECT_TRUE(z != nullptr);

  axom::IndexType numNodes = mesh.getNumberOfNodes();
  for(axom::IndexType inode = 0; inode < numNodes; ++inode)
  {
    EXPECT_NEAR(x[inode], x_expected[inode], axom::numeric_limits<double>::epsilon());
    EXPECT_NEAR(y[inode], y_expected[inode], axom::numeric_limits<double>::epsilon());
    EXPECT_NEAR(z[inode], z_expected[inode], axom::numeric_limits<double>::epsilon());
  }
}

#ifdef AXOM_DATA_DIR
//------------------------------------------------------------------------------
TEST(quest_pro_e_reader, cup_pro_e)
{
  constexpr int NUM_NODES = 171;
  constexpr int NUM_TETS = 574;
  constexpr double EPS = axom::numeric_limits<double>::epsilon();

  // STEP 0: Get Pro/E cup example file for testing
  std::string cup = fs::joinPath(AXOM_DATA_DIR, "quest/cup.proe");

  // STEP 1: create a Pro/E reader and read-in the mesh data
  axom::quest::ProEReader reader;
  reader.setFileName(cup);
  int status = reader.read();
  EXPECT_EQ(status, 0);

  // STEP 2: reading the Pro/E mesh data into a axom::mint::Mesh
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> mesh(3, axom::mint::TET);
  reader.getMesh(&mesh);

  // STEP 3: ensure the mesh is what is expected
  EXPECT_EQ(mesh.getNumberOfCells(), NUM_TETS);
  EXPECT_EQ(mesh.getNumberOfNodes(), NUM_NODES);

  const double* x = mesh.getCoordinateArray(axom::mint::X_COORDINATE);
  const double* y = mesh.getCoordinateArray(axom::mint::Y_COORDINATE);
  const double* z = mesh.getCoordinateArray(axom::mint::Z_COORDINATE);
  EXPECT_TRUE(x != nullptr);
  EXPECT_TRUE(y != nullptr);
  EXPECT_TRUE(z != nullptr);

  // STEP 4: Verify a few tetrahedra are as expected

  // Node IDs for the tetrahedron
  axom::IndexType nodeIDs[4];

  // Check first tetrahedron
  mesh.getCellNodeIDs(0, nodeIDs);

  // Note: Node IDs adjusted to start at 0
  EXPECT_EQ(nodeIDs[0], 10);
  EXPECT_EQ(nodeIDs[1], 91);
  EXPECT_EQ(nodeIDs[2], 90);
  EXPECT_EQ(nodeIDs[3], 88);

  EXPECT_NEAR(x[nodeIDs[0]], -71.842646, EPS);
  EXPECT_NEAR(x[nodeIDs[1]], 18.091620, EPS);
  EXPECT_NEAR(x[nodeIDs[2]], 18.091620, EPS);
  EXPECT_NEAR(x[nodeIDs[3]], 16.637112, EPS);

  EXPECT_NEAR(y[nodeIDs[0]], -142.432949, EPS);
  EXPECT_NEAR(y[nodeIDs[1]], -100.389772, EPS);
  EXPECT_NEAR(y[nodeIDs[2]], -158.496833, EPS);
  EXPECT_NEAR(y[nodeIDs[3]], -116.282419, EPS);

  EXPECT_NEAR(z[nodeIDs[0]], 0.000000, EPS);
  EXPECT_NEAR(z[nodeIDs[1]], 0.000000, EPS);
  EXPECT_NEAR(z[nodeIDs[2]], 0.000000, EPS);
  EXPECT_NEAR(z[nodeIDs[3]], -53.197840, EPS);

  // Check last tetrahedron
  mesh.getCellNodeIDs(NUM_TETS - 1, nodeIDs);

  EXPECT_EQ(nodeIDs[0], 170);
  EXPECT_EQ(nodeIDs[1], 84);
  EXPECT_EQ(nodeIDs[2], 143);
  EXPECT_EQ(nodeIDs[3], 83);

  EXPECT_NEAR(x[nodeIDs[0]], 128.940529, EPS);
  EXPECT_NEAR(x[nodeIDs[1]], 68.261330, EPS);
  EXPECT_NEAR(x[nodeIDs[2]], 121.442679, EPS);
  EXPECT_NEAR(x[nodeIDs[3]], 68.137155, EPS);

  EXPECT_NEAR(y[nodeIDs[0]], -93.930323, EPS);
  EXPECT_NEAR(y[nodeIDs[1]], -41.261753, EPS);
  EXPECT_NEAR(y[nodeIDs[2]], -31.346653, EPS);
  EXPECT_NEAR(y[nodeIDs[3]], -41.466486, EPS);

  EXPECT_NEAR(z[nodeIDs[0]], -145.985725, EPS);
  EXPECT_NEAR(z[nodeIDs[1]], -122.672539, EPS);
  EXPECT_NEAR(z[nodeIDs[2]], -174.327849, EPS);
  EXPECT_NEAR(z[nodeIDs[3]], -173.363017, EPS);

  // Check middle tetrahedron
  mesh.getCellNodeIDs((NUM_TETS / 2) - 1, nodeIDs);

  EXPECT_EQ(nodeIDs[0], 117);
  EXPECT_EQ(nodeIDs[1], 118);
  EXPECT_EQ(nodeIDs[2], 114);
  EXPECT_EQ(nodeIDs[3], 115);

  EXPECT_NEAR(x[nodeIDs[0]], 15.508351, EPS);
  EXPECT_NEAR(x[nodeIDs[1]], -30.194340, EPS);
  EXPECT_NEAR(x[nodeIDs[2]], -35.211560, EPS);
  EXPECT_NEAR(x[nodeIDs[3]], 14.765131, EPS);

  EXPECT_NEAR(y[nodeIDs[0]], 124.048744, EPS);
  EXPECT_NEAR(y[nodeIDs[1]], 109.246291, EPS);
  EXPECT_NEAR(y[nodeIDs[2]], 117.715966, EPS);
  EXPECT_NEAR(y[nodeIDs[3]], 78.384516, EPS);

  EXPECT_NEAR(z[nodeIDs[0]], -179.442905, EPS);
  EXPECT_NEAR(z[nodeIDs[1]], -112.054428, EPS);
  EXPECT_NEAR(z[nodeIDs[2]], -178.748301, EPS);
  EXPECT_NEAR(z[nodeIDs[3]], -163.120569, EPS);

  // Step 5: Dump mesh file
  axom::mint::write_vtk(&mesh, "cup.vtk");
}

//------------------------------------------------------------------------------
TEST(quest_pro_e_reader, cup_pro_e_some)
{
  constexpr int NUM_NODES = 171;
  constexpr int NUM_BBOX_TETS = 52;

  // STEP 0: Get Pro/E cup example file for testing
  std::string cup = fs::joinPath(AXOM_DATA_DIR, "quest/cup.proe");

  // STEP 1: create a Pro/E reader and read-in the mesh data
  axom::quest::ProEReader reader;
  reader.setFileName(cup);
  // STEP 1a: specify a bounding box to include
  // A bounding box that catches some of the points: only tets that
  // fall completely within the bbox are retained.
  axom::quest::ProEReader::BBox3D bbox;
  bbox.addPoint(axom::quest::ProEReader::Point3D {-30, -160, -200});
  bbox.addPoint(axom::quest::ProEReader::Point3D {90, 160, -35});
  reader.setTetPredFromBoundingBox(bbox, false);
  int status = reader.read();
  EXPECT_EQ(status, 0);

  // STEP 2: reading the Pro/E mesh data into a axom::mint::Mesh
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> mesh(3, axom::mint::TET);
  reader.getMesh(&mesh);

  // STEP 3: ensure the mesh is what is expected
  EXPECT_EQ(mesh.getNumberOfCells(), NUM_BBOX_TETS);
  EXPECT_EQ(mesh.getNumberOfNodes(), NUM_NODES);

  // Step 5: Dump mesh file
  axom::mint::write_vtk(&mesh, "exclusivebboxcup.vtk");
}

//------------------------------------------------------------------------------
TEST(quest_pro_e_reader, cup_pro_e_some_incl)
{
  constexpr int NUM_NODES = 171;
  constexpr int NUM_BBOX_INCL_TETS = 298;

  // STEP 0: Get Pro/E cup example file for testing
  std::string cup = fs::joinPath(AXOM_DATA_DIR, "quest/cup.proe");

  // STEP 1: create a Pro/E reader and read-in the mesh data
  axom::quest::ProEReader reader;
  reader.setFileName(cup);
  // STEP 1a: specify a bounding box to include
  // A bounding box that catches some of the points: tets having at
  // least one corner falling within the bbox are retained.
  axom::quest::ProEReader::BBox3D bbox;
  bbox.addPoint(axom::quest::ProEReader::Point3D {-30, -160, -200});
  bbox.addPoint(axom::quest::ProEReader::Point3D {90, 160, -35});
  reader.setTetPredFromBoundingBox(bbox);
  int status = reader.read();
  EXPECT_EQ(status, 0);

  // STEP 2: reading the Pro/E mesh data into a axom::mint::Mesh
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> mesh(3, axom::mint::TET);
  reader.getMesh(&mesh);

  // STEP 3: ensure the mesh is what is expected
  EXPECT_EQ(mesh.getNumberOfCells(), NUM_BBOX_INCL_TETS);
  EXPECT_EQ(mesh.getNumberOfNodes(), NUM_NODES);

  // Step 5: Dump mesh file
  axom::mint::write_vtk(&mesh, "inclusivebboxcup.vtk");
}
#endif

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  return RUN_ALL_TESTS();
}
