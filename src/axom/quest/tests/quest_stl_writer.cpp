// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/io/STLWriter.hpp"
#include "axom/quest/io/STLReader.hpp"
#include "axom/mint/mesh/UniformMesh.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/slic.hpp"
#include "axom/fmt.hpp"

// gtest includes
#include "gtest/gtest.h"

// C/C++ includes
#include <cstdio>
#include <string>
#include <fstream>

// Uncomment this line to write new baseline data to stdout.
// #define AXOM_SET_BASELINES

// namespace aliases
namespace mint = axom::mint;
namespace quest = axom::quest;

namespace testing
{

void writeArray(const axom::Array<double>& vec, const std::string& var_name = "v")
{
  axom::fmt::print("axom::Array<double> {} = {{{}}};\n", var_name, axom::fmt::join(vec, ", "));
}

  /// Convert mesh coordinates into arrays that can be easily compared.
  void getCoordinates(const mint::Mesh &mesh, axom::Array<double> &xc, axom::Array<double> & yc)
  {
    for(axom::IndexType cellId = 0; cellId < mesh.getNumberOfCells(); cellId++)
    {
      EXPECT_EQ(mesh.getNumberOfCellNodes(cellId), 3);

      axom::IndexType nodes[3];
      mesh.getCellNodeIDs(cellId, nodes);
      for(int i = 0; i < 3; i++)
      {
        xc.push_back(mesh.getCoordinateArray(mint::X_COORDINATE)[nodes[i]]);
        yc.push_back(mesh.getCoordinateArray(mint::Y_COORDINATE)[nodes[i]]);
      }
    }
  }

bool compareArrays(const axom::Array<double> &A, const axom::Array<double> &B, double tolerance = 1.e-8)
{
  bool eq = A.size() == B.size();
  if(eq)
  {
    for(axom::IndexType i = 0; i < A.size() && eq; i++)
    {
      eq &= axom::utilities::isNearlyEqual(A[i], B[i], tolerance);
      if(!eq)
      {
        SLIC_ERROR(axom::fmt::format("Difference at index {}: {}, {}", i, A[i], B[i]));
      }
    }
  }
  return eq;
}

struct Test2D
{
  axom::Array<double> baselineXCoordinates()
  {
    return axom::Array<double>{{0, 0.5, 0.5, 0, 0.5, 0, 0.5, 1, 1, 0.5, 1, 0.5, 0, 0.5, 0.5, 0, 0.5, 0, 0.5, 1, 1, 0.5, 1, 0.5}};
  }

  axom::Array<double> baselineYCoordinates()
  {
    return axom::Array<double>{{1, 1, 1.5, 1, 1.5, 1.5, 1, 1, 1.5, 1, 1.5, 1.5, 1.5, 1.5, 2, 1.5, 2, 2, 1.5, 1.5, 2, 1.5, 2, 2}};
  }

  void test(const mint::Mesh &mesh, const std::string &filename, bool binary)
  {
    // Write STL file.
    int result = axom::quest::write_stl(&mesh, filename, binary);
    SLIC_INFO(axom::fmt::format("Writing {} -> {}", filename, (result == 0) ? "ok" : "error"));
    EXPECT_EQ(result, 0);

    // Read file back into memory and check the mesh.
    axom::quest::STLReader reader;
    reader.setFileName(filename);
    result = reader.read();
    SLIC_INFO(axom::fmt::format("Reading {} -> {}", filename, (result == 0) ? "ok" : "error"));
    EXPECT_EQ(result, 0);
    if(result == 0)
    {
//      std::cout << "reader.getNumNodes()=" << reader.getNumNodes() << std::endl;
//      std::cout << "reader.getNumFaces()=" << reader.getNumFaces() << std::endl;

      // Get the file data as a mint mesh.
      mint::UnstructuredMesh<mint::SINGLE_SHAPE> readMesh(3, mint::CellType::TRIANGLE);
      reader.getMesh(&readMesh);
      EXPECT_EQ(readMesh.getNumberOfCells(), 8);

//      std::cout << "readMesh.getNumberOfFaces()=" << readMesh.getNumberOfCells() << std::endl;
      axom::Array<double> xc, yc;
      getCoordinates(readMesh, xc, yc);
#if defined(AXOM_SET_BASELINES)
      // Write out results to put into source code.
      testing::writeArray(xc, "xc");
      testing::writeArray(yc, "yc");
#else
      EXPECT_TRUE(testing::compareArrays(xc, baselineXCoordinates()));
      EXPECT_TRUE(testing::compareArrays(yc, baselineYCoordinates()));
#endif
    }

//    axom::utilities::filesystem::removeFile(filename);
  }
};

} // end namespace testing

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
TEST(quest_stl_writer, uniform2d)
{
  // Make mesh.
  const double lower_bound[] = {0., 1.};
  const double upper_bound[] = {1., 2.};
  constexpr axom::IndexType NI = 3;
  constexpr axom::IndexType NJ = 3;
  mint::UniformMesh mesh(lower_bound, upper_bound, NI, NJ);

  testing::Test2D tester;
  tester.test(mesh, "uniform2d.stl", false);
  tester.test(mesh, "uniform2dB.stl", true);
}

#if 0
  const double x_expected[] = {0.0, 1.0, 0.0};
  const double y_expected[] = {0.0, 0.0, 1.0};
  const double z_expected[] = {0.0, 0.0, 0.0};

  const std::string filename = "triangle.stl";

  // STEP 0: generate a temporary STL file for testing
  generate_stl_file(filename);

  // STEP 1: create an STL reader and read-in the mesh data
  quest::STLReader reader;
  reader.setFileName(filename);
  int status = reader.read();
  EXPECT_EQ(status, 0);

  // STEP 2: reading the STL mesh data into a mint::Mesh
  mint::UnstructuredMesh<mint::SINGLE_SHAPE> mesh(3, mint::TRIANGLE);
  reader.getMesh(&mesh);

  // STEP 3: ensure the mesh is what is expected
  EXPECT_EQ(mesh.getNumberOfCells(), 1);
  EXPECT_EQ(mesh.getNumberOfNodes(), 3);

  const double* x = mesh.getCoordinateArray(mint::X_COORDINATE);
  const double* y = mesh.getCoordinateArray(mint::Y_COORDINATE);
  const double* z = mesh.getCoordinateArray(mint::Z_COORDINATE);
  EXPECT_TRUE(x != nullptr);
  EXPECT_TRUE(y != nullptr);
  EXPECT_TRUE(z != nullptr);

  axom::IndexType numNodes = mesh.getNumberOfNodes();
  for(axom::IndexType inode = 0; inode < numNodes; ++inode)
  {
    EXPECT_NEAR(x[inode], x_expected[inode], axom::numeric_limits<double>::epsilon());
    EXPECT_NEAR(y[inode], y_expected[inode], axom::numeric_limits<double>::epsilon());
    EXPECT_NEAR(z[inode], z_expected[inode], axom::numeric_limits<double>::epsilon());
  }  // END for all nodes

  // STEP 4: remove temporary STL file
  axom::utilities::filesystem::removeFile(filename);
}

//------------------------------------------------------------------------------
TEST(quest_stl_writer, read_stl_external)
{
  constexpr axom::IndexType N_NODES = 3;
  constexpr axom::IndexType N_FACES = 1;
  const double x_expected[] = {0.0, 1.0, 0.0};
  const double y_expected[] = {0.0, 0.0, 1.0};
  const double z_expected[] = {0.0, 0.0, 0.0};

  double xin[] = {-1.0, -1.0, -1.0};
  double yin[] = {-1.0, -1.0, -1.0};
  double zin[] = {-1.0, -1.0, -1.0};

  axom::IndexType conn[] = {-1, -1, -1};

  const std::string filename = "triangle.stl";

  // STEP 0: generate a temporary STL file for testing
  generate_stl_file(filename);

  // STEP 1: create an STL reader and read-in the mesh data
  quest::STLReader reader;
  reader.setFileName(filename);
  int status = reader.read();
  EXPECT_EQ(status, 0);

  // STEP 2: reading the STL mesh data into a mint::Mesh
  mint::UnstructuredMesh<mint::SINGLE_SHAPE> mesh(mint::TRIANGLE, N_FACES, conn, N_NODES, xin, yin, zin);
  EXPECT_EQ(mesh.getNumberOfCells(), N_FACES);
  EXPECT_EQ(mesh.getNumberOfNodes(), N_NODES);

  reader.getMesh(&mesh);

  const double* x = mesh.getCoordinateArray(mint::X_COORDINATE);
  const double* y = mesh.getCoordinateArray(mint::Y_COORDINATE);
  const double* z = mesh.getCoordinateArray(mint::Z_COORDINATE);
  EXPECT_TRUE(x != nullptr);
  EXPECT_TRUE(y != nullptr);
  EXPECT_TRUE(z != nullptr);

  axom::IndexType numNodes = mesh.getNumberOfNodes();
  for(axom::IndexType inode = 0; inode < numNodes; ++inode)
  {
    EXPECT_NEAR(x[inode], x_expected[inode], axom::numeric_limits<double>::epsilon());
    EXPECT_NEAR(y[inode], y_expected[inode], axom::numeric_limits<double>::epsilon());
    EXPECT_NEAR(z[inode], z_expected[inode], axom::numeric_limits<double>::epsilon());
  }  // END for all nodes

  // STEP 4: remove temporary STL file
  axom::utilities::filesystem::removeFile(filename);
}
#endif

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  return RUN_ALL_TESTS();
}
