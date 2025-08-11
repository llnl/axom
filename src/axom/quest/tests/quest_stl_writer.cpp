// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/io/STLWriter.hpp"
#include "axom/quest/io/STLReader.hpp"
#include "axom/mint/mesh/UniformMesh.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/mint/mesh/RectilinearMesh.hpp"
#include "axom/mint/mesh/CurvilinearMesh.hpp"

#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/slic.hpp"
#include "axom/fmt.hpp"

// gtest includes
#include "gtest/gtest.h"

// C/C++ includes
#include <cstdio>
#include <string>
#include <fstream>

// Uncomment the line below to write new baseline data to stdout. Then paste new
// baselines into the various baseline methods.
// #define AXOM_SET_BASELINES

// namespace aliases
namespace mint = axom::mint;
namespace quest = axom::quest;

namespace testing
{

void writeArray(const axom::Array<double> &vec, const std::string &var_name = "v")
{
  axom::fmt::print("axom::Array<double> {} = {{{}}};\n", var_name, axom::fmt::join(vec, ", "));
}

/// Convert mesh coordinates into arrays that can be easily compared.
void getCoordinates(const mint::Mesh &mesh, axom::Array<double> &xc, axom::Array<double> &yc)
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

/// Convert mesh coordinates into arrays that can be easily compared.
void getCoordinates(const mint::Mesh &mesh,
                    axom::Array<double> &xc,
                    axom::Array<double> &yc,
                    axom::Array<double> &zc)
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
      zc.push_back(mesh.getCoordinateArray(mint::Z_COORDINATE)[nodes[i]]);
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

/*!
 * \brief This testing class writes a mesh to STL, reads it back, and checks
 *        that the coordinates used by all the triangles are as-expected.
 */
struct Test2D
{
  axom::Array<double> baselineXCoordinates()
  {
    return axom::Array<double> {{0, 0.5, 0.5, 0, 0.5, 0, 0.5, 1, 1, 0.5, 1, 0.5,
                                 0, 0.5, 0.5, 0, 0.5, 0, 0.5, 1, 1, 0.5, 1, 0.5}};
  }

  axom::Array<double> baselineYCoordinates()
  {
    return axom::Array<double> {{1,   1,   1.5, 1,   1.5, 1.5, 1,   1,   1.5, 1,   1.5, 1.5,
                                 1.5, 1.5, 2,   1.5, 2,   2,   1.5, 1.5, 2,   1.5, 2,   2}};
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
      // Get the file data as a mint mesh.
      mint::UnstructuredMesh<mint::SINGLE_SHAPE> readMesh(3, mint::CellType::TRIANGLE);
      reader.getMesh(&readMesh);
      EXPECT_EQ(readMesh.getNumberOfCells(), 8);

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

    axom::utilities::filesystem::removeFile(filename);
  }
};

/*!
 * \brief This testing class writes a mesh to STL, reads it back, and checks
 *        that the coordinates used by all the triangles are as-expected.
 */
struct Test3D
{
  virtual axom::Array<double> baselineXCoordinates()
  {
    return axom::Array<double> {
      {0,   0,   0,   0,   0,   0,   0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1,   1,   1,   1,   1,   1,
       0,   0,   0,   0,   0,   0,   0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1,   1,   1,   1,   1,   1,
       0,   0.5, 0.5, 0,   0.5, 0,   0.5, 1,   1,   0.5, 1,   0.5, 0,   0.5, 0.5, 0,   0.5, 0,
       0.5, 1,   1,   0.5, 1,   0.5, 0,   0.5, 0.5, 0,   0.5, 0,   0.5, 1,   1,   0.5, 1,   0.5,
       0,   0.5, 0.5, 0,   0.5, 0,   0.5, 1,   1,   0.5, 1,   0.5, 0,   0.5, 0.5, 0,   0.5, 0,
       0.5, 1,   1,   0.5, 1,   0.5, 0,   0.5, 0.5, 0,   0.5, 0,   0.5, 1,   1,   0.5, 1,   0.5,
       0,   0.5, 0.5, 0,   0.5, 0,   0.5, 1,   1,   0.5, 1,   0.5}};
  }

  virtual axom::Array<double> baselineYCoordinates()
  {
    return axom::Array<double> {
      {1,   1,   1.5, 1,   1.5, 1.5, 1,   1,   1.5, 1,   1.5, 1.5, 1,   1,   1.5, 1,   1.5, 1.5,
       1.5, 1.5, 2,   1.5, 2,   2,   1.5, 1.5, 2,   1.5, 2,   2,   1.5, 1.5, 2,   1.5, 2,   2,
       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1.5, 1.5, 1.5, 1.5, 1.5, 1.5,
       1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,
       1,   1,   1.5, 1,   1.5, 1.5, 1,   1,   1.5, 1,   1.5, 1.5, 1.5, 1.5, 2,   1.5, 2,   2,
       1.5, 1.5, 2,   1.5, 2,   2,   1,   1,   1.5, 1,   1.5, 1.5, 1,   1,   1.5, 1,   1.5, 1.5,
       1.5, 1.5, 2,   1.5, 2,   2,   1.5, 1.5, 2,   1.5, 2,   2}};
  }

  virtual axom::Array<double> baselineZCoordinates()
  {
    return axom::Array<double> {
      {2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2,
       2, 3, 3, 2, 3, 2, 2, 2, 3, 2, 3, 3, 2, 2, 3, 2, 3, 3, 2, 2, 3, 2, 3, 3, 2, 2, 3, 2, 3, 3,
       2, 2, 3, 2, 3, 3, 2, 2, 3, 2, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
       2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}};
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
      // Get the file data as a mint mesh.
      mint::UnstructuredMesh<mint::SINGLE_SHAPE> readMesh(3, mint::CellType::TRIANGLE);
      reader.getMesh(&readMesh);
      EXPECT_EQ(readMesh.getNumberOfCells(), 40);

      axom::Array<double> xc, yc, zc;
      getCoordinates(readMesh, xc, yc, zc);
#if defined(AXOM_SET_BASELINES)
      // Write out results to put into source code.
      testing::writeArray(xc, "xc");
      testing::writeArray(yc, "yc");
      testing::writeArray(zc, "zc");
#else
      EXPECT_TRUE(testing::compareArrays(xc, baselineXCoordinates()));
      EXPECT_TRUE(testing::compareArrays(yc, baselineYCoordinates()));
      EXPECT_TRUE(testing::compareArrays(zc, baselineZCoordinates()));
#endif
    }

    axom::utilities::filesystem::removeFile(filename);
  }
};

/*!
 * \brief This testing class writes a mesh to STL, reads it back, and checks
 *        that the coordinates used by all the triangles are as-expected.
 */
struct Test3DUns : public Test3D
{
  virtual axom::Array<double> baselineXCoordinates() override
  {
    return axom::Array<double> {
      {0,   0,   0.5, 0,   0.5, 0.5, 0.5, 0.5, 0,   0.5, 0,   0,   0,   0,   0,   0,   0,   0,
       0.5, 0.5, 1,   0.5, 1,   1,   1,   1,   0.5, 1,   0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
       1,   1,   1,   1,   1,   1,   0,   0,   0.5, 0,   0.5, 0.5, 0,   0.5, 0.5, 0,   0.5, 0,
       0,   0,   0,   0,   0,   0,   0.5, 0.5, 1,   0.5, 1,   1,   0.5, 1,   1,   0.5, 1,   0.5,
       0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1,   1,   1,   1,   1,   1,   0,   0.5, 0.5, 0,   0.5, 0,
       0.5, 1,   1,   0.5, 1,   0.5, 0,   0.5, 0.5, 0,   0.5, 0,   0.5, 1,   1,   0.5, 1,   0.5,
       0,   0.5, 0.5, 0,   0.5, 0,   0.5, 1,   1,   0.5, 1,   0.5}};
  }

  virtual axom::Array<double> baselineYCoordinates() override
  {
    return axom::Array<double> {
      {1,   1.5, 1.5, 1,   1.5, 1,   1,   1,   1,   1,   1,   1,   1,   1,   1.5, 1,   1.5, 1.5,
       1,   1.5, 1.5, 1,   1.5, 1,   1,   1,   1,   1,   1,   1,   1,   1.5, 1.5, 1,   1.5, 1,
       1,   1.5, 1.5, 1,   1.5, 1,   1.5, 2,   2,   1.5, 2,   1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5,
       1.5, 1.5, 2,   1.5, 2,   2,   1.5, 2,   2,   1.5, 2,   1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5,
       1.5, 2,   2,   1.5, 2,   1.5, 1.5, 2,   2,   1.5, 2,   1.5, 2,   2,   2,   2,   2,   2,
       2,   2,   2,   2,   2,   2,   1,   1,   1.5, 1,   1.5, 1.5, 1,   1,   1.5, 1,   1.5, 1.5,
       1.5, 1.5, 2,   1.5, 2,   2,   1.5, 1.5, 2,   1.5, 2,   2}};
  }

  virtual axom::Array<double> baselineZCoordinates() override
  {
    return axom::Array<double> {
      {2, 2, 2, 2, 2, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 2, 3, 2,
       2, 2, 3, 2, 3, 3, 2, 2, 3, 2, 3, 3, 2, 2, 2, 2, 2, 2, 3, 3, 2, 3, 2, 2, 2, 3, 3, 2, 3, 2,
       2, 2, 2, 2, 2, 2, 3, 3, 2, 3, 2, 2, 2, 2, 3, 2, 3, 3, 2, 2, 3, 2, 3, 3, 3, 3, 2, 3, 2, 2,
       3, 3, 2, 3, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}};
  }
};

}  // end namespace testing

//------------------------------------------------------------------------------
// UNIT TESTS
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

TEST(quest_stl_writer, rectilinear2d)
{
  // Make mesh.
  const double x[] = {0., 0.5, 1.};
  const double y[] = {1., 1.5, 2.};
  constexpr axom::IndexType NI = 3;
  constexpr axom::IndexType NJ = 3;
  mint::RectilinearMesh mesh(NI, const_cast<double *>(x), NJ, const_cast<double *>(y));

  testing::Test2D tester;
  tester.test(mesh, "rectilinear2d.stl", false);
  tester.test(mesh, "rectilinear2dB.stl", true);
}

TEST(quest_stl_writer, curvilinear2d)
{
  // Make mesh.
  const double x[] = {0., 0.5, 1., 0., 0.5, 1., 0., 0.5, 1.};
  const double y[] = {1., 1., 1., 1.5, 1.5, 1.5, 2., 2., 2.};
  constexpr axom::IndexType NI = 3;
  constexpr axom::IndexType NJ = 3;
  mint::CurvilinearMesh mesh(NI, const_cast<double *>(x), NJ, const_cast<double *>(y));

  testing::Test2D tester;
  tester.test(mesh, "curvilinear2d.stl", false);
  tester.test(mesh, "curvilinear2dB.stl", true);
}

TEST(quest_stl_writer, unstructured2d)
{
  // Make mesh.
  const double x[] = {0., 0.5, 1., 0., 0.5, 1., 0., 0.5, 1.};
  const double y[] = {1., 1., 1., 1.5, 1.5, 1.5, 2., 2., 2.};
  const axom::IndexType conn[] = {0, 1, 4, 0, 4, 3, 1, 2, 5, 1, 5, 4,
                                  3, 4, 7, 3, 7, 6, 4, 5, 8, 4, 8, 7};
  constexpr axom::IndexType nnodes = 9;
  constexpr axom::IndexType numTriangles = 8;
  mint::UnstructuredMesh<mint::SINGLE_SHAPE> mesh(mint::CellType::TRIANGLE,
                                                  numTriangles,  // ncells
                                                  numTriangles,  // cell_capacity
                                                  const_cast<axom::IndexType *>(conn),
                                                  nnodes,  // nnodes
                                                  nnodes,  // node_capacity
                                                  const_cast<double *>(x),
                                                  const_cast<double *>(y));

  testing::Test2D tester;
  tester.test(mesh, "unstructured2d.stl", false);
  tester.test(mesh, "unstructured2dB.stl", true);
}

TEST(quest_stl_writer, uniform3d)
{
  // Make mesh.
  const double lower_bound[] = {0., 1., 2.};
  const double upper_bound[] = {1., 2., 3.};
  constexpr axom::IndexType NI = 3;
  constexpr axom::IndexType NJ = 3;
  constexpr axom::IndexType NK = 2;
  mint::UniformMesh mesh(lower_bound, upper_bound, NI, NJ, NK);

  testing::Test3D tester;
  tester.test(mesh, "uniform3d.stl", false);
  tester.test(mesh, "uniform3dB.stl", true);
}

TEST(quest_stl_writer, rectilinear3d)
{
  // Make mesh.
  const double x[] = {0., 0.5, 1.};
  const double y[] = {1., 1.5, 2.};
  const double z[] = {2., 3.};
  constexpr axom::IndexType NI = 3;
  constexpr axom::IndexType NJ = 3;
  constexpr axom::IndexType NK = 2;
  mint::RectilinearMesh mesh(NI,
                             const_cast<double *>(x),
                             NJ,
                             const_cast<double *>(y),
                             NK,
                             const_cast<double *>(z));

  testing::Test3D tester;
  tester.test(mesh, "rectilinear3d.stl", false);
  tester.test(mesh, "rectilinear3dB.stl", true);
}

TEST(quest_stl_writer, curvilinear3d)
{
  // Make mesh.
  const double x[] = {0., 0.5, 1., 0., 0.5, 1., 0., 0.5, 1., 0., 0.5, 1., 0., 0.5, 1., 0., 0.5, 1.};
  const double y[] = {1., 1., 1., 1.5, 1.5, 1.5, 2., 2., 2., 1., 1., 1., 1.5, 1.5, 1.5, 2., 2., 2.};
  const double z[] = {2., 2., 2., 2., 2., 2., 2., 2., 2., 3., 3., 3., 3., 3., 3., 3., 3., 3.};
  constexpr axom::IndexType NI = 3;
  constexpr axom::IndexType NJ = 3;
  constexpr axom::IndexType NK = 2;
  mint::CurvilinearMesh mesh(NI,
                             const_cast<double *>(x),
                             NJ,
                             const_cast<double *>(y),
                             NK,
                             const_cast<double *>(z));

  testing::Test3D tester;
  tester.test(mesh, "curvilinear3d.stl", false);
  tester.test(mesh, "curvilinear3dB.stl", true);
}

TEST(quest_stl_writer, unstructured3d)
{
  // Make mesh.
  const double x[] = {0., 0.5, 1., 0., 0.5, 1., 0., 0.5, 1., 0., 0.5, 1., 0., 0.5, 1., 0., 0.5, 1.};
  const double y[] = {1., 1., 1., 1.5, 1.5, 1.5, 2., 2., 2., 1., 1., 1., 1.5, 1.5, 1.5, 2., 2., 2.};
  const double z[] = {2., 2., 2., 2., 2., 2., 2., 2., 2., 3., 3., 3., 3., 3., 3., 3., 3., 3.};
  // clang-format off
  const axom::IndexType conn[] = {
    0,1,4,3,9,10,13,12,
    1,2,5,4,10,11,14,13,
    3,4,7,6,12,13,16,15,
    4,5,8,7,13,14,17,16
  };
  // clang-format on

  constexpr axom::IndexType nnodes = 18;
  constexpr axom::IndexType ncells = 4;
  mint::UnstructuredMesh<mint::SINGLE_SHAPE> mesh(mint::CellType::HEX,
                                                  ncells,  // ncells
                                                  ncells,  // cell_capacity
                                                  const_cast<axom::IndexType *>(conn),
                                                  nnodes,  // nnodes
                                                  nnodes,  // node_capacity
                                                  const_cast<double *>(x),
                                                  const_cast<double *>(y),
                                                  const_cast<double *>(z));
  mesh.initializeFaceConnectivity();

  testing::Test3DUns tester;
  tester.test(mesh, "unstructured3d.stl", false);
  tester.test(mesh, "unstructured3dB.stl", true);
}

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  return RUN_ALL_TESTS();
}
