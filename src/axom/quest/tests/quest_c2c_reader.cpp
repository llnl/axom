// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core/utilities/FileUtilities.hpp"

#ifndef AXOM_USE_C2C
  #error These tests should only be included when Axom is configured with C2C
#endif

#include "axom/quest/io/C2CReader.hpp"
#include "axom/quest/LinearizeCurves.hpp"
#include "axom/slic.hpp"
#include "axom/mint.hpp"
#include "axom/primal.hpp"
#include "axom/fmt.hpp"

// gtest includes
#include "gtest/gtest.h"

// C/C++ includes
#include <cstdio>
#include <string>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <math.h>

// namespace aliases
namespace mint = axom::mint;
namespace primal = axom::primal;
namespace quest = axom::quest;
namespace utilities = axom::utilities;

namespace
{
const std::string C2C_LINE_FILENAME = "test_line.contour";
const std::string C2C_CIRCLE_FILENAME = "test_circle.contour";
const std::string C2C_SQUARE_FILENAME = "test_square.contour";
const std::string C2C_SPLINE_FILENAME = "test_spline.contour";
const std::string C2C_NESTED_ASSEMBLY_FILENAME = "test_nested.assembly";
const std::string C2C_TOP_ASSEMBLY_FILENAME = "test_top.assembly";
}  // end anonymous namespace

/// Writes out a c2c file for a circle
void writeSimpleCircle(const std::string& filename)
{
  std::ofstream c2cFile(filename, std::ios::out);
  c2cFile << "piece = circle(origin=(0cm, 0cm), radius=1cm, start=0deg, end=360deg)" << std::endl;
}

/// Writes out a c2c file for a line
void writeSimpleLine(const std::string& filename)
{
  std::ofstream c2cFile(filename, std::ios::out);
  c2cFile << "piece = line(start=(0cm, 0cm), end=(1cm, 1cm))" << std::endl;
}

/// Writes out a c2c file for a unit square
void writeSquare(const std::string& filename)
{
  std::ofstream c2cFile(filename, std::ios::out);
  c2cFile << "point = start" << std::endl;
  c2cFile << "piece = line(start=(0cm, 0cm), end=(0cm, 1cm))" << std::endl;
  c2cFile << "piece = line()" << std::endl;
  c2cFile << "piece = line(start=(1cm, 1cm), end=(1cm, 0cm))" << std::endl;
  c2cFile << "piece = line(end=start)" << std::endl;
}

/// Writes out a c2c file for a rectangle cut by a cosine wave
void writeSpline(const std::string& filename)
{
  std::vector<std::string> pts;
  // generate a cubic spline approximation to a cosine wave
  // with domain [0, 2*PI] and range defined by a y-offset, amplitude and frequency
  const double offset = 1.;
  const double amplitude = .5;
  const double freq = 3;
  const int NPTS = 2 * freq;
  for(int i = 0; i <= NPTS; ++i)
  {
    double y = 2 * M_PI * static_cast<double>(i) / NPTS;
    double x = offset + amplitude * cos(freq * y);
    pts.emplace_back(axom::fmt::format("{} {}", x, y));
  }

  std::ofstream c2cFile(filename, std::ios::out);
  // output sine wave spline
  c2cFile << "point = spline_start" << std::endl;
  c2cFile << axom::fmt::format(
               "piece = rz(units=cm, spline=cubic, beginTan=-180deg, "
               "endTan=-180deg, rz={})",
               axom::fmt::join(pts.rbegin(), pts.rend(), "\n\t\t"))
          << std::endl;
  c2cFile << "point = spline_end" << std::endl;

  // add straight edges within first quadrant
  c2cFile << "piece = line(end=(0cm,0cm))" << std::endl;
  c2cFile << axom::fmt::format("piece = line(end=(0cm,{}cm))", 2 * M_PI) << std::endl;
  c2cFile << "piece = line(end=spline_start)" << std::endl;
}

/// Writes out a c2c assembly file with the given entries
void writeAssembly(const std::string& filename, const std::vector<std::string>& entries)
{
  std::ofstream c2cFile(filename, std::ios::out);
  for(const auto& entry : entries)
  {
    c2cFile << entry << std::endl;
  }
}

TEST(quest_c2c_reader, unsupported_length_units)
{
  quest::C2CReader reader;
  EXPECT_THROW(reader.setLengthUnit(utilities::LengthUnit::am), std::invalid_argument);
  EXPECT_THROW(reader.setLengthUnit(utilities::LengthUnit::fm), std::invalid_argument);
  EXPECT_THROW(reader.setLengthUnit(utilities::LengthUnit::pm), std::invalid_argument);
  EXPECT_THROW(reader.setLengthUnit(utilities::LengthUnit::dm), std::invalid_argument);
  EXPECT_THROW(reader.setLengthUnit(utilities::LengthUnit::dam), std::invalid_argument);
  EXPECT_THROW(reader.setLengthUnit(utilities::LengthUnit::hm), std::invalid_argument);
  EXPECT_THROW(reader.setLengthUnit(utilities::LengthUnit::nm), std::invalid_argument);
  EXPECT_THROW(reader.setLengthUnit(utilities::LengthUnit::angstrom), std::invalid_argument);
  EXPECT_THROW(reader.setLengthUnit(utilities::LengthUnit::unspecified), std::invalid_argument);
}

TEST(quest_c2c_reader, basic_read)
{
  const std::string fileName = C2C_CIRCLE_FILENAME;
  writeSimpleCircle(fileName);

  quest::C2CReader reader;
  reader.setFileName(fileName);

  ASSERT_EQ(0, reader.read());
  reader.log();
}

TEST(quest_c2c_reader, interpolate_circle)
{
  const std::string fileName = C2C_CIRCLE_FILENAME;
  writeSimpleCircle(fileName);

  quest::C2CReader reader;
  reader.setFileName(fileName);

  ASSERT_EQ(0, reader.read());
  reader.log();

  constexpr int DIM = 2;
  using MeshType = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;
  MeshType* mesh = new MeshType(DIM, mint::SEGMENT);

  const int segmentsPerKnotSpan = 25;
  axom::quest::LinearizeCurves lin;
  lin.getLinearMeshUniform(reader.getCurvesView(), mesh, segmentsPerKnotSpan);

  // The circle is defined by a single NURBS curve with four spans
  SLIC_INFO(axom::fmt::format("Mesh has {} nodes and {} cells",
                              mesh->getNumberOfNodes(),
                              mesh->getNumberOfCells()));
  const int numSpans = 4;
  const int expVerts = numSpans * (segmentsPerKnotSpan + 1);
  EXPECT_EQ(expVerts, mesh->getNumberOfNodes());
  const int expSegments = numSpans * segmentsPerKnotSpan;
  EXPECT_EQ(expSegments, mesh->getNumberOfCells());

  // This is a unit circle; check that its vertices have unit magnitude
  {
    double* x = mesh->getCoordinateArray(mint::X_COORDINATE);
    double* y = mesh->getCoordinateArray(mint::Y_COORDINATE);
    const int numPts = mesh->getNumberOfNodes();
    for(int i = 0; i < numPts; ++i)
    {
      double mag = primal::Vector<double, 2> {x[i], y[i]}.norm();
      EXPECT_DOUBLE_EQ(1., mag);
    }
  }

  mint::write_vtk(mesh, "test_circle.vtk");

  delete mesh;
}

TEST(quest_c2c_reader, read_with_axom_length_unit)
{
  const std::string fileName = C2C_CIRCLE_FILENAME;
  writeSimpleCircle(fileName);

  quest::C2CReader reader;
  reader.setFileName(fileName);
  reader.setLengthUnit(utilities::LengthUnit::mm);

  EXPECT_EQ(0, reader.read());

  constexpr int DIM = 2;
  using MeshType = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;
  MeshType* mesh = new MeshType(DIM, mint::SEGMENT);

  const int segmentsPerKnotSpan = 25;
  axom::quest::LinearizeCurves lin;
  lin.getLinearMeshUniform(reader.getCurvesView(), mesh, segmentsPerKnotSpan);

  double* x = mesh->getCoordinateArray(mint::X_COORDINATE);
  double* y = mesh->getCoordinateArray(mint::Y_COORDINATE);
  const int numPts = mesh->getNumberOfNodes();
  for(int i = 0; i < numPts; ++i)
  {
    double mag = primal::Vector<double, 2> {x[i], y[i]}.norm();
    EXPECT_DOUBLE_EQ(10., mag);
  }

  delete mesh;
}

TEST(quest_c2c_reader, interpolate_square)
{
  const std::string fileName = C2C_SQUARE_FILENAME;
  writeSquare(fileName);

  quest::C2CReader reader;
  reader.setFileName(fileName);

  ASSERT_EQ(0, reader.read());
  reader.log();

  constexpr int DIM = 2;
  using MeshType = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;
  MeshType* mesh = new MeshType(DIM, mint::SEGMENT);

  const int segmentsPerKnotSpan = 10;
  axom::quest::LinearizeCurves lin;
  lin.getLinearMeshUniform(reader.getCurvesView(), mesh, segmentsPerKnotSpan);

  SLIC_INFO(axom::fmt::format("Mesh has {} nodes and {} cells",
                              mesh->getNumberOfNodes(),
                              mesh->getNumberOfCells()));

  const int numPieces = 4;
  const int spansPerPiece = 1;
  const int totalSpans = numPieces * spansPerPiece;
  const int expVerts = totalSpans * (segmentsPerKnotSpan + 1);
  EXPECT_EQ(expVerts, mesh->getNumberOfNodes());

  const int expSegs = totalSpans * segmentsPerKnotSpan;
  EXPECT_EQ(expSegs, mesh->getNumberOfCells());

  mint::write_vtk(mesh, "test_square.vtk");

  delete mesh;
}

TEST(quest_c2c_reader, interpolate_spline)
{
  const std::string fileName = C2C_SPLINE_FILENAME;
  writeSpline(fileName);

  quest::C2CReader reader;
  reader.setFileName(fileName);

  ASSERT_EQ(0, reader.read());
  reader.log();

  constexpr int DIM = 2;
  using MeshType = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;
  MeshType* mesh = new MeshType(DIM, mint::SEGMENT);

  const int segmentsPerKnotSpan = 20;
  axom::quest::LinearizeCurves lin;
  lin.getLinearMeshUniform(reader.getCurvesView(), mesh, segmentsPerKnotSpan);

  SLIC_INFO(axom::fmt::format("Mesh has {} nodes and {} cells",
                              mesh->getNumberOfNodes(),
                              mesh->getNumberOfCells()));

  const int numSpans = 6 + 1 + 1 + 1;
  const int expVerts = numSpans * (segmentsPerKnotSpan + 1);
  EXPECT_EQ(expVerts, mesh->getNumberOfNodes());

  const int expSegs = numSpans * segmentsPerKnotSpan;
  EXPECT_EQ(expSegs, mesh->getNumberOfCells());

  mint::write_vtk(mesh, "test_spline.vtk");

  delete mesh;
}

TEST(quest_c2c_reader, duplicate_point_linear_fails_gracefully)
{
#ifdef AXOM_DATA_DIR
  // This file contains an invalid contour -- reading the files should return non-zero
  const auto fileName =
    axom::utilities::filesystem::joinPath(AXOM_DATA_DIR, "contours/duplicate_point_linear.contour");

  quest::C2CReader reader;
  reader.setFileName(fileName);

  EXPECT_NE(0, reader.read());
  EXPECT_EQ(0, reader.getCurvesView().size());
#else
  GTEST_SKIP() << "AXOM_DATA_DIR not defined";
#endif
}

TEST(quest_c2c_reader, nested_assemblies_are_flattened)
{
  const std::string nestedLineFile = "nested_line.contour";
  const std::string nestedCircleFile = "nested_circle.contour";
  const std::string nestedSquareFile = "nested_square.contour";

  writeSimpleLine(nestedLineFile);
  writeSimpleCircle(nestedCircleFile);
  writeSquare(nestedSquareFile);

  writeAssembly(C2C_NESTED_ASSEMBLY_FILENAME,
                {"pieces = contour(path='nested_circle.contour')",
                 "pieces = contour(path='nested_square.contour')"});
  writeAssembly(C2C_TOP_ASSEMBLY_FILENAME,
                {"pieces = contour(path='nested_line.contour')",
                 "pieces = assembly(path='test_nested.assembly')"});

  quest::C2CReader reader;
  reader.setFileName(C2C_TOP_ASSEMBLY_FILENAME);

  EXPECT_EQ(0, reader.read());
  EXPECT_EQ(6, reader.getCurvesView().size());
}

TEST(quest_c2c_reader, missing_contour_file_fails_gracefully)
{
  quest::C2CReader reader;
  reader.setFileName("missing.contour");

  EXPECT_NE(0, reader.read());
  EXPECT_EQ(0, reader.getCurvesView().size());
}

TEST(quest_c2c_reader, missing_assembly_file_fails_gracefully)
{
  quest::C2CReader reader;
  reader.setFileName("missing.assembly");

  EXPECT_NE(0, reader.read());
  EXPECT_EQ(0, reader.getCurvesView().size());
}

TEST(quest_c2c_reader, missing_nested_assembly_member_fails_transactionally)
{
  writeSimpleCircle(C2C_CIRCLE_FILENAME);
  writeAssembly(C2C_TOP_ASSEMBLY_FILENAME,
                {"pieces = contour(path='test_circle.contour')",
                 "pieces = contour(path='missing_nested.contour')"});

  quest::C2CReader reader;
  reader.setFileName(C2C_TOP_ASSEMBLY_FILENAME);

  EXPECT_NE(0, reader.read());
  EXPECT_EQ(0, reader.getCurvesView().size());
}

TEST(quest_c2c_reader, heroic_roses_black_assembly_reads)
{
#ifdef AXOM_DATA_DIR
  const auto fileName =
    axom::utilities::filesystem::joinPath(AXOM_DATA_DIR, "contours/heroic_roses/c2c/black.assembly");

  quest::C2CReader reader;
  reader.setFileName(fileName);

  EXPECT_EQ(0, reader.read());
  EXPECT_GT(reader.getCurvesView().size(), 30);
#else
  GTEST_SKIP() << "AXOM_DATA_DIR not defined";
#endif
}

//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  return RUN_ALL_TESTS();
}
