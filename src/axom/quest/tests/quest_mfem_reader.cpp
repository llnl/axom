// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#ifndef AXOM_USE_MFEM
  #error These tests should only be included when Axom is configured with MFEM
#endif

#include "axom/quest/io/MFEMReader.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"

// gtest includes
#include "gtest/gtest.h"

// namespace aliases
namespace primal = axom::primal;
namespace quest = axom::quest;

//------------------------------------------------------------------------------
std::string pjoin(const std::string &str) { return str; }

std::string pjoin(const char *str) { return std::string(str); }

template <typename... Args>
std::string pjoin(const std::string &str, Args... args)
{
  return axom::utilities::filesystem::joinPath(str, pjoin(args...));
}

template <typename... Args>
std::string pjoin(const char *str, Args... args)
{
  return axom::utilities::filesystem::joinPath(std::string(str), pjoin(args...));
}
//------------------------------------------------------------------------------

TEST(quest_mfem_reader, read_nurbs_curves)
{
  const std::string fileName =
    pjoin(AXOM_DATA_DIR, "contours", "heroic_roses", "mfem", "blue0.mesh");

  quest::MFEMReader reader;
  reader.setFileName(fileName);

  // Read as 9 NURBS curves
  axom::Array<primal::NURBSCurve<double, 2>> curves;
  EXPECT_EQ(reader.read(curves), 0);
  EXPECT_EQ(curves.size(), 9);
}

TEST(quest_mfem_reader, read_curved_polygon)
{
  const std::string fileName =
    pjoin(AXOM_DATA_DIR, "contours", "heroic_roses", "mfem", "blue0.mesh");

  quest::MFEMReader reader;
  reader.setFileName(fileName);

  // Read as 1 CurvedPolygon with 9 edges
  axom::Array<primal::CurvedPolygon<double, 2>> polys;
  EXPECT_EQ(reader.read(polys), 0);
  EXPECT_EQ(polys.size(), 1);
  EXPECT_EQ(polys[0].numEdges(), 9);

  // Read as CurvedPolygon
  polys.clear();
  const std::string fileNameB =
    pjoin(AXOM_DATA_DIR, "contours", "heroic_roses", "mfem_cp", "black.mesh");
  reader.setFileName(fileNameB);
  EXPECT_EQ(reader.read(polys), 0);
  EXPECT_EQ(polys.size(), 73);
  // Pick some curved polygons and check lengths.
  EXPECT_EQ(polys[0].numEdges(), 11);
  EXPECT_EQ(polys[20].numEdges(), 20);
  EXPECT_EQ(polys[40].numEdges(), 4);
  EXPECT_EQ(polys[72].numEdges(), 2);
}

//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  return RUN_ALL_TESTS();
}
