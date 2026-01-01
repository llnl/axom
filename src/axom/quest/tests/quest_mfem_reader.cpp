// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#if !defined(AXOM_USE_MFEM) || !defined(AXOM_USE_SIDRE)
  #error These tests should only be included when Axom is configured with MFEM and SIDRE
#endif

#include "axom/quest/io/MFEMReader.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"

#include "mfem.hpp"

// gtest includes
#include "gtest/gtest.h"

#include <fstream>
#include <string>

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
  axom::Array<int> attributes;
  EXPECT_EQ(reader.read(curves, attributes), quest::MFEMReader::READ_SUCCESS);
  EXPECT_EQ(curves.size(), 9);
  ASSERT_EQ(attributes.size(), 9);
  for(int i = 0; i < attributes.size(); ++i)
  {
    EXPECT_EQ(attributes[i], 1);
  }
}

TEST(quest_mfem_reader, read_curved_polygon)
{
  const std::string fileName =
    pjoin(AXOM_DATA_DIR, "contours", "heroic_roses", "mfem", "blue0.mesh");

  quest::MFEMReader reader;
  reader.setFileName(fileName);

  // Read as 1 CurvedPolygon with 9 edges
  axom::Array<primal::CurvedPolygon<axom::primal::NURBSCurve<double, 2>>> polys;
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

TEST(quest_mfem_reader, preserves_rational_weights)
{
  const std::string fileName =
    pjoin(AXOM_DATA_DIR, "contours", "heroic_roses", "mfem", "brightgreen_over.mesh");

  quest::MFEMReader reader;
  reader.setFileName(fileName);

  // test read on NURBSCurve array and check that at least one curve is rational
  {
    axom::Array<primal::NURBSCurve<double, 2>> curves;
    EXPECT_EQ(reader.read(curves), quest::MFEMReader::READ_SUCCESS);
    ASSERT_GT(curves.size(), 0);

    bool any_rational = false;
    for(const auto &curve : curves)
    {
      if(curve.isRational())
      {
        any_rational = true;
        break;
      }
    }
    EXPECT_TRUE(any_rational);
  }

  // do the same with the extracted Curved polygon array
  {
    axom::Array<primal::CurvedPolygon<axom::primal::NURBSCurve<double, 2>>> polys;
    EXPECT_EQ(reader.read(polys), quest::MFEMReader::READ_SUCCESS);
    ASSERT_GT(polys.size(), 0);

    bool any_rational = false;
    for(const auto &poly : polys)
    {
      for(const auto &cur : poly.getEdges())
      {
        if(cur.isRational())
        {
          any_rational = true;
          break;
        }
      }
    }
    EXPECT_TRUE(any_rational);
  }
}

TEST(quest_mfem_reader, read_curved_polygon_noncontiguous_attributes)
{
  axom::utilities::filesystem::TempFile tmp_mesh("noncontiguous_attributes", ".mesh");

  // This test uses non-contiguous attributes.
  // For testing, we're setting the y-coordinate to be the same as the attribute
  constexpr int attr10 {10};
  constexpr int attr20 {20};

  {
    mfem::Mesh mesh(/*Dim*/ 1, /*NVert*/ 4, /*NElem*/ 2, /*NBdrElem*/ 0, /*spaceDim*/ 2);
    const double v0[] = {0., static_cast<double>(attr20)};
    const double v1[] = {1., static_cast<double>(attr20)};
    const double v2[] = {0., static_cast<double>(attr10)};
    const double v3[] = {1., static_cast<double>(attr10)};

    mesh.AddVertex(v0);
    mesh.AddVertex(v1);
    mesh.AddVertex(v2);
    mesh.AddVertex(v3);

    mesh.AddSegment(0, 1, attr20);
    mesh.AddSegment(2, 3, attr10);

    mesh.FinalizeTopology(/*generate_bdr*/ true);
    mesh.Finalize(/*refine*/ false, /*fix_orientation*/ true);
    mesh.EnsureNodes();

    std::ofstream ofs(tmp_mesh.getPath());
    ASSERT_TRUE(ofs.good());
    mesh.Print(ofs);
  }

  quest::MFEMReader reader;
  reader.setFileName(tmp_mesh.getPath());

  // check that we can successfully read in the attributes to CurvedPolygon array
  {
    axom::Array<primal::CurvedPolygon<axom::primal::NURBSCurve<double, 2>>> polys;
    axom::Array<int> attributes;
    EXPECT_EQ(reader.read(polys, attributes), quest::MFEMReader::READ_SUCCESS);

    ASSERT_EQ(polys.size(), 2);
    ASSERT_EQ(attributes.size(), 2);

    EXPECT_EQ(polys[0].numEdges(), 1);
    EXPECT_EQ(polys[1].numEdges(), 1);

    // note: the curves are added to a map, so the order is not guaranteed
    // let's check that the attributes and geometry match expectations
    // the y-coordinates of the edges start and end vertex should equal the attribute
    for(int i : {0, 1})
    {
      const auto &edge = polys[i][0];
      if(attributes[i] == attr10)
      {
        EXPECT_EQ(edge[0][1], attr10);
        EXPECT_EQ(edge[1][1], attr10);
      }
      else if(attributes[i] == attr20)
      {
        EXPECT_EQ(edge[0][1], attr20);
        EXPECT_EQ(edge[1][1], attr20);
      }
      else
      {
        FAIL() << "Got unexpected attribute for polygon " << i << ": " << attributes[i] << "\n";
      }
    }

    polys.clear();
    EXPECT_EQ(reader.read(polys), quest::MFEMReader::READ_SUCCESS);
    EXPECT_EQ(polys.size(), 2);
  }

  // check that we can successfully read in the attributes to NURBSCurve array
  {
    axom::Array<primal::NURBSCurve<double, 2>> curves;
    axom::Array<int> curve_attributes;
    EXPECT_EQ(reader.read(curves, curve_attributes), quest::MFEMReader::READ_SUCCESS);
    ASSERT_EQ(curves.size(), 2);
    ASSERT_EQ(curve_attributes.size(), 2);

    for(int i : {0, 1})
    {
      const auto &curve = curves[i];
      if(curve_attributes[i] == attr10)
      {
        EXPECT_EQ(curve[0][1], attr10);
        EXPECT_EQ(curve[1][1], attr10);
      }
      else if(curve_attributes[i] == attr20)
      {
        EXPECT_EQ(curve[0][1], attr20);
        EXPECT_EQ(curve[1][1], attr20);
      }
      else
      {
        FAIL() << "Got unexpected attribute for curve " << i << ": " << curve_attributes[i] << "\n";
      }
    }
  }
}

//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  return RUN_ALL_TESTS();
}
