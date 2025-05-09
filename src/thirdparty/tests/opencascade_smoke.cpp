// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"
#include <iostream>

#include "opencascade/Standard_Version.hxx"

#include "opencascade/Bnd_Box.hxx"
#include "opencascade/BRepBndLib.hxx"
#include "opencascade/BRepBuilderAPI_MakeFace.hxx"
#include "opencascade/BRepBuilderAPI_MakeWire.hxx"
#include "opencascade/BRepBuilderAPI_MakeEdge.hxx"
#include "opencascade/BRepBuilderAPI_MakeShell.hxx"
#include "opencascade/BRepPrimAPI_MakeBox.hxx"
#include "opencascade/gp_Pnt.hxx"
#include "opencascade/IFSelect_ReturnStatus.hxx"
#include "opencascade/TopAbs.hxx"
#include "opencascade/TopoDS_Shape.hxx"
#include "opencascade/TopoDS_Face.hxx"
#include "opencascade/TopExp_Explorer.hxx"

#include "opencascade/STEPControl_Writer.hxx"
#include "opencascade/STEPControl_Reader.hxx"

//------------------------------------------------------------------------------
// Some basic smoke tests for Open Cascade
// including checks that we can read and write STEP files
//------------------------------------------------------------------------------
TEST(opencascade_smoke, check_version)
{
  std::cout << "Using Open Cascade version: "  //
            << OCC_VERSION_MAJOR << "."        //
            << OCC_VERSION_MINOR << "."        //
            << OCC_VERSION_MAINTENANCE << std::endl;

  EXPECT_TRUE(OCC_VERSION_MAJOR >= 0);
  EXPECT_TRUE(OCC_VERSION_MINOR >= 0);
  EXPECT_TRUE(OCC_VERSION_MAINTENANCE >= 0);
}

TEST(opencascade_smoke, make_box)
{
  // Basic test of Open Cascade API -- create a box and check its dimensions
  constexpr int box_length = 200;
  constexpr int box_width = 100;
  constexpr int box_depth = 50;

  TopoDS_Shape box = BRepPrimAPI_MakeBox(box_length, box_width, box_depth).Shape();
  EXPECT_FALSE(box.IsNull());

  Bnd_Box boundingBox;
  BRepBndLib::Add(box, boundingBox);

  Standard_Real xMin, yMin, zMin, xMax, yMax, zMax;
  boundingBox.Get(xMin, yMin, zMin, xMax, yMax, zMax);
  EXPECT_EQ(box_length, static_cast<int>(xMax - xMin));
  EXPECT_EQ(box_width, static_cast<int>(yMax - yMin));
  EXPECT_EQ(box_depth, static_cast<int>(zMax - zMin));
}

TEST(opencascade_smoke, generate_step_file)
{
  TopoDS_Shell shell;

  // Build up a tetrahedron into shell
  {
    // Create vertices
    gp_Pnt p1(0, 0, 0);
    gp_Pnt p2(1, 0, 0);
    gp_Pnt p3(0.5, std::sqrt(3) / 2, 0);
    gp_Pnt p4(0.5, std::sqrt(3) / 6, std::sqrt(2.0 / 3.0));

    // Create edges
    TopoDS_Edge e1 = BRepBuilderAPI_MakeEdge(p1, p2);
    TopoDS_Edge e2 = BRepBuilderAPI_MakeEdge(p2, p3);
    TopoDS_Edge e3 = BRepBuilderAPI_MakeEdge(p3, p1);
    TopoDS_Edge e4 = BRepBuilderAPI_MakeEdge(p1, p4);
    TopoDS_Edge e5 = BRepBuilderAPI_MakeEdge(p2, p4);
    TopoDS_Edge e6 = BRepBuilderAPI_MakeEdge(p3, p4);

    // Create faces
    TopoDS_Wire wire1 = BRepBuilderAPI_MakeWire(e1, e2, e3);
    TopoDS_Wire wire2 = BRepBuilderAPI_MakeWire(e1, e5, e4);
    TopoDS_Wire wire3 = BRepBuilderAPI_MakeWire(e2, e5, e6);
    TopoDS_Wire wire4 = BRepBuilderAPI_MakeWire(e3, e4, e6);

    TopoDS_Face face1 = BRepBuilderAPI_MakeFace(wire1);
    TopoDS_Face face2 = BRepBuilderAPI_MakeFace(wire2);
    TopoDS_Face face3 = BRepBuilderAPI_MakeFace(wire3);
    TopoDS_Face face4 = BRepBuilderAPI_MakeFace(wire4);

    // Add faces to shell
    BRep_Builder builder;
    builder.MakeShell(shell);

    builder.Add(shell, face1);
    builder.Add(shell, face2);
    builder.Add(shell, face3);
    builder.Add(shell, face4);
  }

  // Check the number of faces and export to STEP file
  {
    int faceCount = 0;
    for(TopExp_Explorer expl(shell, TopAbs_FACE); expl.More(); expl.Next())
    {
      faceCount++;
    }
    EXPECT_EQ(4, faceCount);

    STEPControl_Writer writer;
    writer.Transfer(shell, STEPControl_AsIs);
    writer.Write("tetrahedron.stp");
  }

  // Import the tetrahedron as a STEP file and check the number of faces
  {
    STEPControl_Reader stepReader;
    IFSelect_ReturnStatus status = stepReader.ReadFile("tetrahedron.stp");
    ASSERT_EQ(IFSelect_RetDone, status);

    stepReader.TransferRoots();
    TopoDS_Shape shape = stepReader.OneShape();
    EXPECT_FALSE(shape.IsNull());

    int faceCount = 0;
    for(TopExp_Explorer expl(shape, TopAbs_FACE); expl.More(); expl.Next())
    {
      faceCount++;
    }
    EXPECT_EQ(4, faceCount);
  }
}
