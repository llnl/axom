// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/mint.hpp"
#include "axom/slic.hpp"

#include "axom/quest/LinearizeCurves.hpp"
#include "axom/quest/GWNMethods.hpp"
#include "axom/quest/io/MFEMReader.hpp"
#include "axom/quest/util/mesh_helpers.hpp"

#include "gtest/gtest.h"

#include "mfem.hpp"

#include <vector>
#include <math.h>

// For use in the 3D test cases
#ifdef AXOM_USE_OPENCASCADE
  #include "axom/quest/io/STEPReader.hpp"
#endif

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

//------------------------------------------------------------------------------
TEST(quest_winding_number_approximations, mfem_mesh_linearization)
{
  using NURBSCurve2D = axom::primal::NURBSCurve<double, 2>;
  const std::string fileName = pjoin(AXOM_DATA_DIR, "contours", "svg", "mfem_logo_simp.mesh");

  // Read the curves from the MFEM mesh
  axom::quest::MFEMReader mfem_reader;
  mfem_reader.setFileName(fileName);

  axom::Array<NURBSCurve2D> curves;
  const int ret = mfem_reader.read(curves);
  if(ret != 0)
  {
    SLIC_ERROR(axom::fmt::format("Failed to read STEP file '{}'", fileName));
  }

  // Get a linearization of the shape
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> poly_mesh(2, axom::mint::SEGMENT);
  axom::quest::LinearizeCurves lc;
  lc.getLinearMeshUniform(curves.view(), &poly_mesh, 10);

  // Get bounding box of the shape
  // Extract the curves and compute their bounding boxes along the way
  axom::primal::BoundingBox<double, 2> shape_bbox;
  for(const auto &cur : curves)
  {
    shape_bbox.addBox(cur.boundingBox());
  }

  const std::vector<int> resolution {10, 10};
  const int query_order = 1;

  // Generate three query grids and fields
  axom::Array<mfem::DataCollection> dc(0, 3);
  std::string names[] = {"direct", "polyline", "polyline_fast"};
  for(int i = 0; i < 3; ++i)
  {
    dc.emplace_back(axom::fmt::format("gwn_{}", names[i]));
    axom::quest::generate_gwn_query_mesh(dc[i],
                                         shape_bbox,
                                         std::vector<double> {},
                                         std::vector<double> {},
                                         resolution,
                                         query_order);
  }

  // Create tolerance object
  axom::primal::WindingTolerances tol;
  constexpr bool useDirectPolyline = true;

  //// Run three different kinds of GWN query ////

  // Direct
  SLIC_INFO("Testing Direct Evaluation");
  axom::quest::DirectGWN2D gwn_direct {};
  gwn_direct.preprocess(curves);
  gwn_direct.query(dc[0], tol);

  // Linearized
  SLIC_INFO("Testing Direct Evaluation of Triangulation");
  axom::quest::PolylineGWN2D<0> gwn_polyline {};
  gwn_polyline.preprocess(&poly_mesh, useDirectPolyline);
  gwn_polyline.query(dc[1], tol);

  // Linearized, fast approximation
  SLIC_INFO("Testing Fast-Approximate Evaluation of Triangulation");
  axom::quest::PolylineGWN2D<0> gwn_polyline_fast {};
  gwn_polyline_fast.preprocess(&poly_mesh, !useDirectPolyline);
  gwn_polyline_fast.query(dc[2], tol);

  // Compare the in-out values between the three fields
  const auto *query_mesh = dc[0].GetMesh();
  const auto num_query_points = query_mesh->GetNodalFESpace()->GetNDofs();

  auto &inout_direct = *dc[0].GetField("inout");
  auto &inout_polyline = *dc[1].GetField("inout");
  auto &inout_polyline_fast = *dc[2].GetField("inout");

  for(int i = 0; i < num_query_points; ++i)
  {
    EXPECT_EQ(inout_direct[i], inout_polyline[i]);
    EXPECT_EQ(inout_direct[i], inout_polyline_fast[i]);
  }
}

#ifdef AXOM_USE_OPENCASCADE
//------------------------------------------------------------------------------
TEST(quest_winding_number_approximations, step_file_triangulation)
{
  const std::string fileName = pjoin(AXOM_DATA_DIR, "quest", "step", "nut.step");

  // Read the step file
  axom::quest::STEPReader step_reader;
  step_reader.setFileName(fileName);

  constexpr bool validate = false;
  const int ret = step_reader.read(validate);
  if(ret != 0)
  {
    SLIC_ERROR(axom::fmt::format("Failed to read STEP file '{}'", fileName));
  }

  // Get patches from the shape
  const auto patches = step_reader.getPatchArray();

  // Get a triangulation of the shape
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> tri_mesh(3, axom::mint::TRIANGLE);
  step_reader.getTriangleMesh(&tri_mesh, 0.01, 0.5);

  // Get bounding box of the shape
  axom::primal::BoundingBox<double, 3> shape_bbox = step_reader.getBRepBoundingBox();
  const std::vector<int> resolution {5, 5, 5};
  const int query_order = 1;

  // Generate three query grids and fields
  axom::Array<mfem::DataCollection> dc(0, 3);
  std::string names[] = {"direct", "tri", "tri_fast"};
  for(int i = 0; i < 3; ++i)
  {
    dc.emplace_back(axom::fmt::format("gwn_{}", names[i]));
    axom::quest::generate_gwn_query_mesh(dc[i],
                                         shape_bbox,
                                         std::vector<double> {},
                                         std::vector<double> {},
                                         resolution,
                                         query_order);
  }

  // Create tolerance object
  axom::primal::WindingTolerances tol;
  constexpr bool useDirectTriangle = true;

  //// Run three different kinds of GWN query ////

  // Direct
  SLIC_INFO("Testing Direct Evaluation");
  axom::quest::DirectGWN3D gwn_direct {};
  gwn_direct.preprocess(patches);
  gwn_direct.query(dc[0], tol);

  // Triangulated
  SLIC_INFO("Testing Direct Evaluation of Polyline");
  axom::quest::TriangleGWN3D<0> gwn_tri {};
  gwn_tri.preprocess(&tri_mesh, useDirectTriangle);
  gwn_tri.query(dc[1], tol);

  // Triangulated, fast approximation
  SLIC_INFO("Testing Fast-Approximate Evaluation of Polyline");
  axom::quest::TriangleGWN3D<0> gwn_tri_fast {};
  gwn_tri_fast.preprocess(&tri_mesh, !useDirectTriangle);
  gwn_tri_fast.query(dc[2], tol);

  // Compare the in-out values between the three fields
  const auto *query_mesh = dc[0].GetMesh();
  const auto num_query_points = query_mesh->GetNodalFESpace()->GetNDofs();

  auto &inout_direct = *dc[0].GetField("inout");
  auto &inout_tri = *dc[1].GetField("inout");
  auto &inout_tri_fast = *dc[2].GetField("inout");

  for(int i = 0; i < num_query_points; ++i)
  {
    EXPECT_EQ(inout_direct[i], inout_tri[i]);
    EXPECT_EQ(inout_direct[i], inout_tri_fast[i]);
  }
}
#endif

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
