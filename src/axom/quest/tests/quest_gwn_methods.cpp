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
#include "axom/quest/FastApproximateGWN.hpp"
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
TEST(quest_gwn_methods, gwn_moment_data_segment)
{
  using Point2D = axom::primal::Point<double, 2>;
  using Moments = axom::quest::GWNMomentData<double, 2, 2>;

  constexpr double EPS = 1e-14;

  // Combine two segments into a cluster and verify raw moments, centroid, and
  // expansion coefficients. Values are chosen to make exact expectations.
  Moments s1(Point2D {0.0, 0.0}, Point2D {1.0, 0.0});
  Moments s2(Point2D {0.0, 1.0}, Point2D {0.0, 2.0});
  Moments cluster = s1 + s2;

  EXPECT_NEAR(cluster.a, 2.0, EPS);
  EXPECT_NEAR(cluster.ap[0], 0.5, EPS);
  EXPECT_NEAR(cluster.ap[1], 1.5, EPS);

  const auto center = cluster.getCenter();
  EXPECT_NEAR(center[0], 0.25, EPS);
  EXPECT_NEAR(center[1], 0.75, EPS);

  // Raw moments (rm)
  EXPECT_NEAR(cluster.rm[0], -1.0, EPS);  // dy
  EXPECT_NEAR(cluster.rm[1], 1.0, EPS);   // dx
  EXPECT_NEAR(cluster.rm[2], 0.0, EPS);   // 0.5*dy*(x0+x1)
  EXPECT_NEAR(cluster.rm[3], 0.5, EPS);   // 0.5*(x1^2-x0^2)
  EXPECT_NEAR(cluster.rm[4], -1.5, EPS);  // 0.5*(y0^2-y1^2)
  EXPECT_NEAR(cluster.rm[5], 0.0, EPS);   // 0.5*dx*(y0+y1)
  EXPECT_NEAR(cluster.rm[6], 0.0, EPS);
  EXPECT_NEAR(cluster.rm[7], 0.0, EPS);
  EXPECT_NEAR(cluster.rm[8], 0.0, EPS);
  EXPECT_NEAR(cluster.rm[9], -7.0 / 3.0, EPS);
  EXPECT_NEAR(cluster.rm[10], 1.0 / 3.0, EPS);
  EXPECT_NEAR(cluster.rm[11], 0.0, EPS);
  EXPECT_NEAR(cluster.rm[12], 0.0, EPS);
  EXPECT_NEAR(cluster.rm[13], 0.0, EPS);

  // Expansion coefficients (ec)
  EXPECT_NEAR(cluster.ec[0], -1.0, EPS);
  EXPECT_NEAR(cluster.ec[1], 1.0, EPS);
  EXPECT_NEAR(cluster.ec[2], 0.25, EPS);
  EXPECT_NEAR(cluster.ec[3], 0.25, EPS);
  EXPECT_NEAR(cluster.ec[4], -0.75, EPS);
  EXPECT_NEAR(cluster.ec[5], -0.75, EPS);
  EXPECT_NEAR(cluster.ec[6], -1.0 / 32.0, EPS);
  EXPECT_NEAR(cluster.ec[7], -3.0 / 32.0, EPS);
  EXPECT_NEAR(cluster.ec[8], 3.0 / 32.0, EPS);
  EXPECT_NEAR(cluster.ec[9], -121.0 / 96.0, EPS);
  EXPECT_NEAR(cluster.ec[10], 25.0 / 96.0, EPS);
  EXPECT_NEAR(cluster.ec[11], -3.0 / 32.0, EPS);
  EXPECT_NEAR(cluster.ec[12], 27.0 / 32.0, EPS);
  EXPECT_NEAR(cluster.ec[13], 9.0 / 32.0, EPS);

  // Verify that cluster raw moments are sum of segment moments
  for(int i = 0; i < cluster.rm.size(); ++i)
  {
    EXPECT_NEAR(s1.rm[i] + s2.rm[i], cluster.rm[i], EPS);
  }
}

//------------------------------------------------------------------------------
TEST(quest_gwn_methods, gwn_moment_data_triangle)
{
  using Point3D = axom::primal::Point<double, 3>;
  using Tri3D = axom::primal::Triangle<double, 3>;
  using Moments = axom::quest::GWNMomentData<double, 3, 2>;

  constexpr double EPS = 1e-14;

  // Right triangle in the XY plane.
  const Tri3D tri(Point3D {0.0, 0.0, 0.0}, Point3D {1.0, 0.0, 0.0}, Point3D {0.0, 1.0, 0.0});
  const Moments m(tri);

  // Area and centroid accumulation.
  EXPECT_NEAR(m.a, 0.5, EPS);
  EXPECT_NEAR(m.ap[0], 1.0 / 6.0, EPS);
  EXPECT_NEAR(m.ap[1], 1.0 / 6.0, EPS);
  EXPECT_NEAR(m.ap[2], 0.0, EPS);

  const auto center = m.getCenter();
  EXPECT_NEAR(center[0], 1.0 / 3.0, EPS);
  EXPECT_NEAR(center[1], 1.0 / 3.0, EPS);
  EXPECT_NEAR(center[2], 0.0, EPS);

  // Selected raw moments (rm): normal-only and a few higher-order entries.
  EXPECT_NEAR(m.rm[0], 0.0, EPS);
  EXPECT_NEAR(m.rm[1], 0.0, EPS);
  EXPECT_NEAR(m.rm[2], 0.5, EPS);
  EXPECT_NEAR(m.rm[5], 1.0 / 6.0, EPS);
  EXPECT_NEAR(m.rm[8], 1.0 / 6.0, EPS);
  EXPECT_NEAR(m.rm[14], 1.0 / 12.0, EPS);
  EXPECT_NEAR(m.rm[17], 1.0 / 24.0, EPS);
  EXPECT_NEAR(m.rm[23], 1.0 / 24.0, EPS);
  EXPECT_NEAR(m.rm[26], 1.0 / 12.0, EPS);

  // Selected expansion coefficients (ec): first-order terms cancel at centroid.
  EXPECT_NEAR(m.ec[0], 0.0, EPS);
  EXPECT_NEAR(m.ec[1], 0.0, EPS);
  EXPECT_NEAR(m.ec[2], 0.5, EPS);
  EXPECT_NEAR(m.ec[5], 0.0, EPS);
  EXPECT_NEAR(m.ec[8], 0.0, EPS);
  EXPECT_NEAR(m.ec[14], 1.0 / 72.0, EPS);
  EXPECT_NEAR(m.ec[17], -1.0 / 144.0, EPS);
  EXPECT_NEAR(m.ec[23], -1.0 / 144.0, EPS);
  EXPECT_NEAR(m.ec[26], 1.0 / 72.0, EPS);
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
void check_mfem_mesh_linearization()
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

  // Generate query grids and fields
  constexpr int num_queries = 6;
  axom::Array<mfem::DataCollection> dc(0, num_queries);
  for(int i = 0; i < num_queries; ++i)
  {
    dc.emplace_back(axom::fmt::format("gwn_method_{}", i));
    axom::quest::generate_gwn_query_mesh(dc[i],
                                         shape_bbox,
                                         std::vector<double> {},
                                         std::vector<double> {},
                                         resolution,
                                         query_order);
  }

  // Create tolerance object
  axom::primal::WindingTolerances tol;
  constexpr bool useDirectEvaluation = true;
  constexpr bool useMemoization = true;

  //// Run six different kinds of GWN query ////
  // We expect all three fields to return the same values in this case because
  //  of the specific arrangement of query points and linearization.
  // In general, discretizing the shape can result in different GWN values
  //  for query points near to individual curves.

  SLIC_INFO("Testing Curve Evaluation");
  axom::quest::NURBSCurveGWNQuery<ExecSpace> gwn_curves {};
  gwn_curves.preprocess(curves, useDirectEvaluation, !useMemoization);
  gwn_curves.query(dc[0], tol);

  axom::quest::NURBSCurveGWNQuery<ExecSpace> gwn_curves_memoized {};
  gwn_curves_memoized.preprocess(curves, useDirectEvaluation, useMemoization);
  gwn_curves_memoized.query(dc[1], tol);

  axom::quest::NURBSCurveGWNQuery<ExecSpace, 0> gwn_curves_fast {};
  gwn_curves_fast.preprocess(curves, !useDirectEvaluation, !useMemoization);
  gwn_curves_fast.query(dc[2], tol);

  axom::quest::NURBSCurveGWNQuery<ExecSpace, 0> gwn_curves_fast_memoized {};
  gwn_curves_fast_memoized.preprocess(curves, !useDirectEvaluation, useMemoization);
  gwn_curves_fast_memoized.query(dc[3], tol);

  SLIC_INFO("Testing Linearization Evaluation");
  axom::quest::PolylineGWNQuery<ExecSpace> gwn_polyline {};
  gwn_polyline.preprocess(&poly_mesh, useDirectEvaluation);
  gwn_polyline.query(dc[4], tol);

  axom::quest::PolylineGWNQuery<ExecSpace, 0> gwn_polyline_fast {};
  gwn_polyline_fast.preprocess(&poly_mesh, !useDirectEvaluation);
  gwn_polyline_fast.query(dc[5], tol);

  // Compare the in-out values between all fields
  const auto *query_mesh = dc[0].GetMesh();
  const auto num_query_points = query_mesh->GetNodalFESpace()->GetNDofs();

  auto &inout_direct = *dc[0].GetField("inout");
  for(int N = 1; N < num_queries; ++N)
  {
    auto &inout_other = *dc[N].GetField("inout");

    for(int i = 0; i < num_query_points; ++i)
    {
      EXPECT_EQ(inout_direct[i], inout_other[i]);
    }
  }
}

#ifdef AXOM_USE_OPENCASCADE
//------------------------------------------------------------------------------
template <typename ExecSpace>
void check_step_file_triangulation()
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
  constexpr int num_queries = 6;
  axom::Array<mfem::DataCollection> dc(0, num_queries);
  for(int i = 0; i < num_queries; ++i)
  {
    dc.emplace_back(axom::fmt::format("gwn_method_{}", i));
    axom::quest::generate_gwn_query_mesh(dc[i],
                                         shape_bbox,
                                         std::vector<double> {},
                                         std::vector<double> {},
                                         resolution,
                                         query_order);
  }

  // Create tolerance object
  axom::primal::WindingTolerances tol;
  constexpr bool useDirectEvaluation = true;
  constexpr bool useMemoization = true;

  //// Run six different kinds of GWN query ////
  // We expect all three fields to return the same values in this case because
  //  of the specific arrangement of query points and triangulation.
  // In general, triangulating the shape can result in different GWN values
  //  for query points near to individual surfaces.

  // Direct
  SLIC_INFO("Testing Direct Evaluation");
  axom::quest::NURBSPatchGWNQuery<ExecSpace> gwn_patches {};
  gwn_patches.preprocess(patches, useDirectEvaluation, !useMemoization);
  gwn_patches.query(dc[0], tol);

  axom::quest::NURBSPatchGWNQuery<ExecSpace> gwn_patches_memoized {};
  gwn_patches_memoized.preprocess(patches, useDirectEvaluation, useMemoization);
  gwn_patches_memoized.query(dc[1], tol);

  axom::quest::NURBSPatchGWNQuery<ExecSpace, 0> gwn_patches_fast {};
  gwn_patches_fast.preprocess(patches, !useDirectEvaluation, !useMemoization);
  gwn_patches_fast.query(dc[2], tol);

  axom::quest::NURBSPatchGWNQuery<ExecSpace, 0> gwn_patches_fast_memoized {};
  gwn_patches_fast_memoized.preprocess(patches, !useDirectEvaluation, useMemoization);
  gwn_patches_fast_memoized.query(dc[3], tol);

  SLIC_INFO("Testing Linearization Evaluation");
  axom::quest::TriangleGWNQuery<ExecSpace> gwn_triangles {};
  gwn_triangles.preprocess(&tri_mesh, useDirectEvaluation);
  gwn_triangles.query(dc[4], tol);

  axom::quest::TriangleGWNQuery<ExecSpace, 0> gwn_triangles_fast {};
  gwn_triangles_fast.preprocess(&tri_mesh, !useDirectEvaluation);
  gwn_triangles_fast.query(dc[5], tol);

  // Compare the in-out values between the three fields
  const auto *query_mesh = dc[0].GetMesh();
  const auto num_query_points = query_mesh->GetNodalFESpace()->GetNDofs();

  auto &inout_direct = *dc[0].GetField("inout");
  for(int N = 1; N < num_queries; ++N)
  {
    auto &inout_other = *dc[N].GetField("inout");

    for(int i = 0; i < num_query_points; ++i)
    {
      EXPECT_EQ(inout_direct[i], inout_other[i]);
    }
  }
}
#endif

//------------------------------------------------------------------------------
TEST(quest_gwn_methods, mfem_mesh_linearization)
{
  check_mfem_mesh_linearization<axom::SEQ_EXEC>();
}

#if defined AXOM_USE_OPENMP && defined(AXOM_USE_RAJA)
TEST(quest_gwn_methods, mfem_mesh_linearization_omp)
{
  check_mfem_mesh_linearization<axom::OMP_EXEC>();
}
#endif

#ifdef AXOM_USE_OPENCASCADE
TEST(quest_gwn_methods, step_file_triangulation)
{
  check_step_file_triangulation<axom::SEQ_EXEC>();
}
#endif

#if defined(AXOM_USE_OPENCASCADE) && defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
TEST(quest_gwn_methods, step_file_triangulation_omp)
{
  check_step_file_triangulation<axom::OMP_EXEC>();
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
