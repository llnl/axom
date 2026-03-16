// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/core.hpp"
#include "axom/mint.hpp"
#include "axom/primal.hpp"
#include "axom/quest.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <string>

namespace fs = axom::utilities::filesystem;
namespace mint = axom::mint;
namespace primal = axom::primal;
namespace quest = axom::quest;

//------------------------------------------------------------------------------
std::string pjoin(const std::string& str) { return str; }

std::string pjoin(const char* str) { return std::string(str); }

template <typename... Args>
std::string pjoin(const std::string& str, Args... args)
{
  return fs::joinPath(str, pjoin(args...));
}

template <typename... Args>
std::string pjoin(const char* str, Args... args)
{
  return fs::joinPath(std::string(str), pjoin(args...));
}

//------------------------------------------------------------------------------
bool isNearInteger(double value, double eps)
{
  const double nearest = std::round(value);
  return std::abs(value - nearest) <= eps;
}

//------------------------------------------------------------------------------
// Use a pared-down version of the TriangleGWN3D method to directly evaluate
//  the for a triangle mesh
struct MiniTriangleGWN3D
{
  using Point3D = primal::Point<double, 3>;
  using Triangle3D = primal::Triangle<double, 3>;
  using BBox3D = primal::BoundingBox<double, 3>;

  void preprocess(const mint::UnstructuredMesh<mint::SINGLE_SHAPE>& mesh)
  {
    const auto ntris = mesh.getNumberOfCells();
    m_triangles.resize(ntris);

    BBox3D shapeBBox;
    BBox3D* shapeBBoxPtr = &shapeBBox;
    mint::for_all_nodes<axom::SEQ_EXEC, mint::xargs::xyz>(
      &mesh,
      AXOM_LAMBDA(axom::IndexType /*nodeIdx*/, double x, double y, double z) {
        shapeBBoxPtr->addPoint(Point3D {x, y, z});
      });

    m_shapeCenter = shapeBBox.getCentroid();
    const auto longestDim = shapeBBox.getLongestDimension();
    m_scale = shapeBBox.getMax()[longestDim] - shapeBBox.getMin()[longestDim];

    const auto& ctr = m_shapeCenter;
    const auto scl = m_scale;
    auto trisView = m_triangles.view();
    mint::for_all_cells<axom::SEQ_EXEC, mint::xargs::coords>(
      &mesh,
      AXOM_LAMBDA(axom::IndexType cellIdx,
                  const axom::numerics::Matrix<double>& coords,
                  [[maybe_unused]] const axom::IndexType* nodeIds) {
        // clang-format off
      trisView[cellIdx] =
          Triangle3D {Point3D {(coords(0, 0) - ctr[0]) / scl, (coords(1, 0) - ctr[1]) / scl, (coords(2, 0) - ctr[2]) / scl},
                      Point3D {(coords(0, 1) - ctr[0]) / scl, (coords(1, 1) - ctr[1]) / scl, (coords(2, 1) - ctr[2]) / scl},
                      Point3D {(coords(0, 2) - ctr[0]) / scl, (coords(1, 2) - ctr[1]) / scl, (coords(2, 2) - ctr[2]) / scl}};
        // clang-format on
      });
  }

  axom::Array<double> evaluate(const axom::Array<Point3D>& queryArr,
                               const double edge_tol,
                               const double EPS) const
  {
    axom::Array<double> gwnArr(queryArr.size());
    for(int qi = 0; qi < queryArr.size(); ++qi)
    {
      const Point3D qScaled((queryArr[qi].array() - m_shapeCenter.array()) / m_scale);

      double wn = 0.;
      for(const auto& tri : m_triangles)
      {
        wn += primal::winding_number(qScaled, tri, edge_tol, EPS);
      }

      gwnArr[qi] = wn;
    }

    return gwnArr;
  }

  Point3D m_shapeCenter;
  double m_scale {1.0};
  axom::Array<Triangle3D> m_triangles;
};

//------------------------------------------------------------------------------
void runSinglePatchTest(const axom::primal::NURBSPatch<double, 3>& patch, int idx)
{
  // This function evaluates the 2D GWN in the patch parameter space on a grid of points
  //  to verify that the value is near a nonnegative integer, i.e. the trimming curves
  //  form closed, CCW oriented loops

  constexpr int npts = 10;
  double u_pts[npts], v_pts[npts];
  axom::numerics::linspace(patch.getMinKnot_u() - 0.11, patch.getMaxKnot_u() + 0.102, u_pts, npts);
  axom::numerics::linspace(patch.getMinKnot_v() - 0.12, patch.getMaxKnot_v() + 0.101, v_pts, npts);

  axom::Array<axom::primal::Point<double, 2>> query_arr;
  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      query_arr.push_back(axom::primal::Point<double, 2> {u, v});
    }
  }

  auto gwn_arr = axom::primal::winding_number(query_arr, patch.getTrimmingCurves());

  constexpr double integer_eps = 1e-3;
  for(auto wn : gwn_arr)
  {
    EXPECT_NEAR(wn, std::round(wn), integer_eps)
      << axom::fmt::format("Patch {} has GWN not near integers\n", idx);
    EXPECT_GE(std::round(wn), 0) << axom::fmt::format("Patch {} has GWN not near integers\n", idx);
  }
}

//------------------------------------------------------------------------------
void runStepFileTest(const std::string& stepFile, double deflection = 0.1)
{
  const std::string fileName = pjoin(AXOM_DATA_DIR, "quest", "step", stepFile);
  SLIC_INFO(axom::fmt::format("Testing STEP file '{}'", fileName));

  quest::STEPReader stepReader;
  stepReader.setVerbosity(false);
  stepReader.setFileName(fileName);

  constexpr bool validate = false;
  const int ret = stepReader.read(validate);
  if(ret != 0)
  {
    SLIC_ERROR(axom::fmt::format("Failed to read STEP file '{}'", fileName));
  }

  const auto& patches = stepReader.getPatchArray();
  EXPECT_GT(patches.size(), 0) << "No NURBS patches were extracted from '" << fileName << "'";
  if(patches.empty())
  {
    return;
  }

  const auto shapeBbox = stepReader.getBRepBoundingBox().scale(1.1);
  auto bboxMin = shapeBbox.getMin();
  auto bboxMax = shapeBbox.getMax();
  const auto bboxDiag = bboxMax.array() - bboxMin.array();

  axom::Array<primal::Point<double, 3>> query_arr;
  for(const double fx : {0.11, 0.251, 0.52, 0.731, 0.91})
  {
    for(const double fy : {0.09, 0.249, 0.51, 0.75, 0.81})
    {
      for(const double fz : {0.105, 0.25, 0.5, 0.751, 0.901})
      {
        query_arr.push_back(primal::Point<double, 3>({bboxMin[0] + fx * bboxDiag[0],
                                                      bboxMin[1] + fy * bboxDiag[1],
                                                      bboxMin[2] + fz * bboxDiag[2]}));
      }
    }
  }

  const primal::WindingTolerances tol;
  const auto gwn_arr = primal::winding_number(query_arr,
                                              patches,
                                              tol.edge_tol,
                                              tol.ls_tol,
                                              tol.quad_tol,
                                              tol.disk_size,
                                              tol.EPS);

  EXPECT_EQ(gwn_arr.size(), query_arr.size());

  SLIC_INFO("-- Testing NURBS read through GWN --");
  constexpr double integer_eps = 1e-3;
  for(int i = 0; i < gwn_arr.size() && i < query_arr.size(); ++i)
  {
    const double wn = gwn_arr[i];
    EXPECT_NEAR(wn, std::round(wn), integer_eps);
  }

  mint::UnstructuredMesh<mint::SINGLE_SHAPE> triMesh(3, mint::TRIANGLE);
  stepReader.getTriangleMesh(&triMesh, deflection);

  MiniTriangleGWN3D triEval;
  triEval.preprocess(triMesh);

  SLIC_INFO("-- Testing triangulation through GWN --");
  const auto tri_gwn_arr = triEval.evaluate(query_arr, tol.edge_tol, tol.EPS);
  EXPECT_EQ(tri_gwn_arr.size(), query_arr.size());
  for(int i = 0; i < tri_gwn_arr.size() && i < query_arr.size(); ++i)
  {
    const double wn = tri_gwn_arr[i];
    EXPECT_NEAR(wn, std::round(wn), integer_eps);
  }

  // Applying this function uniformly to every patch in the surface can help point out
  //  some issues while debugging, but will also catch cases where trimming curves don't
  //  form closed loops, which the STEPReader isn't expected to correct.

  /*
  SLIC_INFO("-- Testing NURBS read through Patch GWN --");
  for(int i = 0; i < patches.size(); ++i)
  {
    runSinglePatchTest(patches[i], i);
  }
  */
}

//------------------------------------------------------------------------------
// If the STEP file is read properly, then the resulting GWN should be integer-valued
TEST(quest_step_reader, test_tet) { runStepFileTest("tet.step"); }
TEST(quest_step_reader, test_cylinder) { runStepFileTest("sliced_cylinder.step"); }
TEST(quest_step_reader, test_nut) { runStepFileTest("nut.step"); }
TEST(quest_step_reader, test_bearings) { runStepFileTest("bearings.step", 0.05); }
TEST(quest_step_reader, test_bolt_clip) { runStepFileTest("bolt_clip.step"); }

// These tests are more expensive and therefore should not be run in the main testing loop,
//  but should still be checked if something changes in STEPReader.cpp
//TEST(quest_step_reader, test_boxed_sphere) { runStepFileTest("boxed_sphere.step"); }
//TEST(quest_step_reader, test_plate) { runStepFileTest("plate.step", 0.01); }

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
