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

#include "opencascade/BRep_Builder.hxx"
#include "opencascade/BRepBuilderAPI_MakeEdge.hxx"
#include "opencascade/BRepBuilderAPI_MakeFace.hxx"
#include "opencascade/BRepBuilderAPI_MakeWire.hxx"
#include "opencascade/BRepLib.hxx"
#include "opencascade/Geom_BSplineSurface.hxx"
#include "opencascade/Geom_CylindricalSurface.hxx"
#include "opencascade/Geom_Plane.hxx"
#include "opencascade/Geom_RectangularTrimmedSurface.hxx"
#include "opencascade/Geom2d_BSplineCurve.hxx"
#include "opencascade/GeomConvert.hxx"
#include "opencascade/STEPControl_Writer.hxx"
#include "opencascade/TColgp_Array1OfPnt2d.hxx"
#include "opencascade/TColStd_Array1OfInteger.hxx"
#include "opencascade/TColStd_Array1OfReal.hxx"
#include "opencascade/TopoDS_Compound.hxx"
#include "opencascade/TopoDS_Face.hxx"
#include "opencascade/gp_Ax1.hxx"
#include "opencascade/gp_Ax3.hxx"
#include "opencascade/gp_Trsf.hxx"

#include <cmath>
#include <string>

namespace fs = axom::utilities::filesystem;
namespace mint = axom::mint;
namespace primal = axom::primal;
namespace quest = axom::quest;

// Enable code to apply test to larger files
//  and check orientation of trimming curves
#define AXOM_TEST_EXPENSIVE_STEP
#define AXOM_TEST_PATCH_ORIENTATION

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
// Use a pared-down version of the TriangleGWN3D method to directly evaluate
//  the winding number for a triangle mesh
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
      [=] AXOM_HOST_DEVICE(axom::IndexType AXOM_UNUSED_PARAM(nodeIdx), double x, double y, double z) {
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
      [=] AXOM_HOST_DEVICE(axom::IndexType cellIdx,
                           const axom::numerics::Matrix<double>& coords,
                           const axom::IndexType* AXOM_UNUSED_PARAM(nodeIds)) {
        // clang-format off
        trisView[cellIdx] = Triangle3D {Point3D {(coords(0, 0) - ctr[0]) / scl, (coords(1, 0) - ctr[1]) / scl, (coords(2, 0) - ctr[2]) / scl},
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
::testing::AssertionResult checkSinglePatchOrientation(const primal::NURBSPatch<double, 3>& patch,
                                                       int idx)
{
  // Sample the net trim winding, including clockwise hole boundaries. Open trims
  // can have fractional winding, so test the sign of the rounded value and report
  // departures from integers separately as closure diagnostics.

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
  double max_integer_error = 0.;
  for(auto wn : gwn_arr)
  {
    max_integer_error = std::max(max_integer_error, std::abs(wn - std::round(wn)));
  }
  SLIC_WARNING_IF(max_integer_error > integer_eps,
                  axom::fmt::format("Patch {} has fractional trim winding; maximum distance to an "
                                    "integer is {}. Check trimming-curve closure.",
                                    idx,
                                    max_integer_error));

  for(int i = 0; i < gwn_arr.size(); ++i)
  {
    if(!(std::round(gwn_arr[i]) >= 0.))
    {
      return ::testing::AssertionFailure()
        << axom::fmt::format("Patch {} has invalid trim orientation at {}: GWN {}",
                             idx,
                             query_arr[i],
                             gwn_arr[i]);
    }
  }
  return ::testing::AssertionSuccess();
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
  const auto bboxDiag = shapeBbox.range();

  axom::Array<primal::Point<double, 3>> query_arr;
  for(const double fx : {0.11, 0.251, 0.52, 0.731, 0.91})
  {
    for(const double fy : {0.09, 0.249, 0.51, 0.75, 0.81})
    {
      for(const double fz : {0.105, 0.25, 0.5, 0.751, 0.901})
      {
        query_arr.push_back(primal::Point<double, 3>({shapeBbox.getMin()[0] + fx * bboxDiag[0],
                                                      shapeBbox.getMin()[1] + fy * bboxDiag[1],
                                                      shapeBbox.getMin()[2] + fz * bboxDiag[2]}));
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
  ASSERT_EQ(0, stepReader.getTriangleMesh(&triMesh, deflection));

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

#ifdef AXOM_TEST_PATCH_ORIENTATION
  SLIC_INFO("-- Testing NURBS read through Patch GWN --");
  for(int i = 0; i < patches.size(); ++i)
  {
    EXPECT_TRUE(checkSinglePatchOrientation(patches[i], i));
  }
#endif
}

TEST(quest_step_reader, trim_orientation_with_open_boundaries)
{
  using Point2D = primal::Point<double, 2>;
  using Curve = primal::NURBSCurve<double, 2>;
  for(double gap : {0., 0.1})
  {
    SCOPED_TRACE(gap);
    primal::NURBSPatch<double, 3> patch(1, 1);
    for(int i = 0; i < 2; ++i)
    {
      for(int j = 0; j < 2; ++j)
      {
        patch(i, j) = {double(i), double(j), 0.};
      }
    }

    const Point2D corners[] = {{0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}, {0., gap}};
    for(int i = 0; i < 4; ++i)
    {
      patch.addTrimmingCurve(Curve(axom::Array<Point2D> {corners[i], corners[i + 1]}, 1));
    }
    const double winding = primal::winding_number(Point2D {0.5, 0.5}, patch.getTrimmingCurves());
    if(gap == 0.)
    {
      EXPECT_NEAR(winding, 1., 1e-8);
    }
    else
    {
      EXPECT_GT(std::abs(winding - std::round(winding)), 1e-3);
    }
    EXPECT_TRUE(checkSinglePatchOrientation(patch, 0));

    for(auto& curve : patch.getTrimmingCurves())
    {
      curve.reverseOrientation();
    }
    EXPECT_FALSE(checkSinglePatchOrientation(patch, 0));
  }
}

//------------------------------------------------------------------------------
// If the STEP file is read properly, then the resulting GWN should be integer-valued
TEST(quest_step_reader, test_tet) { runStepFileTest("tet.step"); }
TEST(quest_step_reader, test_cylinder) { runStepFileTest("sliced_cylinder.step"); }
TEST(quest_step_reader, test_nut) { runStepFileTest("nut.step"); }
TEST(quest_step_reader, test_bearings) { runStepFileTest("bearings.step", 0.05); }
TEST(quest_step_reader, test_bolt_clip) { runStepFileTest("bolt_clip.step"); }

namespace
{
// Every edge uses [0,1], but its stored B-spline pcurve spans [-1,2].
// B-spline pcurves keep their unused tails through the STEP round trip.
TopoDS_Face makeRestrictedPcurveFace()
{
  Handle(Geom_BSplineSurface) surface = GeomConvert::SurfaceToBSplineSurface(
    new Geom_RectangularTrimmedSurface(new Geom_Plane(gp_Ax3()), -2., 3., -2., 3.));
  const double corners[5][2] = {{0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}, {0., 0.}};
  BRepBuilderAPI_MakeWire wire;
  for(int i = 0; i < 4; ++i)
  {
    const auto& a = corners[i];
    const auto& b = corners[i + 1];
    TColgp_Array1OfPnt2d poles(1, 2);
    poles(1) = gp_Pnt2d(2 * a[0] - b[0], 2 * a[1] - b[1]);
    poles(2) = gp_Pnt2d(2 * b[0] - a[0], 2 * b[1] - a[1]);
    TColStd_Array1OfReal knots(1, 2);
    knots(1) = -1.;
    knots(2) = 2.;
    TColStd_Array1OfInteger multiplicities(1, 2);
    multiplicities(1) = multiplicities(2) = 2;
    Handle(Geom2d_BSplineCurve) curve = new Geom2d_BSplineCurve(poles, knots, multiplicities, 1);
    wire.Add(BRepBuilderAPI_MakeEdge(curve, surface, 0., 1.).Edge());
  }
  TopoDS_Face face = BRepBuilderAPI_MakeFace(surface, wire.Wire(), true).Face();
  BRepLib::BuildCurves3d(face);
  return face;
}

void expectPositiveTrimWinding(const primal::NURBSPatch<double, 3>& patch,
                               double inside_u,
                               double inside_v,
                               double outside_u,
                               double outside_v)
{
  axom::Array<primal::Point<double, 2>> queries;
  queries.push_back({inside_u, inside_v});
  queries.push_back({outside_u, outside_v});
  const auto winding = primal::winding_number(queries, patch.getTrimmingCurves());
  ASSERT_EQ(winding.size(), 2);
  EXPECT_NEAR(winding[0], 1., 1e-8);
  EXPECT_NEAR(winding[1], 0., 1e-8);
  EXPECT_TRUE(patch.isVisible(inside_u, inside_v));
  EXPECT_FALSE(patch.isVisible(outside_u, outside_v));
}
}  // namespace

TEST(quest_step_reader, restricted_pcurves_and_located_reversed_face)
{
  for(bool reversed : {false, true})
  {
    SCOPED_TRACE(reversed);
    auto face = makeRestrictedPcurveFace();
    gp_Trsf transform;
    transform.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(1., 0., 0.)), 0.37);
    transform.SetTranslationPart(gp_Vec(3., -2., 5.));
    face.Location(TopLoc_Location(transform));
    if(reversed)
    {
      face.Reverse();
    }

    fs::TempFile file("axom_step_restricted_", "step");
    STEPControl_Writer writer;
    ASSERT_EQ(writer.Transfer(face, STEPControl_AsIs), IFSelect_RetDone);
    ASSERT_EQ(writer.Write(file.getPath().c_str()), IFSelect_RetDone);
    quest::STEPReader reader;
    reader.setFileName(file.getPath());
    ASSERT_EQ(reader.read(false), 0);
    ASSERT_EQ(reader.numPatches(), 1);
    ASSERT_EQ(reader.getPatchIds().size(), 1);
    EXPECT_EQ(reader.getPatchIds()[0], 0);
    const auto& patch = reader.getPatchArray()[0];
    ASSERT_EQ(patch.getNumTrimmingCurves(), 4);
    ASSERT_EQ(reader.getTrimmingCurveWireIds(0).size(), 4);
    for(int wire_id : reader.getTrimmingCurveWireIds(0))
    {
      EXPECT_EQ(wire_id, 0);
    }
    expectPositiveTrimWinding(patch, 0.5, 0.5, 1.5, 0.5);

    for(const auto& curve : patch.getTrimmingCurves())
    {
      for(int i = 0; i <= 8; ++i)
      {
        const auto uv =
          curve.evaluate(curve.getMinKnot() + (curve.getMaxKnot() - curve.getMinKnot()) * i / 8.);
        EXPECT_GE(uv[0], -1e-9);
        EXPECT_LE(uv[0], 1. + 1e-9);
        EXPECT_GE(uv[1], -1e-9);
        EXPECT_LE(uv[1], 1. + 1e-9);
        gp_Pnt expected(reversed ? 1. - uv[0] : uv[0], uv[1], 0.);
        expected.Transform(transform);
        const auto actual = patch.evaluate(uv[0], uv[1]);
        EXPECT_NEAR(actual[0], expected.X(), 1e-8);
        EXPECT_NEAR(actual[1], expected.Y(), 1e-8);
        EXPECT_NEAR(actual[2], expected.Z(), 1e-8);
      }
    }
    gp_Vec expectedNormal(0., 0., reversed ? -1. : 1.);
    expectedNormal.Transform(transform);
    const auto normal = patch.normal(0.5, 0.5);
    const double dot = normal[0] * expectedNormal.X() + normal[1] * expectedNormal.Y() +
      normal[2] * expectedNormal.Z();
    EXPECT_GT(dot, 0.);
  }
}

TEST(quest_step_reader, periodic_faces_and_both_seam_occurrences)
{
  const double period = 2 * std::acos(-1.);
  Handle(Geom_BSplineSurface) surface = GeomConvert::SurfaceToBSplineSurface(
    new Geom_RectangularTrimmedSurface(new Geom_CylindricalSurface(gp_Ax3(), 2.), 0., period, 0., 4.));
  ASSERT_TRUE(surface->IsUPeriodic());
  TopoDS_Compound compound;
  BRep_Builder builder;
  builder.MakeCompound(compound);
  for(double v : {0., 2.})
  {
    // Both source faces share one periodic surface handle.
    auto face = BRepBuilderAPI_MakeFace(surface, 0., period, v, v + 1., 1e-7).Face();
    if(v > 0.)
    {
      face.Reverse();
    }
    builder.Add(compound, face);
  }
  fs::TempFile file("axom_step_periodic_", "step");
  STEPControl_Writer writer;
  ASSERT_EQ(writer.Transfer(compound, STEPControl_AsIs), IFSelect_RetDone);
  ASSERT_EQ(writer.Write(file.getPath().c_str()), IFSelect_RetDone);
  quest::STEPReader reader;
  reader.setFileName(file.getPath());
  ASSERT_EQ(reader.read(false), 0);
  ASSERT_EQ(reader.numPatches(), 2);
  ASSERT_EQ(reader.getPatchIds().size(), 2);
  for(int i = 0; i < 2; ++i)
  {
    SCOPED_TRACE(i);
    EXPECT_TRUE(reader.patchWasOriginallyPeriodic_u(i));
    EXPECT_FALSE(reader.patchWasOriginallyPeriodic_v(i));
    const auto& patch = reader.getPatchArray()[i];
    ASSERT_EQ(patch.getNumTrimmingCurves(), 4);
    ASSERT_EQ(reader.getTrimmingCurveWireIds(i).size(), 4);
    EXPECT_EQ(reader.getPatchIds()[i], i);
    for(int wire_id : reader.getTrimmingCurveWireIds(i))
    {
      EXPECT_EQ(wire_id, 0);
    }
    const double centerV = 2. * i + 0.5;
    expectPositiveTrimWinding(patch, 0.5 * period, centerV, 0.5 * period, centerV + 1.);
    int lowSeams = 0;
    int highSeams = 0;
    for(const auto& curve : patch.getTrimmingCurves())
    {
      const auto a = curve.evaluate(curve.getMinKnot());
      const auto b = curve.evaluate(curve.getMaxKnot());
      lowSeams += std::abs(a[0]) < 1e-8 && std::abs(b[0]) < 1e-8;
      highSeams += std::abs(a[0] - period) < 1e-8 && std::abs(b[0] - period) < 1e-8;
    }
    EXPECT_EQ(lowSeams, 1);
    EXPECT_EQ(highSeams, 1);
    const auto low = patch.evaluate(patch.getMinKnot_u(), centerV);
    const auto high = patch.evaluate(patch.getMaxKnot_u(), centerV);
    EXPECT_NEAR(primal::squared_distance(low, high), 0., 1e-16);
    const auto center = patch.evaluate(0.5 * period, centerV);
    const auto normal = patch.normal(0.5 * period, centerV);
    const double radialDot = center[0] * normal[0] + center[1] * normal[1];
    EXPECT_GT((i == 0 ? 1. : -1.) * radialDot, 0.);
  }
}

#ifdef AXOM_TEST_EXPENSIVE_STEP
// These tests are more expensive and therefore should not be run in the main testing loop,
//  but should still be checked if something changes in STEPReader.cpp
TEST(quest_step_reader, test_boxed_sphere) { runStepFileTest("boxed_sphere.step"); }
TEST(quest_step_reader, test_plate) { runStepFileTest("plate.step", 0.01); }
#endif

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
