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

#include "opencascade/BRep_Tool.hxx"
#include "opencascade/BRepTools_WireExplorer.hxx"
#include "opencascade/GeomAPI_ProjectPointOnCurve.hxx"
#include "opencascade/STEPControl_Reader.hxx"
#include "opencascade/TopExp_Explorer.hxx"
#include "opencascade/TopoDS.hxx"
#include "opencascade/TopoDS_Wire.hxx"
#include "opencascade/BRep_Builder.hxx"
#include "opencascade/BRepBuilderAPI_MakeEdge.hxx"
#include "opencascade/BRepBuilderAPI_MakeFace.hxx"
#include "opencascade/BRepBuilderAPI_MakeWire.hxx"
#include "opencascade/BRepBuilderAPI_MakeVertex.hxx"
#include "opencascade/BRepCheck_Analyzer.hxx"
#include "opencascade/BRepLib.hxx"
#include "opencascade/Geom_BSplineSurface.hxx"
#include "opencascade/Geom_BSplineCurve.hxx"
#include "opencascade/Geom_CylindricalSurface.hxx"
#include "opencascade/Geom_ConicalSurface.hxx"
#include "opencascade/Geom_Plane.hxx"
#include "opencascade/Geom_RectangularTrimmedSurface.hxx"
#include "opencascade/Geom2d_BSplineCurve.hxx"
#include "opencascade/Geom2d_Line.hxx"
#include "opencascade/GeomConvert.hxx"
#include "opencascade/STEPControl_Writer.hxx"
#include "opencascade/TColgp_Array1OfPnt2d.hxx"
#include "opencascade/TColgp_Array1OfPnt.hxx"
#include "opencascade/TColgp_Array2OfPnt.hxx"
#include "opencascade/TColStd_Array1OfInteger.hxx"
#include "opencascade/TColStd_Array1OfReal.hxx"
#include "opencascade/TopoDS_Compound.hxx"
#include "opencascade/TopoDS_Face.hxx"
#include "opencascade/TopoDS_Vertex.hxx"
#include "opencascade/gp_Dir2d.hxx"
#include "opencascade/gp_Ax1.hxx"
#include "opencascade/gp_Ax3.hxx"
#include "opencascade/gp_Trsf.hxx"

#include <array>
#include <cmath>
#include <limits>
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
void writeStepFile(const TopoDS_Shape& shape, const std::string& path)
{
  STEPControl_Writer writer;
  ASSERT_EQ(writer.Transfer(shape, STEPControl_AsIs), IFSelect_RetDone);
  ASSERT_EQ(writer.Write(path.c_str()), IFSelect_RetDone) << path;
}

// Import the source geometry without going through Axom's conversion and export.
void readStepReference(const std::string& path, STEPControl_Reader& reader)
{
  ASSERT_EQ(reader.ReadFile(path.c_str()), IFSelect_RetDone) << path;
  ASSERT_GT(reader.TransferRoots(), 0) << path;
}

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

// Compare each occurrence with its source edge, including direction and closure.
void expectSingleWireTrimGeometry(const quest::STEPReader& reader,
                                  int patchIndex,
                                  TopoDS_Face sourceFace)
{
  const bool faceReversed = sourceFace.Orientation() == TopAbs_REVERSED;
  sourceFace.Orientation(TopAbs_FORWARD);
  struct SourceEdge
  {
    Handle(Geom_Curve) curve;
    double first, last, tolerance;
    bool reversed;
  };
  axom::Array<SourceEdge> originals;
  for(TopExp_Explorer w(sourceFace, TopAbs_WIRE); w.More(); w.Next())
  {
    for(BRepTools_WireExplorer e(TopoDS::Wire(w.Current()), sourceFace); e.More(); e.Next())
    {
      double a, b;
      const auto curve = BRep_Tool::Curve(e.Current(), a, b);
      ASSERT_FALSE(curve.IsNull());
      originals.push_back(
        {curve,
         a,
         b,
         std::max({1e-7, BRep_Tool::Tolerance(sourceFace), BRep_Tool::Tolerance(e.Current())}),
         (e.Current().Orientation() == TopAbs_REVERSED) != faceReversed});
    }
  }
  const auto& patch = reader.getPatchArray()[patchIndex];
  ASSERT_EQ(originals.size(), patch.getNumTrimmingCurves());

  primal::BoundingBox<double, 2> domain({patch.getMinKnot_u(), patch.getMinKnot_v()},
                                        {patch.getMaxKnot_u(), patch.getMaxKnot_v()});
  domain.expand(1e-7);
  ASSERT_EQ(reader.getTrimmingCurveWireIds(patchIndex).size(), patch.getNumTrimmingCurves());
  for(int i = 0; i < patch.getNumTrimmingCurves(); ++i)
  {
    SCOPED_TRACE(i);
    EXPECT_EQ(reader.getTrimmingCurveWireIds(patchIndex)[i], 0);
    const auto& curve = patch.getTrimmingCurve(i);
    const auto& edge = originals[i];
    constexpr int numSamples = 258;
    constexpr int samplesPerSpan = 18;
    axom::Array<double> parameters(numSamples);
    axom::numerics::linspace(curve.getMinKnot(), curve.getMaxKnot(), parameters.data(), numSamples);
    for(int k = 0; k + 1 < curve.getNumKnots(); ++k)
    {
      const double a = curve.getKnot(k), b = curve.getKnot(k + 1);
      if(b > a)
      {
        double spanSamples[samplesPerSpan];
        axom::numerics::linspace(a, b, spanSamples, samplesPerSpan);
        parameters.append(axom::ArrayView<double>(spanSamples, samplesPerSpan));
      }
    }
    for(double t : parameters)
    {
      const auto uv = curve.evaluate(t);
      EXPECT_TRUE(domain.contains(uv)) << uv << " outside " << domain;
      const auto xyz = patch.evaluate(uv[0], uv[1]);
      const gp_Pnt point(xyz[0], xyz[1], xyz[2]);
      GeomAPI_ProjectPointOnCurve projection(point, edge.curve, edge.first, edge.last);
      double error = std::min(point.Distance(edge.curve->Value(edge.first)),
                              point.Distance(edge.curve->Value(edge.last)));
      if(projection.NbPoints() > 0)
      {
        error = std::min(error, projection.LowerDistance());
      }
      EXPECT_LE(error, edge.tolerance);
    }
    for(bool atEnd : {false, true})
    {
      const auto uv = curve.evaluate(atEnd ? curve.getMaxKnot() : curve.getMinKnot());
      const auto xyz = patch.evaluate(uv[0], uv[1]);
      const auto expected = edge.curve->Value(atEnd != edge.reversed ? edge.last : edge.first);
      EXPECT_LE(gp_Pnt(xyz[0], xyz[1], xyz[2]).Distance(expected), edge.tolerance);
    }
    const auto end = curve.evaluate(faceReversed ? curve.getMinKnot() : curve.getMaxKnot());
    const auto& next = patch.getTrimmingCurve((i + 1) % patch.getNumTrimmingCurves());
    const auto start = next.evaluate(faceReversed ? next.getMaxKnot() : next.getMinKnot());
    EXPECT_LT(primal::squared_distance(end, start), 1e-14);
  }
}

// S(u,v) = (0.004u,v,0) on [0,1]^2. A UV gap of 2e-5 is only 8e-8 in 3D.
// The physical gap stays within Open Cascade's vertex tolerance
// while the UV gap exceeds the diagnostic threshold.
TopoDS_Face makeTrimDiagnosticFace(double topControlV, double gap, bool addHole = false)
{
  TColgp_Array2OfPnt poles(1, 2, 1, 2);
  for(int i = 1; i <= 2; ++i)
  {
    for(int j = 1; j <= 2; ++j)
    {
      poles(i, j) = gp_Pnt(0.004 * (i - 1), j - 1, 0.);
    }
  }
  TColStd_Array1OfReal knots(1, 2);
  knots(1) = 0.;
  knots(2) = 1.;
  TColStd_Array1OfInteger mults(1, 2);
  mults(1) = mults(2) = 2;
  Handle(Geom_BSplineSurface) surface =
    new Geom_BSplineSurface(poles, knots, knots, mults, mults, 1, 1);
  BRep_Builder builder;
  TopoDS_Face face;
  builder.MakeFace(face, surface, 1e-7);
  auto makeWire = [&](double lo, double hi, double midV, double gapUV) {
    const std::array<gp_Pnt2d, 4> corners {gp_Pnt2d(lo, lo),
                                           gp_Pnt2d(hi, lo),
                                           gp_Pnt2d(hi, hi),
                                           gp_Pnt2d(lo, hi)};
    std::array<TopoDS_Vertex, 4> vertices;
    for(int i = 0; i < 4; ++i)
    {
      vertices[i] =
        BRepBuilderAPI_MakeVertex(surface->Value(corners[i].X(), corners[i].Y())).Vertex();
      builder.UpdateVertex(vertices[i], 1e-7);
    }
    TopoDS_Wire wire;
    builder.MakeWire(wire);
    for(int i = 0; i < 4; ++i)
    {
      const int degree = i == 2 ? 2 : 1;
      TColgp_Array1OfPnt2d uv(1, degree + 1);
      uv(1) = corners[i];
      uv(degree + 1) = corners[(i + 1) % 4];
      if(i == 1)
      {
        uv(1).SetX(uv(1).X() + gapUV);
      }
      if(degree == 2)
      {
        uv(2) = gp_Pnt2d((lo + hi) * .5, midV);
      }
      // This affine surface maps the pcurve poles to an exact 3D curve with
      // the same parameterization, so SameRange and SameParameter hold below.
      TColgp_Array1OfPnt xyz(1, degree + 1);
      for(int j = 1; j <= degree + 1; ++j)
      {
        xyz(j) = surface->Value(uv(j).X(), uv(j).Y());
      }
      mults(1) = mults(2) = degree + 1;
      Handle(Geom2d_BSplineCurve) pcurve = new Geom2d_BSplineCurve(uv, knots, mults, degree);
      Handle(Geom_BSplineCurve) spaceCurve = new Geom_BSplineCurve(xyz, knots, mults, degree);
      auto edge =
        BRepBuilderAPI_MakeEdge(spaceCurve, vertices[i], vertices[(i + 1) % 4], 0., 1.).Edge();
      builder.UpdateEdge(edge, pcurve, face, 1e-7);
      builder.Range(edge, face, 0., 1.);
      builder.SameRange(edge, true);
      builder.SameParameter(edge, true);
      builder.Add(wire, edge);
    }
    wire.Closed(true);
    return wire;
  };
  builder.Add(face, makeWire(.2, .8, topControlV, addHole ? 0. : gap));
  if(addHole)
  {
    builder.Add(face, makeWire(.35, .65, .65, gap).Reversed());
  }
  return face;
}

std::string readWithDiagnostics(quest::STEPReader& reader, bool validate)
{
  ::testing::internal::CaptureStdout();
  int status = -1;
  EXPECT_NO_THROW(status = reader.read(validate));
  axom::slic::flushStreams();
  const auto messages = ::testing::internal::GetCapturedStdout();
  EXPECT_EQ(status, 0) << messages;
  return messages;
}

void expectTrimGaps(const quest::STEPReader::NURBSPatch& patch,
                    axom::ArrayView<const int> wire_ids,
                    double expected_uv_gap,
                    double expected_physical_gap)
{
  using Vector3D = primal::Vector<double, 3>;
  constexpr double eps = 1e-12;
  constexpr double source_tol = 1e-7;

  ASSERT_EQ(wire_ids.size(), patch.getNumTrimmingCurves());
  double max_squared_gap = 0.0;
  double max_physical_gap = 0.0;
  // Match endpoints within each wire without assuming traversal order.
  // This avoids reusing STEPReader's consecutive-curve gap calculation.
  for(int i = 0; i < patch.getNumTrimmingCurves(); ++i)
  {
    const auto& curve = patch.getTrimmingCurve(i);
    const auto end = curve.evaluate(curve.getMaxKnot());
    const auto end_point = patch.evaluate(end[0], end[1]);
    double nearest_squared_gap = std::numeric_limits<double>::max();
    double physical_gap = 0.0;
    for(int j = 0; j < patch.getNumTrimmingCurves(); ++j)
    {
      if(wire_ids[i] == wire_ids[j])
      {
        const auto& next = patch.getTrimmingCurve(j);
        const auto start = next.evaluate(next.getMinKnot());
        const double squared_gap = primal::squared_distance(end, start);
        if(squared_gap < nearest_squared_gap)
        {
          nearest_squared_gap = squared_gap;
          physical_gap = Vector3D(end_point, patch.evaluate(start[0], start[1])).norm();
        }
      }
    }
    max_squared_gap = std::max(max_squared_gap, nearest_squared_gap);
    max_physical_gap = std::max(max_physical_gap, physical_gap);
  }
  EXPECT_NEAR(std::sqrt(max_squared_gap), expected_uv_gap, eps);
  EXPECT_NEAR(max_physical_gap, expected_physical_gap, eps);
  EXPECT_LT(max_physical_gap, source_tol);
}

}  // namespace

TEST(quest_step_reader, read_failure)
{
  fs::TempFile empty_file("axom_step_empty_", "step");
  const axom::Array<std::string> paths {empty_file.getPath(), empty_file.getPath() + ".missing"};
  for(const auto& path : paths)
  {
    SCOPED_TRACE(path);
    for(bool validate : {false, true})
    {
      SCOPED_TRACE(validate);
      quest::STEPReader reader;
      reader.setFileName(path);
      EXPECT_NE(reader.read(validate), 0);
      EXPECT_TRUE(reader.getPatchArray().empty());
      EXPECT_TRUE(reader.getPatchIds().empty());
      EXPECT_FALSE(reader.getBRepBoundingBox().isValid());
      mint::UnstructuredMesh<mint::SINGLE_SHAPE> mesh(3, mint::TRIANGLE);
      EXPECT_NE(reader.getTriangleMesh(&mesh), 0);
      EXPECT_EQ(mesh.getNumberOfCells(), 0);
    }
  }
}

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
    ASSERT_NO_FATAL_FAILURE(writeStepFile(face, file.getPath()));
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
    const primal::Vector<double, 3> expected_normal {expectedNormal.X(),
                                                     expectedNormal.Y(),
                                                     expectedNormal.Z()};
    const auto normal = patch.normal(0.5, 0.5);
    EXPECT_GT(normal.dot(expected_normal), 0.);

    // Both triangulation modes must preserve the face's location and orientation.
    // The trim encloses a unit square; the untrimmed support spans [-2,3]^2.
    for(bool trimmed : {true, false})
    {
      SCOPED_TRACE(axom::fmt::format("trimmed={}", trimmed));
      mint::UnstructuredMesh<mint::SINGLE_SHAPE> mesh(3, mint::TRIANGLE);
      ASSERT_EQ(reader.getTriangleMesh(&mesh, 0.1, 0.5, false, trimmed), 0);
      ASSERT_GT(mesh.getNumberOfCells(), 0);
      const auto* patch_ids = mesh.getFieldPtr<int>("patch_index", mint::CELL_CENTERED);
      ASSERT_NE(patch_ids, nullptr);
      double area = 0.0;
      for(axom::IndexType cell = 0; cell < mesh.getNumberOfCells(); ++cell)
      {
        EXPECT_EQ(patch_ids[cell], 0);
        axom::IndexType nodes[3];
        mesh.getCellNodeIDs(cell, nodes);
        primal::Triangle<double, 3> triangle;
        for(int j = 0; j < 3; ++j)
        {
          mesh.getNode(nodes[j], triangle[j].data());
          const auto& point = triangle[j];
          const auto local = gp_Pnt(point[0], point[1], point[2]).Transformed(transform.Inverted());
          EXPECT_NEAR(local.Z(), 0.0, 1e-8);
          EXPECT_GE(local.X(), (trimmed ? 0.0 : -2.0) - 1e-8);
          EXPECT_LE(local.X(), (trimmed ? 1.0 : 3.0) + 1e-8);
          EXPECT_GE(local.Y(), (trimmed ? 0.0 : -2.0) - 1e-8);
          EXPECT_LE(local.Y(), (trimmed ? 1.0 : 3.0) + 1e-8);
        }
        EXPECT_GT(triangle.normal().dot(expected_normal), 0.0);
        area += triangle.area();
      }
      EXPECT_NEAR(area, trimmed ? 1.0 : 25.0, 1e-8);
    }
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
  ASSERT_NO_FATAL_FAILURE(writeStepFile(compound, file.getPath()));
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

// Exercise seam crossings in U and V, on translated domains and reversed faces.
TEST(quest_step_reader, periodic_bspline_seam_crossing)
{
  const double period = 2 * std::acos(-1.);
  for(bool periodic_u : {true, false})
  {
    for(bool reversed : {false, true})
    {
      for(double seam : {-period, 2 * period})
      {
        SCOPED_TRACE(
          axom::fmt::format("periodic_u={} reversed={} seam={}", periodic_u, reversed, seam));
        Handle(Geom_BSplineSurface) surface = GeomConvert::SurfaceToBSplineSurface(
          new Geom_RectangularTrimmedSurface(new Geom_CylindricalSurface(gp_Ax3(), 2.),
                                             0.,
                                             period,
                                             0.,
                                             2.));
        if(!periodic_u)
        {
          surface->ExchangeUV();
        }
        const double u0 = periodic_u ? seam - 0.25 : 0.25;
        const double u1 = periodic_u ? seam + 0.25 : 1.25;
        const double v0 = periodic_u ? 0.25 : seam - 0.25;
        const double v1 = periodic_u ? 1.25 : seam + 0.25;
        auto face = BRepBuilderAPI_MakeFace(surface, u0, u1, v0, v1, 1e-7).Face();
        if(reversed)
        {
          face.Reverse();
        }
        fs::TempFile file("axom_step_seam_", "step");
        ASSERT_NO_FATAL_FAILURE(writeStepFile(face, file.getPath()));
        quest::STEPReader reader;
        reader.setFileName(file.getPath());
        const auto messages = readWithDiagnostics(reader, true);
        EXPECT_EQ(messages.find("leave the exported UV domain"), std::string::npos) << messages;
        EXPECT_EQ(messages.find("UV endpoint gap"), std::string::npos) << messages;
        ASSERT_EQ(reader.numPatches(), 1);
        ASSERT_EQ(reader.getPatchIds().size(), 1);
        EXPECT_EQ(reader.getPatchIds()[0], 0);
        EXPECT_EQ(reader.patchWasOriginallyPeriodic_u(0), periodic_u);
        EXPECT_EQ(reader.patchWasOriginallyPeriodic_v(0), !periodic_u);
        const auto& patch = reader.getPatchArray()[0];
        ASSERT_EQ(patch.getNumTrimmingCurves(), 4);
        ASSERT_EQ(reader.getTrimmingCurveWireIds(0).size(), 4);
        for(int wire_id : reader.getTrimmingCurveWireIds(0))
        {
          EXPECT_EQ(wire_id, 0);
        }
        primal::BoundingBox<double, 2> trimBounds;
        primal::BoundingBox<double, 2> domain({patch.getMinKnot_u(), patch.getMinKnot_v()},
                                              {patch.getMaxKnot_u(), patch.getMaxKnot_v()});
        domain.expand(1e-8);
        for(const auto& curve : patch.getTrimmingCurves())
        {
          constexpr int numSamples = 129;
          double parameters[numSamples];
          axom::numerics::linspace(curve.getMinKnot(), curve.getMaxKnot(), parameters, numSamples);
          for(double t : parameters)
          {
            const auto uv = curve.evaluate(t);
            trimBounds.addPoint(uv);
            EXPECT_TRUE(domain.contains(uv)) << uv << " outside " << domain;
            const auto actual = patch.evaluate(uv[0], uv[1]);
            const double source_u =
              reversed ? patch.getMinKnot_u() + patch.getMaxKnot_u() - uv[0] : uv[0];
            const auto expected = surface->Value(source_u, uv[1]);
            EXPECT_NEAR(actual[0], expected.X(), 1e-8);
            EXPECT_NEAR(actual[1], expected.Y(), 1e-8);
            EXPECT_NEAR(actual[2], expected.Z(), 1e-8);
            EXPECT_GT(actual[0], 1.8);
            EXPECT_GE(actual[2], 0.25 - 1e-8);
            EXPECT_LE(actual[2], 1.25 + 1e-8);
          }
        }
        const auto center = trimBounds.getCentroid();
        expectPositiveTrimWinding(patch,
                                  center[0],
                                  center[1],
                                  periodic_u ? center[0] : trimBounds.getMax()[0] + 0.25,
                                  periodic_u ? trimBounds.getMax()[1] + 0.25 : center[1]);
      }
    }
  }
}

TEST(quest_step_reader, periodic_analytic_cone)
{
  const double period = 2 * std::acos(-1.);
  for(bool reversed : {false, true})
  {
    SCOPED_TRACE(reversed);
    Handle(Geom_ConicalSurface) surface = new Geom_ConicalSurface(gp_Ax3(), 0.35, 3.);
    const double cut = -0.2 * period;
    BRep_Builder builder;
    TopoDS_Face face;
    builder.MakeFace(face, surface, 1e-7);
    const auto bottom = BRepBuilderAPI_MakeVertex(surface->Value(cut, 0.5)).Vertex();
    const auto top = BRepBuilderAPI_MakeVertex(surface->Value(cut, 2.5)).Vertex();
    auto seam = BRepBuilderAPI_MakeEdge(surface->UIso(cut), bottom, top, 0.5, 2.5).Edge();
    builder.UpdateEdge(seam,
                       new Geom2d_Line(gp_Pnt2d(cut + period, 0.), gp_Dir2d(0., 1.)),
                       new Geom2d_Line(gp_Pnt2d(cut, 0.), gp_Dir2d(0., 1.)),
                       face,
                       1e-7);
    builder.Range(seam, face, 0.5, 2.5);

    // Split each rim at its circle's origin and reuse both seam occurrences.
    auto splitRim = [&](double v, const TopoDS_Vertex& vertex) {
      const auto circle = surface->VIso(v);
      const auto marker = BRepBuilderAPI_MakeVertex(surface->Value(0., v)).Vertex();
      const double split = cut + period;
      std::array<TopoDS_Edge, 2> edges {
        BRepBuilderAPI_MakeEdge(circle, vertex, marker, split, period).Edge(),
        BRepBuilderAPI_MakeEdge(circle, marker, vertex, 0., split).Edge()};
      builder.UpdateEdge(edges[0], new Geom2d_Line(gp_Pnt2d(-period, v), gp_Dir2d(1., 0.)), face, 1e-7);
      builder.Range(edges[0], face, split, period);
      builder.UpdateEdge(edges[1], new Geom2d_Line(gp_Pnt2d(0., v), gp_Dir2d(1., 0.)), face, 1e-7);
      builder.Range(edges[1], face, 0., split);
      return edges;
    };
    const auto lower = splitRim(0.5, bottom);
    const auto upper = splitRim(2.5, top);
    BRepBuilderAPI_MakeWire wire;
    wire.Add(lower[0]);
    wire.Add(lower[1]);
    wire.Add(seam);
    wire.Add(TopoDS::Edge(upper[1].Reversed()));
    wire.Add(TopoDS::Edge(upper[0].Reversed()));
    wire.Add(TopoDS::Edge(seam.Reversed()));
    ASSERT_TRUE(wire.IsDone());
    builder.Add(face, wire.Wire());
    ASSERT_TRUE(BRepCheck_Analyzer(face).IsValid());
    if(reversed)
    {
      face.Reverse();
    }
    fs::TempFile file("axom_step_cone_", "step");
    ASSERT_NO_FATAL_FAILURE(writeStepFile(face, file.getPath()));
    STEPControl_Reader reference;
    ASSERT_NO_FATAL_FAILURE(readStepReference(file.getPath(), reference));
    TopExp_Explorer faceExp(reference.OneShape(), TopAbs_FACE);
    ASSERT_TRUE(faceExp.More());

    quest::STEPReader reader;
    reader.setFileName(file.getPath());
    const auto messages = readWithDiagnostics(reader, true);
    EXPECT_EQ(messages.find("leave the exported UV domain"), std::string::npos) << messages;
    EXPECT_EQ(messages.find("UV endpoint gap"), std::string::npos) << messages;
    ASSERT_EQ(reader.numPatches(), 1);
    ASSERT_EQ(reader.getPatchIds().size(), 1);
    EXPECT_EQ(reader.getPatchIds()[0], 0);
    EXPECT_TRUE(reader.patchWasOriginallyPeriodic_u(0));
    EXPECT_FALSE(reader.patchWasOriginallyPeriodic_v(0));
    ASSERT_EQ(reader.getPatchArray()[0].getNumTrimmingCurves(), 6);
    expectSingleWireTrimGeometry(reader, 0, TopoDS::Face(faceExp.Current()));
  }
}

TEST(quest_step_reader, solid56_periodic_face26)
{
  // Gmsh benchmarks/step/Solid_56.step at af3c98576c3d09637b4b6c85f32d7025ea1aa6b8.
  const auto path = pjoin(AXOM_DATA_DIR, "quest", "step", "Solid_56.step");
  if(!fs::pathExists(path))
  {
    GTEST_SKIP() << "Solid_56.step is not available in AXOM_DATA_DIR";
  }
  STEPControl_Reader reference;
  ASSERT_NO_FATAL_FAILURE(readStepReference(path, reference));
  TopExp_Explorer faceExp(reference.OneShape(), TopAbs_FACE);
  for(int i = 0; i < 26 && faceExp.More(); ++i)
  {
    faceExp.Next();
  }
  ASSERT_TRUE(faceExp.More());
  auto sourceFace = TopoDS::Face(faceExp.Current());
  quest::STEPReader reader;
  reader.setFileName(path);
  int status = -1;
  ASSERT_NO_THROW(status = reader.read(false));
  ASSERT_EQ(status, 0);
  ASSERT_EQ(reader.numPatches(), 48);
  ASSERT_EQ(reader.getPatchIds().size(), 48);
  ASSERT_EQ(reader.getPatchIds()[26], 26);
  EXPECT_TRUE(reader.patchWasOriginallyPeriodic_u(26));
  ASSERT_EQ(reader.getPatchArray()[26].getNumTrimmingCurves(), 8);
  expectSingleWireTrimGeometry(reader, 26, sourceFace);
}

TEST(quest_step_reader, stored_periodic_cone_fixtures)
{
  using Vector3D = primal::Vector<double, 3>;
  constexpr double geom_tol = 1e-7;

  // The same wavy conical region, with shifted and native seam cuts.
  for(const char* filename : {"periodic_seam_cone.step", "fixed_seam_offset_cone.step"})
  {
    SCOPED_TRACE(filename);
    const auto path = pjoin(AXOM_DATA_DIR, "quest", "step", filename);
    ASSERT_TRUE(fs::pathExists(path)) << "Missing STEP regression fixture: " << path;
    // Load the source model separately, before Axom's NURBS conversion,
    // and use its 3D edges and tolerances as references for the exported trims below.
    STEPControl_Reader reference;
    ASSERT_NO_FATAL_FAILURE(readStepReference(path, reference));
    ASSERT_TRUE(BRepCheck_Analyzer(reference.OneShape()).IsValid());
    TopExp_Explorer face_exp(reference.OneShape(), TopAbs_FACE);
    ASSERT_TRUE(face_exp.More());
    const auto source_face = TopoDS::Face(face_exp.Current());
    face_exp.Next();
    ASSERT_FALSE(face_exp.More());

    // Both validation modes must preserve the six trim occurrences and their geometry.
    for(bool validate : {false, true})
    {
      SCOPED_TRACE(validate);
      quest::STEPReader reader;
      reader.setFileName(path);
      const auto messages = readWithDiagnostics(reader, validate);
      EXPECT_EQ(messages.find("leave the exported UV domain"), std::string::npos) << messages;
      EXPECT_EQ(messages.find("UV endpoint gap"), std::string::npos) << messages;
      ASSERT_EQ(reader.numPatches(), 1);
      ASSERT_EQ(reader.getPatchIds().size(), 1);
      EXPECT_EQ(reader.getPatchIds()[0], 0);
      EXPECT_TRUE(reader.patchWasOriginallyPeriodic_u(0));
      EXPECT_FALSE(reader.patchWasOriginallyPeriodic_v(0));
      const auto& patch = reader.getPatchArray()[0];
      EXPECT_TRUE(patch.isValidNURBS());
      ASSERT_EQ(patch.getNumTrimmingCurves(), 6);
      // Compare in physical space because the rational patch's U parameter
      // differs from the source cone's angle.
      expectSingleWireTrimGeometry(reader, 0, source_face);
      const double u = (patch.getMinKnot_u() + patch.getMaxKnot_u()) * 0.5;
      const double v = (patch.getMinKnot_v() + patch.getMaxKnot_v()) * 0.5;
      expectPositiveTrimWinding(patch, u, v, u, patch.getMaxKnot_v() + 0.25);
      // The cone has radius 3 at z=0, half-angle 0.35, and an outward normal.
      const auto point = patch.evaluate(u, v);
      const Vector3D radial {point[0], point[1], 0.0};
      EXPECT_NEAR(radial.norm(), 3.0 + point[2] * std::tan(0.35), geom_tol);
      EXPECT_GT(radial.dot(patch.normal(u, v)), 0.0);
    }
  }
}

TEST(quest_step_reader, stored_uv_diagnostic_fixtures)
{
  constexpr double eps = 1e-12;
  constexpr int num_samples = 129;

  // Expected values come from the fixture construction. On the gap fixture's
  // S(u,v) = (0.004u,v,0), a 2e-5 UV gap is 8e-8 in 3D, within its 1e-7 tolerance.
  struct TestCase
  {
    const char* filename;
    double scale_u;
    double max_v;
    double uv_gap;
  };
  const axom::Array<TestCase> test_cases = {{"nonperiodic_trim_excursion.step", 1.0, 1.2, 0.0},
                                            {"trim_gap.step", 0.004, 0.9, 2e-5},
                                            {"closed_control.step", 0.004, 0.9, 0.0}};

  for(const auto& test_case : test_cases)
  {
    SCOPED_TRACE(test_case.filename);
    const auto path = pjoin(AXOM_DATA_DIR, "quest", "step", test_case.filename);
    ASSERT_TRUE(fs::pathExists(path)) << "Missing STEP regression fixture: " << path;
    // Confirm that Open Cascade accepts the source BRep. Axom's exported UV containment
    // and closure checks must still distinguish the excursion, gap, and control.
    STEPControl_Reader reference;
    ASSERT_NO_FATAL_FAILURE(readStepReference(path, reference));
    ASSERT_TRUE(BRepCheck_Analyzer(reference.OneShape()).IsValid());

    for(bool validate : {false, true})
    {
      SCOPED_TRACE(validate);
      quest::STEPReader reader;
      reader.setFileName(path);
      const auto messages = readWithDiagnostics(reader, validate);
      EXPECT_EQ(messages.find("leave the exported UV domain") != std::string::npos,
                validate && test_case.max_v > 1.0)
        << messages;
      EXPECT_EQ(messages.find("UV endpoint gap") != std::string::npos, validate && test_case.uv_gap > 0.0)
        << messages;
      ASSERT_EQ(reader.numPatches(), 1);
      ASSERT_EQ(reader.getPatchIds().size(), 1);
      EXPECT_EQ(reader.getPatchIds()[0], 0);
      EXPECT_FALSE(reader.patchWasOriginallyPeriodic_u(0));
      EXPECT_FALSE(reader.patchWasOriginallyPeriodic_v(0));
      const auto& patch = reader.getPatchArray()[0];
      EXPECT_TRUE(patch.isValidNURBS());
      ASSERT_EQ(patch.getNumTrimmingCurves(), 4);
      const auto wire_ids = reader.getTrimmingCurveWireIds(0);
      ASSERT_EQ(wire_ids.size(), 4);
      for(int wire_id : wire_ids)
      {
        EXPECT_EQ(wire_id, 0);
      }
      EXPECT_DOUBLE_EQ(patch.getMinKnot_u(), 0.0);
      EXPECT_DOUBLE_EQ(patch.getMaxKnot_u(), 1.0);
      EXPECT_DOUBLE_EQ(patch.getMinKnot_v(), 0.0);
      EXPECT_DOUBLE_EQ(patch.getMaxKnot_v(), 1.0);
      for(double u : {0.0, 0.5, 1.0})
      {
        for(double v : {0.0, 0.5, 1.0})
        {
          const auto point = patch.evaluate(u, v);
          EXPECT_NEAR(point[0], test_case.scale_u * u, eps);
          EXPECT_NEAR(point[1], v, eps);
          EXPECT_NEAR(point[2], 0.0, eps);
        }
      }

      // Validation must preserve the excursion and gap. All trim endpoints lie inside
      // the patch domain, including those of the bowed trim.
      const primal::BoundingBox<double, 2> domain({0.0, 0.0}, {1.0, 1.0});
      double max_v = 0.0;
      axom::Array<double> parameters(num_samples);
      for(const auto& curve : patch.getTrimmingCurves())
      {
        EXPECT_TRUE(curve.isValidNURBS());
        EXPECT_TRUE(domain.contains(curve.evaluate(curve.getMinKnot())));
        EXPECT_TRUE(domain.contains(curve.evaluate(curve.getMaxKnot())));
        axom::numerics::linspace(curve.getMinKnot(), curve.getMaxKnot(), parameters.data(), num_samples);
        for(double t : parameters)
        {
          max_v = std::max(max_v, curve.evaluate(t)[1]);
        }
      }
      EXPECT_NEAR(max_v, test_case.max_v, eps);
      expectTrimGaps(patch, wire_ids, test_case.uv_gap, test_case.scale_u * test_case.uv_gap);
    }
  }
}

TEST(quest_step_reader, exported_uv_domain_diagnostics)
{
  for(bool reversed : {false, true})
  {
    // Both curves have control points outside [0,1]^2. Only the second curve leaves it.
    for(double topControlV : {1.1, 1.6})
    {
      SCOPED_TRACE(axom::fmt::format("reversed={} topControlV={}", reversed, topControlV));
      auto face = makeTrimDiagnosticFace(topControlV, 0.);
      if(reversed)
      {
        face.Reverse();
      }
      fs::TempFile file("axom_step_uv_domain_", "step");
      ASSERT_NO_FATAL_FAILURE(writeStepFile(face, file.getPath()));
      for(bool validate : {false, true})
      {
        SCOPED_TRACE(validate);
        quest::STEPReader reader;
        reader.setFileName(file.getPath());
        const auto messages = readWithDiagnostics(reader, validate);
        EXPECT_EQ(messages.find("leave the exported UV domain") != std::string::npos,
                  validate && topControlV == 1.6)
          << messages;
        EXPECT_EQ(messages.find("UV endpoint gap"), std::string::npos) << messages;
        ASSERT_EQ(reader.numPatches(), 1);
        const auto& patch = reader.getPatchArray()[0];
        ASSERT_EQ(patch.getNumTrimmingCurves(), 4);
        EXPECT_EQ(patch.getMinKnot_v(), 0.);
        EXPECT_EQ(patch.getMaxKnot_v(), 1.);
        primal::BoundingBox<double, 2> domain({0., 0.}, {1., 1.});
        double maxV = 0.;
        for(const auto& curve : patch.getTrimmingCurves())
        {
          EXPECT_TRUE(domain.contains(curve.evaluate(curve.getMinKnot())));
          EXPECT_TRUE(domain.contains(curve.evaluate(curve.getMaxKnot())));
          maxV = std::max(maxV, curve.evaluate((curve.getMinKnot() + curve.getMaxKnot()) * .5)[1]);
        }
        EXPECT_NEAR(maxV, (.8 + topControlV) * .5, 1e-12);
      }
    }
  }
}

TEST(quest_step_reader, exported_uv_gap_diagnostics)
{
  for(bool reversed : {false, true})
  {
    for(bool addHole : {false, true})
    {
      for(double gap : {0., 2e-5})
      {
        SCOPED_TRACE(axom::fmt::format("reversed={} addHole={} gap={}", reversed, addHole, gap));
        auto face = makeTrimDiagnosticFace(.8, gap, addHole);
        ASSERT_TRUE(BRepCheck_Analyzer(face).IsValid());
        if(reversed)
        {
          face.Reverse();
        }
        fs::TempFile file("axom_step_uv_gap_", "step");
        ASSERT_NO_FATAL_FAILURE(writeStepFile(face, file.getPath()));
        STEPControl_Reader reference;
        ASSERT_NO_FATAL_FAILURE(readStepReference(file.getPath(), reference));
        ASSERT_TRUE(BRepCheck_Analyzer(reference.OneShape()).IsValid());
        for(bool validate : {false, true})
        {
          SCOPED_TRACE(validate);
          quest::STEPReader reader;
          reader.setFileName(file.getPath());
          const auto messages = readWithDiagnostics(reader, validate);
          EXPECT_EQ(messages.find("UV endpoint gap") != std::string::npos, validate && gap > 0.)
            << messages;
          EXPECT_EQ(messages.find("leave the exported UV domain"), std::string::npos) << messages;
          if(validate && gap > 0.)
          {
            EXPECT_NE(messages.find(axom::fmt::format("[Patch 0 Wire {}] Exported trimming wire",
                                                      addHole ? 1 : 0)),
                      std::string::npos)
              << messages;
          }
          ASSERT_EQ(reader.numPatches(), 1);
          const auto& patch = reader.getPatchArray()[0];
          const auto wireIds = reader.getTrimmingCurveWireIds(0);
          ASSERT_EQ(patch.getNumTrimmingCurves(), addHole ? 8 : 4);
          ASSERT_EQ(wireIds.size(), patch.getNumTrimmingCurves());
          EXPECT_EQ(wireIds[0], 0);
          EXPECT_EQ(wireIds[wireIds.size() - 1], addHole ? 1 : 0);
          expectTrimGaps(patch, wireIds, gap, .004 * gap);
        }
      }
    }
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
