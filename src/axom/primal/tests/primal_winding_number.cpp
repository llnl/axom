// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/primal.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/fmt.hpp"

#include "gtest/gtest.h"

// C++ headers
#include <cmath>
#include <iostream>
#include <fstream>
#include <type_traits>

namespace primal = axom::primal;

TEST(primal_winding_number, simple_cases)
{
  // Test points that are straightforwardly "inside" or "outside" the closed shape
  using Point2D = primal::Point<double, 2>;
  using Triangle = primal::Triangle<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;
  using CPolygon = primal::CurvedPolygon<Bezier>;

  double abs_tol = 1e-8;
  double edge_tol = 1e-8;
  double EPS = primal::PRIMAL_TINY;

  // Simple closed shape with cubic edges
  Point2D top_nodes[] = {Point2D {0.0, 0.0},
                         Point2D {0.0, 1.0},
                         Point2D {-1.0, 1.0},
                         Point2D {-1.0, 0.0}};
  Bezier top_curve(top_nodes, 3);

  Point2D bot_nodes[] = {Point2D {-1.0, 0.0},
                         Point2D {-1.0, -1.0},
                         Point2D {0.0, -1.0},
                         Point2D {0.0, 0.0}};
  Bezier bot_curve(bot_nodes, 3);
  Bezier simple_shape_edges[] = {top_curve, bot_curve};
  CPolygon simple_shape(simple_shape_edges, 2);

  axom::Array<Point2D> inside_points, outside_points;

  // Check interior points
  for(int i = 1; i < 7; i++)
  {
    double offset = std::pow(10, -i);

    inside_points.push_back(Point2D {-0.352, 0.72 - offset});
    inside_points.push_back(Point2D {-0.352 - offset, 0.72});

    outside_points.push_back(Point2D {-0.352, 0.72 + offset});
    outside_points.push_back(Point2D {-0.352 + offset, 0.72});
  }

  // Evaluate with memoized algorithm
  auto inside_gwn = winding_number(inside_points, simple_shape, edge_tol, EPS);
  auto outside_gwn = winding_number(outside_points, simple_shape, edge_tol, EPS);

  for(int i = 0; i < inside_points.size(); ++i)
  {
    EXPECT_NEAR(inside_gwn[i], 1.0, abs_tol);
    EXPECT_NEAR(outside_gwn[i], 0.0, abs_tol);
  }

  // Test containment on non-convex shape, where the query point is outside
  //  the control polygon, but interior to the closed Bezier curve
  Point2D cubic_loop_nodes[] = {Point2D {0.0, 0.0},
                                Point2D {2.0, 1.0},
                                Point2D {-1.0, 1.0},
                                Point2D {1.0, 0.0}};
  Bezier cubic_loop(cubic_loop_nodes, 3);

  EXPECT_NEAR(winding_number(Point2D({0.4, 0.21}), cubic_loop, edge_tol, EPS),
              -0.630526441742,
              abs_tol);

  // Test containment on a 2D triangle
  Triangle tri(Point2D {1.0, -1.0}, Point2D {0.5, 2.0}, Point2D {-2.0, 0.5});
  const bool includeBoundary = true;
  for(double y = -2.0; y < 2.0; y += 0.15)
  {
    const auto q = Point2D {0.0, y};
    if(tri.contains(q))
    {
      EXPECT_EQ(winding_number(q, tri, includeBoundary), 1);
    }
    else
    {
      EXPECT_EQ(winding_number(q, tri, includeBoundary), 0);
    }
  }

  // Reverse the orientation, which flips the winding number
  tri = Triangle(Point2D {1.0, -1.0}, Point2D {2.0, 0.5}, Point2D {0.5, -2.0});
  for(double y = -2.0; y < 2.0; y += 0.15)
  {
    const auto q = Point2D {0.0, y};
    if(tri.contains(q))
    {
      EXPECT_EQ(winding_number(q, tri, includeBoundary), -1);
    }
    else
    {
      EXPECT_EQ(winding_number(q, tri, includeBoundary), 0);
    }
  }
}

TEST(primal_winding_number, closure_edge_cases)
{
  // Tests for when query is on the linear closure
  using Point2D = primal::Point<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;
  using Segment = primal::Segment<double, 2>;

  double abs_tol = 1e-8;
  double edge_tol = 1e-8;
  double EPS = primal::PRIMAL_TINY;

  // Test on linear cases
  Segment linear(Point2D {0.0, 0.0}, Point2D {1.0, 1.0});

  EXPECT_NEAR(winding_number(Point2D({-0.45, -0.45}), linear, edge_tol), 0.0, abs_tol);
  EXPECT_NEAR(winding_number(Point2D({1.45, 1.45}), linear, edge_tol), 0.0, abs_tol);

  // Extra tests if initial and terminal tangent lines are collinear
  Point2D quartic_nodes[] = {Point2D {0.1, 0.0},
                             Point2D {1.0, 0.0},
                             Point2D {0.0, 1.0},
                             Point2D {-1.0, 0.0},
                             Point2D {-0.1, 0.0}};
  Bezier quartic(quartic_nodes, 4);

  // Tangent lines in opposite directions
  EXPECT_NEAR(winding_number(Point2D({0, 0}), quartic, edge_tol, EPS), 0.5, abs_tol);
  EXPECT_NEAR(winding_number(Point2D({2.5, 0}), quartic, edge_tol, EPS), 0.0, abs_tol);
  EXPECT_NEAR(winding_number(Point2D({-2.5, 0}), quartic, edge_tol, EPS), 0.0, abs_tol);

  // Tests a potential issue where the query point is treated as being
  //  on the closure, but not on the edge of the approximating polygon.
  for(int i = 2; i < 15; ++i)
  {
    // In all cases, the winding number should be *near* 0.5.
    //  If the tolerances don't match, we would get an "off-by-0.5" error

    double diff = std::pow(10, -i);
    EXPECT_NEAR(winding_number(Point2D({0, diff}), quartic, 0.5 * diff, EPS), 0.5, 0.1);
    EXPECT_NEAR(winding_number(Point2D({0, diff}), quartic, 1.0 * diff, EPS), 0.5, 0.1);
    EXPECT_NEAR(winding_number(Point2D({0, diff}), quartic, 2.0 * diff, EPS), 0.5, 0.1);

    EXPECT_NEAR(winding_number(Point2D({0, -diff}), quartic, 0.5 * diff, EPS), 0.5, 0.1);
    EXPECT_NEAR(winding_number(Point2D({0, -diff}), quartic, 1.0 * diff, EPS), 0.5, 0.1);
    EXPECT_NEAR(winding_number(Point2D({0, -diff}), quartic, 2.0 * diff, EPS), 0.5, 0.1);
  }

  // Flip the curve vertically
  quartic[2] = Point2D({0.0, -1.0});
  EXPECT_NEAR(winding_number(Point2D({0, 0}), quartic, edge_tol, EPS), -0.5, abs_tol);
  EXPECT_NEAR(winding_number(Point2D({2.5, 0}), quartic, edge_tol, EPS), 0.0, abs_tol);
  EXPECT_NEAR(winding_number(Point2D({-2.5, 0}), quartic, edge_tol, EPS), 0.0, abs_tol);

  // Flip one of the tangent lines
  quartic[1] = Point2D({0.0, 0.0});
  EXPECT_NEAR(winding_number(Point2D({0, 0}), quartic, edge_tol, EPS), -0.5, abs_tol);
  EXPECT_NEAR(winding_number(Point2D({2.5, 0}), quartic, edge_tol, EPS), 0.0, abs_tol);
  EXPECT_NEAR(winding_number(Point2D({-2.5, 0}), quartic, edge_tol, EPS), 0.0, abs_tol);

  // Flip vertically again
  quartic[2] = Point2D({0.0, 1.0});
  EXPECT_NEAR(winding_number(Point2D({0, 0}), quartic, edge_tol, EPS), 0.5, abs_tol);
  EXPECT_NEAR(winding_number(Point2D({2.5, 0}), quartic, edge_tol, EPS), 0.0, abs_tol);
  EXPECT_NEAR(winding_number(Point2D({-2.5, 0}), quartic, edge_tol, EPS), 0.0, abs_tol);
}

TEST(primal_winding_number, corner_cases)
{
  // Tests for when query is identically on either endpoint of the Bezier curve.
  //  Conventionally undefined mathematically, we return the limiting value,
  //  which depends on tangent lines at the query point
  using Point2D = primal::Point<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;

  double abs_tol = 1e-8;
  double edge_tol = 1e-8;
  double EPS = primal::PRIMAL_TINY;

  // Line segment
  Point2D linear_nodes[] = {Point2D {0.0, 0.0}, Point2D {1.0, 1.0}};
  Bezier linear(linear_nodes, 1);

  // Cubic curve
  Point2D cubic_nodes[] = {Point2D {0.0, 0.0},
                           Point2D {0.0, 1.0},
                           Point2D {-1.0, 1.0},
                           Point2D {-1.0, 0.0}};
  Bezier cubic(cubic_nodes, 3);

  EXPECT_NEAR(  // Query on endpoint of linear
    winding_number(Point2D({0.0, 0.0}), linear, edge_tol, EPS),
    0.0,
    abs_tol);
  EXPECT_NEAR(  // Query on endpoint of linear
    winding_number(Point2D({1.0, 1.0}), linear, edge_tol, EPS),
    0.0,
    abs_tol);

  EXPECT_NEAR(  // Query on initial endpoint of cubic
    winding_number(Point2D({-1.0, 0.0}), cubic, edge_tol, EPS),
    0.25,
    abs_tol);
  EXPECT_NEAR(  // Query on terminal endpoint of cubic
    winding_number(Point2D({-1.0, 0.0}), cubic, edge_tol, EPS),
    0.25,
    abs_tol);

  // The query is on the endpoint after one bisection
  EXPECT_NEAR(winding_number(Point2D({-0.5, 0.75}), cubic, edge_tol, EPS), 0.312832962673, abs_tol);

  // Query point on both endpoints
  Point2D closed_cubic_nodes[] = {Point2D {0.0, 0.0},
                                  Point2D {2.0, 1.0},
                                  Point2D {-1.0, 1.0},
                                  Point2D {0.0, 0.0}};
  Bezier cubic_closed(closed_cubic_nodes, 3);
  EXPECT_NEAR(winding_number(Point2D({0.0, 0.0}), cubic_closed, edge_tol, EPS),
              0.301208191175,
              abs_tol);

  // Extra tests if initial and terminal tangent lines are collinear
  Point2D quartic_nodes[] = {Point2D {0.1, 0.0},
                             Point2D {1.0, 0.0},
                             Point2D {0.0, 1.0},
                             Point2D {-1.0, 0.0},
                             Point2D {-0.1, 0.0}};
  Bezier quartic(quartic_nodes, 4);

  EXPECT_NEAR(winding_number(Point2D({0.1, 0}), quartic, edge_tol, EPS), 0.5, abs_tol);
  EXPECT_NEAR(winding_number(Point2D({-0.1, 0}), quartic, edge_tol, EPS), 0.5, abs_tol);

  // Flip the curve vertically
  quartic[2] = Point2D({0.0, -1.0});
  EXPECT_NEAR(winding_number(Point2D({0.1, 0}), quartic, edge_tol, EPS), -0.5, abs_tol);
  EXPECT_NEAR(winding_number(Point2D({-0.1, 0}), quartic, edge_tol, EPS), -0.5, abs_tol);

  // Flip one of the tangent lines
  quartic[1] = Point2D({0.0, 0.0});
  EXPECT_NEAR(winding_number(Point2D({0.1, 0}), quartic, edge_tol, EPS), 0, abs_tol);
  EXPECT_NEAR(winding_number(Point2D({-0.1, 0}), quartic, edge_tol, EPS), -0.5, abs_tol);

  // Flip vertically again
  quartic[2] = Point2D({0.0, 1.0});
  EXPECT_NEAR(winding_number(Point2D({0.1, 0}), quartic, edge_tol, EPS), 0, abs_tol);
  EXPECT_NEAR(winding_number(Point2D({-0.1, 0}), quartic, edge_tol, EPS), 0.5, abs_tol);
}

TEST(primal_winding_number, edge_cases)
{
  // Tests for when query is identically on interior of Bezier curve.
  //   Conventionally undefined mathematically. Uses endpoint formulas
  //   to determine value after some number of bisections
  using Point2D = primal::Point<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;

  double abs_tol = 1e-4;
  double edge_tol = 1e-8;
  double EPS = primal::PRIMAL_TINY;

  // Line segment
  Point2D linear_nodes[] = {Point2D {0.0, 0.0}, Point2D {1.0, 1.0}};
  Bezier linear(linear_nodes, 1);

  // At any point on a line, returns 0
  for(double t = 0.1; t < 1; t += 0.1)
  {
    EXPECT_NEAR(winding_number(Point2D({t, t}), linear, edge_tol, EPS), 0.0, abs_tol);
  }

  // Cubic curve, where query is not on an endpoint after any number of bisections
  Point2D cubic_nodes[] = {Point2D {0.0, 0.0},
                           Point2D {0.0, 1.0},
                           Point2D {-1.0, 1.0},
                           Point2D {-1.0, 0.0}};
  Bezier cubic(cubic_nodes, 3);

  EXPECT_NEAR(winding_number(cubic.evaluate(0.1), cubic, edge_tol, EPS), 0.276676361896, abs_tol);
  EXPECT_NEAR(winding_number(cubic.evaluate(0.4), cubic, edge_tol, EPS), 0.310998033871, abs_tol);
  EXPECT_NEAR(winding_number(cubic.evaluate(0.7), cubic, edge_tol, EPS), 0.305165888012, abs_tol);

  // Cubic curve with internal loop
  Point2D cubic_loop_nodes[] = {Point2D {0.0, 0.0},
                                Point2D {2.0, 1.0},
                                Point2D {-1.0, 1.0},
                                Point2D {1.0, 0.0}};
  Bezier cubic_loop(cubic_loop_nodes, 3);
  EXPECT_NEAR(winding_number(Point2D({0.5, 0.3}), cubic_loop, edge_tol, EPS), 0.327979130377, abs_tol);
  EXPECT_NEAR(winding_number(Point2D({0.5, 0.75}), cubic_loop, edge_tol, EPS), 0.687167046798, abs_tol);
}

TEST(primal_winding_number, degenerate_cases)
{
  // Tests for when Bezier curves are defined with duplicate nodes
  using Point2D = primal::Point<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;

  double abs_tol = 1e-8;
  double edge_tol = 1e-8;
  double EPS = primal::PRIMAL_TINY;

  // Flat curve with anti-parallel tangent lines
  Point2D double_bt_nodes[] = {Point2D {0.0, 1.0},
                               Point2D {0.0, 2.0},
                               Point2D {0.0, 0.0},
                               Point2D {0.0, -2.0},
                               Point2D {0.0, -1.0}};
  Bezier double_bt(double_bt_nodes, 4);

  Point2D linear_nodes[] = {Point2D {0.0, 1.0}, Point2D {0.0, -1.0}};
  Bezier linear(linear_nodes, 1);

  for(double t = -3.0; t <= 3.0; t += 0.1)
  {
    EXPECT_NEAR(winding_number(Point2D({0.0, t}), double_bt, edge_tol, EPS), 0.0, abs_tol);
    EXPECT_NEAR(winding_number(Point2D({1.0, t}), double_bt, edge_tol, EPS),
                winding_number(Point2D({1.0, t}), linear, edge_tol, EPS),
                abs_tol);
  }

  // Check endpoints specifically
  EXPECT_NEAR(winding_number(Point2D({0.0, 1.0}), double_bt, edge_tol, EPS), 0.0, abs_tol);
  EXPECT_NEAR(winding_number(Point2D({0.0, -1.0}), double_bt, edge_tol, EPS), 0.0, abs_tol);

  // empty curve, high order.
  Point2D empty_nodes[] = {Point2D {0.0, 0.0},
                           Point2D {0.0, 0.0},
                           Point2D {0.0, 0.0},
                           Point2D {0.0, 0.0}};
  Bezier empty_curve(empty_nodes, 3);
  axom::Array<Point2D> test_points = {Point2D {0.0, 0.0},
                                      Point2D {0.0, 1.0},
                                      Point2D {1.0, 0.0},
                                      Point2D {1.0, 1.0},
                                      Point2D {0.0, 1.0}};

  for(auto pt : test_points)
  {
    EXPECT_NEAR(winding_number(pt, empty_curve, edge_tol, EPS), 0, abs_tol);
  }

  // Check default empty Bezier curves
  Bezier very_empty_curve(-1);
  for(auto pt : test_points)
  {
    EXPECT_NEAR(winding_number(pt, very_empty_curve, edge_tol, EPS), 0, abs_tol);
  }

  very_empty_curve.setOrder(0);
  for(auto pt : test_points)
  {
    EXPECT_NEAR(winding_number(pt, very_empty_curve, edge_tol, EPS), 0, abs_tol);
  }

  // Cubic curve with many duplicated endpoints
  Point2D cubic_nodes[] = {Point2D {0.0, 0.0},
                           Point2D {0.0, 0.0},
                           Point2D {0.0, 0.0},
                           Point2D {0.0, 1.0},
                           Point2D {-1.0, 1.0},
                           Point2D {-1.0, 0.0},
                           Point2D {-1.0, 0.0},
                           Point2D {-1.0, 0.0}};
  Bezier cubic(cubic_nodes, 7);

  EXPECT_NEAR(  // Query on initial endpoint of cubic
    winding_number(Point2D({-1.0, 0.0}), cubic, edge_tol, EPS),
    0.25,
    abs_tol);
  EXPECT_NEAR(  // Query on terminal endpoint of cubic
    winding_number(Point2D({-1.0, 0.0}), cubic, edge_tol, EPS),
    0.25,
    abs_tol);
}

TEST(primal_winding_number, rational_bezier_winding_number)
{
  using Point2D = primal::Point<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;
  using CPolygon = primal::CurvedPolygon<Bezier>;

  double abs_tol = 1e-8;
  double edge_tol = 0;
  double EPS = 0;

  // Simple quarter circle shape
  Point2D circle_nodes[] = {Point2D {1.0, 0.0}, Point2D {1.0, 1.0}, Point2D {0.0, 1.0}};
  double weights[] = {2.0, 1.0, 1.0};
  Bezier circle_arc(circle_nodes, weights, 2);

  Point2D leg1_nodes[] = {Point2D {0.0, 1.0}, {0.0, 0.0}};
  Bezier leg1(leg1_nodes, 1);

  Point2D leg2_nodes[] = {Point2D {0.0, 0.0}, {1.0, 0.0}};
  Bezier leg2(leg2_nodes, 1);

  CPolygon quarter_circle;
  quarter_circle.addEdge(circle_arc);
  quarter_circle.addEdge(leg1);
  quarter_circle.addEdge(leg2);

  for(double theta = 0.01; theta < 1.5; theta += 0.05)
  {
    for(int i = 1; i < 9; i++)
    {
      const double offset = std::pow(10, -i);
      const double ri = 1.0 - offset;
      const double ro = 1.0 + offset;

      EXPECT_NEAR(winding_number(Point2D({ri * std::cos(theta), ri * std::sin(theta)}),
                                 quarter_circle,
                                 edge_tol,
                                 EPS),
                  1.0,
                  abs_tol);

      EXPECT_NEAR(winding_number(Point2D({ro * std::cos(theta), ro * std::sin(theta)}),
                                 quarter_circle,
                                 edge_tol,
                                 EPS),
                  0.0,
                  abs_tol);
    }
  }
}

TEST(primal_winding_number, nurbs_winding_numbers)
{
  // Define a nurbs curve that represents a circle
  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  PointType data[7] = {PointType {1.0, 0.0},
                       PointType {1.0, 2.0},
                       PointType {-1.0, 2.0},
                       PointType {-1.0, 0.0},
                       PointType {-1.0, -2.0},
                       PointType {1.0, -2.0},
                       PointType {1.0, 0.0}};
  double weights[7] = {1.0, 1. / 3., 1. / 3., 1.0, 1. / 3., 1. / 3., 1.0};

  double knots[11] = {0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0};

  NURBSCurveType circle(data, weights, 7, knots, 11);

  // Check the winding number on a simple grid of points
  for(double x = -2.0; x <= 2.0; x += 0.21)
  {
    for(double y = -2.0; y <= 2.0; y += 0.21)
    {
      PointType query {x, y};
      double gwn = winding_number(query, circle);

      if(x * x + y * y > 1.0)
      {
        EXPECT_DOUBLE_EQ(std::lround(gwn), 0.0);
      }
      else
      {
        EXPECT_DOUBLE_EQ(std::lround(gwn), 1.0);
      }
    }
  }
}

TEST(primal_winding_number, nurbs_patch_gwn_cache_accessors_are_zero_copy)
{
  using Patch = primal::NURBSPatch<double, 3>;
  using Cache = primal::detail::NURBSPatchGWNCache<double>;

  static_assert(std::is_same_v<decltype(std::declval<const Cache&>().getControlPoints()),
                               decltype(std::declval<const Patch&>().getControlPoints())>,
                "NURBSPatchGWNCache::getControlPoints() must match NURBSPatch signature");
  static_assert(std::is_same_v<decltype(std::declval<const Cache&>().getWeights()),
                               decltype(std::declval<const Patch&>().getWeights())>,
                "NURBSPatchGWNCache::getWeights() must match NURBSPatch signature");
  static_assert(std::is_same_v<decltype(std::declval<const Cache&>().getKnots_u()),
                               decltype(std::declval<const Patch&>().getKnots_u())>,
                "NURBSPatchGWNCache::getKnots_u() must match NURBSPatch signature");
  static_assert(std::is_same_v<decltype(std::declval<const Cache&>().getKnots_v()),
                               decltype(std::declval<const Patch&>().getKnots_v())>,
                "NURBSPatchGWNCache::getKnots_v() must match NURBSPatch signature");
  static_assert(std::is_same_v<decltype(std::declval<const Cache&>().getTrimmingCurves()),
                               decltype(std::declval<const Patch&>().getTrimmingCurves())>,
                "NURBSPatchGWNCache::getTrimmingCurves() must match NURBSPatch signature");
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  int result = RUN_ALL_TESTS();

  return result;
}
