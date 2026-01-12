// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/primal.hpp"
#include "axom/slic.hpp"
#include "axom/fmt.hpp"
#include <iostream>

#include "gtest/gtest.h"

// MFEM includes
#ifdef AXOM_USE_MFEM
  #include "mfem.hpp"
#endif

namespace primal = axom::primal;

TEST(primal_integral, evaluate_area_integral)
{
  using Point2D = primal::Point<double, 2>;
  using BCurve = primal::BezierCurve<double, 2>;
  using CPolygon = primal::CurvedPolygon<BCurve>;
  double abs_tol = 1e-10;

  // Quadrature nodes. Should be sufficiently high to pass tests
  int npts = 20;

  // Define anonymous functions for testing
  auto const_integrand = [](Point2D /*x*/) -> double { return 1.0; };
  auto poly_integrand = [](Point2D x) -> double { return x[0] * x[1] * x[1]; };
  auto transc_integrand = [](Point2D x) -> double { return std::sin(x[0] * x[1]); };

  // Test on triangular domain
  Point2D trinodes1[] = {Point2D {0.0, 0.0}, Point2D {1.0, 0.0}};
  BCurve tri1(trinodes1, 1);

  Point2D trinodes2[] = {Point2D {1.0, 0.0}, Point2D {0.0, 1.0}};
  BCurve tri2(trinodes2, 1);

  Point2D trinodes3[] = {Point2D {0.0, 1.0}, Point2D {0.0, 0.0}};
  BCurve tri3(trinodes3, 1);

  BCurve triangle_edges[] = {tri1, tri2, tri3};
  CPolygon triangle(triangle_edges, 3);

  // Compare against hand computed/high-precision calculated values
  EXPECT_NEAR(evaluate_area_integral(triangle, const_integrand, npts), 0.5, abs_tol);
  EXPECT_NEAR(evaluate_area_integral(triangle, poly_integrand, npts), 1.0 / 60.0, abs_tol);
  EXPECT_NEAR(evaluate_area_integral(triangle, transc_integrand, npts), 0.0415181074232, abs_tol);

  // Test on parabolic domain (between f(x) = 1-x^2 and g(x) = x^2-1, shifted to the right 1 unit)
  Point2D paranodes1[] = {Point2D {2.0, 0.0}, Point2D {1.0, 2.0}, Point2D {0.0, 0.0}};
  BCurve para1(paranodes1, 2);

  Point2D paranodes2[] = {Point2D {0.0, 0.0}, Point2D {1.0, -2.0}, Point2D {2.0, 0.0}};
  BCurve para2(paranodes2, 2);

  BCurve parabola_edges[] = {para1, para2};
  CPolygon parabola_shape(parabola_edges, 2);

  // Compare against hand computed/high-precision calculated values
  EXPECT_NEAR(evaluate_area_integral(parabola_shape, const_integrand, npts), 8.0 / 3.0, abs_tol);
  EXPECT_NEAR(evaluate_area_integral(parabola_shape, poly_integrand, npts), 64.0 / 105.0, abs_tol);
  EXPECT_NEAR(evaluate_area_integral(parabola_shape, transc_integrand, npts), 0.0, abs_tol);

  // Ensure compatibility with curved polygons
  BCurve pedges[2] = {para1, para2};
  CPolygon parabola_polygon(pedges, 2);
  EXPECT_NEAR(evaluate_area_integral(parabola_polygon, const_integrand, npts), 8.0 / 3.0, abs_tol);
  EXPECT_NEAR(evaluate_area_integral(parabola_polygon, poly_integrand, npts), 64.0 / 105.0, abs_tol);
  EXPECT_NEAR(evaluate_area_integral(parabola_polygon, transc_integrand, npts), 0.0, abs_tol);

  // Test on a unit square
  Point2D squarenodes1[] = {Point2D {0.0, 0.0}, Point2D {1.0, 0.0}};
  BCurve square1(squarenodes1, 1);

  Point2D squarenodes2[] = {Point2D {1.0, 0.0}, Point2D {1.0, 1.0}};
  BCurve square2(squarenodes2, 1);

  Point2D squarenodes3[] = {Point2D {1.0, 1.0}, Point2D {0.0, 1.0}};
  BCurve square3(squarenodes3, 1);

  Point2D squarenodes4[] = {Point2D {0.0, 1.0}, Point2D {0.0, 0.0}};
  BCurve square4(squarenodes4, 1);

  BCurve square_edges[] = {square1, square2, square3, square4};
  CPolygon square(square_edges, 4);

  EXPECT_NEAR(evaluate_area_integral(square, const_integrand, npts), 1.0, abs_tol);
}

TEST(primal_integral, evaluate_area_integral_aggregate)
{
  using Point2D = primal::Point<double, 2>;
  using BCurve = primal::BezierCurve<double, 2>;
  using CPolygon = primal::CurvedPolygon<BCurve>;
  using ReturnType = primal::Vector<double, 3>;
  double abs_tol = 1e-10;

  // Quadrature nodes. Should be sufficiently high to pass tests
  int npts = 20;

  // Define anonymous function for testing component-wise integration
  //  of an function with a vector return type.
  // This is equivalent to evaluating three separate integrands
  //  without unnecessarily repeating geometric processing.
  auto aggregate_integrand = [](Point2D x) -> ReturnType {
    return ReturnType {1.0, x[0] * x[1] * x[1], std::sin(x[0] * x[1])};
  };

  // Test on triangular domain
  Point2D trinodes1[] = {Point2D {0.0, 0.0}, Point2D {1.0, 0.0}};
  BCurve tri1(trinodes1, 1);

  Point2D trinodes2[] = {Point2D {1.0, 0.0}, Point2D {0.0, 1.0}};
  BCurve tri2(trinodes2, 1);

  Point2D trinodes3[] = {Point2D {0.0, 1.0}, Point2D {0.0, 0.0}};
  BCurve tri3(trinodes3, 1);

  BCurve triangle_edges[] = {tri1, tri2, tri3};
  CPolygon triangle(triangle_edges, 3);

  // Compare against hand computed/high-precision calculated values
  auto observed = evaluate_area_integral(triangle, aggregate_integrand, npts);
  auto expected = ReturnType {0.5, 1.0 / 60.0, 0.0415181074232};
  for(int i = 0; i < 3; ++i)
  {
    EXPECT_NEAR(observed[i], expected[i], abs_tol);
  }

  // Test on parabolic domain (between f(x) = 1-x^2 and g(x) = x^2-1, shifted to the right 1 unit)
  Point2D paranodes1[] = {Point2D {2.0, 0.0}, Point2D {1.0, 2.0}, Point2D {0.0, 0.0}};
  BCurve para1(paranodes1, 2);

  Point2D paranodes2[] = {Point2D {0.0, 0.0}, Point2D {1.0, -2.0}, Point2D {2.0, 0.0}};
  BCurve para2(paranodes2, 2);

  BCurve parabola_edges[] = {para1, para2};
  CPolygon parabola_shape(parabola_edges, 2);

  // Compare against hand computed/high-precision calculated values
  observed = evaluate_area_integral(parabola_shape, aggregate_integrand, npts);
  expected = ReturnType {8.0 / 3.0, 64.0 / 105.0, 0.0};
  for(int i = 0; i < 3; ++i)
  {
    EXPECT_NEAR(observed[i], expected[i], abs_tol);
  }

  // Ensure compatibility with curved polygons
  BCurve pedges[2] = {para1, para2};
  CPolygon parabola_polygon(pedges, 2);
  observed = evaluate_area_integral(parabola_polygon, aggregate_integrand, npts);
  for(int i = 0; i < 3; ++i)
  {
    EXPECT_NEAR(observed[i], expected[i], abs_tol);
  }
}

TEST(primal_integral, evaluate_line_integral_scalar)
{
  using Point2D = primal::Point<double, 2>;
  using BCurve = primal::BezierCurve<double, 2>;
  using CPolygon = primal::CurvedPolygon<BCurve>;
  double abs_tol = 1e-10;

  // Quadrature nodes. Should be sufficiently high to pass tests
  int npts = 30;

  // Define anonymous functions for testing
  auto const_integrand = [](Point2D /*x*/) -> double { return 1.0; };
  auto poly_integrand = [](Point2D x) -> double { return x[0] * x[1] * x[1]; };
  auto transc_integrand = [](Point2D x) -> double { return std::sin(x[0] * x[1]); };

  // Test on single parabolic segment
  Point2D paranodes[] = {Point2D {-1.0, 1.0}, Point2D {0.5, -2.0}, Point2D {2.0, 4.0}};
  BCurve parabola_segment(paranodes, 2);

  // Compare against hand computed/high-precision calculated values.

  // Constant integrand line integral is equivalent to arc-length calculation
  EXPECT_NEAR(evaluate_line_integral(parabola_segment, const_integrand, npts), 6.12572661998, abs_tol);

  EXPECT_NEAR(evaluate_line_integral(parabola_segment, poly_integrand, npts), 37.8010703669, abs_tol);
  EXPECT_NEAR(evaluate_line_integral(parabola_segment, transc_integrand, npts),
              0.495907795678,
              abs_tol);

  // Test on a collection of Bezier curves
  Point2D segnodes1[] = {Point2D {-1.0, -1.0},
                         Point2D {-1.0 / 3.0, 1.0},
                         Point2D {1.0 / 3.0, -1.0},
                         Point2D {1.0, 1.0}};
  BCurve cubic_segment(segnodes1, 3);

  Point2D segnodes2[] = {Point2D {1.0, 1.0}, Point2D {-1.0, 0.0}};
  BCurve linear_segment(segnodes2, 1);

  Point2D segnodes3[] = {Point2D {-1.0, 0.0}, Point2D {-3.0, 1.0}, Point2D {-1.0, 2.0}};
  BCurve quadratic_segment(segnodes3, 2);

  BCurve connected_curve_edges[] = {cubic_segment, linear_segment, quadratic_segment};
  CPolygon connected_curve(connected_curve_edges, 3);

  EXPECT_NEAR(evaluate_line_integral(connected_curve, const_integrand, npts), 8.28968500196, abs_tol);
  EXPECT_NEAR(evaluate_line_integral(connected_curve, poly_integrand, npts), -5.97565740064, abs_tol);
  EXPECT_NEAR(evaluate_line_integral(connected_curve, transc_integrand, npts),
              -0.574992518405,
              abs_tol);

  // Test algorithm on disconnected curves
  BCurve disconnected_curve_edges[] = {cubic_segment, quadratic_segment};
  CPolygon disconnected_curve(disconnected_curve_edges, 2);

  EXPECT_NEAR(evaluate_line_integral(disconnected_curve, const_integrand, npts),
              6.05361702446,
              abs_tol);
  EXPECT_NEAR(evaluate_line_integral(disconnected_curve, poly_integrand, npts),
              -6.34833539689,
              abs_tol);
  EXPECT_NEAR(evaluate_line_integral(disconnected_curve, transc_integrand, npts),
              -0.914161242161,
              abs_tol);
}

TEST(primal_integral, evaluate_line_integral_aggregate)
{
  using Point2D = primal::Point<double, 2>;
  using BCurve = primal::BezierCurve<double, 2>;
  using CPolygon = primal::CurvedPolygon<BCurve>;
  using ReturnType = primal::Vector<double, 3>;
  double abs_tol = 1e-10;

  // Quadrature nodes. Should be sufficiently high to pass tests
  int npts = 30;

  // Define anonymous function for testing component-wise integration
  //  of an function with a vector return type.
  // This is equivalent to evaluating three separate integrands
  //  without unnecessarily repeating geometric processing.
  auto aggregate_integrand = [](Point2D x) -> ReturnType {
    return ReturnType {1.0, x[0] * x[1] * x[1], std::sin(x[0] * x[1])};
  };

  // Test on single parabolic segment
  Point2D paranodes[] = {Point2D {-1.0, 1.0}, Point2D {0.5, -2.0}, Point2D {2.0, 4.0}};
  BCurve parabola_segment(paranodes, 2);

  // Compare against hand computed/high-precision calculated values.
  auto observed = evaluate_line_integral(parabola_segment, aggregate_integrand, npts);
  auto expected = ReturnType {6.12572661998, 37.8010703669, 0.495907795678};
  for(int i = 0; i < 3; ++i)
  {
    EXPECT_NEAR(observed[i], expected[i], abs_tol);
  }

  // Test on a collection of Bezier curves
  Point2D segnodes1[] = {Point2D {-1.0, -1.0},
                         Point2D {-1.0 / 3.0, 1.0},
                         Point2D {1.0 / 3.0, -1.0},
                         Point2D {1.0, 1.0}};
  BCurve cubic_segment(segnodes1, 3);

  Point2D segnodes2[] = {Point2D {1.0, 1.0}, Point2D {-1.0, 0.0}};
  BCurve linear_segment(segnodes2, 1);

  Point2D segnodes3[] = {Point2D {-1.0, 0.0}, Point2D {-3.0, 1.0}, Point2D {-1.0, 2.0}};
  BCurve quadratic_segment(segnodes3, 2);

  BCurve connected_curve_edges[] = {cubic_segment, linear_segment, quadratic_segment};
  CPolygon connected_curve(connected_curve_edges, 3);

  observed = evaluate_line_integral(connected_curve, aggregate_integrand, npts);
  expected = ReturnType {8.28968500196, -5.97565740064, -0.574992518405};
  for(int i = 0; i < 3; ++i)
  {
    EXPECT_NEAR(observed[i], expected[i], abs_tol);
  }

  // Test algorithm on disconnected curves
  BCurve disconnected_curve_edges[] = {cubic_segment, quadratic_segment};
  CPolygon disconnected_curve(disconnected_curve_edges, 2);

  observed = evaluate_line_integral(disconnected_curve, aggregate_integrand, npts);
  expected = ReturnType {6.05361702446, -6.34833539689, -0.914161242161};
  for(int i = 0; i < 3; ++i)
  {
    EXPECT_NEAR(observed[i], expected[i], abs_tol);
  }
}

TEST(primal_integral, evaluate_line_integral_vector)
{
  using Point2D = primal::Point<double, 2>;
  using Vector2D = primal::Vector<double, 2>;
  using BCurve = primal::BezierCurve<double, 2>;
  using CPolygon = primal::CurvedPolygon<BCurve>;
  double abs_tol = 1e-10;

  // Quadrature nodes. Should be sufficiently high to pass tests
  int npts = 30;

  // Test on a single line segment
  auto vec_field = [](Point2D x) -> Vector2D { return Vector2D({x[1] * x[1], 3 * x[0] - 6 * x[1]}); };

  Point2D segnodes[] = {Point2D {3.0, 7.0}, Point2D {0.0, 12.0}};
  BCurve linear_segment(segnodes, 1);

  // Compare against hand computed values
  EXPECT_NEAR(evaluate_vector_line_integral(linear_segment, vec_field, npts), -1079.0 / 2.0, abs_tol);

  // Test on a closed curve
  auto area_field = [](Point2D x) -> Vector2D { return Vector2D({-0.5 * x[1], 0.5 * x[0]}); };
  auto conservative_field = [](Point2D x) -> Vector2D {
    return Vector2D({2 * x[0] * x[1] * x[1], 2 * x[0] * x[0] * x[1]});
  };
  auto winding_field = [](Point2D x) -> Vector2D {
    double denom = 2 * M_PI * (x[0] * x[0] + x[1] * x[1]);
    return Vector2D({-x[1] / denom, x[0] / denom});
  };

  Point2D paranodes1[] = {Point2D {1.0, 0.0}, Point2D {0.0, 2.0}, Point2D {-1.0, 0.0}};
  BCurve para1(paranodes1, 2);

  Point2D paranodes2[] = {Point2D {-1.0, 0.0}, Point2D {0.0, -2.0}, Point2D {1.0, 0.0}};
  BCurve para2(paranodes2, 2);

  BCurve parabola_shape_edges[] = {para1, para2};
  CPolygon parabola_shape(parabola_shape_edges, 2);

  // This vector field calculates the area of the region
  EXPECT_NEAR(evaluate_vector_line_integral(parabola_shape, area_field, npts), 8.0 / 3.0, abs_tol);

  // This vector field is conservative, so it should evaluate to zero
  EXPECT_NEAR(evaluate_vector_line_integral(parabola_shape, conservative_field, npts), 0.0, abs_tol);

  // This vector field is generated by a in/out query, should return 1 (inside)
  EXPECT_NEAR(evaluate_vector_line_integral(parabola_shape, winding_field, npts), 1.0, abs_tol);

  // Test algorithm on disconnected curves
  Point2D paranodes2_shifted[] = {Point2D {-1.0, -1.0}, Point2D {0.0, -3.0}, Point2D {1.0, -1.0}};
  BCurve para2_shift(paranodes2_shifted, 2);

  BCurve disconnected_parabola_edges[] = {para1, para2_shift};
  CPolygon disconnected_parabola_shape(disconnected_parabola_edges, 2);

  EXPECT_NEAR(evaluate_vector_line_integral(disconnected_parabola_shape, area_field, npts),
              11.0 / 3.0,
              abs_tol);
  EXPECT_NEAR(evaluate_vector_line_integral(disconnected_parabola_shape, conservative_field, npts),
              0.0,
              abs_tol);
  EXPECT_NEAR(evaluate_vector_line_integral(disconnected_parabola_shape, winding_field, npts),
              0.75,
              abs_tol);
}

TEST(primal_integral, evaluate_integral_3D)
{
  using Point3D = primal::Point<double, 3>;
  using Vector3D = primal::Vector<double, 3>;
  using BCurve = primal::BezierCurve<double, 3>;
  double abs_tol = 1e-10;

  // Quadrature nodes. Should be sufficiently high to pass tests
  int npts = 30;

  const int order = 3;
  Point3D data[order + 1] = {Point3D {0.6, 1.2, 1.0},
                             Point3D {1.3, 1.6, 1.8},
                             Point3D {2.9, 2.4, 2.3},
                             Point3D {3.2, 3.5, 3.0}};
  BCurve spatial_arc(data, 3);

  auto const_integrand = [](Point3D /*x*/) -> double { return 1.0; };
  auto transc_integrand = [](Point3D x) -> double { return std::sin(x[0] * x[1] * x[2]); };

  auto vector_field = [](Point3D x) -> Vector3D {
    return Vector3D({4 * x[1] * x[1], 8 * x[0] * x[1], 1.0});
  };

  // Test line integral on scalar domain againt values computed with external software
  EXPECT_NEAR(evaluate_line_integral(spatial_arc, const_integrand, npts), 4.09193268998, abs_tol);
  EXPECT_NEAR(evaluate_line_integral(spatial_arc, transc_integrand, npts), 0.515093324547, abs_tol);

  // Test line integral on vector domain againt values computed with external software
  EXPECT_NEAR(evaluate_vector_line_integral(spatial_arc, vector_field, npts), 155.344, abs_tol);
}

TEST(primal_integral, evaluate_integral_rational)
{
  using Point2D = primal::Point<double, 2>;
  using Vector2D = primal::Vector<double, 2>;
  using BCurve = primal::BezierCurve<double, 2>;
  using CPolygon = primal::CurvedPolygon<BCurve>;
  double abs_tol = 1e-10;

  // Quadrature nodes. Should be sufficiently high to pass tests
  int npts = 20;

  // Elliptical arc shape
  Point2D ellipse_nodes[] = {Point2D {2.0, 0.0}, Point2D {2.0, 1.0}, Point2D {0.0, 1.0}};
  double weights[] = {2.0, 1.0, 1.0};
  BCurve ellipse_arc(ellipse_nodes, weights, 2);

  Point2D leg1_nodes[] = {Point2D {0.0, 1.0}, {0.0, 0.0}};
  BCurve leg1(leg1_nodes, 1);

  Point2D leg2_nodes[] = {Point2D {0.0, 0.0}, {2.0, 0.0}};
  BCurve leg2(leg2_nodes, 1);

  CPolygon quarter_ellipse;
  quarter_ellipse.addEdge(ellipse_arc);
  quarter_ellipse.addEdge(leg1);
  quarter_ellipse.addEdge(leg2);

  auto const_integrand = [](Point2D /*x*/) -> double { return 1.0; };
  auto transc_integrand = [](Point2D x) -> double { return std::sin(x[0] * x[1]); };

  auto area_field = [](Point2D x) -> Vector2D { return Vector2D({-0.5 * x[1], 0.5 * x[0]}); };
  auto conservative_field = [](Point2D x) -> Vector2D {
    return Vector2D({2 * x[0] * x[1] * x[1], 2 * x[0] * x[0] * x[1]});
  };

  // Test area integrals with scalar integrand againt values computed with external software
  EXPECT_NEAR(evaluate_area_integral(quarter_ellipse, const_integrand, npts),
              M_PI * 2 * 1 / 4.0,
              abs_tol);
  EXPECT_NEAR(evaluate_area_integral(quarter_ellipse, transc_integrand, npts), 0.472951736306, abs_tol);

  // Test line integral on scalar domain againt values computed with external software
  EXPECT_NEAR(evaluate_line_integral(ellipse_arc, const_integrand, npts), 2.42211205514, abs_tol);
  EXPECT_NEAR(evaluate_line_integral(ellipse_arc, transc_integrand, npts), 1.38837959326, abs_tol);

  // Test line integral on vector domain againt values computed with external software
  EXPECT_NEAR(evaluate_vector_line_integral(ellipse_arc, area_field, npts),
              M_PI * 2 * 1 / 4.0,
              abs_tol);
  EXPECT_NEAR(evaluate_vector_line_integral(quarter_ellipse, conservative_field, npts), 0, abs_tol);
}

TEST(primal_integral, evaluate_integral_nurbs)
{
  using Point2D = primal::Point<double, 2>;
  using Vector2D = primal::Vector<double, 2>;
  using NCurve = primal::NURBSCurve<double, 2>;
  using CPolygon = primal::CurvedPolygon<NCurve>;
  double abs_tol = 1e-10;

  // Quadrature nodes. Should be sufficiently high to pass tests
  int npts = 20;

  // Test integrals with same integrand and curves as `evaluate_integral_rational`,
  //  but insert some knots to make the Bezier extraction more interesting

  // Elliptical arc shape
  Point2D ellipse_nodes[] = {Point2D {2.0, 0.0}, Point2D {2.0, 1.0}, Point2D {0.0, 1.0}};
  double weights[] = {2.0, 1.0, 1.0};
  NCurve ellipse_arc(ellipse_nodes, weights, 3, 2);
  ellipse_arc.insertKnot(0.3, 2);
  ellipse_arc.insertKnot(0.7, 1);

  Point2D leg1_nodes[] = {Point2D {0.0, 1.0}, {0.0, 0.0}};
  NCurve leg1(leg1_nodes, 2, 1);
  leg1.insertKnot(0.4, 1);

  Point2D leg2_nodes[] = {Point2D {0.0, 0.0}, {2.0, 0.0}};
  NCurve leg2(leg2_nodes, 2, 1);
  leg2.insertKnot(0.6, 1);

  CPolygon quarter_ellipse;
  quarter_ellipse.addEdge(ellipse_arc);
  quarter_ellipse.addEdge(leg1);
  quarter_ellipse.addEdge(leg2);

  auto const_integrand = [](Point2D /*x*/) -> double { return 1.0; };
  auto transc_integrand = [](Point2D x) -> double { return std::sin(x[0] * x[1]); };

  auto area_field = [](Point2D x) -> Vector2D { return Vector2D({-0.5 * x[1], 0.5 * x[0]}); };
  auto conservative_field = [](Point2D x) -> Vector2D {
    return Vector2D({2 * x[0] * x[1] * x[1], 2 * x[0] * x[0] * x[1]});
  };

  EXPECT_NEAR(evaluate_area_integral(quarter_ellipse, const_integrand, npts),
              M_PI * 2 * 1 / 4.0,
              abs_tol);
  EXPECT_NEAR(evaluate_area_integral(quarter_ellipse, transc_integrand, npts), 0.472951736306, abs_tol);

  EXPECT_NEAR(evaluate_line_integral(ellipse_arc, const_integrand, npts), 2.42211205514, abs_tol);
  EXPECT_NEAR(evaluate_line_integral(ellipse_arc, transc_integrand, npts), 1.38837959326, abs_tol);

  EXPECT_NEAR(evaluate_vector_line_integral(ellipse_arc, area_field, npts),
              M_PI * 2 * 1 / 4.0,
              abs_tol);
  EXPECT_NEAR(evaluate_vector_line_integral(quarter_ellipse, conservative_field, npts), 0, abs_tol);
}

TEST(primal_integral, evaluate_nurbs_surface_normal)
{
  const int DIM = 3;
  using Point2D = primal::Point<double, 2>;
  using Point3D = primal::Point<double, 3>;
  using NURBSPatchType = primal::NURBSPatch<double, DIM>;

  const int npts_u = 5;
  const int npts_v = 4;

  const int degree_u = 3;
  const int degree_v = 2;

  // clang-format off
	Point3D controlPoints[5 * 4] = {
		Point3D {0, 0, 0}, Point3D {0, 4,  0}, Point3D {0, 8, -3}, Point3D {0, 12, 0},
		Point3D {2, 0, 6}, Point3D {2, 4,  0}, Point3D {2, 8,  0}, Point3D {2, 12, 0},
		Point3D {4, 0, 0}, Point3D {4, 4,  0}, Point3D {4, 8,  3}, Point3D {4, 12, 0},
		Point3D {6, 0, 0}, Point3D {6, 4, -3}, Point3D {6, 8,  0}, Point3D {6, 12, 0},
		Point3D {8, 0, 0}, Point3D {8, 4,  0}, Point3D {8, 8,  0}, Point3D {8, 12, 0} };

	double weights[5 * 4] = {
		1.0, 2.0, 3.0, 2.0,
		2.0, 3.0, 4.0, 3.0,
		3.0, 4.0, 5.0, 4.0,
		4.0, 5.0, 6.0, 5.0,
		5.0, 6.0, 7.0, 6.0 };
  // clang-format on

  NURBSPatchType nPatch(controlPoints, weights, npts_u, npts_v, degree_u, degree_v);
  nPatch.makeTriviallyTrimmed();

  int npts = 20;
  double abs_tol = 1e-10;

  // Calculate the average surface normal using boundary formulation
  auto ueda_formula = nPatch.calculateUntrimmedPatchNormal(npts);

  // Calculate the average surface normal component by component,
  //  i.e. \int_S \hat{n}(x) dA, where \hat{n} is the unit normal at surface point x

  // Because the NURBS surface has discontinuous derivatives at u, v = 0.5,
  //  this integrand does too, and we need to split to get accurate results
  //  using Gaussian quadrature

  NURBSPatchType split_patch[4];
  nPatch.split(0.5, 0.5, split_patch[0], split_patch[1], split_patch[2], split_patch[3]);

  for(int N = 0; N < 3; ++N)
  {
    auto avg_surface_normal_integrand = [&nPatch, &N](Point2D x) -> double {
      return nPatch.normal(x[0], x[1])[N];
    };

    // Integrate over each patch
    double direct_formula = 0.0;
    for(int i = 0; i < 4; ++i)
    {
      direct_formula +=
        evaluate_area_integral(split_patch[i].getTrimmingCurves(), avg_surface_normal_integrand, npts);
    }

    EXPECT_NEAR(ueda_formula[N], direct_formula, abs_tol);
  }
}

TEST(primal_integral, evaluate_integral_nurbs_gwn_cache)
{
  using Point2D = primal::Point<double, 2>;
  using Vector2D = primal::Vector<double, 2>;
  using NCurve = primal::NURBSCurve<double, 2>;
  using NCache = primal::NURBSCurve<double, 2>;
  using CPolygon = primal::CurvedPolygon<NCache>;
  double abs_tol = 1e-10;

  // Quadrature nodes. Should be sufficiently high to pass tests
  int npts = 20;

  // Test integrals with same integrand and curves as `evaluate_integral_rational`,
  //  but with curves added to detail::NURBSCurveGWNCache objects to ensure template compatibility,
  //  even if there isn't a compelling reason to use GWN caches for this purpose

  // Elliptical arc shape
  Point2D ellipse_nodes[] = {Point2D {2.0, 0.0}, Point2D {2.0, 1.0}, Point2D {0.0, 1.0}};
  double weights[] = {2.0, 1.0, 1.0};
  NCurve ellipse_arc(ellipse_nodes, weights, 3, 2);
  ellipse_arc.insertKnot(0.3, 2);
  ellipse_arc.insertKnot(0.7, 1);

  Point2D leg1_nodes[] = {Point2D {0.0, 1.0}, {0.0, 0.0}};
  NCurve leg1(leg1_nodes, 2, 1);
  leg1.insertKnot(0.4, 1);

  Point2D leg2_nodes[] = {Point2D {0.0, 0.0}, {2.0, 0.0}};
  NCurve leg2(leg2_nodes, 2, 1);
  leg2.insertKnot(0.6, 1);

  CPolygon quarter_ellipse;
  quarter_ellipse.addEdge(NCache(ellipse_arc));
  quarter_ellipse.addEdge(NCache(leg1));
  quarter_ellipse.addEdge(NCache(leg2));

  auto const_integrand = [](Point2D /*x*/) -> double { return 1.0; };
  auto transc_integrand = [](Point2D x) -> double { return std::sin(x[0] * x[1]); };

  auto area_field = [](Point2D x) -> Vector2D { return Vector2D({-0.5 * x[1], 0.5 * x[0]}); };
  auto conservative_field = [](Point2D x) -> Vector2D {
    return Vector2D({2 * x[0] * x[1] * x[1], 2 * x[0] * x[0] * x[1]});
  };

  EXPECT_NEAR(evaluate_area_integral(quarter_ellipse, const_integrand, npts),
              M_PI * 2 * 1 / 4.0,
              abs_tol);
  EXPECT_NEAR(evaluate_area_integral(quarter_ellipse, transc_integrand, npts), 0.472951736306, abs_tol);

  EXPECT_NEAR(evaluate_line_integral(ellipse_arc, const_integrand, npts), 2.42211205514, abs_tol);
  EXPECT_NEAR(evaluate_line_integral(ellipse_arc, transc_integrand, npts), 1.38837959326, abs_tol);

  EXPECT_NEAR(evaluate_vector_line_integral(ellipse_arc, area_field, npts),
              M_PI * 2 * 1 / 4.0,
              abs_tol);
  EXPECT_NEAR(evaluate_vector_line_integral(quarter_ellipse, conservative_field, npts), 0, abs_tol);
}

#ifdef AXOM_USE_MFEM
TEST(primal_integral, check_axom_mfem_quadrature_values)
{
  const int N = 200;

  for(int npts = 1; npts <= N; ++npts)
  {
    // Generate the Axom quadrature rule
    axom::numerics::QuadratureRule axom_rule = axom::numerics::get_gauss_legendre(npts);

    // Generate the MFEM quadrature rule
    static mfem::IntegrationRules my_IntRules(0, mfem::Quadrature1D::GaussLegendre);
    const mfem::IntegrationRule& mfem_rule = my_IntRules.Get(mfem::Geometry::SEGMENT, 2 * npts - 1);

    // Check that the nodes and weights are the same between the two rules
    for(int j = 0; j < npts; ++j)
    {
      EXPECT_NEAR(axom_rule.node(j), mfem_rule.IntPoint(j).x, axom::numeric_limits<double>::epsilon());
      EXPECT_NEAR(axom_rule.weight(j),
                  mfem_rule.IntPoint(j).weight,
                  axom::numeric_limits<double>::epsilon());
    }
  }
}
#endif

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  int result = RUN_ALL_TESTS();

  return result;
}
