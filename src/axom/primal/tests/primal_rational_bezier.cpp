// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! 
 * \file primal_rational_bezier.cpp
 * \brief This file tests primal's BezierCurve functionality
 * with respect to rational Bezier curves
 */

#include "axom/config.hpp"

#include "gtest/gtest.h"

#include "axom/primal.hpp"
#include "axom/slic.hpp"
#include "axom/fmt.hpp"

namespace primal = axom::primal;

//------------------------------------------------------------------------------
TEST(primal_rationalbezier, constructor)
{
  const int DIM = 3;
  using CoordType = double;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  {
    SLIC_INFO("Testing default BezierCurve constructor ");
    BezierCurveType bCurve;
    EXPECT_FALSE(bCurve.isRational());
  }

  {
    SLIC_INFO("Testing BezierCurve order constructor ");

    BezierCurveType bCurve(1);
    EXPECT_FALSE(bCurve.isRational());
  }
}

//----------------------------------------------------------------------------------
TEST(primal_rationalbezier, point_array_constructor)
{
  SLIC_INFO("Testing point array constructor");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  PointType controlPoints[2] = {PointType {0.6, 1.2, 1.0}, PointType {0.0, 1.6, 1.8}};
  double weights[2] = {0.5, 1.5};

  BezierCurveType bCurve(controlPoints, weights, 1);

  EXPECT_EQ(1, bCurve.getOrder());
  for(int p = 0; p <= bCurve.getOrder(); ++p)
  {
    auto& pt = bCurve[p];
    auto& wt = bCurve.getWeight(p);
    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_DOUBLE_EQ(controlPoints[p][i], pt[i]);
    }
    EXPECT_DOUBLE_EQ(weights[p], wt);
  }
}

//------------------------------------------------------------------------------
TEST(primal_rationalbezier, evaluate)
{
  SLIC_INFO("Testing Rational Bezier evaluation");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  const int max_order = 3;
  PointType data[max_order + 1] = {PointType {0.6, 1.2, 1.0},
                                   PointType {1.3, 1.6, 1.8},
                                   PointType {2.9, 2.4, 2.3},
                                   PointType {3.2, 3.5, 3.0}};
  double weights[max_order + 1] = {1, 2, 3, 4};

  // clang-format off
  PointType exp_vals[4][3] = {{PointType {0.6, 1.2, 1.0}, PointType {    0.6,     1.2,     1.0}, PointType {0.6, 1.2, 1.0}},
                              {PointType {0.6, 1.2, 1.0}, PointType {16./15., 22./15., 23./15.}, PointType {1.3, 1.6, 1.8}},
                              {PointType {0.6, 1.2, 1.0}, PointType { 1.8125,    1.85,  1.8875}, PointType {2.9, 2.4, 2.3}},
                              {PointType {0.6, 1.2, 1.0}, PointType {  2.365,    2.32,   2.225}, PointType {3.2, 3.5, 3.0}}};
  // clang-format on

  // Test evaluation at various orders, using the first `ord`
  //  points in `data` as control points.
  for(int ord = 0; ord <= max_order; ++ord)
  {
    BezierCurveType curve(data, weights, ord);

    PointType calc_start = curve.evaluate(0.0);
    PointType calc_mid = curve.evaluate(0.5);
    PointType calc_end = curve.evaluate(1.0);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(calc_start[i], exp_vals[ord][0][i], 1e-15);
      EXPECT_NEAR(calc_mid[i], exp_vals[ord][1][i], 1e-15);
      EXPECT_NEAR(calc_end[i], exp_vals[ord][2][i], 1e-15);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_rationalbezier, first_derivative)
{
  SLIC_INFO("Testing Rational Bezier derivative calculation");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  const int max_order = 3;
  PointType data[max_order + 1] = {PointType {0.6, 1.2, 1.0},
                                   PointType {1.3, 1.6, 1.8},
                                   PointType {2.9, 2.4, 2.3},
                                   PointType {3.2, 3.5, 3.0}};
  double weights[max_order + 1] = {1, 2, 3, 4};

  // clang-format off
  VectorType exp_vals[4][3] = {{VectorType {0.0, 0.0, 0.0}, VectorType {    0.0,     0.0,     0.0}, VectorType {    0.0,     0.0,   0.0}},
                               {VectorType {1.4, 0.8, 1.6}, VectorType {28./45., 16./45., 32./45.}, VectorType {   0.35,     0.2,   0.4}},
                               {VectorType {2.8, 1.6, 3.2}, VectorType { 2.2375,    1.15,  1.0625}, VectorType {32./15., 16./15., 2./3.}},
                               {VectorType {4.2, 2.4, 4.8}, VectorType {  2.652,   2.256,    1.62}, VectorType {  0.675,   2.475, 1.575}}};
  // clang-format on

  // Test derivative calculation at various orders, using the first `ord`
  //  points in `data` as control points.
  for(int ord = 0; ord <= max_order; ++ord)
  {
    BezierCurveType curve(data, weights, ord);

    VectorType calc_start = curve.dt(0.0);
    VectorType calc_mid = curve.dt(0.5);
    VectorType calc_end = curve.dt(1.0);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(calc_start[i], exp_vals[ord][0][i], 1e-15);
      EXPECT_NEAR(calc_mid[i], exp_vals[ord][1][i], 1e-15);
      EXPECT_NEAR(calc_end[i], exp_vals[ord][2][i], 1e-15);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_rationalbezier, second_derivative)
{
  SLIC_INFO("Testing Rational Bezier second derivative calculation");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  const int max_order = 3;
  PointType data[max_order + 1] = {PointType {0.6, 1.2, 1.0},
                                   PointType {1.3, 1.6, 1.8},
                                   PointType {2.9, 2.4, 2.3},
                                   PointType {3.2, 3.5, 3.0}};
  double weights[max_order + 1] = {1, 2, 3, 4};

  // clang-format off
  VectorType exp_vals[4][3] = {{VectorType { 0.0,  0.0,   0.0}, VectorType {       0.0,       0.0,        0.0}, VectorType {    0.0,     0.0,      0.0}},
                               {VectorType {-2.8, -1.6,  -3.2}, VectorType {-112./135., -64./135., -128./135.}, VectorType {  -0.35,    -0.2,     -0.4}},
                               {VectorType {-3.0, -2.4, -11.4}, VectorType {    -0.375,      -0.3,     -1.425}, VectorType { -1./9., -8./90., -19./45.}},
                               {VectorType {-0.6, -2.4, -24.6}, VectorType {   -3.8448,    0.3456,     -0.888}, VectorType {-4.0125,  0.4875,   0.3375}}};
  // clang-format on

  // Test derivative calculation at various orders, using the first `ord`
  //  points in `data` as control points.
  for(int ord = 0; ord <= max_order; ++ord)
  {
    BezierCurveType curve(data, weights, ord);

    VectorType calc_start = curve.dtdt(0.0);
    VectorType calc_mid = curve.dtdt(0.5);
    VectorType calc_end = curve.dtdt(1.0);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(calc_start[i], exp_vals[ord][0][i], 1e-14);
      EXPECT_NEAR(calc_mid[i], exp_vals[ord][1][i], 1e-14);
      EXPECT_NEAR(calc_end[i], exp_vals[ord][2][i], 1e-14);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_rationalbezier, split_cubic)
{
  SLIC_INFO("Testing Bezier splitting of a cubic");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  const int order = 3;
  PointType data[order + 1] = {PointType {0.6, 1.2, 1.0},
                               PointType {1.3, 1.6, 1.8},
                               PointType {2.9, 2.4, 2.3},
                               PointType {3.2, 3.5, 3.0}};

  double weights[order + 1] = {1, 2, 3, 4};

  BezierCurveType b2Curve(data, weights, order);

  BezierCurveType b3Curve(order);  // Checks split with order constructor
  BezierCurveType b4Curve;         // Checks split with default constructor
  b2Curve.split(.5, b3Curve, b4Curve);

  PointType b3Nodes[order + 1] = {PointType {0.6, 1.2, 1.0},
                                  PointType {16.0 / 15.0, 22.0 / 15.0, 23.0 / 15.0},
                                  PointType {1.8125, 1.85, 1.8875},
                                  PointType {2.365, 2.32, 2.225}};
  double b3Weights[order + 1] = {1.0, 1.5, 2.0, 2.5};

  PointType b4Nodes[order + 1] = {PointType {2.365, 2.32, 2.225},
                                  PointType {41.0 / 15.0, 79 / 30.0, 2.45},
                                  PointType {43.0 / 14.0, 106.0 / 35.0, 2.7},
                                  PointType {3.2, 3.5, 3.0}};
  double b4Weights[order + 1] = {2.5, 3.0, 3.5, 4.0};

  BezierCurveType b3True(b3Nodes, b3Weights, 3);
  BezierCurveType b4True(b4Nodes, b4Weights, 3);

  for(int i = 0; i < DIM; ++i)
  {
    for(int p = 0; p <= order; ++p)
    {
      EXPECT_NEAR(b3Curve[p][i], b3True[p][i], 1e-15);
      EXPECT_NEAR(b4Curve[p][i], b4True[p][i], 1e-15);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_rationalbezier, split_point)
{
  SLIC_INFO("Testing rational Bezier splitting for order 0");

  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  // Test order-0
  {
    const int order = 0;
    PointType data[order + 1] = {PointType::make_point(0.6, 1.2)};
    double weights[order + 1] = {1.5};

    BezierCurveType b(data, weights, order);

    BezierCurveType c1, c2;
    b.split(0.5, c1, c2);

    EXPECT_NEAR(weights[0], c1.getWeight(0), 1e-15);
    EXPECT_NEAR(weights[0], c2.getWeight(0), 1e-15);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(data[0][i], c1[0][i], 1e-15);
      EXPECT_NEAR(data[0][i], c2[0][i], 1e-15);
    }
  }
}

TEST(primal_rationalbezier, split_linear)
{
  SLIC_INFO("Testing rational Bezier splitting for order 1");

  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  const int order = 1;
  PointType data[order + 1] = {PointType {-1, -5}, PointType {1, 5}};
  double weights[order + 1] = {2, 8};

  BezierCurveType b(data, weights, order);

  {
    BezierCurveType c1, c2;
    b.split(0.5, c1, c2);

    EXPECT_NEAR(-1., c1[0][0], 1e-12);
    EXPECT_NEAR(-5., c1[0][1], 1e-12);
    EXPECT_NEAR(0.6, c1[1][0], 1e-12);
    EXPECT_NEAR(3.0, c1[1][1], 1e-12);

    EXPECT_NEAR(0.6, c2[0][0], 1e-12);
    EXPECT_NEAR(3.0, c2[0][1], 1e-12);
    EXPECT_NEAR(1.0, c2[1][0], 1e-12);
    EXPECT_NEAR(5.0, c2[1][1], 1e-12);
  }

  {
    BezierCurveType c1, c2;
    b.split(0.2, c1, c2);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(-1., c1[0][0], 1e-12);
      EXPECT_NEAR(-5., c1[0][1], 1e-12);
      EXPECT_NEAR(0.0, c1[1][0], 1e-12);
      EXPECT_NEAR(0.0, c1[1][1], 1e-12);

      EXPECT_NEAR(0.0, c2[0][0], 1e-12);
      EXPECT_NEAR(0.0, c2[0][1], 1e-12);
      EXPECT_NEAR(1.0, c2[1][0], 1e-12);
      EXPECT_NEAR(5.0, c2[1][1], 1e-12);
    }
  }
}

TEST(primal_beziercurve, isLinear)
{
  SLIC_INFO("Testing isLinear() on Rational Bezier curves");

  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using SegmentType = primal::Segment<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  // As calculated, linearity should be independent of weights
  {
    const int order = 2;
    auto curve = BezierCurveType(order);
    curve.makeRational();
    EXPECT_TRUE(curve.isLinear());

    // straight line
    curve[0] = PointType {1, 1};
    curve[1] = PointType {2, 2};
    curve[2] = PointType {3, 3};
    EXPECT_TRUE(curve.isLinear());

    // increase middle weight and check linearity with different tolerances
    VectorType v(curve[2], curve[0]);
    auto normal = VectorType {-v[1], v[0]};
    curve[1].array() += 0.005 * normal.array();
    SLIC_INFO("Updated curve: " << curve);

    SegmentType s(curve[2], curve[0]);
    SLIC_INFO("Sq dist: " << primal::squared_distance(curve[1], s));

    for(int i = 1; i < 128; i = i * 2)
    {
      curve.setWeight(1, i);

      // linear for a coarse tolerance
      {
        const double tol = 0.1;
        const double tol_sq = tol * tol;
        EXPECT_TRUE(curve.isLinear(tol_sq));
      }

      // non-linear for a finer tolerance
      {
        const double tol = 0.01;
        const double tol_sq = tol * tol;
        EXPECT_FALSE(curve.isLinear(tol_sq));
      }
    }
  }
}

TEST(primal_beziercurve, reverseOrientation)
{
  SLIC_INFO("Testing reverseOrientation() on Rational Bezier curves");

  {
    const int DIM = 2;
    using CoordType = double;
    using PointType = primal::Point<CoordType, DIM>;
    using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

    // test different orders
    for(int order = 0; order <= 10; ++order)
    {
      // control points for curve monotonically increase
      axom::Array<PointType> pts(order + 1);
      axom::Array<double> weights(order + 1);
      for(int i = 0; i <= order; ++i)
      {
        pts[i] = PointType(i);
        weights[i] = i + 1;
      }
      BezierCurveType curve(pts.data(), weights.data(), order);

      for(int i = 1; i <= order; ++i)
      {
        EXPECT_GT(curve[i][0], curve[i - 1][0]);
      }

      // create a reversed curve and check that it monotonically decreases
      BezierCurveType reversed = curve;
      reversed.reverseOrientation();

      for(int i = 1; i <= order; ++i)
      {
        EXPECT_LT(reversed[i][0], reversed[i - 1][0]);
      }

      // Check that the control points are actually reversed
      for(int i = 0; i <= order; ++i)
      {
        EXPECT_EQ(curve[i], reversed[order - i]);
      }

      // check that reversing again reverts to the original
      BezierCurveType reversedAgain = reversed;
      reversedAgain.reverseOrientation();
      EXPECT_EQ(curve, reversedAgain);
    }
  }
}

TEST(primal_beziercurve, isRational)
{
  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  const double abs_tol = 1e-8;

  const int order = 3;
  PointType data[order + 1] = {PointType {0.6, 1.2, 1.0},
                               PointType {1.3, 1.6, 1.8},
                               PointType {2.9, 2.4, 2.3},
                               PointType {3.2, 3.5, 3.0}};

  BezierCurveType bCurve(data, order);

  BezierCurveType rCurve(data, order);

  EXPECT_FALSE(rCurve.isRational());
  rCurve.makeRational();
  EXPECT_TRUE(rCurve.isRational());

  // Verify that makeRational makes curve trivially rational
  for(int i = 1; i <= order; i++)
  {
    EXPECT_DOUBLE_EQ(rCurve.getWeight(0), rCurve.getWeight(i));
  }

  // Verify that trivially rational Bezier curve is identical
  //  to polynomial Bezier curve
  for(double t = 0; t <= 1; t += 0.5)
  {
    EXPECT_NEAR(bCurve.evaluate(t)[0], rCurve.evaluate(t)[0], abs_tol);
    EXPECT_NEAR(bCurve.evaluate(t)[1], rCurve.evaluate(t)[1], abs_tol);
  }

  rCurve.makeNonrational();
  EXPECT_FALSE(rCurve.isRational());
}

TEST(primal_rationalbezier, rational_intersection)
{
  using Point2D = primal::Point<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;
  constexpr double abs_tol = 1e-8;

  // Intersecting of rational, circular arc shapes
  Point2D bot_nodes[] = {Point2D {1.0, 0.0}, Point2D {1.0, 1.0}, Point2D {0.0, 1.0}};

  Point2D top_nodes[] = {Point2D {1.3, 0.3}, Point2D {0.3, 0.3}, Point2D {0.3, 1.3}};

  double weights[] = {2.0, 1.0, 1.0};

  Bezier bottom_arc(bot_nodes, weights, 2);
  Bezier top_arc(top_nodes, weights, 2);

  axom::Array<double> sp, tp;
  EXPECT_TRUE(intersect(bottom_arc, top_arc, sp, tp));

  EXPECT_NEAR(bottom_arc.evaluate(sp[0])[0], top_arc.evaluate(tp[0])[0], abs_tol);
  EXPECT_NEAR(bottom_arc.evaluate(sp[0])[1], top_arc.evaluate(tp[0])[1], abs_tol);

  EXPECT_NEAR(bottom_arc.evaluate(sp[1])[0], top_arc.evaluate(tp[1])[0], abs_tol);
  EXPECT_NEAR(bottom_arc.evaluate(sp[1])[1], top_arc.evaluate(tp[1])[1], abs_tol);

  // Check intersection of nonrational and rational curve
  Point2D line_nodes[] = {Point2D {1.3, 0}, Point2D {0, 1.3}};
  Bezier line(line_nodes, 1);

  sp.clear();
  tp.clear();
  EXPECT_TRUE(intersect(bottom_arc, line, sp, tp));

  EXPECT_NEAR(bottom_arc.evaluate(sp[0])[0], line.evaluate(tp[0])[0], abs_tol);
  EXPECT_NEAR(bottom_arc.evaluate(sp[0])[1], line.evaluate(tp[0])[1], abs_tol);

  EXPECT_NEAR(bottom_arc.evaluate(sp[1])[0], line.evaluate(tp[1])[0], abs_tol);
  EXPECT_NEAR(bottom_arc.evaluate(sp[1])[1], line.evaluate(tp[1])[1], abs_tol);
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
