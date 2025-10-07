// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! 
 * \file primal_bezier_curve.cpp
 * \brief This file tests primal's Bezier curve functionality
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/operators/squared_distance.hpp"

namespace primal = axom::primal;

//------------------------------------------------------------------------------
TEST(primal_beziercurve, sizing_constructors)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;
  using CoordsVec = BezierCurveType::CoordsVec;

  // Testing default BezierCurve constructor w/ ord -1
  {
    BezierCurveType bCurve;

    constexpr int expOrder = -1;
    EXPECT_EQ(expOrder, bCurve.getOrder());
    EXPECT_EQ(expOrder + 1, bCurve.getControlPoints().size());
    EXPECT_EQ(expOrder + 1, bCurve.getNumControlPoints());
    EXPECT_EQ(CoordsVec(), bCurve.getControlPoints());
    EXPECT_FALSE(bCurve.isRational());
  }

  // Testing BezierCurve constructor w/ different orders
  for(int ord = -1; ord < 4; ++ord)
  {
    BezierCurveType bCurve(ord);
    EXPECT_EQ(ord, bCurve.getOrder());
    EXPECT_EQ(ord + 1, bCurve.getNumControlPoints());
    EXPECT_EQ(ord + 1, static_cast<int>(bCurve.getControlPoints().size()));
    EXPECT_FALSE(bCurve.isRational());
  }
}

//----------------------------------------------------------------------------------
TEST(primal_beziercurve, point_array_constructors_polynomial)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;
  using CoordsVec = typename BezierCurveType::CoordsVec;

  constexpr int order = 1;
  PointType controlPoints[2] = {PointType {0.6, 1.2, 1.0}, PointType {0.0, 1.6, 1.8}};
  CoordsVec pt_array {controlPoints[0], controlPoints[1]};

  auto checkCurve = [&](const BezierCurveType& bCurve) {
    EXPECT_EQ(order, bCurve.getOrder());
    for(int p = 0; p <= bCurve.getOrder(); ++p)
    {
      EXPECT_EQ(controlPoints[p], bCurve[p]);
    }
    EXPECT_FALSE(bCurve.isRational());
  };

  // test constructor from C-array
  {
    BezierCurveType bCurve(controlPoints, order);
    checkCurve(bCurve);
  }

  // test constructor from ArrayView generated from C-array
  {
    BezierCurveType bCurve(axom::ArrayView<PointType>(controlPoints, order + 1), order);
    checkCurve(bCurve);
  }

  // test constructor from axom Array
  {
    BezierCurveType bCurve(pt_array, order);
    checkCurve(bCurve);
  }

  // test constructor from ArrayView generated from axom Array
  {
    BezierCurveType bCurve(pt_array.view(), order);
    checkCurve(bCurve);
  }
}

//----------------------------------------------------------------------------------
TEST(primal_beziercurve, point_array_constructors_rational)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;
  using CoordsVec = typename BezierCurveType::CoordsVec;
  using WeightsVec = typename BezierCurveType::WeightsVec;

  constexpr int order = 2;
  PointType controlPoints[3] = {PointType {0.6, 1.2, 1.0},
                                PointType {0.0, 1.6, 1.8},
                                PointType {2.0, 2.6, 2.8}};
  double weights[3] = {1., .5, 1.};

  CoordsVec pts {controlPoints[0], controlPoints[1], controlPoints[2]};
  WeightsVec wts {weights[0], weights[1], weights[2]};

  auto checkCurve = [&](const BezierCurveType& bCurve) {
    EXPECT_EQ(order, bCurve.getOrder());
    for(int p = 0; p <= bCurve.getOrder(); ++p)
    {
      EXPECT_EQ(controlPoints[p], bCurve[p]);
      EXPECT_EQ(weights[p], bCurve.getWeight(p));
    }
    EXPECT_TRUE(bCurve.isRational());
  };

  // test constructor from C-array
  {
    BezierCurveType bCurve(controlPoints, weights, order);
    checkCurve(bCurve);
  }

  // test constructor from ArrayViews over the C-arrays
  {
    BezierCurveType bCurve(axom::ArrayView<PointType>(controlPoints, order + 1),
                           axom::ArrayView<double>(weights, order + 1),
                           order);
    checkCurve(bCurve);
  }

  // test constructor from axom Arrays
  {
    BezierCurveType bCurve(pts, wts, order);
    checkCurve(bCurve);
  }

  // test constructor from ArrayViews over the axom arrays
  {
    BezierCurveType bCurve(pts.view(), wts.view(), order);
    checkCurve(bCurve);
  }
}

//----------------------------------------------------------------------------------
TEST(primal_beziercurve, set_order)
{
  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO("Test adding control points to empty Bezier curve");

  BezierCurveType bCurve;
  EXPECT_EQ(-1, bCurve.getOrder());
  EXPECT_EQ(0, bCurve.getNumControlPoints());

  const int order = 1;
  PointType controlPoints[2] = {PointType {0.6, 1.2, 1.0}, PointType {0.0, 1.6, 1.8}};

  bCurve.setOrder(order);
  EXPECT_EQ(order, bCurve.getOrder());
  EXPECT_EQ(order + 1, bCurve.getNumControlPoints());

  bCurve[0] = controlPoints[0];
  bCurve[1] = controlPoints[1];

  for(int p = 0; p <= bCurve.getOrder(); ++p)
  {
    EXPECT_EQ(controlPoints[p], bCurve[p]);
  }

  bCurve.clear();
  EXPECT_EQ(-1, bCurve.getOrder());
  EXPECT_EQ(0, bCurve.getNumControlPoints());
  EXPECT_FALSE(bCurve.isRational());
}

//------------------------------------------------------------------------------
TEST(primal_beziercurve, evaluate)
{
  SLIC_INFO("Testing Bezier evaluation");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  const int max_order = 3;
  PointType data[max_order + 1] = {PointType {0.6, 1.2, 1.0},
                                   PointType {1.3, 1.6, 1.8},
                                   PointType {2.9, 2.4, 2.3},
                                   PointType {3.2, 3.5, 3.0}};

  // clang-format off
  PointType exp_vals[4][3] = {{PointType {0.6, 1.2, 1.0}, PointType {  0.6,    1.2,    1.0}, PointType {0.6, 1.2, 1.0}},
                              {PointType {0.6, 1.2, 1.0}, PointType { 0.95,    1.4,    1.4}, PointType {1.3, 1.6, 1.8}},
                              {PointType {0.6, 1.2, 1.0}, PointType {1.525,    1.7,  1.725}, PointType {2.9, 2.4, 2.3}},
                              {PointType {0.6, 1.2, 1.0}, PointType { 2.05, 2.0875, 2.0375}, PointType {3.2, 3.5, 3.0}}};
  // clang-format on

  // Test evaluation at various orders, using the first `ord`
  //  points in `data` as control points.
  for(int ord = 0; ord <= max_order; ++ord)
  {
    BezierCurveType curve(data, ord);

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
TEST(primal_beziercurve_, first_derivatives)
{
  SLIC_INFO("Testing Bezier derivative calculation");

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

  // clang-format off
  VectorType exp_vals[4][3] = {{VectorType {0.0, 0.0, 0.0}, VectorType { 0.0,   0.0,   0.0}, VectorType {0.0, 0.0, 0.0}},
                               {VectorType {0.7, 0.4, 0.8}, VectorType { 0.7,   0.4,   0.8}, VectorType {0.7, 0.4, 0.8}},
                               {VectorType {1.4, 0.8, 1.6}, VectorType { 2.3,   1.2,   1.3}, VectorType {3.2, 1.6, 1.0}},
                               {VectorType {2.1, 1.2, 2.4}, VectorType {3.15, 2.325, 1.875}, VectorType {0.9, 3.3, 2.1}}};
  // clang-format on

  // Test derivative calculation at various orders, using the first `ord`
  //  points in `data` as control points.
  for(int ord = 0; ord <= max_order; ++ord)
  {
    BezierCurveType curve(data, ord);

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
TEST(primal_beziercurve_, second_derivative)
{
  SLIC_INFO("Testing Bezier second derivative calculation");

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

  // clang-format off
  VectorType exp_vals[4][3] = {{VectorType {0.0, 0.0,  0.0}, VectorType { 0.0, 0.0,  0.0}, VectorType { 0.0, 0.0,  0.0}},
                               {VectorType {0.0, 0.0,  0.0}, VectorType { 0.0, 0.0,  0.0}, VectorType { 0.0, 0.0,  0.0}},
                               {VectorType {1.8, 0.8, -0.6}, VectorType { 1.8, 0.8, -0.6}, VectorType { 1.8, 0.8, -0.6}},
                               {VectorType {5.4, 2.4, -1.8}, VectorType {-1.2, 2.1, -0.3}, VectorType {-7.8, 1.8,  1.2}}};
  // clang-format on

  // Test evaluation at various orders, using the first `ord`
  //  points in `data` as control points.
  for(int ord = 0; ord <= max_order; ++ord)
  {
    BezierCurveType curve(data, ord);

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
TEST(primal_beziercurve_, batch_derivatives)
{
  SLIC_INFO("Testing Bezier batched derivative calculations");

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

  const double t = 0.6;
  PointType batch1_val, batch2_val;
  VectorType batch1_dt, batch2_dt, batch2_dtdt;

  // Test the derivative calculations on various orders
  for(int ord = 0; ord <= max_order; ++ord)
  {
    BezierCurveType curve(data, ord);

    curve.evaluateFirstDerivative(t, batch1_val, batch1_dt);
    curve.evaluateSecondDerivative(t, batch2_val, batch2_dt, batch2_dtdt);
    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(curve.evaluate(t)[i], batch1_val[i], 1e-15);
      EXPECT_NEAR(curve.dt(t)[i], batch1_dt[i], 1e-15);

      EXPECT_NEAR(curve.evaluate(t)[i], batch2_val[i], 1e-15);
      EXPECT_NEAR(curve.dt(t)[i], batch2_dt[i], 1e-15);
      EXPECT_NEAR(curve.dtdt(t)[i], batch2_dtdt[i], 1e-15);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziercurve, split_cubic)
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
  BezierCurveType b2Curve(data, order);

  BezierCurveType b3Curve(order);  // Checks split with order constructor
  BezierCurveType b4Curve;         // Checks split with default constructor
  b2Curve.split(.5, b3Curve, b4Curve);

  // clang-format off
  PointType b3Coords[4] = {PointType {  0.6,    1.2,    1.0}, 
                           PointType {  .95,    1.4,    1.4},
                           PointType {1.525,    1.7,  1.725}, 
                           PointType { 2.05, 2.0875, 2.0375}};
                           
                           
  PointType b4Coords[4] = {PointType{  2.05, 2.0875, 2.0375},
                           PointType{ 2.575,  2.475,   2.35},
                           PointType{  3.05,   2.95,   2.65},
                           PointType{   3.2,    3.5,    3.0}};
  // clang-format on        

  BezierCurveType b3True(b3Coords, 3);
  BezierCurveType b4True(b4Coords, 3);
  for(int i = 0; i < DIM; ++i)
  {
    for(int p = 0; p <= order; ++p)
    {
      EXPECT_DOUBLE_EQ(b3Curve[p][i], b3True[p][i]);
      EXPECT_DOUBLE_EQ(b4Curve[p][i], b4True[p][i]);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziercurve, split_point)
{
  SLIC_INFO("Testing Bezier splitting for order 0");

  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  // Test order-0
  {
    const int order = 0;
    PointType data[order + 1] = {PointType::make_point(0.6, 1.2)};
    BezierCurveType b(data, order);

    BezierCurveType c1, c2;
    b.split(0.5, c1, c2);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_DOUBLE_EQ(data[0][i], c1[0][i]);
      EXPECT_DOUBLE_EQ(data[0][i], c2[0][i]);
    }
  }
}

TEST(primal_beziercurve, split_linear)
{
  SLIC_INFO("Testing Bezier splitting for order 1");

  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  const int order = 1;
  PointType data[order + 1] = {PointType {-1, -5}, PointType {1, 5}};
  BezierCurveType b(data, order);

  {
    BezierCurveType c1, c2;
    b.split(0.5, c1, c2);

    EXPECT_DOUBLE_EQ(-1., c1[0][0]);
    EXPECT_DOUBLE_EQ(-5., c1[0][1]);
    EXPECT_DOUBLE_EQ(0., c1[1][0]);
    EXPECT_DOUBLE_EQ(0., c1[1][1]);

    EXPECT_DOUBLE_EQ(0., c2[0][0]);
    EXPECT_DOUBLE_EQ(0., c2[0][1]);
    EXPECT_DOUBLE_EQ(1., c2[1][0]);
    EXPECT_DOUBLE_EQ(5., c2[1][1]);
  }

  {
    BezierCurveType c1, c2;
    const double t = 0.25;
    b.split(0.25, c1, c2);

    PointType interp = PointType::lerp(data[0], data[1], t);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_DOUBLE_EQ(data[0][i], c1[0][i]);
      EXPECT_DOUBLE_EQ(interp[i], c1[1][i]);

      EXPECT_DOUBLE_EQ(interp[i], c2[0][i]);
      EXPECT_DOUBLE_EQ(data[1][i], c2[1][i]);
    }
  }
}

TEST(primal_beziercurve, split_quadratic)
{
  SLIC_INFO("Testing Bezier splitting for order 2");

  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  const double t = .42;
  const int order = 2;

  // Control points for the three levels of the quadratic de Casteljau algorithm
  PointType lev0[3] = {PointType {1.1, 1.1},
                       PointType {5.5, 5.5},
                       PointType {9.9, 2.2}};

  PointType lev1[2] = {PointType::lerp(lev0[0], lev0[1], t),
                       PointType::lerp(lev0[1], lev0[2], t)};

  PointType lev2[1] = {PointType::lerp(lev1[0], lev1[1], t)};

  BezierCurveType b(lev0, order);

  // Define expected control points for curves 1 and 2
  BezierCurveType expC1(order);
  expC1[0] = lev0[0];
  expC1[1] = lev1[0];
  expC1[2] = lev2[0];

  BezierCurveType expC2(order);
  expC2[0] = lev2[0];
  expC2[1] = lev1[1];
  expC2[2] = lev0[2];

  // Split the curve
  BezierCurveType c1, c2;
  b.split(t, c1, c2);

  SLIC_INFO(""                                          //
            << "Original quadratic: " << b              //
            << "\nCurves after splitting at t = " << t  //
            << "\n\t c1: " << c1                        //
            << "\n\t c2: " << c2);

  // Check values
  for(int p = 0; p <= order; ++p)
  {
    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_DOUBLE_EQ(expC1[p][i], c1[p][i]);
      EXPECT_DOUBLE_EQ(expC2[p][i], c2[p][i]);
    }
  }
}

TEST(primal_beziercurve, isLinear)
{
  SLIC_INFO("Testing isLinear() on Bezier curves");

  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using SegmentType = primal::Segment<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  // order 0 -- always true
  {
    const int order = 0;
    auto curve = BezierCurveType(order);
    EXPECT_TRUE(curve.isLinear());

    curve[0] = PointType::make_point(1., 1.);
    EXPECT_TRUE(curve.isLinear());
  }

  // order 1 -- always true
  {
    const int order = 1;
    auto curve = BezierCurveType(order);
    EXPECT_TRUE(curve.isLinear());

    curve[0] = PointType {1., 1.8};
    curve[1] = PointType {-12., 3.5};
    EXPECT_TRUE(curve.isLinear());
  }

  // order 2
  {
    const int order = 2;
    auto curve = BezierCurveType(order);
    EXPECT_TRUE(curve.isLinear());

    // straight line
    curve[0] = PointType {1, 1};
    curve[1] = PointType {2, 2};
    curve[2] = PointType {3, 3};
    EXPECT_TRUE(curve.isLinear());

    // move middle point and check linearity with different tolerances
    VectorType v(curve[2], curve[0]);
    auto normal = VectorType {-v[1], v[0]};
    curve[1].array() += 0.005 * normal.array();
    SLIC_INFO("Updated curve: " << curve);

    SegmentType s(curve[2], curve[0]);
    SLIC_INFO("Sq dist: " << primal::squared_distance(curve[1], s));

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

TEST(primal_beziercurve, reverseOrientation)
{
  SLIC_INFO("Testing reverseOrientation() on Bezier curves");

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
      for(int i = 0; i <= order; ++i)
      {
        pts[i] = PointType(i);
      }
      BezierCurveType curve(pts.data(), order);

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

//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
