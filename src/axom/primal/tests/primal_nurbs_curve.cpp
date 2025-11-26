// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! 
 * \file primal_nurbs_curve.cpp
 * \brief This file tests primal's NURBS curve functionality
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/primal/geometry/NURBSCurve.hpp"
#include <math.h>

namespace primal = axom::primal;

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, default_constructor)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  auto check_curve = [](const NURBSCurveType& cur) {
    constexpr axom::IndexType zero {0};

    EXPECT_EQ(-1, cur.getDegree());
    EXPECT_EQ(zero, cur.getOrder());
    EXPECT_EQ(zero, cur.getControlPoints().size());
    EXPECT_EQ(zero, cur.getKnots().getArray().size());
    EXPECT_FALSE(cur.isValidNURBS());
  };

  // check default constructor
  {
    NURBSCurveType nCurve;
    check_curve(nCurve);
  }

  // explicitly specify degree of -1
  {
    NURBSCurveType nCurve(-1);
    check_curve(nCurve);
  }

  // explicitly specify degree of -1 and 0 points
  {
    NURBSCurveType nCurve(0, -1);
    check_curve(nCurve);
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, sizing_constructors)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  for(int deg = 0; deg < 5; ++deg)
  {
    for(int npts = deg + 1; npts < 10; ++npts)
    {
      NURBSCurveType nCurve(npts, deg);

      EXPECT_TRUE(nCurve.isValidNURBS());

      EXPECT_EQ(deg, nCurve.getDegree());
      EXPECT_EQ(deg + 1, nCurve.getOrder());
      EXPECT_EQ(npts, nCurve.getControlPoints().size());
      EXPECT_EQ(npts + deg + 1, nCurve.getKnots().getArray().size());
      EXPECT_FALSE(nCurve.isRational());
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, bezier_constructors)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  constexpr int order = 2;

  axom::Array<PointType> controlPoints {PointType {0.6, 1.2, 1.0},
                                        PointType {0.0, 1.6, 1.8},
                                        PointType {0.2, 1.4, 2.0}};
  axom::Array<double> weights {.25, .5, .75};

  primal::BezierCurve<double, DIM> nBez(controlPoints, order);
  EXPECT_FALSE(nBez.isRational());

  primal::BezierCurve<double, DIM> rBez(controlPoints, weights, order);
  EXPECT_TRUE(rBez.isRational());

  for(const auto& bez : {nBez, rBez})
  {
    NURBSCurveType cur(bez);
    EXPECT_TRUE(cur.isValidNURBS());
    EXPECT_EQ(cur.isRational(), bez.isRational());

    EXPECT_EQ(order, cur.getDegree());
    EXPECT_EQ(order + 1, cur.getOrder());
    EXPECT_EQ(controlPoints.size(), cur.getControlPoints().size());
    EXPECT_EQ(controlPoints.size() + order + 1, cur.getKnots().getArray().size());

    for(int p = 0; p < controlPoints.size(); ++p)
    {
      EXPECT_EQ(controlPoints[p], cur[p]);
      if(cur.isRational())
      {
        EXPECT_EQ(weights[p], cur.getWeight(p));
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, knotless_array_constructors)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  constexpr int npts = 3;
  constexpr int degree = 1;

  // Construct from C-Style arrays
  PointType controlPoints[npts] = {PointType {0.6, 1.2, 1.0},
                                   PointType {0.0, 1.6, 1.8},
                                   PointType {0.2, 1.4, 2.0}};

  double weights[npts] = {1.0, 2.0, 3.0};

  auto check_curve = [&](const NURBSCurveType& cur, bool expect_rational) {
    EXPECT_TRUE(cur.isValidNURBS());

    EXPECT_EQ(degree, cur.getDegree());
    EXPECT_EQ(degree + 1, cur.getOrder());
    EXPECT_EQ(npts, cur.getControlPoints().size());
    EXPECT_EQ(npts + degree + 1, cur.getKnots().getArray().size());

    EXPECT_EQ(cur.isRational(), expect_rational);

    for(int p = 0; p < cur.getNumControlPoints(); ++p)
    {
      EXPECT_EQ(controlPoints[p], cur[p]);
      if(expect_rational)
      {
        EXPECT_EQ(weights[p], cur.getWeight(p));
      }
    }
  };

  // test constructors over C-arrays
  {
    NURBSCurveType nCurve(controlPoints, npts, degree);
    check_curve(nCurve, false);

    NURBSCurveType wCurve(controlPoints, weights, npts, degree);
    check_curve(wCurve, true);
  }

  // test constructors over ArrayViews over C-arrays
  {
    axom::ArrayView<PointType> cp(controlPoints, npts);
    axom::ArrayView<double> wt(weights, npts);

    NURBSCurveType nCurve(cp, degree);
    check_curve(nCurve, false);

    NURBSCurveType wCurve(cp, wt, degree);
    check_curve(wCurve, true);
  }

  // test over Axom arrays
  {
    axom::Array<PointType> cp;
    cp.assign(std::begin(controlPoints), std::end(controlPoints));

    axom::Array<double> wt;
    wt.assign(std::begin(weights), std::end(weights));

    // test over axom::Arrays
    {
      NURBSCurveType nCurve(cp, degree);
      check_curve(nCurve, false);

      NURBSCurveType wCurve(cp, wt, degree);
      check_curve(wCurve, true);
    }

    // test over ArrayViews from Arrays
    {
      NURBSCurveType nCurve(cp.view(), degree);
      check_curve(nCurve, false);
      NURBSCurveType wCurve(cp.view(), wt.view(), degree);
      check_curve(wCurve, true);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, knotted_array_constructor)
{
  SLIC_INFO("Testing knot array constructor");

  constexpr int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  constexpr int npts = 3;
  constexpr int degree = 1;
  constexpr int nkts = npts + degree + 1;

  // Construct from C-Style arrays
  PointType controlPoints[npts] = {PointType {0.6, 1.2, 1.0},
                                   PointType {0.0, 1.6, 1.8},
                                   PointType {0.2, 1.4, 2.0}};

  double weights[npts] = {1.0, 2.0, 3.0};
  double knots[nkts] = {0.0, 0.0, 0.2, 1.0, 1.0};

  auto check_curve = [&](const NURBSCurveType& cur, bool expect_rational) {
    EXPECT_TRUE(cur.isValidNURBS());

    EXPECT_EQ(degree, cur.getDegree());
    EXPECT_EQ(degree + 1, cur.getOrder());
    EXPECT_EQ(npts, cur.getControlPoints().size());
    EXPECT_EQ(nkts, cur.getKnots().getArray().size());

    EXPECT_EQ(cur.isRational(), expect_rational);

    for(int p = 0; p < cur.getNumControlPoints(); ++p)
    {
      EXPECT_EQ(controlPoints[p], cur[p]);
      if(expect_rational)
      {
        EXPECT_EQ(weights[p], cur.getWeight(p));
      }
    }
  };

  // test constructors over C-arrays
  {
    NURBSCurveType nCurve(controlPoints, npts, knots, nkts);
    NURBSCurveType wCurve(controlPoints, weights, npts, knots, nkts);

    check_curve(nCurve, false);
    check_curve(wCurve, true);
  }

  // test constructors over ArrayViews from C-arrays
  {
    axom::ArrayView<PointType> cp(controlPoints, npts);
    axom::ArrayView<double> wt(weights, npts);
    axom::ArrayView<double> kv(knots, nkts);

    NURBSCurveType nCurve(cp, degree);
    check_curve(nCurve, false);

    NURBSCurveType wCurve(cp, wt, degree);
    check_curve(wCurve, true);
  }

  // Construct from axom::Array
  {
    axom::Array<PointType> cp_arr;
    cp_arr.assign(std::begin(controlPoints), std::end(controlPoints));

    axom::Array<double> wt_arr;
    wt_arr.assign(std::begin(weights), std::end(weights));

    axom::Array<double> kt_arr;
    kt_arr.assign(std::begin(knots), std::end(knots));

    primal::KnotVector<double> knotvec(kt_arr, degree);

    // test using Axom Arrays
    {
      NURBSCurveType nCurve(cp_arr, kt_arr);
      check_curve(nCurve, false);

      NURBSCurveType wCurve(cp_arr, wt_arr, kt_arr);
      check_curve(wCurve, true);
    }

    // test using Axom arrays and KnotVectors
    {
      NURBSCurveType nCurve(cp_arr, knotvec);
      check_curve(nCurve, false);

      NURBSCurveType wCurve(cp_arr, wt_arr, knotvec);
      check_curve(wCurve, true);
    }

    // test using ArrayViews from Arrays as well as KnotVectors
    {
      NURBSCurveType nCurve(cp_arr.view(), knotvec);
      check_curve(nCurve, false);

      NURBSCurveType wCurve(cp_arr.view(), wt_arr.view(), knotvec);
      check_curve(wCurve, true);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, set_degree)
{
  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  SLIC_INFO("Test adding control points to empty NURBS curve");

  NURBSCurveType nCurve;
  EXPECT_EQ(-1, nCurve.getDegree());

  const int degree = 1;
  const int npts = 3;
  PointType controlPoints[3] = {PointType {0.6, 1.2, 1.0},
                                PointType {0.0, 1.6, 1.8},
                                PointType {0.2, 1.4, 2.0}};

  nCurve.setNumControlPoints(npts);
  nCurve.setDegree(degree);

  EXPECT_EQ(nCurve.getDegree(), degree);
  EXPECT_EQ(nCurve.getNumControlPoints(), npts);
  EXPECT_EQ(nCurve.getNumKnots(), npts + degree + 1);

  nCurve[0] = controlPoints[0];
  nCurve[1] = controlPoints[1];
  nCurve[2] = controlPoints[2];

  for(int p = 0; p < npts; ++p)
  {
    EXPECT_EQ(nCurve[p], controlPoints[p]);
  }

  nCurve.clear();
  EXPECT_EQ(-1, nCurve.getDegree());
  EXPECT_FALSE(nCurve.isRational());

  nCurve.setParameters(npts, degree);
  nCurve.makeRational();
  EXPECT_TRUE(nCurve.isRational());

  nCurve.setWeight(0, 2.0);
  EXPECT_EQ(nCurve.getWeight(0), 2.0);
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, evaluate)
{
  SLIC_INFO("Testing NURBS evaluation");

  using CoordType = double;
  using Point1D = primal::Point<CoordType, 1>;
  using Point2D = primal::Point<CoordType, 2>;
  using Point3D = primal::Point<CoordType, 3>;

  using NURBSCurve1D = primal::NURBSCurve<CoordType, 1>;
  using NURBSCurve2D = primal::NURBSCurve<CoordType, 2>;
  using NURBSCurve3D = primal::NURBSCurve<CoordType, 3>;

  const int max_degree = 3;
  Point1D data_1d[max_degree + 1] = {Point1D {0.6}, Point1D {1.3}, Point1D {2.9}, Point1D {3.2}};

  Point2D data_2d[max_degree + 1] = {Point2D {0.6, 1.2},
                                     Point2D {1.3, 1.6},
                                     Point2D {2.9, 2.4},
                                     Point2D {3.2, 3.5}};

  Point3D data_3d[max_degree + 1] = {Point3D {0.6, 1.2, 1.0},
                                     Point3D {1.3, 1.6, 1.8},
                                     Point3D {2.9, 2.4, 2.3},
                                     Point3D {3.2, 3.5, 3.0}};

  double weights[4] = {1.0, 2.0, 3.0, 4.0};

  // clang-format off
  Point3D exp_start_vals[4][4] =  // degree 0                  degree 1                   degree 2                   degree 3
           /* 1 pt */ {{Point3D {0.6, 1.2, 1.0}, Point3D {-999, -999, -999}, Point3D {-999, -999, -999}, Point3D {-999, -999, -999}},
           /* 2 pt */  {Point3D {0.6, 1.2, 1.0}, Point3D { 0.6,  1.2,  1.0}, Point3D {-999, -999, -999}, Point3D {-999, -999, -999}},
           /* 3 pt */  {Point3D {0.6, 1.2, 1.0}, Point3D { 0.6,  1.2,  1.0}, Point3D { 0.6,  1.2,  1.0}, Point3D {-999, -999, -999}},
           /* 4 pt */  {Point3D {0.6, 1.2, 1.0}, Point3D { 0.6,  1.2,  1.0}, Point3D { 0.6,  1.2,  1.0}, Point3D { 0.6,  1.2,  1.0}}};

  Point3D exp_mid_vals[4][4] =  // degree 0                         degree 1                             degree 2                       degree 3
           /* 1 pt */ {{Point3D {0.6, 1.2, 1.0}, Point3D {   -999,     -999,     -999}, Point3D {  -999,  -999,   -999}, Point3D { -999, -999,  -999}},
           /* 2 pt */  {Point3D {1.3, 1.6, 1.8}, Point3D {16./15.,  22./15.,  23./15.}, Point3D {  -999,  -999,   -999}, Point3D { -999, -999,  -999}},
           /* 3 pt */  {Point3D {1.3, 1.6, 1.8}, Point3D {    1.3,      1.6,      1.8}, Point3D {1.8125,  1.85, 1.8875}, Point3D { -999, -999,  -999}},
           /* 4 pt */  {Point3D {2.9, 2.4, 2.3}, Point3D {   2.26,     2.08,      2.1}, Point3D {  2.26,  2.08,    2.1}, Point3D {2.365, 2.32, 2.225}}};

  Point3D exp_end_vals[4][4] =     // degree 0                degree 1                     degree 2                   degree 3
           /* 1 pt */ {{Point3D {0.6, 1.2, 1.0}, Point3D {-999, -999, -999}, Point3D {-999, -999, -999}, Point3D {-999, -999, -999}},  
           /* 2 pt */  {Point3D {1.3, 1.6, 1.8}, Point3D { 1.3,  1.6,  1.8}, Point3D {-999, -999, -999}, Point3D {-999, -999, -999}},
           /* 3 pt */  {Point3D {2.9, 2.4, 2.3}, Point3D { 2.9,  2.4,  2.3}, Point3D { 2.9,  2.4,  2.3}, Point3D {-999, -999, -999}},
           /* 4 pt */  {Point3D {3.2, 3.5, 3.0}, Point3D { 3.2,  3.5,  3.0}, Point3D { 3.2,  3.5,  3.0}, Point3D { 3.2,  3.5,  3.0}}};
  // clang-format on

  // Test evaluation at various spatial dimensions, various degrees,
  //  and various number of control points in `data`.
  for(int npts = 1; npts <= 4; ++npts)
  {
    for(int deg = 0; deg <= npts - 1; ++deg)
    {
      // 1D NURBS Curve
      NURBSCurve1D curve1(data_1d, weights, npts, deg);

      Point1D calc_start1 = curve1.evaluate(0.0);
      Point1D calc_mid1 = curve1.evaluate(0.5);
      Point1D calc_end1 = curve1.evaluate(1.0);

      for(int i = 0; i < 1; ++i)
      {
        EXPECT_NEAR(calc_start1[i], exp_start_vals[npts - 1][deg][i], 1e-15);
        EXPECT_NEAR(calc_mid1[i], exp_mid_vals[npts - 1][deg][i], 1e-15);
        EXPECT_NEAR(calc_end1[i], exp_end_vals[npts - 1][deg][i], 1e-15);
      }

      // 2D NURBS Curve
      NURBSCurve2D curve2(data_2d, weights, npts, deg);

      Point2D calc_start2 = curve2.evaluate(0.0);
      Point2D calc_mid2 = curve2.evaluate(0.5);
      Point2D calc_end2 = curve2.evaluate(1.0);

      for(int i = 0; i < 2; ++i)
      {
        EXPECT_NEAR(calc_start2[i], exp_start_vals[npts - 1][deg][i], 1e-15);
        EXPECT_NEAR(calc_mid2[i], exp_mid_vals[npts - 1][deg][i], 1e-15);
        EXPECT_NEAR(calc_end2[i], exp_end_vals[npts - 1][deg][i], 1e-15);
      }

      // 3D NURBS Curve
      NURBSCurve3D curve3(data_3d, weights, npts, deg);

      Point3D calc_start3 = curve3.evaluate(0.0);
      Point3D calc_mid3 = curve3.evaluate(0.5);
      Point3D calc_end3 = curve3.evaluate(1.0);

      for(int i = 0; i < 3; ++i)
      {
        EXPECT_NEAR(calc_start3[i], exp_start_vals[npts - 1][deg][i], 1e-15);
        EXPECT_NEAR(calc_mid3[i], exp_mid_vals[npts - 1][deg][i], 1e-15);
        EXPECT_NEAR(calc_end3[i], exp_end_vals[npts - 1][deg][i], 1e-15);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, first_derivatives)
{
  SLIC_INFO("Testing NURBS derivative calculation");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  const int max_degree = 3;
  PointType data[max_degree + 1] = {PointType {0.6, 1.2, 1.0},
                                    PointType {1.3, 1.6, 1.8},
                                    PointType {2.9, 2.4, 2.3},
                                    PointType {3.2, 3.5, 3.0}};

  double weights[4] = {1.0, 2.0, 3.0, 4.0};

  // clang-format off
  VectorType exp_start_vals[4][4] =  // degree 0                  degree 1                      degree 2                   degree 3
           /* 1 pt */ {{VectorType {0.0, 0.0, 0.0}, VectorType {-999, -999, -999}, VectorType {-999, -999, -999}, VectorType {-999, -999, -999}},
           /* 2 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType { 1.4,  0.8,  1.6}, VectorType {-999, -999, -999}, VectorType {-999, -999, -999}},
           /* 3 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType { 2.8,  1.6,  3.2}, VectorType { 2.8,  1.6,  3.2}, VectorType {-999, -999, -999}},
           /* 4 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType { 4.2,  2.4,  4.8}, VectorType { 5.6,  3.2,  6.4}, VectorType { 4.2,  2.4,  4.8}}};

  VectorType exp_mid_vals[4][4] =  // degree 0                    degree 1                    degree 2                         degree 3
           /* 1 pt */ {{VectorType {0.0, 0.0, 0.0}, VectorType {   -999,    -999,    -999}, VectorType {  -999,  -999,   -999}, VectorType { -999,  -999, -999}},
           /* 2 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {28./45., 16./45., 32./45.}, VectorType {  -999,  -999,   -999}, VectorType { -999,  -999, -999}},
           /* 3 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {    4.8,     2.4,     1.5}, VectorType {2.2375,  1.15, 1.0625}, VectorType { -999,  -999, -999}},
           /* 4 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {  4.608,   2.304,    1.44}, VectorType { 3.072, 1.536,   0.96}, VectorType {2.652, 2.256, 1.62}}};

  VectorType exp_end_vals[4][4] =  // degree 0                  degree 1                   degree 2                   degree 3
           /* 1 pt */ {{VectorType {0.0, 0.0, 0.0}, VectorType {   -999,    -999,  -999}, VectorType {   -999,    -999,  -999}, VectorType { -999,  -999,  -999}},
           /* 2 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {   0.35,     0.2,   0.4}, VectorType {   -999,    -999,  -999}, VectorType { -999,  -999,  -999}},
           /* 3 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {32./15., 16./15., 2./3.}, VectorType {32./15., 16./15., 2./3.}, VectorType { -999,  -999,  -999}},
           /* 4 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {  0.675,   2.475, 1.575}, VectorType {    0.9,     3.3,   2.1}, VectorType {0.675, 2.475, 1.575}}};
  // clang-format on

  // Test evaluation at various degrees and various number of control
  //  points in `data`.
  for(int npts = 1; npts <= 4; ++npts)
  {
    for(int deg = 0; deg <= npts - 1; ++deg)
    {
      NURBSCurveType curve(data, weights, npts, deg);

      VectorType calc_start = curve.dt(0.0);
      VectorType calc_mid = curve.dt(0.5);
      VectorType calc_end = curve.dt(1.0);

      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(calc_start[i], exp_start_vals[npts - 1][deg][i], 1e-13);
        EXPECT_NEAR(calc_mid[i], exp_mid_vals[npts - 1][deg][i], 1e-13);
        EXPECT_NEAR(calc_end[i], exp_end_vals[npts - 1][deg][i], 1e-13);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, second_derivatives)
{
  SLIC_INFO("Testing NURBS second derivative calculation");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  const int max_degree = 3;
  PointType data[max_degree + 1] = {PointType {0.6, 1.2, 1.0},
                                    PointType {1.3, 1.6, 1.8},
                                    PointType {2.9, 2.4, 2.3},
                                    PointType {3.2, 3.5, 3.0}};

  double weights[4] = {1.0, 2.0, 3.0, 4.0};

  // clang-format off
  VectorType exp_start_vals[4][4] =  // degree 0                  degree 1                   degree 2                   degree 3
           /* 1 pt */ {{VectorType {0.0, 0.0, 0.0},  VectorType { -999,  -999,  -999},  VectorType { -999,  -999,  -999},  VectorType {-999, -999,  -999}},
           /* 2 pt */  {VectorType {0.0, 0.0, 0.0},  VectorType { -2.8,  -1.6,  -3.2},  VectorType { -999,  -999,  -999},  VectorType {-999, -999,  -999}},
           /* 3 pt */  {VectorType {0.0, 0.0, 0.0},  VectorType {-11.2,  -6.4, -12.8},  VectorType { -3.0,  -2.4, -11.4},  VectorType {-999, -999,  -999}},
           /* 4 pt */  {VectorType {0.0, 0.0, 0.0},  VectorType {-25.2, -14.4, -28.8},  VectorType {-34.0, -20.8, -54.8},  VectorType {-0.6, -2.4, -24.6}}};

  VectorType exp_mid_vals[4][4] =  // degree 0                     degree 1                    degree 2                         degree 3
           /* 1 pt */ {{VectorType {0.0, 0.0, 0.0}, VectorType {      -999,       -999,       -999}, VectorType {   -999,   -999,   -999}, VectorType {   -999,   -999,   -999}},
           /* 2 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {-112./135., -64.0/135., -128./135.}, VectorType {   -999,   -999,   -999}, VectorType {   -999,   -999,   -999}},
           /* 3 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {      -9.6,       -4.8,       -3.0}, VectorType { -0.375,   -0.3, -1.425}, VectorType {   -999,   -999,   -999}},
           /* 4 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {  -11.0592,    -5.5296,     -3.456}, VectorType {-5.1712, 9.5744,  6.144}, VectorType {-3.8448, 0.3456, -0.888}}};

  VectorType exp_end_vals[4][4] =  // degree 0                             degree 1                               degree 2                              degree 3
           /* 1 pt */ {{VectorType {0.0, 0.0, 0.0}, VectorType {     -999,     -999,    -999}, VectorType {  -999,    -999,     -999}, VectorType {   -999,   -999,   -999}},
           /* 2 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {    -0.35,     -0.2,    -0.4}, VectorType {  -999,    -999,     -999}, VectorType {   -999,   -999,   -999}},
           /* 3 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {-128./45., -64./45.,  -8./9.}, VectorType {-1./9., -8./90., -19./45.}, VectorType {   -999,   -999,   -999}},
           /* 4 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {  -1.0125,  -3.7125, -2.3625}, VectorType {  -2.9,    -0.5,     -0.3}, VectorType {-4.0125, 0.4875, 0.3375}}};
  // clang-format on

  // Test evaluation at various degrees and various number of control
  //  points in `data`.
  for(int npts = 1; npts <= 4; ++npts)
  {
    for(int deg = 0; deg <= npts - 1; ++deg)
    {
      NURBSCurveType curve(data, weights, npts, deg);

      VectorType calc_start = curve.dtdt(0.0);
      VectorType calc_mid = curve.dtdt(0.5);
      VectorType calc_end = curve.dtdt(1.0);

      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(calc_start[i], exp_start_vals[npts - 1][deg][i], 1e-13);
        EXPECT_NEAR(calc_mid[i], exp_mid_vals[npts - 1][deg][i], 1e-13);
        EXPECT_NEAR(calc_end[i], exp_end_vals[npts - 1][deg][i], 1e-13);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, knot_insertion)
{
  SLIC_INFO("Testing NURBS knot insertion");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  PointType data[4] = {PointType {0.6, 1.2, 1.0},
                       PointType {1.3, 1.6, 1.8},
                       PointType {2.9, 2.4, 2.3},
                       PointType {3.2, 3.5, 3.0}};

  double weights[4] = {1.0, 2.0, 3.0, 4.0};

  NURBSCurveType curve(data, weights, 4, 3);
  NURBSCurveType curve_knots(data, weights, 4, 3);

  // Knot insertion shouldn't change the paramterization or position
  //  of the curve.

  // Insert a knot at 0.5
  auto num_knots = curve_knots.getNumKnots();
  curve_knots.insertKnot(0.5, 3);
  EXPECT_EQ(curve_knots.getNumKnots(), num_knots + 3);
  for(double t = 0.0; t <= 1.0; t += 0.1)
  {
    PointType p = curve.evaluate(t);
    PointType p_knots = curve_knots.evaluate(t);
    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(p[i], p_knots[i], 1e-13);
    }
  }

  num_knots = curve_knots.getNumKnots();
  curve_knots.insertKnot(0.7, 1);
  EXPECT_EQ(curve_knots.getNumKnots(), num_knots + 1);
  for(double t = 0.0; t <= 1.0; t += 0.1)
  {
    PointType p = curve.evaluate(t);
    PointType p_knots = curve_knots.evaluate(t);
    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(p[i], p_knots[i], 1e-13);
    }
  }

  num_knots = curve_knots.getNumKnots();
  curve_knots.insertKnot(0.4, 2);
  EXPECT_EQ(curve_knots.getNumKnots(), num_knots + 2);
  for(double t = 0.0; t <= 1.0; t += 0.1)
  {
    PointType p = curve.evaluate(t);
    PointType p_knots = curve_knots.evaluate(t);
    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(p[i], p_knots[i], 1e-13);
    }
  }

  // Inserting knots shouldn't increase its multiplicty
  //  greater than the degree
  num_knots = curve_knots.getNumKnots();
  curve_knots.insertKnot(0.5, 5);  // Already has degree 3
  curve_knots.insertKnot(0.0, 5);  // Already has degree 4
  curve_knots.insertKnot(1.0, 5);  // Already has degree 4
  EXPECT_EQ(curve_knots.getNumKnots(), num_knots);

  // This method inserts knots with a target degree
  // This wont't change the knot vector since 0.4 is already inserted twice
  curve_knots.insertKnot(0.4, 1);
  curve_knots.insertKnot(0.4, 2);
  curve_knots.insertKnot(0.4, 1);
  curve_knots.insertKnot(0.4, 2);
  EXPECT_EQ(curve_knots.getNumKnots(), num_knots);

  for(double t = 0.0; t <= 1.0; t += 0.1)
  {
    PointType p = curve.evaluate(t);
    PointType p_knots = curve_knots.evaluate(t);
    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(p[i], p_knots[i], 1e-13);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, curve_splitting)
{
  SLIC_INFO("Testing NURBS curve splitting");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  PointType data[4] = {PointType {0.6, 1.2, 1.0},
                       PointType {1.3, 1.6, 1.8},
                       PointType {2.9, 2.4, 2.3},
                       PointType {3.2, 3.5, 3.0}};

  double weights[4] = {1.0, 2.0, 3.0, 4.0};

  for(int deg = 1; deg <= 3; ++deg)
  {
    NURBSCurveType curve(data, weights, 4, deg), dummy1, dummy2;

    // Do some knot insertion to make it interesting
    curve.insertKnot(0.3, 2);
    curve.insertKnot(0.6, 1);
    curve.insertKnot(0.8, 1);

    NURBSCurveType curve1, curve2, curve3;
    curve.split(0.3, curve1, curve2);
    curve2.split(0.6, curve2, curve3);

    for(double t = 0.0; t < 1.0; t += 0.05)
    {
      PointType p = curve.evaluate(t);

      if(t <= 0.3)
      {
        PointType p1 = curve1.evaluate(t);
        for(int i = 0; i < DIM; ++i)
        {
          EXPECT_NEAR(p[i], p1[i], 1e-13);
        }
      }

      if(t >= 0.3 && t <= 0.6)
      {
        PointType p2 = curve2.evaluate(t);
        for(int i = 0; i < DIM; ++i)
        {
          EXPECT_NEAR(p[i], p2[i], 1e-13);
        }
      }

      if(t >= 0.6)
      {
        PointType p3 = curve3.evaluate(t);
        for(int i = 0; i < DIM; ++i)
        {
          EXPECT_NEAR(p[i], p3[i], 1e-13);
        }
      }
    }
  }

  return;
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, bezier_extraction)
{
  SLIC_INFO("Testing NURBS Bezier extraction");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  PointType data[4] = {PointType {0.6, 1.2, 1.0},
                       PointType {1.3, 1.6, 1.8},
                       PointType {2.9, 2.4, 2.3},
                       PointType {3.2, 3.5, 3.0}};

  double weights[4] = {1.0, 2.0, 3.0, 4.0};

  NURBSCurveType curve(data, weights, 4, 3);

  // Do knot insertion, which determines where the Bezier
  //  splitting happens
  curve.insertKnot(0.33, 3);
  curve.insertKnot(0.66, 1);
  curve.insertKnot(0.77, 2);

  auto bezier_list = curve.extractBezier();

  EXPECT_EQ(bezier_list.size(), 4);

  for(double t = 0.0; t < 1.0; t += 0.05)
  {
    PointType p = curve.evaluate(t);

    if(t <= 0.33)
    {
      PointType p1 = bezier_list[0].evaluate(t / 0.33);
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(p[i], p1[i], 1e-13);
      }
    }

    if(t >= 0.33 && t <= 0.66)
    {
      PointType p2 = bezier_list[1].evaluate((t - 0.33) / (0.66 - 0.33));
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(p[i], p2[i], 1e-13);
      }
    }

    if(t >= 0.66 && t <= 0.77)
    {
      PointType p3 = bezier_list[2].evaluate((t - 0.66) / (0.77 - 0.66));
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(p[i], p3[i], 1e-13);
      }
    }

    if(t >= 0.77)
    {
      PointType p4 = bezier_list[3].evaluate((t - 0.77) / (1 - 0.77));
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(p[i], p4[i], 1e-13);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, bezier_extraction_zero)
{
  SLIC_INFO("Testing NURBS Bezier extraction on degree 0 curve");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  PointType data[4] = {PointType {0.6, 1.2, 1.0},
                       PointType {1.3, 1.6, 1.8},
                       PointType {2.9, 2.4, 2.3},
                       PointType {3.2, 3.5, 3.0}};

  double weights[4] = {1.0, 2.0, 3.0, 4.0};

  NURBSCurveType curve(data, weights, 4, 0);

  NURBSCurveType bezier_curve;
  auto bezier_list = curve.extractBezier();

  EXPECT_EQ(bezier_list.size(), curve.getNumControlPoints());

  for(double t = 0.0; t < 1.0; t += 0.05)
  {
    PointType p = curve.evaluate(t);

    if(t < 0.25)
    {
      PointType pt = bezier_list[0].evaluate(t / 0.25);
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(p[i], pt[i], 1e-13);
      }
    }

    if(t > 0.25 && t < 0.5)
    {
      PointType pt = bezier_list[1].evaluate((t - 0.25) / (0.5 - 0.25));
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(p[i], pt[i], 1e-13);
      }
    }

    if(t > 0.5 && t < 0.75)
    {
      PointType pt = bezier_list[2].evaluate((t - 0.5) / (0.75 - 0.5));
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(p[i], pt[i], 1e-13);
      }
    }

    if(t > 0.75)
    {
      PointType pt = bezier_list[3].evaluate((t - 0.75) / (1 - 0.75));
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(p[i], pt[i], 1e-13);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, bezier_extraction_full)
{
  SLIC_INFO("Testing NURBS Bezier extraction on a 'Bezier' curve");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  PointType data[4] = {PointType {0.6, 1.2, 1.0},
                       PointType {1.3, 1.6, 1.8},
                       PointType {2.9, 2.4, 2.3},
                       PointType {3.2, 3.5, 3.0}};

  double weights[4] = {1.0, 2.0, 3.0, 4.0};

  NURBSCurveType curve(data, weights, 4, 3);

  NURBSCurveType bezier_curve;
  auto bezier_list = curve.extractBezier();

  EXPECT_EQ(bezier_list.size(), 1);

  for(double t = 0.0; t < 1.0; t += 0.05)
  {
    PointType p = curve.evaluate(t);
    PointType pt = bezier_list[0].evaluate(t);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(p[i], pt[i], 1e-13);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, nurbs_reverse_orientation)
{
  SLIC_INFO("Testing NURBS reverse orientation");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  PointType data[4] = {PointType {0.6, 1.2, 1.0},
                       PointType {1.3, 1.6, 1.8},
                       PointType {2.9, 2.4, 2.3},
                       PointType {3.2, 3.5, 3.0}};

  double weights[4] = {1.0, 2.0, 3.0, 4.0};

  NURBSCurveType curve(data, weights, 4, 3);

  // Insert some knots to stress test reversal
  curve.insertKnot(0.33, 3);
  curve.insertKnot(0.66, 1);
  curve.insertKnot(0.77, 2);

  NURBSCurveType curve_reversed(curve);
  curve_reversed.reverseOrientation();

  for(double t = 0.0; t < 1.0; t += 0.05)
  {
    PointType p = curve.evaluate(t);
    PointType p_reversed = curve_reversed.evaluate(1.0 - t);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(p[i], p_reversed[i], 1e-13);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, nurbs_knot_normalization)
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
  double scaled_knots[11];

  for(int i = 0; i < 11; ++i)
  {
    scaled_knots[i] = 2.0 * knots[i] - 5.0;
  }

  NURBSCurveType circle(data, weights, 7, knots, 11);
  NURBSCurveType scaled_circle(data, weights, 7, scaled_knots, 11);

  // Evaluate the curve along the parameterization of each, as they should match
  for(double t = 0.0; t <= 1.0; t += 0.1)
  {
    PointType p = circle.evaluate(t);
    PointType p_scaled = scaled_circle.evaluate(2.0 * t - 5.0);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(p[i], p_scaled[i], 1e-13);
    }
  }

  // Re-normalizing the scaled circle should return the original
  scaled_circle.normalize();
  for(double t = 0.0; t <= 1.0; t += 0.1)
  {
    PointType p = circle.evaluate(t);
    PointType p_scaled = scaled_circle.evaluate(t);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(p[i], p_scaled[i], 1e-13);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, circular_arc_constructor)
{
  // Define a nurbs curve that represents a circle
  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  PointType center {1.0, 2.0};
  double radius = 2.3;

  // Check trivial arc (start_theta == end_theta)
  auto invalid_circle =
    NURBSCurveType::make_circular_arc_nurbs(1.0, 1.0, center[0], center[1], radius);
  EXPECT_FALSE(invalid_circle.isValidNURBS());

  // clang-format off
  double start_theta[] = {0.0,     -1.0, 1.0,            2.0,            2.0};
  double end_theta[]   = {2.0*M_PI, 1.0, 1.0 + 2*M_PI/3, 2.0 + 4*M_PI/3, 5.0};
  // clang-format on

  constexpr int npts = 11;
  double t_pts[npts];
  axom::numerics::linspace(0.0, 1.0, t_pts, npts);

  for(int i = 0; i < 5; ++i)
  {
    // Check for counter-clockwise and clockwise arcs
    for(bool is_ccw : {true, false})
    {
      double start_theta_i = start_theta[i];
      double end_theta_i = end_theta[i];

      if(!is_ccw)
      {
        std::swap(start_theta_i, end_theta_i);
      }

      auto circle = NURBSCurveType::make_circular_arc_nurbs(start_theta_i,
                                                            end_theta_i,
                                                            center[0],
                                                            center[1],
                                                            radius);

      // Check the first endpoint of the curve
      PointType start = circle.evaluate(0.0);
      PointType start_ex = PointType {center[0] + radius * std::cos(start_theta_i),
                                      center[1] + radius * std::sin(start_theta_i)};

      EXPECT_NEAR(start[0], start_ex[0], 1e-13);
      EXPECT_NEAR(start[1], start_ex[1], 1e-13);

      // Check the second endpoint of the curve
      PointType end = circle.evaluate(1.0);
      PointType end_ex = PointType {center[0] + radius * std::cos(end_theta_i),
                                    center[1] + radius * std::sin(end_theta_i)};

      EXPECT_NEAR(end[0], end_ex[0], 1e-13);
      EXPECT_NEAR(end[1], end_ex[1], 1e-13);

      // Check the magnitude of the points elsewhere along the curve
      for(int j = 0; j < npts; ++j)
      {
        PointType p = circle.evaluate(t_pts[j]);

        double distance = std::sqrt((p[0] - center[0]) * (p[0] - center[0]) +
                                    (p[1] - center[1]) * (p[1] - center[1]));

        EXPECT_NEAR(distance, radius, 1e-13);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, linear_segment_constructor)
{
  // Define a nurbs curve that represents a circle
  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  constexpr int npts = 11;
  double t_pts[npts];
  axom::numerics::linspace(0.0, 1.0, t_pts, npts);

  PointType start {1.0, 2.0};
  PointType end {3.0, 4.0};
  NURBSCurveType line = NURBSCurveType::make_linear_segment_nurbs(start, end);

  // Check points along the curve
  for(int j = 0; j < npts; ++j)
  {
    PointType p = line.evaluate(t_pts[j]);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(p[i], start[i] + t_pts[j] * (end[i] - start[i]), 1e-13);
    }
  }

  // Check a curve with start == end
  end = start;
  line = NURBSCurveType::make_linear_segment_nurbs(start, end);

  for(int j = 0; j < npts; ++j)
  {
    PointType p = line.evaluate(t_pts[j]);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(p[i], start[i], 1e-13);
    }
  }
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
