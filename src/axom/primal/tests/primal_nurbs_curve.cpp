// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! 
 * \file primal_nurbs_curve.cpp
 * \brief This file tests primal's NURBS curve functionality
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/primal/geometry/NURBSCurve.hpp"
#include "axom/primal/operators/intersect.hpp"
#include <math.h>

namespace primal = axom::primal;

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, is_linear_predicate)
{
  constexpr int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  constexpr double tol = 1e-12;

  // Degree-1 segment
  {
    NURBSCurveType c(2, 1);
    c[0] = PointType {0.0, 0.0};
    c[1] = PointType {1.0, 0.0};
    EXPECT_TRUE(c.isLinear(tol));
  }

  // Degree-2, collinear control polygon
  {
    NURBSCurveType c(3, 2);
    c[0] = PointType {0.0, 0.0};
    c[1] = PointType {0.5, 0.0};
    c[2] = PointType {1.0, 0.0};
    EXPECT_TRUE(c.isLinear(tol));
  }

  // Degree-2, non-collinear interior point
  {
    NURBSCurveType c(3, 2);
    c[0] = PointType {0.0, 0.0};
    c[1] = PointType {0.5, 1e-3};
    c[2] = PointType {1.0, 0.0};
    EXPECT_FALSE(c.isLinear(tol));
  }

  // Rational, collinear should still be linear
  {
    NURBSCurveType c(3, 2);
    c[0] = PointType {0.0, 0.0};
    c[1] = PointType {0.5, 0.0};
    c[2] = PointType {1.0, 0.0};
    c.makeRational();
    c.setWeight(0, 1.0);
    c.setWeight(1, 2.0);
    c.setWeight(2, 0.5);
    EXPECT_TRUE(c.isLinear(tol));
  }

  // Strict mode requires a uniform control-point distribution along the endpoint segment
  {
    // evenly spaced -> strict linear
    NURBSCurveType c(3, 2);
    c[0] = PointType {0.0, 0.0};
    c[1] = PointType {0.5, 0.0};
    c[2] = PointType {1.0, 0.0};
    EXPECT_TRUE(c.isLinear(tol, /*useStrictLinear=*/true));
  }
  {
    // collinear but not evenly spaced -> not strict linear
    NURBSCurveType c(3, 2);
    c[0] = PointType {0.0, 0.0};
    c[1] = PointType {0.2, 0.0};
    c[2] = PointType {1.0, 0.0};
    EXPECT_TRUE(c.isLinear(tol, /*useStrictLinear=*/false));
    EXPECT_FALSE(c.isLinear(tol, /*useStrictLinear=*/true));
  }
  {
    // slight deviation within tolerance should still be linear (non-strict)
    NURBSCurveType c(3, 2);
    c[0] = PointType {0.0, 0.0};
    c[1] = PointType {0.5, 0.5e-6};  // squared distance 2.5e-13 < 1e-12
    c[2] = PointType {1.0, 0.0};
    EXPECT_TRUE(c.isLinear(tol, /*useStrictLinear=*/false));
  }
}

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

    for(axom::IndexType p = 0; p < cur.getNumControlPoints(); ++p)
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

    for(axom::IndexType p = 0; p < cur.getNumControlPoints(); ++p)
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
  // Define a nurbs curve that represents a line segment
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

//------------------------------------------------------------------------------
template <typename NURBS2D, typename NURBS3D, typename VectorType = typename NURBS3D::VectorType>
NURBS3D promoteTo3D(const NURBS2D& input, const VectorType& xvec, const VectorType& yvec)
{
  SLIC_ASSERT(NURBS2D::PointType::DIMENSION == 2);
  SLIC_ASSERT(NURBS3D::PointType::DIMENSION == 3);

  using PointType = typename NURBS3D::PointType;

  NURBS3D output(input.getNumControlPoints(), input.getDegree());
  for(int i = 0; i < input.getNumControlPoints(); ++i)
  {
    const auto& p2 = input[i];
    output[i] = PointType((xvec * p2[0] + yvec * p2[1]).data(), PointType::DIMENSION);
  }

  output.setKnots(input.getKnots());

  if(input.isRational())
  {
    output.makeRational();
    for(int i = 0; i < input.getNumControlPoints(); ++i)
    {
      output.setWeight(i, input.getWeight(i));
    }
  }

  return output;
}

template <typename T>
void checkCircularArcCurvature2D(T tol)
{
  using NURBSCurveType = primal::NURBSCurve<T, 2>;

  const T radius = 4.;
  const auto curve =
    NURBSCurveType::make_circular_arc_nurbs(T(0.15) * T(M_PI), T(1.85) * T(M_PI), T(0.), T(0.), radius);

  constexpr int numSamples = 17;
  for(int i = 0; i < numSamples; ++i)
  {
    const T t = static_cast<T>(i) / static_cast<T>(numSamples - 1);
    EXPECT_NEAR(curve.curvature(t), T(1.) / radius, tol);
  }
}

template <typename T>
void checkCircularArcCurvature3D(T tol)
{
  using NURBSCurve2D = primal::NURBSCurve<T, 2>;
  using NURBSCurve3D = primal::NURBSCurve<T, 3>;
  using Vector3D = primal::Vector<T, 3>;

  const T radius = 4.;
  const auto curve2d =
    NURBSCurve2D::make_circular_arc_nurbs(T(0.1) * T(M_PI), T(1.6) * T(M_PI), T(0.), T(0.), radius);

  const T a0 = T(M_PI) / T(4.);
  const T a1 = a0 + T(M_PI) / T(2.);
  const Vector3D uvec {std::cos(a0), T(0.), -std::sin(a0)};
  const Vector3D vvec {std::cos(a1), T(0.), -std::sin(a1)};
  const auto curve3d = promoteTo3D<NURBSCurve2D, NURBSCurve3D>(curve2d, uvec, vvec);

  constexpr int numSamples = 17;
  for(int i = 0; i < numSamples; ++i)
  {
    const T t = static_cast<T>(i) / static_cast<T>(numSamples - 1);
    EXPECT_NEAR(curve3d.curvature(t), T(1.) / radius, tol);
  }
}

template <typename T>
void checkCircularArcCurvatureDerivative2D(T tol)
{
  using NURBSCurveType = primal::NURBSCurve<T, 2>;

  const T radius = 4.;
  const auto curve =
    NURBSCurveType::make_circular_arc_nurbs(T(0.15) * T(M_PI), T(1.85) * T(M_PI), T(0.), T(0.), radius);

  constexpr int numSamples = 17;
  for(int i = 0; i < numSamples; ++i)
  {
    const T t = static_cast<T>(i) / static_cast<T>(numSamples - 1);
    EXPECT_NEAR(curve.curvatureDerivative(t), T(0.), tol);
  }
}

template <typename T>
void checkCircularArcCurvatureDerivative3D(T tol)
{
  using NURBSCurve2D = primal::NURBSCurve<T, 2>;
  using NURBSCurve3D = primal::NURBSCurve<T, 3>;
  using Vector3D = primal::Vector<T, 3>;

  const T radius = 4.;
  const auto curve2d =
    NURBSCurve2D::make_circular_arc_nurbs(T(0.1) * T(M_PI), T(1.6) * T(M_PI), T(0.), T(0.), radius);

  const T a0 = T(M_PI) / T(4.);
  const T a1 = a0 + T(M_PI) / T(2.);
  const Vector3D uvec {std::cos(a0), T(0.), -std::sin(a0)};
  const Vector3D vvec {std::cos(a1), T(0.), -std::sin(a1)};
  const auto curve3d = promoteTo3D<NURBSCurve2D, NURBSCurve3D>(curve2d, uvec, vvec);

  constexpr int numSamples = 17;
  for(int i = 0; i < numSamples; ++i)
  {
    const T t = static_cast<T>(i) / static_cast<T>(numSamples - 1);
    EXPECT_NEAR(curve3d.curvatureDerivative(t), T(0.), tol);
  }
}

TEST(primal_nurbscurve, curvature2d)
{
  checkCircularArcCurvature2D<float>(1.e-5F);
  checkCircularArcCurvature2D<double>(1.e-10);
}

TEST(primal_nurbscurve, curvature3d)
{
  checkCircularArcCurvature3D<float>(1.e-5F);
  checkCircularArcCurvature3D<double>(1.e-10);
}

TEST(primal_nurbscurve, curvature_derivative2d)
{
  checkCircularArcCurvatureDerivative2D<float>(1.e-5F);
  checkCircularArcCurvatureDerivative2D<double>(1.e-10);
}

TEST(primal_nurbscurve, curvature_derivative3d)
{
  checkCircularArcCurvatureDerivative3D<float>(1.e-5F);
  checkCircularArcCurvatureDerivative3D<double>(1.e-10);
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, curvature_reverse_orientation)
{
  constexpr int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  PointType controlPoints[3] = {PointType {0.0, 0.0}, PointType {0.5, 0.5}, PointType {1.0, 0.0}};
  NURBSCurveType curve(controlPoints, 3, 2);
  NURBSCurveType reversed(curve);
  reversed.reverseOrientation();

  for(const CoordType t : {0.0, 0.25, 0.5, 0.75, 1.0})
  {
    EXPECT_NEAR(curve.curvature(t), -reversed.curvature(1.0 - t), 1e-12);
    EXPECT_NEAR(curve.curvatureDerivative(t), reversed.curvatureDerivative(1.0 - t), 1e-12);
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, quartic_graph_derivatives_and_curvature)
{
  constexpr int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  // This degree-4 nonrational NURBS span exactly represents
  // f(x) = 0.25 * (x + 1) * (x - 2) * (x^2 - 0.25)
  // over x in [-1, 2], parameterized by x(t) = -1 + 3 t.
  auto analyticY = [](CoordType t) {
    return (81.0 / 4.0) * t * t * t * t - (135.0 / 4.0) * t * t * t + (243.0 / 16.0) * t * t -
      (27.0 / 16.0) * t;
  };

  auto analyticYp = [](CoordType t) {
    return 81.0 * t * t * t - (405.0 / 4.0) * t * t + (243.0 / 8.0) * t - (27.0 / 16.0);
  };

  auto analyticYpp = [](CoordType t) { return 243.0 * t * t - (405.0 / 2.0) * t + (243.0 / 8.0); };

  auto analyticYppp = [](CoordType t) { return 486.0 * t - (405.0 / 2.0); };

  auto analyticCurvature = [&](CoordType t) {
    const CoordType yp = analyticYp(t);
    const CoordType ypp = analyticYpp(t);
    const CoordType denom = 9.0 + yp * yp;
    return 3.0 * ypp / std::pow(denom, 1.5);
  };

  auto analyticCurvatureDerivative = [&](CoordType t) {
    const CoordType yp = analyticYp(t);
    const CoordType ypp = analyticYpp(t);
    const CoordType yppp = analyticYppp(t);
    const CoordType denom = 9.0 + yp * yp;
    return 3.0 * yppp / std::pow(denom, 1.5) - 9.0 * yp * ypp * ypp / std::pow(denom, 2.5);
  };

  PointType controlPoints[5] = {PointType {-1.0, 0.0},
                                PointType {-0.25, -27.0 / 64.0},
                                PointType {0.5, 27.0 / 16.0},
                                PointType {1.25, -135.0 / 64.0},
                                PointType {2.0, 0.0}};

  NURBSCurveType curve(controlPoints, 5, 4);

  for(const CoordType t : {0.0, 1.0 / 6.0, 1.0 / 3.0, 0.5, 2.0 / 3.0, 5.0 / 6.0, 1.0})
  {
    const PointType eval = curve.evaluate(t);
    const VectorType Dt = curve.dt(t);
    const VectorType DtDt = curve.dtdt(t);

    PointType evalBatch;
    axom::Array<VectorType> ders;
    curve.evaluateDerivatives(t, 3, evalBatch, ders);

    EXPECT_NEAR(eval[0], -1.0 + 3.0 * t, 1e-13);
    EXPECT_NEAR(eval[1], analyticY(t), 1e-13);
    EXPECT_NEAR(evalBatch[0], eval[0], 1e-13);
    EXPECT_NEAR(evalBatch[1], eval[1], 1e-13);

    EXPECT_NEAR(Dt[0], 3.0, 1e-13);
    EXPECT_NEAR(Dt[1], analyticYp(t), 1e-12);
    EXPECT_NEAR(DtDt[0], 0.0, 1e-13);
    EXPECT_NEAR(DtDt[1], analyticYpp(t), 1e-11);

    ASSERT_EQ(ders.size(), 3);
    EXPECT_NEAR(ders[0][0], 3.0, 1e-13);
    EXPECT_NEAR(ders[0][1], analyticYp(t), 1e-12);
    EXPECT_NEAR(ders[1][0], 0.0, 1e-13);
    EXPECT_NEAR(ders[1][1], analyticYpp(t), 1e-11);
    EXPECT_NEAR(ders[2][0], 0.0, 1e-13);
    EXPECT_NEAR(ders[2][1], analyticYppp(t), 1e-10);

    EXPECT_NEAR(curve.curvature(t), analyticCurvature(t), 1e-12);
    EXPECT_NEAR(curve.curvatureDerivative(t), analyticCurvatureDerivative(t), 1e-11);
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, repeated_knot_kink_behavior)
{
  constexpr int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  // This quadratic nonrational NURBS curve has an internal knot with
  // multiplicity equal to the degree, so it is only C^0 at t = 0.5.
  PointType controlPoints[5] = {PointType {0.0, 0.0},
                                PointType {1.0, 0.0},
                                PointType {1.0, 1.0},
                                PointType {2.0, 1.0},
                                PointType {3.0, 1.0}};
  CoordType knots[8] = {0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0};
  NURBSCurveType curve(controlPoints, 5, knots, 8);

  ASSERT_TRUE(curve.isValidNURBS());

  const PointType pStart = curve.evaluate(0.0);
  const PointType pEnd = curve.evaluate(1.0);
  const VectorType dtStart = curve.dt(0.0);
  const VectorType dtEnd = curve.dt(1.0);
  const VectorType dtdtStart = curve.dtdt(0.0);
  const VectorType dtdtEnd = curve.dtdt(1.0);

  EXPECT_NEAR(pStart[0], 0.0, 1e-14);
  EXPECT_NEAR(pStart[1], 0.0, 1e-14);
  EXPECT_NEAR(dtStart[0], 4.0, 1e-12);
  EXPECT_NEAR(dtStart[1], 0.0, 1e-12);
  EXPECT_NEAR(dtdtStart[0], -8.0, 1e-12);
  EXPECT_NEAR(dtdtStart[1], 8.0, 1e-12);
  EXPECT_NEAR(curve.curvature(0.0), 0.5, 1e-12);

  EXPECT_NEAR(pEnd[0], 3.0, 1e-14);
  EXPECT_NEAR(pEnd[1], 1.0, 1e-14);
  EXPECT_NEAR(dtEnd[0], 4.0, 1e-12);
  EXPECT_NEAR(dtEnd[1], 0.0, 1e-12);
  EXPECT_NEAR(dtdtEnd[0], 0.0, 1e-12);
  EXPECT_NEAR(dtdtEnd[1], 0.0, 1e-12);
  EXPECT_NEAR(curve.curvature(1.0), 0.0, 1e-12);
  EXPECT_NEAR(curve.curvatureDerivative(1.0), 0.0, 1e-12);

  const CoordType tLine = 0.75;
  EXPECT_NEAR(curve.curvature(tLine), 0.0, 1e-12);
  EXPECT_NEAR(curve.curvatureDerivative(tLine), 0.0, 1e-12);

  const CoordType eps = 1.e-4;
  const CoordType tLeft = 0.5 - eps;
  const CoordType tRight = 0.5 + eps;

  const PointType pMid = curve.evaluate(0.5);
  const PointType pLeft = curve.evaluate(tLeft);
  const PointType pRight = curve.evaluate(tRight);

  EXPECT_NEAR(pMid[0], 1.0, 1e-14);
  EXPECT_NEAR(pMid[1], 1.0, 1e-14);
  EXPECT_NEAR(pLeft[0], 1.0 - 4.0 * eps * eps, 1e-12);
  EXPECT_NEAR(pLeft[1], 1.0 - 4.0 * eps + 4.0 * eps * eps, 1e-12);
  EXPECT_NEAR(pRight[0], 1.0 + 4.0 * eps, 1e-12);
  EXPECT_NEAR(pRight[1], 1.0, 1e-12);

  const VectorType dtLeft = curve.dt(tLeft);
  const VectorType dtRight = curve.dt(tRight);
  const VectorType dtdtLeft = curve.dtdt(tLeft);
  const VectorType dtdtRight = curve.dtdt(tRight);

  EXPECT_NEAR(dtLeft[0], 8.0 * eps, 1e-10);
  EXPECT_NEAR(dtLeft[1], 4.0 - 8.0 * eps, 1e-10);
  EXPECT_NEAR(dtRight[0], 4.0, 1e-12);
  EXPECT_NEAR(dtRight[1], 0.0, 1e-12);

  EXPECT_NEAR(dtdtLeft[0], -8.0, 1e-10);
  EXPECT_NEAR(dtdtLeft[1], 8.0, 1e-10);
  EXPECT_NEAR(dtdtRight[0], 0.0, 1e-12);
  EXPECT_NEAR(dtdtRight[1], 0.0, 1e-12);

  // Position is continuous at the knot, but the sided derivatives are not.
  EXPECT_GT((dtLeft - dtRight).norm(), 1.0);
  EXPECT_GT((dtdtLeft - dtdtRight).norm(), 1.0);

  const CoordType kLeft = curve.curvature(tLeft);
  const CoordType kRight = curve.curvature(tRight);
  const CoordType expectedKLeft = 1.0 / (2.0 * std::pow(1.0 - 4.0 * eps + 8.0 * eps * eps, 1.5));
  EXPECT_NEAR(kLeft, expectedKLeft, 1e-9);
  EXPECT_NEAR(kRight, 0.0, 1e-12);
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, nurbscurve_intersections)
{
  constexpr int DIM = 2;
  using CoordType = double;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;
  using Point2D = primal::Point<CoordType, DIM>;

  constexpr int max_degree = 3;

  // Define two nurbs curves in 2D intersecting at one point.
  const Point2D data1_2d[max_degree + 1] = {Point2D {0.6, 1.2},
                                            Point2D {1.3, 1.6},
                                            Point2D {2.9, 2.4},
                                            Point2D {3.2, 3.5}};

  const Point2D data2_2d[max_degree + 1] = {Point2D {0.5, 3.4},
                                            Point2D {1.2, 2.3},
                                            Point2D {2.8, 1.5},
                                            Point2D {3.1, 1.1}};

  constexpr double weights[4] = {1.0, 2.0, 3.0, 4.0};

  constexpr int degree = 3;
  constexpr int npts = 4;

  NURBSCurveType curve1(data1_2d, weights, npts, degree);
  NURBSCurveType curve2(data2_2d, weights, npts, degree);

  Point2D intersection1, intersection2;

  axom::Array<CoordType> p1, p2, q1, q2;
  const bool found = intersect(curve1, curve2, p1, p2);
  const axom::IndexType num_intersections = p1.size();
  EXPECT_TRUE(found && num_intersections == 1 && num_intersections == p2.size());

  for(axom::IndexType j = 0; j < num_intersections; ++j)
  {
    intersection1 = curve1.evaluate(p1[j]);
    intersection2 = curve2.evaluate(p2[j]);

    for(int i = 0; i < DIM; ++i) EXPECT_NEAR(intersection1[i], intersection2[i], 1e-8);
  }

  // Test two curves that do not intersect.
  const Point2D data3_2d[max_degree + 1] = {Point2D {0.5, -3.4},
                                            Point2D {1.2, -2.3},
                                            Point2D {2.8, -1.5},
                                            Point2D {3.1, -1.1}};

  NURBSCurveType curve3(data3_2d, weights, npts, degree);
  const bool not_found = !intersect(curve2, curve3, q1, q2);
  EXPECT_TRUE(not_found && q1.size() == 0 && q2.size() == 0);
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, nurbscurve_circle_intersections)
{
  constexpr int DIM = 2;
  using CoordType = double;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  // Test two circles that intersect at two points.
  const auto circle1 = NURBSCurveType::make_circular_arc_nurbs(0.0, 2.0 * M_PI, 0.0, 0.0, 1.0);
  const auto circle2 = NURBSCurveType::make_circular_arc_nurbs(0.0, 2.0 * M_PI, 1.0, 0.0, 1.0);

  axom::Array<CoordType> p1, p2;
  const bool found = intersect(circle1, circle2, p1, p2);
  EXPECT_TRUE(found && p1.size() == 2 && p2.size() == 2);
}

primal::NURBSCurve<double, 2> make_cubic_shape()
{
  // Open cubic NURBS curve (degree 3), non-rational.
  // 24 control points, 28 knots (n + p + 2: 24 + 3 + 1 = 28).
  using Point2D = primal::Point<double, 2>;

  axom::Array<double> weights(24);
  for(int i = 0; i < 24; ++i) weights[i] = 1.0;

  return primal::NURBSCurve<double, 2>(
    axom::Array<Point2D> {
      Point2D {-5.0, +3.6}, Point2D {-3.4, +3.0}, Point2D {-1.8, +2.2}, Point2D {-0.4, +1.0},
      Point2D {+1.0, -0.4}, Point2D {+2.4, -1.6}, Point2D {+3.6, -2.4}, Point2D {+4.4, -1.0},
      Point2D {+3.6, +0.8}, Point2D {+2.0, +2.0}, Point2D {+0.0, +2.8}, Point2D {-2.0, +2.0},
      Point2D {-3.2, +0.8}, Point2D {-2.6, -0.8}, Point2D {-1.0, -1.6}, Point2D {+0.4, -2.2},
      Point2D {+1.6, -2.6}, Point2D {+0.0, -3.4}, Point2D {-2.2, -3.0}, Point2D {-2.8, -1.4},
      Point2D {-1.0, +0.4}, Point2D {+1.8, +1.8}, Point2D {+3.4, +3.0}, Point2D {+5.0, +3.6}},
    weights,
    axom::Array<double> {0.0,          0.0,          0.0,          0.0,          0.0476190476,
                         0.0952380952, 0.1428571429, 0.1904761905, 0.2380952381, 0.2857142857,
                         0.3333333333, 0.3809523810, 0.4285714286, 0.4761904762, 0.5238095238,
                         0.5714285714, 0.6190476190, 0.6666666667, 0.7142857143, 0.7619047619,
                         0.8095238095, 0.8571428571, 0.9047619048, 0.9523809524, 1.0,
                         1.0,          1.0,          1.0});
}

primal::NURBSCurve<double, 2> make_ellipse_curve()
{
  // Quadratic rational NURBS ellipse (degree 2), 9-control-point 4-arc construction.
  // Center, semi-axes:
  //   center      = (+0.055509, +0.246636)
  //   semi-axes   = (2.506091, 2.506234)   (a == b -> circle in this case)
  //   rotation    = +0.0000 degrees
  // Corner control-point weights are w = cos(pi/4) = sqrt(2)/2.
  using Point2D = primal::Point<double, 2>;
  const double w = 1.0 / std::sqrt(2.0);

  return primal::NURBSCurve<double, 2>(
    axom::Array<Point2D> {Point2D {+2.5615994286, +0.2466357797},
                          Point2D {+2.5615994286, +2.7528699841},
                          Point2D {+0.0555086484, +2.7528699841},
                          Point2D {-2.4505821318, +2.7528699841},
                          Point2D {-2.4505821318, +0.2466357797},
                          Point2D {-2.4505821318, -2.2595984246},
                          Point2D {+0.0555086484, -2.2595984246},
                          Point2D {+2.5615994286, -2.2595984246},
                          Point2D {+2.5615994286, +0.2466357797}},
    axom::Array<double> {1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0},
    axom::Array<double> {0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0});
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, nurbscurve_self_intersections)
{
  constexpr int DIM = 2;
  using CoordType = double;
  using Point2D = primal::Point<double, 2>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  NURBSCurveType curve1 = make_cubic_shape();
  NURBSCurveType curve2 = make_ellipse_curve();

  // Note: This pair of NURBS curves has eight intersections at five unique intersection points
  axom::Array<CoordType> p1, p2;
  const bool found = intersect(curve1, curve2, p1, p2);
  EXPECT_TRUE(found && p1.size() == 8 && p2.size() == 8);

  axom::Array<Point2D> intersections {Point2D {-1.7060766733448747, 2.0292202966044366},
                                      Point2D {2.044376271360025, -1.2782097206437852},
                                      Point2D {1.884074237582652, 1.9604430112255964},
                                      Point2D {-2.1644784828526533, -0.9161966686461217},
                                      Point2D {0.5039968961448462, -2.2191102081421588}};

  primal::Point<double, 2> int1, int2;
  std::array<int, 8> intid = {0, 1, 2, 0, 3, 4, 3, 2};

  for(int i = 0; i < p1.size(); ++i)
  {
    int1 = curve1.evaluate(p1[i]);
    int2 = curve2.evaluate(p2[i]);

    for(int j = 0; j < DIM; ++j)
    {
      EXPECT_NEAR(int1[j], int2[j], 1e-8);
      EXPECT_NEAR(int1[j], intersections[intid[i]][j], 1e-4);
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
