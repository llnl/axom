// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! 
 * \file primal_bezier_patch.cpp
 * \brief This file tests primal's Bezier patch functionality
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/operators/squared_distance.hpp"

namespace primal = axom::primal;

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, sizing_constructors)
{
  using CoordType = double;
  using BezierPatchType = primal::BezierPatch<CoordType, 3>;
  using CoordsMat = BezierPatchType::CoordsMat;
  using PointType = BezierPatchType::PointType;
  using PointArrayView = axom::ArrayView<const PointType, 2>;
  using WeightsArrayView = axom::ArrayView<const double, 2>;

  // testing default BezierPatch constructor
  {
    BezierPatchType bPatch;
    EXPECT_FALSE(bPatch.isRational());

    EXPECT_EQ(-1, bPatch.getOrder_u());
    EXPECT_EQ(-1, bPatch.getOrder_v());

    EXPECT_EQ(0, bPatch.getControlPoints().shape()[0]);
    EXPECT_EQ(0, bPatch.getControlPoints().shape()[1]);

    EXPECT_EQ(CoordsMat(), bPatch.getControlPoints());
  }

  // testing BezierPatch order constructor
  for(int u = 0; u < 5; ++u)
  {
    for(int v = 0; v < 5; ++v)
    {
      for(const BezierPatchType& patch :
          {BezierPatchType(u, v),
           BezierPatchType(PointArrayView(nullptr, {0, 0}), u, v),
           BezierPatchType(PointArrayView(nullptr, {0, 0}), WeightsArrayView(nullptr, {0, 0}), u, v)})
      {
        EXPECT_FALSE(patch.isRational());

        EXPECT_EQ(u, patch.getOrder_u());
        EXPECT_EQ(v, patch.getOrder_v());

        EXPECT_EQ(u + 1, patch.getControlPoints().shape()[0]);
        EXPECT_EQ(v + 1, patch.getControlPoints().shape()[1]);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, array_constructors)
{
  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;

  SLIC_INFO("Testing point array constructor");

  constexpr int order_u = 2;
  constexpr int order_v = 1;

  // clang-format off
  PointType controlPoints[6] = {
    PointType {0.0, 0.0, 1.0}, PointType {0.0, 1.0, 2.0}, PointType {0.0, 2.0, 3.0},
    PointType {1.0, 0.0, 4.0}, PointType {1.0, 1.0, 5.0}, PointType {1.0, 2.0, 6.0}};

  CoordType weights[6] = {0.009, 1.019, 2.029, 
                          3.109, 4.119, 5.129};
  // clang-format on

  auto check_patch = [=](const BezierPatchType& patch, bool expect_rational) {
    EXPECT_EQ(order_u, patch.getOrder_u());
    EXPECT_EQ(order_v, patch.getOrder_v());

    EXPECT_EQ(order_u + 1, patch.getControlPoints().shape()[0]);
    EXPECT_EQ(order_v + 1, patch.getControlPoints().shape()[1]);

    EXPECT_EQ(patch.isRational(), expect_rational);
    if(expect_rational)
    {
      EXPECT_EQ(order_u + 1, patch.getWeights().shape()[0]);
      EXPECT_EQ(order_v + 1, patch.getWeights().shape()[1]);
    }
    else
    {
      EXPECT_EQ(0, patch.getWeights().shape()[0]);
      EXPECT_EQ(0, patch.getWeights().shape()[1]);
    }

    int idx = 0;
    for(int u = 0; u <= patch.getOrder_u(); ++u)
    {
      for(int v = 0; v <= patch.getOrder_v(); ++v, ++idx)
      {
        EXPECT_EQ(patch(u, v), controlPoints[idx]);

        if(expect_rational)
        {
          EXPECT_EQ(patch.getWeight(u, v), weights[idx]);
        }
      }
    }
  };

  // check C-array constructors
  {
    SCOPED_TRACE("Testing C-array constructor, polynomial");
    BezierPatchType nPatch(controlPoints, order_u, order_v);
    check_patch(nPatch, false);

    SCOPED_TRACE("Testing C-array constructor, polynomial w/ null weight");
    BezierPatchType nPatch2(controlPoints, nullptr, order_u, order_v);
    check_patch(nPatch, false);

    SCOPED_TRACE("Testing C-array constructor, rational");
    BezierPatchType rPatch(controlPoints, weights, order_u, order_v);
    check_patch(rPatch, true);
  }

  // check ArrayView constructors (from C-arrays)
  {
    axom::ArrayView<PointType, 2> cp(controlPoints, {order_u + 1, order_v + 1});
    axom::ArrayView<double, 2> w(weights, {order_u + 1, order_v + 1});

    SCOPED_TRACE("Testing ArrayView constructor, polynomial");
    BezierPatchType nPatch(cp, order_u, order_v);
    check_patch(nPatch, false);

    SCOPED_TRACE("Testing ArrayView constructor, rational");
    BezierPatchType rPatch(cp, w, order_u, order_v);
    check_patch(rPatch, true);
  }

  // check 1D Array contructors
  {
    axom::Array<PointType> cp;
    cp.assign(std::begin(controlPoints), std::end(controlPoints));
    axom::Array<double> w;
    w.assign(std::begin(weights), std::end(weights));

    SCOPED_TRACE("Testing Array constructor, polynomial");
    BezierPatchType nPatch(cp, order_u, order_v);
    check_patch(nPatch, false);

    SCOPED_TRACE("Testing Array constructor, rational");
    BezierPatchType rPatch(cp, w, order_u, order_v);
    check_patch(rPatch, true);
  }

  {
    // check 2D array constructors
    axom::Array<PointType, 2> controlPoints_2D(order_u + 1, order_v + 1);
    axom::Array<double, 2> weights_2D(order_u + 1, order_v + 1);
    int idx = 0;
    for(int u = 0; u <= order_u; ++u)
    {
      for(int v = 0; v <= order_v; ++v, ++idx)
      {
        controlPoints_2D(u, v) = controlPoints[idx];
        weights_2D(u, v) = weights[idx];
      }
    }

    {
      SCOPED_TRACE("Testing 2D Array constructor, polynomial");
      BezierPatchType nPatch(controlPoints_2D, order_u, order_v);
      check_patch(nPatch, false);
    }
    {
      SCOPED_TRACE("Testing 2D Array constructor, rational");
      BezierPatchType rPatch(controlPoints_2D, weights_2D, order_u, order_v);
      check_patch(rPatch, true);
    }
    {
      SCOPED_TRACE("Testing 2D ArrayView constructor from Array, polynomial");
      BezierPatchType nPatch(controlPoints_2D.view(), order_u, order_v);
      check_patch(nPatch, false);
    }
    {
      SCOPED_TRACE("Testing 2D Array constructor, rational");
      BezierPatchType rPatch(controlPoints_2D.view(), weights_2D.view(), order_u, order_v);
      check_patch(rPatch, true);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, set_order)
{
  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;

  SLIC_INFO("Test adding control points to an empty Bezier patch");

  BezierPatchType bPatch;
  EXPECT_EQ(bPatch.getOrder_u(), -1);
  EXPECT_EQ(bPatch.getOrder_v(), -1);

  const int order_u = 1;
  const int order_v = 1;

  bPatch.setOrder(order_u, order_v);
  EXPECT_EQ(bPatch.getOrder_u(), order_u);
  EXPECT_EQ(bPatch.getOrder_v(), order_v);

  PointType controlPoints[4] = {PointType {0.0, 0.0, 1.0},
                                PointType {0.0, 1.0, 0.0},
                                PointType {1.0, 0.0, 0.0},
                                PointType {1.0, 1.0, -1.0}};

  bPatch(0, 0) = controlPoints[0];
  bPatch(0, 1) = controlPoints[1];
  bPatch(1, 0) = controlPoints[2];
  bPatch(1, 1) = controlPoints[3];

  int idx = 0;
  for(int p = 0; p <= bPatch.getOrder_u(); ++p)
  {
    for(int q = 0; q <= bPatch.getOrder_v(); ++q, ++idx)
    {
      EXPECT_EQ(controlPoints[idx], bPatch(p, q));
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, make_rational)
{
  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;

  const int order_u = 1;
  const int order_v = 1;

  SLIC_INFO("Test rational property");

  PointType controlPoints[4] = {PointType {0.0, 0.0, 1.0},
                                PointType {0.0, 1.0, 0.0},
                                PointType {1.0, 0.0, 0.0},
                                PointType {1.0, 1.0, -1.0}};

  BezierPatchType rPatch(controlPoints, order_u, order_v);
  EXPECT_FALSE(rPatch.isRational());
  EXPECT_EQ(rPatch.getWeights().size(), 0);

  rPatch.makeRational();
  EXPECT_TRUE(rPatch.isRational());
  EXPECT_EQ(rPatch.getWeights().shape(), rPatch.getControlPoints().shape());

  // makeRational should set all weights to 1
  for(int p = 0; p <= order_u; ++p)
  {
    for(int q = 0; q <= order_v; ++q)
    {
      EXPECT_DOUBLE_EQ(rPatch.getWeight(p, q), 1.0);
    }
  }

  // With all weights 1, the surface should be the same as if unweighted
  BezierPatchType bPatch(controlPoints, order_u, order_v);
  for(double u = 0; u <= 1.0; u += 0.1)
  {
    for(double v = 0; v <= 1.0; v += 0.1)
    {
      auto pt1 = rPatch.evaluate(u, v);
      auto pt2 = bPatch.evaluate(u, v);
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(pt1[i], pt2[i], 1e-10);
      }
    }
  }

  rPatch.makeNonrational();
  EXPECT_FALSE(rPatch.isRational());
  EXPECT_EQ(rPatch.getWeights().size(), 0);
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, evaluate)
{
  SLIC_INFO("Testing bezier Patch evaluation");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;

  const int order_u = 3;
  const int order_v = 2;

  // clang-format off
  PointType controlPoints[(order_u + 1) * (order_v + 1)] = {
    PointType {0, 0, 0}, PointType {0, 4,  0}, PointType {0, 8, -3},
    PointType {2, 0, 6}, PointType {2, 4,  0}, PointType {2, 8,  0},
    PointType {4, 0, 0}, PointType {4, 4,  0}, PointType {4, 8,  3},
    PointType {6, 0, 0}, PointType {6, 4, -3}, PointType {6, 8,  0}};
  // clang-format on

  BezierPatchType bPatch(controlPoints, order_u, order_v);

  // Evaluate the patch at each of the four corners, where the patch interpolates
  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0, 0)[i], controlPoints[0][i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0, 1)[i], controlPoints[2][i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(1, 0)[i], controlPoints[9][i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(1, 1)[i], controlPoints[11][i]);
  }

  // Evaluate the patch at some interior points
  PointType interior_1 {3.0, 4.0, 0.5625};
  PointType interior_2 {1.5, 6.0, -171.0 / 512.0};
  PointType interior_3 {4.5, 2.0, 39.0 / 512.0};

  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0.5, 0.5)[i], interior_1[i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0.25, 0.75)[i], interior_2[i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0.75, 0.25)[i], interior_3[i]);
  }

  // Evaluate calculation is asymmetric across axes. Swap them, and repeat the tests.
  bPatch.swapAxes();
  EXPECT_EQ(bPatch.getOrder_u(), order_v);
  EXPECT_EQ(bPatch.getOrder_v(), order_u);

  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0, 0)[i], controlPoints[0][i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0, 1)[i], controlPoints[9][i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(1, 0)[i], controlPoints[2][i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(1, 1)[i], controlPoints[11][i]);
  }

  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0.5, 0.5)[i], interior_1[i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0.75, 0.25)[i], interior_2[i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0.25, 0.75)[i], interior_3[i]);
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, isocurve)
{
  SLIC_INFO("Testing bezier Patch isocurve conventions");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;

  const int order_u = 3;
  const int order_v = 2;

  // clang-format off
  PointType controlPoints[(order_u + 1) * (order_v + 1)] = {
    PointType {0, 0, 0}, PointType {0, 4,  0}, PointType {0, 8, -3},
    PointType {2, 0, 6}, PointType {2, 4,  0}, PointType {2, 8,  0},
    PointType {4, 0, 0}, PointType {4, 4,  0}, PointType {4, 8,  3},
    PointType {6, 0, 0}, PointType {6, 4, -3}, PointType {6, 8,  0}};
  // clang-format on

  BezierPatchType bPatch(controlPoints, order_u, order_v);

  // isocurve_u should *fix* a value of u, returning a curve parameterized by v
  EXPECT_EQ(bPatch.isocurve_u(0.5).getOrder(), order_v);
  EXPECT_EQ(bPatch.isocurve(0.5, 0).getOrder(), order_v);

  // isocurve_v should *fix* a value of v, returning a curve parameterized by u
  EXPECT_EQ(bPatch.isocurve_v(0.5).getOrder(), order_u);
  EXPECT_EQ(bPatch.isocurve(0.5, 1).getOrder(), order_u);
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, evaluation_degenerate)
{
  SLIC_INFO("Testing Bezier patch evaluation with one degenerate axis");
  // Should reduce to a Bezier curve along the nonempty dimension

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;

  const int order = 3;
  PointType data[order + 1] = {PointType {0.6, 1.2, 1.0},
                               PointType {1.3, 1.6, 1.8},
                               PointType {2.9, 2.4, 2.3},
                               PointType {3.2, 3.5, 3.0}};

  BezierCurveType bCurve(data, order);
  BezierPatchType bPatch(data, order, 0);

  for(double t = 0; t <= 1; t += 0.1)
  {
    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(bCurve.evaluate(t)[i], bPatch.evaluate(t, 0)[i], 1e-10);
      EXPECT_NEAR(bCurve.evaluate(t)[i], bPatch.evaluate(t, 0.5)[i], 1e-10);
      EXPECT_NEAR(bCurve.evaluate(t)[i], bPatch.evaluate(t, 1.0)[i], 1e-10);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, tangent)
{
  SLIC_INFO("Testing bezier Patch tangent calculation");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;

  const int order_u = 3;
  const int order_v = 2;

  // clang-format off
  PointType controlPoints[(order_u + 1) * (order_v + 1)] = {
    PointType {0, 0, 0}, PointType {0, 4,  0}, PointType {0, 8, -3},
    PointType {2, 0, 6}, PointType {2, 4,  0}, PointType {2, 8,  0},
    PointType {4, 0, 0}, PointType {4, 4,  0}, PointType {4, 8,  3},
    PointType {6, 0, 0}, PointType {6, 4, -3}, PointType {6, 8,  0}};
  // clang-format on

  BezierPatchType bPatch(controlPoints, order_u, order_v);

  // Evaluate the tangent at the lower corfour corners,
  //  where the tangents are defined by the control nodes
  VectorType node00_u = VectorType(controlPoints[0], controlPoints[3]);
  VectorType node00_v = VectorType(controlPoints[0], controlPoints[1]);

  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_NEAR(order_u * node00_u[i], bPatch.du(0, 0)[i], 1e-10);
    EXPECT_NEAR(order_v * node00_v[i], bPatch.dv(0, 0)[i], 1e-10);
  }

  // Evaluate the tangent at some interior nodes
  VectorType interior1_u = {6.0, 0.0, 567.0 / 128.0};
  VectorType interior1_v = {0.0, 8.0, -159.0 / 64.0};

  VectorType interior2_u = {6.0, 0.0, -657.0 / 128.0};
  VectorType interior2_v = {0.0, 8.0, -123.0 / 64.0};

  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_NEAR(interior1_u[i], bPatch.du(0.25, 0.75)[i], 1e-10);
    EXPECT_NEAR(interior1_v[i], bPatch.dv(0.25, 0.75)[i], 1e-10);

    EXPECT_NEAR(interior2_u[i], bPatch.du(0.75, 0.25)[i], 1e-10);
    EXPECT_NEAR(interior2_v[i], bPatch.dv(0.75, 0.25)[i], 1e-10);
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, normal)
{
  SLIC_INFO("Testing bezier Patch normal calculation");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;

  const int order_u = 3;
  const int order_v = 2;

  // clang-format off
  PointType controlPoints[(order_u + 1) * (order_v + 1)] = {
    PointType {0, 0, 0}, PointType {0, 4,  0}, PointType {0, 8, -3},
    PointType {2, 0, 6}, PointType {2, 4,  0}, PointType {2, 8,  0},
    PointType {4, 0, 0}, PointType {4, 4,  0}, PointType {4, 8,  3},
    PointType {6, 0, 0}, PointType {6, 4, -3}, PointType {6, 8,  0}};
  // clang-format on

  BezierPatchType bPatch(controlPoints, order_u, order_v);

  // Test on some interior nodes
  VectorType interior_1 = {-279.0 / 16.0, 819.0 / 32.0, 48.0};  // (0.25, 0.25)
  VectorType interior_2 = {-567.0 / 16.0, 477.0 / 32.0, 48.0};  // (0.25, 0.75)
  VectorType interior_3 = {657.0 / 16.0, 369.0 / 32.0, 48.0};   // (0.75, 0.25)
  VectorType interior_4 = {369.0 / 16.0, -513.0 / 32.0, 48.0};  // (0.75, 0.75)

  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_NEAR(interior_1[i], bPatch.normal(0.25, 0.25)[i], 1e-10);
    EXPECT_NEAR(interior_2[i], bPatch.normal(0.25, 0.75)[i], 1e-10);
    EXPECT_NEAR(interior_3[i], bPatch.normal(0.75, 0.25)[i], 1e-10);
    EXPECT_NEAR(interior_4[i], bPatch.normal(0.75, 0.75)[i], 1e-10);
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, rational_evaluation)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;

  constexpr int max_order_u = 3;
  constexpr int order_v = 4;

  // clang-format off
  PointType controlPoints[(max_order_u + 1) * (order_v + 1)] = {
                  PointType {0, 0, 0}, PointType{0, 4,  0}, PointType{0, 8, -3}, PointType{0, 12, 1}, PointType{0, 16, 3},
                  PointType {2, 0, 6}, PointType{2, 4,  5}, PointType{2, 8,  0}, PointType{4, 12, 2}, PointType{2, 16, 2},
                  PointType {4, 0, 0}, PointType{4, 4,  5}, PointType{4, 8,  3}, PointType{2, 12, 3}, PointType{4, 16, 1},
                  PointType {6, 0, 0}, PointType{6, 4, -3}, PointType{6, 8,  0}, PointType{6, 12, 2}, PointType{6, 16, 0}};
  
  double weights[(max_order_u + 1) * (order_v + 1)] = {
                  1.0, 1.0, 1.0, 1.0, 1.0,
                  1.0, 2.0, 3.0, 2.0, 1.0,
                  1.0, 2.0, 3.0, 2.0, 1.0,
                  1.0, 1.0, 1.0, 1.0, 1.0};

  PointType exp_vals[4][3] = {{PointType {0.0, 0.0,   0.0}, PointType {       0.0,         6.4,      -0.8064}, PointType {0.0, 16.0, 3.0}},
                              {PointType {0.4, 0.0,   1.2}, PointType {1561./997.,  6736./997.,   1279./997.}, PointType {1.6, 16.0, 2.2}},
                              {PointType {0.8, 0.0,  1.92}, PointType {  218./91., 8104./1183.,  2579./1183.}, PointType {3.2, 16.0, 1.4}},
                              {PointType {1.2, 0.0, 2.304}, PointType {       3.0, 8104./1183., 11157./4732.}, PointType {4.8, 16.0, 0.6}}};
  // clang-format on

  double u0 = 0.2, u1 = 0.5, u2 = 0.8;
  double v0 = 0.0, v1 = 0.4, v2 = 1.0;

  // For each test order, construct a Bezier patch using the first order_u + 1
  //  rows of control points and weights
  for(int order_u = 0; order_u <= max_order_u; ++order_u)
  {
    BezierPatchType patch(controlPoints, weights, order_u, order_v);

    PointType calc_0 = patch.evaluate(u0, v0);
    PointType calc_1 = patch.evaluate(u1, v1);
    PointType calc_2 = patch.evaluate(u2, v2);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(calc_0[i], exp_vals[order_u][0][i], 1e-15);
      EXPECT_NEAR(calc_1[i], exp_vals[order_u][1][i], 1e-15);
      EXPECT_NEAR(calc_2[i], exp_vals[order_u][2][i], 1e-15);
    }

    patch.swapAxes();

    calc_0 = patch.evaluate(v0, u0);
    calc_1 = patch.evaluate(v1, u1);
    calc_2 = patch.evaluate(v2, u2);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(calc_0[i], exp_vals[order_u][0][i], 1e-15);
      EXPECT_NEAR(calc_1[i], exp_vals[order_u][1][i], 1e-15);
      EXPECT_NEAR(calc_2[i], exp_vals[order_u][2][i], 1e-15);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, rational_first_derivatives)
{
  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;

  const int max_order_u = 3;
  const int order_v = 4;

  // clang-format off
  PointType controlPoints[(max_order_u + 1) * (order_v + 1)] = {
                  PointType {0, 0, 0}, PointType{0, 4,  0}, PointType{0, 8, -3}, PointType{0, 12, 1}, PointType{0, 16, 3},
                  PointType {2, 0, 6}, PointType{2, 4,  5}, PointType{2, 8,  0}, PointType{4, 12, 2}, PointType{2, 16, 2},
                  PointType {4, 0, 0}, PointType{4, 4,  5}, PointType{4, 8,  3}, PointType{2, 12, 3}, PointType{4, 16, 1},
                  PointType {6, 0, 0}, PointType{6, 4, -3}, PointType{6, 8,  0}, PointType{6, 12, 2}, PointType{6, 16, 0}};
  
  double weights[(max_order_u + 1) * (order_v + 1)] = {
                  1.0, 1.0, 1.0, 1.0, 1.0,
                  1.0, 2.0, 3.0, 2.0, 1.0,
                  1.0, 2.0, 3.0, 2.0, 1.0,
                  1.0, 1.0, 1.0, 1.0, 1.0};
                                                                      
  VectorType exp_du[4][3] = {{VectorType {0.0, 0.0,  0.0}, VectorType {             0.0,              0.0,               0.0}, VectorType {0.0, 0.0,  0.0}},
                             {VectorType {2.0, 0.0,  6.0}, VectorType {1951250./994009.,  444000./994009.,  2603726./994009.}, VectorType {2.0, 0.0, -1.0}},
                             {VectorType {4.0, 0.0,  7.2}, VectorType { 301180./107653., 444000./1399489., 4232824./1399489.}, VectorType {4.0, 0.0, -2.0}},
                             {VectorType {6.0, 0.0, 5.76}, VectorType {        330./91.,              0.0,       2523./2366.}, VectorType {6.0, 0.0, -3.0}}};

  VectorType exp_dv[4][3] = {{VectorType {  0.0,  16.0,     0.0}, VectorType {             0.0,               16.0,              -0.064}, VectorType {   0.0,  16.0,     8.0}},
                             {VectorType { 1.28,  19.2,    2.24}, VectorType {1276850./994009.,  12622200./994009.,   -3524050./994009.}, VectorType {-14.08,  28.8,    2.24}},
                             {VectorType {2.048, 21.76,  3.9552}, VectorType {  64650./107653., 16488700./1399489.,  -4637275./1399489.}, VectorType { 4.608, 31.36,  -9.664}},
                             {VectorType {2.304, 23.68, 5.46432}, VectorType {             0.0, 16488700./1399489., -12970725./5597956.}, VectorType { 6.912, 23.68, -11.328}}};
  // clang-format on 

  double u0 = 0.2, u1 = 0.5, u2 = 0.8;
  double v0 = 0.0, v1 = 0.4, v2 = 1.0;

  // For each test order, construct a Bezier patch using the first order_u + 1
  //  rows of control points and weights
  for(int order_u = 0; order_u <= max_order_u; ++order_u)
  {
    BezierPatchType patch(controlPoints, weights, order_u, order_v);

    VectorType calc_du_0 = patch.du(u0, v0);
    VectorType calc_du_1 = patch.du(u1, v1);
    VectorType calc_du_2 = patch.du(u2, v2);

    VectorType calc_dv_0 = patch.dv(u0, v0);
    VectorType calc_dv_1 = patch.dv(u1, v1);
    VectorType calc_dv_2 = patch.dv(u2, v2);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(calc_du_0[i], exp_du[order_u][0][i], 1e-13);
      EXPECT_NEAR(calc_du_1[i], exp_du[order_u][1][i], 1e-13);
      EXPECT_NEAR(calc_du_2[i], exp_du[order_u][2][i], 1e-13);

      EXPECT_NEAR(calc_dv_0[i], exp_dv[order_u][0][i], 1e-13);
      EXPECT_NEAR(calc_dv_1[i], exp_dv[order_u][1][i], 1e-13);
      EXPECT_NEAR(calc_dv_2[i], exp_dv[order_u][2][i], 1e-13);
    }

    // Swap the axes and check that the derivatives are correct
    patch.swapAxes();

    calc_du_0 = patch.dv(v0, u0);
    calc_du_1 = patch.dv(v1, u1);
    calc_du_2 = patch.dv(v2, u2);

    calc_dv_0 = patch.du(v0, u0);
    calc_dv_1 = patch.du(v1, u1);
    calc_dv_2 = patch.du(v2, u2);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(calc_du_0[i], exp_du[order_u][0][i], 1e-13);
      EXPECT_NEAR(calc_du_1[i], exp_du[order_u][1][i], 1e-13);
      EXPECT_NEAR(calc_du_2[i], exp_du[order_u][2][i], 1e-13);

      EXPECT_NEAR(calc_dv_0[i], exp_dv[order_u][0][i], 1e-13);
      EXPECT_NEAR(calc_dv_1[i], exp_dv[order_u][1][i], 1e-13);
      EXPECT_NEAR(calc_dv_2[i], exp_dv[order_u][2][i], 1e-13);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, rational_second_derivatives)
{
  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;

  const int max_order_u = 3;
  const int order_v = 4;

  // clang-format off
  PointType controlPoints[(max_order_u + 1) * (order_v + 1)] = {
                  PointType {0, 0, 0}, PointType{0, 4,  0}, PointType{0, 8, -3}, PointType{0, 12, 1}, PointType{0, 16, 3},
                  PointType {2, 0, 6}, PointType{2, 4,  5}, PointType{2, 8,  0}, PointType{4, 12, 2}, PointType{2, 16, 2},
                  PointType {4, 0, 0}, PointType{4, 4,  5}, PointType{4, 8,  3}, PointType{2, 12, 3}, PointType{4, 16, 1},
                  PointType {6, 0, 0}, PointType{6, 4, -3}, PointType{6, 8,  0}, PointType{6, 12, 2}, PointType{6, 16, 0}};
  
  double weights[(max_order_u + 1) * (order_v + 1)] = {
                  1.0, 1.0, 1.0, 1.0, 1.0,
                  1.0, 2.0, 3.0, 2.0, 1.0,
                  1.0, 2.0, 3.0, 2.0, 1.0,
                  1.0, 1.0, 1.0, 1.0, 1.0};
                                                                    
  // Exact rational values get unwieldy here. For example, exp_dudu[1][1][0] = -2903460000/991026973 exactly.
  VectorType exp_dudu[4][3] = {{VectorType {0.0, 0.0,   0.0}, VectorType {           0.0,             0.0,            0.0}, VectorType {0.0, 0.0, 0.0}},
                               {VectorType {0.0, 0.0,   0.0}, VectorType {-2.92974871432, -0.666653903476, -3.90942365198}, VectorType {0.0, 0.0, 0.0}},
                               {VectorType {0.0, 0.0, -24.0}, VectorType {-2.45334507849, -1.03357131222,     -4.32849911}, VectorType {0.0, 0.0, 0.0}},
                               {VectorType {0.0, 0.0, -50.4}, VectorType {           0.0, -1.90355193931,  -13.2112292415}, VectorType {0.0, 0.0, 0.0}}};

  VectorType exp_dvdv[4][3] = {{VectorType {     0.0,      0.0,        -36.0}, VectorType {            0.0,            0.0,         23.52}, VectorType {     0.0,      0.0,      -24.0}},
                               {VectorType {  -2.048,   -11.52,      -65.984}, VectorType {-0.114397491782, -7.89784961786, 17.3438006414}, VectorType {-166.912,   107.52,    -48.064}},
                               {VectorType {-5.89824, -28.1088,   -93.470976}, VectorType {-0.786421266682, -8.69476156044, 10.1439670088}, VectorType {66.10944, 148.6848, -113.57952}},
                               {VectorType {-8.84736, -44.8512, -116.0229888}, VectorType {            0.0, -8.69476156044, 4.56797689827}, VectorType {54.19008,  44.8512,  -84.39552}}};

  VectorType exp_dudv[4][3] = {{VectorType { 0.0,  0.0,    0.0}, VectorType {           0.0,            0.0,            0.0}, VectorType {  0.0,   0.0,    0.0}},
                               {VectorType { 4.8, 16.0,    6.4}, VectorType {0.882014338474, -4.30534194956, -5.33680771976}, VectorType {-11.2,  16.0,  -10.4}},
                               {VectorType {5.12, 25.6, 12.544}, VectorType {-2.58010892971, -3.12014017951, 0.484714839042}, VectorType {49.28,   6.4, -31.04}},
                               {VectorType {0.96, 28.8, 19.872}, VectorType { -3.6032437554,            0.0,  5.97901269678}, VectorType {-2.88, -28.8,  -0.48}}};
  // clang-format on

  double u0 = 0.2, u1 = 0.5, u2 = 0.8;
  double v0 = 0.0, v1 = 0.4, v2 = 1.0;

  // For each test order, construct a Bezier patch using the first order_u + 1
  //  rows of control points and weights
  for(int order_u = 0; order_u <= max_order_u; ++order_u)
  {
    BezierPatchType patch(controlPoints, weights, order_u, order_v);

    VectorType calc_dudu_0 = patch.dudu(u0, v0);
    VectorType calc_dudu_1 = patch.dudu(u1, v1);
    VectorType calc_dudu_2 = patch.dudu(u2, v2);

    VectorType calc_dvdv_0 = patch.dvdv(u0, v0);
    VectorType calc_dvdv_1 = patch.dvdv(u1, v1);
    VectorType calc_dvdv_2 = patch.dvdv(u2, v2);

    VectorType calc_dudv_0 = patch.dudv(u0, v0);
    VectorType calc_dudv_1 = patch.dudv(u1, v1);
    VectorType calc_dudv_2 = patch.dudv(u2, v2);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(calc_dudu_0[i], exp_dudu[order_u][0][i], 1e-10);
      EXPECT_NEAR(calc_dudu_1[i], exp_dudu[order_u][1][i], 1e-10);
      EXPECT_NEAR(calc_dudu_2[i], exp_dudu[order_u][2][i], 1e-10);

      EXPECT_NEAR(calc_dvdv_0[i], exp_dvdv[order_u][0][i], 1e-10);
      EXPECT_NEAR(calc_dvdv_1[i], exp_dvdv[order_u][1][i], 1e-10);
      EXPECT_NEAR(calc_dvdv_2[i], exp_dvdv[order_u][2][i], 1e-10);

      EXPECT_NEAR(calc_dudv_0[i], exp_dudv[order_u][0][i], 1e-10);
      EXPECT_NEAR(calc_dudv_1[i], exp_dudv[order_u][1][i], 1e-10);
      EXPECT_NEAR(calc_dudv_2[i], exp_dudv[order_u][2][i], 1e-10);
    }

    // Swap the axes and check that the derivatives are correct
    patch.swapAxes();

    calc_dudu_0 = patch.dvdv(v0, u0);
    calc_dudu_1 = patch.dvdv(v1, u1);
    calc_dudu_2 = patch.dvdv(v2, u2);

    calc_dvdv_0 = patch.dudu(v0, u0);
    calc_dvdv_1 = patch.dudu(v1, u1);
    calc_dvdv_2 = patch.dudu(v2, u2);

    calc_dudv_0 = patch.dudv(v0, u0);
    calc_dudv_1 = patch.dudv(v1, u1);
    calc_dudv_2 = patch.dudv(v2, u2);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(calc_dudu_0[i], exp_dudu[order_u][0][i], 1e-10);
      EXPECT_NEAR(calc_dudu_1[i], exp_dudu[order_u][1][i], 1e-10);
      EXPECT_NEAR(calc_dudu_2[i], exp_dudu[order_u][2][i], 1e-10);

      EXPECT_NEAR(calc_dvdv_0[i], exp_dvdv[order_u][0][i], 1e-10);
      EXPECT_NEAR(calc_dvdv_1[i], exp_dvdv[order_u][1][i], 1e-10);
      EXPECT_NEAR(calc_dvdv_2[i], exp_dvdv[order_u][2][i], 1e-10);

      EXPECT_NEAR(calc_dudv_0[i], exp_dudv[order_u][0][i], 1e-10);
      EXPECT_NEAR(calc_dudv_1[i], exp_dudv[order_u][1][i], 1e-10);
      EXPECT_NEAR(calc_dudv_2[i], exp_dudv[order_u][2][i], 1e-10);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, rational_batch_derivatives)
{
  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;

  const int max_order_u = 3;
  const int order_v = 4;

  // clang-format off
  PointType controlPoints[(max_order_u + 1) * (order_v + 1)] = {
                  PointType {0, 0, 0}, PointType{0, 4,  0}, PointType{0, 8, -3}, PointType{0, 12, 1}, PointType{0, 16, 3},
                  PointType {2, 0, 6}, PointType{2, 4,  5}, PointType{2, 8,  0}, PointType{4, 12, 2}, PointType{2, 16, 2},
                  PointType {4, 0, 0}, PointType{4, 4,  5}, PointType{4, 8,  3}, PointType{2, 12, 3}, PointType{4, 16, 1},
                  PointType {6, 0, 0}, PointType{6, 4, -3}, PointType{6, 8,  0}, PointType{6, 12, 2}, PointType{6, 16, 0}};
  
  double weights[(max_order_u + 1) * (order_v + 1)] = {
                  1.0, 1.0, 1.0, 1.0, 1.0,
                  1.0, 2.0, 3.0, 2.0, 1.0,
                  1.0, 2.0, 3.0, 2.0, 1.0,
                  1.0, 1.0, 1.0, 1.0, 1.0};
  // clang-format on

  PointType batch1_val, batch2_val, batch3_val;
  VectorType batch1_du, batch1_dv;
  VectorType batch2_du, batch2_dv, batch2_dudv;
  VectorType batch3_du, batch3_dv, batch3_dudu, batch3_dvdv, batch3_dudv;

  const double u = 0.4, v = 0.6;

  // For each test order, construct a Bezier patch using the first order_u + 1
  //  rows of control points and weights
  for(int order_u = 0; order_u <= max_order_u; ++order_u)
  {
    BezierPatchType patch(controlPoints, weights, order_u, order_v);

    patch.evaluateFirstDerivatives(u, v, batch1_val, batch1_du, batch1_dv);
    patch.evaluate_linear_derivatives(u, v, batch2_val, batch2_du, batch2_dv, batch2_dudv);
    patch.evaluateSecondDerivatives(u,
                                    v,
                                    batch3_val,
                                    batch3_du,
                                    batch3_dv,
                                    batch3_dudu,
                                    batch3_dvdv,
                                    batch3_dudv);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(patch.evaluate(u, v)[i], batch1_val[i], 1e-13);
      EXPECT_NEAR(patch.du(u, v)[i], batch1_du[i], 1e-13);
      EXPECT_NEAR(patch.dv(u, v)[i], batch1_dv[i], 1e-13);

      EXPECT_NEAR(patch.evaluate(u, v)[i], batch2_val[i], 1e-13);
      EXPECT_NEAR(patch.du(u, v)[i], batch2_du[i], 1e-13);
      EXPECT_NEAR(patch.dv(u, v)[i], batch2_dv[i], 1e-13);
      EXPECT_NEAR(patch.dudv(u, v)[i], batch2_dudv[i], 1e-13);

      EXPECT_NEAR(patch.evaluate(u, v)[i], batch3_val[i], 1e-13);
      EXPECT_NEAR(patch.du(u, v)[i], batch2_du[i], 1e-13);
      EXPECT_NEAR(patch.dv(u, v)[i], batch2_dv[i], 1e-13);
      EXPECT_NEAR(patch.dudu(u, v)[i], batch3_dudu[i], 1e-13);
      EXPECT_NEAR(patch.dvdv(u, v)[i], batch3_dvdv[i], 1e-13);
      EXPECT_NEAR(patch.dudv(u, v)[i], batch3_dudv[i], 1e-13);
    }

    // Swap the axes and check that the derivatives are correct
    patch.swapAxes();

    patch.evaluateFirstDerivatives(u, v, batch1_val, batch1_du, batch1_dv);
    patch.evaluate_linear_derivatives(u, v, batch2_val, batch2_du, batch2_dv, batch2_dudv);
    patch.evaluateSecondDerivatives(u,
                                    v,
                                    batch3_val,
                                    batch3_du,
                                    batch3_dv,
                                    batch3_dudu,
                                    batch3_dvdv,
                                    batch3_dudv);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(patch.evaluate(u, v)[i], batch1_val[i], 1e-13);
      EXPECT_NEAR(patch.du(u, v)[i], batch1_du[i], 1e-13);
      EXPECT_NEAR(patch.dv(u, v)[i], batch1_dv[i], 1e-13);

      EXPECT_NEAR(patch.evaluate(u, v)[i], batch2_val[i], 1e-13);
      EXPECT_NEAR(patch.du(u, v)[i], batch2_du[i], 1e-13);
      EXPECT_NEAR(patch.dv(u, v)[i], batch2_dv[i], 1e-13);
      EXPECT_NEAR(patch.dudv(u, v)[i], batch2_dudv[i], 1e-13);

      EXPECT_NEAR(patch.evaluate(u, v)[i], batch3_val[i], 1e-13);
      EXPECT_NEAR(patch.du(u, v)[i], batch3_du[i], 1e-13);
      EXPECT_NEAR(patch.dv(u, v)[i], batch3_dv[i], 1e-13);
      EXPECT_NEAR(patch.dudu(u, v)[i], batch3_dudu[i], 1e-13);
      EXPECT_NEAR(patch.dvdv(u, v)[i], batch3_dvdv[i], 1e-13);
      EXPECT_NEAR(patch.dudv(u, v)[i], batch3_dudv[i], 1e-13);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, split_degenerate)
{
  SLIC_INFO("Testing bezier Patch splitting for order 0 surface");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;

  const int order_u = 0;
  const int order_v = 0;

  PointType data[(order_u + 1) * (order_v + 1)] = {PointType::make_point(0.6, 1.2, 0.2)};

  BezierPatchType p(data, 0, 0);

  BezierPatchType p1, p2, p3, p4;

  p.split(0.5, 0.5, p1, p2, p3, p4);

  for(int i = 0; i < DIM; i++)
  {
    EXPECT_DOUBLE_EQ(data[0][i], p1(0, 0)[i]);
    EXPECT_DOUBLE_EQ(data[0][i], p2(0, 0)[i]);
    EXPECT_DOUBLE_EQ(data[0][i], p3(0, 0)[i]);
    EXPECT_DOUBLE_EQ(data[0][i], p4(0, 0)[i]);
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, split_curve)
{
  SLIC_INFO("Testing Bezier patch split with one degenerate axis");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;

  const int order = 3;
  PointType data[order + 1] = {PointType {0.6, 1.2, 1.0},
                               PointType {1.3, 1.6, 1.8},
                               PointType {2.9, 2.4, 2.3},
                               PointType {3.2, 3.5, 3.0}};

  BezierCurveType bCurve(data, order);
  BezierPatchType bPatch(data, order, 0);

  BezierCurveType c1, c2;
  BezierPatchType p1, p2, p3, p4;

  bCurve.split(0.5, c1, c2);
  bPatch.split(0.5, 0.5, p1, p2, p3, p4);

  // Components split along the u-axis p1/p3 vs p2/p4 should equal curves
  // Components split along the v-axis p1/p2 vs p3/p4 should equal each other
  for(double u = 0; u <= 1; u += 0.1)
  {
    for(double v = 0; v <= 1; v += 0.1)
    {
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(c1.evaluate(u)[i], p1.evaluate(u, v)[i], 1e-10);
        EXPECT_NEAR(c1.evaluate(u)[i], p3.evaluate(u, v)[i], 1e-10);

        EXPECT_NEAR(c2.evaluate(u)[i], p2.evaluate(u, v)[i], 1e-10);
        EXPECT_NEAR(c2.evaluate(u)[i], p4.evaluate(u, v)[i], 1e-10);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, split_plane)
{
  SLIC_INFO("Testing Bezier patch evaluation with one degenerate axis");
  // Should reduce to a Bezier curve along the nonempty dimension

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;

  const int order_u = 1;
  const int order_v = 1;

  // clang-format off
  PointType controlPoints[(order_u + 1) * (order_v + 1)] = {
    PointType {1.0, 0.0, 0.0}, PointType {0.0, 1.0, 0.0},
    PointType {0.0, 0.0, 1.0}, PointType {0.0, 0.0, 1.0}};
  // clang-format on

  BezierPatchType p0(controlPoints, order_u, order_v), p_arr[4];

  p0.split(0.5, 0.5, p_arr[0], p_arr[1], p_arr[2], p_arr[3]);

  // Ensure that after splitting, each point remains on the same hyperplane
  for(double u = 0.0; u <= 1.0; u += 0.1)
  {
    for(double v = 0.0; v <= 1.0; v += 0.1)
    {
      for(int n = 0; n < 4; ++n)
      {
        auto pt = p_arr[n].evaluate(u, v);
        EXPECT_NEAR(pt[0] + pt[1] + pt[2], 1.0, 1e-10);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, split_patch)
{
  SLIC_INFO("Testing Bezier patch evaluation with one degenerate axis");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;

  const int order_u = 3;
  const int order_v = 2;

  // clang-format off
  PointType controlPoints[(order_u + 1) * (order_v + 1)] = {
    PointType {0, 0, 0}, PointType {0, 4,  0}, PointType {0, 8, -3},
    PointType {2, 0, 6}, PointType {2, 4,  0}, PointType {2, 8,  0},
    PointType {4, 0, 0}, PointType {4, 4,  0}, PointType {4, 8,  3},
    PointType {6, 0, 0}, PointType {6, 4, -3}, PointType {6, 8,  0}};
  // clang-format on

  BezierPatchType p0(controlPoints, order_u, order_v), p1, p2, p3, p4;
  p0.split(0.5, 0.5, p1, p2, p3, p4);

  // Check that evaluation the appropriate edges line up
  for(double t = 0; t <= 1; t += 0.1)
  {
    for(int i = 0; i < DIM; i++)
    {
      EXPECT_NEAR(p1.evaluate(1, t)[i], p2.evaluate(0, t)[i], 1e-10);
      EXPECT_NEAR(p3.evaluate(1, t)[i], p4.evaluate(0, t)[i], 1e-10);

      EXPECT_NEAR(p1.evaluate(t, 1)[i], p3.evaluate(t, 0)[i], 1e-10);
      EXPECT_NEAR(p2.evaluate(t, 1)[i], p4.evaluate(t, 0)[i], 1e-10);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, isPlanar)
{
  SLIC_INFO("Testing Bezier patch evaluation with one degenerate axis");
  // Should reduce to a Bezier curve along the nonempty dimension

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;

  // order (0, 0) -- always true
  {
    const int order_u = 0;
    const int order_v = 0;
    BezierPatchType patch(order_u, order_v);
    EXPECT_TRUE(patch.isPlanar());

    patch(0, 0) = PointType::make_point(1.0, 1.0, 1.0);
    EXPECT_TRUE(patch.isPlanar());
  }

  // order (0, 1), (1, 0) -- always true
  {
    PointType data[2] = {PointType {0.0, 1.0, 2.0}, PointType {2.0, 1.0, 0.0}};
    BezierPatchType patch1(data, 0, 1), patch2(data, 1, 0);

    EXPECT_TRUE(patch1.isPlanar());
    EXPECT_TRUE(patch2.isPlanar());
  }

  // order (1, 1) -- sometimes true
  {
    // clang-format off
    PointType data_planar[4] = {
      PointType {1.0, 0.0, 0.0}, PointType {0.0, 1.0, 0.0},
      PointType {0.0, 0.0, 1.0}, PointType {-1.0, 1.0, 1.0}};

    PointType data_nonplanar[4] = {
      PointType {1.0, 0.0, 0.0}, PointType {0.0, 1.0, 0.0},
      PointType {0.0, 0.0, 1.0}, PointType {0.0, 0.0, 2.0}};
    // clang-format on

    EXPECT_TRUE(BezierPatchType(data_planar, 1, 1).isPlanar());
    EXPECT_FALSE(BezierPatchType(data_nonplanar, 1, 1).isPlanar());
  }

  // order (2, 2)
  {
    // clang-format off
    PointType data[9] = {
      PointType {-1.0, -1.0, 0.0}, PointType {-1.0, 0.0, 0.0}, PointType {-1.0, 1.0, 0.0},  
      PointType { 0.0, -1.0, 0.0}, PointType { 0.0, 0.0, 0.0}, PointType { 0.0, 1.0, 0.0},
      PointType { 1.0, -1.0, 0.0}, PointType { 1.0, 0.0, 0.0}, PointType { 1.0, 1.0, 0.0}};
    // clang-format on

    const int ord_u = 2;
    const int ord_v = 2;

    BezierPatchType patch(data, ord_u, ord_v);

    // Begins as planar
    EXPECT_TRUE(patch.isPlanar());

    // Move middle point and check linearity with diffrent tolerances
    patch(1, 1).array() += 0.005 * VectorType({0.0, 0.0, 1.0}).array();

    // Linear for a coarse tolerance
    {
      const double tol = 0.01;
      const double tol_sq = tol * tol;
      EXPECT_TRUE(patch.isPlanar(tol_sq));
    }

    // Non-Linear for a finer tolerance
    {
      const double tol = 0.001;
      const double tol_sq = tol * tol;
      EXPECT_FALSE(patch.isPlanar(tol_sq));
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, reverse_orientation)
{
  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;

  SLIC_INFO("Testing reverseOrientation(axis) on Bezier patches");

  const int order_u = 3;
  const int order_v = 2;

  // clang-format off
  PointType controlPoints[(order_u + 1) * (order_v + 1)] = {
    PointType {0, 0, 0}, PointType {0, 4,  0}, PointType {0, 8, -3},
    PointType {2, 0, 6}, PointType {2, 4,  0}, PointType {2, 8,  0},
    PointType {4, 0, 0}, PointType {4, 4,  0}, PointType {4, 8,  3},
    PointType {6, 0, 0}, PointType {6, 4, -3}, PointType {6, 8,  0}};
  // clang-format on

  BezierPatchType original(controlPoints, order_u, order_v);
  BezierPatchType reversed(controlPoints, order_u, order_v);

  // Reverse along the u-axis
  reversed.reverseOrientation(0);
  for(double u = 0; u <= 1.0; u += 0.1)
  {
    for(double v = 0; v <= 1.0; v += 0.1)
    {
      PointType o_pt = original.evaluate(u, v);
      PointType r_pt = reversed.evaluate(1 - u, v);

      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(o_pt[i], r_pt[i], 1e-10);
      }
    }
  }

  // Reverse along the u-axis again, should return to original
  reversed.reverseOrientation(0);
  for(double u = 0; u <= 1.0; u += 0.1)
  {
    for(double v = 0; v <= 1.0; v += 0.1)
    {
      PointType o_pt = original.evaluate(u, v);
      PointType r_pt = reversed.evaluate(u, v);

      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(o_pt[i], r_pt[i], 1e-10);
      }
    }
  }

  // Reverse along the v-axis
  reversed.reverseOrientation(1);
  for(double u = 0; u <= 1.0; u += 0.1)
  {
    for(double v = 0; v <= 1.0; v += 0.1)
    {
      PointType o_pt = original.evaluate(u, v);
      PointType r_pt = reversed.evaluate(u, 1 - v);

      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(o_pt[i], r_pt[i], 1e-10);
      }
    }
  }

  // Reverse along the u-axis again
  reversed.reverseOrientation(0);
  for(double u = 0; u <= 1.0; u += 0.1)
  {
    for(double v = 0; v <= 1.0; v += 0.1)
    {
      PointType o_pt = original.evaluate(u, v);
      PointType r_pt = reversed.evaluate(1 - u, 1 - v);

      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(o_pt[i], r_pt[i], 1e-10);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, rational_evaluation_split)
{
  // Test behavior with a spherical Bezier patch
  SLIC_INFO("Testing Bezier patch evaluation with rational weights");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;

  const int ord_u = 3;
  const int ord_v = 3;

  // clang-format off
  PointType data[16] = {  
    PointType {0, 0,  1}, PointType {0, 0,  1}, PointType { 0, 0,  1}, PointType { 0, 0,  1},
    PointType {2, 0,  1}, PointType {2, 4,  1}, PointType {-2, 4,  1}, PointType {-2, 0,  1},
    PointType {2, 0, -1}, PointType {2, 4, -1}, PointType {-2, 4, -1}, PointType {-2, 0, -1},
    PointType {0, 0, -1}, PointType {0, 0, -1}, PointType { 0, 0, -1}, PointType { 0, 0, -1}};

  double weights[16] = {1.0,     1.0/3.0, 1.0/3.0, 1.0,
                        1.0/3.0, 1.0/9.0, 1.0/9.0, 1.0/3.0,
                        1.0/3.0, 1.0/9.0, 1.0/9.0, 1.0/3.0,
                        1.0,     1.0/3.0, 1.0/3.0, 1.0};
  // clang-format on

  BezierPatchType hemisphere(data, weights, ord_u, ord_v);

  // Verify that evaluation points are on the sphere
  for(double u = 0; u <= 1.0; u += 0.1)
  {
    for(double v = 0; v <= 1.0; v += 0.1)
    {
      PointType pt = hemisphere.evaluate(u, v);
      EXPECT_NEAR(pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2], 1.0, 1e-10);
    }
  }

  BezierPatchType patches[4];
  hemisphere.split(0.5, 0.5, patches[0], patches[1], patches[2], patches[3]);

  // Verify that evaluation points are on still on the sphere for each subpatch
  for(double u = 0; u <= 1.0; u += 0.1)
  {
    for(double v = 0; v <= 1.0; v += 0.1)
    {
      for(int n = 0; n < 4; ++n)
      {
        PointType pt = patches[n].evaluate(u, v);
        EXPECT_NEAR(pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2], 1.0, 1e-10);
      }
    }
  }

  // Do it again
  BezierPatchType sub_patches[4];
  patches[0].split(0.5, 0.5, sub_patches[0], sub_patches[1], sub_patches[2], sub_patches[3]);

  for(double u = 0; u <= 1.0; u += 0.1)
  {
    for(double v = 0; v <= 1.0; v += 0.1)
    {
      for(int n = 0; n < 4; ++n)
      {
        PointType pt = sub_patches[n].evaluate(u, v);
        EXPECT_NEAR(pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2], 1.0, 1e-10);
      }
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
