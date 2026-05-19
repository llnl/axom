// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file primal_bezier_triangle.cpp
 * \brief This file tests primal's Bezier triangle functionality
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/primal/geometry/BezierTriangle.hpp"

#include <array>
#include <iterator>

namespace primal = axom::primal;

//------------------------------------------------------------------------------
TEST(primal_beziertriangle, sizing_constructors)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using BezierTriangleType = primal::BezierTriangle<CoordType, DIM>;
  using CoordsVec = BezierTriangleType::CoordsVec;

  // testing default BezierTriangle constructor
  {
    BezierTriangleType bTri;
    EXPECT_FALSE(bTri.isRational());

    EXPECT_EQ(-1, bTri.getOrder());
    EXPECT_EQ(0, bTri.getControlPoints().size());
    EXPECT_EQ(CoordsVec(), bTri.getControlPoints());
  }

  // testing BezierTriangle order constructor
  for(int ord = -1; ord < 5; ++ord)
  {
    BezierTriangleType bTri(ord);
    EXPECT_FALSE(bTri.isRational());

    EXPECT_EQ(ord, bTri.getOrder());
    EXPECT_EQ(BezierTriangleType::triSize(ord), bTri.getControlPoints().size());
    EXPECT_EQ(0, bTri.getWeights().size());
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziertriangle, array_constructors)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierTriangleType = primal::BezierTriangle<CoordType, DIM>;

  SLIC_INFO("Testing BezierTriangle array constructors");

  constexpr int ord = 2;
  constexpr int npts = (ord + 1) * (ord + 2) / 2;

  PointType controlPoints[npts];
  CoordType weights[npts];

  for(int i = 0; i <= ord; ++i)
  {
    for(int j = 0; j <= ord - i; ++j)
    {
      const int idx = BezierTriangleType::triIndex(ord, i, j);
      controlPoints[idx] = PointType {static_cast<CoordType>(i),
                                      static_cast<CoordType>(j),
                                      static_cast<CoordType>(i + 2 * j)};
      weights[idx] = static_cast<CoordType>(0.25 + idx);
    }
  }

  auto check_triangle = [&](const BezierTriangleType& tri, bool expect_rational) {
    EXPECT_EQ(ord, tri.getOrder());
    EXPECT_EQ(npts, tri.getControlPoints().size());
    EXPECT_EQ(expect_rational, tri.isRational());
    EXPECT_EQ(expect_rational ? npts : 0, tri.getWeights().size());

    for(int i = 0; i <= ord; ++i)
    {
      for(int j = 0; j <= ord - i; ++j)
      {
        const int idx = BezierTriangleType::triIndex(ord, i, j);
        EXPECT_EQ(tri(i, j), controlPoints[idx]);
        if(expect_rational)
        {
          EXPECT_EQ(tri.getWeights()[idx], weights[idx]);
        }
      }
    }
  };

  // check C-array constructors
  {
    SCOPED_TRACE("Testing C-array constructor, polynomial");
    BezierTriangleType nTri(controlPoints, ord);
    check_triangle(nTri, false);

    SCOPED_TRACE("Testing C-array constructor, polynomial w/ null weight");
    BezierTriangleType nTri2(controlPoints, static_cast<const CoordType*>(nullptr), ord);
    check_triangle(nTri2, false);

    SCOPED_TRACE("Testing C-array constructor, rational");
    BezierTriangleType rTri(controlPoints, weights, ord);
    check_triangle(rTri, true);
  }

  // check ArrayView constructors (from C-arrays)
  {
    axom::ArrayView<PointType> cp(controlPoints, npts);
    axom::ArrayView<CoordType> w(weights, npts);

    SCOPED_TRACE("Testing ArrayView constructor, polynomial");
    BezierTriangleType nTri(cp, axom::ArrayView<CoordType>(nullptr, 0), ord);
    check_triangle(nTri, false);

    SCOPED_TRACE("Testing ArrayView constructor, rational");
    BezierTriangleType rTri(cp, w, ord);
    check_triangle(rTri, true);
  }

  // check 1D Array constructors
  {
    axom::Array<PointType> cp;
    cp.assign(std::begin(controlPoints), std::end(controlPoints));

    axom::Array<CoordType> w;
    w.assign(std::begin(weights), std::end(weights));

    SCOPED_TRACE("Testing Array constructor, polynomial");
    BezierTriangleType nTri(cp, ord);
    check_triangle(nTri, false);

    SCOPED_TRACE("Testing Array constructor, rational");
    BezierTriangleType rTri(cp, w, ord);
    check_triangle(rTri, true);
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziertriangle, set_order)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using BezierTriangleType = primal::BezierTriangle<CoordType, DIM>;

  BezierTriangleType bTri;
  EXPECT_EQ(-1, bTri.getOrder());
  EXPECT_EQ(0, bTri.getControlPoints().size());
  EXPECT_FALSE(bTri.isRational());

  constexpr int ord = 2;
  bTri.setOrder(ord);
  EXPECT_EQ(ord, bTri.getOrder());
  EXPECT_EQ(BezierTriangleType::triSize(ord), bTri.getControlPoints().size());
  EXPECT_FALSE(bTri.isRational());

  // Setting the order should resize weights in rational triangles
  BezierTriangleType rTri(1);
  rTri.getWeights().resize(rTri.getControlPoints().size());
  for(int i = 0; i < rTri.getWeights().size(); ++i)
  {
    rTri.getWeights()[i] = 1.0;
  }
  EXPECT_TRUE(rTri.isRational());

  constexpr int ord2 = 3;
  rTri.setOrder(ord2);
  EXPECT_EQ(ord2, rTri.getOrder());
  EXPECT_EQ(BezierTriangleType::triSize(ord2), rTri.getControlPoints().size());
  EXPECT_EQ(BezierTriangleType::triSize(ord2), rTri.getWeights().size());
}

//------------------------------------------------------------------------------
TEST(primal_beziertriangle, evaluate_linear)
{
  using CoordType = double;
  using BTri = primal::BezierTriangle<CoordType, 3>;
  using PointType = BTri::PointType;
  using VectorType = BTri::VectorType;

  constexpr int ord = 1;
  constexpr CoordType eps = 1e-14;

  BTri tri(ord);

  const PointType p00 {1.0, 2.0, 3.0};
  const PointType p01 {2.0, -1.0, 0.5};
  const PointType p10 {0.25, 1.5, -2.0};
  tri(0, 0) = p00;
  tri(0, 1) = p01;
  tri(1, 0) = p10;

  EXPECT_EQ(tri.evaluate(0.0, 0.0), p00);
  EXPECT_EQ(tri.evaluate(0.0, 1.0), p01);
  EXPECT_EQ(tri.evaluate(1.0, 0.0), p10);

  const CoordType u = 0.3;
  const CoordType v = 0.2;
  const PointType expected = PointType {p00[0] + u * (p10[0] - p00[0]) + v * (p01[0] - p00[0]),
                                        p00[1] + u * (p10[1] - p00[1]) + v * (p01[1] - p00[1]),
                                        p00[2] + u * (p10[2] - p00[2]) + v * (p01[2] - p00[2])};

  const PointType eval = tri.evaluate(u, v);
  for(int d = 0; d < 3; ++d)
  {
    EXPECT_NEAR(eval[d], expected[d], eps);
  }

  PointType e1;
  VectorType Du, Dv;
  tri.evaluateFirstDerivatives(u, v, e1, Du, Dv);
  for(int d = 0; d < 3; ++d)
  {
    EXPECT_NEAR(e1[d], expected[d], eps);
    EXPECT_NEAR(Du[d], (p10[d] - p00[d]), eps);
    EXPECT_NEAR(Dv[d], (p01[d] - p00[d]), eps);
  }

  VectorType DuDu, DvDv, DuDv;
  tri.evaluateSecondDerivatives(u, v, e1, Du, Dv, DuDu, DvDv, DuDv);
  for(int d = 0; d < 3; ++d)
  {
    EXPECT_NEAR(DuDu[d], 0.0, eps);
    EXPECT_NEAR(DvDv[d], 0.0, eps);
    EXPECT_NEAR(DuDv[d], 0.0, eps);
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziertriangle, evaluate_quadratic_second_derivatives_constant)
{
  using CoordType = double;
  using BTri = primal::BezierTriangle<CoordType, 3>;
  using PointType = BTri::PointType;
  using VectorType = BTri::VectorType;

  constexpr int ord = 2;
  constexpr CoordType eps = 1e-12;

  BTri tri(ord);

  // Populate control net with deterministic values
  for(int i = 0; i <= ord; ++i)
  {
    for(int j = 0; j <= ord - i; ++j)
    {
      tri(i, j) = PointType {static_cast<CoordType>(i),
                             static_cast<CoordType>(j),
                             static_cast<CoordType>(i * i + j)};
    }
  }

  // Vertex interpolation
  EXPECT_EQ(tri.evaluate(0.0, 0.0), tri(0, 0));
  EXPECT_EQ(tri.evaluate(0.0, 1.0), tri(0, ord));
  EXPECT_EQ(tri.evaluate(1.0, 0.0), tri(ord, 0));

  auto derivs = [&](CoordType u, CoordType v) {
    PointType e;
    VectorType Du, Dv, DuDu, DvDv, DuDv;
    tri.evaluateSecondDerivatives(u, v, e, Du, Dv, DuDu, DvDv, DuDv);
    return std::array<VectorType, 3> {DuDu, DvDv, DuDv};
  };

  const auto d0 = derivs(0.2, 0.3);
  const auto d1 = derivs(0.6, 0.1);
  for(int d = 0; d < 3; ++d)
  {
    EXPECT_NEAR(d0[0][d], d1[0][d], eps);
    EXPECT_NEAR(d0[1][d], d1[1][d], eps);
    EXPECT_NEAR(d0[2][d], d1[2][d], eps);
  }
}

TEST(primal_beziertriangle, rational_triangles)
{
  using CoordType = double;
  using BTri = primal::BezierTriangle<CoordType, 3>;
  using PointType = BTri::PointType;
  using VectorType = BTri::VectorType;

  constexpr int ord = 2;
  constexpr CoordType eps = 1e-12;

  BTri poly(ord);
  BTri rat(ord);
  rat.getWeights().resize(rat.getControlPoints().size());
  for(int i = 0; i < rat.getWeights().size(); ++i)
  {
    rat.getWeights()[i] = 1.0;
  }

  // Populate a shared control net
  poly(0, 0) = PointType {0.0, 0.0, 0.0};
  poly(0, 1) = PointType {1.0, 0.0, 0.0};
  poly(0, 2) = PointType {2.0, 0.0, 0.0};
  poly(1, 0) = PointType {0.0, 1.0, 0.0};
  poly(1, 1) = PointType {1.0, 1.0, 1.0};
  poly(2, 0) = PointType {0.0, 2.0, 0.0};

  rat(0, 0) = poly(0, 0);
  rat(0, 1) = poly(0, 1);
  rat(0, 2) = poly(0, 2);
  rat(1, 0) = poly(1, 0);
  rat(1, 1) = poly(1, 1);
  rat(2, 0) = poly(2, 0);

  auto check_at = [&](CoordType u, CoordType v) {
    // value
    const PointType p0 = poly.evaluate(u, v);
    const PointType p1 = rat.evaluate(u, v);
    for(int d = 0; d < 3; ++d)
    {
      EXPECT_NEAR(p0[d], p1[d], eps);
    }

    // first derivatives
    PointType e0, e1;
    VectorType du0, dv0, du1, dv1;
    poly.evaluateFirstDerivatives(u, v, e0, du0, dv0);
    rat.evaluateFirstDerivatives(u, v, e1, du1, dv1);
    for(int d = 0; d < 3; ++d)
    {
      EXPECT_NEAR(e0[d], e1[d], eps);
      EXPECT_NEAR(du0[d], du1[d], eps);
      EXPECT_NEAR(dv0[d], dv1[d], eps);
    }

    // second derivatives
    VectorType duu0, dvv0, duv0, duu1, dvv1, duv1;
    poly.evaluateSecondDerivatives(u, v, e0, du0, dv0, duu0, dvv0, duv0);
    rat.evaluateSecondDerivatives(u, v, e1, du1, dv1, duu1, dvv1, duv1);
    for(int d = 0; d < 3; ++d)
    {
      EXPECT_NEAR(e0[d], e1[d], eps);
      EXPECT_NEAR(duu0[d], duu1[d], eps);
      EXPECT_NEAR(dvv0[d], dvv1[d], eps);
      EXPECT_NEAR(duv0[d], duv1[d], eps);
    }
  };

  check_at(0.2, 0.3);
  check_at(0.6, 0.1);
}
