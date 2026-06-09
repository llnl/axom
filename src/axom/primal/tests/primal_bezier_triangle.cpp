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
#include "axom/primal/geometry/GregoryTriangle.hpp"
#include "axom/primal/geometry/GregoryPatch.hpp"

#include <array>

namespace primal = axom::primal;

namespace
{
template <typename BTri>
void fillControlNet(BTri& tri)
{
  // Fills the BezierTriangle with sample data
  using CoordType = typename BTri::PointType::CoordType;
  using PointType = typename BTri::PointType;

  const int ord = tri.getOrder();
  for(int i = 0; i <= ord; ++i)
  {
    for(int j = 0; j <= ord - i; ++j)
    {
      const auto ii = static_cast<CoordType>(i);
      const auto jj = static_cast<CoordType>(j);
      tri(i, j) = PointType {ii, jj, 10.0 * ii - 3.0 * jj};
    }
  }
}

template <typename BTri>
void makeRational(BTri& tri)
{
  using CoordType = typename BTri::PointType::CoordType;

  tri.getWeights().resize(tri.getControlPoints().size());
  for(int k = 0; k < tri.getWeights().size(); ++k)
  {
    // Positive, non-uniform weights
    tri.getWeights()[k] = static_cast<CoordType>(0.5 + 0.25 * k);
  }
  EXPECT_TRUE(tri.isRational());
}
}  // namespace

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
TEST(primal_beziertriangle, parameter_convention_matches_barycentric)
{
  using CoordType = double;
  using BTri = primal::BezierTriangle<CoordType, 3>;
  using PointType = BTri::PointType;
  using Barycentric = BTri::Barycentric;

  constexpr CoordType eps = 1e-14;

  BTri tri(1);
  const PointType pA {1.0, 2.0, 3.0};     // (i,j) = (0,0)
  const PointType pB {4.0, -1.0, 0.5};    // (i,j) = (0,1)
  const PointType pC {-2.0, 0.25, 7.0};   // (i,j) = (1,0)

  tri(0, 0) = pA;
  tri(0, 1) = pB;
  tri(1, 0) = pC;

  // Corner mapping implied by the (u0,v0)->barycentric convention
  EXPECT_EQ(tri.evaluate(0.0, 0.0), pA);  // bary = (1,0,0)
  EXPECT_EQ(tri.evaluate(0.0, 1.0), pB);  // bary = (0,1,0)
  EXPECT_EQ(tri.evaluate(1.0, 0.0), pC);  // bary = (0,0,1)

  // Interior point check: evaluate(u0,v0) == lambda0*A + lambda1*B + lambda2*C,
  // with (lambda0,lambda1,lambda2) = (1-u0-v0, v0, u0)
  const CoordType u0 = 0.3;
  const CoordType v0 = 0.2;
  const Barycentric bary {1.0 - u0 - v0, v0, u0};

  const PointType expected = BTri::triInterpolate(pA, pB, pC, bary);
  const PointType eval = tri.evaluate(u0, v0);
  for(int d = 0; d < 3; ++d)
  {
    EXPECT_NEAR(eval[d], expected[d], eps);
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziertriangle, finite_difference_first_derivatives)
{
  using CoordType = double;
  using BTri = primal::BezierTriangle<CoordType, 3>;
  using PointType = BTri::PointType;
  using VectorType = BTri::VectorType;

  constexpr CoordType h = 1e-7;
  constexpr CoordType tol = 5e-6;

  BTri tri(4);
  fillControlNet(tri);

  const CoordType u0 = 0.23;
  const CoordType v0 = 0.31;

  PointType eval;
  VectorType Du, Dv;
  tri.evaluateFirstDerivatives(u0, v0, eval, Du, Dv);

  const PointType fp_u = tri.evaluate(u0 + h, v0);
  const PointType fm_u = tri.evaluate(u0 - h, v0);
  const PointType fp_v = tri.evaluate(u0, v0 + h);
  const PointType fm_v = tri.evaluate(u0, v0 - h);

  const VectorType Du_fd(fp_u, fm_u);
  const VectorType Dv_fd(fp_v, fm_v);

  for(int d = 0; d < 3; ++d)
  {
    EXPECT_NEAR(Du[d], Du_fd[d] / (2 * h), tol);
    EXPECT_NEAR(Dv[d], Dv_fd[d] / (2 * h), tol);
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziertriangle, gregorytriangle_finite_difference_first_derivatives)
{
  using CoordType = double;
  using PointType = primal::Point<CoordType, 3>;
  using VectorType = primal::Vector<CoordType, 3>;
  using GTri = primal::GregoryTriangle<CoordType>;

  constexpr CoordType h = 1e-7;
  constexpr CoordType tol = 5e-5;

  const std::array<PointType, 3> corners = {
    PointType {0.0, 0.0, 0.0},
    PointType {1.0, 0.0, 0.2},
    PointType {0.0, 1.0, -0.1},
  };

  const std::array<VectorType, 3> normals = {
    VectorType {0.0, 0.0, 1.0},
    VectorType {0.0, 0.0, 1.0},
    VectorType {0.0, 0.0, 1.0},
  };

  GTri gtri(axom::ArrayView<const PointType>(corners.data(), 3),
            axom::ArrayView<const VectorType>(normals.data(), 3));

  const CoordType u0 = 0.21;
  const CoordType v0 = 0.27;

  PointType eval;
  VectorType Du, Dv;
  gtri.evaluateFirstDerivatives(u0, v0, eval, Du, Dv);

  const PointType fp_u = gtri.evaluate(u0 + h, v0);
  const PointType fm_u = gtri.evaluate(u0 - h, v0);
  const PointType fp_v = gtri.evaluate(u0, v0 + h);
  const PointType fm_v = gtri.evaluate(u0, v0 - h);

  const VectorType Du_fd(fp_u, fm_u);
  const VectorType Dv_fd(fp_v, fm_v);

  for(int d = 0; d < 3; ++d)
  {
    EXPECT_NEAR(Du[d], Du_fd[d] / (2 * h), tol);
    EXPECT_NEAR(Dv[d], Dv_fd[d] / (2 * h), tol);
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziertriangle, gregorytriangle_parameter_convention_matches_beziertriangle)
{
  using CoordType = double;
  using PointType = primal::Point<CoordType, 3>;
  using VectorType = primal::Vector<CoordType, 3>;
  using GTri = primal::GregoryTriangle<CoordType>;

  constexpr CoordType eps = 1e-12;

  const std::array<PointType, 3> corners = {
    PointType {0.0, 0.0, 0.0},
    PointType {1.0, 0.0, 0.0},
    PointType {0.0, 1.0, 0.0},
  };

  const std::array<VectorType, 3> normals = {
    VectorType {0.0, 0.0, 1.0},
    VectorType {0.0, 0.0, 1.0},
    VectorType {0.0, 0.0, 1.0},
  };

  GTri gtri(axom::ArrayView<const PointType>(corners.data(), 3),
            axom::ArrayView<const VectorType>(normals.data(), 3));

  // Corner mapping should match BezierTriangle::evaluate convention
  EXPECT_TRUE(gtri.evaluate(0.0, 0.0).isNearlyEqual(corners[0], eps));
  EXPECT_TRUE(gtri.evaluate(0.0, 1.0).isNearlyEqual(corners[1], eps));
  EXPECT_TRUE(gtri.evaluate(1.0, 0.0).isNearlyEqual(corners[2], eps));

  // Edge mapping: verify straight edges for this planar configuration
  auto check_near = [&](const PointType& a, const PointType& b) {
    EXPECT_NEAR(a[0], b[0], eps);
    EXPECT_NEAR(a[1], b[1], eps);
    EXPECT_NEAR(a[2], b[2], eps);
  };

  const CoordType tvals[] = {0.0, 0.2, 0.6, 1.0};

  // edge u0 = 0: from corner0 -> corner1 as v0 goes 0->1
  for(const CoordType t : tvals)
  {
    const PointType expected = PointType::lerp(corners[0], corners[1], t);
    check_near(gtri.evaluate(0.0, t), expected);
  }

  // edge v0 = 0: from corner0 -> corner2 as u0 goes 0->1
  for(const CoordType t : tvals)
  {
    const PointType expected = PointType::lerp(corners[0], corners[2], t);
    check_near(gtri.evaluate(t, 0.0), expected);
  }

  // edge u0 + v0 = 1: from corner1 -> corner2 as u0 goes 0->1 (v0 goes 1->0)
  for(const CoordType t : tvals)
  {
    const PointType expected = PointType::lerp(corners[1], corners[2], t);
    check_near(gtri.evaluate(t, 1.0 - t), expected);
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

//------------------------------------------------------------------------------
TEST(primal_beziertriangle, edges)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using BTri = primal::BezierTriangle<CoordType, DIM>;
  using PointType = BTri::PointType;

  constexpr int ord = 3;
  BTri poly(ord);
  BTri rat(ord);
  rat.getWeights().resize(rat.getControlPoints().size());

  for(int i = 0; i <= ord; ++i)
  {
    for(int j = 0; j <= ord - i; ++j)
    {
      const int idx = BTri::triIndex(ord, i, j);
      const auto ii = static_cast<CoordType>(i);
      const auto jj = static_cast<CoordType>(j);
      poly(i, j) = PointType {ii, jj, 100. * ii + jj};
      rat(i, j) = poly(i, j);
      rat.getWeights()[idx] = 0.25 + idx;
    }
  }

  // edge 0: u + v = 1, from (0,1) to (ord,0)
  {
    const auto e0 = poly.getEdge(0);
    EXPECT_EQ(ord, e0.getOrder());
    EXPECT_FALSE(e0.isRational());
    EXPECT_EQ(ord + 1, e0.getNumControlPoints());
    EXPECT_EQ(poly(0, ord), e0.getInitPoint());
    EXPECT_EQ(poly(ord, 0), e0.getEndPoint());
  }

  // edge 1: v = 0, from (ord,0) to (0,0)
  {
    const auto e1 = poly.getEdge(1);
    EXPECT_EQ(ord, e1.getOrder());
    EXPECT_FALSE(e1.isRational());
    EXPECT_EQ(ord + 1, e1.getNumControlPoints());
    EXPECT_EQ(poly(ord, 0), e1.getInitPoint());
    EXPECT_EQ(poly(0, 0), e1.getEndPoint());
  }

  // edge 2: u = 0, from (0,0) to (0,ord)
  {
    const auto e2 = poly.getEdge(2);
    EXPECT_EQ(ord, e2.getOrder());
    EXPECT_FALSE(e2.isRational());
    EXPECT_EQ(ord + 1, e2.getNumControlPoints());
    EXPECT_EQ(poly(0, 0), e2.getInitPoint());
    EXPECT_EQ(poly(0, ord), e2.getEndPoint());
  }

  // weight propagation to rational edge curves
  {
    const auto e0 = rat.getEdge(0);
    EXPECT_TRUE(e0.isRational());
    EXPECT_EQ(ord, e0.getOrder());
    ASSERT_EQ(ord + 1, e0.getWeights().size());

    const auto e1 = rat.getEdge(1);
    EXPECT_TRUE(e1.isRational());
    EXPECT_EQ(ord, e1.getOrder());
    ASSERT_EQ(ord + 1, e1.getWeights().size());

    const auto e2 = rat.getEdge(2);
    EXPECT_TRUE(e2.isRational());
    EXPECT_EQ(ord, e2.getOrder());
    ASSERT_EQ(ord + 1, e2.getWeights().size());
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziertriangle, split_interior_polynomial)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using BTri = primal::BezierTriangle<CoordType, DIM>;

  constexpr int ord = 3;
  BTri tri(ord);

  fillControlNet(tri);

  const CoordType u = 0.2;
  const CoordType v = 0.3;

  BTri t0, t1, t2;
  tri.split(u, v, t0, t1, t2);

  EXPECT_EQ(ord, t0.getOrder());
  EXPECT_EQ(ord, t1.getOrder());
  EXPECT_EQ(ord, t2.getOrder());
  EXPECT_FALSE(t0.isRational());
  EXPECT_FALSE(t1.isRational());
  EXPECT_FALSE(t2.isRational());

  // Vertex mapping checks
  const auto pA = tri.evaluate(0.0, 0.0);
  const auto pB = tri.evaluate(0.0, 1.0);
  const auto pC = tri.evaluate(1.0, 0.0);
  const auto pQ = tri.evaluate(u, v);

  constexpr double tol = 1e-10;
  for(int i = 0; i < DIM; ++i)
  {
    // t0 -> (B, C, Q)
    EXPECT_NEAR(pB[i], t0.evaluate(0.0, 0.0)[i], tol);
    EXPECT_NEAR(pC[i], t0.evaluate(0.0, 1.0)[i], tol);
    EXPECT_NEAR(pQ[i], t0.evaluate(1.0, 0.0)[i], tol);

    // t1 -> (C, A, Q)
    EXPECT_NEAR(pC[i], t1.evaluate(0.0, 0.0)[i], tol);
    EXPECT_NEAR(pA[i], t1.evaluate(0.0, 1.0)[i], tol);
    EXPECT_NEAR(pQ[i], t1.evaluate(1.0, 0.0)[i], tol);

    // t2 -> (A, B, Q)
    EXPECT_NEAR(pA[i], t2.evaluate(0.0, 0.0)[i], tol);
    EXPECT_NEAR(pB[i], t2.evaluate(0.0, 1.0)[i], tol);
    EXPECT_NEAR(pQ[i], t2.evaluate(1.0, 0.0)[i], tol);
  }

  // Interior point checks via affine parameter mapping
  const CoordType s = 0.2;
  const CoordType t = 0.1;

  // t0 is over vertices (0,1), (1,0), (u,v)
  {
    const CoordType u2 = s * u + t;
    const CoordType v2 = 1.0 + s * (v - 1.0) - t;
    const auto expected = tri.evaluate(u2, v2);
    const auto actual = t0.evaluate(s, t);
    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_NEAR(expected[d], actual[d], tol);
    }
  }

  // t1 is over vertices (1,0), (0,0), (u,v)
  {
    const CoordType u3 = 1.0 + s * (u - 1.0) - t;
    const CoordType v3 = s * v;
    const auto expected = tri.evaluate(u3, v3);
    const auto actual = t1.evaluate(s, t);
    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_NEAR(expected[d], actual[d], tol);
    }
  }

  // t2 is over vertices (0,0), (0,1), (u,v)
  {
    const CoordType u1 = s * u;
    const CoordType v1 = t + s * v;
    const auto expected = tri.evaluate(u1, v1);
    const auto actual = t2.evaluate(s, t);
    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_NEAR(expected[d], actual[d], tol);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziertriangle, split_interior_rational)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using BTri = primal::BezierTriangle<CoordType, DIM>;

  constexpr int ord = 3;
  BTri tri(ord);
  fillControlNet(tri);
  makeRational(tri);

  const CoordType u = 0.2;
  const CoordType v = 0.3;

  BTri t0, t1, t2;
  tri.split(u, v, t0, t1, t2);

  EXPECT_EQ(ord, t0.getOrder());
  EXPECT_EQ(ord, t1.getOrder());
  EXPECT_EQ(ord, t2.getOrder());
  EXPECT_TRUE(t0.isRational());
  EXPECT_TRUE(t1.isRational());
  EXPECT_TRUE(t2.isRational());
  EXPECT_EQ(BTri::triSize(ord), t0.getWeights().size());
  EXPECT_EQ(BTri::triSize(ord), t1.getWeights().size());
  EXPECT_EQ(BTri::triSize(ord), t2.getWeights().size());

  // Vertex mapping checks
  const auto pA = tri.evaluate(0.0, 0.0);
  const auto pB = tri.evaluate(0.0, 1.0);
  const auto pC = tri.evaluate(1.0, 0.0);
  const auto pQ = tri.evaluate(u, v);

  constexpr double tol = 1e-10;
  for(int i = 0; i < DIM; ++i)
  {
    // t0 -> (B, C, Q)
    EXPECT_NEAR(pB[i], t0.evaluate(0.0, 0.0)[i], tol);
    EXPECT_NEAR(pC[i], t0.evaluate(0.0, 1.0)[i], tol);
    EXPECT_NEAR(pQ[i], t0.evaluate(1.0, 0.0)[i], tol);

    // t1 -> (C, A, Q)
    EXPECT_NEAR(pC[i], t1.evaluate(0.0, 0.0)[i], tol);
    EXPECT_NEAR(pA[i], t1.evaluate(0.0, 1.0)[i], tol);
    EXPECT_NEAR(pQ[i], t1.evaluate(1.0, 0.0)[i], tol);

    // t2 -> (A, B, Q)
    EXPECT_NEAR(pA[i], t2.evaluate(0.0, 0.0)[i], tol);
    EXPECT_NEAR(pB[i], t2.evaluate(0.0, 1.0)[i], tol);
    EXPECT_NEAR(pQ[i], t2.evaluate(1.0, 0.0)[i], tol);
  }

  // Interior point checks via affine parameter mapping
  const CoordType s = 0.2;
  const CoordType t = 0.1;

  // t0 is over vertices (0,1), (1,0), (u,v)
  {
    const CoordType u2 = s * u + t;
    const CoordType v2 = 1.0 + s * (v - 1.0) - t;
    const auto expected = tri.evaluate(u2, v2);
    const auto actual = t0.evaluate(s, t);
    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_NEAR(expected[d], actual[d], tol);
    }
  }

  // t1 is over vertices (1,0), (0,0), (u,v)
  {
    const CoordType u3 = 1.0 + s * (u - 1.0) - t;
    const CoordType v3 = s * v;
    const auto expected = tri.evaluate(u3, v3);
    const auto actual = t1.evaluate(s, t);
    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_NEAR(expected[d], actual[d], tol);
    }
  }

  // t2 is over vertices (0,0), (0,1), (u,v)
  {
    const CoordType u1 = s * u;
    const CoordType v1 = t + s * v;
    const auto expected = tri.evaluate(u1, v1);
    const auto actual = t2.evaluate(s, t);
    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_NEAR(expected[d], actual[d], tol);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziertriangle, split_edge_polynomial)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using BTri = primal::BezierTriangle<CoordType, DIM>;
  using PointType = BTri::PointType;

  constexpr int ord = 3;
  BTri tri(ord);

  fillControlNet(tri);

  const auto pA = tri.evaluate(0.0, 0.0);
  const auto pB = tri.evaluate(0.0, 1.0);
  const auto pC = tri.evaluate(1.0, 0.0);

  axom::Array<PointType> vertices {pA, pB, pC};

  const CoordType s = 0.35;

  BTri t0, t1;

  for(int i = 0; i < 3; ++i)
  {
    tri.split(i, s, t0, t1);

    EXPECT_EQ(ord, t0.getOrder());
    EXPECT_EQ(ord, t1.getOrder());
    EXPECT_FALSE(t0.isRational());
    EXPECT_FALSE(t1.isRational());

    // Vertex mapping checks
    const auto pQ = tri.getEdge(i).evaluate(s);

    for(int N = 0; N < DIM; ++N)
    {
      EXPECT_NEAR(vertices[(i + 0) % 3][N], t0.evaluate(0.0, 0.0)[N], 1e-10);
      EXPECT_NEAR(vertices[(i + 1) % 3][N], t0.evaluate(0.0, 1.0)[N], 1e-10);

      EXPECT_NEAR(vertices[(i + 2) % 3][N], t1.evaluate(0.0, 0.0)[N], 1e-10);
      EXPECT_NEAR(vertices[(i + 0) % 3][N], t1.evaluate(0.0, 1.0)[N], 1e-10);

      // Subtriangles agree at the last vertex
      EXPECT_NEAR(pQ[N], t0.evaluate(1.0, 0.0)[N], 1e-10);
      EXPECT_NEAR(pQ[N], t1.evaluate(1.0, 0.0)[N], 1e-10);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziertriangle, split_edge_rational)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using BTri = primal::BezierTriangle<CoordType, DIM>;
  using PointType = BTri::PointType;

  constexpr int ord = 3;
  BTri tri(ord);
  fillControlNet(tri);
  makeRational(tri);

  const auto pA = tri.evaluate(0.0, 0.0);
  const auto pB = tri.evaluate(0.0, 1.0);
  const auto pC = tri.evaluate(1.0, 0.0);

  axom::Array<PointType> vertices {pA, pB, pC};

  const CoordType s = 0.35;

  BTri t0, t1;

  for(int i = 0; i < 3; ++i)
  {
    tri.split(i, s, t0, t1);

    EXPECT_EQ(ord, t0.getOrder());
    EXPECT_EQ(ord, t1.getOrder());
    EXPECT_TRUE(t0.isRational());
    EXPECT_TRUE(t1.isRational());
    EXPECT_EQ(BTri::triSize(ord), t0.getWeights().size());
    EXPECT_EQ(BTri::triSize(ord), t1.getWeights().size());

    // Vertex mapping checks
    const auto pQ = tri.getEdge(i).evaluate(s);

    for(int N = 0; N < DIM; ++N)
    {
      EXPECT_NEAR(vertices[(i + 0) % 3][N], t0.evaluate(0.0, 0.0)[N], 1e-10);
      EXPECT_NEAR(vertices[(i + 1) % 3][N], t0.evaluate(0.0, 1.0)[N], 1e-10);

      EXPECT_NEAR(vertices[(i + 2) % 3][N], t1.evaluate(0.0, 0.0)[N], 1e-10);
      EXPECT_NEAR(vertices[(i + 0) % 3][N], t1.evaluate(0.0, 1.0)[N], 1e-10);

      // Subtriangles agree at the last vertex
      EXPECT_NEAR(pQ[N], t0.evaluate(1.0, 0.0)[N], 1e-10);
      EXPECT_NEAR(pQ[N], t1.evaluate(1.0, 0.0)[N], 1e-10);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziertriangle, split_fourway_polynomial)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using BTri = primal::BezierTriangle<CoordType, DIM>;

  constexpr int ord = 3;
  BTri tri(ord);

  fillControlNet(tri);

  const CoordType s0 = 0.25;  // edge0: B->C
  const CoordType s1 = 0.6;   // edge1: C->A
  const CoordType s2 = 0.4;   // edge2: A->B

  BTri t1, t2, t3, t4;
  tri.split(s0, s1, s2, t1, t2, t3, t4);

  EXPECT_EQ(ord, t1.getOrder());
  EXPECT_EQ(ord, t2.getOrder());
  EXPECT_EQ(ord, t3.getOrder());
  EXPECT_EQ(ord, t4.getOrder());
  EXPECT_FALSE(t1.isRational());
  EXPECT_FALSE(t2.isRational());
  EXPECT_FALSE(t3.isRational());
  EXPECT_FALSE(t4.isRational());

  const auto pA = tri.evaluate(0.0, 0.0);
  const auto pB = tri.evaluate(0.0, 1.0);
  const auto pC = tri.evaluate(1.0, 0.0);

  const auto pP0 = tri.getEdge(0).evaluate(s0);
  const auto pP1 = tri.getEdge(1).evaluate(s1);
  const auto pP2 = tri.getEdge(2).evaluate(s2);

  constexpr double tol = 1e-10;
  for(int d = 0; d < DIM; ++d)
  {
    // t1 = Tri(A, P2, P1)
    EXPECT_NEAR(pA[d], t1.evaluate(0.0, 0.0)[d], tol);
    EXPECT_NEAR(pP2[d], t1.evaluate(0.0, 1.0)[d], tol);
    EXPECT_NEAR(pP1[d], t1.evaluate(1.0, 0.0)[d], tol);

    // t2 = Tri(B, P0, P2)
    EXPECT_NEAR(pB[d], t2.evaluate(0.0, 0.0)[d], tol);
    EXPECT_NEAR(pP0[d], t2.evaluate(0.0, 1.0)[d], tol);
    EXPECT_NEAR(pP2[d], t2.evaluate(1.0, 0.0)[d], tol);

    // t3 = Tri(C, P1, P0)
    EXPECT_NEAR(pC[d], t3.evaluate(0.0, 0.0)[d], tol);
    EXPECT_NEAR(pP1[d], t3.evaluate(0.0, 1.0)[d], tol);
    EXPECT_NEAR(pP0[d], t3.evaluate(1.0, 0.0)[d], tol);

    // t4 = Tri(P1, P0, P2)
    EXPECT_NEAR(pP1[d], t4.evaluate(0.0, 0.0)[d], tol);
    EXPECT_NEAR(pP0[d], t4.evaluate(0.0, 1.0)[d], tol);
    EXPECT_NEAR(pP2[d], t4.evaluate(1.0, 0.0)[d], tol);
  }

  // Shared interior edges agree (same geometry and orientation for the chosen outputs)
  const CoordType s = 0.37;
  const auto eP0P2_from_t2 = t2.getEdge(0).evaluate(s);  // P0 -> P2
  const auto eP0P2_from_t4 = t4.getEdge(0).evaluate(s);  // P0 -> P2

  const auto eP2P1_from_t1 = t1.getEdge(0).evaluate(s);  // P2 -> P1
  const auto eP2P1_from_t4 = t4.getEdge(1).evaluate(s);  // P2 -> P1

  const auto eP1P0_from_t3 = t3.getEdge(0).evaluate(s);  // P1 -> P0
  const auto eP1P0_from_t4 = t4.getEdge(2).evaluate(s);  // P1 -> P0

  for(int d = 0; d < DIM; ++d)
  {
    EXPECT_NEAR(eP0P2_from_t2[d], eP0P2_from_t4[d], tol);
    EXPECT_NEAR(eP2P1_from_t1[d], eP2P1_from_t4[d], tol);
    EXPECT_NEAR(eP1P0_from_t3[d], eP1P0_from_t4[d], tol);
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziertriangle, split_fourway_rational)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using BTri = primal::BezierTriangle<CoordType, DIM>;

  constexpr int ord = 3;
  BTri tri(ord);
  fillControlNet(tri);
  makeRational(tri);

  const CoordType s0 = 0.25;  // edge0: B->C
  const CoordType s1 = 0.6;   // edge1: C->A
  const CoordType s2 = 0.4;   // edge2: A->B

  BTri t1, t2, t3, t4;
  tri.split(s0, s1, s2, t1, t2, t3, t4);

  EXPECT_EQ(ord, t1.getOrder());
  EXPECT_EQ(ord, t2.getOrder());
  EXPECT_EQ(ord, t3.getOrder());
  EXPECT_EQ(ord, t4.getOrder());
  EXPECT_TRUE(t1.isRational());
  EXPECT_TRUE(t2.isRational());
  EXPECT_TRUE(t3.isRational());
  EXPECT_TRUE(t4.isRational());

  const auto pA = tri.evaluate(0.0, 0.0);
  const auto pB = tri.evaluate(0.0, 1.0);
  const auto pC = tri.evaluate(1.0, 0.0);

  const auto pP0 = tri.getEdge(0).evaluate(s0);
  const auto pP1 = tri.getEdge(1).evaluate(s1);
  const auto pP2 = tri.getEdge(2).evaluate(s2);

  constexpr double tol = 1e-10;
  for(int d = 0; d < DIM; ++d)
  {
    // t1 = Tri(A, P2, P1)
    EXPECT_NEAR(pA[d], t1.evaluate(0.0, 0.0)[d], tol);
    EXPECT_NEAR(pP2[d], t1.evaluate(0.0, 1.0)[d], tol);
    EXPECT_NEAR(pP1[d], t1.evaluate(1.0, 0.0)[d], tol);

    // t2 = Tri(B, P0, P2)
    EXPECT_NEAR(pB[d], t2.evaluate(0.0, 0.0)[d], tol);
    EXPECT_NEAR(pP0[d], t2.evaluate(0.0, 1.0)[d], tol);
    EXPECT_NEAR(pP2[d], t2.evaluate(1.0, 0.0)[d], tol);

    // t3 = Tri(C, P1, P0)
    EXPECT_NEAR(pC[d], t3.evaluate(0.0, 0.0)[d], tol);
    EXPECT_NEAR(pP1[d], t3.evaluate(0.0, 1.0)[d], tol);
    EXPECT_NEAR(pP0[d], t3.evaluate(1.0, 0.0)[d], tol);

    // t4 = Tri(P1, P0, P2)
    EXPECT_NEAR(pP1[d], t4.evaluate(0.0, 0.0)[d], tol);
    EXPECT_NEAR(pP0[d], t4.evaluate(0.0, 1.0)[d], tol);
    EXPECT_NEAR(pP2[d], t4.evaluate(1.0, 0.0)[d], tol);
  }

  // Shared interior edges agree
  const CoordType s = 0.37;
  const auto eP0P2_from_t2 = t2.getEdge(0).evaluate(s);  // P0 -> P2
  const auto eP0P2_from_t4 = t4.getEdge(0).evaluate(s);  // P0 -> P2

  const auto eP2P1_from_t1 = t1.getEdge(0).evaluate(s);  // P2 -> P1
  const auto eP2P1_from_t4 = t4.getEdge(1).evaluate(s);  // P2 -> P1

  const auto eP1P0_from_t3 = t3.getEdge(0).evaluate(s);  // P1 -> P0
  const auto eP1P0_from_t4 = t4.getEdge(2).evaluate(s);  // P1 -> P0

  for(int d = 0; d < DIM; ++d)
  {
    EXPECT_NEAR(eP0P2_from_t2[d], eP0P2_from_t4[d], tol);
    EXPECT_NEAR(eP2P1_from_t1[d], eP2P1_from_t4[d], tol);
    EXPECT_NEAR(eP1P0_from_t3[d], eP1P0_from_t4[d], tol);
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziertriangle, uniformSplit_fourway_polynomial)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using BTri = primal::BezierTriangle<CoordType, DIM>;

  constexpr int ord = 3;
  BTri tri(ord);

  fillControlNet(tri);

  BTri t1, t2, t3, t4;
  tri.uniformSplit(t1, t2, t3, t4);

  // Compare against the general (non-optimized) split at s=0.5.
  BTri r1, r2, r3, r4;
  tri.split(0.5, 0.5, 0.5, r1, r2, r3, r4);

  const auto pP0 = tri.getEdge(0).evaluate(0.5);
  const auto pP1 = tri.getEdge(1).evaluate(0.5);
  const auto pP2 = tri.getEdge(2).evaluate(0.5);

  constexpr double tol = 1e-10;
  for(int d = 0; d < DIM; ++d)
  {
    EXPECT_NEAR(pP1[d], t4.evaluate(0.0, 0.0)[d], tol);
    EXPECT_NEAR(pP0[d], t4.evaluate(0.0, 1.0)[d], tol);
    EXPECT_NEAR(pP2[d], t4.evaluate(1.0, 0.0)[d], tol);
  }

  const CoordType s = 0.23;
  const CoordType t = 0.17;
  const auto p11 = t1.evaluate(s, t);
  const auto p12 = r1.evaluate(s, t);
  const auto p21 = t2.evaluate(s, t);
  const auto p22 = r2.evaluate(s, t);
  const auto p31 = t3.evaluate(s, t);
  const auto p32 = r3.evaluate(s, t);
  const auto p41 = t4.evaluate(s, t);
  const auto p42 = r4.evaluate(s, t);

  for(int d = 0; d < DIM; ++d)
  {
    EXPECT_NEAR(p11[d], p12[d], 1e-12);
    EXPECT_NEAR(p21[d], p22[d], 1e-12);
    EXPECT_NEAR(p31[d], p32[d], 1e-12);
    EXPECT_NEAR(p41[d], p42[d], 1e-12);
  }

  // Also confirm vertex ordering matches the general split for s=0.5
  for(int d = 0; d < DIM; ++d)
  {
    EXPECT_NEAR(t1.evaluate(0.0, 0.0)[d], r1.evaluate(0.0, 0.0)[d], tol);
    EXPECT_NEAR(t1.evaluate(0.0, 1.0)[d], r1.evaluate(0.0, 1.0)[d], tol);
    EXPECT_NEAR(t1.evaluate(1.0, 0.0)[d], r1.evaluate(1.0, 0.0)[d], tol);

    EXPECT_NEAR(t2.evaluate(0.0, 0.0)[d], r2.evaluate(0.0, 0.0)[d], tol);
    EXPECT_NEAR(t2.evaluate(0.0, 1.0)[d], r2.evaluate(0.0, 1.0)[d], tol);
    EXPECT_NEAR(t2.evaluate(1.0, 0.0)[d], r2.evaluate(1.0, 0.0)[d], tol);

    EXPECT_NEAR(t3.evaluate(0.0, 0.0)[d], r3.evaluate(0.0, 0.0)[d], tol);
    EXPECT_NEAR(t3.evaluate(0.0, 1.0)[d], r3.evaluate(0.0, 1.0)[d], tol);
    EXPECT_NEAR(t3.evaluate(1.0, 0.0)[d], r3.evaluate(1.0, 0.0)[d], tol);

    EXPECT_NEAR(t4.evaluate(0.0, 0.0)[d], r4.evaluate(0.0, 0.0)[d], tol);
    EXPECT_NEAR(t4.evaluate(0.0, 1.0)[d], r4.evaluate(0.0, 1.0)[d], tol);
    EXPECT_NEAR(t4.evaluate(1.0, 0.0)[d], r4.evaluate(1.0, 0.0)[d], tol);
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziertriangle, uniformSplit_fourway_rational)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using BTri = primal::BezierTriangle<CoordType, DIM>;

  constexpr int ord = 3;
  BTri tri(ord);
  fillControlNet(tri);
  makeRational(tri);

  BTri t1, t2, t3, t4;
  tri.uniformSplit(t1, t2, t3, t4);

  EXPECT_TRUE(t1.isRational());
  EXPECT_TRUE(t2.isRational());
  EXPECT_TRUE(t3.isRational());
  EXPECT_TRUE(t4.isRational());

  // Compare against the general (non-optimized) split at s=0.5.
  BTri r1, r2, r3, r4;
  tri.split(0.5, 0.5, 0.5, r1, r2, r3, r4);

  EXPECT_TRUE(r1.isRational());
  EXPECT_TRUE(r2.isRational());
  EXPECT_TRUE(r3.isRational());
  EXPECT_TRUE(r4.isRational());

  const CoordType s = 0.23;
  const CoordType t = 0.17;
  const auto p11 = t1.evaluate(s, t);
  const auto p12 = r1.evaluate(s, t);
  const auto p21 = t2.evaluate(s, t);
  const auto p22 = r2.evaluate(s, t);
  const auto p31 = t3.evaluate(s, t);
  const auto p32 = r3.evaluate(s, t);
  const auto p41 = t4.evaluate(s, t);
  const auto p42 = r4.evaluate(s, t);

  constexpr double tol = 1e-12;
  for(int d = 0; d < DIM; ++d)
  {
    EXPECT_NEAR(p11[d], p12[d], tol);
    EXPECT_NEAR(p21[d], p22[d], tol);
    EXPECT_NEAR(p31[d], p32[d], tol);
    EXPECT_NEAR(p41[d], p42[d], tol);
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
