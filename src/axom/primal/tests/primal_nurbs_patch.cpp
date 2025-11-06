// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! 
 * \file primal_nurbs_patch.cpp
 * \brief This file tests primal's NURBS patch functionality
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/fmt.hpp"

#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/geometry/NURBSPatch.hpp"
#include "axom/primal/operators/squared_distance.hpp"

#include "axom/core/numerics/matvecops.hpp"

namespace primal = axom::primal;

//------------------------------------------------------------------------------
TEST(primal_nurbspatch, sizing_constructors)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;

  auto check_patch =
    [=](const NURBSPatchType& nPatch, int deg_u, int deg_v, int npts_u, int npts_v, bool expect_empty) {
      EXPECT_EQ(!expect_empty, nPatch.getKnots_u().isValid());
      EXPECT_EQ(!expect_empty, nPatch.getKnots_v().isValid());
      EXPECT_EQ(!expect_empty, nPatch.isValidNURBS());

      EXPECT_EQ(deg_u, nPatch.getDegree_u());
      EXPECT_EQ(deg_v, nPatch.getDegree_v());

      EXPECT_EQ(deg_u + 1, nPatch.getOrder_u());
      EXPECT_EQ(deg_v + 1, nPatch.getOrder_v());

      EXPECT_EQ(npts_u, nPatch.getControlPoints().shape()[0]);
      EXPECT_EQ(npts_v, nPatch.getControlPoints().shape()[1]);
      EXPECT_EQ(npts_u * npts_v, nPatch.getControlPoints().size());

      EXPECT_EQ(npts_u + deg_u + 1, nPatch.getKnots_u().getArray().size());
      EXPECT_EQ(npts_v + deg_v + 1, nPatch.getKnots_v().getArray().size());

      EXPECT_FALSE(nPatch.isRational());
    };

  {
    SCOPED_TRACE("Default NURBS Patch constructor ");
    NURBSPatchType nPatch;
    check_patch(nPatch, -1, -1, 0, 0, true);
  }

  {
    SCOPED_TRACE("Empty NURBS Patch constructor ");
    NURBSPatchType nPatch(-1, -1);
    check_patch(nPatch, -1, -1, 0, 0, true);
  }

  for(int deg_u = 0; deg_u < 5; ++deg_u)
  {
    for(int deg_v = 0; deg_v < 5; ++deg_v)
    {
      {
        SCOPED_TRACE(
          axom::fmt::format("NURBS Patch constructor with deg_u={}, deg_v={}", deg_u, deg_v));
        NURBSPatchType patch = NURBSPatchType(deg_u, deg_v);
        check_patch(patch, deg_u, deg_v, deg_u + 1, deg_v + 1, false);
      }

      for(int npts_u = deg_u + 1; npts_u < deg_u + 5; ++npts_u)
      {
        for(int npts_v = deg_v + 1; npts_v < deg_v + 5; ++npts_v)
        {
          SCOPED_TRACE(axom::fmt::format(
            "NURBS Patch constructor with npts_u={}, npts_v={}, deg_u={}, deg_v={}",
            npts_u,
            npts_v,
            deg_u,
            deg_v));
          NURBSPatchType nPatch(npts_u, npts_v, deg_u, deg_v);
          check_patch(nPatch, deg_u, deg_v, npts_u, npts_v, false);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbspatch, bezier_constructors)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;

  constexpr int ord_u = 2;
  constexpr int ord_v = 3;

  // clang-format off
  axom::Array<PointType> controlPoints {
    PointType {0.0, 0.0, 0.0}, PointType {0.0, 1.0,  1.0}, PointType {0.0, 2.0,  2.0},
    PointType {1.0, 0.0, 3.0}, PointType {1.0, 1.0,  4.0}, PointType {1.0, 2.0,  5.0},
    PointType {2.0, 0.0, 6.0}, PointType {2.0, 1.0,  7.0}, PointType {2.0, 2.0,  8.0},
    PointType {3.0, 0.0, 9.0}, PointType {3.0, 1.0, 10.0}, PointType {3.0, 2.0, 11.0}};

  axom::Array<double> weights { 0.009,  1.019,  2.029, 
                                3.109,  4.119,  5.129, 
                                6.209,  7.219,  8.229, 
                                9.309, 10.319, 11.329};
  // clang-format on

  primal::BezierPatch<double, DIM> nBez(controlPoints, ord_u, ord_v);
  EXPECT_FALSE(nBez.isRational());

  primal::BezierPatch<double, DIM> rBez(controlPoints, weights, ord_u, ord_v);
  EXPECT_TRUE(rBez.isRational());

  for(const auto& bez : {nBez, rBez})
  {
    NURBSPatchType patch(bez);
    EXPECT_TRUE(patch.isValidNURBS());
    EXPECT_EQ(patch.isRational(), bez.isRational());

    EXPECT_EQ(ord_u, patch.getDegree_u());
    EXPECT_EQ(ord_v, patch.getDegree_v());
    EXPECT_EQ(ord_u + 1, patch.getOrder_u());
    EXPECT_EQ(ord_v + 1, patch.getOrder_v());
    EXPECT_EQ(ord_u + 1, patch.getNumControlPoints_u());
    EXPECT_EQ(ord_v + 1, patch.getNumControlPoints_v());
    EXPECT_EQ(controlPoints.size(), patch.getControlPoints().size());
    EXPECT_EQ(bez.getControlPoints().size(), patch.getControlPoints().size());

    for(int p = 0; p < ord_u + 1; ++p)
    {
      for(int q = 0; q < ord_v + 1; ++q)
      {
        EXPECT_EQ(patch(p, q), bez(p, q));
        if(patch.isRational())
        {
          EXPECT_EQ(patch.getWeight(p, q), bez.getWeight(p, q));
        }
      }
    }

    constexpr int nkts_u = 2 * (ord_u + 1);
    constexpr int nkts_v = 2 * (ord_v + 1);

    EXPECT_EQ(nkts_u, patch.getNumKnots_u());
    EXPECT_EQ(nkts_v, patch.getNumKnots_v());

    for(int k_u = 0; k_u < ord_u + 1; ++k_u)
    {
      patch.getKnots_u()[k_u] = 0.;
      patch.getKnots_u()[nkts_u - k_u - 1] = 1.;
    }

    for(int k_v = 0; k_v < ord_v + 1; ++k_v)
    {
      patch.getKnots_v()[k_v] = 0.;
      patch.getKnots_v()[nkts_v - k_v - 1] = 1.;
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbspatch, knotless_array_constructors)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;

  SLIC_INFO("Testing point array constructor");

  constexpr int degree_u = 2;
  constexpr int degree_v = 1;

  constexpr int npts_u = 3;
  constexpr int npts_v = 4;

  // clang-format off
  // Construct from C-style arrays
  PointType controlPoints[12] = {
    PointType {0.0, 0.0, 0.0}, PointType {0.0, 1.0,  1.0}, PointType {0.0, 2.0, 2.0},
    PointType {1.0, 0.0, 3.0}, PointType {1.0, 1.0,  4.0}, PointType {1.0, 2.0, 5.0},
    PointType {2.0, 0.0, 6.0}, PointType {2.0, 1.0,  7.0}, PointType {2.0, 2.0, 8.0},
    PointType {3.0, 0.0, 9.0}, PointType {3.0, 1.0, 10.0}, PointType {3.0, 2.0, 11.0}};

  CoordType weights[12] = { 0.009,  1.019,  2.029, 
                            3.109,  4.119,  5.129, 
                            6.209,  7.219,  8.229, 
                            9.309, 10.319, 11.329};
  // clang-format on

  auto check_patch =
    [=](const NURBSPatchType& patch, int deg_u, int deg_v, int npts_u, int npts_v, bool expect_rational) {
      EXPECT_EQ(deg_u, patch.getDegree_u());
      EXPECT_EQ(deg_v, patch.getDegree_v());

      EXPECT_EQ(npts_u, patch.getControlPoints().shape()[0]);
      EXPECT_EQ(npts_v, patch.getControlPoints().shape()[1]);
      EXPECT_EQ(npts_u * npts_v, patch.getControlPoints().size());

      EXPECT_EQ(npts_u + deg_u + 1, patch.getKnots_u().getArray().size());
      EXPECT_EQ(npts_v + deg_v + 1, patch.getKnots_v().getArray().size());

      if(expect_rational)
      {
        EXPECT_EQ(npts_u, patch.getWeights().shape()[0]);
        EXPECT_EQ(npts_v, patch.getWeights().shape()[1]);
        EXPECT_EQ(npts_u * npts_v, patch.getWeights().size());
      }
      else
      {
        EXPECT_TRUE(patch.getWeights().empty());
      }

      int idx = 0;
      for(int u = 0; u < npts_u; ++u)
      {
        for(int v = 0; v < npts_v; ++v, ++idx)
        {
          EXPECT_EQ(patch(u, v), controlPoints[idx]);
          if(expect_rational)
          {
            EXPECT_EQ(patch.getWeight(u, v), weights[idx]);
          }
        }
      }
    };

  // test C-array constructors
  {
    NURBSPatchType nPatch(controlPoints, npts_u, npts_v, degree_u, degree_v);
    check_patch(nPatch, degree_u, degree_v, npts_u, npts_v, false);

    NURBSPatchType wPatch(controlPoints, weights, npts_u, npts_v, degree_u, degree_v);
    check_patch(wPatch, degree_u, degree_v, npts_u, npts_v, true);
  }

  // test 1D axom::Array constructors
  {
    axom::Array<PointType> cp;
    cp.assign(std::begin(controlPoints), std::end(controlPoints));

    axom::Array<double> w;
    w.assign(std::begin(weights), std::end(weights));

    NURBSPatchType nPatch(cp, npts_u, npts_v, degree_u, degree_v);
    check_patch(nPatch, degree_u, degree_v, npts_u, npts_v, false);

    NURBSPatchType wPatch(cp, w, npts_u, npts_v, degree_u, degree_v);
    check_patch(wPatch, degree_u, degree_v, npts_u, npts_v, true);
  }

  // test 2D axom::Array constructors
  {
    axom::Array<PointType, 2> controlPointsArray2D(npts_u, npts_v);
    axom::Array<double, 2> weightsArray2D(npts_u, npts_v);

    int idx = 0;
    for(int p = 0; p < npts_u; ++p)
    {
      for(int q = 0; q < npts_v; ++q, ++idx)
      {
        controlPointsArray2D(p, q) = controlPoints[idx];
        weightsArray2D(p, q) = weights[idx];
      }
    }

    NURBSPatchType nPatch(controlPointsArray2D, degree_u, degree_v);
    check_patch(nPatch, degree_u, degree_v, npts_u, npts_v, false);

    NURBSPatchType wPatch(controlPointsArray2D, weightsArray2D, degree_u, degree_v);
    check_patch(wPatch, degree_u, degree_v, npts_u, npts_v, true);

    NURBSPatchType nPatch_vw(controlPointsArray2D.view(), degree_u, degree_v);
    check_patch(nPatch_vw, degree_u, degree_v, npts_u, npts_v, false);

    NURBSPatchType wPatch_vw(controlPointsArray2D.view(), weightsArray2D.view(), degree_u, degree_v);
    check_patch(wPatch_vw, degree_u, degree_v, npts_u, npts_v, true);
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbspatch, knot_array_constructor)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;
  using KnotVectorType = primal::KnotVector<CoordType>;

  SLIC_INFO("Testing knot array constructor");

  constexpr int degree_u = 1;
  constexpr int degree_v = 1;

  constexpr int npts_u = 3;
  constexpr int npts_v = 3;

  // clang-format off
  PointType controlPoints[npts_u * npts_v] = {
    PointType {0.0, 0.0, 1.0}, PointType {0.0, 1.0,  0.0}, PointType {0.0, 2.0, 0.0},
    PointType {1.0, 0.0, 0.0}, PointType {1.0, 1.0, -1.0}, PointType {1.0, 2.0, 0.0},
    PointType {2.0, 0.0, 0.0}, PointType {2.0, 1.0,  0.0}, PointType {2.0, 2.0, 1.0}};

  CoordType weights[npts_u * npts_v] = {1.0, 2.0, 3.0, 
                                        2.0, 3.0, 4.0, 
                                        3.0, 4.0, 5.0};
  // clang-format on

  double knots_u[npts_u + degree_u + 1] = {0.0, 0.0, 0.5, 1.0, 1.0};
  double knots_v[npts_v + degree_v + 1] = {0.0, 0.0, 0.5, 1.0, 1.0};

  auto check_patch =
    [=](const NURBSPatchType& patch, int deg_u, int deg_v, int npts_u, int npts_v, bool expect_rational) {
      EXPECT_EQ(deg_u, patch.getDegree_u());
      EXPECT_EQ(deg_v, patch.getDegree_v());

      EXPECT_EQ(npts_u, patch.getControlPoints().shape()[0]);
      EXPECT_EQ(npts_v, patch.getControlPoints().shape()[1]);
      EXPECT_EQ(npts_u * npts_v, patch.getControlPoints().size());

      EXPECT_EQ(npts_u + deg_u + 1, patch.getKnots_u().getArray().size());
      EXPECT_EQ(npts_v + deg_v + 1, patch.getKnots_v().getArray().size());

      if(expect_rational)
      {
        EXPECT_EQ(npts_u, patch.getWeights().shape()[0]);
        EXPECT_EQ(npts_v, patch.getWeights().shape()[1]);
        EXPECT_EQ(npts_u * npts_v, patch.getWeights().size());
      }
      else
      {
        EXPECT_TRUE(patch.getWeights().empty());
      }

      int idx = 0;
      for(int u = 0; u < npts_u; ++u)
      {
        for(int v = 0; v < npts_v; ++v, ++idx)
        {
          EXPECT_EQ(patch(u, v), controlPoints[idx]);
          if(expect_rational)
          {
            EXPECT_EQ(patch.getWeight(u, v), weights[idx]);
          }
        }
      }
    };

  // test from C-style arrays
  {
    NURBSPatchType nPatch(controlPoints,
                          npts_u,
                          npts_v,
                          knots_u,
                          npts_u + degree_u + 1,
                          knots_v,
                          npts_v + degree_v + 1);
    check_patch(nPatch, degree_u, degree_v, npts_u, npts_v, false);

    NURBSPatchType wPatch(controlPoints,
                          weights,
                          npts_u,
                          npts_v,
                          knots_u,
                          npts_u + degree_u + 1,
                          knots_v,
                          npts_v + degree_v + 1);
    check_patch(wPatch, degree_u, degree_v, npts_u, npts_v, true);
  }

  // test with 1D axom::Arrays
  {
    axom::Array<PointType> cp_arr;
    cp_arr.assign(std::begin(controlPoints), std::end(controlPoints));
    axom::Array<double> w_arr;
    w_arr.assign(std::begin(weights), std::end(weights));

    axom::Array<double> knot_arr_u;
    knot_arr_u.assign(std::begin(knots_u), std::end(knots_u));
    axom::Array<double> knot_arr_v;
    knot_arr_v.assign(std::begin(knots_v), std::end(knots_v));

    KnotVectorType knotvec_u(knot_arr_u, degree_u);
    KnotVectorType knotvec_v(knot_arr_v, degree_v);

    {
      SCOPED_TRACE("Testing 1D array constructors with 1D array of knots");

      NURBSPatchType nPatch(cp_arr, npts_u, npts_v, knot_arr_u, knot_arr_v);
      check_patch(nPatch, degree_u, degree_v, npts_u, npts_v, false);

      NURBSPatchType wPatch(cp_arr, w_arr, npts_u, npts_v, knot_arr_u, knot_arr_v);
      check_patch(wPatch, degree_u, degree_v, npts_u, npts_v, true);
    }
    {
      SCOPED_TRACE("Testing 1D array constructors with KnotVectors");

      NURBSPatchType nPatch(cp_arr, npts_u, npts_v, knotvec_u, knotvec_v);
      check_patch(nPatch, degree_u, degree_v, npts_u, npts_v, false);

      NURBSPatchType wPatch(cp_arr, w_arr, npts_u, npts_v, knotvec_u, knotvec_v);
      check_patch(wPatch, degree_u, degree_v, npts_u, npts_v, true);
    }
  }

  // test with 2D axom::Arrays
  {
    axom::Array<PointType, 2> cp_arr_2D(npts_u, npts_v);
    axom::Array<double, 2> w_arr_2D(npts_u, npts_v);

    int idx = 0;
    for(int p = 0; p < npts_u; ++p)
    {
      for(int q = 0; q < npts_v; ++q, ++idx)
      {
        cp_arr_2D(p, q) = controlPoints[idx];
        w_arr_2D(p, q) = weights[idx];
      }
    }

    axom::Array<double> knot_arr_u;
    knot_arr_u.assign(std::begin(knots_u), std::end(knots_u));
    axom::Array<double> knot_arr_v;
    knot_arr_v.assign(std::begin(knots_v), std::end(knots_v));

    KnotVectorType knotvec_u(knot_arr_u, degree_u);
    KnotVectorType knotvec_v(knot_arr_v, degree_v);

    {
      SCOPED_TRACE("Testing 2D array constructors with 1D array of knots");
      NURBSPatchType nPatch(cp_arr_2D, knot_arr_u, knot_arr_v);
      check_patch(nPatch, degree_u, degree_v, npts_u, npts_v, false);

      NURBSPatchType wPatch(cp_arr_2D, w_arr_2D, knot_arr_u, knot_arr_v);
      check_patch(wPatch, degree_u, degree_v, npts_u, npts_v, true);
    }

    {
      SCOPED_TRACE("Testing 2D array constructors with KnotVector");
      NURBSPatchType nPatch(cp_arr_2D, knotvec_u, knotvec_v);
      check_patch(nPatch, degree_u, degree_v, npts_u, npts_v, false);

      NURBSPatchType wPatch(cp_arr_2D, w_arr_2D, knotvec_u, knotvec_v);
      check_patch(wPatch, degree_u, degree_v, npts_u, npts_v, true);
    }

    {
      SCOPED_TRACE("Testing 2D ArrayView constructors with KnotVector");
      NURBSPatchType nPatch(cp_arr_2D.view(), knotvec_u, knotvec_v);
      check_patch(nPatch, degree_u, degree_v, npts_u, npts_v, false);

      NURBSPatchType wPatch(cp_arr_2D.view(), w_arr_2D.view(), knotvec_u, knotvec_v);
      check_patch(wPatch, degree_u, degree_v, npts_u, npts_v, true);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbspatch, set_degree)
{
  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;

  SLIC_INFO("Test adding control points to an empty NURBS patch");

  NURBSPatchType nPatch;
  EXPECT_EQ(nPatch.getDegree_u(), -1);
  EXPECT_EQ(nPatch.getDegree_v(), -1);

  const int degree_u = 1;
  const int degree_v = 1;

  const int npts_u = 3;
  const int npts_v = 3;

  nPatch.setNumControlPoints(npts_u, npts_v);
  nPatch.setDegree(degree_u, degree_v);

  EXPECT_EQ(nPatch.getNumControlPoints_u(), npts_u);
  EXPECT_EQ(nPatch.getNumControlPoints_v(), npts_v);

  EXPECT_EQ(nPatch.getDegree_u(), degree_u);
  EXPECT_EQ(nPatch.getDegree_v(), degree_v);

  EXPECT_EQ(nPatch.getNumKnots_u(), degree_u + npts_u + 1);
  EXPECT_EQ(nPatch.getNumKnots_v(), degree_v + npts_v + 1);

  // clang-format off
  PointType controlPoints[9] = {
    PointType {0.0, 0.0, 1.0}, PointType {0.0, 1.0,  0.0}, PointType {0.0, 2.0, 0.0},
    PointType {1.0, 0.0, 0.0}, PointType {1.0, 1.0, -1.0}, PointType {1.0, 2.0, 0.0},
    PointType {2.0, 0.0, 0.0}, PointType {2.0, 1.0,  0.0}, PointType {2.0, 2.0, 1.0}};

  nPatch(0, 0) = controlPoints[0]; nPatch(0, 1) = controlPoints[1]; nPatch(0, 2) = controlPoints[2];
  nPatch(1, 0) = controlPoints[3]; nPatch(1, 1) = controlPoints[4]; nPatch(1, 2) = controlPoints[5];
  nPatch(2, 0) = controlPoints[6]; nPatch(2, 1) = controlPoints[7]; nPatch(2, 2) = controlPoints[8];
  // clang-format on

  for(int p = 0; p < npts_u; ++p)
  {
    for(int q = 0; q < npts_v; ++q)
    {
      auto& pt = nPatch(p, q);
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_DOUBLE_EQ(controlPoints[p * npts_u + q][i], pt[i]);
      }
    }
  }

  nPatch.clear();
  EXPECT_EQ(nPatch.getDegree_u(), -1);
  EXPECT_EQ(nPatch.getDegree_v(), -1);
  EXPECT_FALSE(nPatch.isRational());

  nPatch.setParameters(npts_u, npts_v, degree_u, degree_v);
  nPatch.makeRational();
  EXPECT_TRUE(nPatch.isRational());

  nPatch.setWeight(0, 0, 2.0);
  EXPECT_DOUBLE_EQ(2.0, nPatch.getWeight(0, 0));
}

//------------------------------------------------------------------------------
TEST(primal_nurbspatch, isocurve_evaluate)
{
  SLIC_INFO("Testing NURBS Patch isocurve evaluation");
  // Test that the isocurves of a NURBS patch are correct,
  //  which we will use to test other evaluation routines

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;

  const int degree_u = 3;
  const int degree_v = 2;

  const int npts_u = 4;
  const int npts_v = 5;

  // clang-format off
  PointType controlPoints[4 * 5] = {
    PointType {0, 0, 0}, PointType {0, 4,  0}, PointType {0, 8, -3}, PointType {0, 12, 0}, PointType {0, 16, 0},
    PointType {2, 0, 6}, PointType {2, 4,  0}, PointType {2, 8,  0}, PointType {2, 12, 0}, PointType {2, 16, 0},
    PointType {4, 0, 0}, PointType {4, 4,  0}, PointType {4, 8,  3}, PointType {4, 12, 0}, PointType {4, 16, 0},
    PointType {6, 0, 0}, PointType {6, 4, -3}, PointType {6, 8,  0}, PointType {6, 12, 0}, PointType {6, 16, 0}};
  
  double weights[4 * 5] = {
    1.0, 2.0, 3.0, 2.0, 1.0,
    2.0, 3.0, 4.0, 3.0, 2.0,
    3.0, 4.0, 5.0, 4.0, 3.0,
    4.0, 5.0, 6.0, 5.0, 4.0};
  // clang-format on

  NURBSPatchType nPatch(controlPoints, weights, npts_u, npts_v, degree_u, degree_v);

  // Test isocurves formed at u/v = 0, 0.5, 0.7
  double isocurve_params[3] = {0.0, 0.5, 0.7};

  // clang-format off
  PointType isocurves_u[3][5] = {
    {PointType {0, 0, 0}, PointType {0, 4, 0}, PointType {0, 8, -3}, PointType {0, 12, 0}, PointType {0, 16, 0}},
    {PointType {3.6, 0, 1.8}, PointType {24./7., 4, -15./28.}, PointType {10./3., 8, 1}, PointType {24./7., 12, 0}, PointType {3.6, 16, 0}},
    {PointType {714./155., 0, 567./775.}, PointType {924./205., 4, -1029./820.}, PointType {378./85., 8, 531./425.}, PointType {924./205., 12, 0}, PointType {714./155., 16, 0}}
    };

  double weights_u[3][5] = {
    {1.0, 2.0, 3.0, 2.0, 1.0},
    {2.5, 3.5, 4.5, 3.5, 2.5},
    {3.1, 4.1, 5.1, 4.1, 3.1}
    };

  PointType isocurves_v[3][4] = {
    {PointType {0, 0, 0}, PointType {2, 0, 6}, PointType {4, 0, 0}, PointType {6, 0, 0}},
    {PointType {0, 8, -27./11.}, PointType {2, 8, 0}, PointType {4, 8, 45./19.}, PointType {6, 8, -15./46.}},
    {PointType {0, 4784./479., -729./479.}, PointType {2, 6868./679., 0}, PointType {4, 5968./586., 810./586.}, PointType {6, 11036./1079., 0}}
    };
  
  double weights_v[3][4] = {
    {1.0, 2.0, 3.0, 4.0},
    {2.75, 3.75, 4.75, 5.75},
    {2.395, 3.395, 4.395, 5.395}
    };
  // clang-format on

  for(int i = 0; i < 3; ++i)
  {
    auto isocurve_u = nPatch.isocurve_u(isocurve_params[i]);
    auto isocurve_v = nPatch.isocurve_v(isocurve_params[i]);

    for(int j = 0; j < 5; ++j)
    {
      auto pt_u = isocurve_u[j];
      auto weight_u = isocurve_u.getWeight(j);

      EXPECT_NEAR(weight_u, weights_u[i][j], 1e-10);
      for(int k = 0; k < DIM; ++k)
      {
        EXPECT_NEAR(pt_u[k], isocurves_u[i][j][k], 1e-10);
      }
    }

    for(int j = 0; j < 4; ++j)
    {
      auto pt_v = isocurve_v[j];
      auto weight_v = isocurve_v.getWeight(j);

      EXPECT_NEAR(weight_v, weights_v[i][j], 1e-10);
      for(int k = 0; k < DIM; ++k)
      {
        EXPECT_NEAR(pt_v[k], isocurves_v[i][j][k], 1e-10);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbspatch, surface_evaluate)
{
  SLIC_INFO("Testing NURBS Patch surface evaluation");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;

  const int npts_u = 5;
  const int npts_v = 4;

  const int degree_u = 3;
  const int degree_v = 2;

  // clang-format off
  PointType controlPoints[5 * 4] = {
    PointType {0, 0, 0}, PointType {0, 4,  0}, PointType {0, 8, -3}, PointType {0, 12, 0},
    PointType {2, 0, 6}, PointType {2, 4,  0}, PointType {2, 8,  0}, PointType {2, 12, 0},
    PointType {4, 0, 0}, PointType {4, 4,  0}, PointType {4, 8,  3}, PointType {4, 12, 0},
    PointType {6, 0, 0}, PointType {6, 4, -3}, PointType {6, 8,  0}, PointType {6, 12, 0},
    PointType {8, 0, 0}, PointType {8, 4,  0}, PointType {8, 8,  0}, PointType {8, 12, 0}};

  double weights[5 * 4] = {
    1.0, 2.0, 3.0, 2.0,
    2.0, 3.0, 4.0, 3.0,
    3.0, 4.0, 5.0, 4.0,
    4.0, 5.0, 6.0, 5.0,
    5.0, 6.0, 7.0, 6.0};
  // clang-format on

  NURBSPatchType nPatch(controlPoints, weights, npts_u, npts_v, degree_u, degree_v);

  // isocurve_u should *fix* a value of u, returning a curve parameterized by v
  EXPECT_EQ(nPatch.isocurve_u(0.5).getDegree(), degree_v);
  EXPECT_EQ(nPatch.isocurve(0.5, 0).getDegree(), degree_v);

  // isocurve_v should *fix* a value of v, returning a curve parameterized by u
  EXPECT_EQ(nPatch.isocurve_v(0.5).getDegree(), degree_u);
  EXPECT_EQ(nPatch.isocurve(0.5, 1).getDegree(), degree_u);

  // Loop over the parameter space of the surface,
  //  and check that `evalaute` matches the results of the two isocurve methods

  constexpr int npts = 11;
  double u_pts[npts], v_pts[npts];
  axom::numerics::linspace(0.0, 1.0, u_pts, npts);
  axom::numerics::linspace(0.0, 1.0, v_pts, npts);

  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      auto pt = nPatch.evaluate(u, v);
      auto pt_u = nPatch.isocurve_u(u).evaluate(v);
      auto pt_v = nPatch.isocurve_v(v).evaluate(u);

      for(int N = 0; N < DIM; ++N)
      {
        EXPECT_NEAR(pt[N], pt_u[N], 1e-10);
        EXPECT_NEAR(pt[N], pt_v[N], 1e-10);
        EXPECT_NEAR(pt_u[N], pt_v[N], 1e-10);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbspatch, first_second_derivatives)
{
  SLIC_INFO("Testing NURBS Patch derivative evaluation");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;

  const int npts_u = 5;
  const int npts_v = 4;

  const int degree_u = 3;
  const int degree_v = 2;

  // clang-format off
  PointType controlPoints[5 * 4] = {
    PointType {0, 0, 0}, PointType {0, 4,  0}, PointType {0, 8, -3}, PointType {0, 12, 0},
    PointType {2, 0, 6}, PointType {2, 4,  0}, PointType {2, 8,  0}, PointType {2, 12, 0},
    PointType {4, 0, 0}, PointType {4, 4,  0}, PointType {4, 8,  3}, PointType {4, 12, 0},
    PointType {6, 0, 0}, PointType {6, 4, -3}, PointType {6, 8,  0}, PointType {6, 12, 0},
    PointType {8, 0, 0}, PointType {8, 4,  0}, PointType {8, 8,  0}, PointType {8, 12, 0}};

  double weights[5 * 4] = {
    1.0, 2.0, 3.0, 2.0,
    2.0, 3.0, 4.0, 3.0,
    3.0, 4.0, 5.0, 4.0,
    4.0, 5.0, 6.0, 5.0,
    5.0, 6.0, 7.0, 6.0};
  // clang-format on

  NURBSPatchType nPatch(controlPoints, weights, npts_u, npts_v, degree_u, degree_v);

  // Loop over the parameter space of the surface,
  //  and check that `evalauteDerivatives` matches the results of the two isocurve methods
  constexpr int npts = 11;
  double u_pts[npts], v_pts[npts];
  axom::numerics::linspace(0.0, 1.0, u_pts, npts);
  axom::numerics::linspace(0.0, 1.0, v_pts, npts);

  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      axom::Array<VectorType, 2> ders;
      nPatch.evaluateDerivatives(u, v, 2, ders);

      auto pt = nPatch.evaluate(u, v);
      auto pt_u = nPatch.isocurve_v(v).dt(u);
      auto pt_v = nPatch.isocurve_u(u).dt(v);
      auto pt_uu = nPatch.isocurve_v(v).dtdt(u);
      auto pt_vv = nPatch.isocurve_u(u).dtdt(v);

      for(int N = 0; N < DIM; ++N)
      {
        EXPECT_NEAR(pt[N], ders[0][0][N], 1e-10);
        EXPECT_NEAR(pt_u[N], ders[1][0][N], 1e-10);
        EXPECT_NEAR(pt_v[N], ders[0][1][N], 1e-10);
        EXPECT_NEAR(pt_uu[N], ders[2][0][N], 1e-10);
        EXPECT_NEAR(pt_vv[N], ders[0][2][N], 1e-10);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbspatch, knot_insertion)
{
  SLIC_INFO("Testing NURBS Patch knot insertion");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;

  const int npts_u = 5;
  const int npts_v = 4;

  const int degree_u = 3;
  const int degree_v = 2;

  // clang-format off
  PointType controlPoints[5 * 4] = {
    PointType {0, 0, 0}, PointType {0, 4,  0}, PointType {0, 8, -3}, PointType {0, 12, 0},
    PointType {2, 0, 6}, PointType {2, 4,  0}, PointType {2, 8,  0}, PointType {2, 12, 0},
    PointType {4, 0, 0}, PointType {4, 4,  0}, PointType {4, 8,  3}, PointType {4, 12, 0},
    PointType {6, 0, 0}, PointType {6, 4, -3}, PointType {6, 8,  0}, PointType {6, 12, 0},
    PointType {8, 0, 0}, PointType {8, 4,  0}, PointType {8, 8,  0}, PointType {8, 12, 0}};

  double weights[5 * 4] = {
    1.0, 2.0, 3.0, 2.0,
    2.0, 3.0, 4.0, 3.0,
    3.0, 4.0, 5.0, 4.0,
    4.0, 5.0, 6.0, 5.0,
    5.0, 6.0, 7.0, 6.0};
  // clang-format on

  NURBSPatchType nPatch(controlPoints, weights, npts_u, npts_v, degree_u, degree_v);

  NURBSPatchType nPatchExtraKnots(controlPoints, weights, npts_u, npts_v, degree_u, degree_v);

  // Insert knots in the u direction
  nPatchExtraKnots.insertKnot_u(0.3, 2);
  nPatchExtraKnots.insertKnot_u(0.5, 1);
  nPatchExtraKnots.insertKnot_u(0.7, 3);

  // Insert a knot in the v direction
  nPatchExtraKnots.insertKnot_v(0.4, 1);
  nPatchExtraKnots.insertKnot_v(0.6, 2);
  nPatchExtraKnots.insertKnot_v(0.8, 3);

  constexpr int npts = 11;
  double u_pts[npts], v_pts[npts];
  axom::numerics::linspace(0.0, 1.0, u_pts, npts);
  axom::numerics::linspace(0.0, 1.0, v_pts, npts);

  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      auto pt1 = nPatch.evaluate(u, v);
      auto pt2 = nPatchExtraKnots.evaluate(u, v);

      for(int N = 0; N < DIM; ++N)
      {
        EXPECT_NEAR(pt1[N], pt2[N], 1e-10);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbspatch, patch_split)
{
  SLIC_INFO("Testing NURBS Patch splitting");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;

  const int npts_u = 5;
  const int npts_v = 4;

  const int degree_u = 3;
  const int degree_v = 2;

  // clang-format off
  PointType controlPoints[5 * 4] = {
    PointType {0, 0, 0}, PointType {0, 4,  0}, PointType {0, 8, -3}, PointType {0, 12, 0},
    PointType {2, 0, 6}, PointType {2, 4,  0}, PointType {2, 8,  0}, PointType {2, 12, 0},
    PointType {4, 0, 0}, PointType {4, 4,  0}, PointType {4, 8,  3}, PointType {4, 12, 0},
    PointType {6, 0, 0}, PointType {6, 4, -3}, PointType {6, 8,  0}, PointType {6, 12, 0},
    PointType {8, 0, 0}, PointType {8, 4,  0}, PointType {8, 8,  0}, PointType {8, 12, 0}};

  double weights[5 * 4] = {
    1.0, 2.0, 3.0, 2.0,
    2.0, 3.0, 4.0, 3.0,
    3.0, 4.0, 5.0, 4.0,
    4.0, 5.0, 6.0, 5.0,
    5.0, 6.0, 7.0, 6.0};
  // clang-format on

  NURBSPatchType nPatch(controlPoints, weights, npts_u, npts_v, degree_u, degree_v);

  NURBSPatchType subpatch1, subpatch2;

  constexpr int npts = 11;
  double u_pts[npts], v_pts[npts];
  axom::numerics::linspace(0.0, 1.0, u_pts, npts);
  axom::numerics::linspace(0.0, 1.0, v_pts, npts);

  double split_vals[3] = {0.3, 0.5, 0.7};
  for(double val : split_vals)
  {
    nPatch.split_u(val, subpatch1, subpatch2);

    for(auto u : u_pts)
    {
      for(auto v : v_pts)
      {
        auto pt1 = nPatch.evaluate(u, v);

        if(u <= val)
        {
          auto pt2 = subpatch1.evaluate(u, v);
          for(int N = 0; N < DIM; ++N)
          {
            EXPECT_NEAR(pt1[N], pt2[N], 1e-10);
          }
        }
        else
        {
          auto pt2 = subpatch2.evaluate(u, v);
          for(int N = 0; N < DIM; ++N)
          {
            EXPECT_NEAR(pt1[N], pt2[N], 1e-10);
          }
        }
      }
    }

    nPatch.split_v(val, subpatch1, subpatch2);

    for(auto u : u_pts)
    {
      for(auto v : v_pts)
      {
        auto pt1 = nPatch.evaluate(u, v);

        if(v <= val)
        {
          auto pt2 = subpatch1.evaluate(u, v);
          for(int N = 0; N < DIM; ++N)
          {
            EXPECT_NEAR(pt1[N], pt2[N], 1e-10);
          }
        }
        else
        {
          auto pt2 = subpatch2.evaluate(u, v);
          for(int N = 0; N < DIM; ++N)
          {
            EXPECT_NEAR(pt1[N], pt2[N], 1e-10);
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbspatch, patch_clip)
{
  SLIC_INFO("Testing NURBS Patch clipping");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;

  const int npts_u = 5;
  const int npts_v = 4;

  const int degree_u = 3;
  const int degree_v = 2;

  // clang-format off
  PointType controlPoints[5 * 4] = {
    PointType {0, 0, 0}, PointType {0, 4,  0}, PointType {0, 8, -3}, PointType {0, 12, 0},
    PointType {2, 0, 6}, PointType {2, 4,  0}, PointType {2, 8,  0}, PointType {2, 12, 0},
    PointType {4, 0, 0}, PointType {4, 4,  0}, PointType {4, 8,  3}, PointType {4, 12, 0},
    PointType {6, 0, 0}, PointType {6, 4, -3}, PointType {6, 8,  0}, PointType {6, 12, 0},
    PointType {8, 0, 0}, PointType {8, 4,  0}, PointType {8, 8,  0}, PointType {8, 12, 0}};

  double weights[5 * 4] = {
    1.0, 2.0, 3.0, 2.0,
    2.0, 3.0, 4.0, 3.0,
    3.0, 4.0, 5.0, 4.0,
    4.0, 5.0, 6.0, 5.0,
    5.0, 6.0, 7.0, 6.0};
  // clang-format on

  NURBSPatchType nPatch(controlPoints, weights, npts_u, npts_v, degree_u, degree_v);
  NURBSPatchType subPatch(nPatch);

  // Clip the patch to the region 0.3 <= u <= 0.7, 0.4 <= v <= 0.6
  subPatch.clip(0.3, 0.7, 0.4, 0.6);

  EXPECT_NEAR(subPatch.getMinKnot_u(), 0.3, 1e-10);
  EXPECT_NEAR(subPatch.getMaxKnot_u(), 0.7, 1e-10);
  EXPECT_NEAR(subPatch.getMinKnot_v(), 0.4, 1e-10);
  EXPECT_NEAR(subPatch.getMaxKnot_v(), 0.6, 1e-10);

  constexpr int npts = 11;
  double u_pts[npts], v_pts[npts];
  axom::numerics::linspace(0.3, 0.7, u_pts, npts);
  axom::numerics::linspace(0.4, 0.6, v_pts, npts);

  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      auto pt1 = nPatch.evaluate(u, v);
      auto pt2 = subPatch.evaluate(u, v);

      for(int N = 0; N < DIM; ++N)
      {
        EXPECT_NEAR(pt1[N], pt2[N], 1e-10);
      }
    }
  }
}
//------------------------------------------------------------------------------
TEST(primal_nurbspatch, nurbs_parameter_space_scaling)
{
  SLIC_INFO("Testing NURBS Patch parameter space expansion");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;

  const int npts_u = 5;
  const int npts_v = 4;

  const int degree_u = 3;
  const int degree_v = 2;

  // clang-format off
  PointType controlPoints[5 * 4] = {
    PointType {0, 0, 0}, PointType {0, 4,  0}, PointType {0, 8, -3}, PointType {0, 12, 0},
    PointType {2, 0, 3}, PointType {2, 4,  0}, PointType {2, 8,  0}, PointType {2, 12, 0},
    PointType {4, 0, 0}, PointType {4, 4,  0}, PointType {4, 8,  3}, PointType {4, 12, 0},
    PointType {6, 0, 0}, PointType {6, 4, -3}, PointType {6, 8,  0}, PointType {6, 12, 0},
    PointType {8, 0, 0}, PointType {8, 4,  0}, PointType {8, 8,  0}, PointType {8, 12, 0}};

  double weights[5 * 4] = {
    1.0, 2.0, 3.0, 2.0,
    2.0, 3.0, 4.0, 3.0,
    3.0, 4.0, 5.0, 4.0,
    4.0, 5.0, 6.0, 5.0,
    5.0, 6.0, 7.0, 6.0};
  // clang-format on

  NURBSPatchType nPatchUntrimmed(controlPoints, weights, npts_u, npts_v, degree_u, degree_v);
  NURBSPatchType nPatchTrimmed(controlPoints, weights, npts_u, npts_v, degree_u, degree_v);
  nPatchTrimmed.addTrimmingCurve(
    primal::NURBSCurve<CoordType, 2>::make_circular_arc_nurbs(0.0, 2.0 * M_PI, 0.5, 0.5, 0.25));

  NURBSPatchType supPatchOriginallyUntrimmed(nPatchUntrimmed);
  NURBSPatchType supPatchUntrimmed(nPatchUntrimmed);
  NURBSPatchType supPatchTrimmed(nPatchTrimmed);

  // Expand the parameter space of the patch
  constexpr double scaleFactor = 1.05;

  // All patches have the same geometry, but different trimming curves
  constexpr bool removeTrimmingCurves = true;
  supPatchOriginallyUntrimmed.scaleParameterSpace(scaleFactor);
  supPatchUntrimmed.scaleParameterSpace(scaleFactor, removeTrimmingCurves);
  supPatchTrimmed.scaleParameterSpace(scaleFactor);

  // Both should be trimmed after this procedure UNLESS the flag is set
  EXPECT_TRUE(supPatchOriginallyUntrimmed.isTrimmed());
  EXPECT_FALSE(supPatchUntrimmed.isTrimmed());
  EXPECT_TRUE(supPatchTrimmed.isTrimmed());

  double min_u = nPatchUntrimmed.getMinKnot_u();
  double max_u = nPatchUntrimmed.getMaxKnot_u();

  double min_v = nPatchUntrimmed.getMinKnot_v();
  double max_v = nPatchUntrimmed.getMaxKnot_v();

  // Check the parameter space of the superpatches' knots
  for(auto& the_patch : {supPatchOriginallyUntrimmed, supPatchUntrimmed, supPatchTrimmed})
  {
    EXPECT_NEAR(the_patch.getMinKnot_u(), min_u - (scaleFactor - 1.0), 1e-10);
    EXPECT_NEAR(the_patch.getMaxKnot_u(), max_u + (scaleFactor - 1.0), 1e-10);
    EXPECT_NEAR(the_patch.getMinKnot_v(), min_v - (scaleFactor - 1.0), 1e-10);
    EXPECT_NEAR(the_patch.getMaxKnot_v(), max_v + (scaleFactor - 1.0), 1e-10);
  }

  // Check that the patches are equal in the original parameter space
  constexpr int npts = 15;
  double u_pts[npts], v_pts[npts];
  axom::numerics::linspace(min_u - (scaleFactor - 1.0), max_u + (scaleFactor - 1.0), u_pts, npts);
  axom::numerics::linspace(min_v - (scaleFactor - 1.0), max_v + (scaleFactor - 1.0), v_pts, npts);

  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      auto orig_pt = nPatchUntrimmed.evaluate(axom::utilities::clampVal(u, min_u, max_u),
                                              axom::utilities::clampVal(v, min_v, max_v));

      auto ext_pt1 = supPatchOriginallyUntrimmed.evaluate(u, v);
      auto ext_pt2 = supPatchUntrimmed.evaluate(u, v);
      auto ext_pt3 = supPatchTrimmed.evaluate(u, v);

      // If the point is in the original parameter space, the two patches should be equal
      if((u >= min_u) && (u <= max_u) && (v >= min_v) && (v <= max_v))
      {
        for(int N = 0; N < DIM; ++N)
        {
          EXPECT_NEAR(orig_pt[N], ext_pt1[N], 1e-10);
          EXPECT_NEAR(orig_pt[N], ext_pt2[N], 1e-10);
          EXPECT_NEAR(orig_pt[N], ext_pt3[N], 1e-10);
        }

        // Visibility on the original parameters should be unchanged after the patch is extended
        EXPECT_EQ(nPatchTrimmed.isVisible(u, v), supPatchTrimmed.isVisible(u, v));
        EXPECT_TRUE(supPatchOriginallyUntrimmed.isVisible(u, v));
        EXPECT_TRUE(supPatchUntrimmed.isVisible(u, v));
      }

      // If not, the points should be "nearby"
      else
      {
        // Check that the points are within a certain distance of each other
        EXPECT_LT(squared_distance(orig_pt, ext_pt1), 6.0 * 6.0);
        EXPECT_LT(squared_distance(orig_pt, ext_pt2), 6.0 * 6.0);
        EXPECT_LT(squared_distance(orig_pt, ext_pt3), 6.0 * 6.0);

        // Only the flagged patch should be visible in the extended parameters
        EXPECT_TRUE(supPatchUntrimmed.isVisible(u, v));
        EXPECT_FALSE(supPatchTrimmed.isVisible(u, v));
        EXPECT_FALSE(supPatchOriginallyUntrimmed.isVisible(u, v));
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbspatch, bezier_extraction)
{
  SLIC_INFO("Testing NURBS Patch Bezier extraction");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;

  const int npts_u = 5;
  const int npts_v = 4;

  const int degree_u = 4;
  const int degree_v = 3;

  // clang-format off
  PointType controlPoints[5 * 4] = {
    PointType {0, 0, 0}, PointType {0, 4,  0}, PointType {0, 8, -3}, PointType {0, 12, 0},
    PointType {2, 0, 6}, PointType {2, 4,  0}, PointType {2, 8,  0}, PointType {2, 12, 0},
    PointType {4, 0, 0}, PointType {4, 4,  0}, PointType {4, 8,  3}, PointType {4, 12, 0},
    PointType {6, 0, 0}, PointType {6, 4, -3}, PointType {6, 8,  0}, PointType {6, 12, 0},
    PointType {8, 0, 0}, PointType {8, 4,  0}, PointType {8, 8,  0}, PointType {8, 12, 0}};

  double weights[5 * 4] = {
    1.0, 2.0, 3.0, 2.0,
    2.0, 3.0, 4.0, 3.0,
    3.0, 4.0, 5.0, 4.0,
    4.0, 5.0, 6.0, 5.0,
    5.0, 6.0, 7.0, 6.0};
  // clang-format on

  NURBSPatchType nPatch(controlPoints, weights, npts_u, npts_v, degree_u, degree_v);

  // Do knot insertion, which determines where the Bezier splitting happens
  nPatch.insertKnot_u(0.33, 3);
  nPatch.insertKnot_u(0.66, 1);
  nPatch.insertKnot_u(0.77, 2);

  nPatch.insertKnot_v(0.25, 2);
  nPatch.insertKnot_v(0.5, 1);
  nPatch.insertKnot_v(0.75, 3);

  auto bezier_list = nPatch.extractBezier();

  EXPECT_EQ(bezier_list.size(), 16);

  constexpr int npts = 11;
  double u_pts[npts], v_pts[npts];

  double u_ranges[5] = {0, 0.33, 0.66, 0.77, 1};
  double v_ranges[5] = {0, 0.25, 0.5, 0.75, 1};

  // bezier_list is ordered lexicographically by v, then u
  for(int i = 0; i < 4; ++i)
  {
    for(int j = 0; j < 4; ++j)
    {
      auto& bPatch = bezier_list[i * 4 + j];

      // Loop over the parameter space of each Bezier patch
      axom::numerics::linspace(u_ranges[i], u_ranges[i + 1], u_pts, npts);
      axom::numerics::linspace(v_ranges[j], v_ranges[j + 1], v_pts, npts);

      for(auto u : u_pts)
      {
        for(auto v : v_pts)
        {
          auto pt1 = nPatch.evaluate(u, v);
          auto pt2 = bPatch.evaluate((u - u_ranges[i]) / (u_ranges[i + 1] - u_ranges[i]),
                                     (v - v_ranges[j]) / (v_ranges[j + 1] - v_ranges[j]));

          for(int N = 0; N < DIM; ++N)
          {
            EXPECT_NEAR(pt1[N], pt2[N], 1e-10);
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbspatch, evaluation_degenerate)
{
  SLIC_INFO("Testing NURBS patch evaluation with one degenerate axis");
  // Should reduce to a Bezier curve along the nonempty dimension

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;

  const int degree = 3;
  PointType data[degree + 1] = {PointType {0.6, 1.2, 1.0},
                                PointType {1.3, 1.6, 1.8},
                                PointType {2.9, 2.4, 2.3},
                                PointType {3.2, 3.5, 3.0}};

  NURBSCurveType nCurve(data, degree + 1, degree);
  NURBSPatchType nPatch(data, degree + 1, 1, degree, 0);

  constexpr int npts = 11;
  double t_pts[npts];
  axom::numerics::linspace(0.0, 1.0, t_pts, npts);

  for(auto t : t_pts)
  {
    for(int N = 0; N < DIM; ++N)
    {
      EXPECT_NEAR(nCurve.evaluate(t)[N], nPatch.evaluate(t, 0)[N], 1e-10);
      EXPECT_NEAR(nCurve.evaluate(t)[N], nPatch.evaluate(t, 0.5)[N], 1e-10);
      EXPECT_NEAR(nCurve.evaluate(t)[N], nPatch.evaluate(t, 1.0)[N], 1e-10);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbspatch, extract_degenerate)
{
  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;

  SLIC_INFO("Testing Bezier extraction on degenerate surface (order 0)");

  const int npts_u = 3;
  const int npts_v = 3;

  // Construct from C-style arrays
  PointType controlPoints[9] = {PointType {0.0, 0.0, 1.0},
                                PointType {0.0, 1.0, 0.0},
                                PointType {0.0, 2.0, 0.0},
                                PointType {1.0, 0.0, 0.0},
                                PointType {1.0, 1.0, -1.0},
                                PointType {1.0, 2.0, 0.0},
                                PointType {2.0, 0.0, 0.0},
                                PointType {2.0, 1.0, 0.0},
                                PointType {2.0, 2.0, 1.0}};

  CoordType weights[9] = {1.0, 2.0, 3.0, 2.0, 3.0, 4.0, 3.0, 4.0, 5.0};

  constexpr int npts = 11;
  double u_pts[npts], v_pts[npts];

  double u_ranges[5];
  double v_ranges[5];

  for(int degree_v = 1; degree_v <= 2; ++degree_v)
  {
    // Degenerate in u, full in v (degrees 1, 2)
    NURBSPatchType nPatch_u(controlPoints, weights, npts_u, npts_v, 0, degree_v);
    auto bezier_list_u = nPatch_u.extractBezier();

    const auto u_spans = nPatch_u.getKnots_u().getNumKnotSpans();
    const auto v_spans = nPatch_u.getKnots_v().getNumKnotSpans();

    EXPECT_EQ(bezier_list_u.size(), npts_u * v_spans);

    axom::numerics::linspace(0.0, 1.0, u_ranges, u_spans + 1);
    axom::numerics::linspace(0.0, 1.0, v_ranges, v_spans + 1);

    // bezier_list is ordered lexicographically by v, then u
    for(int i = 0; i < u_spans; ++i)
    {
      for(int j = 0; j < v_spans; ++j)
      {
        auto& bPatch = bezier_list_u[i * (v_spans) + j];

        // Loop over the parameter space of each Bezier patch
        axom::numerics::linspace(u_ranges[i], u_ranges[i + 1], u_pts, npts);
        axom::numerics::linspace(v_ranges[j], v_ranges[j + 1], v_pts, npts);

        // Order 0 degree curves are discontinuous, so don't check
        //  the boundaries in parameter space
        for(int ui = 1; ui < npts - 1; ++ui)
        {
          for(int vi = 1; vi < npts - 1; ++vi)
          {
            auto pt1 = nPatch_u.evaluate(u_pts[ui], v_pts[vi]);
            auto pt2 = bPatch.evaluate((u_pts[ui] - u_ranges[i]) / (u_ranges[i + 1] - u_ranges[i]),
                                       (v_pts[vi] - v_ranges[j]) / (v_ranges[j + 1] - v_ranges[j]));

            for(int N = 0; N < DIM; ++N)
            {
              EXPECT_NEAR(pt1[N], pt2[N], 1e-10);
            }
          }
        }
      }
    }
  }

  for(int degree_u = 1; degree_u <= 2; ++degree_u)
  {
    // Degenerate in v, full in u (degree 1, 2)
    NURBSPatchType nPatch_v(controlPoints, weights, npts_u, npts_v, degree_u, 0);
    auto bezier_list_v = nPatch_v.extractBezier();

    const auto u_spans = nPatch_v.getKnots_u().getNumKnotSpans();
    const auto v_spans = nPatch_v.getKnots_v().getNumKnotSpans();

    EXPECT_EQ(bezier_list_v.size(), npts_v * u_spans);

    axom::numerics::linspace(0.0, 1.0, u_ranges, u_spans + 1);
    axom::numerics::linspace(0.0, 1.0, v_ranges, v_spans + 1);

    // bezier_list is ordered lexicographically by v, then u
    for(int i = 0; i < u_spans; ++i)
    {
      for(int j = 0; j < v_spans; ++j)
      {
        auto& bPatch = bezier_list_v[i * (v_spans) + j];

        // Loop over the parameter space of each Bezier patch
        axom::numerics::linspace(u_ranges[i], u_ranges[i + 1], u_pts, npts);
        axom::numerics::linspace(v_ranges[j], v_ranges[j + 1], v_pts, npts);

        // Order 0 degree curves are discontinuous, so don't check
        //  the boundaries in parameter space
        for(int ui = 1; ui < npts - 1; ++ui)
        {
          for(int vi = 1; vi < npts - 1; ++vi)
          {
            auto pt1 = nPatch_v.evaluate(u_pts[ui], v_pts[vi]);
            auto pt2 = bPatch.evaluate((u_pts[ui] - u_ranges[i]) / (u_ranges[i + 1] - u_ranges[i]),
                                       (v_pts[vi] - v_ranges[j]) / (v_ranges[j + 1] - v_ranges[j]));

            for(int N = 0; N < DIM; ++N)
            {
              EXPECT_NEAR(pt1[N], pt2[N], 1e-10);
            }
          }
        }
      }
    }
  }

  // Degenerate in both u and v
  NURBSPatchType nPatch_uv(controlPoints, npts_u, npts_v, 0, 0);
  auto bezier_list_uv = nPatch_uv.extractBezier();

  EXPECT_EQ(bezier_list_uv.size(), npts_u * npts_v);

  axom::numerics::linspace(0.0, 1.0, u_pts, npts);
  axom::numerics::linspace(0.0, 1.0, v_pts, npts);

  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      auto pt1 = nPatch_uv.evaluate(u, v);

      if(u < 1.0 / 3.0 && v < 1.0 / 3.0)
      {
        auto pt2 = bezier_list_uv[0].evaluate(u, v);
        for(int N = 0; N < DIM; ++N)
        {
          EXPECT_NEAR(pt1[N], pt2[N], 1e-10);
        }
      }
      else if(u < 1.0 / 3.0 && 1.0 / 3.0 < v && v < 2.0 / 3.0)
      {
        auto pt2 = bezier_list_uv[1].evaluate(u, v);
        for(int N = 0; N < DIM; ++N)
        {
          EXPECT_NEAR(pt1[N], pt2[N], 1e-10);
        }
      }
      else if(u < 1.0 / 3.0 && v > 2.0 / 3.0)
      {
        auto pt2 = bezier_list_uv[2].evaluate(u, v);
        for(int N = 0; N < DIM; ++N)
        {
          EXPECT_NEAR(pt1[N], pt2[N], 1e-10);
        }
      }
      else if(1.0 / 3.0 < u && u < 2.0 / 3.0 && v < 1.0 / 3.0)
      {
        auto pt2 = bezier_list_uv[3].evaluate(u, v);
        for(int N = 0; N < DIM; ++N)
        {
          EXPECT_NEAR(pt1[N], pt2[N], 1e-10);
        }
      }
      else if(1.0 / 3.0 < u && u < 2.0 / 3.0 && 1.0 / 3.0 < v && v < 2.0 / 3.0)
      {
        auto pt2 = bezier_list_uv[4].evaluate(u, v);
        for(int N = 0; N < DIM; ++N)
        {
          EXPECT_NEAR(pt1[N], pt2[N], 1e-10);
        }
      }
      else if(1.0 / 3.0 < u && u < 2.0 / 3.0 && v > 2.0 / 3.0)
      {
        auto pt2 = bezier_list_uv[5].evaluate(u, v);
        for(int N = 0; N < DIM; ++N)
        {
          EXPECT_NEAR(pt1[N], pt2[N], 1e-10);
        }
      }
      else if(u > 2.0 / 3.0 && v < 1.0 / 3.0)
      {
        auto pt2 = bezier_list_uv[6].evaluate(u, v);
        for(int N = 0; N < DIM; ++N)
        {
          EXPECT_NEAR(pt1[N], pt2[N], 1e-10);
        }
      }
      else if(u > 2.0 / 3.0 && 1.0 / 3.0 < v && v < 2.0 / 3.0)
      {
        auto pt2 = bezier_list_uv[7].evaluate(u, v);
        for(int N = 0; N < DIM; ++N)
        {
          EXPECT_NEAR(pt1[N], pt2[N], 1e-10);
        }
      }
      else
      {
        auto pt2 = bezier_list_uv[8].evaluate(u, v);
        for(int N = 0; N < DIM; ++N)
        {
          EXPECT_NEAR(pt1[N], pt2[N], 1e-10);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbspatch, reverse_orientation)
{
  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;

  const int npts_u = 5;
  const int npts_v = 4;

  const int degree_u = 4;
  const int degree_v = 3;

  // clang-format off
  PointType controlPoints[5 * 4] = {
    PointType {0, 0, 0}, PointType {0, 4,  0}, PointType {0, 8, -3}, PointType {0, 12, 0},
    PointType {2, 0, 6}, PointType {2, 4,  0}, PointType {2, 8,  0}, PointType {2, 12, 0},
    PointType {4, 0, 0}, PointType {4, 4,  0}, PointType {4, 8,  3}, PointType {4, 12, 0},
    PointType {6, 0, 0}, PointType {6, 4, -3}, PointType {6, 8,  0}, PointType {6, 12, 0},
    PointType {8, 0, 0}, PointType {8, 4,  0}, PointType {8, 8,  0}, PointType {8, 12, 0}};

  double weights[5 * 4] = {
    1.0, 2.0, 3.0, 2.0,
    2.0, 3.0, 4.0, 3.0,
    3.0, 4.0, 5.0, 4.0,
    4.0, 5.0, 6.0, 5.0,
    5.0, 6.0, 7.0, 6.0};
  // clang-format on

  NURBSPatchType original(controlPoints, weights, npts_u, npts_v, degree_u, degree_v);

  // Add some knots to the patch in the u and v directions to make it more interesting
  original.insertKnot_u(0.3, 2);
  original.insertKnot_u(0.5, 1);

  original.insertKnot_v(0.4, 1);
  original.insertKnot_v(0.6, 2);

  double min_u = 1.0, max_u = 2.0;
  double min_v = -1.0, max_v = 0.5;

  original.rescale_u(min_u, max_u);
  original.rescale_v(min_v, max_v);

  NURBSPatchType reversed(original);

  // Reverse along the u-axis
  reversed.reverseOrientation(0);

  constexpr int npts = 11;
  double u_pts[npts], v_pts[npts];
  axom::numerics::linspace(min_u, max_u, u_pts, npts);
  axom::numerics::linspace(min_v, max_v, v_pts, npts);

  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      PointType o_pt = original.evaluate(u, v);
      PointType r_pt = reversed.evaluate(min_u + max_u - u, v);

      for(int N = 0; N < DIM; ++N)
      {
        EXPECT_NEAR(o_pt[N], r_pt[N], 1e-10);
      }
    }
  }

  // Reverse along the u-axis again, should return to original
  reversed.reverseOrientation(0);
  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      PointType o_pt = original.evaluate(u, v);
      PointType r_pt = reversed.evaluate(u, v);

      for(int N = 0; N < DIM; ++N)
      {
        EXPECT_NEAR(o_pt[N], r_pt[N], 1e-10);
      }
    }
  }

  // Reverse along the v-axis
  reversed.reverseOrientation(1);
  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      PointType o_pt = original.evaluate(u, v);
      PointType r_pt = reversed.evaluate(u, min_v + max_v - v);

      for(int N = 0; N < DIM; ++N)
      {
        EXPECT_NEAR(o_pt[N], r_pt[N], 1e-10);
      }
    }
  }

  // Reverse along the u-axis again
  reversed.reverseOrientation(0);
  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      PointType o_pt = original.evaluate(u, v);
      PointType r_pt = reversed.evaluate(min_u + max_u - u, min_v + max_v - v);

      for(int N = 0; N < DIM; ++N)
      {
        EXPECT_NEAR(o_pt[N], r_pt[N], 1e-10);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbspatch, swap_axes)
{
  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;

  const int npts_u = 5;
  const int npts_v = 4;

  const int degree_u = 4;
  const int degree_v = 3;

  // clang-format off
  PointType controlPoints[5 * 4] = {
    PointType {0, 0, 0}, PointType {0, 4,  0}, PointType {0, 8, -3}, PointType {0, 12, 0},
    PointType {2, 0, 6}, PointType {2, 4,  0}, PointType {2, 8,  0}, PointType {2, 12, 0},
    PointType {4, 0, 0}, PointType {4, 4,  0}, PointType {4, 8,  3}, PointType {4, 12, 0},
    PointType {6, 0, 0}, PointType {6, 4, -3}, PointType {6, 8,  0}, PointType {6, 12, 0},
    PointType {8, 0, 0}, PointType {8, 4,  0}, PointType {8, 8,  0}, PointType {8, 12, 0}};

  double weights[5 * 4] = {
    1.0, 2.0, 3.0, 2.0,
    2.0, 3.0, 4.0, 3.0,
    3.0, 4.0, 5.0, 4.0,
    4.0, 5.0, 6.0, 5.0,
    5.0, 6.0, 7.0, 6.0};
  // clang-format on

  NURBSPatchType original(controlPoints, weights, npts_u, npts_v, degree_u, degree_v);
  NURBSPatchType swapped(controlPoints, weights, npts_u, npts_v, degree_u, degree_v);

  // Swap the u and v axes
  swapped.swapAxes();

  constexpr int npts = 11;
  double u_pts[npts], v_pts[npts];
  axom::numerics::linspace(0.0, 1.0, u_pts, npts);
  axom::numerics::linspace(0.0, 1.0, v_pts, npts);

  for(auto u : u_pts)
  {
    for(auto v : v_pts)
    {
      PointType o_pt = original.evaluate(u, v);
      PointType s_pt = swapped.evaluate(v, u);

      for(int N = 0; N < DIM; ++N)
      {
        EXPECT_NEAR(o_pt[N], s_pt[N], 1e-10);
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
