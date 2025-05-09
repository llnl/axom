// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! 
 * \file primal_compute_moments.cpp
 * \brief This file tests primal's functionality related to computing moments
 */

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"
#include "axom/primal/operators/compute_moments.hpp"
#include "axom/primal/operators/detail/compute_moments_impl.hpp"

namespace primal = axom::primal;

const double EPS = 2e-15;

//------------------------------------------------------------------------------
TEST(primal_compute_moments, sector_area_cubic)
{
  const int DIM = 2;
  using T = double;
  using PointType = primal::Point<T, DIM>;
  using BezierCurveType = primal::BezierCurve<T, DIM>;
  using axom::utilities::isNearlyEqual;

  {
    SLIC_INFO("Testing Bezier sector area calculation for a cubic");
    const int order = 3;
    PointType data[order + 1] = {PointType {0.6, 1.2},
                                 PointType {1.3, 1.6},
                                 PointType {2.9, 2.4},
                                 PointType {3.2, 3.5}};

    BezierCurveType bCurve(data, order);
    const T area = primal::sector_area(bCurve);

    EXPECT_NEAR(.1455, area, EPS);
  }
}

//------------------------------------------------------------------------------
TEST(primal_compute_moments, sector_moment_cubic)
{
  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  {
    SLIC_INFO("Testing Bezier sector moment calculation for a cubic");
    const int order = 3;
    PointType data[order + 1] = {PointType {0.6, 1.2},
                                 PointType {1.3, 1.6},
                                 PointType {2.9, 2.4},
                                 PointType {3.2, 3.5}};

    BezierCurveType bCurve(data, order);
    PointType M = primal::sector_centroid(bCurve);
    EXPECT_NEAR(-.429321428571429, M[0], EPS);
    EXPECT_NEAR(-.354010714285715, M[1], EPS);
  }
}

//------------------------------------------------------------------------------
TEST(primal_compute_moments, sector_area_point)
{
  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  {
    SLIC_INFO("Testing Bezier sector area calculation for a point");
    const int order = 0;
    PointType data[order + 1] = {PointType {0.6, 1.2}};

    BezierCurveType bCurve(data, order);
    EXPECT_DOUBLE_EQ(0., primal::sector_area(bCurve));
  }
}

//------------------------------------------------------------------------------
TEST(primal_compute_moments, sector_moment_point)
{
  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  {
    SLIC_INFO("Testing Bezier sector moment calculation for a point");
    const int order = 0;
    PointType data[order + 1] = {PointType {0.6, 1.2}};

    BezierCurveType bCurve(data, order);
    PointType M = primal::sector_centroid(bCurve);
    EXPECT_DOUBLE_EQ(M[0], 0.0);
    EXPECT_DOUBLE_EQ(M[1], 0.0);
  }
}

//------------------------------------------------------------------------------
TEST(primal_compute_moments, sector_weights)
{
  SLIC_INFO("Testing weights for BezierCurve::sectorArea()");

  // NOTE: Expected weights are provided in the reference paper [Ueda99]
  // See doxygen comment for primal::sector_area(BezierCurve)

  using CoordType = double;
  primal::detail::MemoizedSectorAreaWeights<CoordType> memoizedSectorWeights;

  // order 1
  {
    const int ord = 1;
    auto weights = memoizedSectorWeights.getWeights(ord);

    double binomInv = 1. / axom::utilities::binomialCoefficient(2, 1);
    axom::numerics::Matrix<CoordType> exp(ord + 1, ord + 1);
    // clang-format off
    exp(0,0) =  0; exp(0,1) =  1;
    exp(1,0) = -1; exp(1,1) =  0;
    // clang-format on

    for(int i = 0; i <= ord; ++i)
    {
      for(int j = 0; j <= ord; ++j)
      {
        EXPECT_DOUBLE_EQ(exp(i, j) * binomInv, weights(i, j));
      }
    }
  }

  // order 2
  {
    const int ord = 2;
    auto weights = memoizedSectorWeights.getWeights(ord);

    double binomInv = 1. / axom::utilities::binomialCoefficient(4, 2);
    axom::numerics::Matrix<CoordType> exp(ord + 1, ord + 1);
    // clang-format off
    exp(0,0) =  0; exp(0,1) =  2; exp(0,2) =  1;
    exp(1,0) = -2; exp(1,1) =  0; exp(1,2) =  2;
    exp(2,0) = -1; exp(2,1) = -2; exp(2,2) =  0;
    // clang-format on

    for(int i = 0; i <= ord; ++i)
    {
      for(int j = 0; j <= ord; ++j)
      {
        EXPECT_DOUBLE_EQ(exp(i, j) * binomInv, weights(i, j));
      }
    }
  }

  // order 3
  {
    const int ord = 3;
    auto weights = memoizedSectorWeights.getWeights(ord);

    double binomInv = 1. / axom::utilities::binomialCoefficient(6, 3);
    axom::numerics::Matrix<CoordType> exp(ord + 1, ord + 1);
    // clang-format off
    exp(0,0) =  0; exp(0,1) =  6; exp(0,2) =  3; exp(0,3) =  1;
    exp(1,0) = -6; exp(1,1) =  0; exp(1,2) =  3; exp(1,3) =  3;
    exp(2,0) = -3; exp(2,1) = -3; exp(2,2) =  0; exp(2,3) =  6;
    exp(3,0) = -1; exp(3,1) = -3; exp(3,2) = -6; exp(3,3) =  0;
    // clang-format on

    for(int i = 0; i <= ord; ++i)
    {
      for(int j = 0; j <= ord; ++j)
      {
        EXPECT_DOUBLE_EQ(exp(i, j) * binomInv, weights(i, j));
      }
    }
  }

  // order 4
  {
    const int ord = 4;
    auto weights = memoizedSectorWeights.getWeights(ord);

    double binomInv = 1. / axom::utilities::binomialCoefficient(8, 4);
    axom::numerics::Matrix<CoordType> exp(ord + 1, ord + 1);
    // clang-format off
    exp(0,0) =  0; exp(0,1) = 20; exp(0,2) = 10; exp(0,3) =  4; exp(0,4) =  1;
    exp(1,0) =-20; exp(1,1) =  0; exp(1,2) =  8; exp(1,3) =  8; exp(1,4) =  4;
    exp(2,0) =-10; exp(2,1) = -8; exp(2,2) =  0; exp(2,3) =  8; exp(2,4) = 10;
    exp(3,0) = -4; exp(3,1) = -8; exp(3,2) = -8; exp(3,3) =  0; exp(3,4) = 20;
    exp(4,0) = -1; exp(4,1) = -4; exp(4,2) =-10; exp(4,3) =-20; exp(4,4) =  0;
    // clang-format on

    for(int i = 0; i <= ord; ++i)
    {
      for(int j = 0; j <= ord; ++j)
      {
        EXPECT_DOUBLE_EQ(exp(i, j) * binomInv, weights(i, j));
      }
    }
  }

  // order 5
  {
    const int ord = 5;
    auto weights = memoizedSectorWeights.getWeights(ord);

    double binomInv = 1. / axom::utilities::binomialCoefficient(10, 5);
    axom::numerics::Matrix<CoordType> exp(ord + 1, ord + 1);
    // clang-format off
    exp(0,0) =  0; exp(0,1) = 70; exp(0,2) = 35; exp(0,3) = 15; exp(0,4) =  5; exp(0,5) =  1;
    exp(1,0) =-70; exp(1,1) =  0; exp(1,2) = 25; exp(1,3) = 25; exp(1,4) = 15; exp(1,5) =  5;
    exp(2,0) =-35; exp(2,1) =-25; exp(2,2) =  0; exp(2,3) = 20; exp(2,4) = 25; exp(2,5) = 15;
    exp(3,0) =-15; exp(3,1) =-25; exp(3,2) =-20; exp(3,3) =  0; exp(3,4) = 25; exp(3,5) = 35;
    exp(4,0) = -5; exp(4,1) =-15; exp(4,2) =-25; exp(4,3) =-25; exp(4,4) =  0; exp(4,5) = 70;
    exp(5,0) = -1; exp(5,1) = -5; exp(5,2) =-15; exp(5,3) =-35; exp(5,4) =-70; exp(5,5) =  0;
    // clang-format on

    for(int i = 0; i <= ord; ++i)
    {
      for(int j = 0; j <= ord; ++j)
      {
        EXPECT_DOUBLE_EQ(exp(i, j) * binomInv, weights(i, j));
      }
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
