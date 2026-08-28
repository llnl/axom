// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
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
#include "axom/primal/operators/evaluate_integral_curve.hpp"
#include "axom/primal/operators/detail/compute_moments_impl.hpp"

#include <utility>

namespace primal = axom::primal;

constexpr double EPS = 2e-15;

namespace
{
using Point2D = primal::Point<double, 2>;
using BezierCurve2D = primal::BezierCurve<double, 2>;
using CurvedPolygon2D = primal::CurvedPolygon<BezierCurve2D>;

struct RawPlanarMoments
{
  double m00 {};
  double m10 {};
  double m01 {};
};

CurvedPolygon2D make_quadratic_curve_loop()
{
  return CurvedPolygon2D(axom::Array<BezierCurve2D> {
    BezierCurve2D(axom::Array<Point2D> {Point2D {2., 0.}, Point2D {1., 2.}, Point2D {0., 0.}}, 2),
    BezierCurve2D(axom::Array<Point2D> {Point2D {0., 0.}, Point2D {1., -2.}, Point2D {2., 0.}}, 2)});
}

CurvedPolygon2D make_cubic_curve_loop()
{
  return CurvedPolygon2D(axom::Array<BezierCurve2D> {
    BezierCurve2D(axom::Array<Point2D> {Point2D {0.6, 1.2},
                                        Point2D {1.3, 1.6},
                                        Point2D {2.9, 2.4},
                                        Point2D {3.2, 3.5}},
                  3),
    BezierCurve2D(axom::Array<Point2D> {Point2D {3.2, 3.5}, Point2D {0.6, 1.2}}, 1)});
}

RawPlanarMoments compute_ueda_raw_moments(const CurvedPolygon2D& polygon)
{
  const double area = primal::area(polygon);
  const Point2D centroid = primal::centroid(polygon);

  RawPlanarMoments moments;
  moments.m00 = area;
  moments.m10 = area * centroid[0];
  moments.m01 = area * centroid[1];
  return moments;
}

RawPlanarMoments compute_spectral_raw_moments(const CurvedPolygon2D& polygon, int gauss_points)
{
  auto integrate_spectral = [&polygon, gauss_points](auto&& integrand) {
    return primal::evaluate_area_integral(polygon,
                                          std::forward<decltype(integrand)>(integrand),
                                          gauss_points,
                                          gauss_points);
  };

  RawPlanarMoments moments;
  moments.m00 = integrate_spectral([](Point2D) { return 1.; });
  moments.m10 = integrate_spectral([](const Point2D& x) { return x[0]; });
  moments.m01 = integrate_spectral([](const Point2D& x) { return x[1]; });
  return moments;
}

void expect_raw_moments_near(const RawPlanarMoments& actual,
                             const RawPlanarMoments& expected,
                             double tol)
{
  EXPECT_NEAR(actual.m00, expected.m00, tol);
  EXPECT_NEAR(actual.m10, expected.m10, tol);
  EXPECT_NEAR(actual.m01, expected.m01, tol);
}

void expect_raw_moments_negated(const RawPlanarMoments& actual,
                                const RawPlanarMoments& expected,
                                double tol)
{
  EXPECT_NEAR(actual.m00, -expected.m00, tol);
  EXPECT_NEAR(actual.m10, -expected.m10, tol);
  EXPECT_NEAR(actual.m01, -expected.m01, tol);
}
}  // namespace

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

    EXPECT_NEAR(-.1455, area, EPS);
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
    EXPECT_NEAR(.429321428571429, M[0], EPS);
    EXPECT_NEAR(.354010714285715, M[1], EPS);
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
TEST(primal_compute_moments, spectral_matches_ueda_for_polynomial_bezier_regions)
{
  // Cross-check the closed-form moments against our quadrature on polynomial Bezier regions
  const CurvedPolygon2D quadratic_region = make_quadratic_curve_loop();
  const CurvedPolygon2D cubic_region = make_cubic_curve_loop();

  expect_raw_moments_near(compute_spectral_raw_moments(quadratic_region, 8),
                          compute_ueda_raw_moments(quadratic_region),
                          1e-12);
  expect_raw_moments_near(compute_spectral_raw_moments(cubic_region, 20),
                          compute_ueda_raw_moments(cubic_region),
                          1e-10);
}

//------------------------------------------------------------------------------
TEST(primal_compute_moments, reversed_orientation_negates_raw_moments)
{
  // Reversing a closed boundary changes the signed measure but not the centroid
  const CurvedPolygon2D region = make_cubic_curve_loop();
  CurvedPolygon2D reversed_region = region;
  reversed_region.reverseOrientation();

  const RawPlanarMoments moments = compute_ueda_raw_moments(region);
  const RawPlanarMoments reversed_moments = compute_ueda_raw_moments(reversed_region);
  expect_raw_moments_negated(reversed_moments, moments, 1e-14);

  const Point2D centroid = primal::centroid(region);
  const Point2D reversed_centroid = primal::centroid(reversed_region);
  EXPECT_NEAR(reversed_centroid[0], centroid[0], 1e-14);
  EXPECT_NEAR(reversed_centroid[1], centroid[1], 1e-14);
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
