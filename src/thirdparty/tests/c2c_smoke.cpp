// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"  // for gtest

#include "c2c/config.hpp"
#include "c2c/Literals.hpp"
#include "c2c/Point.hpp"
#include "c2c/EllipticalArc.hpp"
#include "c2c/NURBS.hpp"

#include <cmath>

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(c2c_smoke, check_version)
{
  std::cout << "Using c2c version: "     //
            << C2C_VERSION_MAJOR << "."  //
            << C2C_VERSION_MINOR << "."  //
            << C2C_VERSION_PATCH << std::endl;

  EXPECT_TRUE(C2C_VERSION_MAJOR >= 0);
  EXPECT_TRUE(C2C_VERSION_MINOR >= 0);
  EXPECT_TRUE(C2C_VERSION_PATCH >= 0);
}

TEST(c2c_smoke, check_code)
{
  double EPS = 1.e-8;
  using namespace c2c::literals;

  c2c::Length rRadius = 20.0_cm;
  c2c::Length zRadius = 150.0_mm;  // note: different unit
  c2c::Point origin {100.0_cm, 200.0_cm};
  c2c::RotationalDirection dir = c2c::RotationalDirection::IncreasingAngle;
  c2c::EllipticalArc arc {origin, rRadius, zRadius, 0.0_deg, 90.0_deg, dir};

  // Helper lambda to check if two c2c points are within tolerance EPS of each other
  auto c2cPointsAreNear = [=](const c2c::Point& p1, const c2c::Point& p2) {
    auto cm = c2c::LengthUnit::cm;
    return  //
      std::abs(p1.getR().convertValue(cm) - p2.getR().convertValue(cm)) < EPS &&
      std::abs(p1.getZ().convertValue(cm) - p2.getZ().convertValue(cm)) < EPS;
  };

  // Check properties of the arc's end points
  {
    EXPECT_TRUE(c2cPointsAreNear(c2c::Point {100.0_cm, 215.0_cm}, arc.getStart()));
    EXPECT_TRUE(c2cPointsAreNear(c2c::Point {120.0_cm, 200.0_cm}, arc.getEnd()));
  }

  // Check conversion to NURBS
  {
    // These tests assume the EllipticalArc is converted to a simple
    // endpoint-interpolating quadratic NURBS curve with no internal knots

    c2c::NURBSData curve = toNurbs(arc, c2c::LengthUnit::cm);

    // Check order
    // NOTE: The order of C2C NURBS curves is one higher than the polynomial order
    const int quadratic_c2c_order = 3;
    EXPECT_EQ(quadratic_c2c_order, static_cast<int>(curve.order));

    // Check knots
    EXPECT_EQ(6, static_cast<int>(curve.knots.size()));
    EXPECT_DOUBLE_EQ(0, curve.knots[0]);
    EXPECT_DOUBLE_EQ(0, curve.knots[1]);
    EXPECT_DOUBLE_EQ(0, curve.knots[2]);
    EXPECT_DOUBLE_EQ(1, curve.knots[3]);
    EXPECT_DOUBLE_EQ(1, curve.knots[4]);
    EXPECT_DOUBLE_EQ(1, curve.knots[5]);

    // check control points
    EXPECT_EQ(3, static_cast<int>(curve.controlPoints.size()));
    EXPECT_TRUE(c2cPointsAreNear(arc.getStart(), curve.controlPoints[0]));
    EXPECT_TRUE(c2cPointsAreNear(arc.getEnd(), curve.controlPoints[2]));
    c2c::Point midPoint {arc.getEnd().getR(), arc.getStart().getZ()};
    EXPECT_TRUE(c2cPointsAreNear(midPoint, curve.controlPoints[1]));

    // check weights
    EXPECT_NEAR(1., curve.weights[0], EPS);
    EXPECT_NEAR(std::sqrt(2) / 2., curve.weights[1], EPS);
    EXPECT_NEAR(1., curve.weights[2], EPS);
  }
}
