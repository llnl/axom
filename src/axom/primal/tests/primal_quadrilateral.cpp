// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"

#include "gtest/gtest.h"

namespace primal = axom::primal;

TEST(primal_quadrilateral, construct_from_points)
{
  constexpr int DIM = 2;
  using CoordinateType = double;
  using PointType = primal::Point<CoordinateType, DIM>;
  using QuadrilateralType = primal::Quadrilateral<CoordinateType, DIM>;

  PointType A {-1.0, 0.1};
  PointType B {-0.1, 0.5};
  PointType C {2.0, 0.0};
  PointType D {0.0, 1.0};

  QuadrilateralType quad {A, B, C, D};

  EXPECT_EQ(A, quad[0]);
  EXPECT_EQ(B, quad[1]);
  EXPECT_EQ(C, quad[2]);
  EXPECT_EQ(D, quad[3]);
}

TEST(primal_quadrilateral, construct_from_array_view)
{
  constexpr int DIM = 3;
  using CoordinateType = double;
  using PointType = primal::Point<CoordinateType, DIM>;
  using QuadrilateralType = primal::Quadrilateral<CoordinateType, DIM>;

  // This is a concave quadrilateral
  PointType quadPoints[QuadrilateralType::NUM_QUAD_VERTS] {PointType {-1.0, 0.1, 0.0},
                                                           PointType {-0.1, 0.5, 0.0},
                                                           PointType {2.0, 0.0, 1.0},
                                                           PointType {0.0, 1.0, 1.0}};

  QuadrilateralType quad {axom::ArrayView<PointType> {quadPoints, QuadrilateralType::NUM_QUAD_VERTS}};

  for(int i = 0; i < QuadrilateralType::NUM_QUAD_VERTS; ++i)
  {
    EXPECT_EQ(quadPoints[i], quad[i]);
  }
}

TEST(primal_quadrilateral, area_2D)
{
  constexpr int DIM = 2;
  constexpr double EPS = 1e-12;
  using CoordinateType = double;
  using PointType = primal::Point<CoordinateType, DIM>;
  using QuadrilateralType = primal::Quadrilateral<CoordinateType, DIM>;

  // This is a concave quadrilateral
  QuadrilateralType quad {PointType {-1.0, 0.1},
                          PointType {-0.1, 0.5},
                          PointType {2.0, 0.0},
                          PointType {0.0, 1.0}};

  const CoordinateType signedArea = quad.signedArea();
  const CoordinateType area = quad.area();

  EXPECT_NEAR(0.755, signedArea, EPS);
  EXPECT_NEAR(0.755, area, EPS);
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  int result = RUN_ALL_TESTS();
  return result;
}
