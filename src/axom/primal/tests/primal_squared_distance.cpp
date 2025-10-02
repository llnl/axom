// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core/NumericLimits.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"

#include <cmath>

namespace primal = axom::primal;

//------------------------------------------------------------------------------
TEST(primal_squared_distance, array_to_array)
{
  double A[5];
  double B[5];
  int N;
  double d, expect_d;

  {
    SCOPED_TRACE("1D");
    A[0] = 3.25;
    B[0] = 6.6;
    N = 1;
    expect_d = 11.2225;
    d = primal::squared_distance(A, B, N);
    EXPECT_DOUBLE_EQ(d, expect_d);
  }

  {
    SCOPED_TRACE("2D");
    A[0] = 0;
    A[1] = 0;
    B[0] = 1.5;
    B[1] = 1.5;
    N = 2;
    expect_d = 4.5;
    d = primal::squared_distance(A, B, N);
    EXPECT_DOUBLE_EQ(d, expect_d);
  }

  {
    SCOPED_TRACE("3D");
    A[0] = 6;
    A[1] = 2.3;
    A[2] = -1;
    B[0] = 8;
    B[1] = 1.5;
    B[2] = 0;
    N = 3;
    expect_d = 5.64;
    d = primal::squared_distance(A, B, N);
    EXPECT_DOUBLE_EQ(d, expect_d);
  }
}

//------------------------------------------------------------------------------
TEST(primal_squared_distance, point_to_point)
{
  primal::Point<double, 2> A = primal::Point<double, 2>::make_point(0.0, 0.0);
  primal::Point<double, 2> B = primal::Point<double, 2>::make_point(1.5, 1.5);

  double dist = primal::squared_distance(A, B);
  EXPECT_EQ(4.5, dist);

  double dist2 = primal::squared_distance(B, A);
  EXPECT_EQ(4.5, dist2);
}

//------------------------------------------------------------------------------
TEST(primal_squared_distance, point_to_triangle)
{
  // STEP 0: Setup triangle ABC in 3D
  primal::Point<double, 3> A = primal::Point<double, 3>::make_point(0.0, 0.0, 0.0);
  primal::Point<double, 3> B = primal::Point<double, 3>::make_point(1.5, 1.5, 0.0);
  primal::Point<double, 3> C = primal::Point<double, 3>::make_point(2.5, 0.0, 0.0);
  primal::Triangle<double, 3> tri(A, B, C);

  double dist = 0.0;

  // STEP 1: Check distance of triangle nodes to the triangle
  dist = primal::squared_distance(A, tri);
  EXPECT_DOUBLE_EQ(0.0f, dist);

  dist = primal::squared_distance(B, tri);
  EXPECT_DOUBLE_EQ(0.0f, dist);

  dist = primal::squared_distance(C, tri);
  EXPECT_DOUBLE_EQ(0.0f, dist);

  // STEP 2: Check distance of edge midpoints
  primal::Point<double, 3> mid_of_AB = primal::Point<double, 3>::midpoint(A, B);
  dist = primal::squared_distance(mid_of_AB, tri);
  EXPECT_DOUBLE_EQ(0.0f, dist);

  primal::Point<double, 3> mid_of_BC = primal::Point<double, 3>::midpoint(B, C);
  dist = primal::squared_distance(mid_of_BC, tri);
  EXPECT_DOUBLE_EQ(0.0f, dist);

  primal::Point<double, 3> mid_of_AC = primal::Point<double, 3>::midpoint(A, C);
  dist = primal::squared_distance(mid_of_AC, tri);
  EXPECT_DOUBLE_EQ(0.0f, dist);

  // STEP 3: Check distance of point q, inside ABC and orthogonally
  // projected along z+
  primal::Point<double, 3> q = primal::Point<double, 3>::make_point(1.5, 0.5, 0.5);
  dist = primal::squared_distance(q, tri);
  EXPECT_DOUBLE_EQ(0.25f, dist);

  // STEP 4: Check distance of point xA, within the voronoi region of A
  primal::Point<double, 3> xA = primal::Point<double, 3>::make_point(-1.0, -1.0, 0.0);
  dist = primal::squared_distance(xA, tri);
  EXPECT_DOUBLE_EQ(primal::squared_distance(xA, A), dist);

  // STEP 5: Check distance of point xB, within the voronoi region B
  primal::Point<double, 3> xB = primal::Point<double, 3>::make_point(1.5, 2.0, 0.0);
  dist = primal::squared_distance(xB, tri);
  EXPECT_DOUBLE_EQ(primal::squared_distance(xB, B), dist);

  // STEP 6: Check distance of point xC within the voronoi region C
  primal::Point<double, 3> xC = primal::Point<double, 3>::make_point(3.0, -1.0, 0.5);
  dist = primal::squared_distance(xC, tri);
  EXPECT_DOUBLE_EQ(primal::squared_distance(xC, C), dist);
}

//------------------------------------------------------------------------------
TEST(primal_squared_distance, point_to_segment)
{
  primal::Point<double, 2> A(0.0);
  primal::Point<double, 2> B(1.0, 1);
  primal::Segment<double, 2> S(A, B);

  primal::Point<double, 2> q0 = primal::Point<double, 2>::make_point(0.5, 10);
  double dist = primal::squared_distance(q0, S);
  EXPECT_DOUBLE_EQ(100.0f, dist);

  // STEP 0: check source point
  dist = primal::squared_distance(A, S);
  EXPECT_DOUBLE_EQ(0.0f, dist);

  // STEP 1: check target point
  dist = primal::squared_distance(B, S);
  EXPECT_DOUBLE_EQ(0.0f, dist);

  // STEP 2: check midpoint
  primal::Point<double, 2> Q = primal::Point<double, 2>::midpoint(A, B);
  dist = primal::squared_distance(Q, S);
  EXPECT_DOUBLE_EQ(0.0f, dist);

  // STEP 3: check projection to target point
  primal::Point<double, 2> q1(2.0);
  dist = primal::squared_distance(q1, S);
  EXPECT_DOUBLE_EQ(primal::squared_distance(q1, B), dist);

  // STEP 4: check projection source point
  primal::Point<double, 2> q2(-1.0);
  dist = primal::squared_distance(q2, S);
  EXPECT_DOUBLE_EQ(primal::squared_distance(q2, A), dist);
}

//------------------------------------------------------------------------------
TEST(primal_squared_distance, point_to_bbox)
{
  constexpr int DIM = 3;
  constexpr double EPS = 1e-12;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QBBox = primal::BoundingBox<CoordType, DIM>;

  const QBBox cube(QPoint {-1, -1, -1}, QPoint {1, 1, 1});
  const QBBox empty;

  for(int i = -1; i <= 1; ++i)
  {
    for(int j = -1; j <= 1; ++j)
    {
      for(int k = -1; k <= 1; ++k)
      {
        const QPoint pt {i * 3., j * 3., k * 3.};
        if(i == 0 && j == 0 && k == 0)
        {
          EXPECT_NEAR(0., primal::squared_distance(pt, cube), EPS);
        }
        else
        {
          // if a coordinate is outside the bounding box,
          // it adds 4 == (3-1)^2 units to the squared distance
          const double sqsum = (i == 0 ? 0 : 4) + (j == 0 ? 0 : 4) + (k == 0 ? 0 : 4);
          EXPECT_NEAR(sqsum, primal::squared_distance(pt, cube), EPS);
        }

        EXPECT_EQ(axom::numeric_limits<double>::max(), squared_distance(pt, empty));
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_squared_distance, bbox_to_bbox)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QVector = primal::Vector<CoordType, DIM>;
  using QBBox = primal::BoundingBox<CoordType, DIM>;

  const QBBox middle(QPoint {0.0, 0.0, 0.0}, QPoint {1.0, 1.0, 1.0});
  const QBBox north(QPoint {0., 10., 0.}, QPoint {1., 11., 1.});
  const QBBox east(QPoint {10., 0., 0.}, QPoint {11., 1., 1.});
  const QBBox northeast(QPoint {10., 10., 0.}, QPoint {11., 11., 1.});
  const QBBox northeastup(QPoint {10., 10., 10.}, QPoint {11., 11., 11.});
  QBBox touching(middle);
  touching.shift(QVector::make_vector(0.9, 0.99, 0.999));

  EXPECT_DOUBLE_EQ(squared_distance(middle, north), 81.);
  EXPECT_DOUBLE_EQ(squared_distance(middle, east), 81.);
  EXPECT_DOUBLE_EQ(squared_distance(middle, northeast), 162.);
  EXPECT_DOUBLE_EQ(squared_distance(middle, northeastup), 243.);
  EXPECT_DOUBLE_EQ(squared_distance(middle, touching), 0.);

  // check that squared distances for empty/invalid boxes is max double
  const QBBox empty;
  EXPECT_EQ(axom::numeric_limits<double>::max(), squared_distance(middle, empty));
  EXPECT_EQ(axom::numeric_limits<double>::max(), squared_distance(empty, middle));
  EXPECT_EQ(axom::numeric_limits<double>::max(), squared_distance(empty, empty));
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
