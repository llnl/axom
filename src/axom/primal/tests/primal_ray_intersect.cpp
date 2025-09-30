// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/primal/operators/detail/intersect_ray_impl.hpp"

#include "axom/core/numerics/Matrix.hpp"
#include "axom/core/numerics/matvecops.hpp"
#include "axom/slic.hpp"

#include "gtest/gtest.h"

// C/C++ includes
#include <cmath>

// namespace aliases
namespace primal = axom::primal;
namespace numerics = axom::numerics;

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(primal_ray_intersect, ray_aabb_intersection_2D)
{
  constexpr int N = 3;
  constexpr double LO = 0.0f;
  constexpr double HI = 1.0f;

  using namespace primal::detail;

// Defines the test box:
#define TEST_BOX2D LO, HI, LO, HI

  double x[N];
  numerics::linspace(LO, HI, x, N);

  // test BOTTOM (-y)
  {
    const double n0[] = {0.0, 1.0};
    const double y0 = -1.0f;
    for(int i = 0; i < N; ++i)
    {
      double tmin = 0.0;
      double tmax = 0.0;

      const double sx = x[i];
      const double sy = y0;
      const double nx = n0[0];
      const double ny = n0[1];
      EXPECT_TRUE(intersect_ray(sx, nx, sy, ny, TEST_BOX2D, tmin, tmax));
      EXPECT_FALSE(intersect_ray(sx, -nx, sy, -ny, TEST_BOX2D, tmin, tmax));
    }
  }

  // test RIGHT (+x)
  {
    const double n0[] = {-1.0, 0.0};
    const double x0 = 2.0f;
    for(int i = 0; i < N; ++i)
    {
      double tmin = 0.0;
      double tmax = 0.0;

      const double sx = x0;
      const double sy = x[i];
      const double nx = n0[0];
      const double ny = n0[1];
      EXPECT_TRUE(intersect_ray(sx, nx, sy, ny, TEST_BOX2D, tmin, tmax));
      EXPECT_FALSE(intersect_ray(sx, -nx, sy, -ny, TEST_BOX2D, tmin, tmax));
    }
  }

  // test TOP (+y)
  {
    const double n0[] = {0.0, -1.0};
    const double y0 = 2.0f;
    for(int i = 0; i < N; ++i)
    {
      double tmin = 0.0;
      double tmax = 0.0;

      const double sx = x[i];
      const double sy = y0;
      const double nx = n0[0];
      const double ny = n0[1];
      EXPECT_TRUE(intersect_ray(sx, nx, sy, ny, TEST_BOX2D, tmin, tmax));
      EXPECT_FALSE(intersect_ray(sx, -nx, sy, -ny, TEST_BOX2D, tmin, tmax));
    }
  }

  // test LEFT (-x)
  {
    const double n0[] = {1.0, 0.0};
    const double x0 = -1.0f;
    for(int i = 0; i < N; ++i)
    {
      double tmin = 0.0;
      double tmax = 0.0;

      const double sx = x0;
      const double sy = x[i];
      const double nx = n0[0];
      const double ny = n0[1];
      EXPECT_TRUE(intersect_ray(sx, nx, sy, ny, TEST_BOX2D, tmin, tmax));
      EXPECT_FALSE(intersect_ray(sx, -nx, sy, -ny, TEST_BOX2D, tmin, tmax));
    }
  }

  // Test a bunch of rays emitted from the box center
  constexpr int NUM_ANGLES = 20;
  double angles[NUM_ANGLES];
  numerics::linspace(0.0, 360.0, angles, NUM_ANGLES);

  numerics::Matrix<double> rotation_matrix(2, 2);

  constexpr double PI_OVER_180 = M_PI / 180.0;
  const double xc[] = {0.5, 0.5};
  const double e1[] = {1.0, 0.0};
  for(int i = 0; i < NUM_ANGLES; ++i)
  {
    const double t = angles[i] * PI_OVER_180;
    const double cost = cos(t);
    const double sint = sin(t);

    rotation_matrix(0, 0) = cost;
    rotation_matrix(0, 1) = -sint;
    rotation_matrix(1, 0) = sint;
    rotation_matrix(1, 1) = cost;

    double n[2];
    numerics::matrix_vector_multiply(rotation_matrix, e1, n);

    double tmin = 0.0;
    double tmax = 0.0;
    const double x0 = xc[0];
    const double y0 = xc[1];
    const double nx = n[0];
    const double ny = n[1];
    EXPECT_TRUE(intersect_ray(x0, nx, y0, ny, TEST_BOX2D, tmin, tmax));
    EXPECT_TRUE(intersect_ray(x0, -nx, y0, -ny, TEST_BOX2D, tmin, tmax));
  }

#undef TEST_BOX2D
}

//------------------------------------------------------------------------------
TEST(primal_ray_intersect, ray_aabb_intersection_3D)
{
  constexpr int N = 3;
  constexpr double LO = 0.0f;
  constexpr double HI = 1.0f;

  using namespace primal::detail;

// Defines the test box:
#define TEST_BOX3D LO, HI, LO, HI, LO, HI

  double x[N];
  numerics::linspace(LO, HI, x, N);

  // test BOTTOM (-z)
  {
    const double n0[] = {0.0, 0.0, 1.0};
    const double z0 = -1.0f;

    for(int i = 0; i < N; ++i)
    {
      for(int j = 0; j < N; ++j)
      {
        double tmin = 0.0;
        double tmax = 0.0;

        const double sx = x[i];
        const double sy = x[j];
        const double sz = z0;
        const double nx = n0[0];
        const double ny = n0[1];
        const double nz = n0[2];
        EXPECT_TRUE(intersect_ray(sx, nx, sy, ny, sz, nz, TEST_BOX3D, tmin, tmax));
        EXPECT_FALSE(intersect_ray(sx, -nx, sy, -ny, sz, -nz, TEST_BOX3D, tmin, tmax));
      }  // END for all j
    }  // END for all i
  }

  // test TOP (+z)
  {
    const double n0[] = {0.0, 0.0, -1.0};
    const double z0 = 2.0f;

    for(int i = 0; i < N; ++i)
    {
      for(int j = 0; j < N; ++j)
      {
        double tmin = 0.0;
        double tmax = 0.0;

        const double sx = x[i];
        const double sy = x[j];
        const double sz = z0;
        const double nx = n0[0];
        const double ny = n0[1];
        const double nz = n0[2];
        EXPECT_TRUE(intersect_ray(sx, nx, sy, ny, sz, nz, TEST_BOX3D, tmin, tmax));
        EXPECT_FALSE(intersect_ray(sx, -nx, sy, -ny, sz, -nz, TEST_BOX3D, tmin, tmax));
      }  // END for all j
    }  // END for all i
  }

  // test LEFT (-y)
  {
    const double n0[] = {0.0, 1.0, 0.0};
    const double y0 = -1.0f;

    for(int i = 0; i < N; ++i)
    {
      for(int j = 0; j < N; ++j)
      {
        double tmin = 0.0;
        double tmax = 0.0;

        const double sx = x[i];
        const double sy = y0;
        const double sz = x[j];
        const double nx = n0[0];
        const double ny = n0[1];
        const double nz = n0[2];
        EXPECT_TRUE(intersect_ray(sx, nx, sy, ny, sz, nz, TEST_BOX3D, tmin, tmax));
        EXPECT_FALSE(intersect_ray(sx, -nx, sy, -ny, sz, -nz, TEST_BOX3D, tmin, tmax));
      }  // END for all j
    }  // END for all i
  }

  // test RIGHT (+y)
  {
    const double n0[] = {0.0, -1.0, 0.0};
    const double y0 = 2.0f;

    for(int i = 0; i < N; ++i)
    {
      for(int j = 0; j < N; ++j)
      {
        double tmin = 0.0;
        double tmax = 0.0;

        const double sx = x[i];
        const double sy = y0;
        const double sz = x[j];
        const double nx = n0[0];
        const double ny = n0[1];
        const double nz = n0[2];
        EXPECT_TRUE(intersect_ray(sx, nx, sy, ny, sz, nz, TEST_BOX3D, tmin, tmax));
        EXPECT_FALSE(intersect_ray(sx, -nx, sy, -ny, sz, -nz, TEST_BOX3D, tmin, tmax));
      }  // END for all j
    }  // END for all i
  }

  // test BACK (-x)
  {
    const double n0[] = {1.0, 0.0, 0.0};
    const double x0 = -1.0f;

    for(int i = 0; i < N; ++i)
    {
      for(int j = 0; j < N; ++j)
      {
        double tmin = 0.0;
        double tmax = 0.0;

        const double sx = x0;
        const double sy = x[i];
        const double sz = x[j];
        const double nx = n0[0];
        const double ny = n0[1];
        const double nz = n0[2];
        EXPECT_TRUE(intersect_ray(sx, nx, sy, ny, sz, nz, TEST_BOX3D, tmin, tmax));
        EXPECT_FALSE(intersect_ray(sx, -nx, sy, -ny, sz, -nz, TEST_BOX3D, tmin, tmax));
      }  // END for all j
    }  // END for all i
  }

  // test FRONT (+x)
  {
    const double n0[] = {-1.0, 0.0, 0.0};
    const double x0 = 2.0f;

    for(int i = 0; i < N; ++i)
    {
      for(int j = 0; j < N; ++j)
      {
        double tmin = 0.0;
        double tmax = 0.0;

        const double sx = x0;
        const double sy = x[i];
        const double sz = x[j];
        const double nx = n0[0];
        const double ny = n0[1];
        const double nz = n0[2];
        EXPECT_TRUE(intersect_ray(sx, nx, sy, ny, sz, nz, TEST_BOX3D, tmin, tmax));
        EXPECT_FALSE(intersect_ray(sx, -nx, sy, -ny, sz, -nz, TEST_BOX3D, tmin, tmax));
      }  // END for all j
    }  // END for all i
  }

  // Test a bunch of rays emitted from the box center
  constexpr int NUM_ANGLES = 20;
  double angles[NUM_ANGLES];
  numerics::linspace(0.0, 360.0, angles, NUM_ANGLES);

  numerics::Matrix<double> Rx(3, 3);
  Rx(0, 0) = 1.0;
  numerics::Matrix<double> Ry(3, 3);
  Ry(1, 1) = 1.0;
  numerics::Matrix<double> Rz(3, 3);
  Rz(2, 2) = 1.0;

  constexpr double PI_OVER_180 = M_PI / 180.0;
  const double xc[] = {0.5, 0.5, 0.5};
  const double e1[] = {1.0, 0.0, 0.0};
  const double e2[] = {0.0, 1.0, 0.0};
  const double e3[] = {0.0, 0.0, 1.0};
  double n[3];
  for(int i = 0; i < NUM_ANGLES; ++i)
  {
    const double t = angles[i] * PI_OVER_180;
    const double cost = cos(t);
    const double sint = sin(t);

    double tmin = 0.0;
    double tmax = 0.0;

    double nx = 0.0;
    double ny = 0.0;
    double nz = 0.0;

    const double x0 = xc[0];
    const double y0 = xc[1];
    const double z0 = xc[2];

    // Update Rx
    Rx(1, 1) = cost;
    Rx(1, 2) = -sint;
    Rx(2, 1) = sint;
    Rx(2, 2) = cost;

    numerics::matrix_vector_multiply(Rx, e1, n);
    nx = n[0];
    ny = n[1];
    nz = n[2];
    EXPECT_TRUE(intersect_ray(x0, nx, y0, ny, z0, nz, TEST_BOX3D, tmin, tmax));

    // Update Ry
    Ry(0, 0) = cost;
    Ry(0, 2) = sint;
    Ry(2, 0) = -sint;
    Ry(2, 2) = cost;

    numerics::matrix_vector_multiply(Ry, e2, n);
    nx = n[0];
    ny = n[1];
    nz = n[2];
    EXPECT_TRUE(intersect_ray(x0, nx, y0, ny, z0, nz, TEST_BOX3D, tmin, tmax));

    // Update Rz
    Rz(0, 0) = cost;
    Rz(0, 1) = -sint;
    Rz(1, 0) = sint;
    Rz(1, 1) = cost;

    numerics::matrix_vector_multiply(Rz, e3, n);
    nx = n[0];
    ny = n[1];
    nz = n[2];
    EXPECT_TRUE(intersect_ray(x0, nx, y0, ny, z0, nz, TEST_BOX3D, tmin, tmax));
  }

#undef TEST_BOX3D
}

TEST(primal_ray_intersect, ray_segment_edge_cases)
{
  using namespace primal::detail;

  using PointType = primal::Point<double, 2>;
  using SegmentType = primal::Segment<double, 2>;
  using VectorType = primal::Vector<double, 2>;
  using RayType = primal::Ray<double, 2>;

  RayType the_ray(PointType({0.0, 0.1}), VectorType({1, 2}));
  double ray_param = -1, seg_param = -1;
  const double EPS = 1e-8;

  // Check segments oriented in same direction as ray
  //  and in opposite direction
  SegmentType intersecting_segs[] = {SegmentType(PointType({1, 2.1}), PointType({2, 4.1})),
                                     SegmentType(PointType({0, 0.1}), PointType({1, 2.1})),
                                     SegmentType(PointType({-1, -1.9}), PointType({1, 2.1})),
                                     SegmentType(PointType({-1, -1.9}), PointType({0, 0.1})),
                                     SegmentType(PointType({2, 4.1}), PointType({1, 2.1})),
                                     SegmentType(PointType({1, 2.1}), PointType({0, 0.1})),
                                     SegmentType(PointType({1, 2.1}), PointType({-1, -1.9})),
                                     SegmentType(PointType({0, 0.1}), PointType({-1, -1.9})),
                                     SegmentType(PointType({1, 2.1}), PointType({1, 2.1}))};

  for(auto& seg : intersecting_segs)
  {
    EXPECT_TRUE(intersect_ray(the_ray, seg, ray_param, seg_param, EPS));
    EXPECT_NEAR(the_ray.at(ray_param)[0], seg.at(seg_param)[0], EPS);
    EXPECT_NEAR(the_ray.at(ray_param)[1], seg.at(seg_param)[1], EPS);
  }

  // Check that parallel, but not colinear segments do not intersect
  SegmentType not_intersecting_segs[] = {SegmentType(PointType({-2, -3.9}), PointType({-1, -1.9})),
                                         SegmentType(PointType({-2, -2.9}), PointType({-1, -0.9})),
                                         SegmentType(PointType({-2, -2.9}), PointType({1, 3.1})),
                                         SegmentType(PointType({-1, -1.9}), PointType({-2, -3.9})),
                                         SegmentType(PointType({-1, -0.9}), PointType({-2, -2.9})),
                                         SegmentType(PointType({1, 3.1}), PointType({-2, -2.9}))};

  for(auto& seg : not_intersecting_segs)
  {
    EXPECT_FALSE(intersect_ray(the_ray, seg, ray_param, seg_param, EPS));
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
