// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/primal/geometry/Sphere.hpp"
#include "axom/primal/geometry/Point.hpp"

#include "axom/slic.hpp"
#include "gtest/gtest.h"

namespace primal = axom::primal;

//------------------------------------------------------------------------------
// INTERNAL HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
template <typename PointType>
void check_array(const PointType& a, const PointType& b)
{
  for(int i = 0; i < PointType::dimension(); ++i)
  {
    EXPECT_DOUBLE_EQ(a[i], b[i]);
  }
}

//------------------------------------------------------------------------------
template <int NDIMS>
void check_constructor()
{
  using SphereType = primal::Sphere<double, NDIMS>;
  using PointType = primal::Point<double, NDIMS>;

  const double default_radius = 1.0;
  const double prescribed_radius = 5.0;

  const PointType default_center = PointType::zero();
  const PointType prescribed_center = PointType::ones();

  SphereType S0;
  EXPECT_DOUBLE_EQ(S0.getRadius(), default_radius);
  check_array(S0.getCenter(), default_center);

  SphereType S1(default_radius);
  EXPECT_DOUBLE_EQ(S1.getRadius(), default_radius);
  check_array(S1.getCenter(), default_center);

  SphereType S2(prescribed_radius);
  EXPECT_DOUBLE_EQ(S2.getRadius(), prescribed_radius);
  check_array(S2.getCenter(), default_center);

  SphereType S4(prescribed_center, prescribed_radius);
  EXPECT_DOUBLE_EQ(S4.getRadius(), prescribed_radius);
  check_array(S4.getCenter(), prescribed_center);

  SphereType S5(prescribed_center.data(), prescribed_radius);
  EXPECT_DOUBLE_EQ(S5.getRadius(), prescribed_radius);
  check_array(S5.getCenter(), prescribed_center);
}

//------------------------------------------------------------------------------
template <int NDIMS>
void check_signed_distance_and_orientation()
{
  using PointType = primal::Point<double, NDIMS>;
  using SphereType = primal::Sphere<double, NDIMS>;

  // test default sphere with radius 1.0 centered @ origin
  SphereType sphere;
  const double radius = sphere.getRadius();
  const auto& center = sphere.getCenter();

  // test sphere center
  double signed_distance = sphere.computeSignedDistance(center);
  EXPECT_DOUBLE_EQ(signed_distance, (-1.0) * radius);
  EXPECT_EQ(sphere.getOrientation(center), primal::ON_NEGATIVE_SIDE);

  for(int j = 0; j < NDIMS; ++j)
  {
    PointType x = center;

    // SHIFT RIGHT
    // shift right from center, but still within the sphere
    x[j] = center[j] + 0.5 * radius;
    signed_distance = sphere.computeSignedDistance(x);
    EXPECT_DOUBLE_EQ(signed_distance, (-0.5) * radius);
    EXPECT_EQ(sphere.getOrientation(x), primal::ON_NEGATIVE_SIDE);

    // shift right all the way on the sphere boundary
    x[j] = center[j] + radius;
    signed_distance = sphere.computeSignedDistance(x);
    EXPECT_DOUBLE_EQ(signed_distance, 0.0);
    EXPECT_EQ(sphere.getOrientation(x), primal::ON_BOUNDARY);

    // shift right outside sphere
    x[j] = center[j] + 2 * radius;
    signed_distance = sphere.computeSignedDistance(x);
    EXPECT_DOUBLE_EQ(signed_distance, radius);
    EXPECT_EQ(sphere.getOrientation(x), primal::ON_POSITIVE_SIDE);

    // SHIFT LEFT
    // shift left from center, but still within the sphere
    x[j] = center[j] - 0.5 * radius;
    signed_distance = sphere.computeSignedDistance(x);
    EXPECT_DOUBLE_EQ(signed_distance, (-0.5) * radius);
    EXPECT_EQ(sphere.getOrientation(x), primal::ON_NEGATIVE_SIDE);

    // shift left all the way on the sphere boundary
    x[j] = center[j] - radius;
    signed_distance = sphere.computeSignedDistance(x);
    EXPECT_DOUBLE_EQ(signed_distance, 0.0);
    EXPECT_EQ(sphere.getOrientation(x), primal::ON_BOUNDARY);

    // shift left outside sphere
    x[j] = center[j] - 2 * radius;
    signed_distance = sphere.computeSignedDistance(x);
    EXPECT_DOUBLE_EQ(signed_distance, radius);
    EXPECT_EQ(sphere.getOrientation(x), primal::ON_POSITIVE_SIDE);

  }  // END for all dimensions
}

//------------------------------------------------------------------------------
template <int NDIMS>
void check_sphere_intersection()
{
  using PointType = primal::Point<double, NDIMS>;
  using SphereType = primal::Sphere<double, NDIMS>;

  PointType center {0.0, 0.0, 0.0};

  // STEP 0: test fully overlapping
  SphereType S0;
  EXPECT_TRUE(S0.intersectsWith(S0));

  // STEP 1: test inter-penetrating
  center[0] = 0.5;
  SphereType S1(center);
  EXPECT_TRUE(S0.intersectsWith(S1));

  // STEP 2: test abutting
  center[0] = 2.0;
  SphereType S2(center);
  EXPECT_TRUE(S0.intersectsWith(S2));

  // STEP 3: test non-intersecting
  center[0] = 4.0;
  SphereType S3(center);
  EXPECT_FALSE(S0.intersectsWith(S3));
}

//------------------------------------------------------------------------------
template <int NDIMS>
void check_sphere_containment()
{
  using PointType = primal::Point<double, NDIMS>;
  using SphereType = primal::Sphere<double, NDIMS>;

  PointType center {3.0, 4.0, 0.0};
  double radius = 5.0;

  SphereType S0(center, radius);
  SphereType S1(PointType {1.5, 2.0, 0.0}, radius / 2);
  SphereType S2(S1.getCenter(), radius);

  double tol = 1e-12;

  // STEP 0: test fully containing
  EXPECT_TRUE(S0.contains(S0, -tol));
  EXPECT_TRUE(S0.contains(S1, -tol));

  // STEP 1: test barely not containing.
  EXPECT_FALSE(S0.contains(S0, tol));
  EXPECT_FALSE(S0.contains(S1, tol));

  // STEP 2: test partial containment.
  EXPECT_FALSE(S0.contains(S2));
}

//------------------------------------------------------------------------------
template <int NDIMS>
void check_copy_constructor()
{
  using PointType = primal::Point<double, NDIMS>;
  using SphereType = primal::Sphere<double, NDIMS>;

  const double MAGIC_NUM = 42;
  PointType center(MAGIC_NUM);
  double radius = MAGIC_NUM;

  SphereType s1(center, radius);
  SphereType s2(s1);

  const auto& c1 = s1.getCenter();
  const auto& c2 = s2.getCenter();
  for(int i = 0; i < NDIMS; ++i)
  {
    EXPECT_DOUBLE_EQ(c1[i], c2[i]);
    EXPECT_DOUBLE_EQ(c2[i], MAGIC_NUM);
  }

  EXPECT_DOUBLE_EQ(s1.getRadius(), s2.getRadius());
  EXPECT_DOUBLE_EQ(s2.getRadius(), MAGIC_NUM);
}

//------------------------------------------------------------------------------
template <int NDIMS>
void check_assignment_operator()
{
  using PointType = primal::Point<double, NDIMS>;
  using SphereType = primal::Sphere<double, NDIMS>;

  const double MAGIC_NUM = 42;
  PointType center(MAGIC_NUM);
  double radius = MAGIC_NUM;

  SphereType s1(center, radius);
  SphereType s2;

  /* test assignment */
  s2 = s1;

  const auto& c1 = s1.getCenter();
  const auto& c2 = s2.getCenter();
  for(int i = 0; i < NDIMS; ++i)
  {
    EXPECT_DOUBLE_EQ(c1[i], c2[i]);
    EXPECT_DOUBLE_EQ(c2[i], MAGIC_NUM);
  }

  EXPECT_DOUBLE_EQ(s1.getRadius(), s2.getRadius());
  EXPECT_DOUBLE_EQ(s2.getRadius(), MAGIC_NUM);
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(primal_sphere, constructor)
{
  check_constructor<2>();
  check_constructor<3>();
}

//------------------------------------------------------------------------------
TEST(primal_sphere, copy_constructor)
{
  check_copy_constructor<2>();
  check_copy_constructor<3>();
}

//------------------------------------------------------------------------------
TEST(primal_sphere, assignment_operator)
{
  check_assignment_operator<2>();
  check_assignment_operator<3>();
}

//------------------------------------------------------------------------------
TEST(primal_sphere, signed_distance_and_orientation)
{
  check_signed_distance_and_orientation<2>();
  check_signed_distance_and_orientation<3>();
}

//------------------------------------------------------------------------------
TEST(primal_sphere, sphere_sphere_intersection)
{
  check_sphere_intersection<2>();
  check_sphere_intersection<3>();
}

//------------------------------------------------------------------------------
TEST(primal_sphere, sphere_sphere_containment)
{
  check_sphere_containment<2>();
  check_sphere_containment<3>();
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
