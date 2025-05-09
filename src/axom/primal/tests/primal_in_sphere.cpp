// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"
#include "axom/primal/geometry/Sphere.hpp"
#include "axom/primal/geometry/OrientationResult.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/operators/in_sphere.hpp"

#include "axom/fmt.hpp"

namespace primal = axom::primal;

TEST(primal_in_sphere, test_in_sphere_2d)
{
  const int DIM = 2;
  using PointType = primal::Point<double, DIM>;
  using TriangleType = primal::Triangle<double, DIM>;

  // Define the three vertices of a triangle
  PointType p0 {0, 0};
  PointType p1 {1, 0};
  PointType p2 {0, 1};
  TriangleType tri(p0, p1, p2);

  // Test some points that are inside the circumcircle
  {
    // inside triangle
    PointType q1 {0.1, 0.1};
    EXPECT_TRUE(in_sphere(q1, p0, p1, p2));
    EXPECT_TRUE(in_sphere(q1, tri));

    // outside triangle
    PointType q2 {0.78, 0.6};
    EXPECT_TRUE(in_sphere(q2, p0, p1, p2));
    EXPECT_TRUE(in_sphere(q2, tri));
  }

  // Test some points that are on the circumcircle
  {
    PointType q1 {1, 1};
    EXPECT_FALSE(in_sphere(q1, p0, p1, p2));
    EXPECT_FALSE(in_sphere(q1, tri));

    PointType q2 {0, 1};
    EXPECT_FALSE(in_sphere(q2, p0, p1, p2));
    EXPECT_FALSE(in_sphere(q2, tri));
  }

  // Test some points that are outside the circumcircle
  {
    PointType q1 {1.1, 0};
    EXPECT_FALSE(in_sphere(q1, p0, p1, p2));
    EXPECT_FALSE(in_sphere(q1, tri));

    PointType q2 {-5, -10};
    EXPECT_FALSE(in_sphere(q2, p0, p1, p2));
    EXPECT_FALSE(in_sphere(q2, tri));
  }
}

TEST(primal_in_sphere, test_in_sphere_3d)
{
  const int DIM = 3;
  using PointType = primal::Point<double, DIM>;
  using TetrahedronType = primal::Tetrahedron<double, DIM>;

  // Define the four vertices of a tetrahedron
  PointType p0 {-1, -1, 1};
  PointType p1 {1, -1, -1};
  PointType p2 {-1, 1, -1};
  PointType p3 {1, 1, 1};
  TetrahedronType tet(p0, p1, p2, p3);

  // Test some points that are inside the circumsphere
  {
    PointType q1 {0.5, 0.5, 0.5};
    EXPECT_TRUE(in_sphere(q1, p0, p1, p2, p3));
    EXPECT_TRUE(in_sphere(q1, tet));

    PointType q2 {0., 0., 0.};
    EXPECT_TRUE(in_sphere(q2, p0, p1, p2, p3));
    EXPECT_TRUE(in_sphere(q2, tet));
  }

  // Test some points that are on the circumsphere
  {
    PointType q1 {1, 1, 1};
    EXPECT_FALSE(in_sphere(q1, p0, p1, p2, p3));
    EXPECT_FALSE(in_sphere(q1, tet));

    PointType q2 {-1, 1, 1};
    EXPECT_FALSE(in_sphere(q2, p0, p1, p2, p3));
    EXPECT_FALSE(in_sphere(q2, tet));
  }

  // Test some points that are outside the circumsphere
  {
    PointType q1 {1.1, 1, 1};
    EXPECT_FALSE(in_sphere(q1, p0, p1, p2, p3));
    EXPECT_FALSE(in_sphere(q1, tet));

    PointType q2 {-1.1, 1, 1};
    EXPECT_FALSE(in_sphere(q2, p0, p1, p2, p3));
    EXPECT_FALSE(in_sphere(q2, tet));
  }
}

TEST(primal_in_sphere, compare_to_sphere_orientation)
{
  const int DIM = 3;
  using PointType = primal::Point<double, DIM>;
  using TetrahedronType = primal::Tetrahedron<double, DIM>;
  using SphereType = TetrahedronType::SphereType;

  std::vector<PointType> queryPoints = {PointType {2, 0, 0},
                                        PointType {0.359547, 0.731502, 0.370846}};

  // define some constants for regular tetrahedron
  const double c[3] = {std::sqrt(2) / 3, std::sqrt(6) / 3, 1. / 3};
  std::vector<TetrahedronType> tets = {// A regular tetrahedron centered at the origin
                                       TetrahedronType {PointType {0, 0, 1},
                                                        PointType {2 * c[0], 0, -c[2]},
                                                        PointType {-c[0], -c[1], -c[2]},
                                                        PointType {-c[0], c[1], -c[2]}},

                                       TetrahedronType {PointType {0.942188, 0.211634, 0.395159},
                                                        PointType {0.9105, 0.114406, 0.83865},
                                                        PointType {0.855718, 0.660449, 0.911065},
                                                        PointType {0.393771, 0.572654, 0.471824}}};

  std::vector<primal::OrientationResult> orientations = {
    primal::ON_POSITIVE_SIDE,
    primal::ON_POSITIVE_SIDE,
  };

  const int numChecks = queryPoints.size();
  for(int i = 0; i < numChecks; ++i)
  {
    const auto& query = queryPoints[i];
    const auto& tet = tets[i];
    const auto& exp_orient = orientations[i];
    const SphereType sphere = tet.circumsphere();

    SLIC_INFO(
      axom::fmt::format("Testing point {} against circumsphere {} "
                        "\n\t for tet {}  with signed volume {}"
                        "\n\t signed distance to sphere boundary is {}",
                        query,
                        sphere,
                        tet,
                        tet.signedVolume(),
                        sphere.computeSignedDistance(query)));

    // check that each vertex is on the sphere boundary
    for(int i = 0; i <= DIM; ++i)
    {
      EXPECT_TRUE(sphere.getOrientation(tet[i]) == primal::ON_BOUNDARY);
    }

    // check that the orientation matches expectations
    auto res = sphere.getOrientation(query);
    EXPECT_EQ(exp_orient, res);

    // check that primal::in_sphere() matches sphere.getOrientation()
    //   in_sphere() is true when the point is inside; false on the boundary and outside
    const bool oper_in_sphere = in_sphere(query, tet);
    const bool geom_in_sphere = (res == primal::ON_NEGATIVE_SIDE);
    EXPECT_EQ(oper_in_sphere, geom_in_sphere);
  }
}

TEST(primal_in_sphere, bounding_box_in_sphere)
{
  constexpr int DIM = 2;
  using PointType = primal::Point<double, DIM>;
  using SphereType = primal::Sphere<double, DIM>;
  using BoxType = primal::BoundingBox<double, DIM>;

  SphereType circle(PointType {1.0, 2.0}, 1.3);

  // Circle contains box: All 4 corners inside
  {
    BoxType box(PointType {circle.getCenter()[0] - 0.5, circle.getCenter()[1] - 0.5},
                PointType {circle.getCenter()[0] + 0.5, circle.getCenter()[1] + 0.5});
    EXPECT_TRUE(primal::in_sphere(box, circle));
  }

  // One corner outside
  {
    BoxType box(PointType {0, 1}, PointType {circle.getCenter()[0], circle.getCenter()[1]});
    EXPECT_FALSE(primal::in_sphere(box, circle));
  }

  // Two corners outside
  {
    BoxType box(PointType {0, 1}, PointType {2, 2});
    EXPECT_FALSE(primal::in_sphere(box, circle));
  }

  // Three corners outside
  {
    BoxType box(PointType {0, 1}, PointType {2.2, 2.6});
    EXPECT_FALSE(primal::in_sphere(box, circle));
  }

  // All four corners outside
  {
    BoxType box(PointType {0, 1}, PointType {3, 3});
    EXPECT_FALSE(primal::in_sphere(box, circle));
  }
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  int result = RUN_ALL_TESTS();
  return result;
}
