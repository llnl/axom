// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Plane.hpp"
#include "axom/primal/operators/orientation.hpp"

TEST(primal_orientation, orient3D)
{
  namespace primal = axom::primal;

  using Point3 = primal::Point<double, 3>;
  using Vector3 = primal::Vector<double, 3>;
  using BaryPoint = primal::Point<double, 3>;
  using Tri = primal::Triangle<double, 3>;
  using Plane3 = primal::Plane<double, 3>;

  // Test a few triangles
  for(const auto& tri : {Tri(Point3 {1, 0, 0}, Point3 {0, 1, 0}, Point3 {0, 0, 1}),
                         Tri(Point3 {0, 0, 0}, Point3 {1, 0, 0}, Point3 {1, 1, 0}),
                         Tri(Point3 {0, 0, 0}, Point3 {1.5, 1.5, 0}, Point3 {2.5, 0, 0}),
                         Tri(Point3 {-2, -3, -4}, Point3 {3, 4, 5}, Point3 {-10, 20, -30})})
  {
    // Specify test points relative to triangle via barycentric coordinates
    for(const auto& baryPoint : {BaryPoint {1, 0, 0},
                                 BaryPoint {0, 1, 0},
                                 BaryPoint {0, 0, 1},
                                 BaryPoint {.5, .5, 0},
                                 BaryPoint {.5, 0, .5},
                                 BaryPoint {0, .5, .5},
                                 BaryPoint {1. / 3, 1. / 3, 1. / 3},
                                 BaryPoint {1.25, -12, 11.75}})
    {
      // Sum of barycentric coords should be 1
      EXPECT_DOUBLE_EQ(1., baryPoint[0] + baryPoint[1] + baryPoint[2]);

      // Convert baryPoint to physical space and get normal
      auto phys = tri.baryToPhysical(baryPoint);
      auto normal = tri.normal();

      // check orientation of a few offset points
      // Without offset, the point should be on the same plane
      EXPECT_EQ(primal::ON_BOUNDARY, primal::orientation(phys, tri));

      // Offset along negative normal should have negative orientation
      EXPECT_EQ(primal::ON_NEGATIVE_SIDE, primal::orientation(phys - normal, tri));
      EXPECT_EQ(primal::ON_NEGATIVE_SIDE, primal::orientation(phys - 0.25 * normal, tri));

      // Offset along positive normal should have positive orientation
      EXPECT_EQ(primal::ON_POSITIVE_SIDE, primal::orientation(phys + normal, tri));
      EXPECT_EQ(primal::ON_POSITIVE_SIDE, primal::orientation(phys + 0.25 * normal, tri));

      // check that orientation is equivalent to half-space definition
      {
        EXPECT_LT(normal.dot(Vector3(phys, phys - normal)), 0.);
        EXPECT_LT(normal.dot(Vector3(phys, phys - 0.25 * normal)), 0.);
        EXPECT_GT(normal.dot(Vector3(phys, phys + normal)), 0.);
        EXPECT_GT(normal.dot(Vector3(phys, phys + 0.25 * normal)), 0.);
      }

      // check equivalence to Plane orientation
      {
        const auto plane = Plane3(normal, tri[0]);
        for(const auto& pt :
            {phys, phys - normal, phys - 2.5 * normal, phys + normal, phys + 0.25 * normal})
        {
          EXPECT_EQ(plane.getOrientation(pt), primal::orientation(pt, tri));
        }
      }

      // Orientation on flipped triangle should be flipped
      {
        const auto flipped_triangle = Tri(tri[0], tri[2], tri[1]);
        EXPECT_EQ(primal::ON_BOUNDARY, primal::orientation(phys, flipped_triangle));
        EXPECT_EQ(primal::ON_POSITIVE_SIDE, primal::orientation(phys - normal, flipped_triangle));
        EXPECT_EQ(primal::ON_NEGATIVE_SIDE, primal::orientation(phys + normal, flipped_triangle));
      }

      // check overload with explicit tolerances
      {
        constexpr double TOL = 1e-2;
        constexpr double smallOff = 1e-5;
        constexpr double largeOff = 1e-1;
        const auto unitNormal = normal.unitVector();

        EXPECT_EQ(primal::ON_BOUNDARY, primal::orientation(phys - smallOff * unitNormal, tri, TOL));
        EXPECT_EQ(primal::ON_NEGATIVE_SIDE,
                  primal::orientation(phys - largeOff * unitNormal, tri, TOL));

        EXPECT_EQ(primal::ON_BOUNDARY, primal::orientation(phys + smallOff * unitNormal, tri, TOL));

        EXPECT_EQ(primal::ON_POSITIVE_SIDE,
                  primal::orientation(phys + largeOff * unitNormal, tri, TOL));
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_orientation, orient2D)
{
  namespace primal = axom::primal;

  using Point2 = primal::Point<double, 2>;
  using Vector2 = primal::Vector<double, 2>;
  using Segment = primal::Segment<double, 2>;
  using Plane2 = primal::Plane<double, 2>;

  // Test a few segments
  for(const auto& seg : {Segment(Point2 {0, 0}, Point2 {1, 1}),
                         Segment(Point2 {1, 0}, Point2 {0, 1}),
                         Segment(Point2 {0, 0}, Point2 {1, 0}),
                         Segment(Point2 {0, 2.5}, Point2 {2.5, 0}),
                         Segment(Point2 {-1, 0}, Point2 {0, -1}),
                         Segment(Point2 {-2, -3}, Point2 {3, 4})})
  {
    // Specify test points on segment
    for(const auto& phys :
        {seg.at(0.), seg.at(1.), seg.at(0.5), seg.at(0.33), seg.at(-2.53), seg.at(3.14)})
    {
      auto normal = seg.normal();

      // check orientation of a few offset points
      // Without offset, the point should be on the same plane
      EXPECT_EQ(primal::ON_BOUNDARY, primal::orientation(phys, seg));

      // Offset along negative normal should have negative orientation
      EXPECT_EQ(primal::ON_NEGATIVE_SIDE, primal::orientation(phys - normal, seg));
      EXPECT_EQ(primal::ON_NEGATIVE_SIDE, primal::orientation(phys - 0.25 * normal, seg));

      // Offset along positive normal should have positive orientation
      EXPECT_EQ(primal::ON_POSITIVE_SIDE, primal::orientation(phys + normal, seg));
      EXPECT_EQ(primal::ON_POSITIVE_SIDE, primal::orientation(phys + 0.25 * normal, seg));

      // check that orientation is equivalent to half-space definition
      {
        EXPECT_LT(normal.dot(Vector2(phys, phys - normal)), 0.);
        EXPECT_LT(normal.dot(Vector2(phys, phys - 0.25 * normal)), 0.);
        EXPECT_GT(normal.dot(Vector2(phys, phys + normal)), 0.);
        EXPECT_GT(normal.dot(Vector2(phys, phys + 0.25 * normal)), 0.);
      }

      // check equivalence to Plane orientation
      {
        const auto plane = Plane2(normal, seg[0]);
        for(const auto& pt :
            {phys, phys - normal, phys - 2.5 * normal, phys + normal, phys + 0.25 * normal})
        {
          EXPECT_EQ(plane.getOrientation(pt), primal::orientation(pt, seg));
        }
      }

      // Orientation on flipped triangle should be flipped
      {
        const auto flipped_segment = Segment(seg[1], seg[0]);
        EXPECT_EQ(primal::ON_BOUNDARY, primal::orientation(phys, flipped_segment));
        EXPECT_EQ(primal::ON_POSITIVE_SIDE, primal::orientation(phys - normal, flipped_segment));
        EXPECT_EQ(primal::ON_NEGATIVE_SIDE, primal::orientation(phys + normal, flipped_segment));
      }

      // check overload with explicit tolerances
      {
        constexpr double TOL = 1e-2;
        constexpr double smallOff = 1e-5;
        constexpr double largeOff = 1e-1;
        const auto unitNormal = normal.unitVector();

        EXPECT_EQ(primal::ON_BOUNDARY, primal::orientation(phys - smallOff * unitNormal, seg, TOL));
        EXPECT_EQ(primal::ON_NEGATIVE_SIDE,
                  primal::orientation(phys - largeOff * unitNormal, seg, TOL));

        EXPECT_EQ(primal::ON_BOUNDARY, primal::orientation(phys + smallOff * unitNormal, seg, TOL));
        EXPECT_EQ(primal::ON_POSITIVE_SIDE,
                  primal::orientation(phys + largeOff * unitNormal, seg, TOL));
      }
    }
  }
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
