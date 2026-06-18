// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
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
TEST(primal_orientation, determinant_helpers)
{
  namespace primal = axom::primal;

  using Point2 = primal::Point<double, 2>;
  using Point3 = primal::Point<double, 3>;

  EXPECT_GT(primal::orientation_determinant(Point2 {0., 0.}, Point2 {1., 0.}, Point2 {0., 1.}), 0.);
  EXPECT_LT(primal::orientation_determinant(Point2 {0., 0.}, Point2 {0., 1.}, Point2 {1., 0.}), 0.);

  EXPECT_GT(primal::orientation_determinant(Point3 {0., 0., 0.},
                                            Point3 {1., 0., 0.},
                                            Point3 {0., 1., 0.},
                                            Point3 {0., 0., 1.}),
            0.);
  EXPECT_LT(primal::orientation_determinant(Point3 {0., 0., 0.},
                                            Point3 {0., 1., 0.},
                                            Point3 {1., 0., 0.},
                                            Point3 {0., 0., 1.}),
            0.);
}

//------------------------------------------------------------------------------
// Characterizes the double-precision orientation predicate near degeneracy.
//
// These tests document where the bare double-precision sign is and is not
// reliable (see the precision note in detail/predicate_determinants.hpp).
// They intentionally pin the current behavior so a future exact adaptive
// predicate can be swapped in and validated against the same cases.
//------------------------------------------------------------------------------
TEST(primal_orientation, degeneracy_characterization_2d)
{
  namespace primal = axom::primal;
  using Point2 = primal::Point<double, 2>;
  using Segment = primal::Segment<double, 2>;

  // Exactly collinear points: the orientation determinant is exactly zero,
  // so the classifier reports ON_BOUNDARY for any non-negative tolerance.
  {
    Point2 a {0., 0.}, b {1., 1.}, c {2., 2.};  // all on the line y = x
    const double det = primal::orientation_determinant(a, b, c);
    EXPECT_EQ(det, 0.);
    Segment seg2(a, b);
    EXPECT_EQ(primal::orientation(c, seg2), primal::ON_BOUNDARY);
  }

  // Clearly non-degenerate points classify with a definite (non-boundary) sign.
  {
    Point2 a {0., 0.}, b {1., 0.};
    Segment seg2(a, b);
    EXPECT_NE(primal::orientation(Point2 {0.5, 0.5}, seg2), primal::ON_BOUNDARY);
    EXPECT_NE(primal::orientation(Point2 {0.5, -0.5}, seg2), primal::ON_BOUNDARY);
    // and the two sides receive opposite classifications
    EXPECT_NE(primal::orientation(Point2 {0.5, 0.5}, seg2),
              primal::orientation(Point2 {0.5, -0.5}, seg2));
  }

  // The determinant is not normalized: for a fixed shape, translating the points
  // far from the origin and scaling them up inflates the determinant magnitude.
  // This is why a fixed absolute EPS is inadequate for large coordinates
  // and why callers (e.g. quest::Delaunay) use a scale-aware tolerance.
  {
    Point2 a {0., 0.}, b {1., 0.}, c {0.5, 1.};
    const double small_det = primal::orientation_determinant(a, b, c);

    const double S = 1e6;
    Point2 A {S * a[0], S * a[1]}, B {S * b[0], S * b[1]}, C {S * c[0], S * c[1]};
    const double big_det = primal::orientation_determinant(A, B, C);

    // 2D orientation determinant has degree 2 in the coordinates,
    // so scaling by S multiplies it by ~S^2.
    EXPECT_GT(axom::utilities::abs(big_det), axom::utilities::abs(small_det));
    EXPECT_NEAR(big_det, small_det * S * S, axom::utilities::abs(big_det) * 1e-9);
  }
}

TEST(primal_orientation, degeneracy_characterization_3d)
{
  namespace primal = axom::primal;
  using Point3 = primal::Point<double, 3>;
  using Tri = primal::Triangle<double, 3>;

  // Exactly coplanar query point: determinant is exactly zero -> ON_BOUNDARY.
  {
    Tri tri(Point3 {0., 0., 0.}, Point3 {1., 0., 0.}, Point3 {0., 1., 0.});  // z = 0 plane
    Point3 coplanar {0.25, 0.25, 0.};
    const double det = primal::orientation_determinant(coplanar, tri[0], tri[1], tri[2]);
    EXPECT_EQ(det, 0.);
    EXPECT_EQ(primal::orientation(coplanar, tri), primal::ON_BOUNDARY);
  }

  // Points clearly above / below the plane get opposite, non-boundary signs.
  {
    Tri tri(Point3 {0., 0., 0.}, Point3 {1., 0., 0.}, Point3 {0., 1., 0.});
    const int above = primal::orientation(Point3 {0.25, 0.25, 1.}, tri);
    const int below = primal::orientation(Point3 {0.25, 0.25, -1.}, tri);
    EXPECT_NE(above, primal::ON_BOUNDARY);
    EXPECT_NE(below, primal::ON_BOUNDARY);
    EXPECT_NE(above, below);
  }

  // Degree-3 scaling of the 3D orientation determinant: scaling coordinates by S
  // multiplies the determinant by ~S^3
  {
    Point3 q {0.25, 0.25, 0.5}, p0 {0., 0., 0.}, p1 {1., 0., 0.}, p2 {0., 1., 0.};
    const double small_det = primal::orientation_determinant(q, p0, p1, p2);
    ASSERT_NE(small_det, 0.);

    const double S = 1.0e4;
    auto scaled = [S](const Point3& p) { return Point3 {S * p[0], S * p[1], S * p[2]}; };
    const double big_det =
      primal::orientation_determinant(scaled(q), scaled(p0), scaled(p1), scaled(p2));

    EXPECT_NEAR(big_det, small_det * S * S * S, axom::utilities::abs(big_det) * 1e-9);
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
