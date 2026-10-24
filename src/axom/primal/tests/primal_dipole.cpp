// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/core.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/geometry/Dipole.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/operators/winding_number.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <iostream>

namespace primal = axom::primal;

// The GWN for a curve is the accumulation of dipole components summed across its length.
// We approximate those components by linearizing the curve, and summing dipoles generated at each.
// This is uniquely accurate in 2D, as the curve and the linearization have the same GWN at
//  all points outside a bounding box of each.
TEST(primal_dipole, winding_number_from_linearized_bezier_curve)
{
  using Point2D = primal::Point<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;
  using Segment = primal::Segment<double, 2>;
  using Dipole = primal::Dipole<double, 2>;

  constexpr int num_segments = 128;
  constexpr double edge_tol = 1.e-8;
  constexpr double abs_tol = 1.e-8;

  // This cubic follows a counterclockwise-oriented arc. Segment::normal()
  // gives its right normal, which is the outward normal for this orientation.
  const Point2D control_points[] = {Point2D {-1.5, 0.0},
                                    Point2D {-1.0, 1.0},
                                    Point2D {1.0, 1.0},
                                    Point2D {1.0, 0.5}};
  const Bezier curve(control_points, 3);

  axom::Array<Dipole> dipoles;
  dipoles.reserve(num_segments);
  for(int i = 0; i < num_segments; ++i)
  {
    const auto p0 = curve.evaluate(static_cast<double>(i) / num_segments);
    const auto p1 = curve.evaluate(static_cast<double>(i + 1) / num_segments);
    const Segment segment(p0, p1);
    dipoles.emplace_back(segment.at(0.5), segment.normal(), segment.length());
  }

  const Point2D query {4.0, -4.0};
  const double curve_gwn = primal::winding_number(query, curve, edge_tol);
  const double dipole_gwn = primal::winding_number(query, dipoles);

  double summed_dipole_gwn = 0.0;
  for(const auto& dipole : dipoles)
  {
    summed_dipole_gwn += primal::winding_number(query, dipole);
  }

  EXPECT_DOUBLE_EQ(dipole_gwn, summed_dipole_gwn);
  EXPECT_NEAR(dipole_gwn, curve_gwn, abs_tol);
}

// The GWN for a surface is the accumulation of dipole components summed across its area.
// We approximate those components by triangulating the surface, and summing dipoles
//  defined by the centroid and normal of each triangle.
TEST(primal_dipole, winding_number_from_triangulated_bezier_patch)
{
  using Point3D = primal::Point<double, 3>;
  using BezierPatch = primal::BezierPatch<double, 3>;
  using Triangle = primal::Triangle<double, 3>;
  using Dipole = primal::Dipole<double, 3>;

  constexpr int num_divisions = 32;
  constexpr double abs_tol = 1.e-6;

  // Each pair of triangles is ordered consistently
  //  with the patch's u-by-v orientation.
  const Point3D control_points[] = {Point3D {-1.0, -1.0, 1.0},
                                    Point3D {-1.0, 1.0, 0.0},
                                    Point3D {1.0, -1.0, 0.0},
                                    Point3D {1.0, 1.0, -1.0}};
  const BezierPatch patch(control_points, 1, 1);

  axom::Array<Dipole> dipoles;
  dipoles.reserve(2 * num_divisions * num_divisions);
  for(int i = 0; i < num_divisions; ++i)
  {
    for(int j = 0; j < num_divisions; ++j)
    {
      const double u0 = static_cast<double>(i) / num_divisions;
      const double u1 = static_cast<double>(i + 1) / num_divisions;
      const double v0 = static_cast<double>(j) / num_divisions;
      const double v1 = static_cast<double>(j + 1) / num_divisions;

      const auto p00 = patch.evaluate(u0, v0);
      const auto p10 = patch.evaluate(u1, v0);
      const auto p01 = patch.evaluate(u0, v1);
      const auto p11 = patch.evaluate(u1, v1);
      const Triangle first_triangle(p00, p10, p11);
      const Triangle second_triangle(p00, p11, p01);

      for(const Triangle& triangle : {first_triangle, second_triangle})
      {
        dipoles.emplace_back(triangle.centroid(), triangle.normal(), triangle.area());
      }
    }
  }

  const Point3D query {5.0, 5.0, -5.0};
  const double patch_gwn = primal::winding_number(query, patch);
  const double dipole_gwn = primal::winding_number(query, dipoles);

  double summed_dipole_gwn = 0.0;
  for(const auto& dipole : dipoles)
  {
    summed_dipole_gwn += primal::winding_number(query, dipole);
  }

  EXPECT_DOUBLE_EQ(dipole_gwn, summed_dipole_gwn);
  EXPECT_NEAR(dipole_gwn, patch_gwn, abs_tol);
}
