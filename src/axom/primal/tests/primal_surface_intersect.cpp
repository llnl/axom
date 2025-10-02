// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file primal_surface_intersect.cpp
 * \brief This file tests surface (bilinear/patch/NURBS) intersection routines
 */

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/operators/intersect.hpp"

#include "axom/core/numerics/matvecops.hpp"

#include <cmath>
#include <iomanip>
#include <algorithm>
#include <sstream>

namespace primal = axom::primal;

/*
 * Helper function to compute the intersections of a Bezier patch and a ray 
 * and check that their intersection points match our expectations.
 * Patch parameters are stored in \a exp_u, \a exp_v and \a exp_t.
 * Intersections are computed within tolerance \a eps and our checks use \a test_eps.
 *
 * Param \a shouldPrintIntersections is used for debugging and for generating
 * the initial array of expected intersections.
 */
template <typename CoordType>
void checkIntersections(const primal::Ray<CoordType, 3>& ray,
                        const primal::BezierPatch<CoordType, 3>& patch,
                        const axom::Array<CoordType>& exp_t,
                        const axom::Array<CoordType>& exp_u,
                        const axom::Array<CoordType>& exp_v,
                        double eps,
                        double test_eps,
                        bool isHalfOpen = false,
                        bool shouldPrintIntersections = false)
{
  constexpr int DIM = 3;
  using Array = axom::Array<CoordType>;

  // Check validity of input data.
  // They should have the same size
  EXPECT_EQ(exp_u.size(), exp_v.size());
  EXPECT_EQ(exp_u.size(), exp_t.size());

  const int num_exp_intersections = static_cast<int>(exp_u.size());
  const bool exp_intersect = (num_exp_intersections > 0);

  // Intersect the ray and the patch, intersection parameters will be
  // in arrays (u, v) and t, for the patch and ray, respectively
  Array u, v, t;
  bool ray_intersects = intersect(ray, patch, t, u, v, eps, eps, isHalfOpen);
  EXPECT_EQ(exp_intersect, ray_intersects);
  EXPECT_EQ(u.size(), v.size());
  EXPECT_EQ(u.size(), t.size());

  // check that we found the expected number of intersection points
  const int num_actual_intersections = static_cast<int>(u.size());
  EXPECT_EQ(num_exp_intersections, num_actual_intersections);

  // check that the evaluated intersection points are identical
  for(int i = 0; i < num_actual_intersections; ++i)
  {
    auto p1 = ray.at(t[i]);
    auto p2 = patch.evaluate(u[i], v[i]);

    EXPECT_NEAR(0., primal::squared_distance(p1, p2), test_eps);

    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_NEAR(p1[d], p2[d], test_eps);
    }
  }

  if(shouldPrintIntersections)
  {
    std::stringstream sstr;

    sstr << "Intersections for ray and patch: "
         << "\n\t" << ray << "\n\t" << patch;

    sstr << "\ns (" << u.size() << "): ";
    for(auto i = 0; i < u.size(); ++i)
    {
      sstr << std::setprecision(16) << "(" << u[i] << "," << v[i] << "),";
    }

    sstr << "\nt (" << t.size() << "): ";
    for(auto i = 0; i < t.size(); ++i)
    {
      sstr << std::setprecision(16) << t[i] << ",";
    }

    SLIC_INFO(sstr.str());
  }

  for(int i = 0; i < num_actual_intersections; ++i)
  {
    EXPECT_NEAR(exp_u[i], u[i], test_eps);
    EXPECT_NEAR(exp_v[i], v[i], test_eps);
    EXPECT_NEAR(exp_t[i], t[i], test_eps);

    if(shouldPrintIntersections)
    {
      SLIC_INFO("\t" << i << ": {u:" << u[i] << ", v:" << v[i] << std::setprecision(16)
                     << ", t:" << t[i] << ", u_actual:" << exp_u[i] << ", v_actual:" << exp_v[i]
                     << ", t_actual:" << exp_t[i] << "}");
    }
  }
}

/*
 * Helper function to compute the intersections of a NURBS patch and a ray 
 * and check that their intersection points match our expectations.
 * Patch parameters are stored in \a exp_u, \a exp_v and \a exp_t.
 * Intersections are computed within tolerance \a eps and our checks use \a test_eps.
 *
 * Param \a shouldPrintIntersections is used for debugging and for generating
 * the initial array of expected intersections.
 */
template <typename CoordType>
void checkIntersections(const primal::Ray<CoordType, 3>& ray,
                        const primal::NURBSPatch<CoordType, 3>& patch,
                        const axom::Array<CoordType>& exp_t,
                        const axom::Array<CoordType>& exp_u,
                        const axom::Array<CoordType>& exp_v,
                        double eps,
                        double test_eps,
                        bool countUntrimmed = true,
                        bool isHalfOpen = false,
                        bool shouldPrintIntersections = false)
{
  constexpr int DIM = 3;
  using Array = axom::Array<CoordType>;

  // Check validity of input data.
  // They should have the same size
  EXPECT_EQ(exp_u.size(), exp_v.size());
  EXPECT_EQ(exp_u.size(), exp_t.size());

  const int num_exp_intersections = static_cast<int>(exp_u.size());
  const bool exp_intersect = (num_exp_intersections > 0);

  // Intersect the ray and the patch, intersection parameters will be
  // in arrays (u, v) and t, for the patch and ray, respectively
  Array u, v, t;
  bool ray_intersects = intersect(ray, patch, t, u, v, eps, eps, countUntrimmed, isHalfOpen);
  EXPECT_EQ(exp_intersect, ray_intersects);
  EXPECT_EQ(u.size(), v.size());
  EXPECT_EQ(u.size(), t.size());

  // check that we found the expected number of intersection points
  const int num_actual_intersections = static_cast<int>(u.size());
  EXPECT_EQ(num_exp_intersections, num_actual_intersections);

  // check that the evaluated intersection points are identical
  for(int i = 0; i < num_actual_intersections; ++i)
  {
    auto p1 = ray.at(t[i]);
    auto p2 = patch.evaluate(u[i], v[i]);

    EXPECT_NEAR(0., primal::squared_distance(p1, p2), test_eps);

    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_NEAR(p1[d], p2[d], test_eps);
    }
  }

  if(shouldPrintIntersections)
  {
    std::stringstream sstr;

    sstr << "Intersections for ray and patch: "
         << "\n\t" << ray << "\n\t" << patch;

    sstr << "\ns (" << u.size() << "): ";
    for(auto i = 0; i < u.size(); ++i)
    {
      sstr << std::setprecision(16) << "(" << u[i] << "," << v[i] << "),";
    }

    sstr << "\nt (" << t.size() << "): ";
    for(auto i = 0; i < t.size(); ++i)
    {
      sstr << std::setprecision(16) << t[i] << ",";
    }

    SLIC_INFO(sstr.str());
  }

  for(int i = 0; i < num_actual_intersections; ++i)
  {
    EXPECT_NEAR(exp_u[i], u[i], test_eps);
    EXPECT_NEAR(exp_v[i], v[i], test_eps);
    EXPECT_NEAR(exp_t[i], t[i], test_eps);

    if(shouldPrintIntersections)
    {
      SLIC_INFO("\t" << i << ": {u:" << u[i] << ", v:" << v[i] << std::setprecision(16)
                     << ", t:" << t[i] << ", u_actual:" << exp_u[i] << ", v_actual:" << exp_v[i]
                     << ", t_actual:" << exp_t[i] << "}");
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_surface_inter, bilinear_intersect)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;
  using RayType = primal::Ray<CoordType, DIM>;

  const double eps = 1E-10;
  const double eps_test = 1E-8;

  SLIC_INFO("primal: testing bilinear patch intersection");

  // Set control points
  BezierPatchType bilinear_patch(1, 1);
  bilinear_patch(0, 0) = PointType({-1.0, 1.0, 1.0});
  bilinear_patch(1, 0) = PointType({-1.0, -1.0, 2.0});
  bilinear_patch(1, 1) = PointType({1.0, -1.0, 1.0});
  bilinear_patch(0, 1) = PointType({1.0, 1.0, 2.0});

  // Ray with single intersection
  PointType ray_origin({0.0, 0.0, 1.75});
  VectorType ray_direction({1.0, 1.0, 0.0});
  RayType ray(ray_origin, ray_direction);

  checkIntersections(ray, bilinear_patch, {1.0}, {0.146446609407}, {0.853553390593}, eps, eps_test);

  // Ray with no intersections
  ray_direction = VectorType({1.0, -1.0, 0.0});
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test);

  // Ray with two intersections
  ray_origin = PointType({-1.0, -1.0, 1.75});
  ray_direction = VectorType({1.0, 1.0, 0.0});
  ray = RayType(ray_origin, ray_direction);

  checkIntersections(ray,
                     bilinear_patch,
                     {0.414213562373, 2.41421356237},
                     {0.853553390593, 0.146446609407},
                     {0.146446609407, 0.853553390593},
                     eps,
                     eps_test);

  // Ray with infinitely many intersections
  ray_origin = PointType({-2.0, 0.0, 1.5});
  ray_direction = VectorType({1.0, 0.0, 0.0});
  ray = RayType(ray_origin, ray_direction);

  checkIntersections(ray, bilinear_patch, {2.0}, {0.5}, {0.5}, eps, eps_test);

  // Ray with no intersections on line with infinitely many intersections
  ray_origin = PointType({2.0, 0.0, 1.5});
  ray_direction = VectorType({1.0, 0.0, 0.0});
  ray = RayType(ray_origin, ray_direction);

  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test);

  // Rays with infinitely many intersections, but the origin is on the patch
  ray_origin = PointType({0.4, 0.0, 1.5});
  ray_direction = VectorType({1.0, 0.0, 0.0});
  ray = RayType(ray_origin, ray_direction);

  checkIntersections(ray, bilinear_patch, {0.3}, {0.5}, {0.85}, eps, eps_test);

  ray_direction = VectorType({-1.0, 0.0, 0.0});
  ray = RayType(ray_origin, ray_direction);

  checkIntersections(ray, bilinear_patch, {0.7}, {0.5}, {0.35}, eps, eps_test);
}

//------------------------------------------------------------------------------
TEST(primal_surface_inter, bilinear_boundary_treatment)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;
  using RayType = primal::Ray<CoordType, DIM>;

  const double eps = 1E-10;
  const double eps_test = 1E-8;
  const bool isHalfOpen = true;

  SLIC_INFO("primal: testing bilinear patch intersection");

  // Set control points
  BezierPatchType bilinear_patch(1, 1);
  bilinear_patch(0, 0) = PointType({-1.0, 1.0, 1.0});
  bilinear_patch(1, 0) = PointType({-1.0, -1.0, 2.0});
  bilinear_patch(1, 1) = PointType({1.0, -1.0, 1.0});
  bilinear_patch(0, 1) = PointType({1.0, 1.0, 2.0});

  // Don't count intersections on the u=1 or v=1 isocurves
  PointType ray_origin({0.0, 0.0, 3.0});
  VectorType ray_direction;
  RayType ray(ray_origin, ray_direction);

  // Check both options for isHalfOpen

  ray_direction = VectorType(ray_origin, bilinear_patch.evaluate(0.0, 1.0));
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test, isHalfOpen);
  checkIntersections(ray, bilinear_patch, {sqrt(3.0)}, {0.0}, {1.0}, eps, eps_test, !isHalfOpen);

  ray_direction = VectorType(ray_origin, bilinear_patch.evaluate(0.5, 1.0));
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test, isHalfOpen);
  checkIntersections(ray, bilinear_patch, {sqrt(13.0) / 2.0}, {0.5}, {1.0}, eps, eps_test, !isHalfOpen);

  ray_direction = VectorType(ray_origin, bilinear_patch.evaluate(1.0, 0.5));
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test, isHalfOpen);
  checkIntersections(ray, bilinear_patch, {sqrt(13.0) / 2.0}, {1.0}, {0.5}, eps, eps_test, !isHalfOpen);

  ray_direction = VectorType(ray_origin, bilinear_patch.evaluate(1.0, 0.0));
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test, isHalfOpen);
  checkIntersections(ray, bilinear_patch, {sqrt(3.0)}, {1.0}, {0.0}, eps, eps_test, !isHalfOpen);

  ray_direction = VectorType(ray_origin, bilinear_patch.evaluate(1.0, 0.5));
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test, isHalfOpen);
  checkIntersections(ray, bilinear_patch, {sqrt(13.0) / 2.0}, {1.0}, {0.5}, eps, eps_test, !isHalfOpen);

  ray_direction = VectorType(ray_origin, bilinear_patch.evaluate(1.0, 1.0));
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test, isHalfOpen);
  checkIntersections(ray, bilinear_patch, {sqrt(6.0)}, {1.0}, {1.0}, eps, eps_test, !isHalfOpen);

  // These should record an intersection with both options

  ray_direction = VectorType(ray_origin, bilinear_patch.evaluate(0.0, 0.0));
  ray = RayType(ray_origin, ray_direction);
  for(bool option : {isHalfOpen, !isHalfOpen})
  {
    checkIntersections(ray, bilinear_patch, {sqrt(6.0)}, {0.0}, {0.0}, eps, eps_test, option);
  }

  ray_direction = VectorType(ray_origin, bilinear_patch.evaluate(0.5, 0.0));
  ray = RayType(ray_origin, ray_direction);
  for(bool option : {isHalfOpen, !isHalfOpen})
  {
    checkIntersections(ray, bilinear_patch, {sqrt(13.0) / 2.0}, {0.5}, {0.0}, eps, eps_test, option);
  }

  ray_direction = VectorType(ray_origin, bilinear_patch.evaluate(0.0, 0.5));
  ray = RayType(ray_origin, ray_direction);
  for(bool option : {isHalfOpen, !isHalfOpen})
  {
    checkIntersections(ray, bilinear_patch, {sqrt(13.0) / 2.0}, {0.0}, {0.5}, eps, eps_test, option);
  }
}

//------------------------------------------------------------------------------
TEST(primal_surface_inter, difficult_garp_case)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;
  using RayType = primal::Ray<CoordType, DIM>;

  const double eps = 1E-10;
  const double eps_test = 1E-8;

  SLIC_INFO("primal: testing bilinear patch intersection");

  // Set control points
  BezierPatchType bilinear_patch(1, 1);
  bilinear_patch(0, 0) = PointType({-2.0, 1.0, 2.0});
  bilinear_patch(1, 0) = PointType({-1.0, -1.0, 1.0});
  bilinear_patch(1, 1) = PointType({1.0, -1.0, 1.0});
  bilinear_patch(0, 1) = PointType({2.0, 1.0, 1.0});

  VectorType ray_direction({-1.0, -2.0, 0.0});
  // The first step of the GARP algorithm is to solve a quadratic equation a + bt + ct^2 = 0,
  //   and this configuration of patch + ray direction is such that c = 0

  // Ray with single intersection
  PointType ray_origin({0.0, 1.0, 1.25});
  RayType ray(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {sqrt(20. / 9.)}, {2. / 3.}, {0.25}, eps, eps_test);

  // Ray with no intersections
  ray_origin = PointType({0.0, 1.0, 1.75});
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test);

  bilinear_patch(0, 0) = PointType({-1.0, 1.0, 1.0});
  bilinear_patch(1, 0) = PointType({-1.0, -1.0, 2.0});
  bilinear_patch(1, 1) = PointType({1.0, -1.0, 1.0});
  bilinear_patch(0, 1) = PointType({1.0, 1.0, 2.0});

  // Double roots in the quadratic
  ray_origin = PointType({2.0, 0.0, 1.5});
  ray_direction = VectorType({-1.0, 0.0, 0.0});
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {2.0}, {0.5}, {0.5}, eps, eps_test);

  ray_origin = PointType({0.0, 2.0, 1.5});
  ray_direction = VectorType({0.0, -1.0, 0.0});
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {2.0}, {0.5}, {0.5}, eps, eps_test);

  ray_origin = PointType({0.0, 0.0, 1.75});
  ray_direction = VectorType({1.0, 0.0, 0.0});
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test);

  // Give patch a degeneracy at the point of intersection
  //  at which there are infinitely many parameters of intersection
  bilinear_patch(1, 1) = PointType({-1.0, -1.0, 2.0});

  ray_origin = PointType({-2.0, 0.0, 2.0});
  ray_direction = VectorType({1.0, -1.0, 0.0});
  ray = RayType(ray_origin, ray_direction);

  // Current behavior is to return a single point of intersection,
  //  as in the above cases with infinitely many intersections
  checkIntersections(ray, bilinear_patch, {sqrt(2)}, {1.0}, {0.5}, eps, eps_test);

  // Set ray origin to the point of degeneracy
  ray_origin = PointType({-1.0, -1.0, 2.0});
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {0.0}, {1.0}, {0.5}, eps, eps_test);

  // Set ray origin past the point of degeneracy (no intersections)
  ray_origin = PointType({0.0, -2.0, 2.0});
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test);
}

//------------------------------------------------------------------------------
TEST(primal_surface_inter, flat_bilinear_intersect)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;
  using RayType = primal::Ray<CoordType, DIM>;

  const double eps = 1E-10;
  const double eps_test = 1E-8;

  SLIC_INFO("primal: testing bilinear patch intersection");

  // Set control points
  BezierPatchType bilinear_patch(1, 1);
  bilinear_patch(0, 0) = PointType({-2.0, 1.0, 1.0});
  bilinear_patch(1, 0) = PointType({-1.0, -1.0, 1.0});
  bilinear_patch(1, 1) = PointType({1.0, -1.0, 1.0});
  bilinear_patch(0, 1) = PointType({2.0, 1.0, 1.0});

  // Ray with single intersection
  PointType ray_origin({0.0, 0.0, 2.0});
  VectorType ray_direction({1, 0.5, -1.0});
  RayType ray(ray_origin, ray_direction);

  checkIntersections(ray, bilinear_patch, {1.5}, {0.25}, {11. / 14.}, eps, eps_test);

  // Ray with a single intersection that is coplanar with an isocurve
  ray_direction = VectorType({1.0, 0.0, -1.0});
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {sqrt(2)}, {0.5}, {5. / 6.}, eps, eps_test);

  // Ray with no intersections
  ray_direction = VectorType({1.0, -1.0, -0.5});
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test);

  // Ray with no intersections that is parallel to the patch
  ray_direction = VectorType({1.0, 1.0, 0.0});
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test);

  // Ray with no intersections that is parallel to the patch along an isocurve
  ray_direction = VectorType({1.0, 0.0, 0.0});
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test);

  // Ray with no intersections on a line with infinitely many intersections
  ray_origin = PointType({-2.0, 0.5, 1});
  ray_direction = VectorType({-1.0, 0.0, 0.0});
  ray = RayType(ray_origin, ray_direction);

  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test);

  // Ray with infinitely many intersections
  ray_origin = PointType({-2.0, 0.0, 1.0});
  ray_direction = VectorType({1.0, 0.0, 0.0});
  ray = RayType(ray_origin, ray_direction);

  // Ray that is coplanar with a flat patch
  //  For ease of implementation, always return false
  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test);
}

//------------------------------------------------------------------------------
TEST(primal_surface_inter, flat_selfintersect_bilinear_intersect)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;
  using RayType = primal::Ray<CoordType, DIM>;

  const double eps = 1E-10;
  const double eps_test = 1E-8;

  SLIC_INFO("primal: testing bilinear patch intersection");

  // Set control points for square, hourglass patch
  BezierPatchType bilinear_patch(1, 1);
  bilinear_patch(0, 0) = PointType({-1.0, 1.0, 1.0});
  bilinear_patch(1, 0) = PointType({1.0, -1.0, 1.0});
  bilinear_patch(1, 1) = PointType({-1.0, -1.0, 1.0});
  bilinear_patch(0, 1) = PointType({1.0, 1.0, 1.0});

  // Ray with single intersection at the overlap point
  PointType ray_origin({0.0, 0.0, 2.0});
  VectorType ray_direction({0.0, 0.0, -1.0});
  RayType ray(ray_origin, ray_direction);

  checkIntersections(ray, bilinear_patch, {1.0}, {0.5}, {0.5}, eps, eps_test);

  // Ray with single intersection not at the overlap point
  ray_origin = PointType({0.0, 0.5, 2.0});
  ray = RayType(ray_origin, ray_direction);

  checkIntersections(ray, bilinear_patch, {1.0}, {0.25}, {0.5}, eps, eps_test);

  // Ray with no intersections
  ray_origin = PointType({0.5, 0.0, 2.0});
  ray = RayType(ray_origin, ray_direction);

  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test);

  // Change the patch to have a true self-overlapping region
  bilinear_patch(1, 1) = PointType({-1.0, 0.0, 1.0});

  // Has two intersections in parameter space
  ray_origin = PointType({-0.1, 0.25, 2.0});
  ray = RayType(ray_origin, ray_direction);

  checkIntersections(ray, bilinear_patch, {1.0, 1.0}, {0.6, 5. / 12.}, {0.75, 0.2}, eps, eps_test);

  // Is "tangent" to the overlap, resulting in a single intersection
  ray_origin = PointType({0.0, 0.25, 2.0});
  ray = RayType(ray_origin, ray_direction);

  checkIntersections(ray, bilinear_patch, {1.0}, {0.5}, {0.5}, eps, eps_test);
}

//------------------------------------------------------------------------------
TEST(primal_surface_inter, bezier_surface_intersect)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;
  using RayType = primal::Ray<CoordType, DIM>;

  double rt2 = sqrt(2), rt3 = sqrt(3), rt6 = sqrt(6);

  // Define the nodes and weights for one of six rational, biquartic Bezier patches
  //  that compose the unit sphere. This patch is centered at the +z axis.
  // Nodes and weights obtained from the technical report
  // "Tiling the Sphere with Rational Bezier Patches",
  //  James E. Cobb, University of Utah, 1988

  // clang-format off
  // These are the *homogeneous* coordinates
  axom::Array<PointType> node_data = {
    PointType {4*(1-rt3),     4*(1-rt3),     4*(rt3-1)}, PointType {rt2*(rt3-4),            -rt2, rt2*(4-rt3)}, PointType {4*(1-2*rt3)/3,   0, 4*(2*rt3-1)/3}, PointType {rt2*(rt3-4),           rt2,   rt2*(4-rt3)}, PointType {4*(1-rt3),     4*(rt3-1),     4*(rt3-1)},
    PointType {     -rt2, rt2*(rt3 - 4), rt2*(4 - rt3)}, PointType {(2-3*rt3)/2,     (2-3*rt3)/2,   (rt3+6)/2}, PointType {rt2*(2*rt3-7)/3, 0,       5*rt6/3}, PointType {(2-3*rt3)/2,   (3*rt3-2)/2,     (rt3+6)/2}, PointType {     -rt2,   rt2*(4-rt3),   rt2*(4-rt3)},
    PointType {        0, 4*(1-2*rt3)/3, 4*(2*rt3-1)/3}, PointType {          0, rt2*(2*rt3-7)/3,     5*rt6/3}, PointType {0,               0,   4*(5-rt3)/3}, PointType {          0, rt2*(7-2*rt3)/3,     5*rt6/3}, PointType {        0, 4*(2*rt3-1)/3, 4*(2*rt3-1)/3},
    PointType {      rt2, rt2*(rt3 - 4), rt2*(4 - rt3)}, PointType {(3*rt3-2)/2,     (2-3*rt3)/2,   (rt3+6)/2}, PointType {rt2*(7-2*rt3)/3, 0,       5*rt6/3}, PointType {(3*rt3-2)/2,   (3*rt3-2)/2,     (rt3+6)/2}, PointType {      rt2,   rt2*(4-rt3),   rt2*(4-rt3)},
    PointType {4*(rt3-1),     4*(1-rt3),     4*(rt3-1)}, PointType {rt2*(4-rt3),            -rt2, rt2*(4-rt3)}, PointType {4*(2*rt3-1)/3,   0, 4*(2*rt3-1)/3}, PointType {rt2*(4-rt3),           rt2,   rt2*(4-rt3)}, PointType {4*(rt3-1),     4*(rt3-1),     4*(rt3-1)}};

  axom::Array<double> weight_data = {
         4*(3-rt3), rt2*(3*rt3-2),   4*(5-rt3)/3, rt2*(3*rt3-2),     4*(3-rt3),
     rt2*(3*rt3-2),     (rt3+6)/2, rt2*(rt3+6)/3,     (rt3+6)/2, rt2*(3*rt3-2),
       4*(5-rt3)/3, rt2*(rt3+6)/3, 4*(5*rt3-1)/9, rt2*(rt3+6)/3,   4*(5-rt3)/3,
     rt2*(3*rt3-2),     (rt3+6)/2, rt2*(rt3+6)/3,     (rt3+6)/2, rt2*(3*rt3-2),
         4*(3-rt3), rt2*(3*rt3-2),   4*(5-rt3)/3, rt2*(3*rt3-2),     4*(3-rt3)};
  // clang-format on

  BezierPatchType sphere_face_patch(node_data, weight_data, 4, 4);
  for(int i = 0; i < 5; ++i)
  {
    for(int j = 0; j < 5; ++j)
    {
      sphere_face_patch(i, j).array() /= sphere_face_patch.getWeight(i, j);
    }
  }

  // Intersections with arbitrary patches aren't recorded
  //  with as much precision as the base bilinear case
  const double eps = 1E-5;
  const double eps_test = 1E-5;
  const bool isHalfOpen = true;

  PointType ray_origin({0.0, 0.0, 0.0});

  // Try points which will be interior to all subdivisions
  double u_params[8], v_params[8];
  axom::numerics::linspace(0.0, 1.0, u_params, 8);
  axom::numerics::linspace(0.0, 1.0, v_params, 8);

  for(int i = 0; i < 8; ++i)
  {
    for(int j = 0; j < 8; ++j)
    {
      VectorType ray_direction(ray_origin, sphere_face_patch.evaluate(u_params[i], v_params[j]));
      RayType ray(ray_origin, ray_direction);

      // Points on the edge should not be recorded
      if(i == 9 || j == 9)
      {
        // Check both settings of isHalfOpen
        checkIntersections(ray, sphere_face_patch, {}, {}, {}, eps, eps_test, isHalfOpen);

        checkIntersections(ray,
                           sphere_face_patch,
                           {1.0},
                           {u_params[i]},
                           {v_params[j]},
                           eps,
                           eps_test,
                           !isHalfOpen);
      }
      else
      {
        checkIntersections(ray, sphere_face_patch, {1.0}, {u_params[i]}, {v_params[j]}, eps, eps_test);
      }
    }
  }

  // Use different parameter values so that intersections
  //  are on the boundary of subdivisions
  axom::numerics::linspace(0.0, 1.0, u_params, 5);
  axom::numerics::linspace(0.0, 1.0, v_params, 5);

  for(int i = 0; i < 5; ++i)
  {
    for(int j = 0; j < 5; ++j)
    {
      VectorType ray_direction(ray_origin, sphere_face_patch.evaluate(u_params[i], v_params[j]));
      RayType ray(ray_origin, ray_direction);

      // Points on the edge should not be recorded
      if(i == 8 || j == 8)
      {
        checkIntersections(ray, sphere_face_patch, {}, {}, {}, eps, eps_test, isHalfOpen);

        checkIntersections(ray,
                           sphere_face_patch,
                           {1.0},
                           {u_params[i]},
                           {v_params[j]},
                           eps,
                           eps_test,
                           !isHalfOpen);
      }
      else
      {
        checkIntersections(ray, sphere_face_patch, {1.0}, {u_params[i]}, {v_params[j]}, eps, eps_test);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_surface_inter, NURBS_surface_intersect)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;
  using RayType = primal::Ray<CoordType, DIM>;

  // Represent the sphere with a single NURBS patch
  axom::Array<CoordType> knotvec_u = {-2.0, -2.0, -2.0, -2.0, 1.0, 1.0, 1.0, 1.0};
  axom::Array<CoordType> knotvec_v = {0.0, 0.0, 0.0, 0.0, 1.5, 1.5, 1.5, 3.0, 3.0, 3.0, 3.0};

  // clang-format off
  axom::Array<PointType> node_data = {
    PointType {0, 0,  1}, PointType {0, 0,  1}, PointType { 0, 0,  1}, PointType { 0, 0,  1}, PointType { 0,  0,  1}, PointType {0,  0,  1}, PointType {0, 0,  1},
    PointType {2, 0,  1}, PointType {2, 4,  1}, PointType {-2, 4,  1}, PointType {-2, 0,  1}, PointType {-2, -4,  1}, PointType {2, -4,  1}, PointType {2, 0,  1},
    PointType {2, 0, -1}, PointType {2, 4, -1}, PointType {-2, 4, -1}, PointType {-2, 0, -1}, PointType {-2, -4, -1}, PointType {2, -4, -1}, PointType {2, 0, -1},
    PointType {0, 0, -1}, PointType {0, 0, -1}, PointType { 0, 0, -1}, PointType { 0, 0, -1}, PointType { 0,  0, -1}, PointType {0,  0, -1}, PointType {0, 0, -1}};

  axom::Array<CoordType> weight_data = {
        1.0, 1.0/3.0, 1.0/3.0,     1.0, 1.0/3.0, 1.0/3.0,     1.0,
    1.0/3.0, 1.0/9.0, 1.0/9.0, 1.0/3.0, 1.0/9.0, 1.0/9.0, 1.0/3.0,
    1.0/3.0, 1.0/9.0, 1.0/9.0, 1.0/3.0, 1.0/9.0, 1.0/9.0, 1.0/3.0,
        1.0, 1.0/3.0, 1.0/3.0,     1.0, 1.0/3.0, 1.0/3.0,     1.0};
  // clang-format on

  NURBSPatchType sphere_patch(node_data, weight_data, 4, 7, knotvec_u, knotvec_v);

  // Add extra knots to make the Bezier extraction more interesting
  sphere_patch.insertKnot_u(0.0, 1);
  sphere_patch.insertKnot_u(0.5, 2);
  sphere_patch.insertKnot_v(2.0, 3);

  // Test some parameter values that are at boundaries, knot values,
  //  subdivision boundaries, and interior points to subdivisions.
  // Also, skip the points on the u-edges of the patch,
  //  as they're degenerate in physical space.
  double params_u[6] = {-1.0, -0.9, 0.0, 0.3, 0.5, 0.75};
  double params_v[8] = {0.0, 0.3, 1.0, 1.6, 2.0, 2.5, 2.7, 3.0};

  // Intersections with arbitrary patches aren't recorded
  //  with as much precision as the base bilinear case
  const double eps = 1E-5;
  const double eps_test = 1E-5;
  const bool countUntrimmed = true;
  const bool isHalfOpen = true;

  PointType ray_origin({0.0, 0.0, 0.0});

  for(int i = 0; i < 6; ++i)
  {
    for(int j = 0; j < 8; ++j)
    {
      VectorType ray_direction(ray_origin, sphere_patch.evaluate(params_u[i], params_v[j]));
      RayType ray(ray_origin, ray_direction);

      // The sphere meets itself at the v-edges
      if(j == 0 || j == 7)
      {
        // Once if the surface is half-open
        checkIntersections(ray,
                           sphere_patch,
                           {1.0},
                           {params_u[i]},
                           {0.0},
                           eps,
                           eps_test,
                           !countUntrimmed,
                           isHalfOpen);

        // Twice if the surface is not half-open
        checkIntersections(ray,
                           sphere_patch,
                           {1.0, 1.0},
                           {params_u[i], params_u[i]},
                           {0.0, 3.0},
                           eps,
                           eps_test,
                           !countUntrimmed,
                           !isHalfOpen);
      }
      else
      {
        checkIntersections(ray, sphere_patch, {1.0}, {params_u[i]}, {params_v[j]}, eps, eps_test);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_surface_inter, trimmed_surface_intersect)
{
  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using RayType = primal::Ray<CoordType, DIM>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;

  using ParameterPointType = primal::Point<CoordType, 2>;
  using TrimmingCurveType = primal::NURBSCurve<CoordType, 2>;

  const int degree_u = 1;
  const int degree_v = 1;

  const int npts_u = 3;
  const int npts_v = 3;

  // Construct from 1D axom::Array
  axom::Array<PointType> controlPointsArray({PointType {0.0, 0.0, 1.0},
                                             PointType {0.0, 1.0, 0.0},
                                             PointType {0.0, 2.0, 0.0},
                                             PointType {1.0, 0.0, 0.0},
                                             PointType {1.0, 1.0, -1.0},
                                             PointType {1.0, 2.0, 0.0},
                                             PointType {2.0, 0.0, 0.0},
                                             PointType {2.0, 1.0, 0.0},
                                             PointType {2.0, 2.0, 1.0}});

  NURBSPatchType nPatch(controlPointsArray, npts_u, npts_v, degree_u, degree_v);

  // Add a simple trimming curve from a NURBS curve
  axom::Array<ParameterPointType> trimmingCurveControlPoints {ParameterPointType {0.25, 0.25},
                                                              ParameterPointType {0.75, 0.25},
                                                              ParameterPointType {0.75, 0.75},
                                                              ParameterPointType {0.25, 0.75},
                                                              ParameterPointType {0.25, 0.25}};

  TrimmingCurveType trimmingCurve(trimmingCurveControlPoints, 2);
  nPatch.addTrimmingCurve(trimmingCurve);

  // Check intersection with a ray that intersects the patch at
  //  the parameter point (0.5, 0.5)
  PointType ray_origin({1.0, 1.0, -3.0});
  VectorType ray_direction(PointType {1.0, 1.0, -3.0}, nPatch.evaluate(0.5, 0.5));
  RayType ray(ray_origin, ray_direction);

  const double eps = 1E-5;
  const double eps_test = 1E-5;
  const bool countUntrimmed = true;
  const bool isHalfOpen = true;

  // Intersects the trimmed surface
  checkIntersections(ray, nPatch, {2.0}, {0.5}, {0.5}, eps, eps_test, !countUntrimmed, isHalfOpen);

  // Also intersects the untrimmed surface
  checkIntersections(ray, nPatch, {2.0}, {0.5}, {0.5}, eps, eps_test, countUntrimmed, isHalfOpen);

  // Check intersections with a ray that doesn't intersect the trimmed patch
  ray_direction = VectorType(PointType {1.0, 1.0, -3.0}, nPatch.evaluate(0.8, 0.8));
  ray = RayType(ray_origin, ray_direction);

  // Doesn't intersect the trimmed surface
  checkIntersections(ray, nPatch, {}, {}, {}, eps, eps_test, !countUntrimmed, isHalfOpen);
  checkIntersections(ray, nPatch, {}, {}, {}, eps, eps_test, !countUntrimmed, isHalfOpen);

  // Does intersect the untrimmed surface
  checkIntersections(ray, nPatch, {3.3105890}, {0.8}, {0.8}, eps, eps_test, countUntrimmed, isHalfOpen);
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();

  return result;
}
