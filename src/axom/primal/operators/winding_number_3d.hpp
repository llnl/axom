// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file winding_number_3d.hpp
 *
 * \brief Consists of methods to compute the generalized winding number (GWN) 
 *        for points with respect to various 3D geometric objects.
 */

#ifndef AXOM_PRIMAL_WINDING_NUMBER_3D_HPP_
#define AXOM_PRIMAL_WINDING_NUMBER_3D_HPP_

// Axom includes
#include "axom/core.hpp"
#include "axom/config.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/Polyhedron.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/geometry/NURBSPatch.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"

#include "axom/primal/operators/detail/winding_number_3d_impl.hpp"
#include "axom/primal/operators/intersect.hpp"

// C++ includes
#include <cmath>

// MFEM includes
#ifdef AXOM_USE_MFEM
  #include "mfem.hpp"
#endif

namespace axom
{
namespace primal
{
/*!
 * \brief Computes the GWN for a 3D point wrt a 3D triangle
 *
 * \param [in] query The query point to test
 * \param [in] tri The 3D Triangle object
 * \param [in] isOnFace An optional return parameter if the point is on the triangle
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the GWN as the solid angle modulo 4pi using the formula from 
 *  Oosterom, Strackee, "The Solid Angle of a Plane Triangle" 
 *  IEEE Transactions on Biomedical Engineering, Vol BME-30, No. 2, February 1983
 * with extra adjustments if the triangle takes up a full octant
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 3>& q,
                      const Triangle<T, 3>& tri,
                      bool& isOnFace,
                      const double edge_tol = 1e-8,
                      const double EPS = 1e-8)
{
  using Vec3 = Vector<T, 3>;

  if(tri.area() == 0)
  {
    return 0;
  }

  const Vec3 a = tri[0] - q;
  const Vec3 b = tri[1] - q;
  const Vec3 c = tri[2] - q;

  // Compute norms. Possibly return early
  const double a_norm = a.norm();
  const double b_norm = b.norm();
  const double c_norm = c.norm();

  if(a_norm < edge_tol || b_norm < edge_tol || c_norm < edge_tol)
  {
    return 0;
  }

  const double num = Vec3::scalar_triple_product(a, b, c);
  if(axom::utilities::isNearlyEqual(num, 0.0, EPS))
  {
    isOnFace = true;
    return 0;
  }

  const double denom = a_norm * b_norm * c_norm + a_norm * b.dot(c) +
    b_norm * a.dot(c) + c_norm * a.dot(b);

  // Handle direct cases where argument to atan is undefined
  if(axom::utilities::isNearlyEqual(denom, 0.0, EPS))
  {
    return (num > 0) ? 0.25 : -0.25;
  }

  // Note: denom==0 and num==0 handled above
  if(denom > 0)
  {
    return 0.5 * M_1_PI * atan(num / denom);
  }
  else
  {
    return (num > 0) ? 0.5 * M_1_PI * atan(num / denom) + 0.5
                     : 0.5 * M_1_PI * atan(num / denom) - 0.5;
  }
}

/*!
 * \brief Computes the GWN for a 3D point wrt a 3D triangle
 *
 * \param [in] query The query point to test
 * \param [in] tri The 3D Triangle object
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the GWN for the triangle without an additional return parameter
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 3>& q,
                      const Triangle<T, 3>& tri,
                      const double edge_tol = 1e-8,
                      const double EPS = 1e-8)
{
  bool isOnFace = false;
  return winding_number(q, tri, isOnFace, edge_tol, EPS);
}

/*!
 * \brief Computes the GWN for a 3D point wrt a 3D planar polygon
 *
 * \param [in] query The query point to test
 * \param [in] poly The Polygon object
 * \param [in] isOnFace Return variable to show if the point is on the polygon
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * \pre Assumes the polygon is planar. Otherwise, a meaningless value is returned.
 * 
 * Triangulates the polygon and computes the triangular GWN for each component
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 3>& q,
                      const Polygon<T, 3>& poly,
                      bool& isOnFace,
                      const double edge_tol = 1e-8,
                      const double EPS = 1e-8)
{
  const int num_verts = poly.numVertices();
  if(num_verts < 3)
  {
    return 0;
  }

  double wn = 0.0;
  for(int i = 0; i < num_verts - 2; ++i)
  {
    wn += winding_number(q,
                         Triangle<T, 3>(poly[0], poly[i + 1], poly[i + 2]),
                         isOnFace,
                         edge_tol,
                         EPS);
  }

  return wn;
}

/*!
 * \brief Computes the GWN for a 3D point wrt a 3D planar polygon
 *
 * \param [in] query The query point to test
 * \param [in] poly The Polygon object
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * \pre Assumes the polygon is planar. Otherwise, a meaningless value is returned.
 * 
 * Computes the GWN for the polygon without an additional return parameter
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 3>& q,
                      const Polygon<T, 3>& poly,
                      const double edge_tol = 1e-8,
                      const double EPS = 1e-8)
{
  bool isOnFace = false;
  return winding_number(q, poly, isOnFace, edge_tol, EPS);
}

/*!
 * \brief Computes the winding number for a 3D point wrt a 3D convex polyhedron
 *
 * \param [in] query The query point to test
 * \param [in] poly The Polyhedron object
 * \param [in] includeBoundary If true, points on the boundary are considered interior.
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * \pre Expects the polyhedron to be convex and closed so that the returned value is an integer.
 * 
 * Computes the faces of the polyhedron and computes the GWN for each.
 * The sum is then rounded to the nearest integer, as the shape is assumed to be closed.
 * 
 * \return The integer winding number.
 */
template <typename T>
int winding_number(const Point<T, 3>& query,
                   const Polyhedron<T, 3>& poly,
                   bool includeBoundary = false,
                   double edge_tol = 1e-8,
                   double EPS = 1e-8)
{
  SLIC_ASSERT(poly.hasNeighbors());
  const int num_verts = poly.numVertices();

  axom::Array<int> faces(num_verts * num_verts), face_size(2 * num_verts),
    face_offset(2 * num_verts);
  int face_count;

  poly.getFaces(faces.data(), face_size.data(), face_offset.data(), face_count);

  bool isOnFace = false;
  double wn = 0;
  for(int i = 0; i < face_count; ++i)
  {
    const int N = face_size[i];
    const int i_offset = face_offset[i];
    Polygon<T, 3> the_face(N);
    for(int j = 0; j < N; ++j)
    {
      the_face.addVertex(poly[faces[i_offset + j]]);
    }

    wn += winding_number(query, the_face, isOnFace, edge_tol, EPS);

    if(isOnFace)
    {
      return includeBoundary;
    }
  }

  return std::lround(wn);
}

#ifdef AXOM_USE_MFEM

/*
 * \brief Computes the GWN for a 3D point wrt a 3D Bezier patch
 *
 * \param [in] query The query point to test
 * \param [in] bPatch The Bezier patch object
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] quad_tol The maximum relative error allowed in the quadrature
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * \param [in] depth The current recursive depth
 * 
 * Computes the generalized winding number for a Bezier patch using Stokes theorem.
 *
 * \note Warning: This algorithm is only tested to high accuracy for queries within
 *  1e-5 of the surface. Otherwise, it will return less accurate results.
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 3>& query,
                      const BezierPatch<T, 3>& bPatch,
                      const double edge_tol = 1e-8,
                      const double quad_tol = 1e-8,
                      const double EPS = 1e-8,
                      const int depth = 0)
{
  const int ord_u = bPatch.getOrder_u();
  const int ord_v = bPatch.getOrder_v();
  const bool patchIsRational = bPatch.isRational();
  const double edge_tol_sq = edge_tol * edge_tol;

  // Fix the number of quadrature nodes arbitrarily, but high enough
  //  to `catch` near singularities for refinement
  constexpr int quad_npts = 30;

  // Early return if the patch is approximately polygonal.
  //  Very slight variations in curvature requires small EPS tolerance
  constexpr int MAX_DEPTH = 10;
  if(depth >= MAX_DEPTH || bPatch.isPolygonal(EPS))
  {
    return winding_number(
      query,
      Polygon<T, 3>(axom::Array<Point<T, 3>>(
        {bPatch(0, 0), bPatch(ord_u, 0), bPatch(ord_u, ord_v), bPatch(0, ord_v)})),
      edge_tol,
      PRIMAL_TINY);
  }

  // Use a specific kind of recursion if we are within tol of an endpoint.
  //  Split the surface closer to the corner, assume smallest patch is polygonal,
  //  and set a new edge_tol so corners of the new patch aren't marked as coincident
  constexpr double edge_offset = 0.01;
  if(squared_distance(query, bPatch(0, 0)) <= edge_tol_sq)
  {
    BezierPatch<T, 3> p1, p2, p3, p4;
    bPatch.split(0.0 + edge_offset, 0.0 + edge_offset, p1, p2, p3, p4);
    double new_edge_tol = 0.5 *
      sqrt(axom::utilities::min(
        squared_distance(query, bPatch.evaluate(0.0, 0.0 + edge_offset)),
        squared_distance(query, bPatch.evaluate(0.0 + edge_offset, 0.0))));
    new_edge_tol = axom::utilities::min(new_edge_tol, edge_tol);

    return winding_number(query, p2, new_edge_tol, quad_tol, EPS, depth + 1) +
      winding_number(query, p3, new_edge_tol, quad_tol, EPS, depth + 1) +
      winding_number(query, p4, new_edge_tol, quad_tol, EPS, depth + 1);
  }
  if(squared_distance(query, bPatch(ord_u, 0)) <= edge_tol_sq)
  {
    BezierPatch<T, 3> p1, p2, p3, p4;
    bPatch.split(1.0 - edge_offset, 0.0 + edge_offset, p1, p2, p3, p4);
    double new_edge_tol = 0.5 *
      sqrt(axom::utilities::min(
        squared_distance(query, bPatch.evaluate(1.0, 0.0 + edge_offset)),
        squared_distance(query, bPatch.evaluate(1.0 - edge_offset, 0.0))));
    new_edge_tol = axom::utilities::min(new_edge_tol, edge_tol);

    return winding_number(query, p1, new_edge_tol, quad_tol, EPS, depth + 1) +
      winding_number(query, p3, new_edge_tol, quad_tol, EPS, depth + 1) +
      winding_number(query, p4, new_edge_tol, quad_tol, EPS, depth + 1);
  }
  if(squared_distance(query, bPatch(0, ord_v)) <= edge_tol_sq)
  {
    BezierPatch<T, 3> p1, p2, p3, p4;
    bPatch.split(0.0 + edge_offset, 1.0 - edge_offset, p1, p2, p3, p4);
    double new_edge_tol = 0.5 *
      sqrt(axom::utilities::min(
        squared_distance(query, bPatch.evaluate(0.0 + edge_offset, 1.0)),
        squared_distance(query, bPatch.evaluate(0.0, 1.0 - edge_offset))));
    new_edge_tol = axom::utilities::min(new_edge_tol, edge_tol);

    return winding_number(query, p1, new_edge_tol, quad_tol, EPS, depth + 1) +
      winding_number(query, p2, new_edge_tol, quad_tol, EPS, depth + 1) +
      winding_number(query, p4, new_edge_tol, quad_tol, EPS, depth + 1);
  }
  if(squared_distance(query, bPatch(ord_u, ord_v)) <= edge_tol_sq)
  {
    BezierPatch<T, 3> p1, p2, p3, p4;
    bPatch.split(1.0 - edge_offset, 1.0 - edge_offset, p1, p2, p3, p4);
    double new_edge_tol = 0.5 *
      sqrt(axom::utilities::min(
        squared_distance(query, bPatch.evaluate(1.0, 1.0 - edge_offset)),
        squared_distance(query, bPatch.evaluate(1.0 - edge_offset, 1.0))));
    new_edge_tol = axom::utilities::min(new_edge_tol, edge_tol);

    return winding_number(query, p1, new_edge_tol, quad_tol, EPS, depth + 1) +
      winding_number(query, p2, new_edge_tol, quad_tol, EPS, depth + 1) +
      winding_number(query, p3, new_edge_tol, quad_tol, EPS, depth + 1);
  }

  /* 
   * To use Stokes theorem, we need to identify a separating plane between
   * `query` and the surface, guaranteed through a bounding box.
   * If it does, need to do geometric refinement: Splitting and rotating the curve
   * until we can guarantee this.
   */
  CurvedPolygon<T, 3> boundingPoly(4);

  // Define vector fields whose curl gives us the winding number
  detail::DiscontinuityAxis field_direction;

  // Check an axis-aligned bounding box (most surfaces satisfy this condition)
  BoundingBox<T, 3> bBox(bPatch.boundingBox().expand(edge_tol));
  const bool exterior_x =
    bBox.getMin()[0] > query[0] || query[0] > bBox.getMax()[0];
  const bool exterior_y =
    bBox.getMin()[1] > query[1] || query[1] > bBox.getMax()[1];
  const bool exterior_z =
    bBox.getMin()[2] > query[2] || query[2] > bBox.getMax()[2];

  if(exterior_y || exterior_z)
  {
    field_direction = detail::DiscontinuityAxis::x;
  }
  else if(exterior_x || exterior_z)
  {
    field_direction = detail::DiscontinuityAxis::y;
  }
  else if(exterior_x || exterior_y)
  {
    field_direction = detail::DiscontinuityAxis::z;
  }
  else
  {
    // Next, check an oriented bounding box.
    // If we are interior to the oriented bounding box, then we
    //  cannot guarantee a separating plane, and need geometric refinement.
    OrientedBoundingBox<T, 3> oBox(bPatch.orientedBoundingBox().expand(edge_tol));
    if(oBox.contains(query))
    {
      BezierPatch<T, 3> p1, p2, p3, p4;
      bPatch.split(0.5, 0.5, p1, p2, p3, p4);
      return winding_number(query, p1, edge_tol, quad_tol, EPS, depth + 1) +
        winding_number(query, p2, edge_tol, quad_tol, EPS, depth + 1) +
        winding_number(query, p3, edge_tol, quad_tol, EPS, depth + 1) +
        winding_number(query, p4, edge_tol, quad_tol, EPS, depth + 1);
    }

    // Otherwise, we can apply a rotation to a z-aligned field.
    field_direction = detail::DiscontinuityAxis::rotated;

    // Lambda to generate a 3D rotation matrix from an angle and axis
    // Formulation from https://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
    auto angleAxisRotMatrix =
      [](double theta, const Vector<T, 3>& axis) -> numerics::Matrix<T> {
      const auto unitized = axis.unitVector();
      const double x = unitized[0], y = unitized[1], z = unitized[2];
      const double c = cos(theta), s = sin(theta), C = 1 - c;

      auto matx = numerics::Matrix<T>::zeros(3, 3);

      matx(0, 0) = x * x * C + c;
      matx(0, 1) = x * y * C - z * s;
      matx(0, 2) = x * z * C + y * s;

      matx(1, 0) = y * x * C + z * s;
      matx(1, 1) = y * y * C + c;
      matx(1, 2) = y * z * C - x * s;

      matx(2, 0) = z * x * C - y * s;
      matx(2, 1) = z * y * C + x * s;
      matx(2, 2) = z * z * C + c;

      return matx;
    };

    // Lambda to rotate the input point using the provided rotation matrix
    auto my_rotate_point = [&query](const numerics::Matrix<T>& matx,
                                    const Point<T, 3> input) -> Point<T, 3> {
      Vector<T, 3> shifted(query, input);
      Vector<T, 3> rotated;
      numerics::matrix_vector_multiply(matx, shifted.data(), rotated.data());
      return Point<T, 3>(
        {rotated[0] + query[0], rotated[1] + query[1], rotated[2] + query[2]});
    };

    // Find vector from query to the bounding box
    Point<T, 3> closest = closest_point(query, oBox);
    Vector<T, 3> v0 = Vector<T, 3>(query, closest).unitVector();

    // Find the direction of a ray perpendicular to that
    Vector<T, 3> v1;
    if(axom::utilities::isNearlyEqual(v0[0], v0[1], EPS))
    {
      v1 = Vector<T, 3>({v0[2], v0[2], -v0[0] - v0[1]}).unitVector();
    }
    else
    {
      v1 = Vector<T, 3>({-v0[1] - v0[2], v0[0], v0[0]}).unitVector();
    }

    // Rotate v0 around v1 until it is perpendicular to the plane spanned by k and v1
    double ang = (v0[2] < 0 ? 1.0 : -1.0) *
      acos(axom::utilities::clampVal(
        -(v0[0] * v1[1] - v0[1] * v1[0]) / sqrt(v1[0] * v1[0] + v1[1] * v1[1]),
        -1.0,
        1.0));
    auto rotator = angleAxisRotMatrix(ang, v1);

    // Collect rotated curves into the curved Polygon
    // Set up the (0, v) and (1, v) isocurves, rotated
    boundingPoly[0].setOrder(ord_v);
    boundingPoly[2].setOrder(ord_v);
    if(patchIsRational)
    {
      boundingPoly[0].makeRational();
      boundingPoly[2].makeRational();
    }
    for(int q = 0; q <= ord_v; ++q)
    {
      boundingPoly[0][q] = my_rotate_point(rotator, bPatch(ord_u, q));
      boundingPoly[2][q] = my_rotate_point(rotator, bPatch(0, ord_v - q));

      if(patchIsRational)
      {
        boundingPoly[0].setWeight(q, bPatch.getWeight(ord_u, q));
        boundingPoly[2].setWeight(q, bPatch.getWeight(0, ord_v - q));
      }
    }

    // Set up the (u, 0) and (u, 1) isocurves
    boundingPoly[1].setOrder(ord_u);
    boundingPoly[3].setOrder(ord_u);
    if(patchIsRational)
    {
      boundingPoly[1].makeRational();
      boundingPoly[3].makeRational();
    }
    for(int p = 0; p <= ord_u; ++p)
    {
      boundingPoly[1][p] = my_rotate_point(rotator, bPatch(ord_u - p, ord_v));
      boundingPoly[3][p] = my_rotate_point(rotator, bPatch(p, 0));

      if(patchIsRational)
      {
        boundingPoly[1].setWeight(p, bPatch.getWeight(ord_u - p, ord_v));
        boundingPoly[3].setWeight(p, bPatch.getWeight(p, 0));
      }
    }
  }

  // Set up the polygon if we don't need to do any rotation or splitting.
  if(field_direction != detail::DiscontinuityAxis::rotated)
  {
    //  Add the relevant bounding curves to the patch.
    boundingPoly[0] = bPatch.isocurve_u(0);
    boundingPoly[0].reverseOrientation();

    boundingPoly[1] = bPatch.isocurve_v(1);
    boundingPoly[1].reverseOrientation();

    boundingPoly[2] = bPatch.isocurve_u(1);
    boundingPoly[3] = bPatch.isocurve_v(0);
  }

  // Iterate over the edges of the bounding curved polygon, add up the results
  double wn = 0;
  for(int n = 0; n < 4; ++n)
  {
    wn += detail::stokes_winding_number(query,
                                        boundingPoly[n],
                                        field_direction,
                                        quad_npts,
                                        quad_tol);
  }

  return wn;
}

template <typename T>
double winding_number_casting(const Point<T, 3>& query,
                              const BezierPatch<T, 3>& bPatch,
                              const double edge_tol = 1e-8,
                              const double quad_tol = 1e-8,
                              const double EPS = 1e-8)
{
  auto wn_split =
    winding_number_casting_split(query, bPatch, edge_tol, quad_tol, EPS);
  return wn_split.first + wn_split.second;
}

/* This is the main workhorse of the new GWN algorithm!! */
template <typename T>
std::pair<double, double> winding_number_casting_split(
  const Point<T, 3>& query,
  const NURBSPatch<T, 3>& nPatch,
  const Vector<T, 3>& discontinuity_direction,
  const double edge_tol = 1e-8,
  const double quad_tol = 1e-8,
  const double EPS = 1e-8,
  const int depth = 0)
{
  const double edge_tol_sq = edge_tol * edge_tol;

  // Fix the number of quadrature nodes arbitrarily
  constexpr int quad_npts = 100;

  // The first is the GWN from stokes, the second is the jump condition
  std::pair<double, double> wn_split = {0.0, 0.0};

  /* 
   * To use Stokes theorem, we need to identify either a line containing the
   * query that does not intersect the surface, or one that intersects the *interior*
   * of the surface at known locations.
   */

  // Lambda to generate a 3D rotation matrix from an angle and axis
  // Formulation from https://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
  auto angleAxisRotMatrix = [](double theta,
                               const Vector<T, 3>& axis) -> numerics::Matrix<T> {
    const auto unitized = axis.unitVector();
    const double x = unitized[0], y = unitized[1], z = unitized[2];
    const double c = cos(theta), s = sin(theta), C = 1 - c;

    auto matx = numerics::Matrix<T>::zeros(3, 3);

    matx(0, 0) = x * x * C + c;
    matx(0, 1) = x * y * C - z * s;
    matx(0, 2) = x * z * C + y * s;

    matx(1, 0) = y * x * C + z * s;
    matx(1, 1) = y * y * C + c;
    matx(1, 2) = y * z * C - x * s;

    matx(2, 0) = z * x * C - y * s;
    matx(2, 1) = z * y * C + x * s;
    matx(2, 2) = z * z * C + c;

    return matx;
  };

  // Lambda to rotate the input point using the provided rotation matrix
  auto my_rotate_point = [&query](const numerics::Matrix<T>& matx,
                                  const Point<T, 3> input) -> Point<T, 3> {
    Vector<T, 3> shifted(query, input);
    Vector<T, 3> rotated;
    numerics::matrix_vector_multiply(matx, shifted.data(), rotated.data());
    return Point<T, 3>(
      {rotated[0] + query[0], rotated[1] + query[1], rotated[2] + query[2]});
  };

  // Lambda to generate a random orthogonal normal vector to the input
  auto random_orthogonal = [](const Vector<T, 3>& input) -> Vector<T, 3> {
    // Pick a random direction orthogonal to the surface normal
    //  at the center of the disk
    Vector<T, 3> some_perp(
      {input[1] - input[2], input[2] - input[0], input[0] - input[1]});
    double theta = axom::utilities::random_real(0.0, 2.0 * M_PI);
    Vector<T, 3> new_direction = std::cos(theta) * some_perp.unitVector() +
      std::sin(theta) * Vector3D::cross_product(input, some_perp).unitVector();

    return new_direction;
  };

  // Lambda to generate an entirely random unit vector
  auto random_unit = []() -> Vector<T, 3> {
    double theta = axom::utilities::random_real(0.0, 2 * M_PI);
    double u = axom::utilities::random_real(-1.0, 1.0);
    return Vector<T, 3> {sin(theta) * sqrt(1 - u * u),
                         cos(theta) * sqrt(1 - u * u),
                         u};
  };

  // Rotation matrix for the patch
  numerics::Matrix<T> rotator;

  // Prefer to work with a trimmed patch
  NURBSPatch<T, 3> rotatedPatch = nPatch;
  if(!rotatedPatch.isTrimmed()) rotatedPatch.makeSimpleTrimmed();

  BoundingBox<T, 3> bBox = nPatch.boundingBox();
  OrientedBoundingBox<T, 3> oBox = nPatch.orientedBoundingBox();
  auto characteristic_length = bBox.range().norm();

  bBox.expand(0.01 * characteristic_length);
  oBox.expand(0.01 * characteristic_length);

  // Define vector fields whose curl gives us the winding number
  detail::DiscontinuityAxis field_direction;

  // Case 1: Exterior without rotations
  if(!bBox.contains(query))
  {
    const bool exterior_x =
      bBox.getMin()[0] > query[0] || query[0] > bBox.getMax()[0];
    const bool exterior_y =
      bBox.getMin()[1] > query[1] || query[1] > bBox.getMax()[1];
    const bool exterior_z =
      bBox.getMin()[2] > query[2] || query[2] > bBox.getMax()[2];

    if(exterior_y || exterior_z)
    {
      field_direction = detail::DiscontinuityAxis::x;
    }
    else if(exterior_x || exterior_z)
    {
      field_direction = detail::DiscontinuityAxis::y;
    }
    else if(exterior_x || exterior_y)
    {
      field_direction = detail::DiscontinuityAxis::z;
    }
  }
  // Case 1.5: Exterior with rotation
  else if(!oBox.contains(query))
  {
    /* The following steps rotate the patch until the OBB is /not/ 
       directly above or below the query point */
    field_direction = detail::DiscontinuityAxis::rotated;

    // Find vector from query to the bounding box
    Point<T, 3> closest = closest_point(query, oBox);
    Vector<T, 3> v0 = Vector<T, 3>(query, closest).unitVector();

    // Find the direction of a ray perpendicular to that
    Vector<T, 3> v1;
    if(axom::utilities::isNearlyEqual(v0[0], v0[1], EPS))
    {
      v1 = Vector<T, 3>({v0[2], v0[2], -v0[0] - v0[1]}).unitVector();
    }
    else
    {
      v1 = Vector<T, 3>({-v0[1] - v0[2], v0[0], v0[0]}).unitVector();
    }

    // Rotate v0 around v1 until it is perpendicular to the plane spanned by k and v1
    double ang = (v0[2] < 0 ? 1.0 : -1.0) *
      acos(axom::utilities::clampVal(
        -(v0[0] * v1[1] - v0[1] * v1[0]) / sqrt(v1[0] * v1[0] + v1[1] * v1[1]),
        -1.0,
        1.0));
    rotator = angleAxisRotMatrix(ang, v1);
  }
  // Case 2: Cast a ray, record intersections
  else
  {
    field_direction = detail::DiscontinuityAxis::rotated;
    Line<T, 3> discontinuity_axis(query, discontinuity_direction);

    T patch_knot_size = axom::utilities::max(
      rotatedPatch.getKnots_u()[rotatedPatch.getNumKnots_u() - 1] -
        rotatedPatch.getKnots_u()[0],
      rotatedPatch.getKnots_v()[rotatedPatch.getNumKnots_v() - 1] -
        rotatedPatch.getKnots_v()[0]);

    // Tolerance for what counts as "close to a boundary" in parameter space
    T disk_radius = 0.1 * patch_knot_size;

    // Compute intersections with the *untrimmed and extrapolated* patch
    axom::Array<T> up, vp, tp;
    bool isHalfOpen = false, isTrimmed = false;
    bool success = intersect(discontinuity_axis,
                             rotatedPatch,
                             tp,
                             up,
                             vp,
                             1e-10,  // This is a good heuristic value for accuracy
                             EPS,
                             isHalfOpen);

    // for(int kk = 0; kk < up.size(); ++kk)
    // {
    // Point<T, 3> intersection_point = rotatedPatch.evaluate(up[kk], vp[kk]);
    // Point<T, 3> line_point = discontinuity_axis.at(tp[kk]);
    // std::cout << intersection_point << " " << line_point << std::endl;
    // std::cout << std::endl;
    // }

    if(!success)
    {
      // If too many intersections with the untrimmed patch are recorded, then either

      // 1. The ray is parallel to a the surface
      if(up.size() >= 20)
      {
        // Can't handle this case. Need to start over with a different cast direction
        auto new_direction = random_orthogonal(discontinuity_direction);

        // std::cout << "Recasting (very derogatory) at " << depth << std::endl;
        // std::cout << "Degenerate intersection point" << std::endl;
        return winding_number_casting_split(query,
                                            rotatedPatch,
                                            new_direction,
                                            edge_tol,
                                            quad_tol,
                                            EPS,
                                            depth + 1);
      }

      // 2. The surface is degenerate, and so all intersections were recorded
      //  with the same t parameter, and pruned by `intersect()`.

      // Treating this case requires some implementation I don't have yet,
      //  namely clipping trimming curves along u/v isocurves

      // For now, if *any* of the reported intersection points are coincident with the query,
      //  just return a zero
      for(int i = 0; i < up.size(); ++i)
      {
        Point<T, 3> intersection_point = rotatedPatch.evaluate(up[i], vp[i]);
        if(squared_distance(query, intersection_point) <= edge_tol_sq)
        {
          wn_split.first = 0.0;
          wn_split.second = 0.0;
          return wn_split;
        }
      }

      // Otherwise, we can cast again with a different direction
      {
        auto new_direction = random_orthogonal(discontinuity_direction);

        // std::cout << "Degenerate intersection point" << std::endl;
        // std::cout << "Recasting (derogatory) at " << depth << std::endl;
        return winding_number_casting_split(query,
                                            rotatedPatch,
                                            new_direction,
                                            edge_tol,
                                            quad_tol,
                                            EPS);
      }
    }

    // Account for each discontinuity in the integrand on the *untrimmed and extrapolated* surface

    // If no intersection, then nothing to account for

    // Otherwise, account for each discontinuity analytically or through disk subdivision
    for(int i = 0; i < up.size(); ++i)
    {
      // Check for surface degeneracies or tangencies
      Vector<T, 3> the_normal = rotatedPatch.normal(up[i], vp[i]);
      bool bad_intersection =
        axom::utilities::isNearlyEqual(the_normal.norm(), 0.0, EPS) ||
        axom::utilities::isNearlyEqual(
          the_normal.unitVector().dot(discontinuity_direction),
          0.0,
          EPS);

      // Check for surface coincidence
      Point<T, 3> intersection_point = rotatedPatch.evaluate(up[i], vp[i]);
      bool isOnSurface =
        squared_distance(query, intersection_point) <= edge_tol_sq;

      if(bad_intersection && !isOnSurface)
      {
        // If a far away ray intersects the surface at a tangent/cusp,
        //  can recast and try again
        auto new_direction = random_unit();

        // std::cout << "Recasting: " << the_normal.norm() << ", "
        // << the_normal.unitVector().dot(discontinuity_direction)
        // << std::endl;
        // std::cout << "Recasting at " << depth << std::endl;
        return winding_number_casting_split(query,
                                            rotatedPatch,
                                            new_direction,
                                            edge_tol,
                                            quad_tol,
                                            EPS,
                                            depth + 1);
      }

      if(isOnSurface)
      {
        // If the query point is on the surface, we need to consider a smaller disk
        //  to ensure its winding number is known to be near-zero
        disk_radius *= 0.1;
      }

      // This method accomplishes 2 tasks:
      //  > Determines if the disk of "safe" radius is *entirely* inside or outside trimming curves
      //  > If the disk intersects the trimming curves, performs disk subdivision
      bool isDiskInside, isDiskOutside, ignoreInteriorDisk = true;
      NURBSPatch<T, 3> disk_patch;

      {
        // rotatedPatch.printTrimmingCurves(
        //   "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\graphical_"
        //   "abstract\\trimming_curves\\original.txt");

        // AXOM_ANNOTATE_SCOPE("DISK_SPLIT");
        rotatedPatch.diskSplit(up[i],
                               vp[i],
                               disk_radius,
                               rotatedPatch,
                               disk_patch,
                               isDiskInside,
                               isDiskOutside,
                               ignoreInteriorDisk);

        // disk_patch.printTrimmingCurves(
        //   "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\graphical_"
        //   "abstract\\trimming_curves\\disk.txt");
        // rotatedPatch.printTrimmingCurves(
        //   "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\graphical_"
        //   "abstract\\trimming_curves\\remaining.txt");

        // axom::Array<T> up_1, vp_1, tp_1;
        // bool isHalfOpen = false, isTrimmed = false;
        // intersect(discontinuity_axis,
        //           disk_patch,
        //           tp_1,
        //           up_1,
        //           vp_1,
        //           1e-8,  // This is a good heuristic value for accuracy
        //           EPS,
        //           isHalfOpen);

        // for(int kk = 0; kk < up_1.size(); ++kk)
        // {
        //   Point<T, 3> intersection_point = disk_patch.evaluate(up_1[kk], vp_1[kk]);
        //   Point<T, 3> line_point = discontinuity_axis.at(tp_1[kk]);
        //   std::cout << intersection_point << std::endl;
        //   std::cout << line_point << std::endl;
        // }

        // std::cout << disk_patch.getKnots_u() << std::endl;
        // std::cout << disk_patch.getKnots_v() << std::endl;
        // disk_patch.clip(up[i] - 2 * disk_radius,
        //                 up[i] + 2 * disk_radius,
        //                 vp[i] - 2 * disk_radius,
        //                 vp[i] + 2 * disk_radius);
        // std::cout << disk_patch.getKnots_u() << std::endl;
        // std::cout << disk_patch.getKnots_v() << std::endl;

        // axom::Array<T> up_2, vp_2, tp_2;
        // intersect(discontinuity_axis,
        //           disk_patch,
        //           tp_2,
        //           up_2,
        //           vp_2,
        //           1e-8,  // This is a good heuristic value for accuracy
        //           EPS,
        //           isHalfOpen);

        // for(int kk = 0; kk < up_2.size(); ++kk)
        // {
        //   Point<T, 3> intersection_point = disk_patch.evaluate(up_2[kk], vp_2[kk]);
        //   Point<T, 3> line_point = discontinuity_axis.at(tp_2[kk]);
        //   std::cout << intersection_point << std::endl;
        //   std::cout << line_point << std::endl;
        // }

        // std::cout << "===========" << std::endl;
        // std::cout << up_1 << std::endl;
        // std::cout << up_2 << std::endl;
        // std::cout << "===========" << std::endl;
        // std::cout << vp_1 << std::endl;
        // std::cout << vp_2 << std::endl;
        // std::cout << "===========" << std::endl;
        // std::cout << discontinuity_axis.at(tp_e[0]) << std::endl;
        // std::cout << disk_patch.evaluate(up_e[0], vp_e[0]) << std::endl;
        // std::cout << rotatedPatch.evaluate(up_e[0], vp_e[0]) << std::endl;
        // std::cout << rotatedPatch.evaluate(up[i], vp[i]) << std::endl;
        // std::cout << disk_patch.evaluate(up[i], vp[i]) << std::endl;
        // disk_patch.clip(up[i] - 2 * disk_radius,
        //                 up[i] + 2 * disk_radius,
        //                 vp[i] - 2 * disk_radius,
        //                 vp[i] + 2 * disk_radius);
        // std::cout << disk_patch.evaluate(up[i], vp[i]) << std::endl;
        // std::cout << rotatedPatch.evaluate(up[i], vp[i]) << std::endl;
        // std::cout << disk_patch.evaluate(up[i], vp[i]) << std::endl;
        // std::cout << up[i] << ": " << up_e << std::endl;
        // std::cout << vp[i] << ": " << vp_e << std::endl;
        // std::cout << std::endl;
      }

      // If the query point is on the surface, the contribution of the disk is near-zero,
      //  and we only need to puncture the larger surface to proceed
      if(isOnSurface)
      {
        wn_split.first += 0.0;
        wn_split.second += 0.0;
      }
      // If the disk overlaps with a trimming curve
      else if(!isDiskInside && !isDiskOutside)
      {
        auto new_direction = random_orthogonal(the_normal);
        // std::cout << "SUBCAST " << new_direction[0] << ", " << new_direction[1]
        //           << ", " << new_direction[2] << std::endl;
        // We compute the contribution of the disk directly,
        //  but with a different direction to avoid repeated subdivision
        // std::cout << "Disk Subdivision" << std::endl;
        // if(depth < 3)
        {
          // std::cout << "Disk dividing at " << depth << std::endl;
          wn_split.first += winding_number_casting(query,
                                                   disk_patch,
                                                   new_direction,
                                                   edge_tol,
                                                   quad_tol,
                                                   EPS,
                                                   depth + 1);
        }
      }
      // If the disk is entirely inside or outside, the jump condition is known
      else if(isDiskOutside)
      {
        wn_split.second += 0.0;
      }
      else if(isDiskInside)
      {
        Vector<T, 3> the_direction =
          Vector<T, 3>(query, intersection_point).unitVector();

        // Do a dot product to see what side of the surface the intersection is on
        wn_split.second += std::copysign(0.5, the_normal.dot(the_direction));
      }
    }

    // Rotate the patch so that the discontinuity direction is aligned with the z-axis
    Vector<T, 3> axis = {discontinuity_direction[1],
                         -discontinuity_direction[0],
                         0.0};

    double ang =
      acos(axom::utilities::clampVal(discontinuity_direction[2], -1.0, 1.0));

    rotator = angleAxisRotMatrix(ang, axis);
  }

  if(field_direction == detail::DiscontinuityAxis::rotated)
  {
    // The trimming curves for rotatedPatch have been changed as needed,
    //  but we need to rotate the control points
    auto patch_shape = rotatedPatch.getControlPoints().shape();
    for(int i = 0; i < patch_shape[0]; ++i)
    {
      for(int j = 0; j < patch_shape[1]; ++j)
      {
        rotatedPatch(i, j) = my_rotate_point(rotator, nPatch(i, j));
      }
    }
  }

  wn_split.first += detail::stokes_winding_number(query,
                                                  rotatedPatch,
                                                  field_direction,
                                                  quad_npts,
                                                  quad_tol);

  return wn_split;
}

template <typename T>
double winding_number_casting(const Point<T, 3>& query,
                              const NURBSPatch<T, 3>& nPatch,
                              const double edge_tol = 1e-8,
                              const double quad_tol = 1e-8,
                              const double EPS = 1e-8,
                              const int depth = 0)
{
  double theta = axom::utilities::random_real(0.0, 2 * M_PI);
  double u = axom::utilities::random_real(-1.0, 1.0);
  auto cast_direction =
    Vector<T, 3> {sin(theta) * sqrt(1 - u * u), cos(theta) * sqrt(1 - u * u), u};

  auto wn_split = winding_number_casting_split(query,
                                               nPatch,
                                               cast_direction,
                                               edge_tol,
                                               quad_tol,
                                               EPS,
                                               depth);
  return wn_split.first + wn_split.second;
}

template <typename T>
double winding_number_casting(const Point<T, 3>& query,
                              const NURBSPatch<T, 3>& nPatch,
                              const Vector<T, 3>& discontinuity_direction,
                              const double edge_tol = 1e-8,
                              const double quad_tol = 1e-8,
                              const double EPS = 1e-8,
                              const int depth = 0)
{
  auto wn_split = winding_number_casting_split(query,
                                               nPatch,
                                               discontinuity_direction,
                                               edge_tol,
                                               quad_tol,
                                               EPS,
                                               depth);
  return wn_split.first + wn_split.second;
}

template <typename T>
std::pair<double, double> winding_number_casting_split(
  const Point<T, 3>& query,
  const BezierPatch<T, 3>& bPatch,
  const double edge_tol = 1e-8,
  const double quad_tol = 1e-8,
  const double EPS = 1e-8)
{
  const int ord_u = bPatch.getOrder_u();
  const int ord_v = bPatch.getOrder_v();
  const bool patchIsRational = bPatch.isRational();
  const double edge_tol_sq = edge_tol * edge_tol;

  // Fix the number of quadrature nodes arbitrarily, but high enough
  //  to `catch` near singularities for refinement
  constexpr int quad_npts = 50;

  // The first is the GWN from stokes, the second is the jump condition
  std::pair<double, double> wn_split = {0.0, 0.0};

  /* 
   * To use Stokes theorem, we need to identify either a line containing the
   * query that does not intersect the surface, or one that intersects the *interior*
   * of the surface a known number of times.
   */
  CurvedPolygon<T, 3> boundingPoly(4);

  // Lambda to generate a 3D rotation matrix from an angle and axis
  // Formulation from https://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
  auto angleAxisRotMatrix = [](double theta,
                               const Vector<T, 3>& axis) -> numerics::Matrix<T> {
    const auto unitized = axis.unitVector();
    const double x = unitized[0], y = unitized[1], z = unitized[2];
    const double c = cos(theta), s = sin(theta), C = 1 - c;

    auto matx = numerics::Matrix<T>::zeros(3, 3);

    matx(0, 0) = x * x * C + c;
    matx(0, 1) = x * y * C - z * s;
    matx(0, 2) = x * z * C + y * s;

    matx(1, 0) = y * x * C + z * s;
    matx(1, 1) = y * y * C + c;
    matx(1, 2) = y * z * C - x * s;

    matx(2, 0) = z * x * C - y * s;
    matx(2, 1) = z * y * C + x * s;
    matx(2, 2) = z * z * C + c;

    return matx;
  };

  // Lambda to rotate the input point using the provided rotation matrix
  auto my_rotate_point = [&query](const numerics::Matrix<T>& matx,
                                  const Point<T, 3> input) -> Point<T, 3> {
    Vector<T, 3> shifted(query, input);
    Vector<T, 3> rotated;
    numerics::matrix_vector_multiply(matx, shifted.data(), rotated.data());
    return Point<T, 3>(
      {rotated[0] + query[0], rotated[1] + query[1], rotated[2] + query[2]});
  };

  // Rotation matrix for the patch
  numerics::Matrix<T> rotator;

  OrientedBoundingBox<T, 3> oBox(bPatch.orientedBoundingBox().expand(edge_tol));
  if(!oBox.contains(query))  // Only do casting while debugging
  {
    /* The following steps rotate the patch until the OBB is /not/ 
       directly above or below the query point */

    // Find vector from query to the bounding box
    Point<T, 3> closest = closest_point(query, oBox);
    Vector<T, 3> v0 = Vector<T, 3>(query, closest).unitVector();

    // Find the direction of a ray perpendicular to that
    Vector<T, 3> v1;
    if(axom::utilities::isNearlyEqual(v0[0], v0[1], EPS))
    {
      v1 = Vector<T, 3>({v0[2], v0[2], -v0[0] - v0[1]}).unitVector();
    }
    else
    {
      v1 = Vector<T, 3>({-v0[1] - v0[2], v0[0], v0[0]}).unitVector();
    }

    // Rotate v0 around v1 until it is perpendicular to the plane spanned by k and v1
    double ang = (v0[2] < 0 ? 1.0 : -1.0) *
      acos(axom::utilities::clampVal(
        -(v0[0] * v1[1] - v0[1] * v1[0]) / sqrt(v1[0] * v1[0] + v1[1] * v1[1]),
        -1.0,
        1.0));
    rotator = angleAxisRotMatrix(ang, v1);
  }
  else
  {
    /* The following steps cast a ray and count /interior/ intersections with 
       the patch. If the intersection is tangent or on a boundary, try again */

    // Initial cast with a z-aligned field (i.e., no rotation)
    Vector<T, 3> singularity_direction;
    bool queryOnSurface = false;

    while(true)
    {
      wn_split.second = 0.0;
      bool retry_cast = false;

      // Pick a random cast direction
      double theta = axom::utilities::random_real(0.0, 2 * M_PI);
      double u = axom::utilities::random_real(-1.0, 1.0);
      singularity_direction = Vector<T, 3> {sin(theta) * sqrt(1 - u * u),
                                            cos(theta) * sqrt(1 - u * u),
                                            u};

      constexpr double buffer = 1e-3;

      // Compute intersections with the patch
      axom::Array<T> up, vp, tp;
      Line<T, 3> singularity_axis(query, singularity_direction);
      intersect(singularity_axis, bPatch, tp, up, vp, buffer);

      // If the line doesn't intersect the surface, then we can dodge it
      if(up.size() == 0)
      {
        wn_split.second += 0.0;
      }
      // Otherwise, we need to adjust for the contribution of a small removed
      // disk, depending on the orientation
      else
      {
        for(int i = 0; i < up.size(); ++i)
        {
          bool intersectionOnBoundary = false;
          bool badIntersection = false;

          Point<T, 3> intersect_point = bPatch.evaluate(up[i], vp[i]);

          Vector<T, 3> the_direction =
            Vector<T, 3>(query, intersect_point).unitVector();
          Vector<T, 3> the_normal = bPatch.normal(up[i], vp[i]).unitVector();

          // Do a dot product between the normal and the cast direction
          //  to see what side of the surface the intersection is on
          double surf_orientation = the_normal.dot(the_direction);
          double cast_orientation = singularity_direction.dot(the_direction);

          queryOnSurface =
            squared_distance(query, intersect_point) <= edge_tol_sq;
          badIntersection =
            axom::utilities::isNearlyEqual(surf_orientation, 0.0, 1e-3);
          intersectionOnBoundary = up[i] < buffer || up[i] > 1 - buffer ||
            vp[i] < buffer || vp[i] > 1 - buffer;

          if(queryOnSurface && !intersectionOnBoundary)
          {
            wn_split.second += 0;  // Do nothing for the jump condition
          }
          else if(!badIntersection && !queryOnSurface && !intersectionOnBoundary)
          {
            // Account for the jump condition analytically
            wn_split.second += std::copysign(0.5, surf_orientation);
          }
          else  // Brute force the rest of the split
          {
            BezierPatch<T, 3> coincident_patch(bPatch);
            BezierPatch<T, 3> subpatches[4];
            int n_subpatches = 0;

            // Split the patch around the intersection point
            // clang-format off
            if(up[i] - 0.1 > 0.0)
              coincident_patch.split_u(up[i] - 0.1, subpatches[n_subpatches++], coincident_patch);
            if(vp[i] - 0.1 > 0.0)
              coincident_patch.split_v(vp[i] - 0.1, subpatches[n_subpatches++], coincident_patch);
            if(up[i] + 0.1 < 1.0)
              coincident_patch.split_u(0.2 / (1.0 - (up[i] - 0.1) ), coincident_patch, subpatches[n_subpatches++]);
            if(vp[i] + 0.1 < 1.0)
              coincident_patch.split_v(0.2 / (1.0 - (vp[i] - 0.1) ), coincident_patch, subpatches[n_subpatches++]);
            // clang-format on

            for(int j = 0; j < n_subpatches; ++j)
            {
              auto wn_subpatch = winding_number_casting(query,
                                                        subpatches[j],
                                                        edge_tol,
                                                        quad_tol,
                                                        EPS);
              wn_split.first += wn_subpatch;
            }

            wn_split.first +=
              detail::surface_winding_number(query, coincident_patch, 10, quad_tol);
            return wn_split;
          }
        }
      }

      if(retry_cast)
      {
        continue;
      }

      // Define a rotation such that the cast ray is vertical
      Vector<T, 3> axis = {singularity_direction[1],
                           -singularity_direction[0],
                           0.0};
      double ang =
        acos(axom::utilities::clampVal(singularity_direction[2], -1.0, 1.0));

      rotator = angleAxisRotMatrix(ang, axis);

      break;
    }
  }

  /* With the rotation defined, 
      collect rotated curves into the curved Polygon */

  // Set up the (0, v) and (1, v) isocurves, rotated
  boundingPoly[0].setOrder(ord_v);
  boundingPoly[2].setOrder(ord_v);
  if(patchIsRational)
  {
    boundingPoly[0].makeRational();
    boundingPoly[2].makeRational();
  }
  for(int q = 0; q <= ord_v; ++q)
  {
    boundingPoly[0][q] = my_rotate_point(rotator, bPatch(ord_u, q));
    boundingPoly[2][q] = my_rotate_point(rotator, bPatch(0, ord_v - q));

    if(patchIsRational)
    {
      boundingPoly[0].setWeight(q, bPatch.getWeight(ord_u, q));
      boundingPoly[2].setWeight(q, bPatch.getWeight(0, ord_v - q));
    }
  }

  // Set up the (u, 0) and (u, 1) isocurves
  boundingPoly[1].setOrder(ord_u);
  boundingPoly[3].setOrder(ord_u);
  if(patchIsRational)
  {
    boundingPoly[1].makeRational();
    boundingPoly[3].makeRational();
  }
  for(int p = 0; p <= ord_u; ++p)
  {
    boundingPoly[1][p] = my_rotate_point(rotator, bPatch(ord_u - p, ord_v));
    boundingPoly[3][p] = my_rotate_point(rotator, bPatch(p, 0));

    if(patchIsRational)
    {
      boundingPoly[1].setWeight(p, bPatch.getWeight(ord_u - p, ord_v));
      boundingPoly[3].setWeight(p, bPatch.getWeight(p, 0));
    }
  }

  // Iterate over the edges of the bounding curved polygon, add up the results
  // std::cout << std::setprecision(15);
  // std::cout << query << std::endl;

  for(int n = 0; n < 4; ++n)
  {
    // std::cout << boundingPoly[n] << std::endl;
    auto this_val =
      detail::stokes_winding_number(query,
                                    boundingPoly[n],
                                    detail::DiscontinuityAxis::rotated,
                                    quad_npts,
                                    quad_tol);
    wn_split.first += this_val;
  }

  return wn_split;
}

template <typename T>
double winding_number_direct(const Point<T, 3>& query,
                             const BezierPatch<T, 3>& bPatch,
                             const double edge_tol = 1e-8,
                             const double quad_tol = 1e-8,
                             const double EPS = 1e-8)
{
  // Compute the winding number with a direct 2D surface integral
//   return detail::surface_winding_number(query, bPatch, 100, quad_tol, false);
  return detail::surface_winding_number(query, bPatch, 100);
}

template <typename T>
double winding_number_casting(const Point<T, 3>& query,
                              const NURBSPatchData<T>& nPatchData,
                              int& case_code,
                              int& integrated_trimming_curves,
                              const double edge_tol = 1e-8,
                              const double quad_tol = 1e-8,
                              const double EPS = 1e-8,
                              const int depth = 0)
{
  double theta = axom::utilities::random_real(0.0, 2 * M_PI);
  double u = axom::utilities::random_real(-1.0, 1.0);
  auto cast_direction =
    Vector<T, 3> {sin(theta) * sqrt(1 - u * u), cos(theta) * sqrt(1 - u * u), u};
  cast_direction = Vector<T, 3>({1, 0.5, 1}).unitVector();
  auto wn_split = winding_number_casting_split(query,
                                               nPatchData,
                                               cast_direction,
                                               case_code,
                                               integrated_trimming_curves,
                                               edge_tol,
                                               quad_tol,
                                               EPS,
                                               depth);
  return wn_split.first + wn_split.second;
}

template <typename T>
std::pair<double, double> winding_number_casting_split(
  const Point<T, 3>& query,
  const NURBSPatchData<T>& nPatchData,
  const Vector<T, 3>& discontinuity_direction,
  int& case_code,
  int& integrated_trimming_curves,
  const double edge_tol = 1e-8,
  const double quad_tol = 1e-8,
  const double EPS = 1e-8,
  const int depth = 0)
{
  const double edge_tol_sq = edge_tol * edge_tol;

  // Fix the number of quadrature nodes arbitrarily
  constexpr int quad_npts = 15;

  // The first is the GWN from stokes, the second is the jump condition
  std::pair<double, double> wn_split = {0.0, 0.0};

  /* 
   * To use Stokes theorem, we need to identify either a line containing the
   * query that does not intersect the surface, or one that intersects the *interior*
   * of the surface at known locations.
   */

  // Lambda to generate a 3D rotation matrix from an angle and axis
  // Formulation from https://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
  auto angleAxisRotMatrix = [](double theta,
                               const Vector<T, 3>& axis) -> numerics::Matrix<T> {
    const auto unitized = axis.unitVector();
    const double x = unitized[0], y = unitized[1], z = unitized[2];
    const double c = cos(theta), s = sin(theta), C = 1 - c;

    auto matx = numerics::Matrix<T>::zeros(3, 3);

    matx(0, 0) = x * x * C + c;
    matx(0, 1) = x * y * C - z * s;
    matx(0, 2) = x * z * C + y * s;

    matx(1, 0) = y * x * C + z * s;
    matx(1, 1) = y * y * C + c;
    matx(1, 2) = y * z * C - x * s;

    matx(2, 0) = z * x * C - y * s;
    matx(2, 1) = z * y * C + x * s;
    matx(2, 2) = z * z * C + c;

    return matx;
  };

  // Lambda to generate a random orthogonal normal vector to the input
  auto random_orthogonal = [](const Vector<T, 3>& input) -> Vector<T, 3> {
    // Pick a random direction orthogonal to the surface normal
    //  at the center of the disk
    Vector<T, 3> some_perp(
      {input[1] - input[2], input[2] - input[0], input[0] - input[1]});
    double theta = axom::utilities::random_real(0.0, 2.0 * M_PI);
    Vector<T, 3> new_direction = std::cos(theta) * some_perp.unitVector() +
      std::sin(theta) * Vector3D::cross_product(input, some_perp).unitVector();

    return new_direction;
  };

  // Lambda to generate an entirely random unit vector
  auto random_unit = []() -> Vector<T, 3> {
    double theta = axom::utilities::random_real(0.0, 2 * M_PI);
    double u = axom::utilities::random_real(-1.0, 1.0);
    return Vector<T, 3> {sin(theta) * sqrt(1 - u * u),
                         cos(theta) * sqrt(1 - u * u),
                         u};
  };

  // Rotation matrix for the patch
  numerics::Matrix<T> rotator;

  // Prefer to work with a trimmed patch
  NURBSPatch<T, 3> integralPatch = nPatchData.patch;

  // Define vector fields whose curl gives us the winding number
  detail::DiscontinuityAxis field_direction;
  bool extraTrimming = false;

  auto bBox = BoundingBox<T, 3>(nPatchData.bbox);
  auto characteristic_length = bBox.range().norm();
  bBox.expand(0.01 * characteristic_length);
  auto oBox =
    OrientedBoundingBox<T, 3>(nPatchData.obox).expand(0.01 * characteristic_length);

  // Case 1: Exterior without rotations
  if(!bBox.contains(query))
  {
    case_code = 0;
    integrated_trimming_curves = integralPatch.getNumTrimmingCurves();

    const bool exterior_x =
      bBox.getMin()[0] > query[0] || query[0] > bBox.getMax()[0];
    const bool exterior_y =
      bBox.getMin()[1] > query[1] || query[1] > bBox.getMax()[1];
    const bool exterior_z =
      bBox.getMin()[2] > query[2] || query[2] > bBox.getMax()[2];

    if(exterior_x || exterior_y)
    {
      field_direction = detail::DiscontinuityAxis::z;
    }
    else if(exterior_y || exterior_z)
    {
      field_direction = detail::DiscontinuityAxis::x;
    }
    else if(exterior_x || exterior_z)
    {
      field_direction = detail::DiscontinuityAxis::y;
    }
  }
  // Case 1.5: Exterior with rotation
  else if(!oBox.contains(query))
  {
    case_code = 1;
    integrated_trimming_curves = integralPatch.getNumTrimmingCurves();

    /* The following steps rotate the patch until the OBB is /not/ 
       directly above or below the query point */
    field_direction = detail::DiscontinuityAxis::rotated;

    // Find vector from query to the bounding box
    Point<T, 3> closest = closest_point(query, oBox);
    Vector<T, 3> v0 = Vector<T, 3>(query, closest).unitVector();

    // Find the direction of a ray perpendicular to that
    Vector<T, 3> v1;
    if(std::abs(v0[2]) > std::abs(v0[0]))
    {
      v1 = Vector<T, 3>({v0[2], v0[2], -v0[0] - v0[1]}).unitVector();
    }
    else
    {
      v1 = Vector<T, 3>({-v0[1] - v0[2], v0[0], v0[0]}).unitVector();
    }

    // Rotate v0 around v1 until it is perpendicular to the plane spanned by k and v1
    double ang = (v0[2] < 0 ? 1.0 : -1.0) *
      acos(axom::utilities::clampVal(
        -(v0[0] * v1[1] - v0[1] * v1[0]) / sqrt(v1[0] * v1[0] + v1[1] * v1[1]),
        -1.0,
        1.0));
    rotator = angleAxisRotMatrix(ang, v1);
  }
  // Case 2: Cast a ray, record intersections
  else
  {
    case_code = 2;
    integrated_trimming_curves = integralPatch.getNumTrimmingCurves();

    // wn_split.first += 0.0;
    // return wn_split;

    field_direction = detail::DiscontinuityAxis::rotated;
    Line<T, 3> discontinuity_axis(query, discontinuity_direction);

    T patch_knot_size = axom::utilities::max(
      integralPatch.getKnots_u()[integralPatch.getNumKnots_u() - 1] -
        integralPatch.getKnots_u()[0],
      integralPatch.getKnots_v()[integralPatch.getNumKnots_v() - 1] -
        integralPatch.getKnots_v()[0]);

    // Tolerance for what counts as "close to a boundary" in parameter space
    T disk_radius = 0.1 * patch_knot_size;

    // Compute intersections with the *untrimmed and extrapolated* patch
    axom::Array<T> up, vp, tp;
    bool isHalfOpen = false, isTrimmed = false;
    bool success = intersect(discontinuity_axis,
                             nPatchData,
                             tp,
                             up,
                             vp,
                             1e-10,  // This is a good heuristic value for accuracy
                             EPS,
                             isHalfOpen);

    if(!success)
    {
      // If too many intersections with the untrimmed patch are recorded, then either

      // 1. The ray is parallel to a the surface
      if(up.size() >= 20)
      {
        // Can't handle this case. Need to start over with a different cast direction
        auto new_direction = random_orthogonal(discontinuity_direction);

        // std::cout << "Recasting (very derogatory) at " << depth << std::endl;
        // std::cout << "Degenerate intersection point" << std::endl;
        return winding_number_casting_split(query,
                                            nPatchData.patch,
                                            new_direction,
                                            edge_tol,
                                            quad_tol,
                                            EPS,
                                            depth + 1);
      }

      // 2. The surface is degenerate, and so all intersections were recorded
      //  with the same t parameter, and pruned by `intersect()`.

      // Treating this case requires some implementation I don't have yet,
      //  namely clipping trimming curves along u/v isocurves

      // For now, if *any* of the reported intersection points are coincident with the query,
      //  just return a zero
      for(int i = 0; i < up.size(); ++i)
      {
        Point<T, 3> intersection_point = nPatchData.patch.evaluate(up[i], vp[i]);
        if(squared_distance(query, intersection_point) <= edge_tol_sq)
        {
          wn_split.first = 0.0;
          wn_split.second = 0.0;
          return wn_split;
        }
      }

      // Otherwise, we can cast again with a different direction
      {
        auto new_direction = random_orthogonal(discontinuity_direction);

        // std::cout << "Degenerate intersection point" << std::endl;
        // std::cout << "Recasting (derogatory) at " << depth << std::endl;
        return winding_number_casting_split(query,
                                            nPatchData.patch,
                                            new_direction,
                                            edge_tol,
                                            quad_tol,
                                            EPS);
      }
    }

    // Account for each discontinuity in the integrand on the *untrimmed and extrapolated* surface

    // If no intersection, then nothing to account for

    // Otherwise, account for each discontinuity analytically or through disk subdivision
    for(int i = 0; i < up.size(); ++i)
    {
      // Check for surface degeneracies or tangencies
      Vector<T, 3> the_normal = nPatchData.patch.normal(up[i], vp[i]);
      // std::cout << the_normal.norm() << " " << the_normal.unitVector().dot(discontinuity_direction) << std::endl;
      bool bad_intersection =
        axom::utilities::isNearlyEqual(the_normal.norm(), 0.0, EPS) ||
        axom::utilities::isNearlyEqual(
          the_normal.unitVector().dot(discontinuity_direction),
          0.0,
          EPS);

      // Check for surface coincidence
      Point<T, 3> intersection_point = nPatchData.patch.evaluate(up[i], vp[i]);
      bool isOnSurface =
        squared_distance(query, intersection_point) <= edge_tol_sq;

      if(bad_intersection && !isOnSurface)
      {
        // If a far away ray intersects the surface at a tangent/cusp,
        //  can recast and try again
        auto new_direction = random_unit();

        // std::cout << "Recasting: " << the_normal.norm() << ", "
        // << the_normal.unitVector().dot(discontinuity_direction)
        // << std::endl;
        // std::cout << "Recasting at " << depth << std::endl;
        return winding_number_casting_split(query,
                                            nPatchData.patch,
                                            new_direction,
                                            edge_tol,
                                            quad_tol,
                                            EPS,
                                            depth + 1);
      }

      if(isOnSurface)
      {
        // If the query point is on the surface, we need to consider a smaller disk
        //  to ensure its winding number is known to be near-zero
        disk_radius *= 0.1;
      }

      // This method accomplishes 2 tasks:
      //  > Determines if the disk of "safe" radius is *entirely* inside or outside trimming curves
      //  > If the disk intersects the trimming curves, performs disk subdivision
      bool isDiskInside, isDiskOutside, ignoreInteriorDisk = true;
      NURBSPatch<T, 3> disk_patch;
      int old_num_trim = integralPatch.getNumTrimmingCurves();

      {
        // AXOM_ANNOTATE_SCOPE("DISK_SPLIT");
        integralPatch.diskSplit(up[i],
                                vp[i],
                                disk_radius,
                                integralPatch,
                                disk_patch,
                                isDiskInside,
                                isDiskOutside,
                                ignoreInteriorDisk);

        // nPatchData.patch.printTrimmingCurves(
        //   "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\graphical_"
        //   "abstract\\trimming_curves\\original.txt");
        // disk_patch.printTrimmingCurves(
        //   "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\graphical_"
        //   "abstract\\trimming_curves\\disk.txt");
        // integralPatch.printTrimmingCurves(
        //   "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\graphical_"
        //   "abstract\\trimming_curves\\remaining.txt");

        // Extra trimming is applied if the disk is NOT inside and if the disk is NOT outside,
        // and if the disk is NOT ignored while inside
        extraTrimming = extraTrimming || (!isDiskInside && !isDiskOutside) ||
          (isDiskInside && !ignoreInteriorDisk);

        // --caliper report, counts
        // --caliper counts
      }

      if(extraTrimming)
      {
        case_code = 3;
        integrated_trimming_curves += disk_patch.getNumTrimmingCurves() +
          (integralPatch.getNumTrimmingCurves() - old_num_trim);
      }

      // If the query point is on the surface, the contribution of the disk is near-zero,
      //  and we only need to puncture the larger surface to proceed
      if(isOnSurface)
      {
        wn_split.first += 0.0;
        wn_split.second += 0.0;
      }
      // If the disk overlaps with a trimming curve
      else if(!isDiskInside && !isDiskOutside)
      {
        auto new_direction = random_orthogonal(the_normal);
        std::cout << std::setprecision(15) << std::endl;
        std::cout << new_direction[0] << "," << std::endl;
        std::cout << new_direction[1] << "," << std::endl;
        std::cout << new_direction[2] << "};" << std::endl;
        new_direction =
          Vector<T, 3> {0.183876962693649, -0.850254833796789, 0.493209874389953};

        // We compute the contribution of the disk directly,
        //  but with a different direction to avoid repeated subdivision
        // std::cout << "Disk Subdivision" << std::endl;
        // if(depth < 3)
        {
          // std::cout << "Disk dividing at " << depth << std::endl;
          wn_split.first += winding_number_casting(query,
                                                   disk_patch,
                                                   new_direction,
                                                   edge_tol,
                                                   quad_tol,
                                                   EPS,
                                                   depth + 1);
        }
      }
      // If the disk is entirely inside or outside, the jump condition is known
      else if(isDiskOutside)
      {
        wn_split.second += 0.0;
      }
      else if(isDiskInside)
      {
        Vector<T, 3> the_direction =
          Vector<T, 3>(query, intersection_point).unitVector();

        // Do a dot product to see what side of the surface the intersection is on
        wn_split.second += std::copysign(0.5, the_normal.dot(the_direction));
      }
    }

    // Rotate the patch so that the discontinuity direction is aligned with the z-axis
    Vector<T, 3> axis = {discontinuity_direction[1],
                         -discontinuity_direction[0],
                         0.0};

    double ang =
      acos(axom::utilities::clampVal(discontinuity_direction[2], -1.0, 1.0));

    rotator = angleAxisRotMatrix(ang, axis);

    // std::cout << wn_split.second << std::endl;
    // std::cout << std::endl;
  }

  if(extraTrimming)
  {
    // Can't use cached quadrature rules, cause it's unclear which ones to use

    //  Rotate it if we need to
    if(field_direction == detail::DiscontinuityAxis::rotated)
    {
      // The trimming curves for rotatedPatch have been changed as needed,
      //  but we need to rotate the control points
      auto patch_shape = integralPatch.getControlPoints().shape();
      for(int i = 0; i < patch_shape[0]; ++i)
      {
        for(int j = 0; j < patch_shape[1]; ++j)
        {
          integralPatch(i, j) =
            detail::rotate_point(rotator, query, nPatchData.patch(i, j));
        }
      }
    }

    wn_split.first += detail::stokes_winding_number(query,
                                                    integralPatch,
                                                    field_direction,
                                                    quad_npts,
                                                    quad_tol);
  }
  else
  {
    // It's easier if we don't need to rotate the patch
    if(field_direction != detail::DiscontinuityAxis::rotated)
    {
      wn_split.first += detail::stokes_winding_number_cached(query,
                                                             nPatchData,
                                                             field_direction,
                                                             quad_npts,
                                                             quad_tol);
    }
    else
    {
      wn_split.first += detail::stokes_winding_number_cached_rotated(query,
                                                                     nPatchData,
                                                                     rotator,
                                                                     quad_npts,
                                                                     quad_tol);
    }
  }

  return wn_split;
}

template <typename T>
double simple_coincident_wn(const Point<T, 3>& query,
                            const NURBSPatch<T, 3>& nPatch,
                            const Vector<T, 3>& discontinuity_direction,
                            const double quad_tol = 1e-8)
{
  // Fix the number of quadrature nodes arbitrarily
  constexpr int quad_npts = 15;

  double wn = 0.0;

  /* 
   * To use Stokes theorem, we need to identify either a line containing the
   * query that does not intersect the surface, or one that intersects the *interior*
   * of the surface at known locations.
   */

  // Lambda to generate a 3D rotation matrix from an angle and axis
  // Formulation from https://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
  auto angleAxisRotMatrix = [](double theta,
                               const Vector<T, 3>& axis) -> numerics::Matrix<T> {
    const auto unitized = axis.unitVector();
    const double x = unitized[0], y = unitized[1], z = unitized[2];
    const double c = cos(theta), s = sin(theta), C = 1 - c;

    auto matx = numerics::Matrix<T>::zeros(3, 3);

    matx(0, 0) = x * x * C + c;
    matx(0, 1) = x * y * C - z * s;
    matx(0, 2) = x * z * C + y * s;

    matx(1, 0) = y * x * C + z * s;
    matx(1, 1) = y * y * C + c;
    matx(1, 2) = y * z * C - x * s;

    matx(2, 0) = z * x * C - y * s;
    matx(2, 1) = z * y * C + x * s;
    matx(2, 2) = z * z * C + c;

    return matx;
  };

  // Rotation matrix for the patch
  numerics::Matrix<T> rotator;

  // Prefer to work with a trimmed patch
  NURBSPatch<T, 3> integralPatch = nPatch;

  // Define vector fields whose curl gives us the winding number
  detail::DiscontinuityAxis field_direction;
  bool extraTrimming = false;

  auto bBox = nPatch.boundingBox();
  auto characteristic_length = bBox.range().norm();
  bBox.expand(0.01 * characteristic_length);
  auto oBox = nPatch.orientedBoundingBox().expand(0.01 * characteristic_length);

  // Case 2: Cast a ray, record intersections
  {
    field_direction = detail::DiscontinuityAxis::rotated;
    Line<T, 3> discontinuity_axis(query, discontinuity_direction);

    T patch_knot_size = axom::utilities::max(
      integralPatch.getKnots_u()[integralPatch.getNumKnots_u() - 1] -
        integralPatch.getKnots_u()[0],
      integralPatch.getKnots_v()[integralPatch.getNumKnots_v() - 1] -
        integralPatch.getKnots_v()[0]);

    // Rotate the patch so that the discontinuity direction is aligned with the z-axis
    Vector<T, 3> axis = {discontinuity_direction[1],
                         -discontinuity_direction[0],
                         0.0};

    double ang =
      acos(axom::utilities::clampVal(discontinuity_direction[2], -1.0, 1.0));

    rotator = angleAxisRotMatrix(ang, axis);
  }

  //  Rotate it if we need to
  if(field_direction == detail::DiscontinuityAxis::rotated)
  {
    // The trimming curves for rotatedPatch have been changed as needed,
    //  but we need to rotate the control points
    auto patch_shape = integralPatch.getControlPoints().shape();
    for(int i = 0; i < patch_shape[0]; ++i)
    {
      for(int j = 0; j < patch_shape[1]; ++j)
      {
        integralPatch(i, j) = detail::rotate_point(rotator, query, nPatch(i, j));
      }
    }
  }

  wn += detail::stokes_winding_number(query,
                                      integralPatch,
                                      field_direction,
                                      quad_npts,
                                      0);

  return wn;
}

template <typename T>
double sphere_winding_number_casting(const Point<T, 3>& query,
                                     const NURBSPatchData<T>& nPatchData,
                                     int& case_code,
                                     int& integrated_trimming_curves,
                                     const double edge_tol = 1e-8,
                                     const double quad_tol = 1e-8,
                                     const double EPS = 1e-8,
                                     const int depth = 0)
{
  const double edge_tol_sq = edge_tol * edge_tol;

  // We know the surface is parts of a sphere, so pick the discontinuity direction advantageously
  const Vector<T, 3>& discontinuity_direction = {0, 0, 1};

  // Fix the number of quadrature nodes arbitrarily
  constexpr int quad_npts = 20;

  // The first is the GWN from stokes, the second is the jump condition
  double wn = 0.0;

  /* 
   * To use Stokes theorem, we need to identify either a line containing the
   * query that does not intersect the surface, or one that intersects the *interior*
   * of the surface at known locations.
   */

  // Lambda to generate a 3D rotation matrix from an angle and axis
  // Formulation from https://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
  auto angleAxisRotMatrix = [](double theta,
                               const Vector<T, 3>& axis) -> numerics::Matrix<T> {
    const auto unitized = axis.unitVector();
    const double x = unitized[0], y = unitized[1], z = unitized[2];
    const double c = cos(theta), s = sin(theta), C = 1 - c;

    auto matx = numerics::Matrix<T>::zeros(3, 3);

    matx(0, 0) = x * x * C + c;
    matx(0, 1) = x * y * C - z * s;
    matx(0, 2) = x * z * C + y * s;

    matx(1, 0) = y * x * C + z * s;
    matx(1, 1) = y * y * C + c;
    matx(1, 2) = y * z * C - x * s;

    matx(2, 0) = z * x * C - y * s;
    matx(2, 1) = z * y * C + x * s;
    matx(2, 2) = z * z * C + c;

    return matx;
  };

  // Lambda to generate a random orthogonal normal vector to the input
  auto random_orthogonal = [](const Vector<T, 3>& input) -> Vector<T, 3> {
    // Pick a random direction orthogonal to the surface normal
    //  at the center of the disk
    Vector<T, 3> some_perp(
      {input[1] - input[2], input[2] - input[0], input[0] - input[1]});
    double theta = axom::utilities::random_real(0.0, 2.0 * M_PI);
    Vector<T, 3> new_direction = std::cos(theta) * some_perp.unitVector() +
      std::sin(theta) * Vector3D::cross_product(input, some_perp).unitVector();

    return new_direction;
  };

  // Lambda to generate an entirely random unit vector
  auto random_unit = []() -> Vector<T, 3> {
    double theta = axom::utilities::random_real(0.0, 2 * M_PI);
    double u = axom::utilities::random_real(-1.0, 1.0);
    return Vector<T, 3> {sin(theta) * sqrt(1 - u * u),
                         cos(theta) * sqrt(1 - u * u),
                         u};
  };

  // Rotation matrix for the patch
  numerics::Matrix<T> rotator;

  // Prefer to work with a trimmed patch
  NURBSPatch<T, 3> integralPatch = nPatchData.patch;

  // Define vector fields whose curl gives us the winding number
  detail::DiscontinuityAxis field_direction;
  bool extraTrimming = false;

  auto bBox = BoundingBox<T, 3>(nPatchData.bbox);
  auto characteristic_length = bBox.range().norm();
  bBox.expand(0.01 * characteristic_length);
  auto oBox =
    OrientedBoundingBox<T, 3>(nPatchData.obox).expand(0.01 * characteristic_length);

  // Case 2: Cast a ray, record intersections
  if(true)
  {
    case_code = 2;
    integrated_trimming_curves = integralPatch.getNumTrimmingCurves();

    // wn_split.first += 0.0;
    // return wn_split;

    field_direction = detail::DiscontinuityAxis::rotated;
    Line<T, 3> discontinuity_axis(query, discontinuity_direction);

    // Compute intersections with the *untrimmed and extrapolated* patch
    axom::Array<T> up, vp, tp;
    bool isHalfOpen = false, isTrimmed = false;
    bool success = intersect(discontinuity_axis,
                             nPatchData,
                             tp,
                             up,
                             vp,
                             1e-10,  // This is a good heuristic value for accuracy
                             EPS,
                             isHalfOpen);

    // Account for each discontinuity in the integrand on the *untrimmed and extrapolated* surface

    // If no intersection, then nothing to account for

    // Otherwise, account for each discontinuity analytically or through disk subdivision
    for(int i = 0; i < up.size(); ++i)
    {
      // Check for surface degeneracies or tangencies
      Vector<T, 3> the_normal = nPatchData.patch.normal(up[i], vp[i]);

      // Check for surface coincidence
      Point<T, 3> intersection_point = nPatchData.patch.evaluate(up[i], vp[i]);
      bool isOnSurface =
        squared_distance(query, intersection_point) <= edge_tol_sq;

      // If the query point is on the surface, the contribution of the disk is near-zero,
      //  and we only need to puncture the larger surface to proceed
      if(isOnSurface)
      {
        wn += 0.0;
      }
      else
      {
        Vector<T, 3> the_direction =
          Vector<T, 3>(query, intersection_point).unitVector();

        // Do a dot product to see what side of the surface the intersection is on
        wn += std::copysign(0.5, the_normal.dot(the_direction));
      }
    }

    // Rotate the patch so that the discontinuity direction is aligned with the z-axis
    Vector<T, 3> axis = {discontinuity_direction[1],
                         -discontinuity_direction[0],
                         0.0};

    double ang =
      acos(axom::utilities::clampVal(discontinuity_direction[2], -1.0, 1.0));

    rotator = angleAxisRotMatrix(ang, axis);
  }

  if(true)
  {
    {
      wn += detail::stokes_winding_number_cached_rotated_adaptless(query,
                                                                   nPatchData,
                                                                   rotator,
                                                                   quad_npts,
                                                                   quad_tol);
    }
  }

  return wn;
}
#endif

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_WINDING_NUMBER_3D_HPP_
