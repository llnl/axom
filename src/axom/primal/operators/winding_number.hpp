// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file winding_number.hpp
 *
 * \brief Consists of methods to compute the generalized winding number (GWN) 
 *        for points with respect to various geometric objects.
 */

#ifndef AXOM_PRIMAL_WINDING_NUMBER_HPP_
#define AXOM_PRIMAL_WINDING_NUMBER_HPP_

// Axom includes
#include "axom/core.hpp"
#include "axom/config.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/Polyhedron.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/NURBSCurve.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/geometry/NURBSPatch.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"

#include "axom/primal/operators/detail/winding_number_2d_impl.hpp"
#include "axom/primal/operators/detail/winding_number_3d_impl.hpp"

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
/*
 * \brief Compute the GWN for a 2D point wrt a 2D line segment
 *
 * \param [in] q The query point to test
 * \param [in] s The line segment
 * \param [in] edge_tol The tolerance at which a point is on the line
 *
 * \return The GWN
 */
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const Segment<T, 2>& s,
                      double edge_tol = 1e-8)
{
  bool dummy_isOnEdge = false;
  return detail::linear_winding_number(q, s[0], s[1], dummy_isOnEdge, edge_tol);
}

/*
 * \brief Compute the winding number for a 2D point wrt a 2D triangle
 *
 * \param [in] q The query point to test
 * \param [in] tri The triangle
 * \param [in] includeBoundary If true, points on the boundary are considered interior.
 * \param [in] edge_tol The tolerance at which a point is on the line
 *
 * The triangle is assumed to be closed, so the winding number is an integer
 * 
 * \return The integer winding number
 */
template <typename T>
int winding_number(const Point<T, 2>& q,
                   const Triangle<T, 2>& tri,
                   bool includeBoundary = false,
                   double edge_tol = 1e-8)
{
  return winding_number(
    q,
    Polygon<T, 2>(axom::Array<Point<T, 2>>({tri[0], tri[1], tri[2]})),
    includeBoundary,
    edge_tol);
}

/*!
 * \brief Computes the winding number for a 2D point wrt a 2D polygon
 *
 * \param [in] R The query point to test
 * \param [in] P The Polygon object to test for containment
 * \param [out] isOnEdge An optional return parameter if the point is on the boundary
 * \param [in] includeBoundary If true, points on the boundary are considered interior
 * \param [in] edge_tol The distance at which a point is considered on the boundary
 * 
 * \return The integer winding number
 */
template <typename T>
int winding_number(const Point<T, 2>& R,
                   const Polygon<T, 2>& P,
                   bool& isOnEdge,
                   bool includeBoundary = false,
                   double edge_tol = 1e-8)
{
  return detail::polygon_winding_number(R, P, isOnEdge, includeBoundary, edge_tol);
}

/*!
 * \brief Computes the winding number for a 2D point wrt a 2D polygon
 *
 * \param [in] R The query point to test
 * \param [in] P The Polygon object to test for containment
 * \param [in] includeBoundary If true, points on the boundary are considered interior
 * \param [in] edge_tol The distance at which a point is considered on the boundary
 * 
 * Computes the integer winding number for a polygon without an additional
 *  return parameter for whether the point is on the boundary.
 * 
 * \return The integer winding number
 */
template <typename T>
int winding_number(const Point<T, 2>& R,
                   const Polygon<T, 2>& P,
                   bool includeBoundary = false,
                   double edge_tol = 1e-8)
{
  bool isOnEdge = false;
  return detail::polygon_winding_number(R, P, isOnEdge, includeBoundary, edge_tol);
}

/*!
 * \brief Computes the GWN for a 2D point wrt a 2D NURBS curve
 *
 * \param [in] query The query point to test
 * \param [in] n The NURBS curve object 
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const NURBSCurve<T, 2>& n,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  bool dummy_isOnCurve = false;
  return detail::nurbs_winding_number(q, n, dummy_isOnCurve, edge_tol, EPS);
}

/*!
 * \brief Computes the GWN for a 2D point wrt a 2D NURBS curve
 *
 * \param [in] query The query point to test
 * \param [in] bezier The Bezier curve object 
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const BezierCurve<T, 2>& bezier,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  bool dummy_isOnCurve = false;
  return detail::bezier_winding_number(q, bezier, dummy_isOnCurve, edge_tol, EPS);
}

/*!
 * \brief Computes the GWN for a 2D point wrt to a 2D curved polygon
 *
 * \param [in] query The query point to test
 * \param [in] cpoly The CurvedPolygon object
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the GWN for the curved polygon by summing the GWN for each curved edge
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const CurvedPolygon<T, 2>& cpoly,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  bool dummy_isOnCurve = false;

  double ret_val = 0.0;
  for(int i = 0; i < cpoly.numEdges(); i++)
  {
    ret_val +=
      detail::bezier_winding_number(q, cpoly[i], dummy_isOnCurve, edge_tol, EPS);
  }

  return ret_val;
}

/*!
 * \brief Computes the GWN for a 2D point wrt to a collection of 2D Bezier curves
 *
 * \param [in] query The query point to test
 * \param [in] carray The array of Bezier curves
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Sums the GWN at `query` for each curved edge
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const axom::Array<BezierCurve<T, 2>>& carray,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  bool dummy_isOnCurve = false;
  double ret_val = 0.0;
  for(int i = 0; i < carray.size(); i++)
  {
    ret_val +=
      detail::bezier_winding_number(q, carray[i], dummy_isOnCurve, edge_tol, EPS);
  }

  return ret_val;
}

/*!
 * \brief Computes the GWN for a 2D point wrt to a collection of 2D NURBS curves
 *
 * \param [in] query The query point to test
 * \param [in] narray The array of NURBS curves
 * \param [in] edge_tol The physical distance level at which objects are considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Sums the GWN at `query` for each curved edge
 * 
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const axom::Array<NURBSCurve<T, 2>>& narray,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  bool dummy_isOnCurve = false;
  double ret_val = 0.0;
  for(int i = 0; i < narray.size(); i++)
  {
    ret_val +=
      detail::nurbs_winding_number(q, narray[i], dummy_isOnCurve, edge_tol, EPS);
  }

  return ret_val;
}

//@{
//! @name Winding number operations between 3D points and primitives

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

// #ifdef AXOM_USE_MFEM

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
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 3>& query,
                      const BezierPatch<T, 3>& bPatch,
                      const double edge_tol = 1e-8,
                      const double quad_tol = 1e-8,
                      const double EPS = 1e-8)
{
  NURBSPatch<T, 3> nPatch(bPatch);
  nPatch.makeTriviallyTrimmed();
  return winding_number(query, nPatch, edge_tol, quad_tol, EPS);
}

/*
 * \brief Computes the GWN for a 3D point wrt a 3D NURBS patch
 *
 * \param [in] query The query point to test
 * \param [in] nPatch The NURBS patch object
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] quad_tol The maximum relative error allowed in the quadrature
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * \param [in] depth The current recursive depth
 * 
 * Computes the generalized winding number for a NURBS patch using Stokes theorem.
 *
 * \return The GWN.
 */
template <typename T>
double winding_number(const Point<T, 3>& query,
                      const NURBSPatch<T, 3>& nPatch,
                      const double edge_tol = 1e-8,
                      const double quad_tol = 1e-8,
                      const double EPS = 1e-8,
                      const int depth = 0)
{
  const double edge_tol_sq = edge_tol * edge_tol;

  // Fix the number of quadrature points arbitrarily
  constexpr int quad_npts = 15;

  // Store the winding number
  double the_gwn = 0.0;

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
  auto rotate_point = [&query](const numerics::Matrix<T>& matx,
                               const Point<T, 3> input) -> Point<T, 3> {
    Vector<T, 3> shifted(query, input);
    Vector<T, 3> rotated;
    numerics::matrix_vector_multiply(matx, shifted.data(), rotated.data());
    return Point<T, 3>(
      {rotated[0] + query[0], rotated[1] + query[1], rotated[2] + query[2]});
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

  // Make the patch trivially trimmed, if necessary
  NURBSPatch<T, 3> nPatchTrimmedMore(nPatch);
  if(!nPatch.isTrimmed())
  {
    nPatchTrimmedMore.makeTriviallyTrimmed();
  }

  // Define vector fields whose curl gives us the winding number
  detail::DiscontinuityAxis field_direction;

  // Generate slightly expanded bounding boxes
  auto bBox = nPatch.boundingBox();
  auto oBox = nPatch.orientedBoundingBox();

  auto characteristic_length = bBox.range().norm();

  bBox.expand(0.01 * characteristic_length);
  oBox.expand(0.01 * characteristic_length);

  // Case 1: Exterior without rotations
  if(!bBox.contains(query))
  {
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
  else
  {
    field_direction = detail::DiscontinuityAxis::rotated;
    Vector<T, 3> discontinuity_direction = random_unit();
    Line<T, 3> discontinuity_axis(query, discontinuity_direction);

    T patch_knot_size =
      axom::utilities::max(nPatch.getMaxKnot_u() - nPatch.getMinKnot_u(),
                           nPatch.getMaxKnot_v() - nPatch.getMinKnot_v());

    // Tolerance for what counts as "close to a boundary" in parameter space
    T disk_radius = 0.1 * patch_knot_size;

    // Compute intersections with the *untrimmed and extrapolated* patch
    axom::Array<T> up, vp, tp;
    bool isHalfOpen = false, isTrimmed = false;
    bool success = intersect(discontinuity_axis,
                             nPatchTrimmedMore,
                             tp,
                             up,
                             vp,
                             1e-10,  // This is a good heuristic value for accuracy
                             EPS,
                             isHalfOpen);

    if(!success)
    {
      // This is the part where we haven't figured it out exactly
    }

    // If no intersections are recorded, then nothing extra to account for

    // Otherwise, account for each discontinuity analytically,
    //  or recursively through disk subdivision
    for(int i = 0; i < up.size(); ++i)
    {
      // Compute the intersection point on the surface
      Point<T, 3> the_point( nPatch.evaluate(up[i], vp[i]) );
      Vector<T, 3> the_normal = nPatch.normal(up[i], vp[i]);

      // Check for bad intersections, i.e.,
      //  > There normal is poorly defined (cusp)
      //  > The normal is tangent to the axis of discontinuity
      bool bad_intersection =
        axom::utilities::isNearlyEqual(the_normal.norm(), 0.0, EPS) ||
        axom::utilities::isNearlyEqual(
          the_normal.unitVector().dot(discontinuity_direction),
          0.0,
          EPS);

      bool isOnSurface = squared_distance(query, the_point) <= edge_tol_sq;

      if(bad_intersection && !isOnSurface)
      {
        // If a non-coincident ray intersects the surface at a tangent/cusp,
        //  can recast and try again
        return winding_number(query, nPatch, edge_tol, quad_tol, EPS, depth + 1);
      }

      if(isOnSurface)
      {
        // If the query point is on the surface, then shrink the disk
        //  to ensure its winding number is known to be near-zero
        disk_radius = 0.1 * disk_radius;
      }

      // Consider a disk around the intersection point via NURBSPatch::diskSplit.
      //   If the disk intersects any trimming curves, need to do disk subdivision.
      //   If not, we can compute the winding number without changing the trimming curvse
      const bool ignoreInteriorDisk = true, clipDisk = true;
      bool isDiskInside, isDiskOutside;
      NURBSPatch<T, 3> the_disk;
      nPatchTrimmedMore.diskSplit(up[i],
                                  vp[i],
                                  disk_radius,
                                  the_disk,
                                  nPatchTrimmedMore,
                                  isDiskInside,
                                  isDiskOutside,
                                  ignoreInteriorDisk,
                                  clipDisk);

      if(isOnSurface)
      {
        // If the query point is on the surface, the contribution of the disk is near-zero
        //  and we only needed to puncture the larger surface to proceed
        continue;
      }
      else if(!isDiskInside && !isDiskOutside)
      {
        // If the disk overlapped with the trimming curves, evaluate the winding number for the disk
        the_gwn +=
          winding_number(query, the_disk, edge_tol, quad_tol, EPS, depth + 1);
      }
      else if(isDiskOutside)
      {
        // If the disk is entirely outside the trimming curves, can just look at the boundary
        continue;
      }
      else if(isDiskInside)
      {
        // If the disk is entirely inside the trimming curves,
        //  need to account for the scalar field discontinuity
        auto the_direction = Vector<T, 3>(query, the_point).unitVector();
        the_gwn += std::copysign(0.5, the_normal.dot(the_direction));
      }
    }

    // Rotate the patch so that the discontinuity is aligned with the z-axis
    Vector<T, 3> axis = {discontinuity_direction[1],
                         -discontinuity_direction[0],
                         0};

    double ang =
      std::acos(axom::utilities::clampVal(discontinuity_direction[2], -1.0, 1.0));

    rotator = angleAxisRotMatrix(ang, axis);
  }

  if(field_direction == detail::DiscontinuityAxis::rotated)
  {
    // The trimming curves for rotatedPatch have been changed as needed,
    //  but we need to rotate the control points
    auto patch_shape = nPatchTrimmedMore.getControlPoints().shape();
    for(int i = 0; i < patch_shape[0]; ++i)
    {
      for(int j = 0; j < patch_shape[1]; ++j)
      {
        nPatchTrimmedMore(i, j) = rotate_point(rotator, nPatch(i, j));
      }
    }
  }

  the_gwn += detail::stokes_winding_number_evaluate(query,
                                                    nPatchTrimmedMore,
                                                    field_direction,
                                                    quad_npts,
                                                    quad_tol);

  return the_gwn;
}
// #endif
//@}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_WINDING_NUMBER_H_
