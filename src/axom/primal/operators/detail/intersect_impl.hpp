// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file intersect_impl.hpp
 *
 * This file provides helper functions for testing whether
 * geometric primitives intersect
 */

#ifndef AXOM_PRIMAL_INTERSECT_IMPL_HPP_
#define AXOM_PRIMAL_INTERSECT_IMPL_HPP_

#include "axom/core/Macros.hpp"
#include "axom/core/numerics/Determinants.hpp"
#include "axom/core/utilities/Utilities.hpp"

#include "axom/primal/operators/detail/fuzzy_comparators.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"
#include "axom/primal/geometry/Plane.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/Ray.hpp"
#include "axom/primal/geometry/Line.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"

namespace axom
{
namespace primal
{
namespace detail
{
//--------------------TYPE ALIASES AND FUNCTION DECLARATIONS ------------------

using Vector3 = primal::Vector<double, 3>;
using Point3 = primal::Point<double, 3>;
using Triangle3 = primal::Triangle<double, 3>;
using Triangle2 = primal::Triangle<double, 2>;
using Point2 = primal::Point<double, 2>;

AXOM_HOST_DEVICE
bool nonzeroSignMatch(double x, double y, double z, double EPS = 1E-12);

AXOM_HOST_DEVICE
bool twoZeros(double x, double y, double z, double EPS = 1E-12);

AXOM_HOST_DEVICE
bool oneZeroOthersMatch(double x, double y, double z, double EPS = 1E-12);

AXOM_HOST_DEVICE
int countZeros(double x, double y, double z, double EPS = 1E-12);

AXOM_HOST_DEVICE
double twoDcross(const Point2& A, const Point2& B, const Point2& C);

/*!
 * This function finds where p1 lies in relation to the vertices of t2
 * and calls either checkEdge() or checkVertex().
 * \return status true iff triangle p1 q1 r1 intersects triangle p2 q2 r2.
 */
AXOM_HOST_DEVICE
bool intersectPermuted2DTriangles(const Point2& p1,
                                  const Point2& q1,
                                  const Point2& r1,
                                  const Point2& p2,
                                  const Point2& q2,
                                  const Point2& r2,
                                  bool includeBoundary,
                                  double EPS);

/*!
 * Triangle 1 vertices have been permuted to CCW: permute t2 to CCW
 * and call worker function to test for intersection.
 *
 * q1 and r1 both lie in the negative half-space defined by t2; p1 lies in
 * t2's plane or in its positive half-space.
 * The sign of dp2, dq2, and dr2 indicates whether the associated vertex
 * of t2 lies in the positive or negative half-space defined by t1.
 */
AXOM_HOST_DEVICE
bool intersectOnePermutedTriangle(const Point3& p1,
                                  const Point3& q1,
                                  const Point3& r1,
                                  const Point3& p2,
                                  const Point3& q2,
                                  const Point3& r2,
                                  double dp2,
                                  double dq2,
                                  double dr2,
                                  Vector3& normal,
                                  bool includeBoundary,
                                  double EPS);

AXOM_HOST_DEVICE
bool intersectTwoPermutedTriangles(const Point3& p1,
                                   const Point3& q1,
                                   const Point3& r1,
                                   const Point3& p2,
                                   const Point3& q2,
                                   const Point3& r2,
                                   bool includeBoundary,
                                   double EPS);

/*!
 * Project (nearly) coplanar triangles 1 and 2 on an axis; call 2D worker
 * function to test for intersection.
 */
AXOM_HOST_DEVICE
bool intersectCoplanar3DTriangles(const Point3& p1,
                                  const Point3& q1,
                                  const Point3& r1,
                                  const Point3& p2,
                                  const Point3& q2,
                                  const Point3& r2,
                                  Vector3 normal,
                                  bool includeBoundary,
                                  double EPS);

/*!
 * \brief Orients triangle vertices so both triangles are CCW; calls worker.
 * \return status true iff t1 and t2 intersect
 *
 * Determine triangle orientation, then call the worker function with
 * vertices from t1 and t2 permuted to ensure CCW orientation.
 */
AXOM_HOST_DEVICE
bool TriangleIntersection2D(const Triangle2& t1, const Triangle2& t2, bool includeBoundary, double EPS);

//------------------------------ IMPLEMENTATIONS ------------------------------

/*! @{ @name 3D triangle-triangle intersection */

/*!
 * \brief Tests if 3D Triangles t1 and t2 intersect.
 * \return status true iff t1 intersects with t2, otherwise, false.
 *
 * This algorithm is modeled after Devillers and Guigue (2002).  It computes
 * the line of intersection L of the triangles' planes and finds the segments
 * S1, S2 where t1 and t2 intersect L.  If those segments intersect, the
 * triangles must intersect.  Note that edge and vertex intersections are
 * reported as hits: an edge or vertex intersecting any part of another
 * triangle causes a return value of true.  This is consistent with the paper.
 *
 * Coplanar triangles are handled with a decision tree, testing the
 * relative position of triangle vertices.

 * Olivier Devillers and Phillipe Guigue, Faster Triangle-Triangle Intersection
 * Tests, RR-4488, INRIA (2002).  https://hal.inria.fr/inria-00072100/
 */
template <typename T>
AXOM_HOST_DEVICE bool intersect_tri3D_tri3D(const Triangle<T, 3>& t1,
                                            const Triangle<T, 3>& t2,
                                            bool includeBoundary,
                                            double EPS)
{
  SLIC_CHECK_MSG(!t1.degenerate(), "\n\n WARNING \n\n Triangle " << t1 << " is degenerate");
  SLIC_CHECK_MSG(!t2.degenerate(), "\n\n WARNING \n\n Triangle " << t2 << " is degenerate");

  // Step 1: Check if all the vertices of triangle 1 lie on the same side of
  // the plane created by triangle 2:

  // Vector3 t2Normal = Vector3::cross_product(Vector3(t2[2], t2[0]),
  //                                           Vector3(t2[2], t2[1]));
  Vector3 t2Normal = t2.normal().unitVector();
  const double dp1 = (t1[0] - t2[2]).dot(t2Normal);
  const double dq1 = (t1[1] - t2[2]).dot(t2Normal);
  const double dr1 = (t1[2] - t2[2]).dot(t2Normal);

  if(nonzeroSignMatch(dp1, dq1, dr1, EPS))
  {
    return false;
  }

  if(!includeBoundary && (twoZeros(dp1, dq1, dr1, EPS) || oneZeroOthersMatch(dp1, dq1, dr1, EPS)))
  {
    return false;
  }

  // Step 2: Check if all the vertices of triangle 2 lie on the same side of
  // the plane created by triangle 1:

  // Vector3 t1Normal = Vector3::cross_product(Vector3(t1[0], t1[1]),
  //                                           Vector3(t1[0], t1[2]));
  Vector3 t1Normal = t1.normal().unitVector();
  const double dp2 = (t2[0] - t1[2]).dot(t1Normal);
  const double dq2 = (t2[1] - t1[2]).dot(t1Normal);
  const double dr2 = (t2[2] - t1[2]).dot(t1Normal);

  if(nonzeroSignMatch(dp2, dq2, dr2, EPS))
  {
    return false;
  }

  if(!includeBoundary && (twoZeros(dp2, dq2, dr2, EPS) || oneZeroOthersMatch(dp2, dq2, dr2, EPS)))
  {
    return false;
  }

  /* Note: Because we know that all the vertices either triangle do not
     lay on the same side of the plane formed by the other triangle, we
     know that for each triangle, exactly 1 out of 3 points exists on one
     side of the plane formed by the other triangle.


     Step 3: We apply a circular permutation of triangle 1 such that its
     first point is the only point on the triangle that lies on one side of
     the plane formed by triangle 2 (with the other 2 on the other side),
     while handling the special case of one of the vertices lying on the
     plane formed by triangle 2.  We then perform a swap operation on the
     second and third points of triangle 2 to map the first point of
     triangle 1 to the positive halfspace formed by triangle 2's plane.
   */

  // compare the signs to create a convenient permutation of the vertices
  // of triangle 1

  if(isGt(dp1, 0.0, EPS))
  {
    if(isGt(dq1, 0.0, EPS))
    {
      return intersectOnePermutedTriangle(t1[2],
                                          t1[0],
                                          t1[1],
                                          t2[0],
                                          t2[2],
                                          t2[1],
                                          dp2,
                                          dr2,
                                          dq2,
                                          t1Normal,
                                          includeBoundary,
                                          EPS);
    }
    else if(isGt(dr1, 0.0, EPS))
    {
      return intersectOnePermutedTriangle(t1[1],
                                          t1[2],
                                          t1[0],
                                          t2[0],
                                          t2[2],
                                          t2[1],
                                          dp2,
                                          dr2,
                                          dq2,
                                          t1Normal,
                                          includeBoundary,
                                          EPS);
    }
    else
    {
      return intersectOnePermutedTriangle(t1[0],
                                          t1[1],
                                          t1[2],
                                          t2[0],
                                          t2[1],
                                          t2[2],
                                          dp2,
                                          dq2,
                                          dr2,
                                          t1Normal,
                                          includeBoundary,
                                          EPS);
    }
  }
  else if(isLt(dp1, 0.0, EPS))
  {
    if(isLt(dq1, 0.0, EPS))
    {
      return intersectOnePermutedTriangle(t1[2],
                                          t1[0],
                                          t1[1],
                                          t2[0],
                                          t2[1],
                                          t2[2],
                                          dp2,
                                          dq2,
                                          dr2,
                                          t1Normal,
                                          includeBoundary,
                                          EPS);
    }
    else if(isLt(dr1, 0.0, EPS))
    {
      return intersectOnePermutedTriangle(t1[1],
                                          t1[2],
                                          t1[0],
                                          t2[0],
                                          t2[1],
                                          t2[2],
                                          dp2,
                                          dq2,
                                          dr2,
                                          t1Normal,
                                          includeBoundary,
                                          EPS);
    }
    else
    {
      return intersectOnePermutedTriangle(t1[0],
                                          t1[1],
                                          t1[2],
                                          t2[0],
                                          t2[2],
                                          t2[1],
                                          dp2,
                                          dr2,
                                          dq2,
                                          t1Normal,
                                          includeBoundary,
                                          EPS);
    }
  }
  else  //dp1 ~= 0
  {
    if(isLt(dq1, 0.0, EPS))
    {
      if(isGeq(dr1, 0.0, EPS))
      {
        return intersectOnePermutedTriangle(t1[1],
                                            t1[2],
                                            t1[0],
                                            t2[0],
                                            t2[2],
                                            t2[1],
                                            dp2,
                                            dr2,
                                            dq2,
                                            t1Normal,
                                            includeBoundary,
                                            EPS);
      }
      else
      {
        return intersectOnePermutedTriangle(t1[0],
                                            t1[1],
                                            t1[2],
                                            t2[0],
                                            t2[1],
                                            t2[2],
                                            dp2,
                                            dq2,
                                            dr2,
                                            t1Normal,
                                            includeBoundary,
                                            EPS);
      }
    }
    else if(isGt(dq1, 0.0, EPS))
    {
      if(isGt(dr1, 0.0, EPS))
      {
        return intersectOnePermutedTriangle(t1[0],
                                            t1[1],
                                            t1[2],
                                            t2[0],
                                            t2[2],
                                            t2[1],
                                            dp2,
                                            dr2,
                                            dq2,
                                            t1Normal,
                                            includeBoundary,
                                            EPS);
      }
      else
      {
        return intersectOnePermutedTriangle(t1[1],
                                            t1[2],
                                            t1[0],
                                            t2[0],
                                            t2[1],
                                            t2[2],
                                            dp2,
                                            dq2,
                                            dr2,
                                            t1Normal,
                                            includeBoundary,
                                            EPS);
      }
    }
    else
    {
      if(isGt(dr1, 0.0, EPS))
      {
        return intersectOnePermutedTriangle(t1[2],
                                            t1[0],
                                            t1[1],
                                            t2[0],
                                            t2[1],
                                            t2[2],
                                            dp2,
                                            dq2,
                                            dr2,
                                            t1Normal,
                                            includeBoundary,
                                            EPS);
      }
      else if(isLt(dr1, 0.0, EPS))
      {
        return intersectOnePermutedTriangle(t1[2],
                                            t1[0],
                                            t1[1],
                                            t2[0],
                                            t2[2],
                                            t2[1],
                                            dp2,
                                            dr2,
                                            dq2,
                                            t1Normal,
                                            includeBoundary,
                                            EPS);
      }
      else
      {
        return intersectCoplanar3DTriangles(t1[0],
                                            t1[1],
                                            t1[2],
                                            t2[0],
                                            t2[1],
                                            t2[2],
                                            t1Normal,
                                            includeBoundary,
                                            EPS);
      }
    }
  }
}

/*!
 * \brief Tests for general 3D triangle-triangle intersection.
 * \return status true iff triangles 1 and 2 intersect.
 *
 * Previous tests have ruled out cases where planes of triangles 1 and 2
 * are parallel or identical, as well as cases where plane 1 does not
 * intersect triangle 2 (and vice versa).  The vertices have been permuted
 * so p1 is across plane 2 from q1 and r1, and p2 is across plane 1
 * from q2 and r2.
 *
 * The core of the Devillers and Guigue's method examines the line l0
 * defined by the intersection of
 * planes 1 and 2.  Assume triangles 1 and 2 both intersect this line, which
 * they must if they intersect each other.  Then, the intersection of triangle
 * 1 with l0 is a segment s1, and the intersection of triangle 2 with l0 is
 * a segment s2.  If s1 and s2 overlap, triangles 1 and 2 intersect.
 *
 * This function implements Equation 1 from Devillers and Guigue (2002), p.8
 * using a hint from p.10 that greatly simplifies the computation.
 *
 * Helper function for T-T intersect.
 */
AXOM_HOST_DEVICE
inline bool intersectTwoPermutedTriangles(const Point3& p1,
                                          const Point3& q1,
                                          const Point3& r1,
                                          const Point3& p2,
                                          const Point3& q2,
                                          const Point3& r2,
                                          bool includeBoundary,
                                          double EPS)
{
  /* Step 5: From step's 1 through 4, we now have two triangles that,
     if intersecting, have a line that intersects segments p1r1, p1q1,
     p2q2, and p2r2.  We check if these two intervals overlap:
   */
  const bool bdr = includeBoundary;

  return isLpeq(Vector3(q1, q2).dot(Triangle3(q1, p2, p1).normal()), 0.0, bdr, EPS) &&
    isLpeq(Vector3(p1, r2).dot(Triangle3(p1, p2, r1).normal()), 0.0, bdr, EPS);
}

/*! @} */

/*! @{ @name 2D triangle-triangle intersection */
/*
 * \brief Tests if 2D Triangles t1 and t2 intersect.
 * \return status true iff t1 intersects with t2, otherwise, false.
 *
 * Note that edge and vertex intersections are reported as hits: an edge or
 * vertex intersecting any part of another triangle causes a return value of
 * true.  This is consistent with the paper.
 */
template <typename T>
bool intersect_tri2D_tri2D(const primal::Triangle<T, 2>& t1,
                           const primal::Triangle<T, 2>& t2,
                           bool includeBoundary,
                           double EPS)
{
  SLIC_CHECK_MSG(!t1.degenerate(), "\n\n WARNING \n\n Triangle " << t1 << " is degenerate");
  SLIC_CHECK_MSG(!t2.degenerate(), "\n\n WARNING \n\n Triangle " << t2 << " is degenerate");

  return TriangleIntersection2D(t1, t2, includeBoundary, EPS);
}

/*!
 * \brief Check for 2D triangle-edge intersection, given p1 close to r2p2.
 * \return status true iff coplanar CCW triangles 1 and 2 intersect.
 */
AXOM_HOST_DEVICE
inline bool checkEdge(const Point2& p1,
                      const Point2& q1,
                      const Point2& r1,
                      const Point2& p2,
                      const Point2& r2,
                      bool includeBoundary,
                      double EPS);

/*!
 * \brief Check for 2D triangle-edge intersection, given p1 close to r2.
 * \return status true iff coplanar CCW triangles 1 and 2 intersect.
 */
AXOM_HOST_DEVICE
inline bool checkVertex(const Point2& p1,
                        const Point2& q1,
                        const Point2& r1,
                        const Point2& p2,
                        const Point2& q2,
                        const Point2& r2,
                        bool includeBoundary,
                        double EPS);

/*!
 * \brief Compute cross product of two 2D vectors as if they were 3D.
 * \return Cross product of A C and B C.
 *
 * This function treats three Point2 values as corners of a 3D triangle with
 * zero Z-coordinate.  Thus we can calculate the cross product of A C with
 * B C using only the k-hat term, since the other terms go to zero.  A
 * positive value indicates CCW orientation.
 *
 * \note The result is equal to twice the signed area of a 2D triangle
 * with vertices (A,B,C) (in CCW order).
 */
AXOM_HOST_DEVICE
inline double twoDcross(const Point2& A, const Point2& B, const Point2& C)
{
  return (((A[0] - C[0]) * (B[1] - C[1]) - (A[1] - C[1]) * (B[0] - C[0])));
}

/*!
 * \brief Check if x, y, and z all have the same sign.
 */
AXOM_HOST_DEVICE
inline bool nonzeroSignMatch(double x, double y, double z, double EPS)
{
  return !(axom::utilities::isNearlyEqual(x, 0., EPS)) &&
    !(axom::utilities::isNearlyEqual(y, 0., EPS)) && !(axom::utilities::isNearlyEqual(z, 0., EPS)) &&
    (0 < x) - (x < 0) == (0 < y) - (y < 0) && (0 < x) - (x < 0) == (0 < z) - (z < 0);
}

/*!
 * \brief Check if two of x, y, and z are near zero.
 */
AXOM_HOST_DEVICE
inline bool twoZeros(double x, double y, double z, double EPS)
{
  return countZeros(x, y, z, EPS) == 2;
}

/*!
 * \brief Check if one of x, y, and z is near zero and the others' signs match.
 */
AXOM_HOST_DEVICE
inline bool oneZeroOthersMatch(double x, double y, double z, double EPS)
{
  namespace util = axom::utilities;
  return countZeros(x, y, z, EPS) == 1 &&
    ((util::isNearlyEqual(x, 0.0, EPS) && isGt(y * z, 0.0, EPS)) ||
     (util::isNearlyEqual(y, 0.0, EPS) && isGt(z * x, 0.0, EPS)) ||
     (util::isNearlyEqual(z, 0.0, EPS) && isGt(x * y, 0.0, EPS)));
}

/*!
 * \brief Count the number of arguments near zero.
 */
AXOM_HOST_DEVICE
inline int countZeros(double x, double y, double z, double EPS)
{
  return (int)axom::utilities::isNearlyEqual(x, 0.0, EPS) +
    (int)axom::utilities::isNearlyEqual(y, 0.0, EPS) +
    (int)axom::utilities::isNearlyEqual(z, 0.0, EPS);
}

/*! @} */

/*!
 * \brief Helper function to find disjoint projections for the AABB-triangle
 * test
 * \param d0 The first value defining the test interval
 * \param d1 The second value defining the test interval
 * \param d2 The third value defining the test interval
 * \param r Radius of projection
 * \return True of the intervals are disjoint, false otherwise
 */
AXOM_HOST_DEVICE bool intervalsDisjoint(double d0, double d1, double d2, double r);
AXOM_HOST_DEVICE bool crossEdgesDisjoint(double d0, double d1, double r);

/*! @{ @name Triangle-bbox intersection */

/*!
 * \brief Determines if a triangle and a bounding box intersect
 *        (but does not find the intersections)
 * \param [in] tri user-supplied triangle (with three vertices).
 * \param [in] bb user-supplied axis aligned bounding box.
 * \return true iff tri intersects with bb, otherwise, false.
 */
template <typename T>
AXOM_HOST_DEVICE bool intersect_tri_bbox(const primal::Triangle<T, 3>& tri,
                                         const primal::BoundingBox<T, 3>& bb)
{
  // Note: Algorithm is derived from the one presented in chapter 5.2.9 of
  //   Real Time Collision Detection book by Christer Ericson
  // based on Akenine-Moller algorithm (Journal of Graphics Tools)
  //
  // It uses the Separating Axis Theorem to look for disjoint projections
  // along various axes associated with Faces and Edges of the AABB and
  // triangle.
  // There are 9 tests for the cross products of edges
  //           3 tests for the AABB face normals
  //           1 test for the triangle face normal
  // We use early termination if we find a separating axis between the shapes

  using PointType = typename BoundingBox<T, 3>::PointType;
  using VectorType = typename BoundingBox<T, 3>::VectorType;

  // Check for early return -- nothing intersects an empty/invalid bounding box
  if(!bb.isValid())
  {
    return false;
  }

  // Extent: vector center to max corner of BB
  const VectorType e = 0.5 * bb.range();

  // Make the AABB center the origin by moving the triangle vertices
  const PointType center = bb.getMin() + e;
  const VectorType v[3] = {tri[0] - center, tri[1] - center, tri[2] - center};

  // Create the edge vectors of the triangle
  const VectorType f[3] = {v[1] - v[0], v[2] - v[1], v[0] - v[2]};

  /* clang-format off */

  // Test cross products of edges between triangle edge vectors f and cube normals (9 tests)
  // -- using separating axis theorem on the cross product of edges of triangle and face normals of AABB
  // Each test involves three cross products, two of which have the same value
  // The commented parameters highlights this symmetry.
#define XEDGE_R( _0, _1, _I )      e[ _0 ] * axom::utilities::abs( f[ _I ][ _1 ]) + e[ _1 ] * axom::utilities::abs(f[ _I ][ _0 ])
#define XEDGE_S( _0, _1, _V, _F ) -v[ _V ][ _0 ] * f[ _F ][ _1 ] + v[ _V ][ _1 ] * f[ _F ][ _0 ]

  if ( crossEdgesDisjoint(/*XEDGE_S(1,2,0,0),*/ XEDGE_S(1,2,1,0),   XEDGE_S(1,2,2,0),   XEDGE_R(1,2,0)) ||
       crossEdgesDisjoint(  XEDGE_S(1,2,0,1),/* XEDGE_S(1,2,1,1),*/ XEDGE_S(1,2,2,1),   XEDGE_R(1,2,1)) ||
       crossEdgesDisjoint(  XEDGE_S(1,2,0,2),   XEDGE_S(1,2,1,2),/* XEDGE_S(1,2,2,2),*/ XEDGE_R(1,2,2)) ||
       crossEdgesDisjoint(/*XEDGE_S(2,0,0,0),*/ XEDGE_S(2,0,1,0),   XEDGE_S(2,0,2,0),   XEDGE_R(0,2,0)) ||
       crossEdgesDisjoint(  XEDGE_S(2,0,0,1),/* XEDGE_S(2,0,1,1),*/ XEDGE_S(2,0,2,1),   XEDGE_R(0,2,1)) ||
       crossEdgesDisjoint(  XEDGE_S(2,0,0,2),   XEDGE_S(2,0,1,2),/* XEDGE_S(2,0,2,2),*/ XEDGE_R(0,2,2)) ||
       crossEdgesDisjoint(/*XEDGE_S(0,1,0,0),*/ XEDGE_S(0,1,1,0),   XEDGE_S(0,1,2,0),   XEDGE_R(0,1,0)) ||
       crossEdgesDisjoint(  XEDGE_S(0,1,0,1),/* XEDGE_S(0,1,1,1),*/ XEDGE_S(0,1,2,1),   XEDGE_R(0,1,1)) ||
       crossEdgesDisjoint(  XEDGE_S(0,1,0,2),   XEDGE_S(0,1,1,2),/* XEDGE_S(0,1,2,2),*/ XEDGE_R(0,1,2)) )
  {
    return false;
  }
  /* clang-format on */

#undef XEDGE_R
#undef XEDEG_S

  /// Test face normals of bounding box (3 tests)
  if(intervalsDisjoint(v[0][0], v[1][0], v[2][0], e[0]) ||
     intervalsDisjoint(v[0][1], v[1][1], v[2][1], e[1]) ||
     intervalsDisjoint(v[0][2], v[1][2], v[2][2], e[2]))
  {
    return false;
  }

  /// Final test -- face normal of triangle's plane
  const VectorType planeNormal = VectorType::cross_product(f[0], f[1]);
  const double planeDist = planeNormal.dot(VectorType(tri[0]));

  const double r = e[0] * axom::utilities::abs(planeNormal[0])  //
    + e[1] * axom::utilities::abs(planeNormal[1])               //
    + e[2] * axom::utilities::abs(planeNormal[2]);
  const double s = planeNormal.dot(VectorType(center)) - planeDist;

  return axom::utilities::abs(s) <= r;
}

// ------------------------------------------------------------------------------

/*!
 * \brief Helper function for Triangle/BoundingBox intersection test
 */
AXOM_HOST_DEVICE inline bool crossEdgesDisjoint(double d0, double d1, double r)
{
  return axom::utilities::max(-axom::utilities::max(d0, d1), axom::utilities::min(d0, d1)) > r;
}

/*! @} */

/*! @{ @name Triangle-ray intersection */

/*!
 * \brief Tests if 3D triangle tri intersects with 3D ray R.
 * \param [in] tri The input triangle
 * \param [in] R The input ray
 * \param [out] t Intersection point of tri and R, w.r.t. parametrization of R
 * \param [out] p Intersection point of tri and R, in **un-normalized**
 *   barycentric coordinates relative to tri.  To normalize, divide each
 *   component by the sum of the components.
 * \note If there is an intersection, the intersection point pt is:
 *                     pt = R.at(t)
 * \note If R is coplanar with tri, this routine will report a miss.  If R's
 *       origin lies within tri, this routine will report a miss.
 * \return status true iff tri intersects with R, otherwise, false.
 *
 * This algorithm is modeled after Woop, Benthin, and Wald (2013).  It
 * transforms the ray onto the unit z vector and applies the same transform
 * to the triangle's vertices.  This transform simplifies the calculation of
 * barycentric coordinates of the positive z-axis's intersection with the
 * triangle.  If any of these coordinates are less than zero, the ray misses.
 * Note that for efficiency, the barycentric coordinates of the intersection
 * are not normalized.  To normalize, divide each coordinate by the sum of
 * the coordinates.
 *
 * If any of the barycentric coordinates are equal to zero, more care is
 * needed to check if the ray hits the edge or misses.
 *
 * Sven Woop, Carsten Benthin, Ingo Wald, "Watertight Ray/Triangle
 * Intersection," Journal of Computer Graphics Techniques (JCGT), vol. 2,
 * no. 1, 65–82, 2013  http://jcgt.org/published/0002/01/05/
 */
template <typename T>
AXOM_HOST_DEVICE bool intersect_tri_ray(const Triangle<T, 3>& tri,
                                        const Ray<T, 3>& R,
                                        T& t,
                                        Point<double, 3>& p)
{
  // Ray origins inside of the triangle are considered a miss.
  // This is a good thing, as pointed out by Matt Larsen in January 2017,
  // because of a common technique in ray chasing through a mesh:
  //  1. Cast the ray until it intersects a triangle (record the hit)
  //  2. Starting from the hit point, re-cast the ray to find the next hit
  // Origin-is-miss avoids finding the already-hit triangle as the next
  // hit.
  //
  // Coplanar rays, evidenced by det == 0, are also considered as a miss.
  // I (Arlie Capps, Jan. 2017) don't understand the motivation at this
  // point, but I'll accept this for now.

  using NumArray = NumericArray<T, 3>;
  const T zero = T();

  //find out dimension where ray direction is maximal
  int kx, ky, kz;

  NumArray r = axom::abs(R.direction().array());

  //z-direction largest
  if((r[2] >= r[0]) && (r[2] >= r[1]))
  {
    kz = 2;
  }
  //y direction largest
  else if((r[1] >= r[0]) && (r[1] >= r[2]))
  {
    kz = 1;
  }
  //x direction largest
  else
  {
    kz = 0;
  }

  //assign other dimensions of the ray
  kx = (kz + 1) % 3;
  ky = (kz + 2) % 3;

  //if necessary swap  ky and kx to preserve triangle winding
  if(R.direction()[kz] < zero)
  {
    axom::utilities::swap(kx, ky);
  }

  //calculate shear constants
  NumericArray<T, 3> shear(1.0f / R.direction()[kz], 3);
  shear[0] *= R.direction()[kx];
  shear[1] *= R.direction()[ky];

  //A,B,C are the triangle vertices offset to Ray's origin
  NumArray A = tri[0].array() - R.origin().array();
  NumArray B = tri[1].array() - R.origin().array();
  NumArray C = tri[2].array() - R.origin().array();

  //shear and scale the vertices
  const T Ax = A[kx] - shear[0] * A[kz];
  const T Ay = A[ky] - shear[1] * A[kz];
  const T Bx = B[kx] - shear[0] * B[kz];
  const T By = B[ky] - shear[1] * B[kz];
  const T Cx = C[kx] - shear[0] * C[kz];
  const T Cy = C[ky] - shear[1] * C[kz];

  //scaled barycentric coordinates
  p[0] = Cx * By - Cy * Bx;
  p[1] = Ax * Cy - Ay * Cx;
  p[2] = Bx * Ay - By * Ax;
  const T& U = p[0];
  const T& V = p[1];
  const T& W = p[2];

  //edge testing
  if((U < zero || V < zero || W < zero) && (U > zero || V > zero || W > zero))
  {
    return false;
  }

  //calculate determinant
  const T det = U + V + W;

  if(det == zero)
  {
    return false;
  }

  //calculate scaled z-coordinates of the vertices and use them to calculate hit
  // distance
  const T Az = shear[2] * A[kz];
  const T Bz = shear[2] * B[kz];
  const T Cz = shear[2] * C[kz];
  t = (U * Az + V * Bz + W * Cz);  // save the parameter of the intersection w.r.t. ray
                                   // R

  //make sure hit is in correct direction
  if(((t < zero) && !(det < zero)) || ((det < zero) && !(t < zero)))
  {
    return false;
  }

  t /= det;
  return true;
}

/*! @} */

/*!
 * \brief Tests if 3D triangle tri intersects with 3D ray S.
 * \param [in] tri The input triangle
 * \param [in] S The input segment
 * \param [out] t Intersection point of tri and S, w.r.t. parametrization of S
 * \param [out] p Intersection point of tri and S, in barycentric coordinates relative to tri
 * \return status true iff tri intersects with R, otherwise, false.
 *
 * This routine uses intersect_tri_ray(), which see.
 */
template <typename T>
AXOM_HOST_DEVICE bool intersect_tri_segment(const Triangle<T, 3>& tri,
                                            const Segment<T, 3>& S,
                                            T& t,
                                            Point<double, 3>& p)
{
  Ray<T, 3> r(S.source(), Vector3(S.source(), S.target()));

  //Ray-triangle intersection does not check endpoints, so we explicitly check here
  if(tri.contains(S.source()))
  {
    t = 0;
    p = tri.physToBarycentric(S.source());
    return true;
  }
  if(tri.contains(S.target()))
  {
    t = 1;
    p = tri.physToBarycentric(S.target());
    return true;
  }

  // The triangle only intersects the segment if it intersects the ray defined by one
  // of its endpoints and the direction defined by its two endpoints.
  // We can parametrize the line as:  r.origin() + t * r.direction()
  // Values of the parameter t between 0 and the length of the segment correspond
  // to points on the segment.
  // Note: if intersect_tri_ray() is true, t must be greater than zero
  if(intersect_tri_ray(tri, r, t, p))
  {
    t = t / static_cast<T>(S.length());
    return t <= 1;
  }
  return false;
}

/*!
 * \brief Determines if a 1D OBB intersects a 1D OBB.
 * \param [in] b1 A 1D OrientedBoundingBox
 * \param [in] b2 A 1D OrientedBoundingBox
 * \return true iff b1 intersects with b2, otherwise, false.
 */
template <typename T>
bool intersect_obb1D_obb1D(const OrientedBoundingBox<T, 1>& b1, const OrientedBoundingBox<T, 1>& b2)
{
  T c1 = b1.getCentroid()[0];
  T c2 = b2.getCentroid()[0];

  T e1 = b1.getExtents()[0];
  T e2 = b2.getExtents()[0];

  if(c1 + e1 > c2 - e2)
  {
    return true;
  }
  if(c2 + e2 > c1 - e1)
  {
    return true;
  }

  return false;
}

/*!
 * \brief Determines if a 2D OBB intersects a 2D OBB.
 * \param [in] b1 A 2D OrientedBoundingBox
 * \param [in] b2 A 2D OrientedBoundingBox
 * \return true iff b1 intersects with b2, otherwise, false.
 */
template <typename T>
bool intersect_obb2D_obb2D(const OrientedBoundingBox<T, 2>& b1, const OrientedBoundingBox<T, 2>& b2)
{
  Vector<T, 2> c1(b1.getCentroid());
  Vector<T, 2> c2(b2.getCentroid());

  Vector<T, 2> e1 = b1.getExtents();
  Vector<T, 2> e2 = b2.getExtents();

  const Vector<T, 2>* u1 = b1.getAxes();
  const Vector<T, 2>* u2 = b2.getAxes();

  Vector<T, 2> d = c2 - c1;

  for(int i = 0; i < 2; ++i)
  {
    if(utilities::abs<T>(d.dot(u1[i])) > e1[i] + utilities::abs<T>((e2[0] * u2[0]).dot(u1[i])) +
         utilities::abs<T>((e2[1] * u2[1]).dot(u1[i])))
    {
      return false;
    }
  }

  for(int i = 0; i < 2; ++i)
  {
    if(utilities::abs<T>(d.dot(u2[i])) > e2[i] + utilities::abs<T>((e1[0] * u1[0]).dot(u2[i])) +
         utilities::abs<T>((e1[1] * u1[1]).dot(u2[i])))
    {
      return false;
    }
  }

  return true;
}

/*!
 * \brief Determines if a 3D OBB intersects a 3D OBB.
 * \param [in] b1 A 3D OrientedBoundingBox
 * \param [in] b2 A 3D OrientedBoundingBox
 * \param [in] EPS error tolerance for intersection
 * \return true iff b1 intersects with b2, otherwise, false.
 */
template <typename T>
bool intersect_obb3D_obb3D(const OrientedBoundingBox<T, 3>& b1,
                           const OrientedBoundingBox<T, 3>& b2,
                           double EPS)
{
  Vector<T, 3> d = Vector<T, 3>(b1.getCentroid()) - Vector<T, 3>(b2.getCentroid());

  Vector<T, 3> e1 = b1.getExtents();
  Vector<T, 3> e2 = b2.getExtents();

  const Vector<T, 3>* u1 = b1.getAxes();
  const Vector<T, 3>* u2 = b2.getAxes();

  // compute r and r^T here:
  Vector<T, 3> r[3];
  Vector<T, 3> rt[3];
  for(int i = 0; i < 3; ++i)
  {
    for(int j = 0; j < 3; ++j)
    {
      r[i][j] = utilities::abs<T>(u1[i].dot(u2[j]));
      rt[i][j] = utilities::abs<T>(u1[j].dot(u2[i]));
    }
  }

  // check for separating planes parallel to faces
  for(int i = 0; i < 3; ++i)
  {
    if(utilities::abs<T>(d.dot(u1[i])) > e1[i] + e2.dot(r[i]) + EPS)
    {
      return false;
    }
  }

  for(int i = 0; i < 3; ++i)
  {
    if(utilities::abs<T>(d.dot(u2[i])) > e2[i] + e1.dot(rt[i]) + EPS)
    {
      return false;
    }
  }

  // check for separating planes with normals parallel to cross product of edges
  for(int i = 0; i < 3; ++i)
  {
    for(int j = 0; j < 3; ++j)
    {
      T left = utilities::abs<T>(d.dot(u1[(i + 2) % 3]) * r[(i + 1) % 3][j] -
                                 (d.dot(u1[(i + 1) % 3]) * r[(i + 2) % 3][j]));
      T right = (e1[(i + 1) % 3] * r[(i + 2) % 3][j] + e1[(i + 2) % 3] * r[(i + 1) % 3][j]);
      right += (e2[(i + 1) % 3] * r[i][(j + 2) % 3] + e2[(i + 2) % 3] * r[i][(j + 1) % 3]);
      if(left > right + EPS)
      {
        return false;
      }
    }
  }

  // didn't find a separating anything
  return true;
}

/*!
 * \brief Determines if a 3D plane intersects a 3D bounding box.
 *        By default (i.e. \a checkOverlaps is false), checks if |s| <= r, 
 *        where \a s is the distance of the bounding box center to the plane,
 *        and \a r is the projected radius of the bounding box along the line
 *        parallel to the plane normal and going through the box center.
 *        If \a checkOverlaps is true, checks if |s| < r,
 *        where the bounding box overlaps both half spaces of the plane.
 * \param [in] p A 3D plane
 * \param [in] bb A 3D bounding box
 * \param [in] checkOverlaps If true, checks if bounding box overlaps both 
 *             halfspaces of the plane. Otherwise, overlap of both halfspaces 
 *             is not guaranteed. Default is false.
 * \param [in] EPS tolerance parameter for determining if \a s
 *             is just within min/max of \a r.
 * \return true iff plane intersects with bounding box, otherwise, false.
 */
template <typename T>
AXOM_HOST_DEVICE bool intersect_plane_bbox(const Plane<T, 3>& p,
                                           const BoundingBox<T, 3>& bb,
                                           bool checkOverlaps = false,
                                           double EPS = 1E-12)
{
  using VectorType = Vector<T, 3>;

  // Check for early return -- planes cannot intersect invalid bounding boxes
  if(!bb.isValid())
  {
    return false;
  }

  const auto centroid = bb.getCentroid();
  const VectorType e = bb.getMax() - centroid;

  const auto& N = p.getNormal();
  const T r = e[0] * utilities::abs<T>(N[0])  //
    + e[1] * utilities::abs<T>(N[1])          //
    + e[2] * utilities::abs<T>(N[2]);

  const T s = p.signedDistance(centroid);

  return checkOverlaps ? isLt(utilities::abs<T>(s), r, EPS) : isLeq(utilities::abs<T>(s), r, EPS);
}

/*!
 * \brief Determines if a plane intersects a segment.
 * \param [in] plane A plane
 * \param [in] seg A segment
 * \param [out] t Intersection point of plane and seg, w.r.t. 
 *   parametrization of seg
 * \param [in] EPS tolerance parameter for determining if 0.0 <= t <= 1.0
 * \return true iff plane intersects with segment, otherwise, false.
 */
template <typename T, int DIM>
AXOM_HOST_DEVICE bool intersect_plane_seg(const Plane<T, DIM>& plane,
                                          const Segment<T, DIM>& seg,
                                          T& t,
                                          double EPS = 1E-12)
{
  using VectorType = Vector<T, DIM>;

  VectorType ab(seg.source(), seg.target());
  VectorType normal = plane.getNormal();

  t = (plane.getOffset() - normal.dot(VectorType(seg.source()))) / (normal.dot(ab));

  if(isGeq(t, 0.0, EPS) && isLeq(t, 1.0, EPS))
  {
    return true;
  }

  return false;
}

/*!
 * \brief Determines if a 2D sphere intersects a bounding box.
 * \param [in] circle A 2D sphere
 * \param [in] bb A 2D bounding box
 * \param [in] EPS Tolerance for determining if intersection nearly occurs
 * \return true iff sphere intersects a bounding box, false otherwise.
 */
template <typename T>
bool intersect_circle_bbox(const Sphere<T, 2>& circle, const BoundingBox<T, 2>& bbox, double EPS = 1E-12)
{
  auto center = circle.getCenter();
  auto radius = circle.getRadius();
  T dx = axom::utilities::clampVal(center[0], bbox.getMin()[0], bbox.getMax()[0]);
  T dy = axom::utilities::clampVal(center[1], bbox.getMin()[1], bbox.getMax()[1]);

  // Find the distance from the center to the closest point on the bounding box,
  //  and check if that distance is less than the radius (plus a tolerance)
  if((center[0] - dx) * (center[0] - dx) + (center[1] - dy) * (center[1] - dy) <=
     radius * radius + EPS)
  {
    return true;
  }

  return false;
}

AXOM_HOST_DEVICE
inline bool intersectOnePermutedTriangle(const Point3& p1,
                                         const Point3& q1,
                                         const Point3& r1,
                                         const Point3& p2,
                                         const Point3& q2,
                                         const Point3& r2,
                                         double dp2,
                                         double dq2,
                                         double dr2,
                                         Vector3& normal,
                                         bool includeBoundary,
                                         double EPS)
{
  /* Step 4: repeat Step 3, except doing it for triangle 2
     instead of triangle 1 */
  if(isGt(dp2, 0.0, EPS))
  {
    if(isGt(dq2, 0.0, EPS))
    {
      return intersectTwoPermutedTriangles(p1, r1, q1, r2, p2, q2, includeBoundary, EPS);
    }
    else if(isGt(dr2, 0.0, EPS))
    {
      return intersectTwoPermutedTriangles(p1, r1, q1, q2, r2, p2, includeBoundary, EPS);
    }
    else
    {
      return intersectTwoPermutedTriangles(p1, q1, r1, p2, q2, r2, includeBoundary, EPS);
    }
  }
  else if(isLt(dp2, 0.0, EPS))
  {
    if(isLt(dq2, 0.0, EPS))
    {
      return intersectTwoPermutedTriangles(p1, q1, r1, r2, p2, q2, includeBoundary, EPS);
    }
    else if(isLt(dr2, 0.0, EPS))
    {
      return intersectTwoPermutedTriangles(p1, q1, r1, q2, r2, p2, includeBoundary, EPS);
    }
    else
    {
      return intersectTwoPermutedTriangles(p1, r1, q1, p2, q2, r2, includeBoundary, EPS);
    }
  }
  else
  {
    if(isLt(dq2, 0.0, EPS))
    {
      if(isGeq(dr2, 0.0, EPS))
      {
        return intersectTwoPermutedTriangles(p1, r1, q1, q2, r2, p2, includeBoundary, EPS);
      }
      else
      {
        return intersectTwoPermutedTriangles(p1, q1, r1, p2, q2, r2, includeBoundary, EPS);
      }
    }
    else if(isGt(dq2, 0.0, EPS))
    {
      if(isGt(dr2, 0.0, EPS))
      {
        return intersectTwoPermutedTriangles(p1, r1, q1, p2, q2, r2, includeBoundary, EPS);
      }
      else
      {
        return intersectTwoPermutedTriangles(p1, q1, r1, q2, r2, p2, includeBoundary, EPS);
      }
    }
    else
    {
      if(isGt(dr2, 0.0, EPS))
      {
        return intersectTwoPermutedTriangles(p1, q1, r1, r2, p2, q2, includeBoundary, EPS);
      }
      else if(isLt(dr2, 0.0, EPS))
      {
        return intersectTwoPermutedTriangles(p1, r1, q1, r2, p2, q2, includeBoundary, EPS);
      }
      else
      {
        return intersectCoplanar3DTriangles(p1, q1, r1, p2, q2, r2, normal, includeBoundary, EPS);
      }
    }
  }
}

AXOM_HOST_DEVICE
inline bool intersectCoplanar3DTriangles(const Point3& p1,
                                         const Point3& q1,
                                         const Point3& r1,
                                         const Point3& p2,
                                         const Point3& q2,
                                         const Point3& r2,
                                         Vector3 normal,
                                         bool includeBoundary,
                                         double EPS)
{
  /* Co-planar triangles are projected onto the axis that maximizes their
     area and the 2d intersection used to check if they intersect.
   */

  //find triangle with maximum area:
  for(int i = 0; i < 3; i++)
  {
    normal[i] = axom::utilities::abs(normal[i]);
  }

  if((isGt(normal[0], normal[2], EPS)) && (isGeq(normal[0], normal[1], EPS)))
  {
    //if x projection area greatest, project on YZ and return 2D checker

    const Triangle2 t1_2da = Triangle2(Point2::make_point(q1[2], q1[1]),
                                       Point2::make_point(p1[2], p1[1]),
                                       Point2::make_point(r1[2], r1[1]));

    const Triangle2 t2_2da = Triangle2(Point2::make_point(q2[2], q2[1]),
                                       Point2::make_point(p2[2], p2[1]),
                                       Point2::make_point(r2[2], r2[1]));

    return TriangleIntersection2D(t1_2da, t2_2da, includeBoundary, EPS);
  }
  else if(isGt(normal[1], normal[2], EPS) && isGeq(normal[1], normal[0], EPS))
  {
    //if y projection area greatest, project on XZ and return 2D checker
    const Triangle2 t1_2da = Triangle2(Point2::make_point(q1[0], q1[2]),
                                       Point2::make_point(p1[0], p1[2]),
                                       Point2::make_point(r1[0], r1[2]));

    const Triangle2 t2_2da = Triangle2(Point2::make_point(q2[0], q2[2]),
                                       Point2::make_point(p2[0], p2[2]),
                                       Point2::make_point(r2[0], r2[2]));

    return TriangleIntersection2D(t1_2da, t2_2da, includeBoundary, EPS);
  }

  //if z projection area greatest, project on XY and return 2D checker
  const Triangle2 t1_2da = Triangle2(Point2::make_point(p1[0], p1[1]),
                                     Point2::make_point(q1[0], q1[1]),
                                     Point2::make_point(r1[0], r1[1]));

  const Triangle2 t2_2da = Triangle2(Point2::make_point(p2[0], p2[1]),
                                     Point2::make_point(q2[0], q2[1]),
                                     Point2::make_point(r2[0], r2[1]));

  return TriangleIntersection2D(t1_2da, t2_2da, includeBoundary, EPS);
}

AXOM_HOST_DEVICE
inline bool TriangleIntersection2D(const Triangle2& t1,
                                   const Triangle2& t2,
                                   bool includeBoundary,
                                   double EPS)
{
  if(isLt(twoDcross(t1[0], t1[1], t1[2]), 0.0, EPS))
  {
    if((isLt(twoDcross(t2[0], t2[1], t2[2]), 0.0, EPS)))
    {
      return intersectPermuted2DTriangles(t1[0], t1[2], t1[1], t2[0], t2[2], t2[1], includeBoundary, EPS);
    }
    else
    {
      return intersectPermuted2DTriangles(t1[0], t1[2], t1[1], t2[0], t2[1], t2[2], includeBoundary, EPS);
    }
  }
  else
  {
    if(isLt(twoDcross(t2[0], t2[1], t2[2]), 0.0, EPS))
    {
      return intersectPermuted2DTriangles(t1[0], t1[1], t1[2], t2[0], t2[2], t2[1], includeBoundary, EPS);
    }
    else
    {
      return intersectPermuted2DTriangles(t1[0], t1[1], t1[2], t2[0], t2[1], t2[2], includeBoundary, EPS);
    }
  }
}

AXOM_HOST_DEVICE
inline bool intersectPermuted2DTriangles(const Point2& p1,
                                         const Point2& q1,
                                         const Point2& r1,
                                         const Point2& p2,
                                         const Point2& q2,
                                         const Point2& r2,
                                         bool includeBoundary,
                                         double EPS)
{
  // Step 2: Orient triangle 2 to be counter clockwise and break the problem
  // into two generic cases (where we test the vertex for intersection or the
  // edges).
  //
  // See paper at https://hal.inria.fr/inria-00072100/document for more details

  if(isGpeq(twoDcross(p2, q2, p1), 0.0, includeBoundary, EPS))
  {
    if(isGpeq(twoDcross(q2, r2, p1), 0.0, includeBoundary, EPS))
    {
      if(isGpeq(twoDcross(r2, p2, p1), 0.0, includeBoundary, EPS))
      {
        return true;
      }
      else
      {
        return checkEdge(p1, q1, r1, p2, r2, includeBoundary, EPS);  //T1 clockwise
      }
    }
    else
    {
      if(isGpeq(twoDcross(r2, p2, p1), 0.0, includeBoundary, EPS))
      {
        //5 region decomposition with p1 in the +-- region
        return checkEdge(p1, q1, r1, r2, q2, includeBoundary, EPS);
      }
      else
      {
        return checkVertex(p1, q1, r1, p2, q2, r2, includeBoundary, EPS);
      }
    }
  }
  else
  {
    if(isGpeq(twoDcross(q2, r2, p1), 0.0, includeBoundary, EPS))
    {
      if(isGpeq(twoDcross(r2, p2, p1), 0.0, includeBoundary, EPS))
      {
        //four region decomposition.  ++- region
        return checkEdge(p1, q1, r1, q2, p2, includeBoundary, EPS);
      }
      else
      {
        return checkVertex(p1, q1, r1, q2, r2, p2, includeBoundary, EPS);
      }
    }
    else
    {
      return checkVertex(p1, q1, r1, r2, p2, q2, includeBoundary, EPS);
    }
  }
}

AXOM_HOST_DEVICE
inline bool checkEdge(const Point2& p1,
                      const Point2& q1,
                      const Point2& r1,
                      const Point2& p2,
                      const Point2& r2,
                      bool includeBoundary,
                      double EPS)
{
  if(isGpeq(twoDcross(r2, p2, q1), 0.0, includeBoundary, EPS))
  {
    if(isGpeq(twoDcross(r2, p1, q1), 0.0, includeBoundary, EPS))
    {
      if(isGpeq(twoDcross(p1, p2, q1), 0.0, includeBoundary, EPS))
      {
        return true;
      }
      else
      {
        if(isGpeq(twoDcross(p1, p2, r1), 0.0, includeBoundary, EPS) &&
           isGpeq(twoDcross(q1, r1, p2), 0.0, includeBoundary, EPS))
        {
          return true;
        }
        else
        {
          return false;
        }
      }
    }
    else
    {
      return false;
    }
  }
  else
  {
    if(isGpeq(twoDcross(r2, p2, r1), 0.0, includeBoundary, EPS) &&
       isGpeq(twoDcross(q1, r1, r2), 0.0, includeBoundary, EPS) &&
       isGpeq(twoDcross(p1, p2, r1), 0.0, includeBoundary, EPS))
    {
      return true;
    }
    else
    {
      return false;
    }
  }
}

AXOM_HOST_DEVICE
inline bool checkVertex(const Point2& p1,
                        const Point2& q1,
                        const Point2& r1,
                        const Point2& p2,
                        const Point2& q2,
                        const Point2& r2,
                        bool includeBoundary,
                        double EPS)
{
  // The tests `isGpeq(twoDcross(...))` are checking the orientation
  // of the triangle defined by its three arguments (CCW vs CW)
  // Note: Comments in this function refer to regions
  // in Figure 8 of the paper: https://hal.inria.fr/inria-00072100/document

  // clang-format off
  if (isGpeq(twoDcross(r2, p2, q1), 0.0, includeBoundary, EPS))                // q1 is in {R_22, R_23, R_24, R_25}
  {
    if (isGpeq(twoDcross(q2, r2, q1), 0.0, includeBoundary, EPS))                  // q1 is in {R_23, R_24} or possibly R_22
    {
      if (isGpeq(twoDcross(p1, p2, q1), 0.0, includeBoundary, EPS))                    // q1 in in R_23  or possible R_22
      {
        if (isLpeq(twoDcross(p1, q2, q1), 0.0, includeBoundary, EPS))
        {
          return true;                                                                       // q1 is in R_23  -- intersect
        }
        else
        {
          return false;                                                                      // q1 in in R_22 -- no intersection
        }
      }
      else                                                                             // q1 is in R_24
      {
        if (isGpeq(twoDcross(p1, p2, r1), 0.0, includeBoundary, EPS) &&
            isLpeq(twoDcross(q1, p2, r1), 0.0, includeBoundary, EPS))
        {
          return true;                                                                     // r1 is in R_23 -- intersect!
        }
        else
        {
          return false;                                                                    // no intersection
        }
      }
    }
    else                                                                           // q1 is in {R_22, R_25}
    {
      if (isLpeq(twoDcross(p1, q2, q1), 0.0, includeBoundary, EPS) &&
          isGpeq(twoDcross(q2, r2, r1), 0.0, includeBoundary, EPS) &&
          isGpeq(twoDcross(q1, r1, q2), 0.0, includeBoundary, EPS))
      {
        return true;
      }
      else
      {
        return false;
      }
    }
  }
  else                                                                         // q1 is in R_21
  {
    if (isGpeq(twoDcross(r2, p2, r1), 0.0, includeBoundary, EPS))
    {
      if (isGpeq(twoDcross(q1, r1, r2), 0.0, includeBoundary, EPS))
      {
        if (isGpeq(twoDcross(r1, p1, p2), 0.0, includeBoundary, EPS))
        {
          return true;
        }
        else
        {
          return false;
        }
      }
      else
      {
        if (isGpeq(twoDcross(q1, r1, q2), 0.0, includeBoundary, EPS) &&
            isGpeq(twoDcross(q2, r2, r1), 0.0, includeBoundary, EPS))
        {
          return true;
        }
        else
        {
          return false;
        }
      }
    }
    else
    {
      return false;
    }
  }
  // clang-format on
}

AXOM_HOST_DEVICE inline bool intervalsDisjoint(double d0, double d1, double d2, double r)
{
  if(d1 < d0)
  {
    axom::utilities::swap(d1, d0);  // d0 < d1
  }
  if(d2 > d1)
  {
    axom::utilities::swap(d2, d1);  // d1 is max(d0,d1,d2)
  }
  else if(d2 < d0)
  {
    axom::utilities::swap(d2, d0);  // d0 is min(d0,d1,d2)
  }
  SLIC_ASSERT(d0 <= d1 && d0 <= d2);
  SLIC_ASSERT(d1 >= d0 && d1 >= d2);

  return d1 < -r || d0 > r;
}

AXOM_SUPPRESS_HD_WARN
template <typename T>
AXOM_HOST_DEVICE bool intersect_plane_tet3d(const Plane<T, 3>& p,
                                            const Tetrahedron<T, 3>& tet,
                                            Polygon<T, 3>& intersection)
{
  intersection.clear();

  T distances[4];
  distances[0] = p.signedDistance(tet[0]);
  distances[1] = p.signedDistance(tet[1]);
  distances[2] = p.signedDistance(tet[2]);
  distances[3] = p.signedDistance(tet[3]);

  constexpr T zero {0};
  int gt[4];
  gt[0] = distances[0] > zero;
  gt[1] = distances[1] > zero;
  gt[2] = distances[2] > zero;
  gt[3] = distances[3] > zero;
  int caseNumber = (gt[3] << 3) | (gt[2] << 2) | (gt[1] << 1) | gt[0];

  // Cases 0, 15 indicate all points were on one side of the zero. That would
  // normally not produce geometry. We may have zeroes and some bigger negatives
  // so try picking a better case.
  if(caseNumber == 0 || caseNumber == 15)
  {
    int lt[4];
    lt[0] = distances[0] < zero;
    lt[1] = distances[1] < zero;
    lt[2] = distances[2] < zero;
    lt[3] = distances[3] < zero;
    caseNumber = (~((lt[3] << 3) | (lt[2] << 2) | (lt[1] << 1) | lt[0])) & 0xf;
  }

  // clang-format off
  // Edge lookup (pairs of point ids)
  const std::uint8_t edges[] = {
   /* 0 */ // None
   /* 1 */ 0, 2, 0, 1, 0, 3,
   /* 2 */ 0, 1, 1, 2, 1, 3,
   /* 3 */ 0, 2, 1, 2, 1, 3, 0, 3,
   /* 4 */ 1, 2, 0, 2, 2, 3,
   /* 5 */ 1, 2, 0, 1, 0, 3, 2, 3,
   /* 6 */ 0, 1, 0, 2, 2, 3, 1, 3,
   /* 7 */ 2, 3, 1, 3, 0, 3,
   /* 8 */ 0, 3, 1, 3, 2, 3,
   /* 9 */ 0, 2, 0, 1, 1, 3, 2, 3,
   /* 10 */ 0, 1, 1, 2, 2, 3, 0, 3,
   /* 11 */ 0, 2, 1, 2, 2, 3,
   /* 12 */ 1, 2, 0, 2, 0, 3, 1, 3,
   /* 13 */ 1, 2, 0, 1, 1, 3,
   /* 14 */ 0, 1, 0, 2, 0, 3,
   /* 15 */ // None
  };
  const std::uint8_t edges_offsets[] = {0, 0, 6, 12, 20, 26, 34, 42, 48, 54, 62, 70, 76, 84, 90, 96, 96};
  // clang-format on

  // Add new vertices to the polygon according to the case.
  const int n = static_cast<int>(edges_offsets[caseNumber + 1] - edges_offsets[caseNumber]) >> 1;
  for(int i = 0; i < n; i++)
  {
    const auto eOffset = static_cast<int>(edges_offsets[caseNumber]) + (i << 1);
    const auto p0 = static_cast<int>(edges[eOffset]);
    const auto p1 = static_cast<int>(edges[eOffset + 1]);

    const double d0 = (distances[p0] < 0.) ? -distances[p0] : distances[p0];
    const double d1 = (distances[p1] < 0.) ? -distances[p1] : distances[p1];

    const auto t = d0 / (d0 + d1);
    intersection.addVertex(Point<T, 3>::lerp(tet[p0], tet[p1], t));
  }

  return caseNumber > 0 && caseNumber < 15;
}

/*! \brief Determines if a line intersects a bilinear patch.
 * \param [in] line The line to intersect with the bilinear patch.
 * \param [in] p0 The first corner of the bilinear patch.
 * \param [in] p1 The second corner in ccw order.
 * \param [in] p2 The third corner.
 * \param [in] p3 The fourth corner.
 * \param [out] t Array to append the t parameters of intersections wrt the ray.
 * \param [out] u Array to append the u parameters of intersections wrt the patch.
 * \param [out] v Array to append the v parameters of intersections wrt the patch.
 * \param [in] EPS The parameter space tolerance for intersection.
 * \param [in] isRay If true, only return intersections with t >= 0.
 *
 * Implements GARP algorithm from Chapter 8 of Ray Tracing Gems (2019)
 *
 * \note A bilinear patch is parameterized in [0 - EPS, 1 + EPS]^2 
 * 
 * \note There can be either 0, 1, 2 discrete intersections, or a continuous segment
 *  of intersections, in which case the method returns the ceneter of this segment as one point.
 * 
 * \note Always returns false if the line is coplanar to a planar polygon
 * 
 * \return true if the line intersects the bilinear patch, otherwise false.
 */
AXOM_HOST_DEVICE
inline bool intersect_line_bilinear_patch(const Line<double, 3>& line,
                                          const Point3& p00,
                                          const Point3& p10,
                                          const Point3& p11,
                                          const Point3& p01,
                                          axom::StaticArray<double, 2>& t,
                                          axom::StaticArray<double, 2>& u,
                                          axom::StaticArray<double, 2>& v,
                                          double EPS = 1e-8,
                                          bool isRay = false)
{
  Vector3 e10(p00, p10);
  Vector3 e11(p10, p11);
  Vector3 e00(p00, p01);

  Vector3 qn = Vector3::cross_product(e10, Vector3(p11, p01));

  Vector3 q00(line.origin(), p00), q10(line.origin(), p10);

  // Solve a quadratic to find the parameters u0 of the B(u0, v) isocurves
  //  that are closest to the line
  double au = Vector3::scalar_triple_product(q00, line.direction(), e00);
  double cu = Vector3::dot_product(qn, line.direction());
  double bu = Vector3::scalar_triple_product(q10, line.direction(), e11) - au - cu;

  // Rescale the coefficients to avoid (some) numerical issues
  double su =
    axom::utilities::max(axom::utilities::abs(au),
                         axom::utilities::max(axom::utilities::abs(bu), axom::utilities::abs(cu)));
  au /= su;
  bu /= su;
  cu /= su;

  // Decide what to do based on the coefficients
  if(!axom::utilities::isNearlyEqual(su, 0.0, primal::PRIMAL_TINY))
  {
    double u1, u2;

    // Decide what to do based on the discriminant
    double det = bu * bu - 4 * au * cu;
    if(det < 0)  // No solutions
    {
      return false;
    }
    else if(axom::utilities::isNearlyEqual(cu, 0.0, EPS))  // Quadratic is a line
    {
      u1 = -au / bu;
      u2 = -1;
    }
    else if(axom::utilities::isNearlyEqual(det, 0.0, EPS))  // One repeated solution
    {
      u1 = -bu / (2 * cu);
      u2 = -1;
    }
    else
    {
      u1 = 0.5 * (-bu - std::copysign(std::sqrt(det), bu));
      u2 = au / u1;
      u1 /= cu;
    }

    // Find the point on the isocurve that is closest to the ray
    for(auto u0 : {u1, u2})
    {
      if(u0 < -EPS || u0 > 1.0 + EPS)
      {
        continue;
      }

      Vector3 pa = (1 - u0) * q00 + u0 * q10;
      Vector3 pb = (1 - u0) * e00 + u0 * e11;  // actually stores pb - pa
      Vector3 n = Vector3::cross_product(line.direction(), pb);
      det = Vector3::dot_product(n, n);

      // Need a separate tolerance for this for the case of small patches
      if(!axom::utilities::isNearlyEqual(det, 0.0, primal::PRIMAL_TINY))
      {
        n = Vector3::cross_product(n, pa);
        double t0 = Vector3::dot_product(n, pb) / det;
        double v0 = Vector3::dot_product(n, line.direction()) / det;
        if(-EPS <= v0 && v0 <= 1.0 + EPS)
        {
          if(t0 >= -EPS || !isRay)
          {
            t.push_back(t0);
            u.push_back(u0);
            v.push_back(v0);
          }
        }
      }
      else  // Ray is parallel to the line segment pa + v * (pb - pa)
      {
        // Determine if the line is colinear to the segment
        double cross = Vector3::cross_product(pa, line.direction()).norm();
        if(axom::utilities::isNearlyEqual(cross, 0.0, EPS))
        {
          // Parameters of intersection are non-unique,
          //  so take the center of the segment the intersection
          //  (this avoids any inclusion issues at the boundary)
          const double t1 = Vector3::dot_product(pa, line.direction());
          const double t2 = Vector3::dot_product(pa + pb, line.direction());

          if(!isRay || axom::utilities::min(t1, t2) > 0.0)
          {
            // Always an intersection in this case
            t.push_back(0.5 * (t1 + t2));
            v.push_back(0.5);
            u.push_back(u0);

            return true;
          }
          else if(t1 * t2 <= 0)
          {
            // Means the origin is inside the segment

            // Switch based on orientation of segment
            if(t1 == t2)
            {
              t.push_back(0.5 * t1);
              v.push_back(0.5);
            }
            else if(t1 < t2)
            {
              t.push_back(0.5 * t2);
              v.push_back((t1 - 0.5 * t2) / (t1 - t2));
            }
            else
            {
              t.push_back(0.5 * t1);
              v.push_back(-0.5 * t1 / (t2 - t1));
            }
            u.push_back(u0);

            return true;
          }
          else
          {
            // No intersection
            return false;
          }
        }
      }
    }
  }
  else
  {
    // Switch to finding B(u, v0) isocurves instead
    Vector3 e01(p01, p11);
    Vector3 qm = Vector3::cross_product(e00, -e11);
    Vector3 q01(line.origin(), p01);

    // Find the analogous coefficients for B(u, v0) isocurves
    double av = Vector3::scalar_triple_product(q00, line.direction(), e10);
    double cv = Vector3::dot_product(qm, line.direction());
    double bv = Vector3::scalar_triple_product(q01, line.direction(), e01) - av - cv;

    double sv =
      axom::utilities::max(axom::utilities::abs(av),
                           axom::utilities::max(axom::utilities::abs(bv), axom::utilities::abs(cv)));
    av /= sv;
    bv /= sv;
    cv /= sv;

    // If these are also all zero, then the ray is coplanar to a polygonal patch
    if(axom::utilities::isNearlyEqual(sv, 0.0, primal::PRIMAL_TINY))
    {
      // This case indicates the ray is tangent to the surface,
      //  which we don't count as an intersection
      return false;
    }

    double v1, v2;

    // Decide what to do based on the discriminant
    double det = bv * bv - 4 * av * cv;
    if(det < 0)  // No solutions
    {
      return false;
    }
    else if(axom::utilities::isNearlyEqual(cv, 0.0, EPS))  // Quadratic is a line
    {
      v1 = -av / bv;
      v2 = -1;
    }
    else if(axom::utilities::isNearlyEqual(det, 0.0, EPS))  // One repeated solution
    {
      v1 = -bv / (2 * cv);
      v2 = -1;
    }
    else
    {
      v1 = 0.5 * (-bv - std::copysign(std::sqrt(det), bv));
      v2 = av / v1;
      v1 /= cv;
    }

    // Find the point on the isocurve that is closest to the ray
    for(auto v0 : {v1, v2})
    {
      if(v0 < -EPS || v0 > 1.0 + EPS)
      {
        continue;
      }

      Vector3 pa = (1.0 - v0) * q00 + v0 * q01;
      Vector3 pb = (1.0 - v0) * e10 + v0 * e01;  // actually stores pb - pa
      Vector3 n = Vector3::cross_product(line.direction(), pb);
      det = Vector3::dot_product(n, n);

      // Need a separate tolerance for this for the case of small patches
      if(!axom::utilities::isNearlyEqual(det, 0.0, primal::PRIMAL_TINY))
      {
        n = Vector3::cross_product(n, pa);
        double t0 = Vector3::dot_product(n, pb) / det;
        double u0 = Vector3::dot_product(n, line.direction()) / det;

        if(-EPS <= u0 && u0 <= 1.0 + EPS)
        {
          if(t0 >= -EPS || !isRay)
          {
            t.push_back(t0);
            u.push_back(u0);
            v.push_back(v0);
          }
        }
      }
      else  // Ray is parallel to the line segment pa + u * (pb - pa)
      {
        // Determine if the line is colinear to the segment
        double cross = Vector3::cross_product(pa, line.direction()).norm();
        if(axom::utilities::isNearlyEqual(cross, 0.0, EPS))
        {
          // Parameters of intersection are non-unique,
          //  so take the center of the segment the intersection
          //  (this avoids any inclusion issues at the boundary)
          const double t1 = Vector3::dot_product(pa, line.direction());
          const double t2 = Vector3::dot_product(pa + pb, line.direction());

          if(!isRay || axom::utilities::min(t1, t2) > 0.0)
          {
            // Always an intersection in this case
            t.push_back(0.5 * (t1 + t2));
            u.push_back(0.5);
            v.push_back(v0);

            return true;
          }
          else if(t1 * t2 <= 0)
          {
            // Means the origin is inside the segment

            // Switch based on orientation of segment
            if(t1 == t2)
            {
              t.push_back(0.5 * t1);
              u.push_back(0.5);
            }
            else if(t1 < t2)
            {
              t.push_back(0.5 * t2);
              u.push_back((t1 - 0.5 * t2) / (t1 - t2));
            }
            else
            {
              t.push_back(0.5 * t1);
              u.push_back(-0.5 * t1 / (t2 - t1));
            }
            v.push_back(v0);

            return true;
          }
          else
          {
            // No intersection
            return false;
          }
        }
      }
    }
  }

  return !t.empty();
}

/*!
 * \brief Examine candidate (t,u,v) values and select the ones that are not
 *        duplicates, storing them in the supplied output arrays.
 *
 * \param [in] tc Candidate t values of intersection points.
 * \param [in] uc Candidate u values of intersection points.
 * \param [in] vc Candidate v values of intersection points.
 * \param [out] t Selected t values of intersection points.
 * \param [out] u Selected u values of intersection points.
 * \param [out] v Selected v values of intersection points.
 * \param [in] EPS The tolerance for intersection (for parameter distances).
 * \param [in] isHalfOpen True if the patch is parameterized in [0,1)^2.
 *
 * \note Moved from intersect for Ray/BezierPatch and templated it so it
 *       supports different candidate and output array types.
 */
template <typename CandidateArrayType, typename ArrayType>
bool select_candidates(const CandidateArrayType& tc,
                       const CandidateArrayType& uc,
                       const CandidateArrayType& vc,
                       ArrayType& t,
                       ArrayType& u,
                       ArrayType& v,
                       double EPS = 1e-8,
                       bool isHalfOpen = false)
{
  using T = typename ArrayType::value_type;

  // Remove duplicates from the (u, v) intersection points
  //  (Note it's not possible for (u_1, v_1) == (u_2, v_2) and t_1 != t_2)
  const double sq_EPS = EPS * EPS;

  // The number of reported intersection points will be small,
  //  so we don't need to fully sort the list
  SLIC_WARNING_IF(tc.size() > 10,
                  "Large number of intersections detected, eliminating "
                  "duplicates may be slow");

  for(int i = 0; i < tc.size(); ++i)
  {
    // Also remove any intersections on the half-interval boundaries
    if(isHalfOpen && (uc[i] >= 1.0 - EPS || vc[i] >= 1.0 - EPS))
    {
      continue;
    }

    Point<T, 2> uv({uc[i], vc[i]});

    bool foundDuplicate = false;
    for(int j = i + 1; !foundDuplicate && j < tc.size(); ++j)
    {
      if(squared_distance(uv, Point<T, 2>({uc[j], vc[j]})) < sq_EPS)
      {
        foundDuplicate = true;
      }
    }

    if(!foundDuplicate)
    {
      t.push_back(tc[i]);
      u.push_back(uc[i]);
      v.push_back(vc[i]);
    }
  }

  return !t.empty();
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif  // AXOM_PRIMAL_INTERSECT_IMPL_HPP_
