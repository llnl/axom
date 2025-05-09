// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file closest_point.hpp
 *
 * \brief Consists of a set of methods that compute the closest point on a
 *  geometric primitive B from another geometric primitive A.
 *
 */

#ifndef AXOM_PRIMAL_CLOSEST_POINT_HPP_
#define AXOM_PRIMAL_CLOSEST_POINT_HPP_

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Sphere.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"
#include "axom/primal/operators/detail/intersect_impl.hpp"
namespace axom
{
namespace primal
{
/*!
 * \brief Computes the closest point from a point, P, to a given segment.
 *
 * \param [in] P the query point
 * \param [in] seg user-supplied segment
 * \param [out] loc location along the line segment
 * \param [in] EPS fuzz factor for equality comparisons
 * \return cp the closest point from a point P and a segment
 *
 * \note loc \f$ \in [0, 1] \f% represents the fraction of the way from A to B,
 *  with 0 corresponding to A and 1 corresponding to B.
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE inline Point<T, NDIMS> closest_point(const Point<T, NDIMS>& P,
                                                      const Segment<T, NDIMS>& seg,
                                                      T* loc,
                                                      double EPS = PRIMAL_TINY)
{
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;

  using detail::isGeq;
  using detail::isLeq;

  constexpr T ZERO {0.};
  constexpr T ONE {1.};

  const PointType& A = seg[0];
  const PointType& B = seg[1];

  const VectorType AB(A, B);

  // Compute length of the projection of AP onto AB
  T t = VectorType(A, P).dot(AB);

  if(isLeq(t, ZERO, EPS))
  {
    if(loc)
    {
      *loc = ZERO;
    }

    return A;
  }
  else
  {
    const T squaredNormAB = AB.squared_norm();

    if(isGeq(t, squaredNormAB, EPS))
    {
      if(loc)
      {
        *loc = ONE;
      }

      return B;
    }
    else if(utilities::isNearlyEqual(squaredNormAB, ZERO, EPS))
    {
      // Segment is degenerate (A and B are collocated),
      // so we can pick either end point. We pick A.
      if(loc)
      {
        *loc = ZERO;
      }

      return A;
    }
    else
    {
      // Normalize t
      t /= squaredNormAB;

      if(loc)
      {
        *loc = t;
      }

      return A + AB * t;
    }
  }
}

/*!
 * \brief Computes the closest point from a point, P, to a given segment.
 *
 * \param [in] P the query point
 * \param [in] seg user-supplied segment
 * \param [in] EPS fuzz factor for equality comparisons
 * \return cp the closest point from a point P and a segment
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE inline Point<T, NDIMS> closest_point(const Point<T, NDIMS>& P,
                                                      const Segment<T, NDIMS>& seg,
                                                      double EPS = PRIMAL_TINY)
{
  T* loc = nullptr;
  return closest_point(P, seg, loc, EPS);
}

/*!
 * \brief Computes the closest point from a point, P, to a given triangle.
 *
 * \param [in] P the query point
 * \param [in] tri user-supplied triangle.
 * \param [out] loc int pointer to store location of closest point (optional).
 * \return cp the closest point from a point P and a triangle.
 *
 * \note If the optional int pointer is supplied for `loc`, the method returns
 *  the location of the closest point, which is illustrated in the schematic
 *  diagram below and encoded as follows:
 * <ul>
 *  <li> loc \f$ \in [0,2] \f$, loc corresponds to the triangle node index </li>
 *  <li> loc \f$ \in [-3,-1] \f$, abs(loc) corresponds to an edge </li>
 *  <li> loc >= 3, loc is on a triangle face </li>
 * </ul>
 *
 * \verbatim
 *
 *            2
 *           /\
 *    (-3)--/  \--(-2)
 *         /    \
 *        /_ _ _ \
 *       0   |    1
 *           |
 *         (-1)
 *
 * \endverbatim
 *
 * \pre NDIMS==2 || NDIMS==3
 *
 * \note Implementation is based on "Real Time Collision Detection,
 *  Chapter 5.1.5 Closest Point on Triangle to Point".
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE inline Point<T, NDIMS> closest_point(const Point<T, NDIMS>& P,
                                                      const Triangle<T, NDIMS>& tri,
                                                      int* loc = nullptr,
                                                      double EPS = PRIMAL_TINY)
{
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;

  using detail::isGeq;
  using detail::isLeq;

  const PointType& A = tri[0];
  const PointType& B = tri[1];
  const PointType& C = tri[2];

  // Check if P in vertex region outside A
  const VectorType ab(A, B);
  const VectorType ac(A, C);
  const VectorType ap(A, P);
  const T d1 = VectorType::dot_product(ab, ap);
  const T d2 = VectorType::dot_product(ac, ap);
  if(isLeq(d1, T(0), EPS) && isLeq(d2, T(0), EPS))
  {
    // A is the closest point
    if(loc != nullptr)
    {
      *loc = 0;
    }

    return A;
  }

  //----------------------------------------------------------------------------
  // Check if P in vertex region outside B
  const VectorType bp(B, P);
  const T d3 = VectorType::dot_product(ab, bp);
  const T d4 = VectorType::dot_product(ac, bp);
  if(isGeq(d3, T(0), EPS) && isLeq(d4, d3, EPS))
  {
    // B is the closest point
    if(loc != nullptr)
    {
      *loc = 1;
    }

    return B;
  }

  //----------------------------------------------------------------------------
  // Check if P in edge region of AB
  const T vc = d1 * d4 - d3 * d2;
  if(isLeq(vc, T(0), EPS) && isGeq(d1, T(0), EPS) && isLeq(d3, T(0), EPS) &&
     !utilities::isNearlyEqual(d1, d3, EPS))  // Additional check for degenerate triangles
  {
    const T v = d1 / (d1 - d3);
    const VectorType v_ab = ab * v;

    if(loc != nullptr)
    {
      *loc = -1;
    }

    return A + v_ab;
  }

  //----------------------------------------------------------------------------
  // Check if P in vertex region outside C
  const VectorType cp(C, P);
  const T d5 = VectorType::dot_product(ab, cp);
  const T d6 = VectorType::dot_product(ac, cp);
  if(isGeq(d6, T(0), EPS) && isLeq(d5, d6, EPS))
  {
    // C is the closest point
    if(loc != nullptr)
    {
      *loc = 2;
    }

    return C;
  }

  //----------------------------------------------------------------------------
  // Check if P in edge region of AC
  const T vb = d5 * d2 - d1 * d6;
  if(isLeq(vb, T(0), EPS) && isGeq(d2, T(0), EPS) && isLeq(d6, T(0), EPS))
  {
    const T w = d2 / (d2 - d6);
    const VectorType w_ac = ac * w;

    if(loc != nullptr)
    {
      *loc = -3;
    }

    return A + w_ac;
  }

  //----------------------------------------------------------------------------
  // Check if P in edge region of BC
  T va = d3 * d6 - d5 * d4;
  if(isLeq(va, T(0), EPS) && isGeq(d4 - d3, T(0), EPS) && isGeq(d5 - d6, T(0), EPS))
  {
    const T w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    const VectorType bc(B, C);
    const VectorType w_bc = bc * w;

    if(loc != nullptr)
    {
      *loc = -2;
    }

    return B + w_bc;
  }

  //----------------------------------------------------------------------------
  // P is inside face region
  const T denom = T(1) / (va + vb + vc);
  const T v = vb * denom;
  const T w = vc * denom;
  const VectorType N = (ab * v) + (ac * w);

  if(loc != nullptr)
  {
    *loc = Triangle<T, NDIMS>::NUM_TRI_VERTS;
  }

  return A + N;
}

/*!
 * \brief Computes the closest point from a point to a given OBB.
 *
 * \param [in] pt the query pt.
 * \param [in] obb user-supplied oriented bounding box.
 * \return cp the closest point from a point pt and an OBB.
 */
template <typename T, int NDIMS>
inline Point<T, NDIMS> closest_point(const Point<T, NDIMS>& pt,
                                     const OrientedBoundingBox<T, NDIMS>& obb)
{
  Vector<T, NDIMS> e = obb.getExtents();
  const Vector<T, NDIMS>* u = obb.getAxes();

  Vector<T, NDIMS> pt_l(obb.toLocal(pt));
  Vector<T, NDIMS> res(obb.getCentroid());

  for(int i = 0; i < NDIMS; i++)
  {
    // since the local coordinates are individually constrained, we can simply
    // choose the "best" local coordinate in each axis direction
    if(pt_l[i] <= e[i] && pt_l[i] >= -e[i])
    {
      res += pt_l[i] * u[i];
    }
    else if(pt_l[i] > e[i])
    {
      res += e[i] * u[i];
    }
    else
    {
      res -= e[i] * u[i];
    }
  }

  return Point<T, NDIMS>(res.array());
}

/*!
 * \brief Computes the closest point from a point, P, to a sphere.
 *
 * \param [in] P the query point
 * \param [in] sphere user-supplied sphere
 * \return cp the closest point on \a sphere to point \a P
 *
 * \note The closest point is uniquely defined everywhere except at the sphere's center.
 * We handle that case by returning the point on the sphere along an arbitrary direction
 * (specifically, the direction determined by \a primal::Vector::unitVector() )
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE inline Point<T, NDIMS> closest_point(const Point<T, NDIMS>& P,
                                                      const Sphere<T, NDIMS>& sphere)
{
  using VectorType = Vector<T, NDIMS>;

  const auto v = VectorType(sphere.getCenter(), P).unitVector();
  return sphere.getCenter() + sphere.getRadius() * v;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_CLOSEST_POINT_HPP_
