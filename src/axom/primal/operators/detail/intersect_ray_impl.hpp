// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_INTERSECT_RAY_IMPL_HPP_
#define AXOM_PRIMAL_INTERSECT_RAY_IMPL_HPP_

// numerics includes
#include "axom/core/numerics/floating_point_limits.hpp"

// primal includes
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Ray.hpp"
#include "axom/primal/geometry/Line.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"

#include <type_traits>  // for std::is_floating_point< T >()

namespace axom
{
namespace primal
{
namespace detail
{
/*!
 * \brief Computes the intersection of the given ray \a R with the segment \a S.
 * 
 * When there is a valid intersection (within tolerance \a EPS), 
 * \a ray_param returns the parametric coordinate of the intersection point along \a R
 * and \a seg_param returns the parametric coordinate of the intersection point along \a S.
 * If the intersection point is nonunique, \a seg_param and \a ray_param return only
 * a single point of intersection at the center of the region.
 * \return status true iff R intersects with S, otherwise, false.
 */
template <typename T>
AXOM_HOST_DEVICE inline bool intersect_ray(const primal::Ray<T, 2>& R,
                                           const primal::Segment<T, 2>& S,
                                           T& ray_param,
                                           T& seg_param,
                                           const double EPS)
{
  AXOM_STATIC_ASSERT(std::is_floating_point<T>::value);

  // Find the parameterized direction of the segment.
  const auto seg_dir = primal::Vector<T, 2>(S.source(), S.target());
  const auto& ray_dir = R.direction();

  // Define lower/upper threshold
  const double tlow = 0.0 - EPS;
  const double thigh = 1.0 + EPS;

  // Equating R(ray_param)=S(seg_param) yields a system of two equations and
  // two unknowns, namely, t0 and t1. We can solve this system directly using Cramer's Rule.
  const double denom = numerics::determinant(ray_dir[0], -seg_dir[0], ray_dir[1], -seg_dir[1]);

  // If denom is (nearly) zero (within a numerical tolerance), the ray and segment are parallel
  if(axom::utilities::isNearlyEqual(denom, 0.0, EPS))
  {
    // Alternate logic for nearly parallel lines
    const primal::Vector<T, 2> col1(R.origin(), S.source());
    double orientation1 = numerics::determinant(col1[0], ray_dir[0], col1[1], ray_dir[1]);

    const primal::Vector<T, 2> col2(R.origin(), S.target());
    double orientation2 = numerics::determinant(col2[0], ray_dir[0], col2[1], ray_dir[1]);

    if(orientation1 * orientation2 > 0.0)
    {
      // The two segment endpoints are on the same side of the ray
      return false;
    }
    else if(axom::utilities::isNearlyEqual(orientation1, 0.0, primal::PRIMAL_TINY) &&
            axom::utilities::isNearlyEqual(orientation2, 0.0, primal::PRIMAL_TINY))
    {
      // Both endpoints are on the ray's line (but maybe not the ray itself)

      // Check orientation of segment relative to ray
      const double t1 = col1[0] * ray_dir[0] + col1[1] * ray_dir[1];
      const double t2 = t1 + seg_dir.dot(ray_dir);

      if(std::min(t1, t2) > 0.0)
      {
        // The origin is outside the segment,
        //  but the ray intersects
        ray_param = 0.5 * (t1 + t2);
        seg_param = 0.5;
      }
      else if(t1 * t2 <= 0)
      {
        // Means the origin is inside the segment

        // Switch based on orientation of segment
        if(t1 == t2)
        {
          ray_param = 0.5 * t1;
          seg_param = 0.5;
        }
        else if(t1 < t2)
        {
          ray_param = 0.5 * t2;
          seg_param = (t1 - 0.5 * t2) / (t1 - t2);
        }
        else
        {
          ray_param = 0.5 * t1;
          seg_param = -0.5 * t1 / (t2 - t1);
        }
      }
      else
      {
        // No intersection
        return false;
      }

      return true;
    }
    else if(axom::utilities::isNearlyEqual(orientation1, 0.0, primal::PRIMAL_TINY))
    {
      // One of the endpoints is on the ray's line (but maybe not the ray itself)
      ray_param = col1[0] * ray_dir[0] + col1[1] * ray_dir[1];
      seg_param = 0.0;

      return true;
    }
    else if(axom::utilities::isNearlyEqual(orientation2, 0.0, primal::PRIMAL_TINY))
    {
      // One of the endpoints is on the ray's line (but maybe not the ray itself)
      ray_param = col2[0] * ray_dir[0] + col2[1] * ray_dir[1];
      seg_param = 1.0;

      return true;
    }

    // If the two endpoints *are* on opposite sides of the ray, then Cramer's rule
    //  will be stable enough to work, since the intersection won't be at infinity.
  }

  // Solve for the ray_param and seg_param directly using Cramer's rule
  const auto sol = S.source().array() - R.origin().array();

  // Note: ray_param is an OUT parameter of this function
  ray_param = numerics::determinant(sol[0], -seg_dir[0], sol[1], -seg_dir[1]) / denom;

  // Note: seg_param is an OUT parameter of this function
  seg_param = numerics::determinant(ray_dir[0], sol[0], ray_dir[1], sol[1]) / denom;

  // Necessary and sufficient criteria for an intersection between
  // ray, R(t0),  and a finite segment S(t1) are:
  // 1. ray_param >= tlow w.r.t. the ray R(ray_param).
  // 2. tlow >= seg_param >= thigh w.r.t. the segment S(seg_param).
  return ((ray_param >= tlow) && (seg_param >= tlow) && (seg_param <= thigh));
}

/*!
 * \brief Helper routine for ray / AABB intersection test.
 *
 * \param [in] x0 coordinate component of the ray origin
 * \param [in] n normal component of the ray direction
 * \param [in] bbmin the AABB min coordinate along a direction
 * \param [in] bbmax the AABB max coordinate along a direction
 *
 * \param [in,out] tmin coordinate of closest intersection point along ray
 * \param [in,out] tmax coordinate of farthest intersection point along ray
 *
 * \param [in] EPS user-supplied tolerance.
 *
 * \return status true if the ray intersects, otherwise false.
 *
 * \note This routine is called by the intersect ray/AABB methods for each spatial dimension.
 */
template <typename T>
AXOM_HOST_DEVICE inline bool intersect_ray_bbox_test(const T& x0,
                                                     const T& n,
                                                     const T& bbmin,
                                                     const T& bbmax,
                                                     T& tmin,
                                                     T& tmax,
                                                     T EPS)
{
  AXOM_STATIC_ASSERT(std::is_floating_point<T>::value);

  constexpr T ZERO = static_cast<T>(0);

  bool status = true;

  if(axom::utilities::isNearlyEqual(n, ZERO, EPS))
  {
    status = (((x0 < bbmin) || (x0 > bbmax)) ? false : true);
  }
  else
  {
    const T invn = static_cast<T>(1.0) / n;
    T t1 = (bbmin - x0) * invn;
    T t2 = (bbmax - x0) * invn;

    if(t1 > t2)
    {
      axom::utilities::swap(t1, t2);
    }

    tmin = axom::utilities::max(tmin, t1);
    tmax = axom::utilities::min(tmax, t2);

    status = ((tmin > tmax) ? false : true);
  }

  return status;
}

/*!
 * \brief Checks if a specified ray intersects with the supplied bounding box.
 *
 * \param [in] x0 x-coordinate of the ray source point
 * \param [in] nx x-component of the ray normal
 * \param [in] y0 y-coordinate of the ray source point
 * \param [in] ny y-component of the ray normal
 * \param [in] xmin x-coordinate of the lower bounding box corner
 * \param [in] xmax x-coordinate of the upper bounding box corner
 * \param [in] ymin y-coordinate of the lower bounding box corner
 * \param [in] ymax y-coordinate of the upper bounding box corner
 * \param [out] t length of intersection point from the ray source point.
 * \param [in] EPS optional tolerance
 *
 * \return status true if the ray intersects the bounding box, otherwise, false.
 */
template <typename T>
AXOM_HOST_DEVICE inline bool intersect_ray(const T& x0,
                                           const T& nx,
                                           const T& y0,
                                           const T& ny,
                                           const T& xmin,
                                           const T& xmax,
                                           const T& ymin,
                                           const T& ymax,
                                           T& t,
                                           T EPS = numerics::floating_point_limits<T>::epsilon())
{
  AXOM_STATIC_ASSERT(std::is_floating_point<T>::value);

  t = 0.0f;
  T tmax = axom::numerics::floating_point_limits<T>::max();

  bool status = true;
  status = status && intersect_ray_bbox_test(x0, nx, xmin, xmax, t, tmax, EPS);
  status = status && intersect_ray_bbox_test(y0, ny, ymin, ymax, t, tmax, EPS);

  return status;
}

/*!
 * \brief Checks if a specified ray intersects with the supplied bounding box.
 *
 * \param [in] x0 x-coordinate of the ray source point
 * \param [in] nx x-component of the ray normal
 * \param [in] y0 y-coordinate of the ray source point
 * \param [in] ny y-component of the ray normal
 * \param [in] z0 z-coordinate of the ray source point
 * \param [in] nz z-component of the ray normal
 * \param [in] xmin x-coordinate of the lower bounding box corner
 * \param [in] xmax x-coordinate of the upper bounding box corner
 * \param [in] ymin y-coordinate of the lower bounding box corner
 * \param [in] ymax y-coordinate of the upper bounding box corner
 * \param [in] zmin z-coordinate of the lower bounding box corner
 * \param [in] zmax z-coordinate of the upper bounding box corner
 * \param [out] t length of intersection point from the ray source point.
 * \param [in] EPS optional tolerance
 *
 * \return status true if the ray intersects the bounding box, otherwise, false.
 */
template <typename T>
AXOM_HOST_DEVICE inline bool intersect_ray(const T& x0,
                                           const T& nx,
                                           const T& y0,
                                           const T& ny,
                                           const T& z0,
                                           const T& nz,
                                           const T& xmin,
                                           const T& xmax,
                                           const T& ymin,
                                           const T& ymax,
                                           const T& zmin,
                                           const T& zmax,
                                           T& t,
                                           T EPS = numerics::floating_point_limits<T>::epsilon())
{
  AXOM_STATIC_ASSERT(std::is_floating_point<T>::value);

  t = 0.0f;
  T tmax = axom::numerics::floating_point_limits<T>::max();

  bool status = true;
  status = status && intersect_ray_bbox_test(x0, nx, xmin, xmax, t, tmax, EPS);
  status = status && intersect_ray_bbox_test(y0, ny, ymin, ymax, t, tmax, EPS);
  status = status && intersect_ray_bbox_test(z0, nz, zmin, zmax, t, tmax, EPS);

  return status;
}

/*!
 * \brief Computes the intersection of the given ray, R, with the Box, bb
 *
 * \param [in] R the specified ray
 * \param [in] bb the user-supplied axis-aligned bounding box
 * \param [in,out] tmin coordinate of closest intersection point along ray
 * \param [in,out] tmax coordinate of farthest intersection point along ray
 * \param [out] ip the intersection point where R intersects bb
 *
 * \return status true iff bb intersects with R, otherwise, false
 *
 * \see primal::Ray
 * \see primal::Segment
 * \see primal::BoundingBox
 *
 * \note Computes Ray Box intersection using the slab method from pg 180 of
 *  Real Time Collision Detection by Christer Ericson.
 */
template <typename T, int DIM>
AXOM_HOST_DEVICE inline bool intersect_ray(const primal::Ray<T, DIM>& R,
                                           const primal::BoundingBox<T, DIM>& bb,
                                           T& tmin,
                                           T& tmax,
                                           T EPS)
{
  AXOM_STATIC_ASSERT(std::is_floating_point<T>::value);

  bool intersects = true;
  for(int d = 0; d < DIM; ++d)
  {
    intersects = intersects &&
      intersect_ray_bbox_test(R.origin()[d],
                              R.direction()[d],
                              bb.getMin()[d],
                              bb.getMax()[d],
                              tmin,
                              tmax,
                              EPS);
  }

  return intersects;
}

template <typename T, int DIM>
AXOM_HOST_DEVICE inline bool intersect_ray(const primal::Ray<T, DIM>& R,
                                           const primal::BoundingBox<T, DIM>& bb,
                                           primal::Point<T, DIM>& ip,
                                           T EPS = numerics::floating_point_limits<T>::epsilon())
{
  AXOM_STATIC_ASSERT(std::is_floating_point<T>::value);

  T tmin = axom::numerics::floating_point_limits<T>::min();
  T tmax = axom::numerics::floating_point_limits<T>::max();

  bool intersects = true;
  for(int d = 0; d < DIM; ++d)
  {
    intersects = intersects &&
      intersect_ray_bbox_test(R.origin()[d],
                              R.direction()[d],
                              bb.getMin()[d],
                              bb.getMax()[d],
                              tmin,
                              tmax,
                              EPS);
  }

  if(intersects)
  {
    ip = R.at(tmin);
  }

  return intersects;
}

template <typename T, int DIM>
AXOM_HOST_DEVICE inline bool intersect_line(const primal::Line<T, DIM>& L,
                                            const primal::BoundingBox<T, DIM>& bb,
                                            primal::Point<T, DIM>& ip,
                                            T EPS = numerics::floating_point_limits<T>::epsilon())
{
  AXOM_STATIC_ASSERT(std::is_floating_point<T>::value);

  T tmin = -axom::numerics::floating_point_limits<T>::max();
  T tmax = axom::numerics::floating_point_limits<T>::max();

  bool intersects = true;
  for(int d = 0; d < DIM; ++d)
  {
    intersects = intersects &&
      intersect_ray_bbox_test(L.origin()[d],
                              L.direction()[d],
                              bb.getMin()[d],
                              bb.getMax()[d],
                              tmin,
                              tmax,
                              EPS);
  }

  if(intersects)
  {
    ip = L.at(tmin);
  }

  return intersects;
}

}  // namespace detail
}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_INTERSECT_RAY_IMPL_HPP_
