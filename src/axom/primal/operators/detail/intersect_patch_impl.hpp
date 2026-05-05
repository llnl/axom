// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file intersect_patch_impl.hpp
 *
 * This file provides helper functions for testing the intersection
 * of rays and Bezier patches
 */

#ifndef AXOM_PRIMAL_INTERSECT_PATCH_IMPL_HPP_
#define AXOM_PRIMAL_INTERSECT_PATCH_IMPL_HPP_

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"

#include "axom/primal/operators/intersect.hpp"
#include "axom/primal/operators/in_polygon.hpp"
#include "axom/primal/operators/detail/intersect_impl.hpp"
#include "axom/primal/operators/detail/intersect_ray_impl.hpp"

#include <vector>

namespace axom
{
namespace primal
{
namespace detail
{
//---------------------------- FUNCTION DECLARATIONS ---------------------------

/*!
 * \brief Recursive function to find the intersections between a line and a Bezier patch
 *
 * \param [in] line The input line
 * \param [in] patch The input patch
 * \param [out] tp Arrays to append parametric coordinates of intersections in \a line
 * \param [out] up Arrays to append parametric coordinates of intersections in \a patch
 * \param [out] vp Arrays to append parametric coordinates of intersections in \a patch
 * \param [in] order_u The order of \a line in the u direction
 * \param [in] order_v The order of \a line in the v direction
 * \param [in] u_offset The offset in parameter space for \a patch in the u direction
 * \param [in] u_scale The scale in parameter space for \a patch in the u direction
 * \param [in] v_offset The offset in parameter space for \a patch in the v direction
 * \param [in] v_scale The scale in parameter space for \a patch in the v direction
 * \param [in] sq_tol Numerical tolerance for physical distances
 * \param [in] EPS Numerical tolerance in parameter space
 * \param [in] isRay True if the line is a ray, i.e., only return nonnegative t values
 * \param [out] success False if an early return was triggered
 *
 * A line can only intersect a Bezier patch if it intersects its bounding box.
 * The base case of the recursion is when we can approximate the patch as
 * parametrically bilinear, where we directly find their intersections. Otherwise,
 * check for intersections recursively after bisecting the patch in each direction.
 *
 * \note This detail function returns all found intersections within EPS of parameter space,
 *  including identical intersections reported by each subdivision. 
 * The calling `intersect` routine should remove duplicates and enforce half-open behavior. 
 * 
 * \note This function assumes that all intersections have multiplicity 
 *  one, i.e. does not find tangencies. 
 *
 * \warning This function returns early if we record excessive intersections.
 *    This implies the patch is degenerate at the point of intersection.
 * 
 * \return False if an early return was triggered (failure). True otherwise
 */
template <typename T>
bool intersect_line_patch(const Line<T, 3> &line,
                          const BezierPatch<T, 3> &patch,
                          axom::Array<T> &tp,
                          axom::Array<T> &up,
                          axom::Array<T> &vp,
                          int order_u,
                          int order_v,
                          double u_offset,
                          double u_scale,
                          double v_offset,
                          double v_scale,
                          double sq_tol,
                          double EPS,
                          bool isRay,
                          bool &success);

//------------------------------ IMPLEMENTATIONS ------------------------------

template <typename T>
bool intersect_line_patch(const Line<T, 3> &line,
                          const BezierPatch<T, 3> &patch,
                          axom::Array<T> &tp,
                          axom::Array<T> &up,
                          axom::Array<T> &vp,
                          int order_u,
                          int order_v,
                          double u_offset,
                          double u_scale,
                          double v_offset,
                          double v_scale,
                          double sq_tol,
                          double EPS,
                          bool isRay,
                          bool &success)
{
  using BPatch = BezierPatch<T, 3>;

  // Early return if we start to record excessive intersections.
  //  This implies the patch is degenerate at the point of intersection.
  if(tp.size() > order_v * order_u)
  {
    success = false;
    return true;
  }

  // Check bounding box to skip the subdivision procedure
  // Expand the box a bit so that intersections near subdivision boundaries are accurately recorded
  if(!intersect(line, patch.boundingBox().expand(10 * EPS)))
  {
    return false;
  }

  if(patch.isBilinear(sq_tol, true))
  {
    // Store candidate intersection points
    axom::StaticArray<T, 2> tc, uc, vc;

    detail::intersect_line_bilinear_patch(line,
                                          patch(0, 0),
                                          patch(order_u, 0),
                                          patch(order_u, order_v),
                                          patch(0, order_v),
                                          tc,
                                          uc,
                                          vc,
                                          EPS,
                                          isRay);

    // Check intersections based on subdivision-scaled tolerances
    for(int i = 0; i < tc.size(); ++i)
    {
      const T t0 = tc[i];
      const T u0 = uc[i];
      const T v0 = vc[i];

      // Use EPS to record points near the boundary of the bilinear approximation
      if(u0 >= -EPS / u_scale && u0 <= 1.0 + EPS / u_scale && v0 >= -EPS / v_scale &&
         v0 <= 1.0 + EPS / v_scale)
      {
        if(t0 >= -EPS || !isRay)
        {
          up.push_back(u_offset + u0 * u_scale);
          vp.push_back(v_offset + v0 * v_scale);
          tp.push_back(t0);
        }
      }
    }
  }
  else
  {
    constexpr double splitVal = 0.5;
    constexpr double scaleFac = 0.5;

    BPatch subpatches[4];
    patch.split(splitVal, splitVal, subpatches[0], subpatches[1], subpatches[2], subpatches[3]);

    u_scale *= scaleFac;
    v_scale *= scaleFac;

    for(int i = 0; i < 4; ++i)
    {
      // If we already found a degenerate intersection, we can skip the rest
      if(!success)
      {
        return !tp.empty();
      }

      // Check all four subpatches even if intersections are found, as we want to find them all
      intersect_line_patch(line,
                           subpatches[i],
                           tp,
                           up,
                           vp,
                           order_u,
                           order_v,
                           u_offset + (i % 2) * u_scale,
                           u_scale,
                           v_offset + (i / 2) * v_scale,
                           v_scale,
                           sq_tol,
                           EPS,
                           isRay,
                           success);
    }
  }

  return !tp.empty();
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif  // AXOM_PRIMAL_INTERSECT_PATCH_IMPL_HPP_
