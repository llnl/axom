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
 * \param [out] tp Parametric coordinates of intersections in \a line
 * \param [out] up Parametric coordinates of intersections in \a patch
 * \param [out] vp Parametric coordinates of intersections in \a patch
 * \param [in] order_u The order of \a line in the u direction
 * \param [in] order_v The order of \a line in the v direction
 * \param [in] u_offset The offset in parameter space for \a patch in the u direction
 * \param [in] u_scale The scale in parameter space for \a patch in the u direction
 * \param [in] v_offset The offset in parameter space for \a patch in the v direction
 * \param [in] v_scale The scale in parameter space for \a patch in the v direction
 * \param [in] sq_tol Numerical tolerance for physical distances
 * \param [in] EPS Numerical tolerance in parameter space
 *
 * A ray can only intersect a Bezier patch if it intersects its bounding box.
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
 * \return True if the line intersects the patch, False otherwise
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
                          bool isRay);

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
                          bool isRay)
{
  using BPatch = BezierPatch<T, 3>;

  // Early return if we start to record excessive intersections.
  //  This implies the patch is degenerate at the point of intersection.
  constexpr int max_intersections = 100;
  if(tp.size() >= max_intersections)
  {
    SLIC_WARNING(
      "Unexpectedly large number of intersections recorded."
      "Returning early to avoid excessive computation.");
    return false;
  }

  // Check bounding box to short-circuit the intersection
  //  Need to expand the box a bit so that intersections near subdivision boundaries
  //  are accurately recorded
  if(!intersect(line, patch.boundingBox().scale(1.5)))
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

    BPatch p1(order_u, order_v), p2(order_u, order_v), p3(order_u, order_v), p4(order_u, order_v);

    patch.split(splitVal, splitVal, p1, p2, p3, p4);
    u_scale *= scaleFac;
    v_scale *= scaleFac;

    // Note: we want to find all intersections, so don't short-circuit
    intersect_line_patch(line,
                         p1,
                         tp,
                         up,
                         vp,
                         order_u,
                         order_v,
                         u_offset,
                         u_scale,
                         v_offset,
                         v_scale,
                         sq_tol,
                         EPS,
                         isRay);

    intersect_line_patch(line,
                         p2,
                         tp,
                         up,
                         vp,
                         order_u,
                         order_v,
                         u_offset + u_scale,
                         u_scale,
                         v_offset,
                         v_scale,
                         sq_tol,
                         EPS,
                         isRay);

    intersect_line_patch(line,
                         p3,
                         tp,
                         up,
                         vp,
                         order_u,
                         order_v,
                         u_offset,
                         u_scale,
                         v_offset + v_scale,
                         v_scale,
                         sq_tol,
                         EPS,
                         isRay);

    intersect_line_patch(line,
                         p4,
                         tp,
                         up,
                         vp,
                         order_u,
                         order_v,
                         u_offset + u_scale,
                         u_scale,
                         v_offset + v_scale,
                         v_scale,
                         sq_tol,
                         EPS,
                         isRay);
  }

  return tp.size() > 0;
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif  // AXOM_PRIMAL_INTERSECT_PATCH_IMPL_HPP_
