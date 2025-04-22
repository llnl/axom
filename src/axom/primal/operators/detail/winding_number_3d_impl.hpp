// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef PRIMAL_WINDING_NUMBER_3D_IMPL_HPP_
#define PRIMAL_WINDING_NUMBER_3D_IMPL_HPP_

// Axom includes
#include "axom/config.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/operators/is_convex.hpp"
#include "axom/primal/operators/squared_distance.hpp"

// C++ includes
#include <math.h>

// MFEM includes
#ifdef AXOM_USE_MFEM
  #include "mfem.hpp"
#endif

namespace axom
{
namespace primal
{
namespace detail
{
/// Type to indicate orientation of singularities relative to surface
enum class SingularityAxis
{
  x,
  y,
  z,
  rotated
};

#ifdef AXOM_USE_MFEM
/*!
 * \brief Evaluates the integral of the "anti-curl" of the GWN integrand
 *        (via Stokes' theorem) at a point wrt to a 3D Bezier curve
 *
 * \param [in] q The query point to test
 * \param [in] curve The BezierCurve object
 * \param [in] ax The axis (relative to query) denoting which anti-curl we use
 * \param [in] npts The number of points used in each Gaussian quadrature
 * \param [in] quad_tol The maximum relative error allowed in each quadrature
 * 
 * Applies a non-adaptive quadrature to a BezierCurve using one of three possible
 * "anti-curl" vector fields, the curl of each of which is equal to <x, y, z>/||x||^3.
 * With the proper "anti-curl" selected, integrating this along a closed curve is equal
 * to evaluating the winding number of a surface with that curve as the boundary.
 * 
 * \note This is only meant to be used for `winding_number<BezierPatch>()`,
 *  and the result does not make sense outside of that context.
 *
 * \return The value of the integral
 */
template <typename T>
double stokes_winding_number(const Point<T, 3>& q,
                             const BezierCurve<T, 3>& curve,
                             const SingularityAxis ax,
                             int npts,
                             double quad_tol)
{
  // Generate the quadrature rules in parameter space
  static mfem::IntegrationRules my_IntRules(0, mfem::Quadrature1D::GaussLegendre);
  const mfem::IntegrationRule& quad_rule = my_IntRules.Get(mfem::Geometry::SEGMENT, 2 * npts - 1);

  double quadrature = 0.0;
  for(int qi = 0; qi < quad_rule.GetNPoints(); ++qi)
  {
    // Get quadrature points in space (shifted by the query)
    const Vector<T, 3> node(q, curve.evaluate(quad_rule.IntPoint(qi).x));
    const Vector<T, 3> node_dt(curve.dt(quad_rule.IntPoint(qi).x));
    const double node_norm = node.norm();

    // Compute one of three vector field line integrals depending on
    //  the orientation of the original surface, indicated through ax.
    switch(ax)
    {
    case(SingularityAxis::x):
      quadrature += quad_rule.IntPoint(qi).weight *
        (node[2] * node[0] * node_dt[1] - node[1] * node[0] * node_dt[2]) /
        (node[1] * node[1] + node[2] * node[2]) / node_norm;
      break;
    case(SingularityAxis::y):
      quadrature += quad_rule.IntPoint(qi).weight *
        (node[0] * node[1] * node_dt[2] - node[2] * node[1] * node_dt[0]) /
        (node[0] * node[0] + node[2] * node[2]) / node_norm;
      break;
    case(SingularityAxis::z):
    case(SingularityAxis::rotated):
      quadrature += quad_rule.IntPoint(qi).weight *
        (node[1] * node[2] * node_dt[0] - node[0] * node[2] * node_dt[1]) /
        (node[0] * node[0] + node[1] * node[1]) / node_norm;
      break;
    }
  }

  // Adaptively refine quadrature over curves if query is not far enough away
  //  from the singularity axis. If rotated, assume you need to adapt.
  bool needs_adapt = false;
  BoundingBox<T, 3> cBox(curve.boundingBox());
  Point<T, 3> centroid = cBox.getCentroid();

  switch(ax)
  {
  case(SingularityAxis::x):
    needs_adapt =
      (q[1] - centroid[1]) * (q[1] - centroid[1]) + (q[2] - centroid[2]) * (q[2] - centroid[2]) <=
      cBox.range().squared_norm();
    break;
  case(SingularityAxis::y):
    needs_adapt =
      (q[0] - centroid[0]) * (q[0] - centroid[0]) + (q[2] - centroid[2]) * (q[2] - centroid[2]) <=
      cBox.range().squared_norm();
    break;
  case(SingularityAxis::z):
    needs_adapt =
      (q[0] - centroid[0]) * (q[0] - centroid[0]) * (q[1] - centroid[1]) * (q[1] - centroid[1]) <=
      cBox.range().squared_norm();
    break;
  case(SingularityAxis::rotated):
    needs_adapt = true;
    break;
  }

  if(needs_adapt)
  {
    return stokes_winding_number_adaptive(q, curve, ax, quad_rule, quadrature, quad_tol);
  }

  return 0.25 * M_1_PI * quadrature;
}
#endif

#ifdef AXOM_USE_MFEM
/*!
 * \brief Accurately evaluates the integral of the "anti-curl" of the GWN integrand
 *        (via Stokes' theorem) at a point wrt to a 3D Bezier curve via recursion
 *
 * \param [in] q The query point to test
 * \param [in] curve The BezierCurve object
 * \param [in] ax The axis (relative to query) denoting which anti-curl we use
 * \param [in] quad_rule The mfem quadrature rule object
 * \param [in] quad_coarse The integral evaluated on the original curve
 * \param [in] quad_tol The maximum relative error allowed in each quadrature
 * \param [in] depth The current recursive depth
 * 
 * Recursively apply quadrature for one of three possible integrals along two halfs 
 * of a curve. The sum of this integral along the subcurves should be equal to to
 * quad_coarse. Otherwise, apply this algorithm to each half recursively.
 *
 * \note This is only meant to be used for `winding_number<BezierPatch>()`,
 *  and the result does not make sense outside of that context.
 * 
 * \return The value of the integral
 */
template <typename T>
double stokes_winding_number_adaptive(const Point<T, 3>& q,
                                      const BezierCurve<T, 3>& curve,
                                      const SingularityAxis ax,
                                      const mfem::IntegrationRule& quad_rule,
                                      const double quad_coarse,
                                      const double quad_tol,
                                      const int depth = 1)
{
  // Split the curve, do the quadrature over both components
  BezierCurve<T, 3> subcurves[2];
  curve.split(0.5, subcurves[0], subcurves[1]);

  double quad_fine[2] = {0.0, 0.0};
  for(int i = 0; i < 2; ++i)
  {
    for(int qi = 0; qi < quad_rule.GetNPoints(); ++qi)
    {
      // Get quad_rulerature points in space (shifted by the query)
      const Vector<T, 3> node(q, subcurves[i].evaluate(quad_rule.IntPoint(qi).x));
      const Vector<T, 3> node_dt(subcurves[i].dt(quad_rule.IntPoint(qi).x));
      const double node_norm = node.norm();

      // Compute one of three vector field line integrals depending on
      //  the orientation of the original surface, indicated through ax.
      switch(ax)
      {
      case(SingularityAxis::x):
        quad_fine[i] += quad_rule.IntPoint(qi).weight *
          (node[2] * node[0] * node_dt[1] - node[1] * node[0] * node_dt[2]) /
          (node[1] * node[1] + node[2] * node[2]) / node_norm;
        break;
      case(SingularityAxis::y):
        quad_fine[i] += quad_rule.IntPoint(qi).weight *
          (node[0] * node[1] * node_dt[2] - node[2] * node[1] * node_dt[0]) /
          (node[0] * node[0] + node[2] * node[2]) / node_norm;
        break;
      case(SingularityAxis::z):
      case(SingularityAxis::rotated):
        quad_fine[i] += quad_rule.IntPoint(qi).weight *
          (node[1] * node[2] * node_dt[0] - node[0] * node[2] * node_dt[1]) /
          (node[0] * node[0] + node[1] * node[1]) / node_norm;
        break;
      }
    }
  }

  constexpr int MAX_DEPTH = 12;
  if(depth >= MAX_DEPTH ||
     axom::utilities::isNearlyEqualRelative(quad_fine[0] + quad_fine[1], quad_coarse, quad_tol, 1e-10))
  {
    return 0.25 * M_1_PI * (quad_fine[0] + quad_fine[1]);
  }
  else
  {
    return stokes_winding_number_adaptive(q,
                                          subcurves[0],
                                          ax,
                                          quad_rule,
                                          quad_fine[0],
                                          quad_tol,
                                          depth + 1) +
      stokes_winding_number_adaptive(q, subcurves[1], ax, quad_rule, quad_fine[1], quad_tol, depth + 1);
  }
}
#endif

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif
