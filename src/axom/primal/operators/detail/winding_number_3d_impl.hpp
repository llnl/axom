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
enum class DiscontinuityAxis
{
  x,
  y,
  z,
  rotated
};

#ifdef AXOM_USE_MFEM
/*!
 * \brief Evaluates the integral of the "anti-curl" of the GWN integrand
 *        (via Stokes' theorem) at a point wrt to the trimming curves of a surface
 *
 * \param [in] query The query point
 * \param [in] patch The NURBSPatch object contianing the trimming curves
 * \param [in] ax The axis (relative to query) denoting which anti-curl we use
 * \param [in] npts The number of points used in each Gaussian quadrature
 * \param [in] quad_tol The maximum relative error allowed in each quadrature
 * 
 * Applies a non-adaptive quadrature to the trimming curves of a NURBS patch using one of three possible
 * "anti-curl" vector fields, the curl of each of which is equal to <x, y, z>/||x||^3.
 * 
 * \note This is only meant to be used for `winding_number<NURBSPatch>()`,
 *  and the result does not make sense outside of that context.
 *
 * \return The value of the integral
 */
template <typename T>
double stokes_winding_number_evaluate(const Point<T, 3>& query,
                                      const NURBSPatch<T, 3>& patch,
                                      const DiscontinuityAxis ax,
                                      int npts,
                                      double quad_tol)
{
  // Generate the quadrature rules in parameter space
  static mfem::IntegrationRules my_IntRules(0, mfem::Quadrature1D::GaussLegendre);
  const mfem::IntegrationRule& quad_rule = my_IntRules.Get(mfem::Geometry::SEGMENT, 2 * npts - 1);

  double quad = 0;
  for(int n = 0; n < patch.getNumTrimmingCurves(); ++n)
  {
    NURBSCurve<T, 2> trimming_curve(patch.getTrimmingCurve(n));
    double quad_coarse = stokes_winding_number_component(query,
                                                         trimming_curve,
                                                         patch,
                                                         ax,
                                                         0,
                                                         0,
                                                         quad_rule,
                                                         quad_tol);

    quad += stokes_winding_number_adaptive(query,
                                           trimming_curve,
                                           quad_rule,
                                           patch,
                                           ax,
                                           0,
                                           0,
                                           quad_coarse,
                                           quad_tol);
  }

  return 0.25 * M_1_PI * quad;
}

template <typename T>
double stokes_winding_number_adaptive(const Point<T, 3>& query,
                                      const NURBSCurve<T, 2>& curve,
                                      const mfem::IntegrationRule& quad_rule,
                                      const NURBSPatch<T, 3>& patch,
                                      const DiscontinuityAxis ax,
                                      const int refinement_level,
                                      const int refinement_index,
                                      const double quad_coarse,
                                      const double quad_tol)
{
  double quad_fine_1 = stokes_winding_number_component(query,
                                                       curve,
                                                       patch,
                                                       ax,
                                                       refinement_level + 1,
                                                       2 * refinement_index,
                                                       quad_rule,
                                                       quad_tol);
  double quad_fine_2 = stokes_winding_number_component(query,
                                                       curve,
                                                       patch,
                                                       ax,
                                                       refinement_level + 1,
                                                       2 * refinement_index + 1,
                                                       quad_rule,
                                                       quad_tol);

  if(refinement_level > 25 ||
     axom::utilities::isNearlyEqualRelative(quad_coarse,
                                            quad_fine_1 + quad_fine_2,
                                            quad_tol,
                                            1e-10))
  {
    return quad_fine_1 + quad_fine_2;
  }

  // Perform the adaptive call over both halves of the curve
  quad_fine_1 = stokes_winding_number_adaptive(query,
                                               curve,
                                               quad_rule,
                                               patch,
                                               ax,
                                               refinement_level + 1,
                                               2 * refinement_index,
                                               quad_fine_1,
                                               quad_tol);

  quad_fine_2 = stokes_winding_number_adaptive(query,
                                               curve,
                                               quad_rule,
                                               patch,
                                               ax,
                                               refinement_level + 1,
                                               2 * refinement_index + 1,
                                               quad_fine_2,
                                               quad_tol);

  return quad_fine_1 + quad_fine_2;
}

template <typename T>
double stokes_winding_number_component(const Point<T, 3>& query,
                                       const NURBSCurve<T, 2>& curve,
                                       const NURBSPatch<T, 3>& patch,
                                       const DiscontinuityAxis ax,
                                       const int refinement_level,
                                       const int refinement_index,
                                       const mfem::IntegrationRule& quad_rule,
                                       const double quad_tol)
{
  double this_quad = 0;
  for(int q = 0; q < quad_rule.GetNPoints(); ++q)
  {
    // Find the right knot span based on refinement level
    double span_length =
      (curve.getMaxKnot() - curve.getMinKnot()) / (1 << refinement_level);
    double span_offset = curve.getMinKnot() + span_length * refinement_index;

    double quad_x = span_offset + span_length * quad_rule.IntPoint(q).x;
    double quad_weight = quad_rule.IntPoint(q).weight * span_length;

    Point<T, 2> c_eval;
    Vector<T, 2> c_Dt;
    curve.evaluate_first_derivative(quad_x, c_eval, c_Dt);

    Point<T, 3> s_eval;
    Vector<T, 3> s_Du, s_Dv;
    patch.evaluateFirstDerivatives(c_eval[0], c_eval[1], s_eval, s_Du, s_Dv);

    // Compute quadrature point and total derivative
    const Vector<T, 3> node(query, s_eval);
    const Vector<T, 3> node_dt(s_Du * c_Dt[0] + s_Dv * c_Dt[1]);

    // Compute one of three vector field line integrals depending on
    //  the orientation of the original surface, indicated through ax
    switch(ax)
    {
    case(DiscontinuityAxis::x):
      this_quad += quad_weight *
        (node[2] * node[0] * node_dt[1] - node[1] * node[0] * node_dt[2]) /
        (node[1] * node[1] + node[2] * node[2]) / node.norm();
      break;
    case(DiscontinuityAxis::y):
      this_quad += quad_weight *
        (node[0] * node[1] * node_dt[2] - node[2] * node[1] * node_dt[0]) /
        (node[0] * node[0] + node[2] * node[2]) / node.norm();
      break;
    case(DiscontinuityAxis::z):
    case(DiscontinuityAxis::rotated):
      this_quad += quad_weight *
        (node[1] * node[2] * node_dt[0] - node[0] * node[2] * node_dt[1]) /
        (node[0] * node[0] + node[1] * node[1]) / node.norm();
      break;
    }
  }

  return this_quad;
}
#endif

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif
