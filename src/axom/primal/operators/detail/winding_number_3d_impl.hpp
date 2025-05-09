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
 * \brief Identify the u/v isoline on which all degenerate intersections occur, 
 *         and "clip out" patches that do not contain this line
 *
 * \param [in] patch The NURBS patch
 * \param [in] up, vp The arrays of intersection coordinates in parameter space 
 * \param [in] clip_radius The width of the strip which is removed in parameter space 
 * \param [out] out_patch1, out_patch2 The patches which are returned on either side of the strip
 * 
 * \note If the relevant isoline occurs within `clip_radius` of a patch edge, 
 *         the correspondong out_patch will be invalid
 * 
 * \return The clipped patch
 */
template <typename T>
void degenerate_surface_processing(const NURBSPatch<T, 3>& patch,
                                   const axom::Array<T> up,
                                   const axom::Array<T> vp,
                                   const T clip_radius,
                                   NURBSPatch<T, 3>& out_patch1,
                                   NURBSPatch<T, 3>& out_patch2)
{
  T mean_u = up[0], var_u = 0.0;
  T mean_v = vp[0], var_v = 0.0;

  // Iterate through the coordinates to identify the correct u/v line
  for(int i = 1; i < up.size(); ++i)
  {
    T new_mean_u = mean_u + (up[i] - mean_u) / static_cast<T>(i + 1);
    T new_mean_v = mean_v + (vp[i] - mean_v) / static_cast<T>(i + 1);

    T new_var_u = var_u + (up[i] - mean_u) * (up[i] - new_mean_u);
    T new_var_v = var_v + (vp[i] - mean_v) * (vp[i] - new_mean_v);

    mean_u = new_mean_u;
    mean_v = new_mean_v;

    var_u = new_var_u;
    var_v = new_var_v;
  }

  NURBSPatch<T, 3> dummy_patch(patch);

  // Indicates a u isocurve
  if(var_u < var_v)
  {
    if(mean_u - clip_radius > patch.getMinKnot_u())
    {
      dummy_patch.split_u(mean_u - clip_radius, out_patch1, dummy_patch);
    }
    else
    {
      out_patch1 = NURBSPatch<T, 3>();
    }

    if(mean_u + clip_radius < patch.getMaxKnot_u())
    {
      dummy_patch.split_u(mean_u + clip_radius, dummy_patch, out_patch2);
    }
    else
    {
      out_patch2 = NURBSPatch<T, 3>();
    }
  }
  else
  {
    if(mean_v - clip_radius > patch.getMinKnot_v())
    {
      dummy_patch.split_v(mean_v - clip_radius, out_patch1, dummy_patch);
    }
    else
    {
      out_patch1 = NURBSPatch<T, 3>();
    }

    if(mean_v + clip_radius < patch.getMaxKnot_v())
    {
      dummy_patch.split_v(mean_v + clip_radius, dummy_patch, out_patch2);
    }
    else
    {
      out_patch2 = NURBSPatch<T, 3>();
    }
  }
}

/*!
 * \brief Evaluates the integral of the "anti-curl" of the GWN integrand
 *        (via Stokes' theorem) at a point wrt to a 3D Bezier curve
 *
 * \param [in] query The query point to test
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
double stokes_winding_number(const Point<T, 3>& query,
                             const BezierCurve<T, 3>& curve,
                             const DiscontinuityAxis ax,
                             int npts,
                             double quad_tol)
{
  // Generate the quadrature rules in parameter space
  static mfem::IntegrationRules my_IntRules(0, mfem::Quadrature1D::GaussLegendre);
  const mfem::IntegrationRule& quad_rule =
    my_IntRules.Get(mfem::Geometry::SEGMENT, 2 * npts - 1);

  double quadrature = 0.0;
  for(int q = 0; q < quad_rule.GetNPoints(); ++q)
  {
    // Get quadrature points in space (shifted by the query)
    const Vector<T, 3> node(query, curve.evaluate(quad_rule.IntPoint(q).x));
    const Vector<T, 3> node_dt(curve.dt(quad_rule.IntPoint(q).x));
    const double node_norm = node.norm();

    // Compute one of three vector field line integrals depending on
    //  the orientation of the original surface, indicated through ax.
    switch(ax)
    {
    case(DiscontinuityAxis::x):
      quadrature += quad_rule.IntPoint(q).weight *
        (node[2] * node[0] * node_dt[1] - node[1] * node[0] * node_dt[2]) /
        (node[1] * node[1] + node[2] * node[2]) / node_norm;
      break;
    case(DiscontinuityAxis::y):
      quadrature += quad_rule.IntPoint(q).weight *
        (node[0] * node[1] * node_dt[2] - node[2] * node[1] * node_dt[0]) /
        (node[0] * node[0] + node[2] * node[2]) / node_norm;
      break;
    case(DiscontinuityAxis::z):
    case(DiscontinuityAxis::rotated):
      quadrature += quad_rule.IntPoint(q).weight *
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
  case(DiscontinuityAxis::x):
    needs_adapt = (query[1] - centroid[1]) * (query[1] - centroid[1]) +
        (query[2] - centroid[2]) * (query[2] - centroid[2]) <=
      cBox.range().squared_norm();
    break;
  case(DiscontinuityAxis::y):
    needs_adapt = (query[0] - centroid[0]) * (query[0] - centroid[0]) +
        (query[2] - centroid[2]) * (query[2] - centroid[2]) <=
      cBox.range().squared_norm();
    break;
  case(DiscontinuityAxis::z):
    needs_adapt = (query[0] - centroid[0]) * (query[0] - centroid[0]) +
        (query[1] - centroid[1]) * (query[1] - centroid[1]) <=
      cBox.range().squared_norm();
    break;
  case(DiscontinuityAxis::rotated):
    needs_adapt = true;
    break;
  }

  if(needs_adapt)
  {
    return stokes_winding_number_adaptive(query,
                                          curve,
                                          ax,
                                          quad_rule,
                                          quadrature,
                                          quad_tol);
  }

  return 0.25 * M_1_PI * quadrature;
}
#endif

#ifdef AXOM_USE_MFEM
/*!
 * \brief Accurately evaluates the integral of the "anti-curl" of the GWN integrand
 *        (via Stokes' theorem) at a point wrt to a 3D Bezier curve via recursion
 *
 * \param [in] query The query point to test
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
double stokes_winding_number_adaptive(const Point<T, 3>& query,
                                      const BezierCurve<T, 3>& curve,
                                      const DiscontinuityAxis ax,
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
    for(int q = 0; q < quad_rule.GetNPoints(); ++q)
    {
      // Get quad_rulerature points in space (shifted by the query)
      const Vector<T, 3> node(query,
                              subcurves[i].evaluate(quad_rule.IntPoint(q).x));
      const Vector<T, 3> node_dt(subcurves[i].dt(quad_rule.IntPoint(q).x));
      const double node_norm = node.norm();

      // Compute one of three vector field line integrals depending on
      //  the orientation of the original surface, indicated through ax.
      switch(ax)
      {
      case(DiscontinuityAxis::x):
        quad_fine[i] += quad_rule.IntPoint(q).weight *
          (node[2] * node[0] * node_dt[1] - node[1] * node[0] * node_dt[2]) /
          (node[1] * node[1] + node[2] * node[2]) / node_norm;
        break;
      case(DiscontinuityAxis::y):
        quad_fine[i] += quad_rule.IntPoint(q).weight *
          (node[0] * node[1] * node_dt[2] - node[2] * node[1] * node_dt[0]) /
          (node[0] * node[0] + node[2] * node[2]) / node_norm;
        break;
      case(DiscontinuityAxis::z):
      case(DiscontinuityAxis::rotated):
        quad_fine[i] += quad_rule.IntPoint(q).weight *
          (node[1] * node[2] * node_dt[0] - node[0] * node[2] * node_dt[1]) /
          (node[0] * node[0] + node[1] * node[1]) / node_norm;
        break;
      }
    }
  }

  constexpr int MAX_DEPTH = 15;
  if(depth >= MAX_DEPTH ||
     axom::utilities::isNearlyEqualRelative(quad_fine[0] + quad_fine[1],
                                            quad_coarse,
                                            quad_tol,
                                            1e-10))
  {
    return 0.25 * M_1_PI * (quad_fine[0] + quad_fine[1]);
  }
  else
  {
    return stokes_winding_number_adaptive(query,
                                          subcurves[0],
                                          ax,
                                          quad_rule,
                                          quad_fine[0],
                                          quad_tol,
                                          depth + 1) +
      stokes_winding_number_adaptive(query,
                                     subcurves[1],
                                     ax,
                                     quad_rule,
                                     quad_fine[1],
                                     quad_tol,
                                     depth + 1);
  }
}

const bool FORCE_ADAPTIVE = true;

template <typename T>
double stokes_winding_number(const Point<T, 3>& query,
                             const NURBSPatch<T, 3>& patch,
                             const DiscontinuityAxis ax,
                             int npts,
                             double quad_tol)
{
  // Generate the quadrature rules in parameter space
  static mfem::IntegrationRules my_IntRules(0, mfem::Quadrature1D::GaussLegendre);
  const mfem::IntegrationRule& quad_rule =
    my_IntRules.Get(mfem::Geometry::SEGMENT, 2 * npts - 1);

  double quad = 0;
  for(int n = 0; n < patch.getNumTrimmingCurves(); ++n)
  {
    NURBSCurve<T, 2> trimming_curve(patch.getTrimmingCurve(n));
    quad += stokes_winding_number_trimming(query,
                                           trimming_curve,
                                           patch,
                                           ax,
                                           quad_rule,
                                           quad_tol);
  }

  return quad;
}

template <typename T>
double stokes_winding_number_trimming(const Point<T, 3>& query,
                                      const NURBSCurve<T, 2>& curve,
                                      const NURBSPatch<T, 3>& patch,
                                      const DiscontinuityAxis ax,
                                      const mfem::IntegrationRule& quad_rule,
                                      const double quad_tol)
{
  // Track an approximate distance of the curve to axis
  //  to determine if we need to adapt
  double min_dist = std::numeric_limits<double>::max();
  double max_dist = 0;

  const T min_knot = curve.getKnot(0);
  const T max_knot = curve.getKnot(curve.getNumKnots() - 1);

  int q_npts = quad_rule.GetNPoints();
  double quad = 0;

  for(int q = 0; q < q_npts; ++q)
  {
    // Get quadrature node shifted to knot span
    T quad_x = quad_rule.IntPoint(q).x * (max_knot - min_knot) + min_knot;
    T quad_weight = quad_rule.IntPoint(q).weight * (max_knot - min_knot);

    // Get quadrature points in space (shifted by the query)
    Point<T, 2> p_eval;
    Vector<T, 2> p_Dt;
    curve.evaluate_first_derivative(quad_x, p_eval, p_Dt);

    Point<T, 3> s_eval;
    Vector<T, 3> s_Du, s_Dv;
    patch.evaluateFirstDerivatives(p_eval[0], p_eval[1], s_eval, s_Du, s_Dv);

    const Vector<T, 3> node(query, s_eval);
    const double node_norm = node.norm();

    // Compute total derivative
    const Vector<T, 3> node_dt(s_Du * p_Dt[0] + s_Dv * p_Dt[1]);

    // Compute one of three vector field line integrals depending on
    //  the orientation of the original surface, indicated through ax.
    double denom;
    switch(ax)
    {
    case(DiscontinuityAxis::x):
      denom = node[1] * node[1] + node[2] * node[2];
      quad += quad_weight *
        (node[2] * node[0] * node_dt[1] - node[1] * node[0] * node_dt[2]) /
        denom / node_norm;
      break;
    case(DiscontinuityAxis::y):
      denom = node[0] * node[0] + node[2] * node[2];
      quad += quad_weight *
        (node[0] * node[1] * node_dt[2] - node[2] * node[1] * node_dt[0]) /
        denom / node_norm;
      break;
    case(DiscontinuityAxis::z):
    case(DiscontinuityAxis::rotated):
      denom = node[0] * node[0] + node[1] * node[1];
      quad += quad_weight *
        (node[1] * node[2] * node_dt[0] - node[0] * node[2] * node_dt[1]) /
        denom / node_norm;
      break;
    }

    min_dist = std::min(min_dist, denom);
    max_dist = std::max(max_dist, denom);
  }

  // Adaptively refine quadrature over curves if the minimum distance
  //  from the axis to the curve is within 10% of the difference between
  //  the maximum and minimum distances (this is a bad heuristic).
  if(FORCE_ADAPTIVE || true)  //min_dist <= 0.1 * (max_dist - min_dist))
  {
    return stokes_winding_number_trimming_adaptive(query,
                                                   curve,
                                                   patch,
                                                   ax,
                                                   quad_rule,
                                                   quad,
                                                   quad_tol);
  }

  return 0.25 * M_1_PI * quad;
}

template <typename T>
double stokes_winding_number_trimming_adaptive(const Point<T, 3>& query,
                                               const NURBSCurve<T, 2>& curve,
                                               const NURBSPatch<T, 3>& patch,
                                               const DiscontinuityAxis ax,
                                               const mfem::IntegrationRule& quad_rule,
                                               const double quad_coarse,
                                               const double quad_tol,
                                               const int depth = 1)
{
  // Split the curve, do the quadrature over both components
  NURBSCurve<T, 2> subcurves[2];

  const T min_knot = curve.getKnot(0);
  const T max_knot = curve.getKnot(curve.getNumKnots() - 1);
  const T mid_knot = 0.5 * (min_knot + max_knot);

  curve.split(mid_knot, subcurves[0], subcurves[1]);

  double quad_fine[2] = {0.0, 0.0};
  for(int i = 0; i < 2; ++i)
  {
    T the_min_knot = (i == 0) ? min_knot : mid_knot;
    T the_max_knot = (i == 0) ? mid_knot : max_knot;

    for(int q = 0; q < quad_rule.GetNPoints(); ++q)
    {
      const T quad_x =
        quad_rule.IntPoint(q).x * (the_max_knot - the_min_knot) + the_min_knot;
      const T quad_weight =
        quad_rule.IntPoint(q).weight * (the_max_knot - the_min_knot);

      // Get quadrature points in space (shifted by the query)
      Point<T, 2> p_eval;
      Vector<T, 2> p_Dt;
      subcurves[i].evaluate_first_derivative(quad_x, p_eval, p_Dt);

      Point<T, 3> s_eval;
      Vector<T, 3> s_Du, s_Dv;
      patch.evaluateFirstDerivatives(p_eval[0], p_eval[1], s_eval, s_Du, s_Dv);

      const Vector<T, 3> node(query, s_eval);
      const double node_norm = node.norm();

      // Compute total derivative
      const Vector<T, 3> node_dt(s_Du * p_Dt[0] + s_Dv * p_Dt[1]);

      // Compute one of three vector field line integrals depending on
      //  the orientation of the original surface, indicated through ax.
      switch(ax)
      {
      case(DiscontinuityAxis::x):
        quad_fine[i] += quad_weight *
          (node[2] * node[0] * node_dt[1] - node[1] * node[0] * node_dt[2]) /
          (node[1] * node[1] + node[2] * node[2]) / node_norm;
        break;
      case(DiscontinuityAxis::y):
        quad_fine[i] += quad_weight *
          (node[0] * node[1] * node_dt[2] - node[2] * node[1] * node_dt[0]) /
          (node[0] * node[0] + node[2] * node[2]) / node_norm;
        break;
      case(DiscontinuityAxis::z):
      case(DiscontinuityAxis::rotated):
        quad_fine[i] += quad_weight *
          (node[1] * node[2] * node_dt[0] - node[0] * node[2] * node_dt[1]) /
          (node[0] * node[0] + node[1] * node[1]) / node_norm;
        break;
      }
    }
  }

  constexpr int MAX_DEPTH = 15;
  if(depth >= MAX_DEPTH ||
     axom::utilities::isNearlyEqualRelative(quad_fine[0] + quad_fine[1],
                                            quad_coarse,
                                            quad_tol,
                                            1e-10))
  {
    return 0.25 * M_1_PI * (quad_fine[0] + quad_fine[1]);
  }
  else
  {
    return stokes_winding_number_trimming_adaptive(query,
                                                   subcurves[0],
                                                   patch,
                                                   ax,
                                                   quad_rule,
                                                   quad_fine[0],
                                                   quad_tol,
                                                   depth + 1) +
      stokes_winding_number_trimming_adaptive(query,
                                              subcurves[1],
                                              patch,
                                              ax,
                                              quad_rule,
                                              quad_fine[1],
                                              quad_tol,
                                              depth + 1);
  }
}

template <typename T>
bool isNearAxisBox(const Point<T, 3>& query,
                   const BoundingBox<T, 3>& bbox,
                   const double beta,
                   const DiscontinuityAxis ax)
{
  auto centroid = bbox.getCentroid();
  Vector<T, 3> c_query(centroid, query);
  Vector<T, 3> box_max(centroid, bbox.getMax());

  double distance_to_axis, bbox_radius;
  switch(ax)
  {
  case(DiscontinuityAxis::x):
    distance_to_axis = c_query[1] * c_query[1] + c_query[2] * c_query[2];
    bbox_radius = box_max[1] * box_max[1] + box_max[2] * box_max[2];
    return distance_to_axis <= beta * bbox_radius;
  case(DiscontinuityAxis::y):
    distance_to_axis = c_query[0] * c_query[0] + c_query[2] * c_query[2];
    bbox_radius = box_max[0] * box_max[0] + box_max[2] * box_max[2];
    return distance_to_axis <= beta * bbox_radius;
  case(DiscontinuityAxis::z):
    distance_to_axis = c_query[0] * c_query[0] + c_query[1] * c_query[1];
    bbox_radius = box_max[0] * box_max[0] + box_max[1] * box_max[1];
    return distance_to_axis <= beta * bbox_radius;
  }

  // Make "true" the default for more reliable quadrature
  return true;
}

template <typename T>
double stokes_winding_number_cached(const Point<T, 3>& query,
                                    const NURBSPatchData<T>& nPatchData,
                                    const DiscontinuityAxis ax,
                                    const int quad_npts,
                                    const double quad_tol)
{
  // Generate the quadrature rules in parameter space
  static mfem::IntegrationRules my_IntRules(0, mfem::Quadrature1D::GaussLegendre);
  const mfem::IntegrationRule& quad_rule =
    my_IntRules.Get(mfem::Geometry::SEGMENT, 2 * quad_npts - 1);

  double quad = 0;
  for(int n = 0; n < nPatchData.patch.getNumTrimmingCurves(); ++n)
  {
    // Get the quadrature points for the curve without any refinement
    auto trimming_curve_data = nPatchData.getQuadratureData(n, quad_rule, 0, 0);
    double quad_coarse =
      single_stokes_cached(query, quad_rule, ax, trimming_curve_data);

    // If near, use adaptive
    if(FORCE_ADAPTIVE || isNearAxisBox(query, trimming_curve_data.bbox, 5.0, ax))
    {
      quad += 0.25 * M_1_PI *
        adaptive_stokes_cached(query,
                               nPatchData,
                               quad_rule,
                               ax,
                               n,
                               0,
                               0,
                               quad_coarse,
                               quad_tol);
    }
    // Otherwise, do some refinement
    else
    {
      quad += 0.25 * M_1_PI * quad_coarse;
    }
  }

  return quad;
}

template <typename T>
double adaptive_stokes_cached(const Point<T, 3>& query,
                              const NURBSPatchData<T>& nPatchData,
                              const mfem::IntegrationRule& quad_rule,
                              const DiscontinuityAxis ax,
                              const int curve_index,
                              const int refinement_level,
                              const int refinement_index,
                              const double quad_coarse,
                              const double quad_tol)
{
  auto trimming_curve_data_1 =
    nPatchData.getQuadratureData(curve_index,
                                 quad_rule,
                                 refinement_level + 1,
                                 2 * refinement_index);
  auto trimming_curve_data_2 =
    nPatchData.getQuadratureData(curve_index,
                                 quad_rule,
                                 refinement_level + 1,
                                 2 * refinement_index + 1);

  double quad_fine_1 =
    single_stokes_cached(query, quad_rule, ax, trimming_curve_data_1);
  double quad_fine_2 =
    single_stokes_cached(query, quad_rule, ax, trimming_curve_data_2);

  if(refinement_level >= 15 ||
     axom::utilities::isNearlyEqualRelative(quad_fine_1 + quad_fine_2,
                                            quad_coarse,
                                            quad_tol,
                                            1e-10))
  {
    return quad_fine_1 + quad_fine_2;
  }
  else
  {
    // If we're not near the axis, we can trust the original value
    if(FORCE_ADAPTIVE ||
       isNearAxisBox(query, trimming_curve_data_1.bbox, 5.0, ax))
      quad_fine_1 = adaptive_stokes_cached(query,
                                           nPatchData,
                                           quad_rule,
                                           ax,
                                           curve_index,
                                           refinement_level + 1,
                                           2 * refinement_index,
                                           quad_fine_1,
                                           quad_tol);

    if(FORCE_ADAPTIVE ||
       isNearAxisBox(query, trimming_curve_data_2.bbox, 5.0, ax))
      quad_fine_2 = adaptive_stokes_cached(query,
                                           nPatchData,
                                           quad_rule,
                                           ax,
                                           curve_index,
                                           refinement_level + 1,
                                           2 * refinement_index + 1,
                                           quad_fine_2,
                                           quad_tol);

    return quad_fine_1 + quad_fine_2;
  }
}

template <typename T>
double single_stokes_cached(const Point<T, 3>& query,
                            const mfem::IntegrationRule& quad_rule,
                            const DiscontinuityAxis ax,
                            const TrimmingCurveQuadratureData<T>& trimming_curve_data)
{
  // Do this without refinement
  double this_quad = 0;
  for(int q = 0; q < quad_rule.GetNPoints(); ++q)
  {
    const Vector<T, 3> node(query,
                            trimming_curve_data.quadrature_points[q].first);
    const Vector<T, 3> node_dt(trimming_curve_data.quadrature_points[q].second);
    const double node_norm = node.norm();

    const double quad_weight =
      quad_rule.IntPoint(q).weight * trimming_curve_data.span_length;

    // Compute one of three vector field line integrals depending on
    //  the orientation of the original surface, indicated through ax.
    switch(ax)
    {
    case(DiscontinuityAxis::x):
      this_quad += quad_weight *
        (node[2] * node[0] * node_dt[1] - node[1] * node[0] * node_dt[2]) /
        (node[1] * node[1] + node[2] * node[2]) / node_norm;
      break;
    case(DiscontinuityAxis::y):
      this_quad += quad_weight *
        (node[0] * node[1] * node_dt[2] - node[2] * node[1] * node_dt[0]) /
        (node[0] * node[0] + node[2] * node[2]) / node_norm;
      break;
    case(DiscontinuityAxis::z):
      this_quad += quad_weight *
        (node[1] * node[2] * node_dt[0] - node[0] * node[2] * node_dt[1]) /
        (node[0] * node[0] + node[1] * node[1]) / node_norm;
      break;
    }
  }

  return this_quad;
}

// --=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--

template <typename T>
Point<T, 3> rotate_point(const numerics::Matrix<T>& matx,
                         const Point<T, 3>& center,
                         const Point<T, 3>& input)
{
  Vector<T, 3> shifted(center, input);
  Vector<T, 3> rotated;
  numerics::matrix_vector_multiply(matx, shifted.data(), rotated.data());
  return Point<T, 3>(
    {rotated[0] + center[0], rotated[1] + center[1], rotated[2] + center[2]});
}

template <typename T>
Point<T, 3> rotate_vector(const numerics::Matrix<T>& matx,
                          const Point<T, 3>& center,
                          const Vector<T, 3>& input)
{
  Vector<T, 3> shifted {input[0] - center[0],
                        input[1] - center[1],
                        input[2] - center[2]};
  Vector<T, 3> rotated;
  numerics::matrix_vector_multiply(matx, shifted.data(), rotated.data());
  return Point<T, 3>(
    {rotated[0] + center[0], rotated[1] + center[1], rotated[2] + center[2]});
}

template <typename T>
bool isNearAxisBoxRotated(const Point<T, 3>& query,
                          const BoundingBox<T, 3>& bbox,
                          const double beta,
                          const axom::numerics::Matrix<T>& rotator)
{
  auto centroid = rotate_point(rotator, query, bbox.getCentroid());
  Vector<T, 3> c_query(centroid, query);
  Vector<T, 3> box_max(centroid, rotate_point(rotator, query, bbox.getMax()));

  double distance_to_axis, bbox_radius;
  distance_to_axis = c_query[0] * c_query[0] + c_query[1] * c_query[1];
  bbox_radius = box_max[0] * box_max[0] + box_max[1] * box_max[1];
  return distance_to_axis <= beta * bbox_radius;
}

template <typename T>
double stokes_winding_number_cached_rotated(const Point<T, 3>& query,
                                            const NURBSPatchData<T>& nPatchData,
                                            const axom::numerics::Matrix<T>& rotator,
                                            const int quad_npts,
                                            const double quad_tol)
{
  // Generate the quadrature rules in parameter space
  static mfem::IntegrationRules my_IntRules(0, mfem::Quadrature1D::GaussLegendre);
  const mfem::IntegrationRule& quad_rule =
    my_IntRules.Get(mfem::Geometry::SEGMENT, 2 * quad_npts - 1);

  double quad = 0;
  for(int n = 0; n < nPatchData.patch.getNumTrimmingCurves(); ++n)
  {
    // Get the quadrature points for the curve without any refinement
    auto trimming_curve_data = nPatchData.getQuadratureData(n, quad_rule, 0, 0);
    double quad_coarse =
      single_stokes_cached_rotated(query, quad_rule, rotator, trimming_curve_data);

    // If far away, use this value
    if(FORCE_ADAPTIVE ||
       isNearAxisBoxRotated(query, trimming_curve_data.bbox, 5.0, rotator))
    {
      quad += 0.25 * M_1_PI *
        adaptive_stokes_cached_rotated(query,
                                       nPatchData,
                                       quad_rule,
                                       rotator,
                                       n,
                                       0,
                                       0,
                                       quad_coarse,
                                       quad_tol);
    }
    // Otherwise, do some refinement
    else
    {
      quad += 0.25 * M_1_PI * quad_coarse;
    }
  }

  return quad;
}

template <typename T>
double stokes_winding_number_cached_rotated_adaptless(
  const Point<T, 3>& query,
  const NURBSPatchData<T>& nPatchData,
  const axom::numerics::Matrix<T>& rotator,
  const int quad_npts,
  const double quad_tol)
{
  // Generate the quadrature rules in parameter space
  static mfem::IntegrationRules my_IntRules(0, mfem::Quadrature1D::GaussLegendre);
  const mfem::IntegrationRule& quad_rule =
    my_IntRules.Get(mfem::Geometry::SEGMENT, 2 * quad_npts - 1);

  double quad = 0;
  for(int n = 0; n < nPatchData.patch.getNumTrimmingCurves(); ++n)
  {
    // Get the quadrature points for the curve without any refinement
    auto trimming_curve_data = nPatchData.getQuadratureData(n, quad_rule, 0, 0);
    double quad_coarse =
      single_stokes_cached_rotated(query, quad_rule, rotator, trimming_curve_data);

    quad += 0.25 * M_1_PI * quad_coarse;
  }

  return quad;
}

template <typename T>
double adaptive_stokes_cached_rotated(const Point<T, 3>& query,
                                      const NURBSPatchData<T>& nPatchData,
                                      const mfem::IntegrationRule& quad_rule,
                                      const axom::numerics::Matrix<T>& rotator,
                                      const int curve_index,
                                      const int refinement_level,
                                      const int refinement_index,
                                      const double quad_coarse,
                                      const double quad_tol)
{
  auto trimming_curve_data_1 =
    nPatchData.getQuadratureData(curve_index,
                                 quad_rule,
                                 refinement_level + 1,
                                 2 * refinement_index);
  auto trimming_curve_data_2 =
    nPatchData.getQuadratureData(curve_index,
                                 quad_rule,
                                 refinement_level + 1,
                                 2 * refinement_index + 1);

  double quad_fine_1 =
    single_stokes_cached_rotated(query, quad_rule, rotator, trimming_curve_data_1);
  double quad_fine_2 =
    single_stokes_cached_rotated(query, quad_rule, rotator, trimming_curve_data_2);

  if(refinement_level >= 15 ||
     axom::utilities::isNearlyEqualRelative(quad_fine_1 + quad_fine_2,
                                            quad_coarse,
                                            quad_tol,
                                            1e-10))
  {
    return quad_fine_1 + quad_fine_2;
  }
  else
  {
    // If we're not near the axis, we can trust the original value
    if(FORCE_ADAPTIVE || isNearAxisBoxRotated(query, trimming_curve_data_1.bbox, 5.0, rotator))
      quad_fine_1 = adaptive_stokes_cached_rotated(query,
                                                   nPatchData,
                                                   quad_rule,
                                                   rotator,
                                                   curve_index,
                                                   refinement_level + 1,
                                                   2 * refinement_index,
                                                   quad_fine_1,
                                                   quad_tol);

    if(FORCE_ADAPTIVE || isNearAxisBoxRotated(query, trimming_curve_data_2.bbox, 5.0, rotator))
      quad_fine_2 = adaptive_stokes_cached_rotated(query,
                                                   nPatchData,
                                                   quad_rule,
                                                   rotator,
                                                   curve_index,
                                                   refinement_level + 1,
                                                   2 * refinement_index + 1,
                                                   quad_fine_2,
                                                   quad_tol);

    return quad_fine_1 + quad_fine_2;
  }
}

template <typename T>
double single_stokes_cached_rotated(
  const Point<T, 3>& query,
  const mfem::IntegrationRule& quad_rule,
  const axom::numerics::Matrix<T>& rotator,
  const TrimmingCurveQuadratureData<T>& trimming_curve_data)
{
  // Do this without refinement
  double this_quad = 0;
  for(int q = 0; q < quad_rule.GetNPoints(); ++q)
  {
    const Vector<T, 3> node(
      query,
      rotate_point(rotator, query, trimming_curve_data.quadrature_points[q].first));
    const Vector<T, 3> node_dt(
      rotate_vector(rotator,
                    Point<T, 3>({0.0, 0.0, 0.0}),
                    trimming_curve_data.quadrature_points[q].second));
    const double node_norm = node.norm();

    const double quad_weight =
      quad_rule.IntPoint(q).weight * trimming_curve_data.span_length;

    // This is the z-aligned axis
    this_quad += quad_weight *
      (node[1] * node[2] * node_dt[0] - node[0] * node[2] * node_dt[1]) /
      (node[0] * node[0] + node[1] * node[1]) / node_norm;
  }

  return this_quad;
}

// --=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--

template <typename T>
double surface_winding_number(const Point<T, 3>& query,
                              const NURBSPatch<T, 3>& patch,
                              int npts)
{
  // Generate the quadrature rules in parameter space
  static mfem::IntegrationRules my_IntRules(0, mfem::Quadrature1D::GaussLegendre);

  const mfem::IntegrationRule& quad_rule =
    my_IntRules.Get(mfem::Geometry::SEGMENT, 2 * npts - 1);

  double quadrature = 0.0;
  for(int qu = 0; qu < quad_rule.GetNPoints(); ++qu)
  {
    for(int qv = 0; qv < quad_rule.GetNPoints(); ++qv)
    {
      Vector3D node(
        query,
        patch.evaluate(quad_rule.IntPoint(qu).x, quad_rule.IntPoint(qv).x));

      // Compute the normal vector
      Vector3D normal =
        patch.normal(quad_rule.IntPoint(qu).x, quad_rule.IntPoint(qv).x);

      quadrature += quad_rule.IntPoint(qu).weight *
        quad_rule.IntPoint(qv).weight * Vector3D::dot_product(node, normal) /
        std::pow(node.norm(), 3);
    }
  }

  return 0.25 * M_1_PI * quadrature;
}

template <typename T>
double surface_winding_number(const Point<T, 3>& query,
                              const BezierPatch<T, 3>& patch,
                              int npts,
                              bool needs_adapt = true)
{
  // Generate the quadrature rules in parameter space
  static mfem::IntegrationRules my_IntRules(0, mfem::Quadrature1D::GaussLegendre);

  const mfem::IntegrationRule& quad_rule =
    my_IntRules.Get(mfem::Geometry::SEGMENT, 2 * npts - 1);

  double quadrature = 0.0;
  for(int qu = 0; qu < quad_rule.GetNPoints(); ++qu)
  {
    for(int qv = 0; qv < quad_rule.GetNPoints(); ++qv)
    {
      Vector3D node(
        query,
        patch.evaluate(quad_rule.IntPoint(qu).x, quad_rule.IntPoint(qv).x));

      // Compute the normal vector
      Vector3D normal =
        patch.normal(quad_rule.IntPoint(qu).x, quad_rule.IntPoint(qv).x);

      quadrature += quad_rule.IntPoint(qu).weight *
        quad_rule.IntPoint(qv).weight * Vector3D::dot_product(node, normal) /
        std::pow(node.norm(), 3);
    }
  }

  return 0.25 * M_1_PI * quadrature;
}

template <typename T>
double surface_winding_number_adaptive(const Point<T, 3>& query,
                                       const BezierPatch<T, 3>& patch,
                                       const mfem::IntegrationRule& quad_rule,
                                       const double quad_coarse,
                                       const double quad_tol,
                                       const int depth = 1)
{
  // Split the patch, do the quadrature over all four components
  BezierPatch<T, 3> subpatches[4];
  patch.split(0.5, 0.5, subpatches[0], subpatches[1], subpatches[2], subpatches[3]);

  double quad_fine[4] = {0.0, 0.0, 0.0, 0.0};
  for(int i = 0; i < 4; ++i)
  {
    for(int qu = 0; qu < quad_rule.GetNPoints(); ++qu)
    {
      for(int qv = 0; qv < quad_rule.GetNPoints(); ++qv)
      {
        Vector3D node(query,
                      subpatches[i].evaluate(quad_rule.IntPoint(qu).x,
                                             quad_rule.IntPoint(qv).x));

        // Compute the normal vector
        Vector3D normal = subpatches[i].normal(quad_rule.IntPoint(qu).x,
                                               quad_rule.IntPoint(qv).x);

        quad_fine[i] += quad_rule.IntPoint(qu).weight *
          quad_rule.IntPoint(qv).weight * Vector3D::dot_product(node, normal) /
          std::pow(node.norm(), 3);
      }
    }
  }

  constexpr int MAX_DEPTH = 12;
  if(depth >= MAX_DEPTH ||
     axom::utilities::isNearlyEqualRelative(
       quad_fine[0] + quad_fine[1] + quad_fine[2] + quad_fine[3],
       quad_coarse,
       quad_tol,
       1e-10))
  {
    return 0.25 * M_1_PI *
      (quad_fine[0] + quad_fine[1] + quad_fine[2] + quad_fine[3]);
  }
  else
  {
    return surface_winding_number_adaptive(query,
                                           subpatches[0],
                                           quad_rule,
                                           quad_fine[0],
                                           quad_tol,
                                           depth + 1) +
      surface_winding_number_adaptive(query,
                                      subpatches[1],
                                      quad_rule,
                                      quad_fine[1],
                                      quad_tol,
                                      depth + 1) +
      surface_winding_number_adaptive(query,
                                      subpatches[2],
                                      quad_rule,
                                      quad_fine[2],
                                      quad_tol,
                                      depth + 1) +
      surface_winding_number_adaptive(query,
                                      subpatches[3],
                                      quad_rule,
                                      quad_fine[3],
                                      quad_tol,
                                      depth + 1);
  }
}
#endif

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif
