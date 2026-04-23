// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
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

#include "axom/core/numerics/transforms.hpp"
#include "axom/core/numerics/quadrature.hpp"

#include "axom/primal/operators/detail/winding_number_3d_memoization.hpp"

// C++ includes
#include <cstdint>
#include <functional>
#include <limits>
#include <math.h>
#include <optional>

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
template <typename NURBSType, typename T>
void degenerate_surface_processing(const NURBSType& nurbs,
                                   const axom::Array<T>& up,
                                   const axom::Array<T>& vp,
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

  NURBSPatch<T, 3> dummy_patch(nurbs.getControlPoints(),
                               nurbs.getWeights(),
                               nurbs.getKnots_u(),
                               nurbs.getKnots_v());
  dummy_patch.setTrimmingCurves(nurbs.getTrimmingCurves());

  // Indicates a u isocurve
  if(var_u < var_v)
  {
    if(mean_u - clip_radius > nurbs.getMinKnot_u())
    {
      dummy_patch.split_u(mean_u - clip_radius, out_patch1, dummy_patch);
    }
    else
    {
      out_patch1 = NURBSPatch<T, 3>();
    }

    if(mean_u + clip_radius < nurbs.getMaxKnot_u())
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
    if(mean_v - clip_radius > nurbs.getMinKnot_v())
    {
      dummy_patch.split_v(mean_v - clip_radius, out_patch1, dummy_patch);
    }
    else
    {
      out_patch1 = NURBSPatch<T, 3>();
    }

    if(mean_v + clip_radius < nurbs.getMaxKnot_v())
    {
      dummy_patch.split_v(mean_v + clip_radius, dummy_patch, out_patch2);
    }
    else
    {
      out_patch2 = NURBSPatch<T, 3>();
    }
  }
}

//! \brief Rotate a point around another point according to the provided rotation matrix
template <typename T>
Point<T, 3> rotate_point(const numerics::Matrix<T>& matx,
                         const Point<T, 3>& center,
                         const Point<T, 3>& input)
{
  Vector<T, 3> shifted(center, input);
  Vector<T, 3> rotated;
  numerics::matrix_vector_multiply(matx, shifted.data(), rotated.data());
  return Point<T, 3>({rotated[0] + center[0], rotated[1] + center[1], rotated[2] + center[2]});
}

//! \brief Rotate a vector around the origin according to the provided rotation matrix
template <typename T>
Vector<T, 3> rotate_vector_origin(const numerics::Matrix<T>& matx, const Vector<T, 3>& input)
{
  Vector<T, 3> shifted {input[0], input[1], input[2]};
  Vector<T, 3> rotated;
  numerics::matrix_vector_multiply(matx, shifted.data(), rotated.data());
  return rotated;
}

/*!
 * \brief Adaptively evaluate the integral of the "anti-curl" of the GWN integrand
 *
 * \param [in] query The query point
 * \param [in] nurbs The NURBSPatchGWNCache object containing the trimming curves
 * \param [in] curve_index The curve on the patch which we want to integrate
 * \param [in] quad_npts The number of quadrature points at each level
 * \param [in] refinement_level The current subdivision levels 
 * \param [in] refinement_index Which subdivision in the level
 * \param [in] ax The axis (relative to query) denoting which anti-curl we use
 * \param [in] rotator This is the rotation matrix to use if ax == DiscontinuityAxis::rotated
 * \param [in] quad_coarse The approximate integral of the curve, 
 *              which should match the sum of the integral over each half
 * \param [in] quad_tol The maximum relative error allowed in each quadrature
 * 
 * \return The value of the integral
 */
template <typename T>
double stokes_gwn_adaptive(const Point<T, 3>& query,
                           const NURBSPatchGWNCache<T>& nurbs,
                           const int curve_index,
                           const int quad_npts,
                           const int refinement_level,
                           const int refinement_index,
                           const DiscontinuityAxis ax,
                           const numerics::Matrix<T>& rotator,
                           const double quad_coarse,
                           const double quad_tol)
{
  const auto& trimming_curve_data_1 = nurbs.getTrimmingCurveQuadratureData(curve_index,
                                                                           quad_npts,
                                                                           refinement_level + 1,
                                                                           2 * refinement_index);
  const auto& trimming_curve_data_2 = nurbs.getTrimmingCurveQuadratureData(curve_index,
                                                                           quad_npts,
                                                                           refinement_level + 1,
                                                                           2 * refinement_index + 1);

  double quad_fine_1 = stokes_gwn_component(query, ax, rotator, trimming_curve_data_1);
  double quad_fine_2 = stokes_gwn_component(query, ax, rotator, trimming_curve_data_2);

  if(refinement_level >= 25 ||
     axom::utilities::isNearlyEqualRelative(quad_fine_1 + quad_fine_2, quad_coarse, quad_tol, 1e-10))
  {
    return quad_fine_1 + quad_fine_2;
  }

  quad_fine_1 = stokes_gwn_adaptive(query,
                                    nurbs,
                                    curve_index,
                                    quad_npts,
                                    refinement_level + 1,
                                    2 * refinement_index,
                                    ax,
                                    rotator,
                                    quad_fine_1,
                                    quad_tol);

  quad_fine_2 = stokes_gwn_adaptive(query,
                                    nurbs,
                                    curve_index,
                                    quad_npts,
                                    refinement_level + 1,
                                    2 * refinement_index + 1,
                                    ax,
                                    rotator,
                                    quad_fine_2,
                                    quad_tol);

  return quad_fine_1 + quad_fine_2;
}

/*!
 * \brief Adaptively evaluate the integral of the "anti-curl" of the GWN integrand
 *
 * \param [in] query The query point
 * \param [in] nurbs The NURBSPatch object contianing the trimming curves
 * \param [in] curve_index The curve on the patch which we want to integrate
 * \param [in] quad_npts The number of quadrature points at each level
 * \param [in] refinement_level The current subdivision levels 
 * \param [in] refinement_index Which subdivision in the level
 * \param [in] ax The axis (relative to query) denoting which anti-curl we use
 * \param [in] rotator This is the rotation matrix to use if ax == DiscontinuityAxis::rotated
 * \param [in] quad_coarse The approximate integral of the curve, 
 *              which should match the sum of the integral over each half
 * \param [in] quad_tol The maximum relative error allowed in each quadrature
 * 
 * \return The value of the integral
 */
template <typename T>
double stokes_gwn_adaptive(const Point<T, 3>& query,
                           const NURBSPatch<T, 3>& nurbs,
                           const int curve_index,
                           const int quad_npts,
                           const int refinement_level,
                           const int refinement_index,
                           const DiscontinuityAxis ax,
                           const numerics::Matrix<T>& rotator,
                           const double quad_coarse,
                           const double quad_tol)
{
  auto trimming_curve_data_1 = TrimmingCurveQuadratureData<T>(nurbs,
                                                              curve_index,
                                                              quad_npts,
                                                              refinement_level + 1,
                                                              2 * refinement_index);
  auto trimming_curve_data_2 = TrimmingCurveQuadratureData<T>(nurbs,
                                                              curve_index,
                                                              quad_npts,
                                                              refinement_level + 1,
                                                              2 * refinement_index + 1);

  double quad_fine_1 = stokes_gwn_component(query, ax, rotator, trimming_curve_data_1);
  double quad_fine_2 = stokes_gwn_component(query, ax, rotator, trimming_curve_data_2);

  if(refinement_level >= 25 ||
     axom::utilities::isNearlyEqualRelative(quad_fine_1 + quad_fine_2, quad_coarse, quad_tol, 1e-10))
  {
    return quad_fine_1 + quad_fine_2;
  }

  quad_fine_1 = stokes_gwn_adaptive(query,
                                    nurbs,
                                    curve_index,
                                    quad_npts,
                                    refinement_level + 1,
                                    2 * refinement_index,
                                    ax,
                                    rotator,
                                    quad_fine_1,
                                    quad_tol);

  quad_fine_2 = stokes_gwn_adaptive(query,
                                    nurbs,
                                    curve_index,
                                    quad_npts,
                                    refinement_level + 1,
                                    2 * refinement_index + 1,
                                    ax,
                                    rotator,
                                    quad_fine_2,
                                    quad_tol);

  return quad_fine_1 + quad_fine_2;
}

/*!
 * \brief Evaluates the integral of the "anti-curl" of the GWN integrand on one curve
 *
 * \param [in] query The query point
 * \param [in] ax The axis (relative to query) denoting which anti-curl we use
 * \param [in] rotator This is the rotation matrix to use if ax == DiscontinuityAxis::rotated
 * \param [in] trimming_curve_data The struct with all quadrature info (nodes and tangents)
 * 
 * \note This is only meant to be used for `stokes_gwn_evaluate()`,
 *  and the result does not make sense outside of that context.
 * 
 * \return The value of the integral
 */
template <typename T>
double stokes_gwn_component(const Point<T, 3>& query,
                            const DiscontinuityAxis ax,
                            const numerics::Matrix<T>& rotator,
                            const TrimmingCurveQuadratureData<T>& trimming_curve_data)
{
  using VectorType = Vector<T, 3>;

  const int num_pts = trimming_curve_data.getNumPoints();
  const auto quad_points = trimming_curve_data.getQuadraturePoints();
  const auto quad_tangents = trimming_curve_data.getQuadratureTangents();
  const auto quad_weights = trimming_curve_data.getQuadratureWeights();

  // Shared loop structure for the different discontinuity axes:
  //  - build node vector (quadrature point relative to query, optionally rotated)
  //  - build tangent vector (optionally rotated)
  //  - accumulate weight * numerator(node, dt) / (denom(node) * ||node||)
  auto accumulate = [&](auto&& make_node, auto&& make_dt, auto&& numer_fn, auto&& denom_fn) -> double {
    double this_quad = 0.0;
    for(int q = 0; q < num_pts; ++q)
    {
      const auto& node = make_node(q);
      const auto& node_dt = make_dt(q);
      const double inv = 1.0 / (denom_fn(node) * node.norm());
      this_quad += quad_weights[q] * numer_fn(node, node_dt) * inv;
    }
    return this_quad;
  };

  auto numer_x = [](const VectorType& node, const VectorType& dt) -> double {
    return node[0] * (node[2] * dt[1] - node[1] * dt[2]);
  };
  auto numer_y = [](const VectorType& node, const VectorType& dt) -> double {
    return node[1] * (node[0] * dt[2] - node[2] * dt[0]);
  };
  auto numer_z = [](const VectorType& node, const VectorType& dt) -> double {
    return node[2] * (node[1] * dt[0] - node[0] * dt[1]);
  };

  auto denom_x = [](const VectorType& node) -> double {
    return (node[1] * node[1]) + (node[2] * node[2]);
  };
  auto denom_y = [](const VectorType& node) -> double {
    return (node[0] * node[0]) + (node[2] * node[2]);
  };
  auto denom_z = [](const VectorType& node) -> double {
    return (node[0] * node[0]) + (node[1] * node[1]);
  };

  // Non-rotated paths share the same node/dt construction.
  auto make_node_unrot = [&](int q) -> VectorType { return quad_points[q] - query; };
  auto make_dt_unrot = [&](int q) -> const VectorType& { return quad_tangents[q]; };

  // The axis is constant for the full loop; split into specialized loops to avoid
  // per-iteration switching/branching overhead.
  if(ax == DiscontinuityAxis::rotated)
  {
    // Rotate query-relative vectors directly; avoid the extra translate in
    // rotate_point(..., query, P) - query.
    const T* r = rotator.data();
    const double r00 = r[0], r01 = r[3], r02 = r[6];
    const double r10 = r[1], r11 = r[4], r12 = r[7];
    const double r20 = r[2], r21 = r[5], r22 = r[8];

    auto make_node_rot = [&](int q) -> VectorType {
      const auto& qp = quad_points[q];
      const double dx = qp[0] - query[0];
      const double dy = qp[1] - query[1];
      const double dz = qp[2] - query[2];

      return VectorType {r00 * dx + r01 * dy + r02 * dz,
                         r10 * dx + r11 * dy + r12 * dz,
                         r20 * dx + r21 * dy + r22 * dz};
    };

    auto make_dt_rot = [&](int q) -> VectorType {
      const auto& qt = quad_tangents[q];

      return VectorType {r00 * qt[0] + r01 * qt[1] + r02 * qt[2],
                         r10 * qt[0] + r11 * qt[1] + r12 * qt[2],
                         r20 * qt[0] + r21 * qt[1] + r22 * qt[2]};
    };

    return accumulate(make_node_rot, make_dt_rot, numer_z, denom_z);
  }

  if(ax == DiscontinuityAxis::x)
  {
    return accumulate(make_node_unrot, make_dt_unrot, numer_x, denom_x);
  }

  if(ax == DiscontinuityAxis::y)
  {
    return accumulate(make_node_unrot, make_dt_unrot, numer_y, denom_y);
  }

  // Default: z-axis field.
  return accumulate(make_node_unrot, make_dt_unrot, numer_z, denom_z);
}

/*!
 * \brief Evaluates the integral of the "anti-curl" of the GWN integrand
 *        (via Stokes' theorem) at a point wrt to the trimming curves of a surface
 *
 * \param [in] query The query point
 * \param [in] nurbs The NURBSPatchGWNCache object contianing the trimming curves
 * \param [in] quad_npts The number of points used in each Gaussian quadrature
 * \param [in] ax The axis (relative to query) denoting which anti-curl we use
 * \param [in] rotator This is the rotation matrix to use if ax == DiscontinuityAxis::rotated
 * \param [in] quad_tol The maximum relative error allowed in each quadrature
 * 
 * Applies a non-adaptive quadrature to the trimming curves of a NURBS patch using one of three possible
 * "anti-curl" vector fields, the curl of each of which is equal to <x, y, z>/||x||^3.
 *
 * \return The value of the integral
 */
template <typename T>
double stokes_gwn_evaluate(const Point<T, 3>& query,
                           const NURBSPatchGWNCache<T>& nurbs,
                           const int quad_npts,
                           DiscontinuityAxis ax,
                           const numerics::Matrix<T>& rotator,
                           const double quad_tol)
{
  constexpr double gwn_modulo = 0.25 * M_1_PI;

  // Can't rotate the patch as pre-processing if working with cached data
  double quad = 0;
  for(int n = 0; n < nurbs.getNumTrimmingCurves(); ++n)
  {
    // Get the quadrature points for the curve on the patch without any refinement
    const auto& trimming_curve_data = nurbs.getTrimmingCurveQuadratureData(n, quad_npts, 0, 0);
    const double quad_coarse = stokes_gwn_component(query, ax, rotator, trimming_curve_data);

    quad += gwn_modulo *
      stokes_gwn_adaptive(query, nurbs, n, quad_npts, 0, 0, ax, rotator, quad_coarse, quad_tol);
  }

  return quad;
}

/*!
 * \brief Overload of stokes_gwn_evaluate if NURBSPatch is not memoized
 *
 * \param [in] query The query point
 * \param [in] nurbs The NURBSPatch object containing the trimming curves
 * \param [in] quad_npts The number of points used in each Gaussian quadrature
 * \param [in] ax The axis (relative to query) denoting which anti-curl we use
 * \param [in] rotator This is the rotation matrix to use if ax == DiscontinuityAxis::rotated
 * \param [in] quad_tol The maximum relative error allowed in each quadrature
 * 
 * If quadrature nodes are not memoized, then we can rotate the patch by its control points
 *  instead of rotating individual quadrature nodes
 *
 * \return The value of the integral
 */
template <typename T>
double stokes_gwn_evaluate(const Point<T, 3>& query,
                           const NURBSPatch<T, 3>& nurbs,
                           const int quad_npts,
                           DiscontinuityAxis ax,
                           const numerics::Matrix<T>& rotator,
                           const double quad_tol)
{
  // Only copy/rotate the patch when needed; otherwise work directly with `nurbs`.
  const NURBSPatch<T, 3>* nurbs_eval = &nurbs;
  std::optional<NURBSPatch<T, 3>> rotated_patch;

  if(ax == DiscontinuityAxis::rotated)
  {
    rotated_patch.emplace(nurbs);
    auto& rotated = *rotated_patch;

    const auto patch_shape = nurbs.getControlPoints().shape();
    for(int i = 0; i < patch_shape[0]; ++i)
    {
      for(int j = 0; j < patch_shape[1]; ++j)
      {
        rotated(i, j) = rotate_point(rotator, query, nurbs(i, j));
      }
    }

    nurbs_eval = &rotated;

    // Change the rotation axis so we don't rotate a second time
    ax = DiscontinuityAxis::z;
  }

  constexpr double gwn_modulo = 0.25 * M_1_PI;

  double quad = 0;
  for(int n = 0; n < nurbs_eval->getNumTrimmingCurves(); ++n)
  {
    // Get the quadrature points for the curve on the *rotated* patch without any refinement
    const auto trimming_curve_data = TrimmingCurveQuadratureData<T>(*nurbs_eval, n, quad_npts, 0, 0);
    const double quad_coarse = stokes_gwn_component(query, ax, rotator, trimming_curve_data);

    quad += gwn_modulo *
      stokes_gwn_adaptive(query, *nurbs_eval, n, quad_npts, 0, 0, ax, rotator, quad_coarse, quad_tol);
  }

  return quad;
}

/*!
 * \brief Computes the GWN for a 3D point wrt a 3D NURBS patch with precomputed data
 *
 * \tparam NURBSType The memoized (NURBSPatchGWNCache) or un-memoized (NURBSPatch) surface type
 * \param [in] query The query point to test
 * \param [in] nurbs The NURBS patch object with precomputed data
 * \param [in] cast_direction The vector which defines the discrete correction term
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] ls_tol The tolerance for the line-surface intersection routine
 * \param [in] quad_tol The maximum relative error allowed in the quadrature
 * \param [in] disk_size The size of extracted disks as a percent of parameter bbox diagonal
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * Computes the generalized winding number for a NURBS patch using Stokes theorem,
 *  along with a correction term determined by a line-surface intersection test
 *
 * \pre Assumes that the NURBS patch is trimmed, and has been slightly extended in 
 *       parameter space so that trimming curves arent't on the boundary of the untrimmed patch
 * \return The GWN.
 */
template <typename NURBSType, typename T>
double nurbs_winding_number(const Point<T, 3>& query,
                            const NURBSType& nurbs,
                            const Vector<T, 3>& cast_direction,
                            const double edge_tol = 1e-8,
                            const double ls_tol = 1e-8,
                            const double quad_tol = 1e-8,
                            const double disk_size = 0.01,
                            const double EPS = 1e-8)
{
  // Skip processing of degenerate surfaces
  if(nurbs.getNumControlPoints_u() <= 1 || nurbs.getNumControlPoints_v() <= 1)
  {
    return 0.0;
  }

  // Also skip processing of surfaces with zero trimming curves
  if(nurbs.getNumTrimmingCurves() == 0)
  {
    return 0.0;
  }

  const double edge_tol_sq = edge_tol * edge_tol;

  // Fix the number of quadrature points arbitrarily
  constexpr int k_quad_npts = 15;
  constexpr double k_bbox_expand_frac = 0.01;
  constexpr double k_degenerate_clip_frac = 0.01;
  constexpr double k_disk_shrink_on_surface_frac = 0.01;
  constexpr int k_max_noncoincident_intersections = 5;
  constexpr int k_max_recast_attempts = 128;

  // Store the winding number
  double the_gwn = 0.0;

  /*
   * To use Stokes theorem, we need to identify either a line containing the
   * query that does not intersect the surface, or one that intersects the *interior*
   * of the surface at known locations.
   */

  auto hash_combine = [](std::uint64_t seed, std::uint64_t value) -> std::uint64_t {
    return seed ^ (value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2));
  };

  auto splitmix64 = [](std::uint64_t value) -> std::uint64_t {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
  };

  auto deterministic_unit = [&](std::uint64_t seed) -> Vector<T, 3> {
    const std::uint64_t theta_bits = splitmix64(seed);
    const std::uint64_t z_bits = splitmix64(theta_bits);
    constexpr double inv_max = 1.0 / static_cast<double>(std::numeric_limits<std::uint64_t>::max());
    const double theta = 2.0 * M_PI * static_cast<double>(theta_bits) * inv_max;
    const double z = 2.0 * static_cast<double>(z_bits) * inv_max - 1.0;
    const double radial_sq = 1.0 - z * z;
    const double radial = std::sqrt(radial_sq > 0.0 ? radial_sq : 0.0);
    return Vector<T, 3> {radial * std::cos(theta), radial * std::sin(theta), z};
  };

  auto recast_seed = [&]() -> std::uint64_t {
    std::uint64_t seed = 0;
    const std::hash<double> hasher {};
    for(int dim = 0; dim < 3; ++dim)
    {
      seed = hash_combine(seed, static_cast<std::uint64_t>(hasher(static_cast<double>(query[dim]))));
      seed =
        hash_combine(seed,
                     static_cast<std::uint64_t>(hasher(static_cast<double>(cast_direction[dim]))));
    }
    return seed;
  };
  const std::uint64_t base_seed = recast_seed();

  // Recasting is implemented as a loop (not recursion) to avoid stack overflows
  // in OpenMP thread stacks on pathological intersections.
  Vector<T, 3> cast_direction_local = cast_direction;
  for(int recast_attempt = 0; recast_attempt < k_max_recast_attempts; ++recast_attempt)
  {
    const bool can_recast = (recast_attempt + 1 < k_max_recast_attempts);

    auto request_recast = [&]() -> bool {
      if(can_recast)
      {
        cast_direction_local =
          deterministic_unit(hash_combine(base_seed, static_cast<std::uint64_t>(recast_attempt + 1)));
        return true;
      }
      return false;
    };

    the_gwn = 0.0;

    // Rotation matrix for the patch
    numerics::Matrix<T> rotator;

    // Lazily allocate space for the patch which contains all surface boundaries,
    //  and any extra trimming curves added by disk extraction.
    // Note: For most query points (relative to a patch bounding box), we exit early
    //       and can avoid making a deep copy of the surface and trimming curves.
    std::optional<NURBSPatch<T, 3>> nurbs_modified;

    // Define vector fields whose curl gives us the winding number
    DiscontinuityAxis field_direction = DiscontinuityAxis::rotated;
    bool extraTrimming = false;

    // Generate slightly expanded bounding boxes
    auto bBox = nurbs.boundingBox();
    auto oBox = nurbs.orientedBoundingBox();

    auto patch_diameter = bBox.range().norm();

    bBox.expand(k_bbox_expand_frac * patch_diameter);
    oBox.expand(k_bbox_expand_frac * patch_diameter);

    // Case 1: Exterior without rotations
    if(!bBox.contains(query))
    {
      double best_dist = -1.0;
      auto consider_axis = [&](double dist, DiscontinuityAxis axis) {
        if(dist > best_dist)
        {
          best_dist = dist;
          field_direction = axis;
        }
      };

      if(query[0] <= bBox.getMin()[0])
      {
        consider_axis(bBox.getMin()[0] - query[0], DiscontinuityAxis::y);
      }
      else if(query[0] >= bBox.getMax()[0])
      {
        consider_axis(query[0] - bBox.getMax()[0], DiscontinuityAxis::y);
      }

      if(query[1] <= bBox.getMin()[1])
      {
        consider_axis(bBox.getMin()[1] - query[1], DiscontinuityAxis::z);
      }
      else if(query[1] >= bBox.getMax()[1])
      {
        consider_axis(query[1] - bBox.getMax()[1], DiscontinuityAxis::z);
      }

      if(query[2] <= bBox.getMin()[2])
      {
        consider_axis(bBox.getMin()[2] - query[2], DiscontinuityAxis::y);
      }
      else if(query[2] >= bBox.getMax()[2])
      {
        consider_axis(query[2] - bBox.getMax()[2], DiscontinuityAxis::x);
      }
    }
    // Case 1.5: Exterior with rotation
    else if(!oBox.contains(query))
    {
      // Rotate the patch until the OBB is not directly above/below the query point.
      field_direction = DiscontinuityAxis::rotated;

      // Find vector from query to the bounding box
      const Point<T, 3> closest = closest_point(query, oBox);
      const Vector<T, 3> v0 = Vector<T, 3>(query, closest).unitVector();

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
      const double ang = (v0[2] < 0 ? 1.0 : -1.0) *
        acos(axom::utilities::clampVal(
          -(v0[0] * v1[1] - v0[1] * v1[0]) / sqrt(v1[0] * v1[0] + v1[1] * v1[1]),
          -1.0,
          1.0));

      rotator = numerics::transforms::axisRotation(ang, v1[0], v1[1], v1[2]);
    }
    else
    {
      field_direction = DiscontinuityAxis::rotated;
      const Line<T, 3> discontinuity_axis(query, cast_direction_local);
      bool should_recast = false;

      // Tolerance for what counts as "close to a boundary" in parameter space
      T disk_radius = disk_size * nurbs.getParameterSpaceDiagonal();

      // Compute intersections with the *untrimmed and extrapolated* patch
      axom::Array<T> up, vp, tp;
      const bool isHalfOpen = false, countUntrimmed = true;

      bool success = true;
      nurbs_modified.emplace(nurbs.getControlPoints(),
                             nurbs.getWeights(),
                             nurbs.getKnots_u(),
                             nurbs.getKnots_v());
      nurbs_modified->setTrimmingCurves(nurbs.getTrimmingCurves());

      intersect(discontinuity_axis,
                *nurbs_modified,
                tp,
                up,
                vp,
                ls_tol,
                EPS,
                countUntrimmed,
                isHalfOpen,
                success);

      if(!success)
      {
        // Look at the intersection points
        int num_noncoincident = 0;
        for(int i = 0; i < tp.size(); ++i)
        {
          const Point<T, 3> the_point(nurbs_modified->evaluate(up[i], vp[i]));
          // If any of the intersection points are coincident with the surface,
          //  then attempt to clip out all degenerate intersections, and retry
          if(squared_distance(query, the_point) <= edge_tol_sq)
          {
            NURBSPatch<T, 3> clipped_patch1, clipped_patch2;

            degenerate_surface_processing(nurbs,
                                          up,
                                          vp,
                                          k_degenerate_clip_frac * disk_radius,
                                          clipped_patch1,
                                          clipped_patch2);

            return nurbs_winding_number(query,
                                        clipped_patch1,
                                        cast_direction_local,
                                        edge_tol,
                                        ls_tol,
                                        quad_tol * 1e-5,
                                        disk_size,
                                        EPS) +
              nurbs_winding_number(query,
                                   clipped_patch2,
                                   cast_direction_local,
                                   edge_tol,
                                   ls_tol,
                                   quad_tol * 1e-5,
                                   disk_size,
                                   EPS);
          }
          else
          {
            num_noncoincident++;
          }

          // If more than 5 (arbitrary) are *not* coincident with the surface,
          //  re-cast and try again. This is to avoid cases where the point *is*
          //  coincident with the surface, but the first recorded point of
          //  intersection is not after multiple re-casts.
          if(num_noncoincident > k_max_noncoincident_intersections && request_recast())
          {
            should_recast = true;
            break;
          }
        }
      }

      if(should_recast)
      {
        continue;
      }

      // If no intersections are recorded, then nothing extra to account for

      // Otherwise, account for each discontinuity analytically,
      //  or recursively through disk subdivision
      for(int i = 0; i < up.size(); ++i)
      {
        // Compute the intersection point on the surface
        const Point<T, 3> the_point(nurbs_modified->evaluate(up[i], vp[i]));
        const Vector<T, 3> the_normal = nurbs_modified->normal(up[i], vp[i]);

        // Check for bad intersections, i.e.,
        //  > There normal is poorly defined (cusp)
        //  > The normal is tangent to the axis of discontinuity
        const bool bad_intersection = axom::utilities::isNearlyEqual(the_normal.norm(), 0.0, EPS) ||
          axom::utilities::isNearlyEqual(the_normal.unitVector().dot(cast_direction_local), 0.0, EPS);

        const bool isOnSurface = squared_distance(query, the_point) <= edge_tol_sq;

        if(bad_intersection && !isOnSurface && request_recast())
        {
          // If a non-coincident ray intersects the surface at a tangent/cusp,
          //  can recast and try again with the memoized patch
          should_recast = true;
          break;
        }

        if(isOnSurface)
        {
          // If the query point is on the surface, then shrink the disk
          //  to ensure its winding number is known to be near-zero
          disk_radius = k_disk_shrink_on_surface_frac * disk_radius;
        }

        // Consider a disk around the intersection point via NURBSPatch::diskSplit.
        //   If the disk intersects any trimming curves, need to do disk subdivision.
        //   If not, we can compute the winding number without changing the trimming curvse
        const bool ignoreInteriorDisk = true;
        bool isDiskInside, isDiskOutside;
        NURBSPatch<T, 3> the_disk;

        nurbs_modified->diskSplit(up[i],
                                  vp[i],
                                  disk_radius,
                                  the_disk,
                                  *nurbs_modified,
                                  isDiskInside,
                                  isDiskOutside,
                                  ignoreInteriorDisk,
                                  disk_radius);

        extraTrimming = extraTrimming || (!isDiskInside && !isDiskOutside) ||
          (isDiskInside && !ignoreInteriorDisk);

        if(isOnSurface)
        {
          // If the query point is on the surface, the contribution of the disk is near-zero
          //  and we only needed to puncture the larger surface to proceed
          continue;
        }
        else if(!isDiskInside && !isDiskOutside)
        {
          // If the disk overlapped with the trimming curves, evaluate the winding number for the disk
          //  with a cast ray that is mostly in the direction of the normal (assuming it's non-zero)
          const auto perturbation =
            deterministic_unit(hash_combine(base_seed, static_cast<std::uint64_t>(i + 1)));
          Vector<T, 3> new_cast_direction = the_disk.normal(up[i], vp[i]);
          new_cast_direction = (new_cast_direction.norm() < EPS)
            ? perturbation
            : (new_cast_direction.unitVector() + 0.1 * perturbation).unitVector();

          the_gwn += nurbs_winding_number(query,
                                          the_disk,
                                          new_cast_direction,
                                          edge_tol,
                                          ls_tol,
                                          quad_tol,
                                          disk_size,
                                          EPS);
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
          const auto the_direction = Vector<T, 3>(query, the_point).unitVector();
          the_gwn += std::copysign(0.5, the_normal.dot(the_direction));
        }
      }

      if(should_recast)
      {
        continue;
      }

      // Rotate the patch so that the discontinuity is aligned with the z-axis
      const double ang = std::acos(axom::utilities::clampVal(cast_direction_local[2], -1.0, 1.0));
      rotator =
        numerics::transforms::axisRotation(ang, cast_direction_local[1], -cast_direction_local[0], 0);
    }

    if(extraTrimming)
    {
      the_gwn +=
        stokes_gwn_evaluate(query, *nurbs_modified, k_quad_npts, field_direction, rotator, quad_tol);
    }
    else
    {
      the_gwn += stokes_gwn_evaluate(query, nurbs, k_quad_npts, field_direction, rotator, quad_tol);
    }

    return the_gwn;
  }

  return the_gwn;
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif
