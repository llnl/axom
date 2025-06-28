// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef PRIMAL_WINDING_NUMBER_3D_DEBUG_HPP_
#define PRIMAL_WINDING_NUMBER_3D_DEBUG_HPP_

// Axom includes
#include "axom/config.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/operators/is_convex.hpp"
#include "axom/primal/operators/squared_distance.hpp"

#include "axom/core/numerics/transforms.hpp"

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
enum class DiscontinuityAxis_debug
{
  x,
  y,
  z,
  rotated
};

// #ifdef AXOM_USE_MFEM
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
void degenerate_surface_processing_debug(const NURBSPatch<T, 3>& patch,
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

/*
 * \brief Computes the GWN for a 3D point wrt a 3D NURBS patch
 *
 * \param [in] query The query point to test
 * \param [in] nPatch The NURBS patch object
 * \param [in] cast_direction The direction of the ray to compute surface intersections
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] ls_tol The tolerance for the line-surface intersection routine
 * \param [in] quad_tol The maximum relative error allowed in the quadrature
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * \param [in] depth The current recursive depth
 * 
 * Computes the generalized winding number for a NURBS patch using Stokes theorem.
 *
 * \pre Assumes that the NURBS patch is trimmed, and has been slightly extended in 
 *       parameter space so that trimming curves arent't on the boundary of the untrimmed patch
 * \return The GWN.
 * 
 * \sa NURBSPatch::makeTriviallyTrimmed
 * \sa NURBSPatch::scaleParameterSpace
 */
template <typename T>
double nurbs_winding_number_debug(const Point<T, 3>& query,
                                  const NURBSPatch<T, 3>& nPatch,
                                  const Vector<T, 3>& cast_direction,
                                  int& node_count,
                                  const double edge_tol = 1e-8,
                                  const double ls_tol = 1e-8,
                                  const double quad_tol = 1e-8,
                                  const double EPS = 1e-8,
                                  const int depth = 0)
{
  // Skip processing of degenerate surfaces
  if(nPatch.getNumControlPoints_u() <= 1 || nPatch.getNumControlPoints_v() <= 1)
  {
    return 0.0;
  }

  // Also skip processing of surfaces with zero trimming curves
  if(nPatch.getNumTrimmingCurves() == 0)
  {
    return 0.0;
  }

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

  // Lambda to rotate the input point using the provided rotation matrix
  auto rotate_point_debug = [&query](const numerics::Matrix<T>& matx,
                                     const Point<T, 3> input) -> Point<T, 3> {
    Vector<T, 3> shifted(query, input);
    Vector<T, 3> rotated;
    numerics::matrix_vector_multiply(matx, shifted.data(), rotated.data());
    return Point<T, 3>({rotated[0] + query[0], rotated[1] + query[1], rotated[2] + query[2]});
  };

  // Lambda to generate an entirely random unit vector
  auto random_unit = []() -> Vector<T, 3> {
    double theta = axom::utilities::random_real(0.0, 2 * M_PI);
    double u = axom::utilities::random_real(-1.0, 1.0);
    return Vector<T, 3> {sin(theta) * sqrt(1 - u * u), cos(theta) * sqrt(1 - u * u), u};
  };

  // Rotation matrix for the patch
  numerics::Matrix<T> rotator;

  // Allocate space for the patch which contains all surface boundaries
  NURBSPatch<T, 3> nPatchWithBoundaries(nPatch), the_disk;

  // Define vector fields whose curl gives us the winding number
  DiscontinuityAxis_debug field_direction;

  // Generate slightly expanded bounding boxes
  auto bBox = nPatch.boundingBox();
  auto oBox = nPatch.orientedBoundingBox();

  auto patch_diameter = bBox.range().norm();

  bBox.expand(0.01 * patch_diameter);
  oBox.expand(0.01 * patch_diameter);

  // Case 1: Exterior without rotations
  if(!bBox.contains(query))
  {
    const bool exterior_x = bBox.getMin()[0] > query[0] || query[0] > bBox.getMax()[0];
    const bool exterior_y = bBox.getMin()[1] > query[1] || query[1] > bBox.getMax()[1];
    const bool exterior_z = bBox.getMin()[2] > query[2] || query[2] > bBox.getMax()[2];

    if(exterior_x || exterior_y)
    {
      field_direction = DiscontinuityAxis_debug::z;
    }
    else if(exterior_y || exterior_z)
    {
      field_direction = DiscontinuityAxis_debug::x;
    }
    else if(exterior_x || exterior_z)
    {
      field_direction = DiscontinuityAxis_debug::y;
    }
  }
  // Case 1.5: Exterior with rotation
  else if(!oBox.contains(query))
  {
    /* The following steps rotate the patch until the OBB is /not/ 
       directly above or below the query point */
    field_direction = DiscontinuityAxis_debug::rotated;

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
    field_direction = DiscontinuityAxis_debug::rotated;
    const Line<T, 3> discontinuity_axis(query, cast_direction);

    // Tolerance for what counts as "close to a boundary" in parameter space
    T disk_radius = 0.01 * nPatch.getParameterSpaceDiagonal();

    // Compute intersections with the *untrimmed and extrapolated* patch
    axom::Array<T> up, vp, tp;
    const bool isHalfOpen = false, countUntrimmed = true;

    bool success = true;
    intersect(discontinuity_axis, nPatch, tp, up, vp, ls_tol, EPS, countUntrimmed, isHalfOpen, success);

    if(!success)
    {
      // Look at the intersection points
      int num_noncoincident = 0;
      for(int i = 0; i < tp.size(); ++i)
      {
        const Point<T, 3> the_point(nPatch.evaluate(up[i], vp[i]));
        // If any of the intersection points are coincident with the surface,
        //  then attempt to clip out all degenerate intersections, and retry
        if(squared_distance(query, the_point) <= edge_tol_sq)
        {
          NURBSPatch<T, 3> clipped_patch1, clipped_patch2;

          degenerate_surface_processing_debug(nPatch,
                                              up,
                                              vp,
                                              0.01 * disk_radius,
                                              clipped_patch1,
                                              clipped_patch2);

          return nurbs_winding_number_debug(query,
                                            clipped_patch1,
                                            cast_direction,
                                            node_count,
                                            edge_tol,
                                            ls_tol,
                                            quad_tol * 1e-5,
                                            EPS,
                                            depth + 1) +
            nurbs_winding_number_debug(query,
                                       clipped_patch2,
                                       cast_direction,
                                       node_count,
                                       edge_tol,
                                       ls_tol,
                                       quad_tol * 1e-5,
                                       EPS,
                                       depth + 1);
        }
        else
        {
          num_noncoincident++;
        }

        // If more than 5 (arbitrary) are *not* coincident with the surface,
        //  re-cast and try again. This is to avoid cases where the point *is*
        //  coincident with the surface, but the first recorded point of
        //  intersection is not after multiple re-casts.
        if(num_noncoincident > 5)
        {
          const auto new_cast_direction = random_unit();
          return nurbs_winding_number_debug(query,
                                            nPatch,
                                            new_cast_direction,
                                            node_count,
                                            edge_tol,
                                            ls_tol,
                                            quad_tol,
                                            EPS,
                                            depth + 1);
        }
      }
    }

    // If no intersections are recorded, then nothing extra to account for

    // Otherwise, account for each discontinuity analytically,
    //  or recursively through disk subdivision
    for(int i = 0; i < up.size(); ++i)
    {
      // Compute the intersection point on the surface
      const Point<T, 3> the_point(nPatch.evaluate(up[i], vp[i]));
      const Vector<T, 3> the_normal = nPatch.normal(up[i], vp[i]);

      // Check for bad intersections, i.e.,
      //  > There normal is poorly defined (cusp)
      //  > The normal is tangent to the axis of discontinuity
      const bool bad_intersection = axom::utilities::isNearlyEqual(the_normal.norm(), 0.0, EPS) ||
        axom::utilities::isNearlyEqual(the_normal.unitVector().dot(cast_direction), 0.0, EPS);

      const bool isOnSurface = squared_distance(query, the_point) <= edge_tol_sq;

      if(bad_intersection && !isOnSurface)
      {
        // If a non-coincident ray intersects the surface at a tangent/cusp,
        //  can recast and try again
        const auto new_cast_direction = random_unit();
        return nurbs_winding_number_debug(query,
                                          nPatch,
                                          new_cast_direction,
                                          node_count,
                                          edge_tol,
                                          ls_tol,
                                          quad_tol,
                                          EPS,
                                          depth + 1);
      }

      if(isOnSurface)
      {
        // If the query point is on the surface, then shrink the disk
        //  to ensure its winding number is known to be near-zero
        disk_radius = 0.01 * disk_radius;
      }

      // Consider a disk around the intersection point via NURBSPatch::diskSplit.
      //   If the disk intersects any trimming curves, need to do disk subdivision.
      //   If not, we can compute the winding number without changing the trimming curvse
      const bool ignoreInteriorDisk = true, clipDisk = true;
      bool isDiskInside, isDiskOutside;
      NURBSPatch<T, 3> the_disk;
      nPatchWithBoundaries.diskSplit(up[i],
                                     vp[i],
                                     disk_radius,
                                     the_disk,
                                     nPatchWithBoundaries,
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
        const auto new_cast_direction = random_unit();
        the_gwn += nurbs_winding_number_debug(query,
                                              the_disk,
                                              new_cast_direction,
                                              node_count,
                                              edge_tol,
                                              ls_tol,
                                              quad_tol,
                                              EPS,
                                              depth + 1);
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

    // Rotate the patch so that the discontinuity is aligned with the z-axis
    const double ang = std::acos(axom::utilities::clampVal(cast_direction[2], -1.0, 1.0));
    rotator = numerics::transforms::axisRotation(ang, cast_direction[1], -cast_direction[0], 0);
  }

  if(field_direction == DiscontinuityAxis_debug::rotated)
  {
    // The trimming curves for rotatedPatch have been changed as needed,
    //  but we need to rotate the control points
    const auto patch_shape = nPatchWithBoundaries.getControlPoints().shape();
    for(int i = 0; i < patch_shape[0]; ++i)
    {
      for(int j = 0; j < patch_shape[1]; ++j)
      {
        nPatchWithBoundaries(i, j) = rotate_point_debug(rotator, nPatchWithBoundaries(i, j));
      }
    }
  }

  the_gwn +=
    stokes_gwn_evaluate(query, nPatchWithBoundaries, field_direction, quad_npts, node_count, quad_tol);

  return the_gwn;
}

/*
 * \brief Computes the GWN for a 3D point wrt a 3D NURBS patch with precomputed data
 *
 * \param [in] query The query point to test
 * \param [in] nPatch The NURBS patch object with precomputed data
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] ls_tol The tolerance for the line-surface intersection routine
 * \param [in] quad_tol The maximum relative error allowed in the quadrature
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * \param [in] depth The current recursive depth
 * 
 * Computes the generalized winding number for a NURBS patch using Stokes theorem.
 *
 * \pre Assumes that the NURBS patch is trimmed, and has been slightly extended in 
 *       parameter space so that trimming curves arent't on the boundary of the untrimmed patch
 * \return The GWN.
 */
template <typename T>
double nurbs_data_winding_number_debug(const Point<T, 3>& query,
                                       const NURBSPatchData<T>& nPatchData,
                                       const Vector<T, 3>& cast_direction,
                                       int& node_count,
                                       int& case_code,
                                       int& integrated_curves,
                                       const double edge_tol = 1e-8,
                                       const double ls_tol = 1e-8,
                                       const double quad_tol = 1e-8,
                                       const double EPS = 1e-8,
                                       const int depth = 0)
{
  // Skip processing of degenerate surfaces
  if(nPatchData.patch.getNumControlPoints_u() <= 1 || nPatchData.patch.getNumControlPoints_v() <= 1)
  {
    return 0.0;
  }

  // Also skip processing of surfaces with zero trimming curves
  if(nPatchData.patch.getNumTrimmingCurves() == 0)
  {
    return 0.0;
  }

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

  // Lambda to rotate the input point using the provided rotation matrix
  auto rotate_point_debug = [&query](const numerics::Matrix<T>& matx,
                                     const Point<T, 3> input) -> Point<T, 3> {
    Vector<T, 3> shifted(query, input);
    Vector<T, 3> rotated;
    numerics::matrix_vector_multiply(matx, shifted.data(), rotated.data());
    return Point<T, 3>({rotated[0] + query[0], rotated[1] + query[1], rotated[2] + query[2]});
  };

  // Lambda to generate an entirely random unit vector
  auto random_unit = []() -> Vector<T, 3> {
    double theta = axom::utilities::random_real(0.0, 2 * M_PI);
    double u = axom::utilities::random_real(-1.0, 1.0);
    return Vector<T, 3> {sin(theta) * sqrt(1 - u * u), cos(theta) * sqrt(1 - u * u), u};
  };

  // Rotation matrix for the patch
  numerics::Matrix<T> rotator;

  // Allocate space for the patch which contains all surface boundaries
  NURBSPatch<T, 3> nPatchWithBoundaries(nPatchData.patch), the_disk;

  // Define vector fields whose curl gives us the winding number
  DiscontinuityAxis_debug field_direction;
  bool extraTrimming = false;

  // Generate slightly expanded bounding boxes
  auto bBox = nPatchData.bbox;
  auto oBox = nPatchData.obox;

  auto patch_diameter = bBox.range().norm();

  bBox.expand(0.01 * patch_diameter);
  oBox.expand(0.01 * patch_diameter);

  // Case 1: Exterior without rotations
  if(!bBox.contains(query))
  {
    case_code = 0;
    integrated_curves = nPatchWithBoundaries.getNumTrimmingCurves();

    const bool exterior_x = bBox.getMin()[0] > query[0] || query[0] > bBox.getMax()[0];
    const bool exterior_y = bBox.getMin()[1] > query[1] || query[1] > bBox.getMax()[1];
    const bool exterior_z = bBox.getMin()[2] > query[2] || query[2] > bBox.getMax()[2];

    if(exterior_x || exterior_y)
    {
      field_direction = DiscontinuityAxis_debug::z;
    }
    else if(exterior_y || exterior_z)
    {
      field_direction = DiscontinuityAxis_debug::x;
    }
    else if(exterior_x || exterior_z)
    {
      field_direction = DiscontinuityAxis_debug::y;
    }
  }
  // Case 1.5: Exterior with rotation
  else if(!oBox.contains(query))
  {
    case_code = 1;
    integrated_curves = nPatchWithBoundaries.getNumTrimmingCurves();

    /* The following steps rotate the patch until the OBB is /not/ 
       directly above or below the query point */
    field_direction = DiscontinuityAxis_debug::rotated;

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
    case_code = 2;
    integrated_curves = nPatchWithBoundaries.getNumTrimmingCurves();

    field_direction = DiscontinuityAxis_debug::rotated;
    const Line<T, 3> discontinuity_axis(query, cast_direction);

    // Tolerance for what counts as "close to a boundary" in parameter space
    T disk_radius = 0.01 * nPatchData.pbox_diag;

    // Compute intersections with the *untrimmed and extrapolated* patch
    axom::Array<T> up, vp, tp;
    const bool isHalfOpen = false, countUntrimmed = true;

    bool success = true;
    intersect(discontinuity_axis,
              nPatchData.patch,
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
        const Point<T, 3> the_point(nPatchData.patch.evaluate(up[i], vp[i]));
        // If any of the intersection points are coincident with the surface,
        //  then attempt to clip out all degenerate intersections, and retry
        if(squared_distance(query, the_point) <= edge_tol_sq)
        {
          NURBSPatch<T, 3> clipped_patch1, clipped_patch2;

          degenerate_surface_processing_debug(nPatchData.patch,
                                              up,
                                              vp,
                                              0.01 * disk_radius,
                                              clipped_patch1,
                                              clipped_patch2);

          return nurbs_winding_number_debug(query,
                                            clipped_patch1,
                                            cast_direction,
                                            node_count,
                                            edge_tol,
                                            ls_tol,
                                            quad_tol * 1e-5,
                                            EPS,
                                            depth + 1) +
            nurbs_winding_number_debug(query,
                                       clipped_patch2,
                                       cast_direction,
                                       node_count,
                                       edge_tol,
                                       ls_tol,
                                       quad_tol * 1e-5,
                                       EPS,
                                       depth + 1);
        }
        else
        {
          num_noncoincident++;
        }

        // If more than 5 (arbitrary) are *not* coincident with the surface,
        //  re-cast and try again. This is to avoid cases where the point *is*
        //  coincident with the surface, but the first recorded point of
        //  intersection is not after multiple re-casts.
        if(num_noncoincident > 5)
        {
          const auto new_cast_direction = random_unit();
          return nurbs_winding_number_debug(query,
                                            nPatchData.patch,
                                            new_cast_direction,
                                            node_count,
                                            edge_tol,
                                            ls_tol,
                                            quad_tol,
                                            EPS,
                                            depth + 1);
        }
      }
    }

    // If no intersections are recorded, then nothing extra to account for

    // Otherwise, account for each discontinuity analytically,
    //  or recursively through disk subdivision
    for(int i = 0; i < up.size(); ++i)
    {
      // Compute the intersection point on the surface
      const Point<T, 3> the_point(nPatchData.patch.evaluate(up[i], vp[i]));
      const Vector<T, 3> the_normal = nPatchData.patch.normal(up[i], vp[i]);

      // Check for bad intersections, i.e.,
      //  > There normal is poorly defined (cusp)
      //  > The normal is tangent to the axis of discontinuity
      const bool bad_intersection = axom::utilities::isNearlyEqual(the_normal.norm(), 0.0, EPS) ||
        axom::utilities::isNearlyEqual(the_normal.unitVector().dot(cast_direction), 0.0, EPS);

      const bool isOnSurface = squared_distance(query, the_point) <= edge_tol_sq;

      if(bad_intersection && !isOnSurface)
      {
        // If a non-coincident ray intersects the surface at a tangent/cusp,
        //  can recast and try again
        const auto new_cast_direction = random_unit();
        return nurbs_winding_number_debug(query,
                                          nPatchData.patch,
                                          new_cast_direction,
                                          node_count,
                                          edge_tol,
                                          ls_tol,
                                          quad_tol,
                                          EPS,
                                          depth + 1);
      }

      if(isOnSurface)
      {
        // If the query point is on the surface, then shrink the disk
        //  to ensure its winding number is known to be near-zero
        disk_radius = 0.01 * disk_radius;
      }

      // Consider a disk around the intersection point via NURBSPatch::diskSplit.
      //   If the disk intersects any trimming curves, need to do disk subdivision.
      //   If not, we can compute the winding number without changing the trimming curvse
      const bool ignoreInteriorDisk = true, clipDisk = true;
      bool isDiskInside, isDiskOutside;
      int old_num_trim = nPatchWithBoundaries.getNumTrimmingCurves();

      NURBSPatch<T, 3> the_disk;
      nPatchWithBoundaries.diskSplit(up[i],
                                     vp[i],
                                     disk_radius,
                                     the_disk,
                                     nPatchWithBoundaries,
                                     isDiskInside,
                                     isDiskOutside,
                                     ignoreInteriorDisk,
                                     clipDisk);

      extraTrimming =
        extraTrimming || (!isDiskInside && !isDiskOutside) || (isDiskInside && !ignoreInteriorDisk);

      if(extraTrimming)
      {
        case_code = 3;
        integrated_curves += the_disk.getNumTrimmingCurves() +
          (nPatchWithBoundaries.getNumTrimmingCurves() - old_num_trim);
      }

      if(isOnSurface)
      {
        // If the query point is on the surface, the contribution of the disk is near-zero
        //  and we only needed to puncture the larger surface to proceed
        continue;
      }
      else if(!isDiskInside && !isDiskOutside)
      {
        // If the disk overlapped with the trimming curves, evaluate the winding number for the disk
        const auto new_cast_direction = random_unit();
        the_gwn += nurbs_winding_number_debug(query,
                                              the_disk,
                                              new_cast_direction,
                                              node_count,
                                              edge_tol,
                                              ls_tol,
                                              quad_tol,
                                              EPS,
                                              depth + 1);
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

    // Rotate the patch so that the discontinuity is aligned with the z-axis
    const double ang = std::acos(axom::utilities::clampVal(cast_direction[2], -1.0, 1.0));
    rotator = numerics::transforms::axisRotation(ang, cast_direction[1], -cast_direction[0], 0);
  }

  if(extraTrimming)
  {
    // Can't use cached quadrature rules
    if(field_direction == DiscontinuityAxis_debug::rotated)
    {
      // The trimming curves for rotatedPatch have been changed as needed,
      //  but we need to rotate the control points
      auto patch_shape = nPatchWithBoundaries.getControlPoints().shape();
      for(int i = 0; i < patch_shape[0]; ++i)
      {
        for(int j = 0; j < patch_shape[1]; ++j)
        {
          nPatchWithBoundaries(i, j) = rotate_point_debug(rotator, nPatchWithBoundaries(i, j));
        }
      }
    }

    the_gwn +=
      stokes_gwn_evaluate(query, nPatchWithBoundaries, field_direction, quad_npts, node_count, quad_tol);
  }
  else
  {
    if(field_direction != DiscontinuityAxis_debug::rotated)
    {
      the_gwn +=
        stokes_gwn_evaluate_cached(query, nPatchData, field_direction, quad_npts, node_count, quad_tol);
    }
    else
    {
      the_gwn +=
        stokes_gwn_evaluate_cached_rotated(query, nPatchData, rotator, quad_npts, node_count, quad_tol);
    }
  }

  return the_gwn;
}

// ================================ ORIGINAL VERSIONS ================================

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
double stokes_gwn_evaluate(const Point<T, 3>& query,
                           const NURBSPatch<T, 3>& patch,
                           const DiscontinuityAxis_debug ax,
                           int npts,
                           int& node_count,
                           double quad_tol)
{
  // Generate the quadrature rules in parameter space
  static mfem::IntegrationRules my_IntRules(0, mfem::Quadrature1D::GaussLegendre);
  const mfem::IntegrationRule& quad_rule = my_IntRules.Get(mfem::Geometry::SEGMENT, 2 * npts - 1);

  double quad = 0;
  for(int n = 0; n < patch.getNumTrimmingCurves(); ++n)
  {
    NURBSCurve<T, 2> trimming_curve(patch.getTrimmingCurve(n));
    double quad_coarse =
      stokes_gwn_component(query, trimming_curve, patch, ax, 0, 0, quad_rule, node_count);

    if(quad_tol >= 0)
    {
      quad += stokes_gwn_adaptive(query,
                                  trimming_curve,
                                  quad_rule,
                                  patch,
                                  ax,
                                  0,
                                  0,
                                  quad_coarse,
                                  node_count,
                                  quad_tol);
    }
    else
    {
      quad += quad_coarse;
    }
  }

  return 0.25 * M_1_PI * quad;
}

template <typename T>
double stokes_gwn_adaptive(const Point<T, 3>& query,
                           const NURBSCurve<T, 2>& curve,
                           const mfem::IntegrationRule& quad_rule,
                           const NURBSPatch<T, 3>& patch,
                           const DiscontinuityAxis_debug ax,
                           const int refinement_level,
                           const int refinement_index,
                           const double quad_coarse,
                           int& node_count,
                           const double quad_tol)
{
  double quad_fine_1 = stokes_gwn_component(query,
                                            curve,
                                            patch,
                                            ax,
                                            refinement_level + 1,
                                            2 * refinement_index,
                                            quad_rule,
                                            node_count);
  double quad_fine_2 = stokes_gwn_component(query,
                                            curve,
                                            patch,
                                            ax,
                                            refinement_level + 1,
                                            2 * refinement_index + 1,
                                            quad_rule,
                                            node_count);

  if(refinement_level > 25 ||
     axom::utilities::isNearlyEqualRelative(quad_coarse, quad_fine_1 + quad_fine_2, quad_tol, 1e-10))
  {
    return quad_fine_1 + quad_fine_2;
  }

  // Perform the adaptive call over both halves of the curve
  quad_fine_1 = stokes_gwn_adaptive(query,
                                    curve,
                                    quad_rule,
                                    patch,
                                    ax,
                                    refinement_level + 1,
                                    2 * refinement_index,
                                    quad_fine_1,
                                    node_count,
                                    quad_tol);

  quad_fine_2 = stokes_gwn_adaptive(query,
                                    curve,
                                    quad_rule,
                                    patch,
                                    ax,
                                    refinement_level + 1,
                                    2 * refinement_index + 1,
                                    quad_fine_2,
                                    node_count,
                                    quad_tol);

  return quad_fine_1 + quad_fine_2;
}

template <typename T>
double stokes_gwn_component(const Point<T, 3>& query,
                            const NURBSCurve<T, 2>& curve,
                            const NURBSPatch<T, 3>& patch,
                            const DiscontinuityAxis_debug ax,
                            const int refinement_level,
                            const int refinement_index,
                            const mfem::IntegrationRule& quad_rule,
                            int& node_count)
{
  node_count += quad_rule.GetNPoints();

  double this_quad = 0;
  for(int q = 0; q < quad_rule.GetNPoints(); ++q)
  {
    // Find the right knot span based on refinement level
    const double span_length = (curve.getMaxKnot() - curve.getMinKnot()) / (1 << refinement_level);
    const double span_offset = curve.getMinKnot() + span_length * refinement_index;

    const double quad_x = span_offset + span_length * quad_rule.IntPoint(q).x;
    const double quad_weight = quad_rule.IntPoint(q).weight * span_length;

    Point<T, 2> c_eval;
    Vector<T, 2> c_Dt;
    curve.evaluateFirstDerivative(quad_x, c_eval, c_Dt);

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
    case(DiscontinuityAxis_debug::x):
      this_quad += quad_weight * (node[2] * node[0] * node_dt[1] - node[1] * node[0] * node_dt[2]) /
        (node[1] * node[1] + node[2] * node[2]) / node.norm();
      break;
    case(DiscontinuityAxis_debug::y):
      this_quad += quad_weight * (node[0] * node[1] * node_dt[2] - node[2] * node[1] * node_dt[0]) /
        (node[0] * node[0] + node[2] * node[2]) / node.norm();
      break;
    case(DiscontinuityAxis_debug::z):
    case(DiscontinuityAxis_debug::rotated):
      this_quad += quad_weight * (node[1] * node[2] * node_dt[0] - node[0] * node[2] * node_dt[1]) /
        (node[0] * node[0] + node[1] * node[1]) / node.norm();
      break;
    }
  }

  return this_quad;
}

// ================================ CACHED VERSIONS ================================

template <typename T>
double stokes_gwn_evaluate_cached(const Point<T, 3>& query,
                                  const NURBSPatchData<T>& nPatchData,
                                  const DiscontinuityAxis_debug ax,
                                  const int quad_npts,
                                  int& node_count,
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
      stokes_gwn_component_cached(query, quad_rule, ax, trimming_curve_data, node_count);

    if(quad_tol >= 0)
    {
      quad += 0.25 * M_1_PI *
        stokes_gwn_adaptive_cached(query,
                                   nPatchData,
                                   quad_rule,
                                   ax,
                                   n,
                                   0,
                                   0,
                                   quad_coarse,
                                   node_count,
                                   quad_tol);
    }
    else
    {
      quad += 0.25 * M_1_PI * quad_coarse;
    }
  }

  return quad;
}

template <typename T>
double stokes_gwn_adaptive_cached(const Point<T, 3>& query,
                                  const NURBSPatchData<T>& nPatchData,
                                  const mfem::IntegrationRule& quad_rule,
                                  const DiscontinuityAxis_debug ax,
                                  const int curve_index,
                                  const int refinement_level,
                                  const int refinement_index,
                                  const double quad_coarse,
                                  int& node_count,
                                  const double quad_tol)
{
  auto trimming_curve_data_1 =
    nPatchData.getQuadratureData(curve_index, quad_rule, refinement_level + 1, 2 * refinement_index);
  auto trimming_curve_data_2 = nPatchData.getQuadratureData(curve_index,
                                                            quad_rule,
                                                            refinement_level + 1,
                                                            2 * refinement_index + 1);

  double quad_fine_1 =
    stokes_gwn_component_cached(query, quad_rule, ax, trimming_curve_data_1, node_count);
  double quad_fine_2 =
    stokes_gwn_component_cached(query, quad_rule, ax, trimming_curve_data_2, node_count);

  if(refinement_level >= 25 ||
     axom::utilities::isNearlyEqualRelative(quad_fine_1 + quad_fine_2, quad_coarse, quad_tol, 1e-10))
  {
    return quad_fine_1 + quad_fine_2;
  }

  quad_fine_1 = stokes_gwn_adaptive_cached(query,
                                           nPatchData,
                                           quad_rule,
                                           ax,
                                           curve_index,
                                           refinement_level + 1,
                                           2 * refinement_index,
                                           quad_fine_1,
                                           node_count,
                                           quad_tol);

  quad_fine_2 = stokes_gwn_adaptive_cached(query,
                                           nPatchData,
                                           quad_rule,
                                           ax,
                                           curve_index,
                                           refinement_level + 1,
                                           2 * refinement_index + 1,
                                           quad_fine_2,
                                           node_count,
                                           quad_tol);

  return quad_fine_1 + quad_fine_2;
}

template <typename T>
double stokes_gwn_component_cached(const Point<T, 3>& query,
                                   const mfem::IntegrationRule& quad_rule,
                                   const DiscontinuityAxis_debug ax,
                                   const TrimmingCurveQuadratureData<T>& trimming_curve_data,
                                   int& node_count)
{
  node_count += quad_rule.GetNPoints();

  // Do this without refinement
  double this_quad = 0;
  for(int q = 0; q < quad_rule.GetNPoints(); ++q)
  {
    const Vector<T, 3> node(query, trimming_curve_data.quadrature_points[q].first);
    const Vector<T, 3> node_dt(trimming_curve_data.quadrature_points[q].second);
    const double node_norm = node.norm();

    const double quad_weight = quad_rule.IntPoint(q).weight * trimming_curve_data.span_length;

    // Compute one of three vector field line integrals depending on
    //  the orientation of the original surface, indicated through ax.
    switch(ax)
    {
    case(DiscontinuityAxis_debug::x):
      this_quad += quad_weight * (node[2] * node[0] * node_dt[1] - node[1] * node[0] * node_dt[2]) /
        (node[1] * node[1] + node[2] * node[2]) / node_norm;
      break;
    case(DiscontinuityAxis_debug::y):
      this_quad += quad_weight * (node[0] * node[1] * node_dt[2] - node[2] * node[1] * node_dt[0]) /
        (node[0] * node[0] + node[2] * node[2]) / node_norm;
      break;
    case(DiscontinuityAxis_debug::z):
      this_quad += quad_weight * (node[1] * node[2] * node_dt[0] - node[0] * node[2] * node_dt[1]) /
        (node[0] * node[0] + node[1] * node[1]) / node_norm;
      break;
    }
  }

  return this_quad;
}

// ================================ CACHED ROTATED VERSIONS ================================

template <typename T>
Point<T, 3> rotate_point_debug(const numerics::Matrix<T>& matx,
                               const Point<T, 3>& center,
                               const Point<T, 3>& input)
{
  Vector<T, 3> shifted(center, input);
  Vector<T, 3> rotated;
  numerics::matrix_vector_multiply(matx, shifted.data(), rotated.data());
  return Point<T, 3>({rotated[0] + center[0], rotated[1] + center[1], rotated[2] + center[2]});
}

template <typename T>
Point<T, 3> rotate_vector_debug(const numerics::Matrix<T>& matx,
                                const Point<T, 3>& center,
                                const Vector<T, 3>& input)
{
  Vector<T, 3> shifted {input[0] - center[0], input[1] - center[1], input[2] - center[2]};
  Vector<T, 3> rotated;
  numerics::matrix_vector_multiply(matx, shifted.data(), rotated.data());
  return Point<T, 3>({rotated[0] + center[0], rotated[1] + center[1], rotated[2] + center[2]});
}

template <typename T>
double stokes_gwn_evaluate_cached_rotated(const Point<T, 3>& query,
                                          const NURBSPatchData<T>& nPatchData,
                                          const axom::numerics::Matrix<T>& rotator,
                                          const int quad_npts,
                                          int& node_count,
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
      stokes_gwn_component_cached_rotated(query, quad_rule, rotator, trimming_curve_data, node_count);

    if(quad_tol >= 0)
    {
      quad += 0.25 * M_1_PI *
        stokes_gwn_adaptive_cached_rotated(query,
                                           nPatchData,
                                           quad_rule,
                                           rotator,
                                           n,
                                           0,
                                           0,
                                           quad_coarse,
                                           node_count,
                                           quad_tol);
    }
    else
    {
      quad += 0.25 * M_1_PI * quad_coarse;
    }
  }

  return quad;
}

template <typename T>
double stokes_gwn_adaptive_cached_rotated(const Point<T, 3>& query,
                                          const NURBSPatchData<T>& nPatchData,
                                          const mfem::IntegrationRule& quad_rule,
                                          const axom::numerics::Matrix<T>& rotator,
                                          const int curve_index,
                                          const int refinement_level,
                                          const int refinement_index,
                                          const double quad_coarse,
                                          int& node_count,
                                          const double quad_tol)
{
  auto trimming_curve_data_1 =
    nPatchData.getQuadratureData(curve_index, quad_rule, refinement_level + 1, 2 * refinement_index);
  auto trimming_curve_data_2 = nPatchData.getQuadratureData(curve_index,
                                                            quad_rule,
                                                            refinement_level + 1,
                                                            2 * refinement_index + 1);

  double quad_fine_1 =
    stokes_gwn_component_cached_rotated(query, quad_rule, rotator, trimming_curve_data_1, node_count);
  double quad_fine_2 =
    stokes_gwn_component_cached_rotated(query, quad_rule, rotator, trimming_curve_data_2, node_count);

  if(refinement_level >= 25 ||
     axom::utilities::isNearlyEqualRelative(quad_fine_1 + quad_fine_2, quad_coarse, quad_tol, 1e-10))
  {
    return quad_fine_1 + quad_fine_2;
  }

  quad_fine_1 = stokes_gwn_adaptive_cached_rotated(query,
                                                   nPatchData,
                                                   quad_rule,
                                                   rotator,
                                                   curve_index,
                                                   refinement_level + 1,
                                                   2 * refinement_index,
                                                   quad_fine_1,
                                                   node_count,
                                                   quad_tol);

  quad_fine_2 = stokes_gwn_adaptive_cached_rotated(query,
                                                   nPatchData,
                                                   quad_rule,
                                                   rotator,
                                                   curve_index,
                                                   refinement_level + 1,
                                                   2 * refinement_index + 1,
                                                   quad_fine_2,
                                                   node_count,
                                                   quad_tol);

  return quad_fine_1 + quad_fine_2;
}

template <typename T>
double stokes_gwn_component_cached_rotated(const Point<T, 3>& query,
                                           const mfem::IntegrationRule& quad_rule,
                                           const axom::numerics::Matrix<T>& rotator,
                                           const TrimmingCurveQuadratureData<T>& trimming_curve_data,
                                           int& node_count)
{
  node_count += quad_rule.GetNPoints();

  // Do this without refinement
  double this_quad = 0;
  for(int q = 0; q < quad_rule.GetNPoints(); ++q)
  {
    const Vector<T, 3> node(
      query,
      rotate_point_debug(rotator, query, trimming_curve_data.quadrature_points[q].first));
    const Vector<T, 3> node_dt(rotate_vector_debug(rotator,
                                                   Point<T, 3>({0.0, 0.0, 0.0}),
                                                   trimming_curve_data.quadrature_points[q].second));
    const double node_norm = node.norm();

    const double quad_weight = quad_rule.IntPoint(q).weight * trimming_curve_data.span_length;

    // This is the z-aligned axis
    this_quad += quad_weight * (node[1] * node[2] * node_dt[0] - node[0] * node[2] * node_dt[1]) /
      (node[0] * node[0] + node[1] * node[1]) / node_norm;
  }

  return this_quad;
}

template <typename T>
double nurbs_winding_number_surface_debug(const Point<T, 3>& query,
                                          const BezierPatch<T, 3>& patch,
                                          int& node_count,
                                          const double quad_tol = 1e-8,
                                          const int depth = 0)
{
  // Generate the quadrature rules in parameter space
  int quad_npts = 15;

  static mfem::IntegrationRules my_IntRules(0, mfem::Quadrature1D::GaussLegendre);
  const mfem::IntegrationRule& quad_rule =
    my_IntRules.Get(mfem::Geometry::SEGMENT, 2 * quad_npts - 1);

  // Get the quadrature on the original patch
  double quad_coarse = surface_winding_number_component(query, patch, quad_rule, node_count);
  if(quad_tol >= 0.0)
  {
    return surface_winding_number_adaptive(query, patch, quad_coarse, quad_rule, quad_tol, node_count);
  }
  else
  {
    // If the quad_tol is negative, we are not using adaptive quadrature
    //  and we can return the result immediately
    return quad_coarse;
  }
}

template <typename T>
double surface_winding_number_adaptive(const Point<T, 3>& query,
                                       const BezierPatch<T, 3>& patch,
                                       double quad_coarse,
                                       const mfem::IntegrationRule& quad_rule,
                                       double quad_tol,
                                       int& node_count)
{
  BezierPatch<T, 3> subpatches[4];
  double sub_gwn[4];

  patch.split(0.5, 0.5, subpatches[0], subpatches[1], subpatches[2], subpatches[3]);
  for(int i = 0; i < 4; ++i)
  {
    // Recursively compute the GWN for each subpatch
    sub_gwn[i] = surface_winding_number_component(query, subpatches[i], quad_rule, node_count);
  }

  double quad_fine = sub_gwn[0] + sub_gwn[1] + sub_gwn[2] + sub_gwn[3];
  if(axom::utilities::isNearlyEqualRelative(quad_coarse, quad_fine, quad_tol, 1e-10))
  {
    return quad_fine;
  }
   
  for(int i = 0; i < 4; ++i)
  {
    sub_gwn[i] =
      surface_winding_number_adaptive(query, subpatches[i], sub_gwn[i], quad_rule, quad_tol, node_count);
  }

  return sub_gwn[0] + sub_gwn[1] + sub_gwn[2] + sub_gwn[3];
}

template <typename T>
double surface_winding_number_component(const Point<T, 3>& query,
                                        const BezierPatch<T, 3>& patch,
                                        const mfem::IntegrationRule& quad_rule,
                                        int& node_count)
{
  node_count += quad_rule.GetNPoints() * quad_rule.GetNPoints();

  double quadrature = 0.0;
  for(int qu = 0; qu < quad_rule.GetNPoints(); ++qu)
  {
    for(int qv = 0; qv < quad_rule.GetNPoints(); ++qv)
    {
      Vector<T, 3> node(query, patch.evaluate(quad_rule.IntPoint(qu).x, quad_rule.IntPoint(qv).x));

      // Compute the normal vector
      Vector<T, 3> normal = patch.normal(quad_rule.IntPoint(qu).x, quad_rule.IntPoint(qv).x);

      quadrature += quad_rule.IntPoint(qu).weight * quad_rule.IntPoint(qv).weight *
        Vector<T, 3>::dot_product(node, normal) / std::pow(node.norm(), 3);
    }
  }

  return 0.25 * M_1_PI * quadrature;
}
// #endif

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif