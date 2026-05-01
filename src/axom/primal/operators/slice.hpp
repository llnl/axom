// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_SLICE_HPP_
#define AXOM_PRIMAL_SLICE_HPP_
#include "axom/core/utilities/Utilities.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/Plane.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"

#include "axom/primal/operators/intersect.hpp"
#include "axom/primal/operators/detail/slice_impl.hpp"

namespace axom
{
namespace primal
{

/*!
 * \brief Slices a 3D tetrahedron with a plane and returns the resulting
 *        polygon.
 *
 * \param [in] tet The tetrahedron to slice
 * \param [in] plane The slicing plane
 *
 * \return The polygon obtained from slicing a tetrahedron with a plane.
 *   When the plane intersects the tetrahedron in a nondegenerate cross
 *   section, the return value is a triangle or quadrilateral. Degenerate
 *   intersections are also represented as polygons: a plane that touches the
 *   tetrahedron at a single vertex returns a one-vertex polygon, and a plane
 *   that intersects the tetrahedron only along one edge returns a two-vertex
 *   polygon.
 *
 * \note For nondegenerate intersections, the polygon vertices are intended to
 *   be ordered so that the polygon normal is aligned with the plane normal.
 */
template <typename T, PolygonArray ARRAY_TYPE = PolygonArray::Dynamic, int MAX_VERTS = DEFAULT_MAX_NUM_VERTICES>
AXOM_HOST_DEVICE primal::Polygon<T, 3, ARRAY_TYPE, MAX_VERTS> slice(const primal::Tetrahedron<T, 3>& tet,
                                                                    const primal::Plane<T, 3>& plane)
{
  return detail::slice_tet_plane<T, ARRAY_TYPE, MAX_VERTS>(tet, plane);
}

}  // namespace primal
}  // namespace axom
#endif
