// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file compute_bounding_box.hpp
 *
 * \brief Consists of functions to create bounding boxes.
 */

#ifndef AXOM_PRIMAL_COMPUTE_BOUNDING_BOX_HPP_
#define AXOM_PRIMAL_COMPUTE_BOUNDING_BOX_HPP_

#include "axom/core/numerics/Matrix.hpp"       // for Matrix
#include "axom/core/Macros.hpp"                // for AXOM_HOST__DEVICE
#include "axom/core/numerics/eigen_solve.hpp"  // for eigen_solve
#include "axom/primal/geometry/Polyhedron.hpp"
#include "axom/primal/geometry/Octahedron.hpp"
#include "axom/primal/geometry/Hexahedron.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Quadrilateral.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"

namespace axom
{
namespace primal
{
/*!
 * \brief Creates a bounding box which contains the collection of passed in
 * points.
 *
 * Call srand() to initialize the random number generator before using this
 * function.
 *
 * \param [in] pts array of points
 * \param [in] n number of points
 * \note if n <= 0, invokes default constructor
 */
template <typename T, int NDIMS>
OrientedBoundingBox<T, NDIMS> compute_oriented_bounding_box(const Point<T, NDIMS> *pts, int n)
{
  return OrientedBoundingBox<T, NDIMS>(pts, n);
}

/*!
 * \brief Creates an oriented bounding box which contains the passed in OBBs.
 *
 * Call srand() to initialize the random number generator before using this
 * function.
 *
 * \param [in] l left obb
 * \param [in] r right obb
 */
template <typename T, int NDIMS>
OrientedBoundingBox<T, NDIMS> merge_boxes(const OrientedBoundingBox<T, NDIMS> &l,
                                          const OrientedBoundingBox<T, NDIMS> &r)
{
  // TODO: See if this initial check can be improved so it's not so costly in
  // cases where it doesn't end up helping
  if(l.contains(r))
  {
    return l;
  }
  if(r.contains(l))
  {
    return r;
  }

  std::vector<Point<T, NDIMS>> lv = l.vertices();
  std::vector<Point<T, NDIMS>> rv = r.vertices();

  const int size = (1 << NDIMS);

  Point<T, NDIMS> pts[2 * size];

  for(int i = 0; i < size; i++)
  {
    pts[i] = lv[i];
    pts[i + size] = rv[i];
  }

  return compute_oriented_bounding_box<T, NDIMS>(pts, 2 * size);
}

/*!
 * \brief Constructor. Creates a bounding box which contains the passed in
 * bounding boxes.
 *
 * \param [in] l left bb
 * \param [in] r right bb
 */

template <typename T, int NDIMS>
BoundingBox<T, NDIMS> merge_boxes(const BoundingBox<T, NDIMS> &l, const BoundingBox<T, NDIMS> &r)
{
  BoundingBox<T, NDIMS> res(l);
  res.addBox(r);
  return res;
}

/*!
 * \brief Creates a bounding box around a Triangle
 * \accelerated
 * \param [in] tri The Triangle
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE BoundingBox<T, NDIMS> compute_bounding_box(const Triangle<T, NDIMS> &tri)
{
  return BoundingBox<T, NDIMS> {tri[0], tri[1], tri[2]};
}

/*!
 * \brief Creates a bounding box around a Quadrilateral
 * \accelerated
 * \param [in] quad The Quadrilateral
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE BoundingBox<T, NDIMS> compute_bounding_box(const Quadrilateral<T, NDIMS> &quad)
{
  return BoundingBox<T, NDIMS> {quad[0], quad[1], quad[2], quad[3]};
}

/*!
 * \brief Creates a bounding box around an Octahedron
 *
 * \param [in] oct The Octahedron
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE BoundingBox<T, NDIMS> compute_bounding_box(const Octahedron<T, NDIMS> &oct)
{
  return BoundingBox<T, NDIMS> {oct[0], oct[1], oct[2], oct[3], oct[4], oct[5]};
}

/*!
 * \brief Creates a bounding box around a Hexahedron
 *
 * \param [in] hex The Hexahedron
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE BoundingBox<T, NDIMS> compute_bounding_box(const Hexahedron<T, NDIMS> &hex)
{
  return BoundingBox<T, NDIMS> {hex[0], hex[1], hex[2], hex[3], hex[4], hex[5], hex[6], hex[7]};
}

/*!
 * \brief Creates a bounding box around a Polyhedron
 *
 * \param [in] poly The Polyhedron
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE BoundingBox<T, NDIMS> compute_bounding_box(const Polyhedron<T, NDIMS> &poly)
{
  BoundingBox<T, NDIMS> res;
  for(int i = 0; i < poly.numVertices(); i++)
  {
    res.addPoint(poly[i]);
  }
  return res;
}

/*!
 * \brief Creates a bounding box around a Tetrahedron
 *
 * \param [in] tet The tetrahedron
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE BoundingBox<T, NDIMS> compute_bounding_box(const Tetrahedron<T, NDIMS> &tet)
{
  return BoundingBox<T, NDIMS> {tet[0], tet[1], tet[2], tet[3]};
}

/*!
 * \brief Creates a bounding box around a Polygon
 *
 * \param [in] poly The polygon
 */
template <typename T, int NDIMS, PolygonArray ARRAY_TYPE = PolygonArray::Dynamic, int MAX_VERTS = DEFAULT_MAX_NUM_VERTICES>
AXOM_HOST_DEVICE BoundingBox<T, NDIMS> compute_bounding_box(
  const Polygon<T, NDIMS, ARRAY_TYPE, MAX_VERTS> &poly)
{
  BoundingBox<T, NDIMS> res;
  for(int i = 0; i < poly.numVertices(); ++i)
  {
    res.addPoint(poly[i]);
  }
  return res;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_COMPUTE_BOUNDING_BOX_HPP_
