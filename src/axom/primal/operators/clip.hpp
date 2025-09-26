// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file clip.hpp
 *
 * \brief Consists of a set of methods to clip a primal primitive against
 *        another primal primitive
 */

#ifndef AXOM_PRIMAL_CLIP_HPP_
#define AXOM_PRIMAL_CLIP_HPP_

#include "axom/core/utilities/Utilities.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Hexahedron.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/Polyhedron.hpp"
#include "axom/primal/geometry/Octahedron.hpp"
#include "axom/primal/geometry/Plane.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"

#include "axom/primal/operators/detail/clip_impl.hpp"

namespace axom
{
namespace primal
{
/*!
 * \brief Clips a 3D triangle against an axis-aligned bounding box in 3D
 *
 * \param [in] tri The triangle to clip
 * \param [in] bbox The bounding box to clip against
 * \return A planar polygon of the triangle clipped against the bounding box.
 *         If the triangle is completely outside the bounding box,
 *         the returned polygon is empty (i.e. it has no vertices).
 *
 * \note Using a specialization of the Sutherland-Hodgeman clipping algorithm
 *       for axis aligned planes
 */
template <typename T>
Polygon<T, 3> clip(const Triangle<T, 3>& tri, const BoundingBox<T, 3>& bbox)
{
  using BoundingBoxType = BoundingBox<T, 3>;
  using PolygonType = Polygon<T, 3>;

  // Use two polygons with pointers for 'back-buffer'-like swapping
  constexpr int MAX_VERTS = 6;
  PolygonType poly[2] = {PolygonType(MAX_VERTS), PolygonType(MAX_VERTS)};
  PolygonType* currentPoly = &poly[0];
  PolygonType* prevPoly = &poly[1];

  // First check if the triangle is contained in the bbox, if not we are empty
  BoundingBoxType triBox {tri[0], tri[1], tri[2]};

  if(!bbox.intersectsWith(triBox))
  {
    return *currentPoly;  // note: currentPoly is empty
  }

  // Set up the initial polygon
  currentPoly->addVertex(tri[0]);
  currentPoly->addVertex(tri[1]);
  currentPoly->addVertex(tri[2]);

  // If all the vertices are contained, we have no work to do
  if(bbox.contains(triBox))
  {
    return *currentPoly;  // Note current poly has the same verts as tri
  }

  // Loop through the planes of the bbox and clip the vertices
  for(int dim = 0; dim < 3; ++dim)
  {
    // Optimization note: we should be able to save some work based on
    // the clipping plane and the triangle's bounding box
    if(triBox.getMin()[dim] < bbox.getMin()[dim])
    {
      axom::utilities::swap(prevPoly, currentPoly);
      detail::clipAxisPlane(prevPoly, currentPoly, 2 * dim + 0, bbox.getMin()[dim]);
    }

    if(triBox.getMax()[dim] > bbox.getMax()[dim])
    {
      axom::utilities::swap(prevPoly, currentPoly);
      detail::clipAxisPlane(prevPoly, currentPoly, 2 * dim + 1, bbox.getMax()[dim]);
    }
  }

  return *currentPoly;
}

/*!
 * \brief Clips a 2D subject polygon against a clip polygon in 2D, returning
 *        their geometric intersection as a polygon
 *
 *  This function clips the subject polygon by the planes obtained from the
 *  clip polygon's edges (normals point inward). Clipping the
 *  subject polygon by each plane gives the polygon above that plane.
 *  Clipping the polygon by a plane involves
 *  finding new vertices at the intersection of the polygon edges and
 *  the plane, and removing vertices from the polygon that are below the
 *  plane.
 *
 *
 * \param [in] subjectPolygon The subject polygon
 * \param [in] clipPolygon The clip polygon
 * \param [in] eps The tolerance for plane point orientation.
 *                 Defaults to 1.e-10.
 * \param [in] tryFixOrientation If true, takes each shape with a negative
 *             signed area and swaps the order of some vertices in that
 *             shape to try to obtain a nonnegative signed area.
 *             Defaults to false.
 *
 * \return A polygon of the subject polygon clipped against the clip polygon.
 *
 * \note Function is based off the Sutherland–Hodgman algorithm.
 *
 * \warning Polygons with static array types must have enough vertices
 *          preallocated for the output polygon. It is mandatory that
 *          MAX_VERTS >= subjectPolygon.numVertices() + clipPolygon.numVertices()
 *          for the output polygon with the largest possible vertex count.
 *          Otherwise, if there is not enough preallocated vertices, output
 *          polygon will have missing vertices.
 *
 * \sa axom::primal::Polygon::addVertex(), axom::StaticArray::push_back()
 *     for behavior when there is not enough preallocated vertices.
 *
 * \warning tryFixOrientation flag does not guarantee the shapes' vertex orders
 *          will be valid. It is the responsiblity of the caller to pass
 *          shapes with a valid vertex order. Otherwise, if the shapes have
 *          invalid vertex orders, the returned Polygon
 *          will have a non-positive and/or unexpected area.
 *
 * \warning If tryFixOrientation flag is false and some of the shapes have
 *          a negative signed area, the returned Polygon
 *          will have a non-positive and/or unexpected area.
 *
 */
AXOM_SUPPRESS_HD_WARN
template <typename T, axom::primal::PolygonArray ARRAY_TYPE, int MAX_VERTS>
AXOM_HOST_DEVICE Polygon<T, 2, ARRAY_TYPE, MAX_VERTS> clip(
  const Polygon<T, 2, ARRAY_TYPE, MAX_VERTS>& subjectPolygon,
  const Polygon<T, 2, ARRAY_TYPE, MAX_VERTS>& clipPolygon,
  double eps = 1.e-10,
  bool tryFixOrientation = false)
{
  return detail::clipPolygonPolygon(subjectPolygon, clipPolygon, eps, tryFixOrientation);
}

/*!
 * \brief Clips a 2D subject polygon against a clip plane in 2D, returning
 *        their geometric intersection as a polygon.
 *
 *  This function makes a clip polygon line segment from the plane and
 *  delegates clipping to the polygon-polygon clipper.
 *
 *
 * \param [in] subjectPolygon The subject polygon
 * \param [in] clipPlane The clip plane
 * \param [in] eps The tolerance for plane point orientation.
 *                 Defaults to 1.e-10.
 * \param [in] tryFixOrientation If true, takes each shape with a negative
 *             signed area and swaps the order of some vertices in that
 *             shape to try to obtain a nonnegative signed area.
 *             Defaults to false.
 *
 * \return A polygon of the subject polygon clipped against the clip plane.
 */
template <typename T, axom::primal::PolygonArray ARRAY_TYPE, int MAX_VERTS>
AXOM_HOST_DEVICE Polygon<T, 2, ARRAY_TYPE, MAX_VERTS> clip(
  const Polygon<T, 2, ARRAY_TYPE, MAX_VERTS>& subjectPolygon,
  const Plane<T, 2>& clipPlane,
  double eps = 1.e-10,
  bool tryFixOrientation = false)
{
  return detail::clipPolygonPlane(subjectPolygon, clipPlane, eps, tryFixOrientation);
}

/*!
 * \brief Clips a 3D hexahedron against a tetrahedron in 3D, returning
 *        the geometric intersection of the hexahedron and the tetrahedron
 *        as a polyhedron
 *
 *  This function clips the hexahedron by the 4 planes obtained from the
 *  tetrahedron's faces (normals point inward). Clipping the
 *  hexahedron/polyhedron by each plane gives the polyhedron above that plane.
 *  Clipping the polyhedron by a plane involves
 *  finding new vertices at the intersection of the polyhedron edges and
 *  the plane, removing vertices from the polyhedron that are below the
 *  plane, and redefining the neighbors for each vertex (a vertex is a
 *  neighbor of another vertex if there is an edge between them).
 *
 *
 * \param [in] hex The hexahedron to clip
 * \param [in] tet The tetrahedron to clip against
 * \param [in] eps The epsilon value
 * \param [in] tryFixOrientation If true, takes each shape with a negative
 *             signed volume and swaps the order of some vertices in that
 *             shape to try to obtain a nonnegative signed volume.
 *             Defaults to false.
 *
 * \return A polyhedron of the hexahedron clipped against the tetrahedron.
 *
 * \note Function is based off clipPolyhedron() in Mike Owen's PolyClipper.
 *
 * \warning tryFixOrientation flag does not guarantee the shapes' vertex orders
 *          will be valid. It is the responsiblity of the caller to pass
 *          shapes with a valid vertex order. Otherwise, if the shapes have
 *          invalid vertex orders, the returned Polyhedron
 *          will have a non-positive and/or unexpected volume.
 *
 * \warning If tryFixOrientation flag is false and some of the shapes have
 *          a negative signed volume, the returned Polyhedron
 *          will have a non-positive and/or unexpected volume.
 *
 */
template <typename T>
AXOM_HOST_DEVICE Polyhedron<T, 3> clip(const Hexahedron<T, 3>& hex,
                                       const Tetrahedron<T, 3>& tet,
                                       double eps = 1.e-10,
                                       bool tryFixOrientation = false)
{
  return detail::clipHexahedron(hex, tet, eps, tryFixOrientation);
}

/*!
 * \brief Clips a 3D hexahedron against a tetrahedron in 3D, returning
 *        the geometric intersection of the hexahedron and the tetrahedron
 *        as a polyhedron
 *
 *  This function clips the hexahedron by the 4 planes obtained from the
 *  tetrahedron's faces (normals point inward). Clipping the
 *  hexahedron/polyhedron by each plane gives the polyhedron above that plane.
 *  Clipping the polyhedron by a plane involves
 *  finding new vertices at the intersection of the polyhedron edges and
 *  the plane, removing vertices from the polyhedron that are below the
 *  plane, and redefining the neighbors for each vertex (a vertex is a
 *  neighbor of another vertex if there is an edge between them).
 *
 *
 * \param [in] tet The tetrahedron to clip against
 * \param [in] hex The hexahedron to clip
 * \param [in] eps The epsilon value
 * \param [in] tryFixOrientation If true, takes each shape with a negative
 *             signed volume and swaps the order of some vertices in that
 *             shape to try to obtain a nonnegative signed volume.
 *             Defaults to false.
 *
 * \return A polyhedron of the hexahedron clipped against the tetrahedron.
 *
 * \note Function is based off clipPolyhedron() in Mike Owen's PolyClipper.
 *
 * \warning tryFixOrientation flag does not guarantee the shapes' vertex orders
 *          will be valid. It is the responsiblity of the caller to pass
 *          shapes with a valid vertex order. Otherwise, if the shapes have
 *          invalid vertex orders, the returned Polyhedron
 *          will have a non-positive and/or unexpected volume.
 *
 * \warning If tryFixOrientation flag is false and some of the shapes have
 *          a negative signed volume, the returned Polyhedron
 *          will have a non-positive and/or unexpected volume.
 *
 */
template <typename T>
AXOM_HOST_DEVICE Polyhedron<T, 3> clip(const Tetrahedron<T, 3>& tet,
                                       const Hexahedron<T, 3>& hex,
                                       double eps = 1.e-10,
                                       bool tryFixOrientation = false)
{
  return clip(hex, tet, eps, tryFixOrientation);
}

/*!
 * \brief Clips a 3D hexahedron against a hexahedron in 3D, returning
 *        the geometric intersection of the hexahedron and the hexahedron
 *        as a polyhedron
 *
 *  This function clips the hexahedron by the 6 planes obtained from the
 *  hexahedron's faces (normals point inward). Clipping the
 *  hexahedron/polyhedron by each plane gives the polyhedron above that plane.
 *  Clipping the polyhedron by a plane involves
 *  finding new vertices at the intersection of the polyhedron edges and
 *  the plane, removing vertices from the polyhedron that are below the
 *  plane, and redefining the neighbors for each vertex (a vertex is a
 *  neighbor of another vertex if there is an edge between them).
 *
 *
 * \param [in] hex1 The hexahedron to clip
 * \param [in] hex2 The hexahedron to clip against
 * \param [in] eps The epsilon value
 * \param [in] tryFixOrientation If true, takes each shape with a negative
 *             signed volume and swaps the order of some vertices in that
 *             shape to try to obtain a nonnegative signed volume.
 *             Defaults to false.
 *
 * \return A polyhedron of the hexahedron clipped against the hexahedron.
 *
 * \note Function is based off clipPolyhedron() in Mike Owen's PolyClipper.
 *
 * \note hex1 and hex2 are assumed to be convex.
 *
 * \warning tryFixOrientation flag does not guarantee the shapes' vertex orders
 *          will be valid. It is the responsiblity of the caller to pass
 *          shapes with a valid vertex order. Otherwise, if the shapes have
 *          invalid vertex orders, the returned Polyhedron
 *          will have a non-positive and/or unexpected volume.
 *
 * \warning If tryFixOrientation flag is false and some of the shapes have
 *          a negative signed volume, the returned Polyhedron
 *          will have a non-positive and/or unexpected volume.
 *
 */
template <typename T>
AXOM_HOST_DEVICE Polyhedron<T, 3> clip(const Hexahedron<T, 3>& hex1,
                                       const Hexahedron<T, 3>& hex2,
                                       double eps = 1.e-10,
                                       bool tryFixOrientation = false)
{
  return detail::clipHexahedron(hex1, hex2, eps, tryFixOrientation);
}

/*!
 * \brief Clips a 3D octahedron against a tetrahedron in 3D, returning
 *        the geometric intersection of the octahedron and the tetrahedron
 *        as a polyhedron
 *
 *  This function clips the octahedron by the 4 planes obtained from the
 *  tetrahedron's faces (normals point inward). Clipping the
 *  octahedron/polyhedron by each plane gives the polyhedron above that plane.
 *  Clipping the polyhedron by a plane involves
 *  finding new vertices at the intersection of the polyhedron edges and
 *  the plane, removing vertices from the polyhedron that are below the
 *  plane, and redefining the neighbors for each vertex (a vertex is a
 *  neighbor of another vertex if there is an edge between them).
 *
 *
 * \param [in] oct The octahedron to clip
 * \param [in] tet The tetrahedron to clip against
 * \param [in] eps The epsilon value
 * \param [in] tryFixOrientation If true, takes each shape with a negative
 *             signed volume and swaps the order of some vertices in that
 *             shape to try to obtain a nonnegative signed volume.
 *             Defaults to false.
 *
 * \return A polyhedron of the octahedron clipped against the tetrahedron.
 *
 * \note Function is based off clipPolyhedron() in Mike Owen's PolyClipper.
 *
 * \warning tryFixOrientation flag does not guarantee the shapes' vertex orders
 *          will be valid. It is the responsiblity of the caller to pass
 *          shapes with a valid vertex order. Otherwise, if the shapes have
 *          invalid vertex orders, the returned Polyhedron
 *          will have a non-positive and/or unexpected volume.
 *
 * \warning If tryFixOrientation flag is false and some of the shapes have
 *          a negative signed volume, the returned Polyhedron
 *          will have a non-positive and/or unexpected volume.
 *
 */
template <typename T>
AXOM_HOST_DEVICE Polyhedron<T, 3> clip(const Octahedron<T, 3>& oct,
                                       const Tetrahedron<T, 3>& tet,
                                       double eps = 1.e-10,
                                       bool tryFixOrientation = false)
{
  return detail::clipOctahedron(oct, tet, eps, tryFixOrientation);
}

/*!
 * \brief Clips a 3D octahedron against a tetrahedron in 3D, returning
 *        the geometric intersection of the octahedron and the tetrahedron
 *        as a polyhedron
 *
 *  This function clips the octahedron by the 4 planes obtained from the
 *  tetrahedron's faces (normals point inward). Clipping the
 *  octahedron/polyhedron by each plane gives the polyhedron above that plane.
 *  Clipping the polyhedron by a plane involves
 *  finding new vertices at the intersection of the polyhedron edges and
 *  the plane, removing vertices from the polyhedron that are below the
 *  plane, and redefining the neighbors for each vertex (a vertex is a
 *  neighbor of another vertex if there is an edge between them).
 *
 *
 * \param [in] tet The tetrahedron to clip against
 * \param [in] oct The octahedron to clip
 * \param [in] tet The tetrahedron to clip against
 * \param [in] eps The epsilon value
 * \param [in] tryFixOrientation If true, takes each shape with a negative
 *             signed volume and swaps the order of some vertices in that
 *             shape to try to obtain a nonnegative signed volume.
 *             Defaults to false.
 *
 * \return A polyhedron of the octahedron clipped against the tetrahedron.
 *
 * \note Function is based off clipPolyhedron() in Mike Owen's PolyClipper.
 *
 * \warning tryFixOrientation flag does not guarantee the shapes' vertex orders
 *          will be valid. It is the responsiblity of the caller to pass
 *          shapes with a valid vertex order. Otherwise, if the shapes have
 *          invalid vertex orders, the returned Polyhedron
 *          will have a non-positive and/or unexpected volume.
 *
 * \warning If tryFixOrientation flag is false and some of the shapes have
 *          a negative signed volume, the returned Polyhedron
 *          will have a non-positive and/or unexpected volume.
 *
 */
template <typename T>
AXOM_HOST_DEVICE Polyhedron<T, 3> clip(const Tetrahedron<T, 3>& tet,
                                       const Octahedron<T, 3>& oct,
                                       double eps = 1.e-10,
                                       bool tryFixOrientation = false)
{
  return clip(oct, tet, eps, tryFixOrientation);
}

/*!
 * \brief Clips a 3D tetrahedron against another tetrahedron in 3D, returning
 *        the geometric intersection as a polyhedron
 *
 *  This function clips the first tetrahedron by the 4 planes obtained from the
 *  second tetrahedron's faces (normals point inward). Clipping the
 *  tetrahedron/polyhedron by each plane gives the polyhedron above that plane.
 *  Clipping the polyhedron by a plane involves
 *  finding new vertices at the intersection of the polyhedron edges and
 *  the plane, removing vertices from the polyhedron that are below the
 *  plane, and redefining the neighbors for each vertex (a vertex is a
 *  neighbor of another vertex if there is an edge between them).
 *
 *
 * \param [in] tet1 The tetrahedron to clip
 * \param [in] tet2 The tetrahedron to clip against
 * \param [in] eps The epsilon value
 * \param [in] tryFixOrientation If true, takes each shape with a negative
 *             signed volume and swaps the order of some vertices in that
 *             shape to try to obtain a nonnegative signed volume.
 *             Defaults to false.
 *
 * \return A polyhedron of the tetrahedron clipped against
 *         the other tetrahedron.
 *
 * \note Function is based off clipPolyhedron() in Mike Owen's PolyClipper.
 *
 * \warning tryFixOrientation flag does not guarantee the shapes' vertex orders
 *          will be valid. It is the responsiblity of the caller to pass
 *          shapes with a valid vertex order. Otherwise, if the shapes have
 *          invalid vertex orders, the returned Polyhedron
 *          will have a non-positive and/or unexpected volume.
 *
 * \warning If tryFixOrientation flag is false and some of the shapes have
 *          a negative signed volume, the returned Polyhedron
 *          will have a non-positive and/or unexpected volume.
 *
 */
template <typename T>
AXOM_HOST_DEVICE Polyhedron<T, 3> clip(const Tetrahedron<T, 3>& tet1,
                                       const Tetrahedron<T, 3>& tet2,
                                       double eps = 1.e-10,
                                       bool tryFixOrientation = false)
{
  return detail::clipTetrahedron(tet1, tet2, eps, tryFixOrientation);
}

/*!
 * \brief Clips a 3D polyhedron against another polyhedron in 3D, returning
 *        the geometric intersection as a polyhedron
 *
 *  This function clips the first polyhedron by the planes obtained from the
 *  second polyhedron's faces (normals point inward). Clipping the
 *  polyhedron by each plane gives the polyhedron above that plane.
 *  Clipping the polyhedron by a plane involves
 *  finding new vertices at the intersection of the polyhedron edges and
 *  the plane, removing vertices from the polyhedron that are below the
 *  plane, and redefining the neighbors for each vertex (a vertex is a
 *  neighbor of another vertex if there is an edge between them).
 *
 *
 * \param [in] tet The tetrahedron to clip
 * \param [in] poly The polyhedron to clip against
 * \param [in] eps The epsilon value
 * \param [in] tryFixOrientation If true, takes each shape with a negative
 *             signed volume and swaps the order of some vertices in that
 *             shape to try to obtain a nonnegative signed volume.
 *             Defaults to false.
 *
 * \return A polyhedron of the clip results
 *
 * \note poly is assumed to be convex. Any non-planar faces in poly will result
 *       in multiple clipping planes for the face.
 *
 * \warning tryFixOrientation flag does not guarantee the shapes' vertex orders
 *          will be valid. It is the responsiblity of the caller to pass
 *          shapes with a valid vertex order. Otherwise, if the shapes have
 *          invalid vertex orders, the returned Polyhedron
 *          will have a non-positive and/or unexpected volume.
 *
 * \warning If tryFixOrientation flag is false and some of the shapes have
 *          a negative signed volume, the returned Polyhedron
 *          will have a non-positive and/or unexpected volume.
 */
template <typename T>
AXOM_HOST_DEVICE Polyhedron<T, 3> clip(const Tetrahedron<T, 3>& tet,
                                       const Polyhedron<T, 3>& poly,
                                       double eps = 1.e-10,
                                       bool tryFixOrientation = false)
{
  return detail::clipTetrahedron(tet, poly, eps, tryFixOrientation);
}

/*!
 * \brief Clips a 3D polyhedron against another polyhedron in 3D, returning
 *        the geometric intersection as a polyhedron
 *
 *  This function clips the first polyhedron by the planes obtained from the
 *  second polyhedron's faces (normals point inward). Clipping the
 *  polyhedron by each plane gives the polyhedron above that plane.
 *  Clipping the polyhedron by a plane involves
 *  finding new vertices at the intersection of the polyhedron edges and
 *  the plane, removing vertices from the polyhedron that are below the
 *  plane, and redefining the neighbors for each vertex (a vertex is a
 *  neighbor of another vertex if there is an edge between them).
 *
 *
 * \param [in] tet The tetrahedron to clip
 * \param [in] poly The polyhedron to clip against
 * \param [in] eps The epsilon value
 * \param [in] tryFixOrientation If true, takes each shape with a negative
 *             signed volume and swaps the order of some vertices in that
 *             shape to try to obtain a nonnegative signed volume.
 *             Defaults to false.
 *
 * \return A polyhedron of the clip results
 *
 * \note poly is assumed to be convex. Any non-planar faces in poly will result
 *       in multiple clipping planes for the face.
 *
 * \warning tryFixOrientation flag does not guarantee the shapes' vertex orders
 *          will be valid. It is the responsiblity of the caller to pass
 *          shapes with a valid vertex order. Otherwise, if the shapes have
 *          invalid vertex orders, the returned Polyhedron
 *          will have a non-positive and/or unexpected volume.
 *
 * \warning If tryFixOrientation flag is false and some of the shapes have
 *          a negative signed volume, the returned Polyhedron
 *          will have a non-positive and/or unexpected volume.
 */
template <typename T>
AXOM_HOST_DEVICE Polyhedron<T, 3> clip(const Polyhedron<T, 3>& poly,
                                       const Tetrahedron<T, 3>& tet,
                                       double eps = 1.e-10,
                                       bool tryFixOrientation = false)
{
  return detail::clipTetrahedron(tet, poly, eps, tryFixOrientation);
}

/*!
 * \brief Clips a 3D hexahedron against a polyhedron in 3D, returning
 *        the geometric intersection as a polyhedron
 *
 *  This function clips the polyhedron by the planes obtained from the
 *  hexahedron's faces (normals point inward). Clipping the
 *  polyhedron by each plane gives the polyhedron above that plane.
 *  Clipping the polyhedron by a plane involves
 *  finding new vertices at the intersection of the polyhedron edges and
 *  the plane, removing vertices from the polyhedron that are below the
 *  plane, and redefining the neighbors for each vertex (a vertex is a
 *  neighbor of another vertex if there is an edge between them).
 *
 * \param [in] hex The hexahedron to clip
 * \param [in] poly The polyhedron to clip against
 * \param [in] eps The epsilon value
 * \param [in] tryFixOrientation If true, takes each shape with a negative
 *             signed volume and swaps the order of some vertices in that
 *             shape to try to obtain a nonnegative signed volume.
 *             Defaults to false.
 *
 * \return A polyhedron of the clip results
 *
 * \note \a hex and \a poly shapes are assumed to be convex. The \a hex shape
 *       is clipped as a polyhedron. Any non-planar faces may give rise to
 *       non-planar clipped faces. The \a poly shape's faces are used as
 *       clipping planes for \a hex. Any non-planar faces in \a poly will
 *       result in multiple clipping planes for the face.
 *
 * \warning tryFixOrientation flag does not guarantee the shapes' vertex orders
 *          will be valid. It is the responsiblity of the caller to pass
 *          shapes with a valid vertex order. Otherwise, if the shapes have
 *          invalid vertex orders, the returned Polyhedron
 *          will have a non-positive and/or unexpected volume.
 *
 * \warning If tryFixOrientation flag is false and some of the shapes have
 *          a negative signed volume, the returned Polyhedron
 *          will have a non-positive and/or unexpected volume.
 */
template <typename T>
AXOM_HOST_DEVICE Polyhedron<T, 3> clip(const Hexahedron<T, 3>& hex,
                                       const Polyhedron<T, 3>& poly,
                                       double eps = 1.e-10,
                                       bool tryFixOrientation = false)
{
  return detail::clipHexahedron(hex, poly, eps, tryFixOrientation);
}

/*!
 * \brief Clips a 3D polyhedron against a hexahedron in 3D, returning
 *        the geometric intersection as a polyhedron
 *
 *  This function clips the polyhedron by the planes obtained from the
 *  hexahedron's faces (normals point inward). Clipping the
 *  polyhedron by each plane gives the polyhedron above that plane.
 *  Clipping the polyhedron by a plane involves
 *  finding new vertices at the intersection of the polyhedron edges and
 *  the plane, removing vertices from the polyhedron that are below the
 *  plane, and redefining the neighbors for each vertex (a vertex is a
 *  neighbor of another vertex if there is an edge between them).
 *
 * \param [in] poly The polyhedron to clip against
 * \param [in] hex The hexahedron to clip
 * \param [in] eps The epsilon value
 * \param [in] tryFixOrientation If true, takes each shape with a negative
 *             signed volume and swaps the order of some vertices in that
 *             shape to try to obtain a nonnegative signed volume.
 *             Defaults to false.
 *
 * \return A polyhedron of the clip results
 *
 * \note \a hex and \a poly shapes are assumed to be convex. Any non-planar
 *       faces may give rise to non-planar clipped faces. The \a hex shape's
 *       faces are used as clipping planes for \a poly. Any non-planar faces
 *       in \a hex will result in multiple clipping planes for the face.
 *
 * \warning tryFixOrientation flag does not guarantee the shapes' vertex orders
 *          will be valid. It is the responsiblity of the caller to pass
 *          shapes with a valid vertex order. Otherwise, if the shapes have
 *          invalid vertex orders, the returned Polyhedron
 *          will have a non-positive and/or unexpected volume.
 *
 * \warning If tryFixOrientation flag is false and some of the shapes have
 *          a negative signed volume, the returned Polyhedron
 *          will have a non-positive and/or unexpected volume.
 */
template <typename T>
AXOM_HOST_DEVICE Polyhedron<T, 3> clip(const Polyhedron<T, 3>& poly,
                                       const Hexahedron<T, 3>& hex,
                                       double eps = 1.e-10,
                                       bool tryFixOrientation = false)
{
  return detail::clipHexahedron(hex, poly, eps, tryFixOrientation);
}

/*!
 * \brief Clips a 3D polyhedron against another polyhedron in 3D, returning
 *        the geometric intersection as a polyhedron
 *
 *  This function clips the first polyhedron by the planes obtained from the
 *  second polyhedron's faces (normals point inward). Clipping the
 *  polyhedron by each plane gives the polyhedron above that plane.
 *  Clipping the polyhedron by a plane involves
 *  finding new vertices at the intersection of the polyhedron edges and
 *  the plane, removing vertices from the polyhedron that are below the
 *  plane, and redefining the neighbors for each vertex (a vertex is a
 *  neighbor of another vertex if there is an edge between them).
 *
 *
 * \param [in] poly1 The polyhedron to clip
 * \param [in] poly2 The polyhedron to clip against
 * \param [in] eps The epsilon value
 *
 * \return A polyhedron of the clip results
 *
 * \note \a poly1 and \a poly2 shapes are assumed to be convex. Any non-planar
 *       faces may give rise to non-planar clipped faces. The \a poly2 shape's
 *       faces are used as clipping planes for \a poly1. Any non-planar faces
 *       in \a poly2 will result in multiple clipping planes for the face.
 */
template <typename T>
AXOM_HOST_DEVICE Polyhedron<T, 3> clip(const Polyhedron<T, 3>& poly1,
                                       const Polyhedron<T, 3>& poly2,
                                       double eps = 1.e-10)
{
  return detail::clipPolyhedronPolyhedron(poly1, poly2, eps);
}

/*!
 * \brief Clips a 3D tetrahedron against the half-space defined by a plane
 *        and returns the resulting polyhedron
 *
 *  This function clips a tetrahedron against the half-space defined by a
 *  plane. This involves finding new vertices at the intersection of the
 *  polyhedron edges and the plane, removing vertices from the polyhedron
 *  that are below the plane, and redefining the neighbors for each vertex
 *  (a vertex is a neighbor of another vertex if there is an edge between
 *  them).
 *
 * \param [in] tet The tetrahedron to clip
 * \param [in] plane The plane defining the half-space used to clip the tetrahedron
 * \param [in] eps The tolerance for plane point orientation
 * \param [in] tryFixOrientation If true and the tetrahedron has a negative
 *             signed volume, swaps the order of some vertices in the
 *             tetrathedron to try to obtain a nonnegative signed volume.
 *             Defaults to false.
 *
 * \return The polyhedron obtained from clipping a tetrahedron against
 *         the half-space defined by a plane.
 *
 * \note Function is based off clipPolyhedron() in Mike Owen's PolyClipper.
 *
 * \warning tryFixOrientation flag does not guarantee the tetrahedron's vertex
 *          order will be valid. It is the responsiblity of the caller to pass
 *          a tetrahedron with a valid vertex order. Otherwise, the returned
 *          polyhedron will have a non-positive and/or unexpected volume.
 *
 * \warning If the tryFixOrientation flag is false and the tetrahedron has
 *          a negative signed volume, the returned polyhedron will have a
 *          non-positive and/or unexpected volume.
 */
template <typename T>
AXOM_HOST_DEVICE Polyhedron<T, 3> clip(const Tetrahedron<T, 3>& tet,
                                       const Plane<T, 3>& plane,
                                       double eps = 1.e-10,
                                       bool tryFixOrientation = false)
{
  return detail::clipTetrahedron(tet, plane, eps, tryFixOrientation);
}

/*!
 * \brief Clips a 3D tetrahedron against the half-space defined by a plane
 *        and returns the resulting polyhedron
 *
 *  This function clips a tetrahedron against the half-space defined by a
 *  plane. This involves finding new vertices at the intersection of the
 *  polyhedron edges and the plane, removing vertices from the polyhedron
 *  that are below the plane, and redefining the neighbors for each vertex
 *  (a vertex is a neighbor of another vertex if there is an edge between
 *  them).
 *
 * \param [in] plane The plane defining the half-space used to clip the tetrahedron
 * \param [in] tet The tetrahedron to clip
 * \param [in] eps The tolerance for plane point orientation
 * \param [in] tryFixOrientation If true and the tetrahedron has a negative
 *             signed volume, swaps the order of some vertices in the
 *             tetrathedron to try to obtain a nonnegative signed volume.
 *             Defaults to false.
 *
 * \return The polyhedron obtained from clipping a tetrahedron against
 *         the half-space defined by a plane.
 *
 * \note Function is based off clipPolyhedron() in Mike Owen's PolyClipper.
 *
 * \warning tryFixOrientation flag does not guarantee the tetrahedron's vertex
 *          order will be valid. It is the responsiblity of the caller to pass
 *          a tetrahedron with a valid vertex order. Otherwise, the returned
 *          polyhedron will have a non-positive and/or unexpected volume.
 *
 * \warning If the tryFixOrientation flag is false and the tetrahedron has
 *          a negative signed volume, the returned polyhedron will have a
 *          non-positive and/or unexpected volume.
 */
template <typename T>
AXOM_HOST_DEVICE Polyhedron<T, 3> clip(const Plane<T, 3>& plane,
                                       const Tetrahedron<T, 3>& tet,
                                       double eps = 1.e-10,
                                       bool tryFixOrientation = false)
{
  return clip(tet, plane, eps, tryFixOrientation);
}

/*!
 * \brief Clips a 3D polyhedron against the half-space defined by a plane
 *        and returns the resulting polyhedron
 *
 *  This function clips a polyhedron against the half-space defined by a
 *  plane. This involves finding new vertices at the intersection of the
 *  polyhedron edges and the plane, removing vertices from the polyhedron
 *  that are below the plane, and redefining the neighbors for each vertex
 *  (a vertex is a neighbor of another vertex if there is an edge between
 *  them).
 *
 * \param [in] plane The plane defining the half-space used to clip the tetrahedron
 * \param [in] poly The polyhedron to clip
 * \param [in] eps The tolerance for plane point orientation
 *
 * \return The polyhedron obtained from clipping a polyhedron against
 *         the half-space defined by a plane.
 */
template <typename T>
AXOM_HOST_DEVICE Polyhedron<T, 3> clip(const Polyhedron<T, 3>& poly,
                                       const Plane<T, 3>& plane,
                                       double eps = 1.e-10)
{
  Polyhedron<T, 3> clipped(poly);
  detail::clipPolyhedron(clipped, plane, eps);
  return clipped;
}

/*!
 * \brief Clips a 3D polyhedron against the half-space defined by a plane
 *        and returns the resulting polyhedron
 *
 *  This function clips a polyhedron against the half-space defined by a
 *  plane. This involves finding new vertices at the intersection of the
 *  polyhedron edges and the plane, removing vertices from the polyhedron
 *  that are below the plane, and redefining the neighbors for each vertex
 *  (a vertex is a neighbor of another vertex if there is an edge between
 *  them).
 *
 * \param [in] poly The polyhedron to clip
 * \param [in] plane The plane defining the half-space used to clip the polyhedron.
 * \param [in] eps The tolerance for plane point orientation
 *
 * \return The polyhedron obtained from clipping a polyhedron against
 *         the half-space defined by a plane.
 */
template <typename T>
AXOM_HOST_DEVICE Polyhedron<T, 3> clip(const Plane<T, 3>& plane,
                                       const Polyhedron<T, 3>& poly,
                                       double eps = 1.e-10)
{
  Polyhedron<T, 3> clipped(poly);
  detail::clipPolyhedron(clipped, plane, eps);
  return clipped;
}

/*!
 * \brief Clips a 3D hexahedron against the half-space defined by a plane
 *        and returns the resulting polyhedron
 *
 *  This function clips a hexahedron against the half-space defined by a
 *  plane. This involves finding new vertices at the intersection of the
 *  hexahedron edges and the plane, removing vertices from the hexahedron
 *  that are below the plane, and redefining the neighbors for each vertex
 *  (a vertex is a neighbor of another vertex if there is an edge between
 *  them).
 *
 * \param [in] hex The hexahedron to clip
 * \param [in] plane The plane defining the half-space used to clip the hexahedron.
 * \param [in] eps The tolerance for plane point orientation
 * \param [in] tryFixOrientation If true, takes each shape with a negative
 *             signed volume and swaps the order of some vertices in that
 *             shape to try to obtain a nonnegative signed volume.
 *             Defaults to false.
 *
 * \return The polyhedron obtained from clipping a hexahedron against
 *         the half-space defined by a plane.
 */
template <typename T>
AXOM_HOST_DEVICE Polyhedron<T, 3> clip(const Hexahedron<T, 3>& hex,
                                       const Plane<T, 3>& plane,
                                       double eps = 1.e-10,
                                       bool tryFixOrientation = false)
{
  // Initialize our polyhedron to clip
  auto poly = Polyhedron<T, 3>::from_primitive(hex, tryFixOrientation);
  return clip(poly, plane, eps);
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_CLIP_HPP_
