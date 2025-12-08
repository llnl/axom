// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file clip_impl.hpp
 *
 * \brief Helper functions for the primal clipping operators
 */

#ifndef AXOM_PRIMAL_CLIP_IMPL_HPP_
#define AXOM_PRIMAL_CLIP_IMPL_HPP_

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Hexahedron.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/Polyhedron.hpp"
#include "axom/primal/geometry/Octahedron.hpp"
#include "axom/primal/geometry/Plane.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"

#include "axom/primal/operators/intersect.hpp"
#include "axom/primal/operators/orientation.hpp"

namespace axom
{
namespace primal
{
namespace detail
{
/*! Returns true when index is even */
inline bool isEven(int index) { return (index & 1) == 0; }

/*!
 * \brief Specialized point plane classifier for axis aligned planes
 *
 * \param [in] pt The plane to classify
 * \param [in] index The index of the axis aligned plane. See below for mapping
 * \param [in] val The plane's coordinate with respect to the given axis
 * \param [in] eps A parameter for thickening width of the plane (default 1e-8)
 *
 * Mapping of index to axis
 * * 0 -> -x axis
 * * 1 -> +x axis
 * * 2 -> -y axis
 * * 3 -> +y axis
 * * 4 -> -z axis
 * * 5 -> +z axis
 *
 * \return An OrientedSide (ON_POSITIVE_SIDE, ON_NEGATIVE_SIDE, ON_BOUNDARY)
 *         value based on the relative orientations of point pt and the
 *         corresponding plane associated with index.
 * \see OrientedPlane enum
 */
template <typename T, int NDIMS>
int classifyPointAxisPlane(const Point<T, NDIMS>& pt, int index, T val, const double eps = 1e-8)
{
  // Note: we are exploiting the fact that the planes are axis aligned
  // So the dot product is +/- the given coordinate.
  // In general, we would need to call distance(pt, plane) here
  const T dist = isEven(index) ? val - pt[index / 2] : pt[index / 2] - val;

  if(dist > eps)
  {
    return ON_POSITIVE_SIDE;
  }
  if(dist < -eps)
  {
    return ON_NEGATIVE_SIDE;
  }

  return ON_BOUNDARY;
}

/*!
 * \brief Finds the clipping intersection point between points a and b.
 *
 * \param [in] a The point behind the plane
 * \param [in] b The point in front of the plane
 * \param [in] index The index of the axis aligned plane.
 * \param [in] val The plane's coordinate with respect to the given axis
 * \return The point between a and b whose corresponding coordinate is val
 *
 * \see classifyPointAxisPlane for description of how index maps to coordinates.
 */
template <typename T, int NDIMS>
Point<T, NDIMS> findIntersectionPoint(const Point<T, NDIMS>& a,
                                      const Point<T, NDIMS>& b,
                                      int index,
                                      T val)
{
  using PointType = Point<T, NDIMS>;

  // Need to find a parameter t for the point pt, such that,
  // * 0 <= t <= 1
  // * pt = a + t * (b-a)
  // * pt[ index/2]  == val

  const int coord = index / 2;
  const T t = (val - a[coord]) / (b[coord] - a[coord]);
  SLIC_ASSERT(0. <= t && t <= 1.);

  const auto ret = PointType::lerp(a, b, t);
  SLIC_ASSERT(classifyPointAxisPlane(ret, index, val) == ON_BOUNDARY);

  return ret;
}

/*!
 * \brief Clips the vertices of the polygon to be behind the plane.
 *
 * This is a specialization of the Sutherland-Hodgeman clipping algorithm
 * for axis-aligned planes
 *
 * \param [in] prevPoly  An input polygon with the vertices to clip
 * \param [out] currentPoly An output polygon whose coordinates are clipped
 *                          against this plane.
 * \param [in] index The index of the axis aligned plane.
 * \param [in] val The plane's coordinate with respect to the given axis
 *
 * \note Algorithm for robust clipping against "thick" planes derived from
 *       Section 8.3 of Christer Ericson's "Real-Time Collision Detection"
 *       and is based on the Sutherland-Hodgeman clipping algorithm.
 *       We are only keeping the "back" polygon, w.r.t. that algorithm.
 * \see classifyPointAxisPlane for description of how index maps to coordinates.
 */
template <typename T, int NDIMS>
void clipAxisPlane(const Polygon<T, NDIMS>* prevPoly, Polygon<T, NDIMS>* currentPoly, int index, T val)
{
  using PointType = Point<T, NDIMS>;

  currentPoly->clear();
  const int numVerts = prevPoly->numVertices();

  if(numVerts == 0)
  {
    return;
  }

  // Initialize point a with the last vertex of the polygon
  const PointType* a = &(*prevPoly)[numVerts - 1];
  int aSide = classifyPointAxisPlane(*a, index, val);

  for(int i = 0; i < numVerts; ++i)
  {
    const PointType* b = &(*prevPoly)[i];
    const int bSide = classifyPointAxisPlane(*b, index, val);

    switch(bSide)
    {
    case ON_POSITIVE_SIDE:
      if(aSide == ON_NEGATIVE_SIDE)
      {
        currentPoly->addVertex(findIntersectionPoint(*a, *b, index, val));
      }
      break;
    case ON_BOUNDARY:
      if(aSide == ON_NEGATIVE_SIDE)
      {
        currentPoly->addVertex(*b);
      }
      break;
    case ON_NEGATIVE_SIDE:
      switch(aSide)
      {
      case ON_POSITIVE_SIDE:
        currentPoly->addVertex(findIntersectionPoint(*a, *b, index, val));
        currentPoly->addVertex(*b);
        break;
      case ON_BOUNDARY:
        currentPoly->addVertex(*a);
        currentPoly->addVertex(*b);
        break;
      case ON_NEGATIVE_SIDE:
        currentPoly->addVertex(*b);
        break;
      }
      break;
    }

    // swap a and b
    a = b;
    aSide = bSide;
  }
}

template <typename T, int NDIMS>
AXOM_HOST_DEVICE void poly_clip_vertices(Polyhedron<T, NDIMS>& poly,
                                         const Plane<T, NDIMS>& plane,
                                         const double eps,
                                         unsigned int& out_clipped)
{
  using SegmentType = Segment<T, NDIMS>;

  // Loop over Polyhedron vertices
  int numVerts = poly.numVertices();
  for(std::int8_t i = 0; i < numVerts; i++)
  {
    int orientation = plane.getOrientation(poly[i], eps);

    // Vertex is under the plane
    if(orientation == ON_NEGATIVE_SIDE)
    {
      // Mark this vertex for removal later
      out_clipped |= 1 << i;

      // Check neighbors for vertex above the plane (edge clipped by plane)
      int numNeighbors = poly.getNumNeighbors(i);
      for(int j = 0; j < numNeighbors; j++)
      {
        std::int8_t neighborIndex = poly.getNeighbors(i)[j];

        int neighborOrientation = plane.getOrientation(poly[neighborIndex], eps);

        // Insert new vertex to polyhedron, where edge intersects plane.
        if(neighborOrientation == ON_POSITIVE_SIDE)
        {
          const int expectedVertexIndex = poly.numVertices();
          AXOM_UNUSED_VAR(expectedVertexIndex);  // silence warning in release configs

          T lerp_val;
          SegmentType seg(poly[i], poly[neighborIndex]);
          intersect(plane, seg, lerp_val);

          int newVertexIndex = poly.addVertex(seg.at(lerp_val));
          // A probable AMD GPU 6.3 compiler bug caused the first assert to fail.
          // The second assert is equivalent because we expect addVertex call to
          // consume the next index and increment the number of vertices.
          // However, when the compilers no longar cause false failures, we should
          // restore the first assert.
          // SLIC_ASSERT(newVertexIndex == expectedVertexIndex);
          SLIC_ASSERT(newVertexIndex == poly.numVertices() - 1);

          poly.addNeighbors(newVertexIndex, {i, neighborIndex});

          // Update current vertex's & neighbor's neighbors with the
          // new vertex
          poly.getNeighbors(i)[j] = newVertexIndex;
          for(int k = 0; k < poly.getNumNeighbors(neighborIndex); k++)
          {
            if(poly.getNeighbors(neighborIndex)[k] == i)
            {
              poly.getNeighbors(neighborIndex)[k] = newVertexIndex;
            }
          }
        }
      }
    }
  }  // end of loop over Polyhedron vertices
}

template <typename T, int NDIMS>
AXOM_HOST_DEVICE void poly_clip_fix_nbrs(Polyhedron<T, NDIMS>& poly,
                                         const Plane<T, NDIMS>& plane,
                                         const int oldVerts,
                                         const double eps,
                                         const unsigned int clipped)
{
  NeighborCollection& poly_nbrs = poly.getNeighbors();
  // Keep copy of old connectivity
  NeighborCollection old_nbrs = poly.getNeighbors();
  for(int i = 0; i < poly.numVertices(); i++)
  {
    // Check clipped created vertices first, then vertices on the plane
    int vIndex = (i + oldVerts) % poly.numVertices();
    int vOrientation = plane.getOrientation(poly[vIndex], eps);

    if(vIndex >= oldVerts || vOrientation == ON_BOUNDARY)
    {
      for(int j = 0; j < poly_nbrs.getNumNeighbors(vIndex); j++)
      {
        int neighborIndex = poly_nbrs[vIndex][j];
        int neighborOrientation = plane.getOrientation(poly[neighborIndex], eps);

        // This neighbor is not newly inserted and below the plane
        if(neighborIndex < oldVerts && neighborOrientation == ON_NEGATIVE_SIDE)
        {
          // Look for 1st vertex along this face not below the plane.
          int iprev = vIndex;
          int inext = neighborIndex;
          int itmp = inext;

          int val = 0;

          while((clipped & (1 << inext)) && (val++ < poly.numVertices()))
          {
            itmp = inext;
            unsigned int next_nbrs = poly_nbrs.getNumNeighbors(inext);
            // Find next vertex along face
            for(unsigned int ni = 0; ni < next_nbrs; ni++)
            {
              if(poly_nbrs[inext][ni] == iprev)
              {
                inext = (ni == 0) ? poly_nbrs[inext][next_nbrs - 1] : poly_nbrs[inext][ni - 1];
                break;
              }
            }

            iprev = itmp;
          }

          // Remove neighbor from list if vertex found was already a neighbor or
          // is the vertex we are currently checking for.
          if(poly_nbrs[vIndex][(j + 1) % poly_nbrs.getNumNeighbors(vIndex)] == inext ||
             inext == vIndex)
          {
            poly_nbrs[vIndex][j] = -1;
          }

          // Otherwise update neighbor lists of vertex found and vertex we are checking for.
          else
          {
            poly_nbrs[vIndex][j] = inext;

            if(inext >= oldVerts)
            {
              poly_nbrs.insertNeighborAtPos(inext, vIndex, 0);
              old_nbrs.insertNeighborAtPos(inext, -1, 0);
            }
            else
            {
              int offset {};
              for(int oi = 0; oi < old_nbrs.getNumNeighbors(inext); oi++)
              {
                if(old_nbrs[inext][oi] == iprev)
                {
                  offset = oi;
                  break;
                }

                // Max offset
                if(oi == old_nbrs.getNumNeighbors(inext) - 1)
                {
                  offset = old_nbrs.getNumNeighbors(inext);
                }
              }
              poly_nbrs.insertNeighborAtPos(inext, vIndex, offset);
              old_nbrs.insertNeighborAtPos(inext, vIndex, offset);
            }
          }
        }
      }  // end of loop over vertex neighbors
    }
  }  // end of loop over Polyhedron vertices
  poly.getNeighbors().pruneNeighbors();
}

template <typename T, int NDIMS>
AXOM_HOST_DEVICE void poly_clip_reindex(Polyhedron<T, NDIMS>& poly, const unsigned int clipped)
{
  // Dictionary for old indices to new indices positions
  std::int8_t newIndices[Polyhedron<T, NDIMS>::MAX_VERTS] = {0};

  Polyhedron<T, NDIMS> old_poly;

  for(int i = 0; i < poly.numVertices(); i++)
  {
    old_poly.addVertex(poly[i]);
    for(int j = 0; j < poly.getNumNeighbors(i); j++)
    {
      old_poly.addNeighbors(i, poly.getNeighbors(i)[j]);
    }
  }

  poly.clear();

  int curIndex = 0;

  for(int i = 0; i < old_poly.numVertices(); i++)
  {
    if(!(clipped & (1 << i)))
    {
      // Non-clipped vertex
      newIndices[i] = curIndex++;
      poly.addVertex(old_poly[i]);
    }
  }

  // Reinsert neighbors into polyhedron
  for(int i = 0; i < old_poly.numVertices(); i++)
  {
    if(!(clipped & (1 << i)))
    {
      for(int j = 0; j < old_poly.getNumNeighbors(i); j++)
      {
        poly.addNeighbors(newIndices[i], {newIndices[old_poly.getNeighbors()[i][j]]});
      }
    }
  }
}

/*!
 * \brief Clips a polyhedron against a half-space defined by a plane
 *
 * \param [inout] poly The polyhedron to clip
 * \param [in] plane The plane defining the half-space used to clip the polyhedron
 * \param [in] eps The tolerance for plane point orientation
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE void clipPolyhedron(Polyhedron<T, NDIMS>& poly,
                                     const Plane<T, NDIMS>& plane,
                                     double eps)
{
  using BoxType = BoundingBox<T, NDIMS>;

  // Check that plane intersects Polyhedron
  if(intersect(plane, BoxType(&poly[0], poly.numVertices()), true, eps))
  {
    int numVerts = poly.numVertices();

    // Each bit value indicates if that Polyhedron vertex is formed from
    // Polyhedron clipping with a plane.
    unsigned int clipped = 0;

    // Clip polyhedron against current plane, generating extra vertices
    // where edges meet the plane.
    poly_clip_vertices(poly, plane, eps, clipped);

    // Adjust connectivity to link up newly-generated vertices.
    poly_clip_fix_nbrs(poly, plane, numVerts, eps, clipped);

    // Reindex polyhedron connectivity by removing vertices on the negative
    // side of the plane.
    poly_clip_reindex(poly, clipped);
  }

  // If entire polyhedron is below a plane (points can be on the plane),
  // it is completely removed.
  else
  {
    bool completeClip = true;

    for(int i = 0; i < poly.numVertices(); i++)
    {
      if(plane.getOrientation(poly[i], eps) == ON_POSITIVE_SIDE)
      {
        completeClip = false;
        break;
      }
    }

    if(completeClip)
    {
      poly.clear();
    }
  }

  return;
}

/*!
 * \brief Clips a polyhedron against an array of planes
 *
 * \param [inout] poly The polyhedron to clip
 * \param [in] planes The array of planes
 * \param [in] eps The tolerance for plane point orientation
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE void clipPolyhedron(Polyhedron<T, NDIMS>& poly,
                                     axom::ArrayView<Plane<T, NDIMS>> planes,
                                     double eps)
{
  using PlaneType = Plane<T, NDIMS>;

  // Clip Polyhedron by each plane
  for(const PlaneType& plane : planes)
  {
    clipPolyhedron(poly, plane, eps);

    if(poly.numVertices() == 0)
    {
      return;
    }
  }

  return;
}

/*!
 * \brief Finds the clipped intersection Polyhedron between Tetrahedron
 *        tet and Polyhedron poly.
 *
 * \param [in] tet The tetrahedron to clip
 * \param [in] poly The polyhedron to clip against
 * \param [in] eps The tolerance for plane point orientation.
 * \param [in] tryFixOrientation Check if the signed volume of each shape is positive.
 * \return The Polyhedron formed from clipping the tetrahedron with a polyhedron.
 *
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE Polyhedron<T, NDIMS> clipTetrahedron(const Tetrahedron<T, NDIMS>& tet,
                                                      const Polyhedron<T, NDIMS>& poly,
                                                      double eps,
                                                      bool tryFixOrientation)
{
  using PlaneType = Plane<T, NDIMS>;
  using PolyhedronType = Polyhedron<T, NDIMS>;

  // Initialize our polyhedron to return
  PolyhedronType clipped = PolyhedronType::from_primitive(tet, tryFixOrientation);

  int numPlanes = 0;
  auto planes = poly.getFaces(numPlanes);
  axom::ArrayView<PlaneType> planesView(planes.data(), numPlanes);
  clipPolyhedron(clipped, planesView, eps);
  return clipped;
}

/*!
 * \brief Finds the clipped intersection Polyhedron between Hexahedron
 *        hex and Polyhedron poly.
 *
 * \param [in] hex The hexahedron to clip
 * \param [in] poly The polyhedron to clip against
 * \param [in] eps The tolerance for plane point orientation.
 * \param [in] tryFixOrientation Check if the signed volume of each shape is positive.
 * \return The Polyhedron formed from clipping the hexahedron with a polyhedron.
 *
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE Polyhedron<T, NDIMS> clipHexahedron(const Hexahedron<T, NDIMS>& hex,
                                                     const Polyhedron<T, NDIMS>& poly,
                                                     double eps,
                                                     bool tryFixOrientation)
{
  using PlaneType = Plane<T, NDIMS>;
  using PolyhedronType = Polyhedron<T, NDIMS>;

  // Initialize our polyhedron to return
  PolyhedronType clipped = PolyhedronType::from_primitive(hex, tryFixOrientation);

  int numPlanes = 0;
  auto planes = poly.getFaces(numPlanes);
  axom::ArrayView<PlaneType> planesView(planes.data(), numPlanes);
  clipPolyhedron(clipped, planesView, eps);
  return clipped;
}

/*!
 * \brief Clips a polyhedron against another polyhedron.
 *
 * \param [in] poly1 The polyhedron to clip
 * \param [in] poly2 The polyhedron to use to clip the first polyhedron.
 * \param [in] eps The tolerance for plane point orientation
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE Polyhedron<T, NDIMS> clipPolyhedronPolyhedron(const Polyhedron<T, NDIMS>& poly1,
                                                               const Polyhedron<T, NDIMS>& poly2,
                                                               double eps)
{
  using PlaneType = Plane<T, NDIMS>;

  // Get the plane faces in poly2.
  int numPlanes = 0;
  auto planes = poly2.getFaces(numPlanes);
  axom::ArrayView<PlaneType> planesView(planes.data(), numPlanes);

  // Clip Polyhedron by each plane in poly2.
  Polyhedron<T, NDIMS> clipped(poly1);
  clipPolyhedron(clipped, planesView, eps);
  return clipped;
}

/*!
 * \brief Finds the clipped intersection Polyhedron between Hexahedron
 *        hex and Tetrahedron tet.
 *
 * \param [in] hex The hexahedron
 * \param [in] tet The tetrahedron
 * \param [in] eps The tolerance for plane point orientation.
 * \param [in] tryFixOrientation Check if the signed volume of each shape is positive.
 * \return The Polyhedron formed from clipping the hexahedron with a tetrahedron.
 *
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE Polyhedron<T, NDIMS> clipHexahedron(const Hexahedron<T, NDIMS>& hex,
                                                     const Tetrahedron<T, NDIMS>& tet,
                                                     double eps,
                                                     bool tryFixOrientation)
{
  using PlaneType = Plane<T, NDIMS>;
  using PolyhedronType = Polyhedron<T, NDIMS>;

  // Initialize our polyhedron to return
  PolyhedronType poly = PolyhedronType::from_primitive(hex, tryFixOrientation);

  // Initialize planes from tetrahedron vertices
  // (Ordering here matters to get the correct winding)
  PlaneType planes[4] = {make_plane(tet[1], tet[3], tet[2]),
                         make_plane(tet[0], tet[2], tet[3]),
                         make_plane(tet[0], tet[3], tet[1]),
                         make_plane(tet[0], tet[1], tet[2])};

  // Adjusts planes in case tetrahedron signed volume is negative
  if(tryFixOrientation)
  {
    PolyhedronType tet_poly = PolyhedronType::from_primitive(tet, tryFixOrientation);
    planes[0] = make_plane(tet_poly[1], tet_poly[3], tet_poly[2]);
    planes[1] = make_plane(tet_poly[0], tet_poly[2], tet_poly[3]);
    planes[2] = make_plane(tet_poly[0], tet_poly[3], tet_poly[1]);
    planes[3] = make_plane(tet_poly[0], tet_poly[1], tet_poly[2]);
  }

  axom::StackArray<IndexType, 1> planeSize = {4};
  axom::ArrayView<PlaneType> planesView(planes, planeSize);

  clipPolyhedron(poly, planesView, eps);
  return poly;
}

/*!
 * \brief This object breaks a hex face into 2 triangles and indicates a point
 *        that can be tested for whether the face is planar.
 */
struct HexFaces
{
  int plane1Ids[3];
  int testId;
  int plane2Ids[3];
};

/*!
 * \brief Get a view containing planes for an input hex (or hex-like) shape.
 *
 * \param shape The input hex shape.
 * \param planes An array of planes large enough to hold the planes we'll make.
 * \param eps A tolerance that indicates how far from the plane a point can be
 *            in order to be planar.
 *
 * \return A view containing the planes that make up the hex shape.
 */
template <typename ShapeType, typename PlaneType>
AXOM_HOST_DEVICE axom::ArrayView<PlaneType> getHexahedronPlanes(const ShapeType& shape,
                                                                PlaneType* planes,
                                                                double eps)
{
  /*
      3-------2
     /|      /|
    / |     / |
   7-------6  |
   |  0----|--1
   | /     | /
   |/      |/
   4-------5

   */
  const HexFaces faces[] = {
    {{1, 2, 0}, 3, {0, 2, 3}},  // Back face
    {{4, 7, 5}, 6, {5, 7, 6}},  // Front face
    {{0, 4, 1}, 5, {4, 5, 1}},  // Bottom face
    {{5, 6, 1}, 2, {1, 6, 2}},  // Right face
    {{2, 6, 3}, 7, {3, 6, 7}},  // Top face
    {{0, 3, 4}, 7, {4, 3, 7}}   // Left face
  };
  int planeCount = 0;
  for(int f = 0; f < 6; f++)
  {
    const HexFaces& face = faces[f];
    // Make a plane with 3 points of the quad face.
    planes[planeCount] =
      make_plane(shape[face.plane1Ids[0]], shape[face.plane1Ids[1]], shape[face.plane1Ids[2]]);
    // Get the signed distance to another point.
    const auto dist = planes[planeCount].signedDistance(shape[face.testId]);
    planeCount++;

    // If the point was far enough off the plane then the face is not
    // planar and we should add another plane for the other triangle
    // in the hex face.
    if(axom::utilities::abs(dist) > eps)
    {
      // The plane was not planar enough. Add the other plane too.
      planes[planeCount] =
        make_plane(shape[face.plane2Ids[0]], shape[face.plane2Ids[1]], shape[face.plane2Ids[2]]);
      planeCount++;
    }
  }

  return axom::ArrayView<PlaneType>(planes, planeCount);
}

/*!
 * \brief Finds the clipped intersection Polyhedron between Hexahedron
 *        hex1 and Hexahedron hex2.
 *
 * \param [in] hex1 The hexahedron to clip.
 * \param [in] hex2 The hexahedron for clipping.
 * \param [in] eps The tolerance for plane point orientation.
 * \param [in] tryFixOrientation Check if the signed volume of each shape is positive.
 *
 * \note hex1 and hex2 are assumed to be convex.
 * \note If hex2's faces are not planar then each quad face will be divided into
 *       2 cutting planes.
 *
 * \return The Polyhedron formed from clipping the hexahedron with a hexahedron.
 *
 */

template <typename T, int NDIMS>
AXOM_HOST_DEVICE Polyhedron<T, NDIMS> clipHexahedron(const Hexahedron<T, NDIMS>& hex1,
                                                     const Hexahedron<T, NDIMS>& hex2,
                                                     double eps,
                                                     bool tryFixOrientation)
{
  using PlaneType = Plane<T, NDIMS>;
  using PolyhedronType = Polyhedron<T, NDIMS>;

  // Initialize our polyhedron to return
  PolyhedronType poly = PolyhedronType::from_primitive(hex1, tryFixOrientation);

  // Get clipping planes from hex2.
  constexpr int MAX_PLANES = 12;
  PlaneType planes[MAX_PLANES];
  axom::ArrayView<PlaneType> planesView;

  // Adjusts planes in case hexhedron signed volume is negative
  if(tryFixOrientation)
  {
    PolyhedronType hex2_poly = PolyhedronType::from_primitive(hex2, tryFixOrientation);
    // Get planes from the Polyhedral version of the plane, which could have
    // reordered the points.
    planesView = getHexahedronPlanes(hex2_poly, planes, eps);
  }
  else
  {
    planesView = getHexahedronPlanes(hex2, planes, eps);
  }

  clipPolyhedron(poly, planesView, eps);

  return poly;
}

/*!
 * \brief Finds the clipped intersection Polyhedron between Octahedron
 *        oct and Tetrahedron tet.
 *
 * \param [in] oct The octahedron
 * \param [in] tet The tetrahedron
 * \param [in] eps The tolerance for plane point orientation.
 * \param [in] tryFixOrientation Check if the signed volume of each shape is positive.
 * \return The Polyhedron formed from clipping the octahedron with a tetrahedron.
 *
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE Polyhedron<T, NDIMS> clipOctahedron(const Octahedron<T, NDIMS>& oct,
                                                     const Tetrahedron<T, NDIMS>& tet,
                                                     double eps,
                                                     bool tryFixOrientation)
{
  using PlaneType = Plane<T, NDIMS>;
  using PolyhedronType = Polyhedron<T, NDIMS>;

  // Initialize our polyhedron to return
  PolyhedronType poly = PolyhedronType::from_primitive(oct, tryFixOrientation);

  // Initialize planes from tetrahedron vertices
  // (Ordering here matters to get the correct winding)
  PlaneType planes[4] = {make_plane(tet[1], tet[3], tet[2]),
                         make_plane(tet[0], tet[2], tet[3]),
                         make_plane(tet[0], tet[3], tet[1]),
                         make_plane(tet[0], tet[1], tet[2])};

  // Adjusts planes in case tetrahedron signed volume is negative
  if(tryFixOrientation)
  {
    PolyhedronType tet_poly = PolyhedronType::from_primitive(tet, tryFixOrientation);
    planes[0] = make_plane(tet_poly[1], tet_poly[3], tet_poly[2]);
    planes[1] = make_plane(tet_poly[0], tet_poly[2], tet_poly[3]);
    planes[2] = make_plane(tet_poly[0], tet_poly[3], tet_poly[1]);
    planes[3] = make_plane(tet_poly[0], tet_poly[1], tet_poly[2]);
  }

  axom::StackArray<IndexType, 1> planeSize = {4};
  axom::ArrayView<PlaneType> planesView(planes, planeSize);

  clipPolyhedron(poly, planesView, eps);
  return poly;
}

/*!
 * \brief Finds the clipped intersection Polyhedron between Tetrahedron
 *        tet1 and Tetrahedron tet2.
 *
 * \param [in] tet1 The tetrahedron to clip
 * \param [in] tet2 The tetrahedron to clip against
 * \param [in] eps The tolerance for plane point orientation.
 * \param [in] tryFixOrientation Check if the signed volume of each shape is positive.
 * \return The Polyhedron formed from clipping the tetrahedron with a tetrahedron.
 *
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE Polyhedron<T, NDIMS> clipTetrahedron(const Tetrahedron<T, NDIMS>& tet1,
                                                      const Tetrahedron<T, NDIMS>& tet2,
                                                      double eps,
                                                      bool tryFixOrientation)
{
  using PlaneType = Plane<T, NDIMS>;
  using PolyhedronType = Polyhedron<T, NDIMS>;

  // Initialize our polyhedron to return
  PolyhedronType poly = PolyhedronType::from_primitive(tet1, tryFixOrientation);

  // Initialize planes from tetrahedron vertices
  // (Ordering here matters to get the correct winding)
  PlaneType planes[4] = {make_plane(tet2[1], tet2[3], tet2[2]),
                         make_plane(tet2[0], tet2[2], tet2[3]),
                         make_plane(tet2[0], tet2[3], tet2[1]),
                         make_plane(tet2[0], tet2[1], tet2[2])};

  // Adjusts planes in case tetrahedron signed volume is negative
  if(tryFixOrientation)
  {
    PolyhedronType tet_poly = PolyhedronType::from_primitive(tet2, tryFixOrientation);
    planes[0] = make_plane(tet_poly[1], tet_poly[3], tet_poly[2]);
    planes[1] = make_plane(tet_poly[0], tet_poly[2], tet_poly[3]);
    planes[2] = make_plane(tet_poly[0], tet_poly[3], tet_poly[1]);
    planes[3] = make_plane(tet_poly[0], tet_poly[1], tet_poly[2]);
  }

  axom::StackArray<IndexType, 1> planeSize = {4};
  axom::ArrayView<PlaneType> planesView(planes, planeSize);

  clipPolyhedron(poly, planesView, eps);
  return poly;
}

/*!
 * \brief Clips a tetrahedron against the half-space defined by a plane
 *        and returns the resulting polyhedron
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
 *         the half-space defined a plane
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
template <typename T, int NDIMS>
AXOM_HOST_DEVICE Polyhedron<T, NDIMS> clipTetrahedron(const Tetrahedron<T, NDIMS>& tet,
                                                      const Plane<T, NDIMS>& plane,
                                                      double eps,
                                                      bool tryFixOrientation)
{
  using PolyhedronType = Polyhedron<T, NDIMS>;
  PolyhedronType poly = PolyhedronType::from_primitive(tet, tryFixOrientation);
  clipPolyhedron(poly, plane, eps);
  return poly;
}

/*!
 * \brief Clips a 2D subject polygon against a clipping plane in 2D, returning
 *        their geometric intersection as a polygon.
 *
 * \param inputList The input polygon to be clipped.
 * \param plane The plane being used for clipping.
 * \param eps The tolerance used for intersection.
 *
 * \return A clipped polygon.
 */
AXOM_SUPPRESS_HD_WARN
template <typename PolygonType, typename PlaneType>
AXOM_HOST_DEVICE PolygonType clipPolygonPlaneSimple(const PolygonType& inputList,
                                                    const PlaneType& plane,
                                                    double eps)
{
  using T = typename PolygonType::PointType::CoordType;
  using SegmentType = Segment<T, 2>;
  using PointType = Point<T, 2>;

  PolygonType outputList;

  const int numVertices = inputList.numVertices();
  for(int iVert = 0; iVert < numVertices; iVert++)
  {
    const int prevVert = ((iVert - 1) == -1) ? (numVertices - 1) : (iVert - 1);
    const PointType& current_point = inputList[iVert];
    const PointType& prev_point = inputList[prevVert];

    T seg_param;
    PointType intersecting_point;
    SegmentType subject_edge(prev_point, current_point);

    const bool intersected = intersect(plane, subject_edge, seg_param, eps);
    if(intersected)
    {
      intersecting_point = subject_edge.at(seg_param);
    }

    const int cur_p_orientation = plane.getOrientation(current_point, eps);
    const int prev_p_orientation = plane.getOrientation(prev_point, eps);

    if(cur_p_orientation == ON_POSITIVE_SIDE)
    {
      if(prev_p_orientation != ON_POSITIVE_SIDE)
      {
        outputList.addVertex(intersecting_point);
      }
      outputList.addVertex(current_point);
    }
    else if(prev_p_orientation == ON_POSITIVE_SIDE)
    {
      // In this branch, there has been a plane crossing from positive to negative.
      if(intersected)
      {
        outputList.addVertex(intersecting_point);
      }
      else
      {
        // If there was not an intersection point then it was close but not
        // within eps tolerance. It would have been really close to current_point.
        // Since there was a plane crossing, add current_point.
        outputList.addVertex(current_point);
      }
    }
  }

  return outputList;
}

/*!
 * \brief Remove duplicate adjacent points in a polygon.
 *
 * \param poly The polygon to filter.
 * \param eps The tolerance for filtering.
 *
 * \note Handles edge cases such as 3 consecutive vertices with orientations
 *       ON_POSITIVE_SIDE, ON_BOUNDARY, ON_POSITIVE where point on the boundary
 *       is added twice.
 *
 * \return A new polygon that has unique points.
 */
AXOM_SUPPRESS_HD_WARN
template <typename PolygonType>
AXOM_HOST_DEVICE PolygonType makeUniquePoints(const PolygonType& poly, double eps)
{
  PolygonType uniqueList;
  const int numVertices = poly.numVertices();
  for(int i = 0; i < numVertices; i++)
  {
    int prevIndex = ((i - 1) == -1) ? (numVertices - 1) : (i - 1);
    const auto& curPoint = poly[i];
    const auto& prevPoint = poly[prevIndex];
    // Check whether curPoint and prevPoint are far enough apart to be different.
    // If so, add the point. NOTE: casts are needed to match eps when the polygon
    // does not use double precision coordinates.
    if(!axom::utilities::isNearlyEqual(static_cast<double>(curPoint[0]),
                                       static_cast<double>(prevPoint[0]),
                                       eps) ||
       !axom::utilities::isNearlyEqual(static_cast<double>(curPoint[1]),
                                       static_cast<double>(prevPoint[1]),
                                       eps))
    {
      uniqueList.addVertex(poly[i]);
    }
  }

  return uniqueList;
}

/*!
 * \brief Clips a 2D subject polygon against a clip polygon in 2D, returning
 *        their geometric intersection as a polygon.
 *
 * \sa axom::primal::clip()
 */
AXOM_SUPPRESS_HD_WARN
template <typename T, axom::primal::PolygonArray ARRAY_TYPE, int MAX_VERTS>
AXOM_HOST_DEVICE Polygon<T, 2, ARRAY_TYPE, MAX_VERTS> clipPolygonPolygon(
  const Polygon<T, 2, ARRAY_TYPE, MAX_VERTS>& subjectPolygon,
  const Polygon<T, 2, ARRAY_TYPE, MAX_VERTS>& clipPolygon,
  double eps = 1.e-10,
  bool tryFixOrientation = false)
{
  SLIC_ASSERT(ARRAY_TYPE == axom::primal::PolygonArray::Dynamic ||
              (ARRAY_TYPE == axom::primal::PolygonArray::Static &&
               MAX_VERTS >= (subjectPolygon.numVertices() + clipPolygon.numVertices())));

  using PlaneType = Plane<T, 2>;
  using PolygonType = Polygon<T, 2, ARRAY_TYPE, MAX_VERTS>;

  PolygonType outputList = subjectPolygon;
  PolygonType planePoints = clipPolygon;

  if(tryFixOrientation)
  {
    if(outputList.signedArea() < 0)
    {
      outputList.reverseOrientation();
    }

    if(planePoints.signedArea() < 0)
    {
      planePoints.reverseOrientation();
    }
  }

  const int numClipEdges = planePoints.numVertices();

  // Iterate through edges of clip polygon, represented as planes
  for(int iEdge = 0; iEdge < numClipEdges; iEdge++)
  {
    PlaneType plane = make_plane(planePoints[iEdge], planePoints[(iEdge + 1) % numClipEdges]);

    outputList = clipPolygonPlaneSimple(outputList, plane, eps);
  }  // end of iteration through edges of clip polygon

  // Remove duplicate points.
  return makeUniquePoints(outputList, eps);
}

/*!
 * \brief Clips a 2D subject polygon against a 2D plane (line), returning
 *        their geometric intersection as a polygon.
 *
 * \sa axom::primal::clip()
 */
template <typename T, axom::primal::PolygonArray ARRAY_TYPE, int MAX_VERTS>
AXOM_HOST_DEVICE Polygon<T, 2, ARRAY_TYPE, MAX_VERTS> clipPolygonPlane(
  const Polygon<T, 2, ARRAY_TYPE, MAX_VERTS>& subjectPolygon,
  const Plane<T, 2>& clipPlane,
  double eps = 1.e-10,
  bool tryFixOrientation = false)
{
  using PolygonType = Polygon<T, 2, ARRAY_TYPE, MAX_VERTS>;

  PolygonType outputList = subjectPolygon;
  if(tryFixOrientation)
  {
    if(outputList.signedArea() < 0)
    {
      outputList.reverseOrientation();
    }
  }

  // Clip the plane.
  return makeUniquePoints(clipPolygonPlaneSimple(outputList, clipPlane, eps), eps);
}

}  // namespace detail
}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_CLIP_IMPL_HPP_
