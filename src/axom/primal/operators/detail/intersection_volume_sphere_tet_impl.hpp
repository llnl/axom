// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/*!
 * \file intersection_volume_sphere_tet_impl.hpp
 *
 * \brief Implements exact sphere-tetrahedron overlap volume.
 *
 * The algorithm follows S. Strobl, A. Formella, and T. Poeschel,
 * "Exact calculation of the overlap volume of spheres and mesh elements,"
 * Journal of Computational Physics 311 (2016), 158-172.
 * https://doi.org/10.1016/j.jcp.2016.02.003
 *
 * The implementation adapts robustness details from the authors' MIT-licensed
 * reference at https://github.com/severinstrobl/overlap, revision
 * 7c5db4ed03653b560eafe0914397205c96e35cc2.
 *
 * MIT License
 *
 * Copyright (c) 2015-2025 Severin Strobl
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "axom/core/NumericLimits.hpp"
#include "axom/core/StackArray.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Sphere.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/operators/compute_bounding_box.hpp"
#include "axom/slic/interface/slic.hpp"

#include <cmath>

namespace axom
{
namespace primal
{
namespace detail
{
namespace sphere_tet
{

using Vector3 = primal::Vector<double, 3>;
using Point3 = primal::Point<double, 3>;
using Tetrahedron3 = primal::Tetrahedron<double, 3>;
using Triangle3 = primal::Triangle<double, 3>;

constexpr int NUM_VERTICES = Tetrahedron3::numVertices();
constexpr int NUM_EDGES = 6;
constexpr int NUM_FACES = 4;
constexpr int NUM_EDGE_ENDPOINTS = 2;
constexpr int NUM_EDGE_FACES = 2;
constexpr int NUM_FACE_VERTICES = 3;
constexpr int NUM_VERTEX_EDGES = 3;
constexpr int NUM_VERTEX_FACES = 3;
constexpr double TINY_EPSILON = 4.0 * axom::numeric_limits<double>::epsilon();
constexpr double LARGE_EPSILON = 1.e-10;

/*!
 * \brief Returns an endpoint vertex for a tetrahedron edge.
 *
 * The topology helpers in this section use the following numbering. Face
 * vertices are ordered so their normals point outward for a positively
 * oriented tetrahedron.
 *
 * \verbatim
 *
 *                    3
 *                    *
 *                e3 /|\ e5
 *                  / | \
 *                 / e4  \
 *                /   |   \
 *             0 *---e2----* 2
 *                \   |   /
 *             e0  \  |  / e1
 *                  \ | /
 *                   \|/
 *                    *
 *                    1
 *
 *   edge 0: (0,1)       face 0: (2,1,0)
 *   edge 1: (1,2)       face 1: (0,1,3)
 *   edge 2: (2,0)       face 2: (1,2,3)
 *   edge 3: (0,3)       face 3: (2,0,3)
 *   edge 4: (1,3)
 *   edge 5: (2,3)
 *
 * A vertex correction numbers its three edge-intersection points locally.
 * Its three local faces connect point pairs (0,1), (0,2), and (1,2):
 *
 *                    2
 *                    *
 *                   / \
 *          face 1  /   \  face 2
 *                 /     \
 *              0 *-------* 1
 *                   face 0
 *
 * \endverbatim
 *
 * \param [in] edge the edge index
 * \param [in] endpoint the local endpoint index
 * \return the tetrahedron vertex index
 * \pre `0 <= edge < NUM_EDGES`
 * \pre `0 <= endpoint < NUM_EDGE_ENDPOINTS`
 */
AXOM_HOST_DEVICE inline int edgeVertex(int edge, int endpoint)
{
  SLIC_ASSERT(edge >= 0 && edge < NUM_EDGES);
  SLIC_ASSERT(endpoint >= 0 && endpoint < NUM_EDGE_ENDPOINTS);

  if(edge < 3)
  {
    return endpoint == 0 ? edge : (edge + 1) % 3;
  }
  return endpoint == 0 ? edge - 3 : 3;
}

/*!
 * \brief Returns an incident face for a tetrahedron edge.
 * \param [in] edge the edge index
 * \param [in] side the local face index
 * \return the tetrahedron face index
 * \pre `0 <= edge < NUM_EDGES`
 * \pre `0 <= side < NUM_EDGE_FACES`
 */
AXOM_HOST_DEVICE inline int edgeFace(int edge, int side)
{
  SLIC_ASSERT(edge >= 0 && edge < NUM_EDGES);
  SLIC_ASSERT(side >= 0 && side < NUM_EDGE_FACES);

  if(side == 0)
  {
    return edge < 3 ? 0 : (edge < 5 ? 1 : 2);
  }

  switch(edge)
  {
  case 0:
    return 1;
  case 1:
  case 4:
    return 2;
  default:
    return 3;
  }
}

/*!
 * \brief Returns an incident edge for a tetrahedron vertex.
 * \param [in] vertex the vertex index
 * \param [in] local_edge the local incident-edge index
 * \return the tetrahedron edge index
 * \pre `0 <= vertex < NUM_VERTICES`
 * \pre `0 <= local_edge < NUM_VERTEX_EDGES`
 */
AXOM_HOST_DEVICE inline int vertexEdge(int vertex, int local_edge)
{
  SLIC_ASSERT(vertex >= 0 && vertex < NUM_VERTICES);
  SLIC_ASSERT(local_edge >= 0 && local_edge < NUM_VERTEX_EDGES);

  switch(vertex)
  {
  case 0:
    return local_edge == 0 ? 0 : (local_edge == 1 ? 2 : 3);
  case 1:
    return local_edge == 0 ? 0 : (local_edge == 1 ? 1 : 4);
  case 2:
    return local_edge == 0 ? 1 : (local_edge == 1 ? 2 : 5);
  default:
    return local_edge + 3;
  }
}

/*!
 * \brief Returns the endpoint index of a vertex on an incident edge.
 * \param [in] vertex the vertex index
 * \param [in] local_edge the local incident-edge index
 * \return the local endpoint index
 * \pre `0 <= vertex < NUM_VERTICES`
 * \pre `0 <= local_edge < NUM_VERTEX_EDGES`
 */
AXOM_HOST_DEVICE inline int vertexEdgeEnd(int vertex, int local_edge)
{
  SLIC_ASSERT(vertex >= 0 && vertex < NUM_VERTICES);
  SLIC_ASSERT(local_edge >= 0 && local_edge < NUM_VERTEX_EDGES);

  if(vertex == 0)
  {
    return local_edge == 1 ? 1 : 0;
  }
  return vertex == 3 || local_edge == 0 ? 1 : 0;
}

/*!
 * \brief Returns an incident face for a tetrahedron vertex.
 * \param [in] vertex the vertex index
 * \param [in] local_face the local incident-face index
 * \return the tetrahedron face index
 * \pre `0 <= vertex < NUM_VERTICES`
 * \pre `0 <= local_face < NUM_VERTEX_FACES`
 */
AXOM_HOST_DEVICE inline int vertexFace(int vertex, int local_face)
{
  SLIC_ASSERT(vertex >= 0 && vertex < NUM_VERTICES);
  SLIC_ASSERT(local_face >= 0 && local_face < NUM_VERTEX_FACES);

  switch(vertex)
  {
  case 0:
    return local_face < 2 ? local_face : 3;
  case 1:
    return local_face;
  case 2:
    return local_face == 0 ? 0 : local_face + 1;
  default:
    return local_face == 0 ? 1 : (local_face == 1 ? 3 : 2);
  }
}

/*!
 * \brief Returns an edge index in the local vertex-correction triangle.
 * \param [in] local_face the local face index
 * \param [in] side the local edge endpoint index
 * \return the vertex-correction point index
 * \pre `0 <= local_face < NUM_VERTEX_FACES`
 * \pre `0 <= side < NUM_EDGE_ENDPOINTS`
 */
AXOM_HOST_DEVICE inline int vertexFaceEdge(int local_face, int side)
{
  SLIC_ASSERT(local_face >= 0 && local_face < NUM_VERTEX_FACES);
  SLIC_ASSERT(side >= 0 && side < NUM_EDGE_ENDPOINTS);

  return side == 0 ? (local_face == 2 ? 1 : 0) : (local_face == 0 ? 1 : 2);
}

/*!
 * \brief Returns a vertex for an oriented tetrahedron face.
 * \param [in] face the face index
 * \param [in] local_vertex the local face-vertex index
 * \return the tetrahedron vertex index
 * \pre `0 <= face < NUM_FACES`
 * \pre `0 <= local_vertex < NUM_FACE_VERTICES`
 */
AXOM_HOST_DEVICE inline int faceVertex(int face, int local_vertex)
{
  SLIC_ASSERT(face >= 0 && face < NUM_FACES);
  SLIC_ASSERT(local_vertex >= 0 && local_vertex < NUM_FACE_VERTICES);

  switch(face)
  {
  case 0:
    return local_vertex == 0 ? 2 : (local_vertex == 1 ? 1 : 0);
  case 1:
    return local_vertex == 2 ? 3 : local_vertex;
  case 2:
    return local_vertex == 0 ? 1 : (local_vertex == 1 ? 2 : 3);
  default:
    return local_vertex == 0 ? 2 : (local_vertex == 1 ? 0 : 3);
  }
}

/*!
 * \brief Normalizes a vector while preserving the zero vector.
 * \param [in] vector the vector
 * \return the normalized vector, or zero when \p vector is zero
 */
AXOM_HOST_DEVICE inline Vector3 normalizedOrZero(const Vector3& vector)
{
  const double norm = vector.norm();
  return norm > 0.0 ? vector / norm : Vector3 {};
}

/*!
 * \brief Computes the angle between normalized vectors.
 * \param [in] lhs the first normalized vector
 * \param [in] rhs the second normalized vector
 * \return the angle in radians
 * \note The clamped dot product keeps roundoff from passing a value outside
 * the domain of `acos()`.
 */
AXOM_HOST_DEVICE inline double angle(const Vector3& lhs, const Vector3& rhs)
{
  return ::acos(axom::utilities::clampVal(lhs.dot(rhs), -1.0, 1.0));
}

/// \brief Stores the representative point and unit normal for a tetrahedron face.
struct Face
{
  Vector3 center;
  Vector3 normal;
};

/// \brief Caches normalized tetrahedron geometry and its derived quantities.
struct TetCache
{
  Tetrahedron3 geometry;
  axom::StackArray<Face, NUM_FACES> faces;
  Vector3 center;
  double volume;
};

/*!
 * \brief Returns a tetrahedron vertex as a vector from the unit-sphere center.
 * \param [in] tet the cached tetrahedron
 * \param [in] vertex the vertex index
 * \return the vertex coordinates as a vector
 */
AXOM_HOST_DEVICE inline Vector3 vertexVector(const TetCache& tet, int vertex)
{
  return Vector3(tet.geometry[vertex]);
}

/*!
 * \brief Constructs an oriented triangular face.
 * \param [in] a the first face vertex
 * \param [in] b the second face vertex
 * \param [in] c the third face vertex
 * \return the face with its center and unit normal initialized
 */
AXOM_HOST_DEVICE inline Face makeFace(const Vector3& a, const Vector3& b, const Vector3& c)
{
  const Triangle3 triangle(Point3(a.array()), Point3(b.array()), Point3(c.array()));
  Face face;
  face.center = Vector3(triangle.centroid());
  face.normal = normalizedOrZero(triangle.normal());
  return face;
}

/*!
 * \brief Constructs a positively oriented tetrahedron cache and its outward faces.
 * \param [in] geometry the tetrahedron geometry
 * \return the initialized tetrahedron cache
 */
AXOM_HOST_DEVICE inline TetCache makeTetCache(Tetrahedron3 geometry)
{
  // Positive orientation makes the face order below point outward.
  geometry.checkAndFixOrientation();

  TetCache tet;
  tet.geometry = geometry;
  const Vector3 v0 = vertexVector(tet, 0);
  const Vector3 v1 = vertexVector(tet, 1);
  const Vector3 v2 = vertexVector(tet, 2);
  const Vector3 v3 = vertexVector(tet, 3);
  tet.faces[0] = makeFace(v2, v1, v0);
  tet.faces[1] = makeFace(v0, v1, v3);
  tet.faces[2] = makeFace(v1, v2, v3);
  tet.faces[3] = makeFace(v2, v0, v3);
  tet.center = (v0 + v1 + v2 + v3) / 4.0;
  tet.volume = geometry.volume();
  return tet;
}

/*!
 * \brief Determines whether a point projects inside a triangular face.
 * \param [in] tet the tetrahedron
 * \param [in] face_index the face index
 * \param [in] point the point
 * \return true when the projected point is inside the face
 */
AXOM_HOST_DEVICE inline bool faceContainsProjection(const TetCache& tet,
                                                    int face_index,
                                                    const Vector3& point)
{
  const Face& face = tet.faces[face_index];
  const Vector3 projected = point - face.normal.dot(point - face.center) * face.normal;
  const Triangle3 triangle(tet.geometry[faceVertex(face_index, 0)],
                           tet.geometry[faceVertex(face_index, 1)],
                           tet.geometry[faceVertex(face_index, 2)]);
  return triangle.contains(Point3(projected.array()));
}

/*!
 * \brief Computes the volume of a unit-sphere cap.
 * \param [in] height the cap height
 * \return the cap volume
 */
AXOM_HOST_DEVICE inline double capVolume(double height)
{
  if(height <= 0.0)
  {
    return 0.0;
  }
  if(height >= 2.0)
  {
    return 4.0 * M_PI / 3.0;
  }
  return (M_PI / 3.0) * height * height * (3.0 - height);
}

/// \brief Stores the ordered parameters of a proper line-sphere intersection.
struct EdgeRoots
{
  bool intersects {false};
  double first {0.0};
  double second {0.0};
};

/*!
 * \brief Intersects a line with the unit sphere.
 * \param [in] base the line origin
 * \param [in] direction the line direction
 * \return the ordered line parameters when two proper intersections exist
 */
AXOM_HOST_DEVICE inline EdgeRoots intersectUnitSphere(const Vector3& base, const Vector3& direction)
{
  const double a = direction.squared_norm();
  if(a == 0.0)
  {
    return {};
  }

  const double b = 2.0 * direction.dot(base);
  const double c = base.squared_norm() - 1.0;
  const double discriminant = b * b - 4.0 * a * c;
  if(discriminant <= 0.0)
  {
    return {};
  }

  // This form of the quadratic formula avoids cancellation in one root;
  // the product of the roots supplies the other root.
  const double q = -0.5 * (b + ::copysign(::sqrt(discriminant), b));
  double first = q / a;
  double second = c / q;
  if(first > second)
  {
    axom::utilities::swap(first, second);
  }

  // The midpoint between roots is the closest point on the line to the
  // sphere center. It remains strictly inside for a proper crossing.
  const Vector3 midpoint = base + (0.5 * (first + second)) * direction;
  if(midpoint.squared_norm() >= 1.0)
  {
    return {};
  }
  return {true, first, second};
}

/*!
 * \brief Computes the principal unit-sphere wedge term.
 * \param [in] distance distance from the sphere center to the wedge edge
 * \param [in] alpha the wedge angle
 * \return the principal wedge volume
 */
AXOM_HOST_DEVICE inline double principalWedge(double distance, double alpha)
{
  const double sin_alpha = ::sin(alpha);
  const double cos_alpha = ::cos(alpha);
  const double a = distance * sin_alpha;
  const double b = ::sqrt(::fabs(1.0 - distance * distance));
  const double c = distance * cos_alpha;
  return (a * b * c) / 3.0 + a * (a * a / 3.0 - 1.0) * ::atan2(b, c) +
    (2.0 / 3.0) * ::atan2(sin_alpha * b, cos_alpha);
}

/*!
 * \brief Computes a regularized unit-sphere wedge volume.
 * \param [in] distance distance from the sphere center to the wedge edge
 * \param [in] alpha the wedge angle
 * \param [in] axial_distance signed distance along the wedge axis
 * \return the wedge volume
 */
AXOM_HOST_DEVICE inline double regularizedWedge(double distance, double alpha, double axial_distance)
{
  const bool obtuse = alpha > 0.5 * M_PI;
  const double acute_alpha = obtuse ? M_PI - alpha : alpha;
  const double principal = principalWedge(distance, acute_alpha);

  if(axial_distance >= 0.0)
  {
    if(obtuse)
    {
      const double height = 1.0 - axial_distance;
      return capVolume(height) - principal;
    }
    return principal;
  }

  const double hemisphere_volume = 2.0 * M_PI / 3.0;
  if(obtuse)
  {
    return hemisphere_volume - principal;
  }
  return hemisphere_volume - (capVolume(1.0 + axial_distance) - principal);
}

/// \brief Stores the representative point and unit normal for a wedge plane.
struct Plane
{
  Vector3 center;
  Vector3 normal;
};

/*!
 * \brief Computes the unit-sphere wedge bounded by two planes.
 * \param [in] first the first bounding plane
 * \param [in] second the second bounding plane
 * \param [in] offset the vector from the sphere center to the planes' intersection line
 * \return the wedge volume
 */
AXOM_HOST_DEVICE inline double generalWedge(const Plane& first,
                                            const Plane& second,
                                            const Vector3& offset)
{
  const double distance = offset.norm();
  // Planes whose intersection line passes through the sphere center bound a
  // sector whose volume depends only on their angle.
  if(distance < TINY_EPSILON)
  {
    return (2.0 / 3.0) * (M_PI - angle(first.normal, second.normal));
  }
  if(distance >= 1.0)
  {
    return 0.0;
  }

  const double first_distance = offset.dot(first.normal);
  const double second_distance = offset.dot(second.normal);
  if(::fabs(first_distance) < TINY_EPSILON || ::fabs(second_distance) < TINY_EPSILON)
  {
    // A plane through the sphere center reduces the general case to the
    // regularized one-plane wedge formula.
    const double alpha = M_PI - angle(first.normal, second.normal);
    const double axial =
      ::fabs(first_distance) > ::fabs(second_distance) ? first_distance : second_distance;
    return regularizedWedge(distance, alpha, axial);
  }

  Vector3 direction = offset / distance;
  if(distance < LARGE_EPSILON)
  {
    // Remove roundoff parallel to the planes' intersection line when the
    // radial direction is poorly conditioned near the sphere center.
    const Vector3 normal_cross =
      normalizedOrZero(Vector3::cross_product(first.normal, second.normal));
    direction = normalizedOrZero(direction - direction.dot(normal_cross) * normal_cross);
  }

  double first_angle = angle(first.normal, direction);
  double second_angle = angle(second.normal, direction);
  const double first_direction = direction.dot(offset - first.center);
  const double second_direction = direction.dot(offset - second.center);

  if(first_distance >= 0.0 && second_distance >= 0.0)
  {
    // The sphere center lies below both oriented planes, so the desired
    // region is the sum of two regularized wedges.
    first_angle = 0.5 * M_PI - ::copysign(first_angle, first_direction);
    second_angle = 0.5 * M_PI - ::copysign(second_angle, second_direction);
    return regularizedWedge(distance, first_angle, first_distance) +
      regularizedWedge(distance, second_angle, second_distance);
  }

  if(first_distance < 0.0 && second_distance < 0.0)
  {
    // The sphere center lies above both planes, so compute the complement of
    // the two opposite wedges.
    first_angle = 0.5 * M_PI + ::copysign(1.0, first_direction) * (first_angle - M_PI);
    second_angle = 0.5 * M_PI + ::copysign(1.0, second_direction) * (second_angle - M_PI);
    return 4.0 * M_PI / 3.0 -
      (regularizedWedge(distance, first_angle, -first_distance) +
       regularizedWedge(distance, second_angle, -second_distance));
  }

  // The center lies between the planes. The desired region is the difference
  // between their regularized wedges.
  first_angle = 0.5 * M_PI -
    ::copysign(1.0, first_direction * first_distance) *
      (first_angle - (first_distance < 0.0 ? M_PI : 0.0));
  second_angle = 0.5 * M_PI -
    ::copysign(1.0, second_direction * second_distance) *
      (second_angle - (second_distance < 0.0 ? M_PI : 0.0));

  const double first_volume = regularizedWedge(distance, first_angle, ::fabs(first_distance));
  const double second_volume = regularizedWedge(distance, second_angle, ::fabs(second_distance));
  return ::fabs(first_volume - second_volume);
}

/// \brief Stores intersected entities and their edge-intersection offsets.
struct IntersectionData
{
  unsigned int vertices {0};
  unsigned int edges {0};
  unsigned int faces {0};
  axom::StackArray<axom::StackArray<Vector3, NUM_EDGE_ENDPOINTS>, NUM_EDGES> edge_points;
};

/*!
 * \brief Determines whether an entity is present in a compact bit mask.
 * \param [in] entities the entity bit mask
 * \param [in] index the entity index
 * \return true when the entity is marked
 */
AXOM_HOST_DEVICE inline bool isMarked(unsigned int entities, int index)
{
  return (entities & (1U << index)) != 0;
}

/*!
 * \brief Determines whether a triangular face intersects the unit sphere.
 * \param [in] tet the tetrahedron
 * \param [in] face the face index
 * \return true when the sphere crosses the face interior
 */
AXOM_HOST_DEVICE inline bool faceIntersectsUnitSphere(const TetCache& tet, int face)
{
  const double distance = tet.faces[face].normal.dot(-tet.faces[face].center);
  return distance * distance < 1.0 && faceContainsProjection(tet, face, Vector3 {});
}

/*!
 * \brief Finds sphere intersections with tetrahedron vertices, edges, and faces.
 * \param [in] tet the tetrahedron in unit-sphere coordinates
 * \return the marked entities and edge intersection points
 */
AXOM_HOST_DEVICE inline IntersectionData findIntersections(const TetCache& tet)
{
  IntersectionData data;

  // Intersect every topological edge and mark its incident entities.
  for(int edge = 0; edge < NUM_EDGES; ++edge)
  {
    const int first_vertex = edgeVertex(edge, 0);
    const int second_vertex = edgeVertex(edge, 1);
    const Vector3 base = vertexVector(tet, first_vertex);
    const Vector3 direction = vertexVector(tet, second_vertex) - base;
    const EdgeRoots roots = intersectUnitSphere(base, direction);
    if(!roots.intersects || roots.first >= 1.0 || roots.second <= 0.0)
    {
      continue;
    }

    data.vertices |= (roots.first < 0.0 ? 1U : 0U) << first_vertex;
    data.vertices |= (roots.second > 1.0 ? 1U : 0U) << second_vertex;
    data.edge_points[edge][0] = roots.first * direction;
    data.edge_points[edge][1] = (roots.second - 1.0) * direction;
    data.edges |= 1U << edge;
    data.faces |= 1U << edgeFace(edge, 0);
    data.faces |= 1U << edgeFace(edge, 1);
  }

  // A vertex correction is valid only when the sphere crosses all three
  // edges incident to that vertex.
  for(int vertex = 0; vertex < NUM_VERTICES; ++vertex)
  {
    if(isMarked(data.vertices, vertex))
    {
      const bool all_edges_intersect = isMarked(data.edges, vertexEdge(vertex, 0)) &&
        isMarked(data.edges, vertexEdge(vertex, 1)) && isMarked(data.edges, vertexEdge(vertex, 2));
      if(!all_edges_intersect)
      {
        data.vertices &= ~(1U << vertex);
      }
    }
  }

  // An edge test does not detect a sphere that crosses only a face interior.
  for(int face = 0; face < NUM_FACES; ++face)
  {
    data.faces |= (faceIntersectsUnitSphere(tet, face) ? 1U : 0U) << face;
  }
  return data;
}

/*!
 * \brief Computes the spherical wedge contribution for an intersected edge.
 * \param [in] tet the tetrahedron
 * \param [in] data the cached sphere intersections
 * \param [in] edge the edge index
 * \return the edge wedge volume
 */
AXOM_HOST_DEVICE inline double edgeWedge(const TetCache& tet, const IntersectionData& data, int edge)
{
  const Face& first_face = tet.faces[edgeFace(edge, 0)];
  const Face& second_face = tet.faces[edgeFace(edge, 1)];
  const Vector3 midpoint = 0.5 *
    ((data.edge_points[edge][0] + vertexVector(tet, edgeVertex(edge, 0))) +
     (data.edge_points[edge][1] + vertexVector(tet, edgeVertex(edge, 1))));
  return generalWedge({first_face.center, first_face.normal},
                      {second_face.center, second_face.normal},
                      midpoint);
}

/*!
 * \brief Computes the spherical-cone correction at an intersected vertex.
 * \param [in] tet the tetrahedron
 * \param [in] data the cached sphere intersections
 * \param [in] vertex the vertex index
 * \return the vertex correction volume
 */
AXOM_HOST_DEVICE inline double vertexCorrection(const TetCache& tet,
                                                const IntersectionData& data,
                                                int vertex)
{
  struct IndexedDistance
  {
    int index;
    double squared_distance;
  };

  // Gather the three edge crossings relative to the selected vertex.
  axom::StackArray<Vector3, NUM_VERTEX_EDGES> relative_points;
  axom::StackArray<Vector3, NUM_VERTEX_EDGES> points;
  axom::StackArray<IndexedDistance, NUM_VERTEX_EDGES> distances;
  const Vector3 vertex_position = vertexVector(tet, vertex);
  for(int local_edge = 0; local_edge < NUM_VERTEX_EDGES; ++local_edge)
  {
    const int edge = vertexEdge(vertex, local_edge);
    relative_points[local_edge] = data.edge_points[edge][vertexEdgeEnd(vertex, local_edge)];
    points[local_edge] = relative_points[local_edge] + vertex_position;
    distances[local_edge] = {local_edge, relative_points[local_edge].squared_norm()};
  }
  // Sort the three crossings from nearest to farthest with a fixed-size
  // sorting network.
  if(distances[1].squared_distance < distances[0].squared_distance)
  {
    axom::utilities::swap(distances[0], distances[1]);
  }
  if(distances[2].squared_distance < distances[1].squared_distance)
  {
    axom::utilities::swap(distances[1], distances[2]);
  }
  if(distances[1].squared_distance < distances[0].squared_distance)
  {
    axom::utilities::swap(distances[0], distances[1]);
  }

  if(distances[1].squared_distance < distances[2].squared_distance * LARGE_EPSILON)
  {
    // A collapsed cone base reduces the corner correction to one edge wedge.
    return edgeWedge(tet, data, vertexEdge(vertex, distances[2].index));
  }

  const Vector3 center = (points[0] + points[1] + points[2]) / 3.0;
  const Triangle3 correction_triangle(Point3(relative_points[0].array()),
                                      Point3(relative_points[1].array()),
                                      Point3(relative_points[2].array()));
  Vector3 normal = normalizedOrZero(correction_triangle.normal());
  // The cone base normal points away from the tetrahedron interior.
  if(normal.dot(tet.center - center) > 0.0)
  {
    normal = -normal;
  }

  const Plane cone_plane {center, normal};
  double segment_correction = 0.0;
  // Each incident tetrahedron face and the cone base bound a spherical wedge
  // that is removed from the cone's cap.
  for(int local_face = 0; local_face < NUM_VERTEX_FACES; ++local_face)
  {
    const Face& face = tet.faces[vertexFace(vertex, local_face)];
    const Vector3 midpoint =
      0.5 * (points[vertexFaceEdge(local_face, 0)] + points[vertexFaceEdge(local_face, 1)]);
    segment_correction += generalWedge(cone_plane, {face.center, -face.normal}, midpoint);
  }

  // The correction combines the straight tetrahedral tip with the spherical
  // cap above it, then removes the three side wedges.
  const double tip_volume =
    ::fabs((-relative_points[2])
             .dot(Vector3::cross_product(relative_points[0] - relative_points[2],
                                         relative_points[1] - relative_points[2]))) /
    6.0;
  const double cap = capVolume(1.0 + normal.dot(-center));
  if(cap < TINY_EPSILON)
  {
    return tip_volume;
  }
  return axom::utilities::max(0.0, tip_volume + cap - segment_correction);
}

/*!
 * \brief Computes unit-sphere overlap with a normalized tetrahedron.
 * \param [in] tet the tetrahedron in unit-sphere coordinates
 * \return the overlap volume
 */
AXOM_HOST_DEVICE inline double unitSphereTetOverlap(const TetCache& tet)
{
  if(tet.volume == 0.0)
  {
    return 0.0;
  }

  // Resolve inexpensive disjoint and fully-contained cases before evaluating
  // analytical boundary terms.
  const auto tet_box = primal::compute_bounding_box(tet.geometry);
  const primal::BoundingBox<double, 3> sphere_box(Point3 {-1.0, -1.0, -1.0}, Point3 {1.0, 1.0, 1.0});
  bool all_vertices_inside = true;
  const primal::Sphere<double, 3> unit_sphere(Point3 {}, 1.0);
  for(int vertex = 0; vertex < NUM_VERTICES; ++vertex)
  {
    all_vertices_inside = all_vertices_inside && unit_sphere.contains(tet.geometry[vertex]);
  }
  if(!tet_box.intersectsWith(sphere_box))
  {
    return 0.0;
  }
  if(all_vertices_inside)
  {
    return tet.volume;
  }

  // Classify the tetrahedron entities touched by the sphere boundary.
  const IntersectionData intersections = findIntersections(tet);
  if(intersections.faces == 0 && tet.geometry.contains(Point3 {}))
  {
    return 4.0 * M_PI / 3.0;
  }
  if(intersections.vertices == 0 && intersections.edges == 0 && intersections.faces == 0)
  {
    return 0.0;
  }

  // Inclusion-exclusion starts with the full sphere, removes the cap outside
  // each face, restores regions removed twice along edges, and removes regions
  // restored too many times around vertices.
  double result = 4.0 * M_PI / 3.0;
  for(int face = 0; face < NUM_FACES; ++face)
  {
    if(isMarked(intersections.faces, face))
    {
      result -= capVolume(1.0 + tet.faces[face].normal.dot(-tet.faces[face].center));
    }
  }
  for(int edge = 0; edge < NUM_EDGES; ++edge)
  {
    if(isMarked(intersections.edges, edge))
    {
      result += edgeWedge(tet, intersections, edge);
    }
  }
  for(int vertex = 0; vertex < NUM_VERTICES; ++vertex)
  {
    if(isMarked(intersections.vertices, vertex))
    {
      result -= vertexCorrection(tet, intersections, vertex);
    }
  }

  // Clamp only roundoff-sized excursions beyond the physical volume range.
  const double maximum = axom::utilities::min(4.0 * M_PI / 3.0, tet.volume);
  const double tolerance = ::sqrt(axom::numeric_limits<double>::epsilon()) * maximum;
  if(result < 0.0 && result > -tolerance)
  {
    result = 0.0;
  }
  if(result > maximum && result - maximum < tolerance)
  {
    result = maximum;
  }
  SLIC_ASSERT(result >= 0.0 && result <= maximum);
  return result;
}

}  // namespace sphere_tet
}  // namespace detail

}  // namespace primal
}  // namespace axom
