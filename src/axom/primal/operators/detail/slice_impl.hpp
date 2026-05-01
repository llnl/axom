// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_SLICE_IMPL_HPP_
#define AXOM_PRIMAL_SLICE_IMPL_HPP_

namespace axom
{
namespace primal
{
namespace detail
{

template <typename PolygonType>
AXOM_HOST_DEVICE bool polygon_has_vertex(const PolygonType& poly,
                                         const typename PolygonType::PointType& pt,
                                         double eps = 1e-10)
{
  for(int i = 0; i < poly.numVertices(); ++i)
  {
    if(poly[i].isNearlyEqual(pt, eps))
    {
      return true;
    }
  }

  return false;
}

template <typename PolygonType>
AXOM_HOST_DEVICE void add_unique_vertex(PolygonType& poly,
                                        const typename PolygonType::PointType& pt,
                                        double eps = 1e-10)
{
  if(!polygon_has_vertex(poly, pt, eps))
  {
    poly.addVertex(pt);
  }
}

/*!
 * \brief Slices a 3D tetrahedron with a plane and returns the resulting
 *        polygon.
 *
 * \param [in] tet The tetrahedron to slice
 * \param [in] plane The slicing plane
 *
 * \return The polygon obtained from slicing a tetrahedron with a plane.
 */
template <typename T, PolygonArray ARRAY_TYPE = PolygonArray::Dynamic, int MAX_VERTS = DEFAULT_MAX_NUM_VERTICES>
AXOM_HOST_DEVICE primal::Polygon<T, 3, ARRAY_TYPE, MAX_VERTS> slice_tet_plane(
  const primal::Tetrahedron<T, 3>& tet,
  const primal::Plane<T, 3>& plane)
{
  Polygon<T, 3, ARRAY_TYPE, MAX_VERTS> intersectionPolygon;
  const int tetEdges[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
  T signedDistances[4];
  bool onPlane[4];

  if(!plane.isValid())
  {
    return intersectionPolygon;
  }

  // Compute the signed distance of each tet vertex to the plane once.
  for(int i = 0; i < 4; ++i)
  {
    signedDistances[i] = plane.signedDistance(tet[i]);
    onPlane[i] = axom::utilities::isNearlyEqual(signedDistances[i], T {0});
  }

  // Check each tet edge for a crossing or an endpoint that lies on the plane.
  for(int i = 0; i < 6; ++i)
  {
    const int v0 = tetEdges[i][0];
    const int v1 = tetEdges[i][1];
    const T sourceDistance = signedDistances[v0];
    const T targetDistance = signedDistances[v1];
    const bool sourceOnPlane = onPlane[v0];
    const bool targetOnPlane = onPlane[v1];

    if(sourceOnPlane && targetOnPlane)
    {
      add_unique_vertex(intersectionPolygon, tet[v0]);
      add_unique_vertex(intersectionPolygon, tet[v1]);
    }
    else if(sourceOnPlane)
    {
      add_unique_vertex(intersectionPolygon, tet[v0]);
    }
    else if(targetOnPlane)
    {
      add_unique_vertex(intersectionPolygon, tet[v1]);
    }
    else if((sourceDistance < T {0} && targetDistance > T {0}) ||
            (sourceDistance > T {0} && targetDistance < T {0}))
    {
      Segment<T, 3> edge(tet[v0], tet[v1]);
      T t {};
      if(primal::intersect(plane, edge, t))
      {
        add_unique_vertex(intersectionPolygon, edge.at(t));
      }
    }
  }
  SLIC_ASSERT(intersectionPolygon.numVertices() <= 4);

  // fix the polygon if it bowties
  if(intersectionPolygon.numVertices() == 4)
  {
    Segment<T, 3> seg1(intersectionPolygon[0], intersectionPolygon[1]);
    Segment<T, 3> seg2(intersectionPolygon[2], intersectionPolygon[3]);
    Point<T, 3> sp;
    if(!primal::intersect(seg1, seg2, sp))
    {
      axom::utilities::swap(intersectionPolygon[2], intersectionPolygon[3]);
    }
  }

  // For nondegenerate slices, orient the polygon so its normal follows the
  // slicing plane normal.
  if(intersectionPolygon.numVertices() >= 3 &&
     intersectionPolygon.normal().dot(plane.getNormal()) < T {0})
  {
    intersectionPolygon.reverseOrientation();
  }

  return intersectionPolygon;
}

}  // namespace detail
}  // namespace primal
}  // namespace axom
#endif
