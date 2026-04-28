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

  // find intersection vertices
  for(int i = 0; i < 4; ++i)
  {
    for(int j = i + 1; j < 4; ++j)
    {
      Segment<T, 3> edge(tet[i], tet[j]);
      T t {};
      if(primal::intersect(plane, edge, t))
      {
        intersectionPolygon.addVertex(edge.at(t));
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
  return intersectionPolygon;
}

}  // namespace detail
}  // namespace primal
}  // namespace axom
#endif
