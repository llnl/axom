// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/*!
 * \file intersection_volume.hpp
 *
 * \brief Consists of a set of methods to find the volume of
 *        intersection (clipping) between a primal primitive and
 *        another primal primitive
 */

#include "axom/primal/geometry/Sphere.hpp"
#include "axom/primal/geometry/Octahedron.hpp"
#include "axom/primal/geometry/Polyhedron.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"

#include "axom/primal/operators/clip.hpp"
#include "axom/primal/operators/detail/intersection_volume_sphere_tet_impl.hpp"

namespace axom
{
namespace primal
{
/*!
 * \brief Finds the exact sphere-tetrahedron intersection volume.
 *
 * The computation normalizes the tetrahedron against a unit sphere and uses
 * inclusion-exclusion over spherical caps, edge wedges, and vertex cones.
 *
 * \tparam T the coordinate type
 * \param [in] sphere the sphere
 * \param [in] tet the tetrahedron
 * \return the unsigned intersection volume
 * \note This overload is available in host and device code.
 */
template <typename T>
AXOM_HOST_DEVICE T intersection_volume(const Sphere<T, 3>& sphere, const Tetrahedron<T, 3>& tet)
{
  const double radius = static_cast<double>(sphere.getRadius());
  if(radius <= 0.0)
  {
    return T {};
  }

  const auto& center = sphere.getCenter();
  axom::StackArray<detail::sphere_tet::Point3, detail::sphere_tet::NUM_VERTICES> vertices;
  for(int vertex = 0; vertex < detail::sphere_tet::NUM_VERTICES; ++vertex)
  {
    vertices[vertex] = detail::sphere_tet::Point3 {
      (static_cast<double>(tet[vertex][0]) - static_cast<double>(center[0])) / radius,
      (static_cast<double>(tet[vertex][1]) - static_cast<double>(center[1])) / radius,
      (static_cast<double>(tet[vertex][2]) - static_cast<double>(center[2])) / radius};
  }
  const detail::sphere_tet::Tetrahedron3 normalizedGeometry {vertices[0],
                                                             vertices[1],
                                                             vertices[2],
                                                             vertices[3]};
  const auto normalizedTet = detail::sphere_tet::makeTetCache(normalizedGeometry);
  const double normalizedVolume = detail::sphere_tet::unitSphereTetOverlap(normalizedTet);
  return static_cast<T>(normalizedVolume * radius * radius * radius);
}

/*!
 * \brief Finds the exact tetrahedron-sphere intersection volume.
 *
 * \tparam T the coordinate type
 * \param [in] tet the tetrahedron
 * \param [in] sphere the sphere
 * \return the unsigned intersection volume
 * \note This overload is available in host and device code.
 */
template <typename T>
AXOM_HOST_DEVICE T intersection_volume(const Tetrahedron<T, 3>& tet, const Sphere<T, 3>& sphere)
{
  return intersection_volume(sphere, tet);
}

/*!
 * \brief Finds the absolute (unsigned) intersection volume between
 *        a hexahedron and a tetrahedron
 *
 * \param [in] hex The hexahedron
 * \param [in] tet The tetrahedron
 * \param [in] eps The tolerance for determining the intersection
 * \param [in] tryFixOrientation If true, takes each shape with a negative
 *             signed volume and swaps the order of some vertices in that
 *             shape to try to obtain a nonnegative signed volume.
 *             Defaults to false.
 *
 * \return Intersection volume between the hexahedron and tetrahedron
 *
 * \warning tryFixOrientation flag does not guarantee the shapes' vertex orders
 *          will be valid. It is the responsiblity of the caller to pass
 *          shapes with a valid vertex order. Otherwise, if the shapes have
 *          invalid vertex orders, the returned volume may be zero
 *          and/or unexpected.
 *
 * \warning If tryFixOrientation flag is false and some of the shapes have
 *          a negative signed volume, the returned volume of intersection
 *          may be zero and/or unexpected.
 *
 */
template <typename T>
AXOM_HOST_DEVICE T intersection_volume(const Hexahedron<T, 3>& hex,
                                       const Tetrahedron<T, 3>& tet,
                                       double eps = 1.e-10,
                                       bool tryFixOrientation = false)
{
  return clip(hex, tet, eps, tryFixOrientation).volume();
}

/*!
 * \brief Finds the absolute (unsigned) intersection volume between
 *        a tetrahedron and a hexahedron
 *
 * \param [in] hex The tetrahedron
 * \param [in] tet The hexahedron
 * \param [in] eps The tolerance for determining the intersection
 * \param [in] tryFixOrientation If true, takes each shape with a negative
 *             signed volume and swaps the order of some vertices in that
 *             shape to try to obtain a nonnegative signed volume.
 *             Defaults to false.
 *
 * \return Intersection volume between the tetrahedron and hexahedron
 *
 * \warning tryFixOrientation flag does not guarantee the shapes' vertex orders
 *          will be valid. It is the responsiblity of the caller to pass
 *          shapes with a valid vertex order. Otherwise, if the shapes have
 *          invalid vertex orders, the returned volume may be zero
 *          and/or unexpected.
 *
 * \warning If tryFixOrientation flag is false and some of the shapes have
 *          a negative signed volume, the returned volume of intersection
 *          may be zero and/or unexpected.
 *
 */
template <typename T>
AXOM_HOST_DEVICE T intersection_volume(const Tetrahedron<T, 3>& tet,
                                       const Hexahedron<T, 3>& hex,
                                       double eps = 1.e-10,
                                       bool tryFixOrientation = false)
{
  return intersection_volume(hex, tet, eps, tryFixOrientation);
}

/*!
 * \brief Finds the absolute (unsigned) intersection volume between
 *        a octahedron and a tetrahedron
 *
 * \param [in] oct The octahedron
 * \param [in] tet The tetrahedron
 * \param [in] eps The tolerance for determining the intersection
 * \param [in] tryFixOrientation If true, takes each shape with a negative
 *             signed volume and swaps the order of some vertices in that
 *             shape to try to obtain a nonnegative signed volume.
 *             Defaults to false.
 *
 * \return Intersection volume between the octahedron and tetrahedron
 *
 * \warning tryFixOrientation flag does not guarantee the shapes' vertex orders
 *          will be valid. It is the responsiblity of the caller to pass
 *          shapes with a valid vertex order. Otherwise, if the shapes have
 *          invalid vertex orders, the returned volume may be zero
 *          and/or unexpected.
 *
 * \warning If tryFixOrientation flag is false and some of the shapes have
 *          a negative signed volume, the returned volume of intersection
 *          may be zero and/or unexpected.
 *
 */
template <typename T>
AXOM_HOST_DEVICE T intersection_volume(const Octahedron<T, 3>& oct,
                                       const Tetrahedron<T, 3>& tet,
                                       double eps = 1.e-10,
                                       bool tryFixOrientation = false)
{
  return clip(oct, tet, eps, tryFixOrientation).volume();
}

/*!
 * \brief Finds the absolute (unsigned) intersection volume between
 *        a tetrahedron and a octahedron
 *
 * \param [in] oct The tetrahedron
 * \param [in] tet The octahedron
 * \param [in] eps The tolerance for determining the intersection
 * \param [in] tryFixOrientation If true, takes each shape with a negative
 *             signed volume and swaps the order of some vertices in that
 *             shape to try to obtain a nonnegative signed volume.
 *             Defaults to false.
 *
 * \return Intersection volume between the tetrahedron and octahedron
 *
 * \warning tryFixOrientation flag does not guarantee the shapes' vertex orders
 *          will be valid. It is the responsiblity of the caller to pass
 *          shapes with a valid vertex order. Otherwise, if the shapes have
 *          invalid vertex orders, the returned volume may be zero
 *          and/or unexpected.
 *
 * \warning If tryFixOrientation flag is false and some of the shapes have
 *          a negative signed volume, the returned volume of intersection
 *          may be zero and/or unexpected.
 *
 */
template <typename T>
AXOM_HOST_DEVICE T intersection_volume(const Tetrahedron<T, 3>& tet,
                                       const Octahedron<T, 3>& oct,
                                       double eps = 1.e-10,
                                       bool tryFixOrientation = false)
{
  return intersection_volume(oct, tet, eps, tryFixOrientation);
}

/*!
 * \brief Finds the absolute (unsigned) intersection volume between
 *        a tetrahedron and another tetrahedron
 *
 * \param [in] tet1 The tetrahedron
 * \param [in] tet2 The other tetrahedron
 * \param [in] eps The tolerance for determining the intersection
 * \param [in] tryFixOrientation If true, takes each shape with a negative
 *             signed volume and swaps the order of some vertices in that
 *             shape to try to obtain a nonnegative signed volume.
 *             Defaults to false.
 *
 * \return Intersection volume between the tetrahedra
 *
 * \warning tryFixOrientation flag does not guarantee the shapes' vertex orders
 *          will be valid. It is the responsiblity of the caller to pass
 *          shapes with a valid vertex order. Otherwise, if the shapes have
 *          invalid vertex orders, the returned volume may be zero
 *          and/or unexpected.
 *
 * \warning If tryFixOrientation flag is false and some of the shapes have
 *          a negative signed volume, the returned volume of intersection
 *          may be zero and/or unexpected.
 *
 */
template <typename T>
AXOM_HOST_DEVICE T intersection_volume(const Tetrahedron<T, 3>& tet1,
                                       const Tetrahedron<T, 3>& tet2,
                                       double eps = 1.e-10,
                                       bool tryFixOrientation = false)
{
  return clip(tet1, tet2, eps, tryFixOrientation).volume();
}

}  // namespace primal
}  // namespace axom
