// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file intersection_volume.hpp
 *
 * \brief Consists of a set of methods to find the volume of
 *        intersection (clipping) between a primal primitive and
 *        another primal primitive
 */

#ifndef AXOM_PRIMAL_INTERSECTION_VOLUME_HPP_
#define AXOM_PRIMAL_INTERSECTION_VOLUME_HPP_

#include "axom/primal/geometry/Tetrahedron.hpp"
#include "axom/primal/geometry/Octahedron.hpp"
#include "axom/primal/geometry/Polyhedron.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"

#include "axom/primal/operators/clip.hpp"

namespace axom
{
namespace primal
{
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

#endif  // AXOM_PRIMAL_INTERSECTION_VOLUME_HPP_
