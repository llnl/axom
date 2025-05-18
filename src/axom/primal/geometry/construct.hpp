// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_CONSTRUCT_HPP_
#define AXOM_PRIMAL_CONSTRUCT_HPP_

#include "axom/config.hpp"
#include "axom/slic.hpp"
#include "axom/core/numerics/Matrix.hpp"
#include "axom/core/numerics/matvecops.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/Polyhedron.hpp"

#include <cmath>

namespace axom
{
namespace primal
{

/*!
 * \brief Create a regular polygon.
 *
 * \param nSides The number of sides in the polygon.
 * \param radius The radius of the polgon from the origin to a point.
 * \param transform An optional transformation matrix, which defaults to the identity
 *                  matrix for the dimension NDIMS.
 *
 * \return A new regular polygon with \a nSides sides and radius \radius,
 *         transformed by the supplied tranformation matrix.
 *
 * \note This is host-only function.
 */
template <typename T, int NDIMS, PolygonArray ARRAY_TYPE = PolygonArray::Dynamic, int MAX_VERTS = DEFAULT_MAX_NUM_VERTICES>
Polygon<T, NDIMS, ARRAY_TYPE, MAX_VERTS> regular_polygon(
  int nSides,
  T radius = T {1},
  const axom::numerics::Matrix<T>& transform = axom::numerics::Matrix<T>::identity(NDIMS))
{
  using PointType = typename axom::primal::Point<T, NDIMS>;

  SLIC_ASSERT(nSides >= 3);
  SLIC_ASSERT(transform.getNumRows() == transform.getNumColumns());

  const double dA = (2. * M_PI) / static_cast<double>(nSides);
  const double a0 = dA * 0.5 - M_PI / 2.;
  double a = a0;
  Polygon<T, NDIMS, ARRAY_TYPE, MAX_VERTS> poly;
  for(int s = 0; s < nSides; s++)
  {
    PointType pt {radius * static_cast<T>(cos(a)), radius * static_cast<T>(sin(a))};

    // Add the transformed point to the polygon.
    poly.addVertex(transform_point(pt, transform));

    a += dA;
  }

  return poly;
}

/*!
 * \brief Construct a regular 3D prism.
 *
 * \param numSides The number of sides on the base polygon.
 * \param radius The radius of the base polygon.
 * \param height The height of the prism.
 * \param transform An optional 4x4 matrix transformation.
 *
 * \return A Polyhedron that represents the regular prism.
 */
template <typename T>
Polyhedron<T, 3> regular_prism(
  int numSides,
  T radius = 1,
  T height = 1,
  const axom::numerics::Matrix<T>& transform = axom::numerics::Matrix<T>::identity(4))
{
  constexpr int MAX_VERTS = Polyhedron<T, 3>::MAX_VERTS;

  SLIC_ASSERT(transform.getNumRows() == 4);
  SLIC_ASSERT(transform.getNumRows() == transform.getNumColumns());

  Polyhedron<T, 3> poly;

  // Add vertices for base and top polygons to the polyhedron.
  std::int8_t verts[MAX_VERTS];
  const T zero {0};
  for(int zi = 0, vIndex = 0; zi < 2; zi++)
  {
    const T z = (zi == 1) ? height : zero;
    // Make a combined transform, xform.
    const auto translation = axom::numerics::transforms::translate<T>(zero, zero, z);
    auto xform = axom::numerics::Matrix<T>::zeros(4, 4);
    axom::numerics::matrix_multiply(transform, translation, xform);

    // Make a polygon and add its points to the polyhedron.
    const auto base = regular_polygon<T, 3>(numSides, radius, xform);
    for(int s = 0; s < numSides; s++)
    {
      verts[vIndex] = static_cast<std::int8_t>(poly.addVertex(base[s]));
      // Increment if we can fit more.
      vIndex = ((vIndex + 1) < MAX_VERTS) ? (vIndex + 1) : vIndex;
    }
  }

  // Define the vertex neighbors.
  for(int s = 0; s < numSides; s++)
  {
    int current = s;
    int prev = (s == 0) ? (numSides - 1) : (s - 1);
    int next = (s == (numSides - 1)) ? 0 : (s + 1);
    int z = s + numSides;
    poly.addNeighbors(current, {verts[next], verts[z], verts[prev]});
  }
  for(int s = 0; s < numSides; s++)
  {
    int current = numSides + s;
    int prev = numSides + ((s == 0) ? (numSides - 1) : (s - 1));
    int next = numSides + ((s == (numSides - 1)) ? 0 : (s + 1));
    int z = s;
    poly.addNeighbors(current, {verts[prev], verts[z], verts[next]});
  }

  return poly;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_CONSTRUCT_HPP_
