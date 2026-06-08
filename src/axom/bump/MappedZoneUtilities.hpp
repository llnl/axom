// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_MAPPED_ZONE_UTILITIES_HPP_
#define AXOM_BUMP_MAPPED_ZONE_UTILITIES_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/core/numerics/Determinants.hpp"
#include "axom/primal.hpp"

namespace axom
{
namespace bump
{
namespace detail
{

/*!
 * \file MappedZoneUtilities.hpp
 *
 * \brief Helper functions for mapping low-order quad and hex zones from
 *        reference-space quadrature coordinates into physical space.
 */

/*!
 * \brief Map a quadrilateral's reference-space point to physical space.
 *
 * \param zone The quadrilateral zone to map.
 * \param coordsetView The coordset view that provides zone vertices.
 * \param u The first reference-space coordinate in [0, 1].
 * \param v The second reference-space coordinate in [0, 1].
 *
 * \return The physical point corresponding to \a (u, v).
 */
template <typename ShapeType, typename CoordsetView>
AXOM_HOST_DEVICE primal::Point<typename CoordsetView::value_type, 2>
mapToPhysicalPoint(const ShapeType& zone, const CoordsetView& coordsetView, double u, double v)
{
  using PointType = primal::Point<typename CoordsetView::value_type, 2>;
  const auto p0 = coordsetView[zone.getId(0)];
  const auto p1 = coordsetView[zone.getId(1)];
  const auto p2 = coordsetView[zone.getId(2)];
  const auto p3 = coordsetView[zone.getId(3)];

  const double n0 = (1.0 - u) * (1.0 - v);
  const double n1 = u * (1.0 - v);
  const double n2 = u * v;
  const double n3 = (1.0 - u) * v;

  PointType pt;
  for(int d = 0; d < 2; ++d)
  {
    pt[d] = n0 * p0[d] + n1 * p1[d] + n2 * p2[d] + n3 * p3[d];
  }
  return pt;
}

/*!
 * \brief Map a hexahedron's reference-space point to physical space.
 *
 * \param zone The hexahedral zone to map.
 * \param coordsetView The coordset view that provides zone vertices.
 * \param u The first reference-space coordinate in [0, 1].
 * \param v The second reference-space coordinate in [0, 1].
 * \param w The third reference-space coordinate in [0, 1].
 *
 * \return The physical point corresponding to \a (u, v, w).
 */
template <typename ShapeType, typename CoordsetView>
AXOM_HOST_DEVICE primal::Point<typename CoordsetView::value_type, 3>
mapToPhysicalPoint(const ShapeType& zone,
                   const CoordsetView& coordsetView,
                   double u,
                   double v,
                   double w)
{
  using PointType = primal::Point<typename CoordsetView::value_type, 3>;
  const auto p0 = coordsetView[zone.getId(0)];
  const auto p1 = coordsetView[zone.getId(1)];
  const auto p2 = coordsetView[zone.getId(2)];
  const auto p3 = coordsetView[zone.getId(3)];
  const auto p4 = coordsetView[zone.getId(4)];
  const auto p5 = coordsetView[zone.getId(5)];
  const auto p6 = coordsetView[zone.getId(6)];
  const auto p7 = coordsetView[zone.getId(7)];

  const double a = 1.0 - u;
  const double b = 1.0 - v;
  const double c = 1.0 - w;

  const double n0 = a * b * c;
  const double n1 = u * b * c;
  const double n2 = u * v * c;
  const double n3 = a * v * c;
  const double n4 = a * b * w;
  const double n5 = u * b * w;
  const double n6 = u * v * w;
  const double n7 = a * v * w;

  PointType pt;
  for(int d = 0; d < 3; ++d)
  {
    pt[d] = n0 * p0[d] + n1 * p1[d] + n2 * p2[d] + n3 * p3[d] + n4 * p4[d] + n5 * p5[d] +
      n6 * p6[d] + n7 * p7[d];
  }
  return pt;
}

/*!
 * \brief Compute the quadrilateral area scale at a reference-space point.
 *
 * \param zone The quadrilateral zone to evaluate.
 * \param coordsetView The coordset view that provides zone vertices.
 * \param u The first reference-space coordinate in [0, 1].
 * \param v The second reference-space coordinate in [0, 1].
 *
 * \return The absolute Jacobian determinant at \a (u, v).
 */
template <typename ShapeType, typename CoordsetView>
AXOM_HOST_DEVICE double computePhysicalMeasureFactor(const ShapeType& zone,
                                                     const CoordsetView& coordsetView,
                                                     double u,
                                                     double v)
{
  using VectorType = primal::Vector<typename CoordsetView::value_type, 2>;

  const auto p0 = coordsetView[zone.getId(0)];
  const auto p1 = coordsetView[zone.getId(1)];
  const auto p2 = coordsetView[zone.getId(2)];
  const auto p3 = coordsetView[zone.getId(3)];

  const double du0 = -(1.0 - v);
  const double du1 = (1.0 - v);
  const double du2 = v;
  const double du3 = -v;

  const double dv0 = -(1.0 - u);
  const double dv1 = -u;
  const double dv2 = u;
  const double dv3 = 1.0 - u;

  VectorType dxdu;
  VectorType dxdv;
  for(int d = 0; d < 2; ++d)
  {
    dxdu[d] = du0 * p0[d] + du1 * p1[d] + du2 * p2[d] + du3 * p3[d];
    dxdv[d] = dv0 * p0[d] + dv1 * p1[d] + dv2 * p2[d] + dv3 * p3[d];
  }

  return axom::utilities::abs(axom::numerics::determinant(dxdu[0], dxdv[0], dxdu[1], dxdv[1]));
}

/*!
 * \brief Compute the hexahedral volume scale at a reference-space point.
 *
 * \param zone The hexahedral zone to evaluate.
 * \param coordsetView The coordset view that provides zone vertices.
 * \param u The first reference-space coordinate in [0, 1].
 * \param v The second reference-space coordinate in [0, 1].
 * \param w The third reference-space coordinate in [0, 1].
 *
 * \return The absolute Jacobian determinant at \a (u, v, w).
 */
template <typename ShapeType, typename CoordsetView>
AXOM_HOST_DEVICE double computePhysicalMeasureFactor(const ShapeType& zone,
                                                     const CoordsetView& coordsetView,
                                                     double u,
                                                     double v,
                                                     double w)
{
  using VectorType = primal::Vector<typename CoordsetView::value_type, 3>;

  const auto p0 = coordsetView[zone.getId(0)];
  const auto p1 = coordsetView[zone.getId(1)];
  const auto p2 = coordsetView[zone.getId(2)];
  const auto p3 = coordsetView[zone.getId(3)];
  const auto p4 = coordsetView[zone.getId(4)];
  const auto p5 = coordsetView[zone.getId(5)];
  const auto p6 = coordsetView[zone.getId(6)];
  const auto p7 = coordsetView[zone.getId(7)];

  const double a = 1.0 - u;
  const double b = 1.0 - v;
  const double c = 1.0 - w;

  const double du0 = -b * c;
  const double du1 = b * c;
  const double du2 = v * c;
  const double du3 = -v * c;
  const double du4 = -b * w;
  const double du5 = b * w;
  const double du6 = v * w;
  const double du7 = -v * w;

  const double dv0 = -a * c;
  const double dv1 = -u * c;
  const double dv2 = u * c;
  const double dv3 = a * c;
  const double dv4 = -a * w;
  const double dv5 = -u * w;
  const double dv6 = u * w;
  const double dv7 = a * w;

  const double dw0 = -a * b;
  const double dw1 = -u * b;
  const double dw2 = -u * v;
  const double dw3 = -a * v;
  const double dw4 = a * b;
  const double dw5 = u * b;
  const double dw6 = u * v;
  const double dw7 = a * v;

  VectorType dxdu;
  VectorType dxdv;
  VectorType dxdw;
  for(int d = 0; d < 3; ++d)
  {
    dxdu[d] = du0 * p0[d] + du1 * p1[d] + du2 * p2[d] + du3 * p3[d] + du4 * p4[d] + du5 * p5[d] +
      du6 * p6[d] + du7 * p7[d];
    dxdv[d] = dv0 * p0[d] + dv1 * p1[d] + dv2 * p2[d] + dv3 * p3[d] + dv4 * p4[d] + dv5 * p5[d] +
      dv6 * p6[d] + dv7 * p7[d];
    dxdw[d] = dw0 * p0[d] + dw1 * p1[d] + dw2 * p2[d] + dw3 * p3[d] + dw4 * p4[d] + dw5 * p5[d] +
      dw6 * p6[d] + dw7 * p7[d];
  }

  return axom::utilities::abs(VectorType::scalar_triple_product(dxdu, dxdv, dxdw));
}

}  // namespace detail
}  // namespace bump
}  // namespace axom

#endif
