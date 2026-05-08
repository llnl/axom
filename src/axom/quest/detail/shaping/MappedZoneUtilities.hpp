// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_MAPPED_ZONE_UTILITIES_HPP_
#define AXOM_QUEST_MAPPED_ZONE_UTILITIES_HPP_

#include "axom/config.hpp"

#include "axom/core.hpp"
#include "axom/core/numerics/Determinants.hpp"
#include "axom/primal.hpp"

/*!
 * \file MappedZoneUtilities.hpp
 *
 * \brief Header-only utilities for evaluating low-order mapped quad/hex zones.
 */

namespace axom
{
namespace quest
{
namespace shaping
{
namespace detail
{

/*!
 * \brief Maps a point in the unit square to a physical quad using bilinear
 *        shape functions.
 *
 * \tparam ShapeType A zone type exposing `getId(i)` for its node ids.
 * \tparam CoordsetView A coordset view whose entries are point-like.
 *
 * \param [in] zone The source zone.
 * \param [in] coordsetView The coordinate storage for the source mesh.
 * \param [in] u The first reference-space coordinate in `[0,1]`.
 * \param [in] v The second reference-space coordinate in `[0,1]`.
 *
 * \return The mapped physical-space point.
 */
template <typename ShapeType, typename CoordsetView>
AXOM_HOST_DEVICE primal::Point<typename CoordsetView::value_type, 2> mapToPhysicalPoint(
  const ShapeType& zone,
  const CoordsetView& coordsetView,
  double u,
  double v)
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
 * \brief Maps a point in the unit cube to a physical hex using trilinear
 *        shape functions.
 *
 * \tparam ShapeType A zone type exposing `getId(i)` for its node ids.
 * \tparam CoordsetView A coordset view whose entries are point-like.
 *
 * \param [in] zone The source zone.
 * \param [in] coordsetView The coordinate storage for the source mesh.
 * \param [in] u The first reference-space coordinate in `[0,1]`.
 * \param [in] v The second reference-space coordinate in `[0,1]`.
 * \param [in] w The third reference-space coordinate in `[0,1]`.
 *
 * \return The mapped physical-space point.
 */
template <typename ShapeType, typename CoordsetView>
AXOM_HOST_DEVICE primal::Point<typename CoordsetView::value_type, 3> mapToPhysicalPoint(
  const ShapeType& zone,
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
 * \brief Evaluates the local physical area scale for a mapped quad.
 *
 * Computes the determinant of the 2x2 Jacobian for the bilinear map from the
 * unit square to the physical zone, returning its absolute value.
 *
 * \tparam ShapeType A zone type exposing `getId(i)` for its node ids.
 * \tparam CoordsetView A coordset view whose entries are point-like.
 *
 * \param [in] zone The source zone.
 * \param [in] coordsetView The coordinate storage for the source mesh.
 * \param [in] u The first reference-space coordinate in `[0,1]`.
 * \param [in] v The second reference-space coordinate in `[0,1]`.
 *
 * \return The local Jacobian area scale `|det(dx/du, dx/dv)|`.
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

  return axom::utilities::abs(
    axom::numerics::determinant(dxdu[0], dxdv[0], dxdu[1], dxdv[1]));
}

/*!
 * \brief Evaluates the local physical volume scale for a mapped hex.
 *
 * Computes the determinant of the 3x3 Jacobian for the trilinear map from the
 * unit cube to the physical zone, returning its absolute value.
 *
 * \tparam ShapeType A zone type exposing `getId(i)` for its node ids.
 * \tparam CoordsetView A coordset view whose entries are point-like.
 *
 * \param [in] zone The source zone.
 * \param [in] coordsetView The coordinate storage for the source mesh.
 * \param [in] u The first reference-space coordinate in `[0,1]`.
 * \param [in] v The second reference-space coordinate in `[0,1]`.
 * \param [in] w The third reference-space coordinate in `[0,1]`.
 *
 * \return The local Jacobian volume scale `|det(dx/du, dx/dv, dx/dw)|`.
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
    dxdu[d] = du0 * p0[d] + du1 * p1[d] + du2 * p2[d] + du3 * p3[d] + du4 * p4[d] +
      du5 * p5[d] + du6 * p6[d] + du7 * p7[d];
    dxdv[d] = dv0 * p0[d] + dv1 * p1[d] + dv2 * p2[d] + dv3 * p3[d] + dv4 * p4[d] +
      dv5 * p5[d] + dv6 * p6[d] + dv7 * p7[d];
    dxdw[d] = dw0 * p0[d] + dw1 * p1[d] + dw2 * p2[d] + dw3 * p3[d] + dw4 * p4[d] +
      dw5 * p5[d] + dw6 * p6[d] + dw7 * p7[d];
  }

  return axom::utilities::abs(VectorType::scalar_triple_product(dxdu, dxdv, dxdw));
}

/*!
 * \brief Returns the tensor-product quadrature point count per zone.
 *
 * \param [in] ruleX The rule in the first logical direction.
 * \param [in] ruleY The rule in the second logical direction.
 * \param [in] ruleZ The rule in the third logical direction.
 * \param [in] dim The logical dimension of the source zones.
 *
 * \return The number of quadrature points generated for one zone.
 */
inline int quadraturePointCount(const numerics::QuadratureRule& ruleX,
                                const numerics::QuadratureRule& ruleY,
                                const numerics::QuadratureRule& ruleZ,
                                int dim)
{
  return dim == 2 ? ruleX.getNumPoints() * ruleY.getNumPoints()
                  : ruleX.getNumPoints() * ruleY.getNumPoints() * ruleZ.getNumPoints();
}

}  // namespace detail
}  // namespace shaping
}  // namespace quest
}  // namespace axom

#endif
