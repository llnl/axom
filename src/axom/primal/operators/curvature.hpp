// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file curvature.hpp
 *
 * \brief Consists of a set of templated routines used to calculate
 *  the curvature, given a curve's derivatives evaluated at a point.
 *
 */

#pragma once

#include "axom/config.hpp"
#include "axom/core/Array.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include "axom/slic/interface/slic.hpp"

namespace axom
{
namespace primal
{

/*! 
 * \brief Evaluates the curvature, given the derivatives.
 *  
 * \param[in] Dt The 1st derivative of the curve.
 * \param[in] DtDt The 2nd derivative of the curve.
 *
 * \return The curvature evaluated at \a u.
 */
template <typename VectorType, typename T = typename VectorType::CoordType>
T curvature(const VectorType& Dt, const VectorType& DtDt)
{
  if constexpr(VectorType::dimension() == 2)
  {
    const double xp = Dt[0];  // x'
    const double yp = Dt[1];  // y'

    const double xpp = DtDt[0];  // x''
    const double ypp = DtDt[1];  // y''

    // This is signed curvature as formulated at:
    // https://en.wikipedia.org/wiki/Curvature#Curvature_of_a_graph
    // k = (x'y'' - y'x'') / pow(x'x' + y'y', 3./2.)
    const double xp2_plus_yp2 = xp * xp + yp * yp;
    SLIC_ASSERT(xp2_plus_yp2 > 0.0);
    return static_cast<T>((xp * ypp - yp * xpp) / (xp2_plus_yp2 * sqrt(xp2_plus_yp2)));
  }
  else
  {
    const double dt_norm = Dt.norm();
    SLIC_ASSERT(dt_norm != 0.0);
    return static_cast<T>(VectorType::cross_product(Dt, DtDt).norm() /
                          (dt_norm * dt_norm * dt_norm));
  }
}

/*!
 * \brief Evaluates the first parameter derivative of curvature using supplied
 *        curve derivatives.
 *
 * \param[in] D1 The 1st derivative of the curve.
 * \param[in] D2 The 2nd derivative of the curve.
 * \param[in] D3 The 3rd derivative of the curve.
 *
 * \return The curvature derivative evaluated with respect to the curve
 *         parameter.
 */
template <typename VectorType, typename T = typename VectorType::CoordType>
T curvatureDerivative(const VectorType& D1, const VectorType& D2, const VectorType& D3)
{
  const double d1_norm = D1.norm();
  SLIC_ASSERT(d1_norm != 0.0);
  const double d1_norm3 = d1_norm * d1_norm * d1_norm;
  const double d1_norm5 = d1_norm3 * d1_norm * d1_norm;

  if constexpr(VectorType::dimension() == 2)
  {
    const double det12 = D1[0] * D2[1] - D1[1] * D2[0];
    const double det13 = D1[0] * D3[1] - D1[1] * D3[0];
    const double d1_dot_d2 = D1.dot(D2);
    return static_cast<T>(det13 / d1_norm3 - 3.0 * det12 * d1_dot_d2 / d1_norm5);
  }
  else
  {
    const auto cross12 = VectorType::cross_product(D1, D2);
    const auto cross13 = VectorType::cross_product(D1, D3);
    const double cross12_norm = cross12.norm();
    const double cross_term = cross12.dot(cross13);
    const double d1_dot_d2 = D1.dot(D2);
    return static_cast<T>(cross_term / (cross12_norm * d1_norm3) -
                          3.0 * cross12_norm * d1_dot_d2 / d1_norm5);
  }
}

}  // namespace primal
}  // namespace axom
