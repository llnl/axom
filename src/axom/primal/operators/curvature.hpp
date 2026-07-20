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
    const T xp = Dt[0];  // x'
    const T yp = Dt[1];  // y'

    const T xpp = DtDt[0];  // x''
    const T ypp = DtDt[1];  // y''

    // This is signed curvature as formulated at:
    // https://en.wikipedia.org/wiki/Curvature#Curvature_of_a_graph
    // k = (x'y'' - y'x'') / pow(x'x' + y'y', 3./2.)
    const T xp2_plus_yp2 = xp * xp + yp * yp;
    return (xp * ypp - yp * xpp) / pow(xp2_plus_yp2, 3. / 2.);
  }
  else
  {
    return VectorType::cross_product(Dt, DtDt).norm() / pow(Dt.norm(), 3.);
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
  const T D1Norm = D1.norm();
  const T D1Norm3 = pow(D1Norm, 3.);
  const T D1Norm5 = pow(D1Norm, 5.);

  if constexpr(VectorType::dimension() == 2)
  {
    const T det12 = D1[0] * D2[1] - D1[1] * D2[0];
    const T det13 = D1[0] * D3[1] - D1[1] * D3[0];
    return det13 / D1Norm3 - 3. * det12 * D1.dot(D2) / D1Norm5;
  }
  else
  {
    const auto cross12 = VectorType::cross_product(D1, D2);
    const auto cross13 = VectorType::cross_product(D1, D3);
    const T cross12Norm = cross12.norm();
    const T crossTerm = cross12.dot(cross13);

    return crossTerm / (cross12Norm * D1Norm3) - 3. * cross12Norm * D1.dot(D2) / D1Norm5;
  }
}

}  // namespace primal
}  // namespace axom
