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

#ifndef AXOM_PRIMAL_CURVATURE_HPP_
#define AXOM_PRIMAL_CURVATURE_HPP_

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
T curvature(const VectorType &Dt, const VectorType &DtDt)
{
  if constexpr (VectorType::dimension() == 2)
  {
    const T xp = Dt[0];    // x'
    const T yp = Dt[1];    // y'

    const T xpp = DtDt[0]; // x''
    const T ypp = DtDt[1]; // y''
  
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
 * \brief Evaluates the curvature derivatives using supplied curve derivatives.
 *  
 * \param[in] d The number of derivatives to compute (1=1st deriv, 2=1st & 2nd derivs) 
 * \param[in] curveDerivs The derivatives, up to order \a d, of the curve.
 * \param[out] ders An array that will contain the curvature derivatives. 
 */ 
template <typename VectorType, typename T = typename VectorType::CoordType>
void curvatureDerivatives(int d,
                          const axom::Array<VectorType> &curveDerivs,
                          axom::Array<T> &ders)
{
  SLIC_ASSERT(d == 1 || d == 2);
  SLIC_ASSERT(curveDerivs.size() == 3);

  ders.resize(d);

  const VectorType& D1 = curveDerivs[0];
  const VectorType& D2 = curveDerivs[1];
  const VectorType& D3 = curveDerivs[2];

#if 0
  // Original 2D only
  const T xp = D1[0];    // x' 
  const T xpp = D2[0];   // x'' 
  const T xppp = D3[0];  // x''' 

  const T yp = D1[1];    // y' 
  const T ypp = D2[1];   // y'' 
  const T yppp = D3[1];  // y''' 
  
  // 1st derivative of curvature. 
  const T xp2_plus_yp2 = xp * xp + yp * yp; 
  const T A = -3. * (xp * ypp - yp * xpp) * 2. * (xp * xpp + yp * ypp); 
  const T B = 2. * pow(xp2_plus_yp2, 5. / 2.); 
  const T C = xp * yppp - yp * xppp; 
  const T D = pow(xp2_plus_yp2, 3. / 2.); 
  ders[0] = A / B + C / D; 
  
  if(d >= 2) 
  { 
    // 2nd derivative of curvature. 
    const T E = 15. * (-yp * xpp + xp * ypp) * 
      pow(2. * xp * xpp + 2. * yp * ypp, 2.) / 
      (4. * pow(xp2_plus_yp2, 7. / 2.)); 
    const T F = 3. * (2. * xp * xpp + 2. * yp * ypp) * 
      (-yp * xppp + xp * yppp) / pow(xp2_plus_yp2, 5. / 2.); 
    const T G = 3. * (-yp * xpp + xp * ypp) * 
      (2. * (xpp * xpp) + 2. * (ypp * ypp) + 2. * xp * xppp + 2. * yp * yppp) / 
      (2. * pow(xp2_plus_yp2, 5. / 2.)); 
    const T H = (-ypp * xppp + xpp * yppp) / pow(xp2_plus_yp2, 3. / 2.); 

    ders[1] = E - F - G + H; 
  }
#else
  const T D1Norm = D1.norm();
  const T D1Norm3 = pow(D1Norm, 3.);
  const T D1Norm5 = pow(D1Norm, 5.);
  const T D1D2Norm = VectorType::cross_product(D1, D2).norm();
  const T D1D3Norm = VectorType::cross_product(D1, D3).norm();

  // 1st derivative of curvature. 
  const T A = -3. * D1D2Norm * 2. * D1.dot(D2);
  const T B = 2. * D1Norm5;
  const T C = D1D3Norm;
  const T D = D1Norm3;
  ders[0] = A / B + C / D; 

  if(d >= 2) 
  { 
    // 2nd derivative of curvature. 
    const T E = 15. * D1D2Norm * pow(2. * D1.dot(D2), 2.) / (4. * pow(D1Norm, 7.));

    const T F = 3. * (2. * D1.dot(D2)) * D1D3Norm / D1Norm5;

    const T G = 3. * D1D2Norm * (D2.squared_norm() * D1.dot(D3)) / D1Norm5;

    const T H = VectorType::cross_product(D2, D3).norm() / D1Norm3;

    ders[1] = E - F - G + H; 
  }
#endif
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_CURVATURE_HPP_
