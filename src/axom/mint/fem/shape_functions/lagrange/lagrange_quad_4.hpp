// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_QUAD4_HPP_
#define MINT_QUAD4_HPP_

// Mint includes
#include "axom/mint/mesh/CellTypes.hpp"
#include "axom/mint/fem/FEBasisTypes.hpp"
#include "axom/mint/fem/shape_functions/ShapeFunction.hpp"

// Slic includes
#include "axom/slic/interface/slic.hpp"

namespace axom
{
namespace mint
{
/*!
 * \brief Lagrange Finite Element definition for Bilinear Quadrilateral.
 *
 * \verbatim
 *
 * quad_4:
 *
 *   3       2
 *   +-------+
 *   |       |
 *   |       |
 *   +-------+
 *   0       1
 *
 * \endverbatim
 *
 * \see ShapeFunction
 */
template <>
class Lagrange<mint::QUAD> : public ShapeFunction<Lagrange<mint::QUAD>>
{
public:
  static CellType getCellType() { return mint::QUAD; }

  static int getType() { return MINT_LAGRANGE_BASIS; }

  static int getNumDofs() { return 4; }

  static int getMaxNewtonIters() { return 16; }

  static int getDimension() { return 2; }

  static double getMin() { return 0; }

  static double getMax() { return 1; }

  static void getCenter(double* center)
  {
    SLIC_ASSERT(center != nullptr);

    center[0] = 0.5;
    center[1] = 0.5;
  }

  static void getCoords(double* coords)
  {
    SLIC_ASSERT(coords != nullptr);

    // node 0
    coords[0] = 0.0;
    coords[1] = 0.0;

    // node 1
    coords[2] = 1.0;
    coords[3] = 0.0;

    // node 2
    coords[4] = 1.0;
    coords[5] = 1.0;

    // node 3
    coords[6] = 0.0;
    coords[7] = 1.0;
  }

  static void computeShape(const double* xr, double* phi)
  {
    SLIC_ASSERT(xr != nullptr);
    SLIC_ASSERT(phi != nullptr);

    const double r = xr[0];
    const double s = xr[1];
    const double rm = 1. - r;
    const double sm = 1. - s;

    phi[0] = rm * sm;
    phi[1] = r * sm;
    phi[2] = r * s;
    phi[3] = rm * s;
  }

  static void computeDerivatives(const double* xr, double* phidot)
  {
    SLIC_ASSERT(xr != nullptr);
    SLIC_ASSERT(phidot != nullptr);

    const double r = xr[0];
    const double s = xr[1];
    const double rm = 1. - r;
    const double sm = 1. - s;

    // r derivatives
    phidot[0] = -sm;
    phidot[1] = sm;
    phidot[2] = s;
    phidot[3] = -s;

    // s derivatives
    phidot[4] = -rm;
    phidot[5] = -r;
    phidot[6] = r;
    phidot[7] = rm;
  }
};

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_QUAD4_HPP_ */
