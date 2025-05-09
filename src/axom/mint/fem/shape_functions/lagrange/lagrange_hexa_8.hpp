// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_HEXA_8_HPP_
#define MINT_HEXA_8_HPP_

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
 * \brief Lagrange Finite Element definition for the Linear Hexahedron
 *
 * \verbatim
 *
 * hexa_8:
 *
 * 7 +---------+ 6
 *   |\        |\
 *   |  \      |  \
 *   | 4 + --------+ 5
 * 3 +---|-----+ 2 |
 *   \   |     \   |
 *    \  |      \  |
 *     \ |       \ |
 *   0  +----------+ 1
 *
 * \endverbatim
 *
 * \see ShapeFunction
 */
template <>
class Lagrange<mint::HEX> : public ShapeFunction<Lagrange<mint::HEX>>
{
public:
  static CellType getCellType() { return mint::HEX; }

  static int getType() { return MINT_LAGRANGE_BASIS; }

  static int getNumDofs() { return 8; }

  static int getMaxNewtonIters() { return 16; }

  static int getDimension() { return 3; }

  static double getMin() { return 0; }

  static double getMax() { return 1; }

  static void getCenter(double* center)
  {
    SLIC_ASSERT(center != nullptr);
    center[0] = center[1] = center[2] = 0.5;
  }

  static void getCoords(double* coords)
  {
    SLIC_ASSERT(coords != nullptr);

    // bottom:
    // node 0
    coords[0] = 0.0;
    coords[1] = 0.0;
    coords[2] = 0.0;

    // node 1
    coords[3] = 1.0;
    coords[4] = 0.0;
    coords[5] = 0.0;

    // node 2
    coords[6] = 1.0;
    coords[7] = 1.0;
    coords[8] = 0.0;

    // node 3
    coords[9] = 0.0;
    coords[10] = 1.0;
    coords[11] = 0.0;

    // top:
    // node 4
    coords[12] = 0.0;
    coords[13] = 0.0;
    coords[14] = 1.0;

    // node 5
    coords[15] = 1.0;
    coords[16] = 0.0;
    coords[17] = 1.0;

    // node 6
    coords[18] = 1.0;
    coords[19] = 1.0;
    coords[20] = 1.0;

    // node 7
    coords[21] = 0.0;
    coords[22] = 1.0;
    coords[23] = 1.0;
  }

  static void computeShape(const double* xr, double* phi)
  {
    SLIC_ASSERT(xr != nullptr);
    SLIC_ASSERT(phi != nullptr);

    const double r = xr[0];
    const double s = xr[1];
    const double t = xr[2];
    const double rm = 1. - r;
    const double sm = 1. - s;
    const double tm = 1. - t;

    const double rm_x_sm = rm * sm;
    const double r_x_sm = r * sm;
    const double r_x_s = r * s;
    const double rm_x_s = rm * s;

    phi[0] = rm_x_sm * tm;
    phi[1] = r_x_sm * tm;
    phi[2] = r_x_s * tm;
    phi[3] = rm_x_s * tm;

    phi[4] = rm_x_sm * t;
    phi[5] = r_x_sm * t;
    phi[6] = r_x_s * t;
    phi[7] = rm_x_s * t;
  }

  static void computeDerivatives(const double* xr, double* phidot)
  {
    SLIC_ASSERT(xr != nullptr);
    SLIC_ASSERT(phidot != nullptr);

    const double r = xr[0];
    const double s = xr[1];
    const double t = xr[2];
    const double rm = 1. - r;
    const double sm = 1. - s;
    const double tm = 1. - t;

    // r derivatives
    phidot[0] = -sm * tm;
    phidot[1] = sm * tm;
    phidot[2] = s * tm;
    phidot[3] = -s * tm;

    phidot[4] = -sm * t;
    phidot[5] = sm * t;
    phidot[6] = s * t;
    phidot[7] = -s * t;

    // s derivatives
    phidot[8] = -rm * tm;
    phidot[9] = -r * tm;
    phidot[10] = r * tm;
    phidot[11] = rm * tm;

    phidot[12] = -rm * t;
    phidot[13] = -r * t;
    phidot[14] = r * t;
    phidot[15] = rm * t;

    // t derivatives
    phidot[16] = -rm * sm;
    phidot[17] = -r * sm;
    phidot[18] = -r * s;
    phidot[19] = -rm * s;

    phidot[20] = rm * sm;
    phidot[21] = r * sm;
    phidot[22] = r * s;
    phidot[23] = rm * s;
  }
};

} /* namespace mint */
} /* namespace axom */
#endif /* MINT_HEXA_8_HPP_ */
