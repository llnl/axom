// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_QUAD9_HPP_
#define MINT_QUAD9_HPP_

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
 * \brief Lagrange Finite Element definition for the Quadratic Quadrilateral.
 *
 * \note The nodes are numbered in the following order: corners(4), edges(4)
 *  and last the centroid (1) as depicted below.
 *
 * \verbatim
 *
 * quad_9:
 *
 *   3   6    2
 *   +---+----+
 *   |        |
 * 7 +   + 8  + 5
 *   |        |
 *   +---+----+
 *   0   4    1
 *
 * \endverbatim
 *
 * \see ShapeFunction
 */
template <>
class Lagrange<mint::QUAD9> : public ShapeFunction<Lagrange<mint::QUAD9>>
{
public:
  static CellType getCellType() { return mint::QUAD9; }

  static int getType() { return MINT_LAGRANGE_BASIS; }

  static int getNumDofs() { return 9; }

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

    // corners:
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

    // edges:
    // node 4
    coords[8] = 0.5;
    coords[9] = 0.0;

    // node 5
    coords[10] = 1.0;
    coords[11] = 0.5;

    // node 6
    coords[12] = 0.5;
    coords[13] = 1.0;

    // node 7
    coords[14] = 0.0;
    coords[15] = 0.5;

    // centroid:
    // node 8
    coords[16] = 0.5;
    coords[17] = 0.5;
  }

  static void computeShape(const double* xr, double* phi)
  {
    SLIC_ASSERT(xr != nullptr);
    SLIC_ASSERT(phi != nullptr);

    const double r = xr[0];
    const double s = xr[1];

    // bar element along r
    const double r1 = (r - 1.) * (2. * r - 1);
    const double r2 = 4. * r * (1. - r);
    const double r3 = r * (2. * r - 1.);

    // bar element along s
    const double s1 = (s - 1.) * (2. * s - 1);
    const double s2 = 4. * s * (1. - s);
    const double s3 = s * (2. * s - 1.);

    phi[0] = r1 * s1;
    phi[1] = r3 * s1;
    phi[2] = r3 * s3;
    phi[3] = r1 * s3;
    phi[4] = r2 * s1;
    phi[5] = r3 * s2;
    phi[6] = r2 * s3;
    phi[7] = r1 * s2;
    phi[8] = r2 * s2;
  }

  static void computeDerivatives(const double* xr, double* phidot)
  {
    SLIC_ASSERT(xr != nullptr);
    SLIC_ASSERT(phidot != nullptr);

    const double r = xr[0];
    const double s = xr[1];

    // bar element along r
    const double r1 = (r - 1.) * (2. * r - 1);
    const double r2 = 4. * r * (1. - r);
    const double r3 = r * (2. * r - 1.);

    // bar element along s
    const double s1 = (s - 1.) * (2. * s - 1);
    const double s2 = 4. * s * (1. - s);
    const double s3 = s * (2. * s - 1.);

    // bar element derivatives along r
    const double dr1 = 4. * r - 3.;
    const double dr2 = 4. - 8. * r;
    const double dr3 = 4. * r - 1.;

    // bar element derivatives along s
    const double ds1 = 4. * s - 3.;
    const double ds2 = 4. - 8. * s;
    const double ds3 = 4. * s - 1.;

    // r derivatives
    phidot[0] = dr1 * s1;
    phidot[1] = dr3 * s1;
    phidot[2] = dr3 * s3;
    phidot[3] = dr1 * s3;
    phidot[4] = dr2 * s1;
    phidot[5] = dr3 * s2;
    phidot[6] = dr2 * s3;
    phidot[7] = dr1 * s2;
    phidot[8] = dr2 * s2;

    // s derivatives
    phidot[9] = r1 * ds1;
    phidot[10] = r3 * ds1;
    phidot[11] = r3 * ds3;
    phidot[12] = r1 * ds3;
    phidot[13] = r2 * ds1;
    phidot[14] = r3 * ds2;
    phidot[15] = r2 * ds3;
    phidot[16] = r1 * ds2;
    phidot[17] = r2 * ds2;
  }
};

} /* namespace mint */
}  // namespace axom

#endif /* MINT_QUAD_9_HPP_ */
