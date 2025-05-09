// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_TETRA_4_HPP_
#define MINT_TETRA_4_HPP_

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
 * \brief Lagrange Finite Element definition for the Linear Tetrahedron
 *
 * \verbatim
 *
 * tetra_4:
 *
 *        3
 *        x
 *       /|\
 *      / | \
 *     /  |  \
 *    /   |   \
 * 0 x----|----x 2
 *    \   |   /
 *     \  |  /
 *      \ | /
 *        x
 *        1
 *
 * \endverbatim
 *
 * \see ShapeFunction
 */
template <>
class Lagrange<mint::TET> : public ShapeFunction<Lagrange<mint::TET>>
{
public:
  static CellType getCellType() { return mint::TET; }

  static int getType() { return MINT_LAGRANGE_BASIS; }

  static int getNumDofs() { return 4; }

  static int getMaxNewtonIters() { return 16; }

  static int getDimension() { return 3; }

  static double getMin() { return 0; }

  static double getMax() { return 1; }

  static void getCenter(double* center)
  {
    SLIC_ASSERT(center != nullptr);
    center[0] = center[1] = center[2] = 0.25;
  }

  static void getCoords(double* coords)
  {
    SLIC_ASSERT(coords != nullptr);

    // node 0
    coords[0] = 0.0;
    coords[1] = 0.0;
    coords[2] = 0.0;

    // node 1
    coords[3] = 1.0;
    coords[4] = 0.0;
    coords[5] = 0.0;

    // node 2
    coords[6] = 0.0;
    coords[7] = 1.0;
    coords[8] = 0.0;

    // node 3
    coords[9] = 0.0;
    coords[10] = 0.0;
    coords[11] = 1.0;
  }

  static void computeShape(const double* xr, double* phi)
  {
    SLIC_ASSERT(xr != nullptr);
    SLIC_ASSERT(phi != nullptr);

    const double r = xr[0];
    const double s = xr[1];
    const double t = xr[2];

    phi[0] = 1 - r - s - t;
    phi[1] = r;
    phi[2] = s;
    phi[3] = t;
  }

  static void computeDerivatives(const double* AXOM_UNUSED_PARAM(xr), double* phidot)
  {
    SLIC_ASSERT(phidot != nullptr);

    // r derivatives
    phidot[0] = -1;
    phidot[1] = 1;
    phidot[2] = 0;
    phidot[3] = 0;

    // s derivatives
    phidot[4] = -1;
    phidot[5] = 0;
    phidot[6] = 1;
    phidot[7] = 0;

    // t derivatives
    phidot[8] = -1;
    phidot[9] = 0;
    phidot[10] = 0;
    phidot[11] = 1;
  }
};

} /* namespace mint */
} /* namespace axom */
#endif /* MINT_TETRA_4_HPP_ */
