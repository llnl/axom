// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_LAGRANGE_HEXA_27_HPP_
#define MINT_LAGRANGE_HEXA_27_HPP_

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
 * \brief Lagrange Finite Element definition for the triquadratic Hexahedron.
 *
 * \note The nodes are numbered in the following order: corners(8),
 *  edges(12), faces(6) and last the centroid (1) as depicted below.
 *
 * \verbatim
 *
 * hexa_27:
 *
 *  Top:
 *
 *    7   14    6
 *    +---+----+
 *    |        |
 * 15 +   + 25 + 13
 *    |        |
 *    +---+----+
 *    4   12    5
 *
 *  Middle:
 *
 *  19   23    18
 *    +---+----+
 *    |        |
 * 20 +   + 26 + 21
 *    |        |
 *    +---+----+
 *   16   22    17
 *
 *
 *  Bottom:
 *
 *    3   10   2
 *    +---+----+
 *    |        |
 * 11 +   + 24 + 9
 *    |        |
 *    +---+----+
 *    0   8    1
 *
 * \endverbatim
 *
 * \see ShapeFunction
 */
template <>
class Lagrange<mint::HEX27> : public ShapeFunction<Lagrange<mint::HEX27>>
{
public:
  static CellType getCellType() { return mint::HEX27; }

  static int getType() { return MINT_LAGRANGE_BASIS; }

  static int getNumDofs() { return 27; }

  static int getMaxNewtonIters() { return 16; }

  static int getDimension() { return 3; }

  static double getMin() { return 0; }

  static double getMax() { return 1; }

  static void getCenter(double* center)
  {
    SLIC_ASSERT(center != nullptr);

    center[0] = 0.5;
    center[1] = 0.5;
    center[2] = 0.5;
  }

  static void getCoords(double* coords)
  {
    SLIC_ASSERT(coords != nullptr);

    // corner nodes
    coords[0] = 0.0;
    coords[1] = 0.0;
    coords[2] = 0.0;  // node 0
    coords[3] = 1.0;
    coords[4] = 0.0;
    coords[5] = 0.0;  // node 1
    coords[6] = 1.0;
    coords[7] = 1.0;
    coords[8] = 0.0;  // node 2
    coords[9] = 0.0;
    coords[10] = 1.0;
    coords[11] = 0.0;  // node 3

    coords[12] = 0.0;
    coords[13] = 0.0;
    coords[14] = 1.0;  // node 4
    coords[15] = 1.0;
    coords[16] = 0.0;
    coords[17] = 1.0;  // node 5
    coords[18] = 1.0;
    coords[19] = 1.0;
    coords[20] = 1.0;  // node 6
    coords[21] = 0.0;
    coords[22] = 1.0;
    coords[23] = 1.0;  // node 7

    // edge nodes
    coords[24] = 0.5;
    coords[25] = 0.0;
    coords[26] = 0.0;  // node 8
    coords[27] = 1.0;
    coords[28] = 0.5;
    coords[29] = 0.0;  // node 9
    coords[30] = 0.5;
    coords[31] = 1.0;
    coords[32] = 0.0;  // node 10
    coords[33] = 0.0;
    coords[34] = 0.5;
    coords[35] = 0.0;  // node 11

    coords[36] = 0.5;
    coords[37] = 0.0;
    coords[38] = 1.0;  // node 12
    coords[39] = 1.0;
    coords[40] = 0.5;
    coords[41] = 1.0;  // node 13
    coords[42] = 0.5;
    coords[43] = 1.0;
    coords[44] = 1.0;  // node 14
    coords[45] = 0.0;
    coords[46] = 0.5;
    coords[47] = 1.0;  // node 15

    coords[48] = 0.0;
    coords[49] = 0.0;
    coords[50] = 0.5;  // node 16
    coords[51] = 1.0;
    coords[52] = 0.0;
    coords[53] = 0.5;  // node 17
    coords[54] = 1.0;
    coords[55] = 1.0;
    coords[56] = 0.5;  // node 18
    coords[57] = 0.0;
    coords[58] = 1.0;
    coords[59] = 0.5;  // node 19

    // face nodes
    coords[60] = 0.0;
    coords[61] = 0.5;
    coords[62] = 0.5;  // node 20
    coords[63] = 1.0;
    coords[64] = 0.5;
    coords[65] = 0.5;  // node 21
    coords[66] = 0.5;
    coords[67] = 0.0;
    coords[68] = 0.5;  // node 22
    coords[69] = 0.5;
    coords[70] = 1.0;
    coords[71] = 0.5;  // node 23
    coords[72] = 0.5;
    coords[73] = 0.5;
    coords[74] = 0.0;  // node 24
    coords[75] = 0.5;
    coords[76] = 0.5;
    coords[77] = 1.0;  // node 25

    // centroid
    coords[78] = 0.5;
    coords[79] = 0.5;
    coords[80] = 0.5;  // node 26
  }

  static void computeShape(const double* xr, double* phi)
  {
    SLIC_ASSERT(xr != nullptr);
    SLIC_ASSERT(phi != nullptr);

    const double r = xr[0];
    const double s = xr[1];
    const double t = xr[2];

    // bar element along r
    const double r1 = (r - 1.) * (2. * r - 1);
    const double r2 = 4. * r * (1. - r);
    const double r3 = r * (2. * r - 1.);

    // bar element along s
    const double s1 = (s - 1.) * (2. * s - 1);
    const double s2 = 4. * s * (1. - s);
    const double s3 = s * (2. * s - 1.);

    // bar element along t
    const double t1 = (t - 1) * (2. * t - 1);
    const double t2 = 4. * t * (1. - t);
    const double t3 = t * (2. * t - 1.);

    // corner nodes
    phi[0] = r1 * s1 * t1;
    phi[1] = r3 * s1 * t1;
    phi[2] = r3 * s3 * t1;
    phi[3] = r1 * s3 * t1;

    phi[4] = r1 * s1 * t3;
    phi[5] = r3 * s1 * t3;
    phi[6] = r3 * s3 * t3;
    phi[7] = r1 * s3 * t3;

    // edge nodes
    phi[8] = r2 * s1 * t1;
    phi[9] = r3 * s2 * t1;
    phi[10] = r2 * s3 * t1;
    phi[11] = r1 * s2 * t1;

    phi[12] = r2 * s1 * t3;
    phi[13] = r3 * s2 * t3;
    phi[14] = r2 * s3 * t3;
    phi[15] = r1 * s2 * t3;

    phi[16] = r1 * s1 * t2;
    phi[17] = r3 * s1 * t2;
    phi[18] = r3 * s3 * t2;
    phi[19] = r1 * s3 * t2;

    // face nodes
    phi[20] = r1 * s2 * t2;
    phi[21] = r3 * s2 * t2;
    phi[22] = r2 * s1 * t2;
    phi[23] = r2 * s3 * t2;
    phi[24] = r2 * s2 * t1;
    phi[25] = r2 * s2 * t3;

    // centroid
    phi[26] = r2 * s2 * t2;
  }

  static void computeDerivatives(const double* xr, double* phidot)
  {
    SLIC_ASSERT(xr != nullptr);
    SLIC_ASSERT(phidot != nullptr);

    const double r = xr[0];
    const double s = xr[1];
    const double t = xr[2];

    // bar element along r
    const double r1 = (r - 1.) * (2. * r - 1);
    const double r2 = 4. * r * (1. - r);
    const double r3 = r * (2. * r - 1.);

    // bar element along s
    const double s1 = (s - 1.) * (2. * s - 1);
    const double s2 = 4. * s * (1. - s);
    const double s3 = s * (2. * s - 1.);

    // bar element along t
    const double t1 = (t - 1) * (2. * t - 1);
    const double t2 = 4. * t * (1. - t);
    const double t3 = t * (2. * t - 1.);

    // bar element derivatives along r
    const double dr1 = 4. * r - 3.;
    const double dr2 = 4. - 8. * r;
    const double dr3 = 4. * r - 1.;

    // bar element derivatives along s
    const double ds1 = 4. * s - 3.;
    const double ds2 = 4. - 8. * s;
    const double ds3 = 4. * s - 1.;

    // bar element derivatives along t
    const double dt1 = 4. * t - 3.;
    const double dt2 = 4. - 8. * t;
    const double dt3 = 4. * t - 1.;

    // r derivatives
    phidot[0] = dr1 * s1 * t1;
    phidot[1] = dr3 * s1 * t1;
    phidot[2] = dr3 * s3 * t1;
    phidot[3] = dr1 * s3 * t1;
    phidot[4] = dr1 * s1 * t3;
    phidot[5] = dr3 * s1 * t3;
    phidot[6] = dr3 * s3 * t3;
    phidot[7] = dr1 * s3 * t3;
    phidot[8] = dr2 * s1 * t1;
    phidot[9] = dr3 * s2 * t1;
    phidot[10] = dr2 * s3 * t1;
    phidot[11] = dr1 * s2 * t1;
    phidot[12] = dr2 * s1 * t3;
    phidot[13] = dr3 * s2 * t3;
    phidot[14] = dr2 * s3 * t3;
    phidot[15] = dr1 * s2 * t3;
    phidot[16] = dr1 * s1 * t2;
    phidot[17] = dr3 * s1 * t2;
    phidot[18] = dr3 * s3 * t2;
    phidot[19] = dr1 * s3 * t2;
    phidot[20] = dr1 * s2 * t2;
    phidot[21] = dr3 * s2 * t2;
    phidot[22] = dr2 * s1 * t2;
    phidot[23] = dr2 * s3 * t2;
    phidot[24] = dr2 * s2 * t1;
    phidot[25] = dr2 * s2 * t3;
    phidot[26] = dr2 * s2 * t2;

    // s derivatives
    phidot[27] = r1 * ds1 * t1;
    phidot[28] = r3 * ds1 * t1;
    phidot[29] = r3 * ds3 * t1;
    phidot[30] = r1 * ds3 * t1;
    phidot[31] = r1 * ds1 * t3;
    phidot[32] = r3 * ds1 * t3;
    phidot[33] = r3 * ds3 * t3;
    phidot[34] = r1 * ds3 * t3;
    phidot[35] = r2 * ds1 * t1;
    phidot[36] = r3 * ds2 * t1;
    phidot[37] = r2 * ds3 * t1;
    phidot[38] = r1 * ds2 * t1;
    phidot[39] = r2 * ds1 * t3;
    phidot[40] = r3 * ds2 * t3;
    phidot[41] = r2 * ds3 * t3;
    phidot[42] = r1 * ds2 * t3;
    phidot[43] = r1 * ds1 * t2;
    phidot[44] = r3 * ds1 * t2;
    phidot[45] = r3 * ds3 * t2;
    phidot[46] = r1 * ds3 * t2;
    phidot[47] = r1 * ds2 * t2;
    phidot[48] = r3 * ds2 * t2;
    phidot[49] = r2 * ds1 * t2;
    phidot[50] = r2 * ds3 * t2;
    phidot[51] = r2 * ds2 * t1;
    phidot[52] = r2 * ds2 * t3;
    phidot[53] = r2 * ds2 * t2;

    // t derivatives
    phidot[54] = r1 * s1 * dt1;
    phidot[55] = r3 * s1 * dt1;
    phidot[56] = r3 * s3 * dt1;
    phidot[57] = r1 * s3 * dt1;
    phidot[58] = r1 * s1 * dt3;
    phidot[59] = r3 * s1 * dt3;
    phidot[60] = r3 * s3 * dt3;
    phidot[61] = r1 * s3 * dt3;
    phidot[62] = r2 * s1 * dt1;
    phidot[63] = r3 * s2 * dt1;
    phidot[64] = r2 * s3 * dt1;
    phidot[65] = r1 * s2 * dt1;
    phidot[66] = r2 * s1 * dt3;
    phidot[67] = r3 * s2 * dt3;
    phidot[68] = r2 * s3 * dt3;
    phidot[69] = r1 * s2 * dt3;
    phidot[70] = r1 * s1 * dt2;
    phidot[71] = r3 * s1 * dt2;
    phidot[72] = r3 * s3 * dt2;
    phidot[73] = r1 * s3 * dt2;
    phidot[74] = r1 * s2 * dt2;
    phidot[75] = r3 * s2 * dt2;
    phidot[76] = r2 * s1 * dt2;
    phidot[77] = r2 * s3 * dt2;
    phidot[78] = r2 * s2 * dt1;
    phidot[79] = r2 * s2 * dt3;
    phidot[80] = r2 * s2 * dt2;
  }
};

} /* namespace mint */
}  // namespace axom

#endif /* MINT_LAGRANGE_HEXA_27_HPP_ */
