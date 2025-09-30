// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_NUMERICS_GAUSS_LEGENDGRE_HPP_
#define AXOM_NUMERICS_GAUSS_LEGENDGRE_HPP_

#include <math.h>

/*!
 * \file gauss_legendre.hpp
 * The functions declared in this header file find the nodes and weights of 
 *   arbitrary order Gauss-Legendre quadrature rules
 */

namespace axom
{
namespace numerics
{

void compute_rule(int np, double* nodes, double* weights);

} /* end namespace numerics */
} /* end namespace axom */

#endif  // AXOM_NUMERICS_GAUSS_LEGENDGRE_HPP_
