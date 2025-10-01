// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_NUMERICS_GAUSS_LEGENDGRE_HPP_
#define AXOM_NUMERICS_GAUSS_LEGENDGRE_HPP_

#include "axom/core/Array.hpp"

/*!
 * \file gauss_legendre.hpp
 * The functions declared in this header file find the nodes and weights of 
 *   arbitrary order Gauss-Legendre quadrature rules
 */

namespace axom
{
namespace numerics
{

//! \brief Class to store quadature nodes and weights
class QuadratureRule
{
  // Define friend functions so rules can only be created via compute_rule() methods
  friend QuadratureRule compute_gauss_legendre(int npts);

public:
  QuadratureRule() = default;

  double node(size_t idx) const { return m_nodes[idx]; };
  double weight(size_t idx) const { return m_weights[idx]; };

private:
  QuadratureRule(int npts) : m_nodes(npts), m_weights(npts) { };

  axom::Array<double> m_nodes;
  axom::Array<double> m_weights;
};

QuadratureRule compute_gauss_legendre(int npts);

const QuadratureRule& get_gauss_legendre(int npts);

} /* end namespace numerics */
} /* end namespace axom */

#endif  // AXOM_NUMERICS_GAUSS_LEGENDGRE_HPP_
