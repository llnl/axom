// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_NUMERICS_QUADRATURE_HPP_
#define AXOM_NUMERICS_QUADRATURE_HPP_

#include "axom/core/Array.hpp"

/*!
 * \file quadrature.hpp
 * The functions declared in this header file find the nodes and weights of 
 *   arbitrary order quadrature rules
 */

namespace axom
{
namespace numerics
{

/*!
 * \class QuadratureRule
 *
 * \brief Stores a fixed array of 1D quadrature points and weights
 */
class QuadratureRule
{
  // Define friend functions so rules can only be created via compute_rule() methods
  friend QuadratureRule compute_gauss_legendre(int npts);

public:
  QuadratureRule() = default;

  //! \brief Accessor for quadrature nodes
  double node(size_t idx) const { return m_nodes[idx]; };

  //! \brief Accessor for quadrature weights
  double weight(size_t idx) const { return m_weights[idx]; };

private:
  //! \brief Private constructor for use in compute_<rule>() methods. Only allocates memory
  QuadratureRule(int npts) : m_nodes(npts), m_weights(npts) { };

  axom::Array<double> m_nodes;
  axom::Array<double> m_weights;
};

/*!
 * \brief Computes a 1D quadrature rule of Gauss-Legendre points 
 *
 * \param [in] npts The number of points in the rule
 * 
 * A Gauss-Legendre rule with \a npts points can exactly integrate
 *  polynomials of order 2 * npts - 1
 *
 * Algorithm adapted from the MFEM implementation in `mfem/fem/intrules.cpp`
 * 
 * \note This method constructs the points by scratch each time, without caching
 * \sa get_gauss_legendre(int)
 *
 * \return The `QuadratureRule` object which contains axom::Array<double>'s of nodes and weights
 */
QuadratureRule compute_gauss_legendre(int npts);

/*!
 * \brief Computes or accesses a precomputed 1D quadrature rule of Gauss-Legendre points 
 *
 * \param [in] npts The number of points in the rule
 * 
 * A Gauss-Legendre rule with \a npts points can exactly integrate
 *  polynomials of order 2 * npts - 1
 *
 * \note If this method has already been called for a given order, it will reuse the same quadrature points
 *  without needing to recompute them
 *
 * \return The `QuadratureRule` object which contains axom::Array<double>'s of nodes and weights
 */
const QuadratureRule& get_gauss_legendre(int npts);

} /* end namespace numerics */
} /* end namespace axom */

#endif  // AXOM_NUMERICS_QUADRATURE_HPP_
