// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_NUMERICS_QUADRATURE_HPP_
#define AXOM_NUMERICS_QUADRATURE_HPP_

#include "axom/core/Array.hpp"
#include "axom/core/memory_management.hpp"

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
 * \brief Stores fixed views to arrays of 1D quadrature points and weights
 */
class QuadratureRule
{
  // Define friend functions so rules can only be created via get_rule() methods
  friend QuadratureRule get_gauss_legendre(int, int);

public:
  //! \brief Accessor for quadrature nodes
  AXOM_HOST_DEVICE
  double node(size_t idx) const { return m_nodes[idx]; };

  //! \brief Accessor for quadrature weights
  AXOM_HOST_DEVICE
  double weight(size_t idx) const { return m_weights[idx]; };

  //! \brief Accessor for the size of the quadrature rule
  AXOM_HOST_DEVICE
  int getNumPoints() const { return static_cast<int>(m_nodes.size()); }

private:
  //! \brief Use a private constructor to avoid creation of an invalid rule
  QuadratureRule(axom::ArrayView<double> nodes, axom::ArrayView<double> weights)
    : m_nodes(nodes)
    , m_weights(weights) { };

  axom::ArrayView<double> m_nodes;
  axom::ArrayView<double> m_weights;
};

/*!
 * \brief Computes a 1D quadrature rule of Gauss-Legendre points 
 *
 * \param [in] npts The number of points in the rule
 * \param [out] nodes The array of 1D nodes
 * \param [out] weights The array of weights
 * 
 * A Gauss-Legendre rule with \a npts points can exactly integrate
 *  polynomials of order 2 * npts - 1
 *
 * Algorithm adapted from the MFEM implementation in `mfem/fem/intrules.cpp`
 * 
 * \note This method constructs the points by scratch each time, without caching
 * \sa get_gauss_legendre(int)
 */
void compute_gauss_legendre_data(int npts,
                                 axom::Array<double>& nodes,
                                 axom::Array<double>& weights,
                                 int allocatorID = axom::getDefaultAllocatorID());

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
 * \return The `QuadratureRule` object which contains axom::ArrayView<double>'s of stored nodes and weights
 */
QuadratureRule get_gauss_legendre(int npts, int allocatorID = axom::getDefaultAllocatorID());

} /* end namespace numerics */
} /* end namespace axom */

#endif  // AXOM_NUMERICS_QUADRATURE_HPP_
