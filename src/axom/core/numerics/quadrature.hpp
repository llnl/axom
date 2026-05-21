// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_NUMERICS_QUADRATURE_HPP_
#define AXOM_NUMERICS_QUADRATURE_HPP_

#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/core/memory_management.hpp"

#include <complex>

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
  friend QuadratureRule get_rational_fejer(axom::ArrayView<const std::complex<double>>, int);

public:
  //! \brief Accessor for the full array of quadrature nodes
  AXOM_HOST_DEVICE
  axom::ArrayView<const double> nodes() const { return axom::ArrayView<const double>(m_nodes); }

  //! \brief Accessor for the full array of quadrature weights
  AXOM_HOST_DEVICE
  axom::ArrayView<const double> weights() const { return axom::ArrayView<const double>(m_weights); }

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

/*!
 * \brief Computes a 1D rational Fejer quadrature rule on [0, 1]
 *
 * \param [in] poles01 The pole sequence in the [0, 1] parameter domain.
 *   Finite poles must lie outside [0, 1]. Infinite poles may be supplied as
 *   `std::complex<double>(inf, 0.)`.
 * \param [out] nodes The array of 1D nodes on [0, 1]
 * \param [out] weights The array of weights on [0, 1]
 *   The returned arrays contain one more entry than the canonicalized pole
 *   sequence used internally. In particular, `m` real/infinite poles produce an
 *   `(m + 1)`-point rule. Non-real poles are canonicalized in conjugate pairs;
 *   if a caller supplies one side of a complex pair, the missing conjugate is
 *   added before the point count is determined.
 *
 * This rule implements the extended rational Fejer construction of
 * Deckers, Mougaida, and Belhadjsalah, "Algorithm 973: Extended Rational
 * Fejer Quadrature Rules Based on Chebyshev Orthogonal Rational Functions"
 * (ACM TOMS, 2017). Internally, Axom first builds the companion rational
 * Chebyshev rule on [-1, 1], following the fixed-pole rational
 * Gauss-Chebyshev / near-best interpolation construction of van Deun,
 * Deckers, Bultheel, and Weideman, "Algorithm 882" (ACM TOMS, 2008),
 * and then rescales the resulting nodes and weights to [0, 1].
 *
 * Fejer quadrature is closely related to Clenshaw-Curtis quadrature: both are
 * Chebyshev-based families of rules. A standard Clenshaw-Curtis rule uses
 * Chebyshev extrema and includes the interval endpoints, while a standard Fejer
 * rule uses interior Chebyshev nodes and excludes the endpoints. The present
 * routine goes one step further and builds a rational Fejer rule, so the weights
 * are adapted to a prescribed pole set rather than to the polynomial-only
 * Chebyshev setting. When all poles are at infinity, this construction reduces
 * to the standard interior Chebyshev / Fejer rule.
 *
 * \note This method constructs the points from scratch each time, without caching.
 * \pre `poles01` is non-empty and all finite real poles lie outside `[0,1]`.
 * \sa get_rational_fejer(axom::ArrayView<const std::complex<double>>, int)
 */
void compute_rational_fejer_data(axom::ArrayView<const std::complex<double>> poles01,
                                 axom::Array<double>& nodes,
                                 axom::Array<double>& weights,
                                 int allocatorID = axom::getDefaultAllocatorID());

/*!
 * \brief Computes or accesses a cached 1D rational Fejer quadrature rule on [0, 1]
 *
 * \param [in] poles01 The pole sequence in the [0, 1] parameter domain.
 *   The returned rule has one more point than the canonicalized pole sequence;
 *   see \c compute_rational_fejer_data() for the point-count convention.
 *
 * \note If this method has already been called for the same pole sequence and allocator,
 *  it will reuse the same quadrature data without recomputing it.
 * \note Although rational Fejer is Chebyshev-based like Clenshaw-Curtis, the
 *  cached rule here is keyed by the explicit pole sequence, since the pole data
 *  changes the exactness space and therefore the resulting weights.
 * \warning This process-wide cache is capped to avoid unbounded growth. Once
 *  the cap is reached, adding a new rule evicts the least recently used cached
 *  rule, which invalidates views returned by the evicted \c QuadratureRule.
 *  Prefer \c compute_rational_fejer_data() for one-off generated pole
 *  sequences or when the caller needs owned arrays.
 * \pre `poles01` is non-empty and all finite real poles lie outside `[0,1]`.
 *
 * \return The `QuadratureRule` object which contains axom::ArrayView<double>'s of stored nodes and weights
 */
QuadratureRule get_rational_fejer(axom::ArrayView<const std::complex<double>> poles01,
                                  int allocatorID = axom::getDefaultAllocatorID());

} /* end namespace numerics */
} /* end namespace axom */

#endif  // AXOM_NUMERICS_QUADRATURE_HPP_
