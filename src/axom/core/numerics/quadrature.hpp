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

#include <cassert>
#include <complex>
#include <utility>

/*!
 * \file quadrature.hpp
 *
 * \brief Numerical quadrature rules for integration
 *
 * This file provides 1D quadrature rules for numerical integration. 
 * It supports two main quadrature families:
 *
 * ## Polynomial Quadrature (Gauss-Legendre)
 *
 * Standard Gauss-Legendre rules integrate smooth polynomial-like functions:
 * - Use for general-purpose integration of smooth functions
 * - N-point rule exactly integrates polynomials of degree 2N-1
 * - See \c get_gauss_legendre() and \c compute_gauss_legendre_data()
 *
 * ## Rational Quadrature (Rational Fejer)
 *
 * Rational Fejer rules adapt to functions with known singularities at prescribed locations called **poles**:
 * - Use when integrands have singularities or sharp features at known locations
 * - Poles are complex values where the function (or its approximation) is singular
 * - Finite real poles must lie outside the integration domain [0, 1]
 * - Complex poles are automatically completed with their conjugate pair
 * - See \c get_rational_fejer() and \c compute_rational_fejer_data()
 *
 * ### When to Use Rational Fejer Quadrature
 *
 * Consider rational quadrature when:
 * - Integrating functions with known singularities outside [0, 1]
 * - Functions have sharp gradients near specific locations (geometric corners)
 * - Standard Gauss-Legendre requires very high orders for accuracy
 * - Integrand is well-approximated by rational functions
 *
 * ### Example: Integrating Near a Singularity
 *
 * \code{.cpp}
 * // Integrate f(x) = 1/(x - 1.5) from 0 to 1
 * // Singularity at x = 1.5 is outside [0, 1]
 *
 * std::vector<std::complex<double>> poles = {std::complex<double>(1.5, 0.0)};
 * auto rule = axom::numerics::get_rational_fejer(poles);
 *
 * double integral = 0.0;
 * for (int i = 0; i < rule.getNumPoints(); ++i) {
 *   integral += rule.weight(i) * (1.0 / (rule.node(i) - 1.5));
 * }
 * // Result: integral ≈ log(0.5) ≈ -0.693 (2 points vs 10+ for Gauss-Legendre)
 * \endcode
 *
 * ### Understanding Poles
 *
 * A **pole** is a location where a function becomes singular (infinite). For
 * rational Fejer quadrature on [0, 1]:
 *
 * - **Finite real poles** must be < 0 or > 1 (outside domain)
 *   - Example: pole at x = 1.5 for function f(x) = 1/(x - 1.5)
 *
 * - **Complex poles** come in conjugate pairs (auto-completed if needed)
 *   - Example: specifying z = 1.5 + 0.2i automatically includes z = 1.5 - 0.2i
 *
 * - **Infinite poles** represent the polynomial limit
 *   - Use std::complex<double>(inf, 0.0) to specify
 *   - All poles at infinity → standard Fejer (Chebyshev-based) rule
 *
 * The number of quadrature points equals (number of canonical poles + 1), where
 * canonical poles include auto-completed conjugates.
 *
 * ### QuadratureRule vs QuadratureRuleView
 *
 * Following Axom's Array/ArrayView pattern:
 * - \c QuadratureRule: Owns node and weight arrays (use for long-term storage)
 * - \c QuadratureRuleView: Non-owning view (lightweight, device-compatible)
 *
 * The \c get_* functions return views over cached data. Use \c .copy() to
 * create an owned copy if needed beyond immediate scope:
 *
 * \code{.cpp}
 * // Get cached view (lightweight)
 * auto view = axom::numerics::get_rational_fejer(poles);
 *
 * // Create owned copy if cache eviction is possible
 * axom::numerics::QuadratureRule owned = view.copy();
 * \endcode
 *
 * ### Cache Behavior
 *
 * Both \c get_gauss_legendre() and \c get_rational_fejer() use LRU caches:
 * - Gauss-Legendre: Unbounded cache (entries never evicted)
 * - Rational Fejer: 65,536 entry cache (LRU eviction when full)
 *
 * When rational Fejer cache is full, eviction invalidates views to evicted
 * rules. Use \c .copy() if you need stable storage beyond immediate use.
 *
 * ### Further Documentation
 *
 * See RATIONAL_QUADRATURE_GUIDE.md for:
 * - Detailed pole selection strategies
 * - Complete worked examples
 * - Performance considerations
 * - Comparison with Gauss-Legendre
 * - Device execution examples
 *
 * ### References
 *
 * - van Deun, Deckers, Bultheel, Weideman, "Algorithm 882: Rational Gauss-Chebyshev
 *   Quadrature", ACM Trans. Math. Softw., 2008 (rational Chebyshev stage)
 * - Deckers, Mougaida, Belhadjsalah, "Algorithm 973: Extended Rational Fejer
 *   Quadrature Rules", ACM Trans. Math. Softw., 2017 (full rational Fejer)
 */

namespace axom
{
namespace numerics
{

class QuadratureRuleView;

/*!
 * \class QuadratureRule
 *
 * \brief Owns arrays for a 1D quadrature rule.
 *
 * `QuadratureRule` stores its own node and weight arrays. Use \c view() when a
 * lightweight non-owning rule is needed, including in device kernels.
 */
class QuadratureRule
{
public:
  /*!
   * \brief Copy nodes and weights into owned arrays
   *
   * \pre `nodes.size() == weights.size()`.
   */
  QuadratureRule(axom::ArrayView<const double> nodes,
                 axom::ArrayView<const double> weights,
                 int allocatorID = axom::getDefaultAllocatorID())
    : m_nodes(nodes, allocatorID)
    , m_weights(weights, allocatorID)
  {
    assert(m_nodes.size() == m_weights.size());
  }

  //! \brief Copy a non-owning rule into owned arrays
  QuadratureRule(const QuadratureRuleView& rule, int allocatorID = axom::getDefaultAllocatorID());

  /*!
   * \brief Take ownership of node and weight arrays
   *
   * \pre `nodes.size() == weights.size()`.
   */
  QuadratureRule(axom::Array<double> nodes, axom::Array<double> weights)
    : m_nodes(std::move(nodes))
    , m_weights(std::move(weights))
  {
    assert(m_nodes.size() == m_weights.size());
  }

  //! \brief Accessor for the full array of quadrature nodes
  axom::ArrayView<const double> nodes() const { return m_nodes.view(); }

  //! \brief Accessor for the full array of quadrature weights
  axom::ArrayView<const double> weights() const { return m_weights.view(); }

  //! \brief Accessor for quadrature nodes
  double node(size_t idx) const { return m_nodes[idx]; };

  //! \brief Accessor for quadrature weights
  double weight(size_t idx) const { return m_weights[idx]; };

  //! \brief Accessor for the size of the quadrature rule
  int getNumPoints() const { return static_cast<int>(m_nodes.size()); }

  //! \brief Return a non-owning view rule backed by this owned rule
  QuadratureRuleView view() const;

private:
  axom::Array<double> m_nodes;
  axom::Array<double> m_weights;
};

/*!
 * \class QuadratureRuleView
 *
 * \brief Non-owning view over a 1D quadrature rule's points and weights.
 *
 * `QuadratureRuleView` is the lightweight, device-capturable representation of
 * a quadrature rule. It does not own the node or weight arrays.
 */
class QuadratureRuleView
{
public:
  AXOM_HOST_DEVICE QuadratureRuleView() = default;

  //! \brief Construct a view over existing node and weight arrays
  AXOM_HOST_DEVICE QuadratureRuleView(axom::ArrayView<const double> nodes,
                                      axom::ArrayView<const double> weights)
    : m_nodes(nodes)
    , m_weights(weights)
  {
    assert(m_nodes.size() == m_weights.size());
  }

  //! \brief Copy this non-owning rule into an owning rule
  QuadratureRule copy(int allocatorID = axom::getDefaultAllocatorID()) const
  {
    return QuadratureRule {*this, allocatorID};
  }

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
  axom::ArrayView<const double> m_nodes;
  axom::ArrayView<const double> m_weights;
};

inline QuadratureRule::QuadratureRule(const QuadratureRuleView& rule, int allocatorID)
  : QuadratureRule(rule.nodes(), rule.weights(), allocatorID)
{ }

inline QuadratureRuleView QuadratureRule::view() const
{
  return QuadratureRuleView {m_nodes.view(), m_weights.view()};
}

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
 * \return A non-owning \c QuadratureRuleView over cached nodes and weights
 */
QuadratureRuleView get_gauss_legendre(int npts, int allocatorID = axom::getDefaultAllocatorID());

/*!
 * \brief Computes a 1D rational Fejer quadrature rule on [0, 1]
 *
 * Rational Fejer quadrature adapts to functions with known singularities or sharp
 * features at prescribed **pole** locations. Unlike polynomial-based Gauss-Legendre,
 * this rule can efficiently integrate functions with near-singularities using far
 * fewer quadrature points.
 *
 * \param [in] poles01 The pole sequence in the [0, 1] parameter domain.
 *   - Finite real poles must lie outside [0, 1] (e.g., 1.5, -0.2, 2.0)
 *   - Complex poles are auto-completed to conjugate pairs
 *   - Infinite poles: use `std::complex<double>(inf, 0.0)`
 *   - The returned rule has `(num_canonical_poles + 1)` points
 *
 * \param [out] nodes The array of 1D nodes on [0, 1]
 * \param [out] weights The array of weights on [0, 1]
 * \param [in] allocatorID Allocator ID for output arrays (default: default allocator)
 *
 * ### Point Count
 *
 * The returned arrays contain one more entry than the canonicalized pole
 * sequence used internally. In particular, `m` real/infinite poles produce an
 * `(m + 1)`-point rule. Non-real poles are canonicalized in conjugate pairs;
 * if a caller supplies one side of a complex pair, the missing conjugate is
 * added before the point count is determined.
 *
 * ### Example: Near-Singularity Integration
 *
 * \code{.cpp}
 * // Integrate f(x) = 1/(x - 1.5) from 0 to 1
 * // Singularity at x = 1.5 is outside [0, 1]
 * axom::Array<std::complex<double>> poles = {std::complex<double>(1.5, 0.0)};
 *
 * axom::Array<double> nodes, weights;
 * axom::numerics::compute_rational_fejer_data(poles, nodes, weights);
 *
 * double integral = 0.0;
 * for (int i = 0; i < nodes.size(); ++i) {
 *   integral += weights[i] / (nodes[i] - 1.5);
 * }
 * // Result: 2 points achieve near machine-precision vs 10+ for Gauss-Legendre
 * \endcode
 *
 * ### Algorithm Details
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
 * In relation to rational Chebyshev quadrature, whose integrands are of the form
 * `f(x)/sqrt(1-x^2)`, rational Fejer modifies the weights to support unweighted
 * integrands of the form `f(x)`.
 * This construction reuses the rational Chebyshev nodes but modifies the weights.
 *
 * Fejer quadrature is closely related to Clenshaw-Curtis quadrature: both are
 * Chebyshev-based families of rules. A standard Clenshaw-Curtis rule uses
 * Chebyshev extrema and includes the interval endpoints, while a standard Fejer
 * rule uses interior Chebyshev nodes and excludes the endpoints. In the rational 
 * case, the weights are adapted to a prescribed pole set. When all poles are at infinity, 
 * this construction reduces to the standard interior Chebyshev / Fejer rule.
 *
 * \note This method constructs the points from scratch each time, without caching.
 *   For repeated use with the same poles, prefer \c get_rational_fejer() which caches computed rules.
 *
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
 * This is the recommended entry point for rational Fejer quadrature when the same
 * pole sequence will be reused. The function maintains a process-wide LRU cache 
 * to avoid recomputation while limiting overall memory consumption.
 *
 * \param [in] poles01 The pole sequence in the [0, 1] parameter domain.
 *   - Finite real poles must lie outside [0, 1]
 *   - Complex poles are auto-completed to conjugate pairs
 *   - Infinite poles: use `std::complex<double>(inf, 0.0)`
 *   - The returned rule has `(num_canonical_poles + 1)` points
 *
 * \param [in] allocatorID Allocator ID for cached arrays (default: default allocator)
 *   - Use unified/device allocators for GPU execution
 *   - Cache key includes allocator ID, so different allocators maintain separate entries
 *
 * ### Caching Behavior
 *
 * - **Cache hit:** Returns a view over previously computed rule (fast O(1) lookup)
 * - **Cache miss:** Computes rule, stores in cache, returns view (O(m^2) construction)
 * - **Cache full:** Evicts least-recently-used entry before inserting new rule
 *
 * ### When to Use .copy()
 *
 * The returned `QuadratureRuleView` points to cached storage. Use `.copy()` if:
 * - Storing the rule for use beyond immediate scope
 * - Uncertain about cache lifetime (e.g., generating many different pole sequences)
 * - Need guaranteed stable pointers to node/weight arrays
 *
 * \code{.cpp}
 * // Safe: immediate use within local scope
 * auto rule = axom::numerics::get_rational_fejer(poles);
 * double result = integrate(rule, f);
 *
 * // Safe: owned copy protects against cache eviction
 * axom::numerics::QuadratureRule owned =
 *   axom::numerics::get_rational_fejer(poles).copy();
 * // ... use owned later, even if cache evicts this entry
 * \endcode
 *
 * \note If this method has already been called for the same pole sequence and allocator,
 *  it will reuse the same quadrature data without recomputing it.
 * \note Cache keys are based on canonical [-1,1] poles, so equivalent pole sequences
 *  (e.g., with/without explicit conjugates) map to the same cached entry.
 *
 * \pre `poles01` is non-empty and all finite real poles lie outside `[0,1]`.
 * \return A non-owning \c QuadratureRuleView over cached nodes and weights
 * \sa compute_rational_fejer_data() for uncached construction
 */
QuadratureRuleView get_rational_fejer(axom::ArrayView<const std::complex<double>> poles01,
                                      int allocatorID = axom::getDefaultAllocatorID());

} /* end namespace numerics */
} /* end namespace axom */

#endif  // AXOM_NUMERICS_QUADRATURE_HPP_
