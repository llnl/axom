// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_NUMERICS_INTERNAL_RATIONAL_QUADRATURE_HPP_
#define AXOM_NUMERICS_INTERNAL_RATIONAL_QUADRATURE_HPP_

#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/core/numerics/quadrature.hpp"

#include <complex>
#include <string>

namespace axom
{
namespace numerics
{
namespace internal
{

using Complex = std::complex<double>;

/*!
 * \note In the rational Fejer implementation, the suffix `_m11` means
 * "minus-one-to-one": the symmetric interval `[-1,1]`. The public Axom API
 * accepts and returns unit-interval (`[0,1]`) rules because downstream curve
 * parameterizations use that domain, but the Deckers / van Deun Cayley
 * transforms, rational Chebyshev phase equation, Blaschke basis, and moment
 * recurrences are naturally stated on `[-1,1]`. Internal helpers and diagnostics
 * keep `_m11` names until the final affine map to `[0,1]` so the implementation
 * stays close to the reference algorithms.
 */

/*!
 * \brief Developer diagnostics for the internal rational Fejer construction.
 *
 * These diagnostics expose the pole canonicalization, the intermediate
 * rational Chebyshev rule, and the final basis/weight assembly used by Axom's
 * `rfejer` implementation. They are intended for direct unit tests and for
 * one-off comparisons against reference implementations.
 */
struct RationalFejerDiagnostics
{
  struct Step
  {
    int step_index {0};
    int coeff_index_before {0};
    int component_count {0};
    int pole_multiplicity_so_far {0};
    std::complex<double> pole_m11 {0.0, 0.0};
    axom::Array<double> integral_coefficients_before;
    axom::Array<double> weighted_row0;
    axom::Array<double> weighted_row1;
    axom::Array<double> projected_row0;
    axom::Array<double> projected_row1;
    axom::Array<double> projected_row0_terms;
    axom::Array<double> projected_row1_terms;
    int projected_row_terms_node_count {0};
    double rhs0 {0.0};
    double rhs1 {0.0};
    double a00 {0.0};
    double a01 {0.0};
    double a10 {0.0};
    double a11 {0.0};
    bool used_pivot {false};
    double factor {0.0};
    double reduced_a11 {0.0};
    double reduced_rhs1 {0.0};
    double solved_x0 {0.0};
    double solved_x1 {0.0};
    axom::Array<double> orthogonal_integrals;
    axom::Array<double> integral_coefficients_after;
  };

  axom::Array<std::complex<double>> canonical_poles_m11;
  axom::Array<std::complex<double>> cayley_poles;
  axom::Array<double> rational_chebyshev_nodes_m11;
  axom::Array<double> rational_chebyshev_weights_m11;
  axom::Array<double> integral_coefficients;
  /// Sampled Fejer correction; final weights equal Chebyshev weights times this array.
  axom::Array<double> weight_correction_m11;
  axom::Array<double> final_weights_m11;
  axom::Array<double> basis_matrix_transpose_m11;
  int basis_matrix_row_count {0};
  int basis_matrix_col_count {0};
  axom::Array<double> nodes_01;
  axom::Array<double> weights_01;
  axom::Array<Step> steps;
};

/*!
 * \brief Compute the rational Chebyshev rule on `[-1,1]`.
 *
 * This is the van Deun et al. rational Chebyshev stage used by the rational Fejer
 * construction. It is exposed internally so tests can validate paper examples
 * and the intermediate node/weight rule before the final Fejer solve.
 *
 * \param [in] poles_m11 Pole sequence in the `[-1,1]` coordinate system.
 * \param [out] nodes Rational Chebyshev nodes on `[-1,1]`.
 * \param [out] weights Rational Chebyshev weights on `[-1,1]`.
 * \param [in] allocatorID Allocator used for the output arrays.
 *
 * \pre `poles_m11` is non-empty and all finite real poles lie outside `[-1,1]`.
 */
void compute_rational_chebyshev_data_m11(axom::ArrayView<const Complex> poles_m11,
                                         axom::Array<double>& nodes,
                                         axom::Array<double>& weights,
                                         int allocatorID = axom::getDefaultAllocatorID());

/*!
 * \brief Compute a fixed-order rational Fejer rule on `[-1,1]`.
 *
 * This helper runs the Deckers et al. rational Fejer construction in the
 * reference `[-1,1]` domain, including pole canonicalization, the intermediate
 * rational Chebyshev stage, and the final rational-basis weight solve.
 *
 * \param [in] poles_m11 Pole sequence in the `[-1,1]` coordinate system.
 * \param [out] nodes Rational Fejer nodes on `[-1,1]`.
 * \param [out] weights Rational Fejer weights on `[-1,1]`.
 * \param [in] allocatorID Allocator used for the output arrays.
 *
 * \pre `poles_m11` is non-empty and all finite real poles lie outside `[-1,1]`.
 */
void compute_rational_fejer_data_m11(axom::ArrayView<const Complex> poles_m11,
                                     axom::Array<double>& nodes,
                                     axom::Array<double>& weights,
                                     int allocatorID = axom::getDefaultAllocatorID());

/*!
 * \brief Populate diagnostics for a rational Fejer rule on `[-1,1]`.
 *
 * The diagnostic data records the same internal quantities used by
 * \c compute_rational_fejer_data_m11(), including canonical poles, Cayley
 * poles, rational Chebyshev nodes/weights, sampled basis columns, basis
 * integral coefficients, and per-step orthogonalization data.
 *
 * \param [in] poles_m11 Pole sequence in the `[-1,1]` coordinate system.
 * \param [out] diagnostics Diagnostic structure populated by the construction.
 * \param [in] allocatorID Allocator used for diagnostic arrays.
 *
 * \note Diagnostics are intended for host-side testing and debugging. The nested arrays
 * in \c RationalFejerDiagnostics should be allocated with a host-accessible allocator.
 *
 * \pre `poles_m11` is non-empty and all finite real poles lie outside `[-1,1]`.
 */
void compute_rational_fejer_diagnostics_m11(axom::ArrayView<const Complex> poles_m11,
                                            RationalFejerDiagnostics& diagnostics,
                                            int allocatorID = axom::getDefaultAllocatorID());

/*!
 * \brief Build the canonical rational Fejer cache key for `[-1,1]` poles.
 *
 * This internal helper is exposed for tests that verify cache-key canonicalization
 * follows the same `[-1,1]` pole pipeline as the rule construction.
 */
std::string make_rational_fejer_cache_key_m11(axom::ArrayView<const Complex> poles_m11,
                                              int allocatorID = axom::getDefaultAllocatorID());

/*!
 * \brief Populate diagnostics for a rational Fejer rule on `[0,1]`.
 *
 * This is the unit-interval diagnostic entry point corresponding to Axom's
 * public rational Fejer rule. It maps the input poles to `[-1,1]`, runs the
 * internal reference-domain construction, and stores both the `[-1,1]`
 * diagnostic arrays and the final `[0,1]` nodes/weights.
 *
 * \param [in] poles01 Pole sequence in the `[0,1]` coordinate system.
 * \param [out] diagnostics Diagnostic structure populated by the construction.
 * \param [in] allocatorID Allocator used for diagnostic arrays.
 *
 * \note Diagnostics are intended for host-side testing and debugging. The nested arrays
 * in \c RationalFejerDiagnostics should be allocated with a host-accessible allocator.
 *
 * \pre `poles01` is non-empty and all finite real poles lie outside `[0,1]`.
 */
void compute_rational_fejer_diagnostics(axom::ArrayView<const Complex> poles01,
                                        RationalFejerDiagnostics& diagnostics,
                                        int allocatorID = axom::getDefaultAllocatorID());

}  // namespace internal
}  // namespace numerics
}  // namespace axom

#endif  // AXOM_NUMERICS_INTERNAL_RATIONAL_QUADRATURE_HPP_
