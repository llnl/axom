// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_NUMERICS_INTERNAL_RATIONAL_QUADRATURE_HPP_
#define AXOM_NUMERICS_INTERNAL_RATIONAL_QUADRATURE_HPP_

#include "axom/core/numerics/quadrature.hpp"

#include <complex>
#include <vector>

namespace axom
{
namespace numerics
{
namespace internal
{

using Complex = std::complex<double>;

/*!
 * \brief Developer diagnostics for the internal rational Fejer construction.
 *
 * These diagnostics expose the pole canonicalization, the intermediate
 * rational-Chebyshev rule, and the final basis/weight assembly used by Axom's
 * `rfejer` implementation. They are intended for direct unit tests and for
 * one-off comparisons against reference implementations such as QuaHOG.
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
    axom::Array<double> basis_coefficients_before;
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
    axom::Array<double> basis_coefficients_after;
  };

  std::vector<std::complex<double>> canonical_poles_m11;
  std::vector<std::complex<double>> cayley_poles;
  axom::Array<double> rcheb_nodes_m11;
  axom::Array<double> rcheb_weights_m11;
  axom::Array<double> basis_coefficients;
  axom::Array<double> basis_expansion_m11;
  axom::Array<double> final_weights_m11;
  axom::Array<double> basis_matrix_transpose_m11;
  int basis_matrix_row_count {0};
  int basis_matrix_col_count {0};
  axom::Array<double> nodes_01;
  axom::Array<double> weights_01;
  std::vector<Step> steps;
};

/*!
 * \brief Direct `rcheb` and `rfejer` entry points on `[-1,1]`.
 *
 * These helpers expose the van Deun / Deckers construction in the same domain
 * used in the source papers so core numerics tests can validate the reference
 * examples without going through the public `[0,1]` wrappers.
 */
void compute_rcheb_data_m11(const std::vector<Complex>& poles_m11,
                            std::vector<double>& nodes,
                            std::vector<double>& weights);

void compute_rational_fejer_data_m11(const std::vector<Complex>& poles_m11,
                                     std::vector<double>& nodes,
                                     std::vector<double>& weights);

void compute_rational_fejer_diagnostics_m11(const std::vector<Complex>& poles_m11,
                                            RationalFejerDiagnostics& diagnostics,
                                            int allocatorID = axom::getDefaultAllocatorID());

/*!
 * \brief Compute rational Fejer diagnostics for a unit-interval pole sequence.
 */
void compute_rational_fejer_diagnostics(
  const std::vector<Complex>& poles01,
  RationalFejerDiagnostics& diagnostics,
  int allocatorID = axom::getDefaultAllocatorID());

}  // namespace internal
}  // namespace numerics
}  // namespace axom

#endif  // AXOM_NUMERICS_INTERNAL_RATIONAL_QUADRATURE_HPP_
