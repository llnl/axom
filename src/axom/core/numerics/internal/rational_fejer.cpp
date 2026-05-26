// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/LRUCache.hpp"
#include "axom/core/numerics/internal/rational_fejer_diagnostics.hpp"
#include "axom/core/numerics/internal/rational_quadrature_common.hpp"

#include <cmath>
#include <mutex>
#include <string>
#include <utility>

namespace axom
{
namespace numerics
{
namespace internal
{

/*!
 * \file rational_fejer.cpp
 *
 * \brief Implements fixed-order extended rational Fejer quadrature.
 *
 * This file implements the rational Fejer rule construction from Deckers,
 * Mougaida, and Belhadjsalah, "Algorithm 973: Extended Rational Fejer
 * Quadrature Rules Based on Chebyshev Orthogonal Rational Functions",
 * ACM Transactions on Mathematical Software 43(4), 2017.
 *
 * Given a prescribed pole sequence, the algorithm canonicalizes the poles,
 * completes missing complex conjugates, and builds the companion rational
 * Chebyshev rule. Rational Fejer keeps those Chebyshev nodes and solves for a
 * correction to the Chebyshev weights so that the final rule integrates
 * unweighted rational basis functions. Public entry points map between Axom's
 * `[0,1]` rule domain and the `[-1,1]` reference domain used by the paper
 * formulas.
 */

/// \brief Computes the asymptotic pole-moment integral used for distant poles.
inline Complex computeAsymptoticPoleMomentIntegral(const Complex& pole, int order, int parity, double tol)
{
  // For poles far from [-1,1], the closed-form logarithmic recurrence loses
  // accuracy. Deckers switches to this inverse-pole expansion and truncates
  // when the omitted tail is below the requested tolerance.
  int j = 1;
  while(std::log(std::pow(std::abs(pole), 2 * j) * (2 * j + order + 3 - parity) *
                 (2 * j + order + 1 - parity) * tol) < 0.0)
  {
    ++j;
  }

  Complex series_sum {0.0, 0.0};
  for(int exponent = 0; exponent <= 2 * j; exponent += 2)
  {
    const double d1 = exponent + order + 3 - parity;
    const double d2 = exponent + order + 1 - parity;
    series_sum += Complex {1.0, 0.0} / (d1 * d2 * std::pow(pole, exponent));
  }

  return 2.0 * (1 - parity) / (order + 1.0) +
    std::pow(-1.0, parity) * 4.0 * order * series_sum / std::pow(pole, 2 - parity);
}

/// \brief Computes exact pole-moment integrals used to seed rational Fejer coefficients.
inline axom::Array<double> computePoleMomentIntegrals(const Pole& pole,
                                                      int multiplicity,
                                                      int component_count,
                                                      double tol)
{
  // These are the non-orthogonal moments for one repeated pole. The Deckers
  // construction uses them as exact seed data before the basis coefficients
  // are orthogonalized against the sampled rational basis.
  axom::Array<Complex> J(multiplicity, multiplicity);
  J.fill(Complex {0.0, 0.0});

  if(pole.isInfinite())
  {
    for(int idx = 1; idx < multiplicity; idx += 2)
    {
      J[idx] = 2.0 / (idx + 2.0);
    }
  }
  else
  {
    const Complex pole_value = pole.value();
    const double c = std::abs(pole_value) - 1.0;
    // Deckers' moment recurrence is numerically split: near the interval,
    // integrate forward from the closed-form logarithmic expression; farther
    // away, use the asymptotic expansion and recurse backward.
    if((c < 1e10 * tol) || (std::log(10.0) * (4 * multiplicity + 50) + 90.0 * std::log(c) < 0.0) ||
       ((multiplicity < 3) && (c < 9.0)))
    {
      J[0] = (pole_value * pole_value - Complex {1.0, 0.0}) *
          std::log((pole_value + Complex {1.0, 0.0}) / (pole_value - Complex {1.0, 0.0})) -
        2.0 * pole_value;

      if(multiplicity > 1)
      {
        Complex F {2.0, 0.0};
        for(int k = 2; k <= multiplicity; ++k)
        {
          const double parity = (k - 1) % 2 == 0 ? 1.0 : -1.0;
          const Complex next =
            (pole_value * pole_value - Complex {1.0, 0.0}) * (1.0 - parity) / (k - 1.0) -
            2.0 * pole_value * J[k - 2] - pole_value * pole_value * F;
          F = J[k - 2];
          J[k - 1] = next;
        }
      }
    }
    else
    {
      const int parity = multiplicity % 2;
      J[multiplicity - 1] =
        computeAsymptoticPoleMomentIntegral(pole_value, multiplicity, parity, tol);
      if(multiplicity > 1)
      {
        Complex F = J[multiplicity - 1];
        J[multiplicity - 2] =
          computeAsymptoticPoleMomentIntegral(pole_value, multiplicity - 1, 1 - parity, tol);
        for(int k = multiplicity - 2; k >= 1; --k)
        {
          const double parity_term = (k + 1) % 2 == 0 ? 1.0 : -1.0;
          const Complex next =
            -((Complex {1.0, 0.0} / (pole_value * pole_value)) - Complex {1.0, 0.0}) *
              (1.0 - parity_term) / (k + 1.0) -
            2.0 * (Complex {1.0, 0.0} / pole_value) * J[k] -
            (Complex {1.0, 0.0} / (pole_value * pole_value)) * F;
          F = J[k];
          J[k - 1] = next;
        }
      }
    }
  }

  if(component_count == 1)
  {
    axom::Array<double> result(multiplicity, multiplicity);
    for(int i = 0; i < multiplicity; ++i)
    {
      result[i] = J[i].real();
    }
    return result;
  }

  axom::Array<double> result(2 * multiplicity, 2 * multiplicity);
  for(int i = 0; i < multiplicity; ++i)
  {
    result[2 * i] = J[i].real();
    result[2 * i + 1] = -J[i].imag();
  }
  return result;
}

/// \brief Stores the rational Fejer basis sampled at the companion Chebyshev nodes.
///
/// The companion rational Chebyshev rule supplies the nodes and seed weights.
/// The Fejer-specific solve below only computes the weight correction sampled in this basis.
class RationalFejerBasis
{
public:
  RationalFejerBasis(axom::ArrayView<const Complex> cayley_poles,
                     QuadratureRule companion_chebyshev_rule)
    : m_companionChebyshevRule(std::move(companion_chebyshev_rule))
    , m_QColumns(buildBasisColumns(cayley_poles, m_companionChebyshevRule.nodes()))
    , m_basisColumnCount(m_companionChebyshevRule.getNumPoints())
    , m_nodeCount(m_companionChebyshevRule.getNumPoints())
  { }

  axom::ArrayView<const double> nodes() const { return m_companionChebyshevRule.nodes(); }

  axom::ArrayView<const double> chebyshevWeights() const
  {
    return m_companionChebyshevRule.weights();
  }

  axom::ArrayView<const double> basisColumns() const { return m_QColumns.view(); }

  int basisColumnCount() const { return m_basisColumnCount; }

  int nodeCount() const { return m_nodeCount; }

  double node(int idx) const { return m_companionChebyshevRule.node(idx); }

  double basisValue(int column, int node) const { return m_QColumns[column * m_nodeCount + node]; }

  axom::Array<double> solveOrthogonalIntegrals(const Pole& pole,
                                               int multiplicity,
                                               int num_basis_columns,
                                               axom::ArrayView<const double> known_integrals,
                                               int component_count,
                                               const RationalFejerStepDiagnosticsRecorder& diagnostics =
                                                 RationalFejerStepDiagnosticsRecorder {}) const
  {
    assert(component_count == 1 || component_count == 2);
    assert(num_basis_columns >= component_count);
    assert(known_integrals.size() >= num_basis_columns);

    const int n = nodeCount();
    const bool has_imaginary_component = component_count == 2;
    const bool diagnostics_enabled = diagnostics.enabled();
    const int first_unknown_column = num_basis_columns - component_count;
    axom::Array<long double> weighted_row0(n, n);
    weighted_row0.fill(0.0L);
    axom::Array<long double> weighted_row1(has_imaginary_component ? n : 0,
                                           has_imaginary_component ? n : 0);
    weighted_row1.fill(0.0L);

    // A real pole contributes one basis component. A non-real pole contributes
    // the real and negative-imaginary parts of the same rational factor, which
    // keeps the eventual quadrature weights real.
    for(int i = 0; i < n; ++i)
    {
      const Complex basis_value = evaluatePoleBasisValue(pole, multiplicity, node(i));

      weighted_row0[i] = static_cast<long double>(chebyshevWeights()[i]) *
        static_cast<long double>(basis_value.real());
      if(has_imaginary_component)
      {
        weighted_row1[i] = static_cast<long double>(chebyshevWeights()[i]) *
          static_cast<long double>(-basis_value.imag());
      }
    }

    diagnostics.recordWeightedRows(weighted_row0, weighted_row1, has_imaginary_component);

    axom::Array<long double> bmat_row0(num_basis_columns, num_basis_columns);
    bmat_row0.fill(0.0L);
    axom::Array<long double> bmat_row1(has_imaginary_component ? num_basis_columns : 0,
                                       has_imaginary_component ? num_basis_columns : 0);
    bmat_row1.fill(0.0L);
    axom::Array<long double> bmat_row0_terms;
    axom::Array<long double> bmat_row1_terms;
    if(diagnostics_enabled)
    {
      bmat_row0_terms = axom::Array<long double>(num_basis_columns * n, num_basis_columns * n);
      if(has_imaginary_component)
      {
        bmat_row1_terms = axom::Array<long double>(num_basis_columns * n, num_basis_columns * n);
      }
    }
    // Project the pole-dependent function against the sampled orthonormal basis.
    // The trailing 1x1 or 2x2 block is then solved for the next basis coefficients.
    for(int col = 0; col < num_basis_columns; ++col)
    {
      long double value0 = 0.0L;
      long double value1 = 0.0L;
      for(int i = 0; i < n; ++i)
      {
        const long double basis_column_value = static_cast<long double>(basisValue(col, i));
        const long double term0 = weighted_row0[i] * basis_column_value;
        value0 += term0;
        if(diagnostics_enabled)
        {
          bmat_row0_terms[col * n + i] = term0;
        }
        if(has_imaginary_component)
        {
          const long double term1 = weighted_row1[i] * basis_column_value;
          value1 += term1;
          if(diagnostics_enabled)
          {
            bmat_row1_terms[col * n + i] = term1;
          }
        }
      }
      bmat_row0[col] = value0;
      if(has_imaginary_component)
      {
        bmat_row1[col] = value1;
      }
    }

    long double rhs0 = static_cast<long double>(known_integrals[first_unknown_column]);
    // Previously solved coefficients already contribute to this projected
    // moment equation; move those known terms to the right-hand side so the
    // trailing 1x1 or 2x2 block can be solved directly.
    for(int col = 0; col < first_unknown_column; ++col)
    {
      rhs0 -= bmat_row0[col] * static_cast<long double>(known_integrals[col]);
    }

    diagnostics.recordProjectedRow0(bmat_row0, bmat_row0_terms, n, rhs0);

    if(component_count == 1)
    {
      return axom::Array<double> {static_cast<double>(rhs0 / bmat_row0[first_unknown_column])};
    }

    long double rhs1 = static_cast<long double>(known_integrals[num_basis_columns - 1]);
    for(int col = 0; col < first_unknown_column; ++col)
    {
      rhs1 -= bmat_row1[col] * static_cast<long double>(known_integrals[col]);
    }

    long double a00 = bmat_row0[first_unknown_column];
    long double a01 = bmat_row0[first_unknown_column + 1];
    long double a10 = bmat_row1[first_unknown_column];
    long double a11 = bmat_row1[first_unknown_column + 1];

    diagnostics.recordProjectedRow1(bmat_row1, bmat_row1_terms, rhs1, a00, a01, a10, a11);

    // The complex-pole case is a tiny 2x2 solve. Pivoting is cheap and prevents
    // the diagnostics path from depending on which projected row happens to
    // carry the larger leading coefficient.
    const bool used_pivot = std::abs(a10) > std::abs(a00);
    if(used_pivot)
    {
      std::swap(a00, a10);
      std::swap(a01, a11);
      std::swap(rhs0, rhs1);
    }

    const long double factor = a10 / a00;
    const long double reduced_a11 = a11 - factor * a01;
    const long double reduced_rhs1 = rhs1 - factor * rhs0;
    const long double x1 = reduced_rhs1 / reduced_a11;
    const long double x0 = (rhs0 - a01 * x1) / a00;
    diagnostics.recordTinySolve(used_pivot, factor, reduced_a11, reduced_rhs1, x0, x1);
    return axom::Array<double> {static_cast<double>(x0), static_cast<double>(x1)};
  }

  axom::Array<double> evaluateExpansionAtChebyshevNodes(axom::ArrayView<const double> coefficients) const
  {
    axom::Array<double> values(nodeCount(), nodeCount());
    values.fill(0.0);
    for(int j = 0; j < coefficients.size(); ++j)
    {
      const double coefficient = coefficients[j];
      if(coefficient == 0.0)
      {
        continue;
      }

      for(int i = 0; i < nodeCount(); ++i)
      {
        values[i] += basisValue(j, i) * coefficient;
      }
    }

    return values;
  }

private:
  static Complex evaluatePoleBasisValue(const Pole& pole, int multiplicity, double node)
  {
    const Complex node_value {node, 0.0};
    if(pole.isInfinite())
    {
      return powInteger(node_value, multiplicity);
    }

    const Complex pole_value = pole.value();
    return powInteger((Complex {1.0, 0.0} - pole_value * node) / (node_value - pole_value),
                      multiplicity);
  }

  static axom::Array<double> buildBasisColumns(axom::ArrayView<const Complex> cayley_poles,
                                               axom::ArrayView<const double> nodes)
  {
    // Sample the orthonormal rational basis at every rational Chebyshev node.
    // Later stages reuse this matrix both when enforcing orthogonality and when
    // rebuilding the final quadrature weights from the recovered coefficients.
    // Layout is column-major by basis function: Q_columns[col * n + node].
    const int n = static_cast<int>(nodes.size());
    axom::Array<double> Q_columns(n * n, n * n);
    Q_columns.fill(0.0);

    for(int i = 0; i < n; ++i)
    {
      Q_columns[i] = 1.0 / std::sqrt(M_PI);
    }

    if(n == 1)
    {
      return Q_columns;
    }

    axom::Array<Complex> z_values(n, n);
    for(int i = 0; i < n; ++i)
    {
      // Nodes x in [-1,1] are lifted to the upper unit circle
      // z = x + i sqrt(1-x^2), where the Blaschke-product formula is evaluated.
      z_values[i] =
        Complex {nodes[i], std::sqrt(axom::utilities::max(0.0, 1.0 - nodes[i] * nodes[i]))};
    }

    axom::Array<Complex> blaschke_products(n, n);
    blaschke_products.fill(Complex {1.0, 0.0});
    axom::Array<Complex> conjugate_blaschke_products(n, n);
    conjugate_blaschke_products.fill(Complex {1.0, 0.0});
    const double normalizer = std::sqrt(1.0 / (2.0 * M_PI));

    for(int i = 0; i < n; ++i)
    {
      const Complex z = z_values[i];
      const double basis_normalizer =
        normalizer * std::sqrt(axom::utilities::max(0.0, 1.0 - std::norm(cayley_poles[0])));
      const Complex value = basis_normalizer *
        (z * conjugate_blaschke_products[i] / (Complex {1.0, 0.0} - cayley_poles[0] * z) +
         Complex {1.0, 0.0} / ((z - cayley_poles[0]) * blaschke_products[i]));
      Q_columns[n + i] = value.real();
    }

    for(int k = 1; k < n - 1; ++k)
    {
      const Complex pole = cayley_poles[k - 1];
      const Complex conjugate_pole = std::conj(pole);
      const double basis_normalizer =
        normalizer * std::sqrt(axom::utilities::max(0.0, 1.0 - std::norm(cayley_poles[k])));
      for(int i = 0; i < n; ++i)
      {
        const Complex z = z_values[i];
        blaschke_products[i] *= (z - pole) / (Complex {1.0, 0.0} - conjugate_pole * z);
        conjugate_blaschke_products[i] *= (z - conjugate_pole) / (Complex {1.0, 0.0} - pole * z);

        const Complex value = basis_normalizer *
          (z * conjugate_blaschke_products[i] / (Complex {1.0, 0.0} - cayley_poles[k] * z) +
           Complex {1.0, 0.0} / ((z - cayley_poles[k]) * blaschke_products[i]));
        Q_columns[(k + 1) * n + i] = value.real();
      }
    }

    return Q_columns;
  }

  QuadratureRule m_companionChebyshevRule;
  axom::Array<double> m_QColumns;
  int m_basisColumnCount {0};
  int m_nodeCount {0};
};

/// \brief Maps canonical poles to the signed Cayley-domain sequence used by the Fejer basis.
inline axom::Array<Complex> makeRationalFejerCayleyPoles(axom::ArrayView<const Pole> canonical_poles)
{
  axom::Array<Complex> cayley_poles(canonical_poles.size(), canonical_poles.size());
  for(axom::IndexType i = 0; i < canonical_poles.size(); ++i)
  {
    // The rational basis is parameterized in the Cayley-transformed pole domain,
    // while the exact pole integrals are organized using the original pole ordering.
    cayley_poles[i] = canonical_poles[i].cayleyTransformPreservingImaginarySign();
  }
  return cayley_poles;
}

/// \brief Prepared pole data shared by Fejer cache lookup, diagnostics, and construction.
///
/// Public callers supply poles in `[0,1]`, while the Deckers construction and
/// cache key are defined from the canonicalized `[-1,1]` pole sequence. Keeping
/// those derived quantities together prevents the cached and uncached paths
/// from drifting apart.
struct RationalFejerPoleData
{
  PoleSequence poles_m11;
  PoleSequence canonical_poles_m11;
  axom::Array<Complex> cayley_poles;
  double pole_tolerance {rationalFejerPoleTolerance()};
  double infinite_pole_threshold {rationalFejerInfinitePoleThreshold()};

  static RationalFejerPoleData from01(axom::ArrayView<const Complex> poles01)
  {
    const PoleSequence poles01_sequence = PoleSequence::fromComplex(poles01);
    poles01_sequence.validate(0.0, 1.0, "[0,1]");
    return fromValidatedM11(poles01_sequence.mapped01ToM11());
  }

  static RationalFejerPoleData fromM11(axom::ArrayView<const Complex> poles_m11)
  {
    return fromM11(PoleSequence::fromComplex(poles_m11));
  }

  static RationalFejerPoleData fromM11(PoleSequence poles_m11)
  {
    poles_m11.validate(-1.0, 1.0, "[-1,1]");
    return fromValidatedM11(std::move(poles_m11));
  }

  std::string cacheKey(int allocatorID) const
  {
    axom::fmt::memory_buffer key;
    axom::fmt::format_to(std::back_inserter(key), "{}|", allocatorID);
    for(const auto& pole : canonical_poles_m11.array())
    {
      axom::fmt::format_to(std::back_inserter(key), "{:a},{:a};", pole.real(), pole.imag());
    }
    return axom::fmt::to_string(key);
  }

private:
  static RationalFejerPoleData fromValidatedM11(PoleSequence poles_m11)
  {
    RationalFejerPoleData data;
    data.poles_m11 = std::move(poles_m11);
    data.canonical_poles_m11 =
      data.poles_m11.canonicalizedAndNormalized(data.pole_tolerance, data.infinite_pole_threshold);
    data.cayley_poles = makeRationalFejerCayleyPoles(data.canonical_poles_m11.view());
    return data;
  }
};

/// \brief Initializes expansion coefficients for the unweighted integration functional.
inline axom::Array<double> initializeUnweightedIntegralCoefficients(int num_poles)
{
  axom::Array<double> integral_coefficients(num_poles + 1, num_poles + 1);
  integral_coefficients.fill(0.0);
  integral_coefficients[0] = 2.0 / std::sqrt(M_PI);
  return integral_coefficients;
}

/// \brief Seeds exact non-orthogonal moments for the next newly encountered pole.
inline int seedNewPoleMomentCoefficients(const PoleSequence& canonical_poles,
                                         const Pole& pole,
                                         int coeff_index,
                                         int component_count,
                                         double pole_tolerance,
                                         axom::Array<double>& integral_coefficients)
{
  const int num_poles = static_cast<int>(canonical_poles.size());
  const int pole_multiplicity_so_far =
    canonical_poles.countMatching(pole, 0, coeff_index, pole_tolerance);

  if(pole_multiplicity_so_far != 1)
  {
    return pole_multiplicity_so_far;
  }

  // Seed all occurrences of a newly encountered pole with the exact
  // non-orthogonal pole moments from the Deckers construction.
  axom::Array<int> repeated_pole_offsets =
    canonical_poles.matchingOffsets(pole, coeff_index - 1, num_poles, pole_tolerance);
  const int multiplicity = static_cast<int>(repeated_pole_offsets.size());

  if(component_count == 2)
  {
    // Complex poles occupy adjacent coefficient slots in real/imaginary
    // order. Interleaving the matching conjugate offsets maps the exact
    // complex moments onto those real-valued slots.
    axom::Array<int> conjugate_offsets =
      canonical_poles.matchingOffsets(pole.conjugate(), coeff_index - 1, num_poles, pole_tolerance);
    repeated_pole_offsets = interleaveOffsets(repeated_pole_offsets, conjugate_offsets);
  }

  const auto pole_integrals =
    computePoleMomentIntegrals(pole, multiplicity, component_count, pole_tolerance);
  for(axom::IndexType i = 0; i < repeated_pole_offsets.size(); ++i)
  {
    integral_coefficients[coeff_index + repeated_pole_offsets[i] - 1] = pole_integrals[i];
  }

  return pole_multiplicity_so_far;
}

/// \brief Solves the unweighted integration functional in the rational Fejer basis.
inline axom::Array<double> solveUnweightedIntegralCoefficients(
  const PoleSequence& canonical_poles,
  const RationalFejerBasis& basis,
  double pole_tolerance,
  RationalFejerDiagnosticsRecorder& diagnostics_recorder)
{
  const int num_poles = static_cast<int>(canonical_poles.size());
  axom::Array<double> integral_coefficients = initializeUnweightedIntegralCoefficients(num_poles);

  int coeff_index = 1;
  while(coeff_index < num_poles + 1)
  {
    const Pole pole = canonical_poles[coeff_index - 1];
    const int component_count = pole.componentCount(pole_tolerance);
    const int pole_multiplicity_so_far = seedNewPoleMomentCoefficients(canonical_poles,
                                                                       pole,
                                                                       coeff_index,
                                                                       component_count,
                                                                       pole_tolerance,
                                                                       integral_coefficients);

    const int num_known_columns = coeff_index + component_count;
    auto known_integrals = integral_coefficients.view().subspan(0, num_known_columns);

    RationalFejerDiagnostics::Step step_diagnostics =
      diagnostics_recorder.beginStep(coeff_index,
                                     component_count,
                                     pole_multiplicity_so_far,
                                     pole,
                                     known_integrals);

    // Enforce orthogonality against the previously determined columns to recover
    // the next one or two rational basis coefficients.
    const RationalFejerStepDiagnosticsRecorder step_recorder =
      diagnostics_recorder.stepRecorder(step_diagnostics);
    const auto orthogonal_integrals = basis.solveOrthogonalIntegrals(pole,
                                                                     pole_multiplicity_so_far,
                                                                     num_known_columns,
                                                                     known_integrals,
                                                                     component_count,
                                                                     step_recorder);
    for(int i = 0; i < component_count; ++i)
    {
      integral_coefficients[coeff_index + i] = orthogonal_integrals[i];
    }

    diagnostics_recorder.finishStep(std::move(step_diagnostics), orthogonal_integrals, known_integrals);

    coeff_index += component_count;
  }

  return integral_coefficients;
}

/// \brief Applies the Fejer weight correction to the companion Chebyshev rule.
inline QuadratureRule applyFejerWeightCorrection(const RationalFejerBasis& basis,
                                                 axom::ArrayView<const double> weight_correction)
{
  axom::Array<double> weights(weight_correction.size(), weight_correction.size());
  const auto chebyshev_weights = basis.chebyshevWeights();
  // Rational Fejer uses the rational Chebyshev nodes. Only the weights change:
  // the sampled unweighted integration functional is a correction factor.
  for(int i = 0; i < weights.size(); ++i)
  {
    weights[i] = chebyshev_weights[i] * weight_correction[i];
  }
  axom::Array<double> nodes(basis.nodes(), axom::getDefaultAllocatorID());
  return QuadratureRule {std::move(nodes), std::move(weights)};
}

/// \brief Computes the companion rational Chebyshev rule that supplies Fejer nodes.
inline QuadratureRule computeCompanionRationalChebyshevRule(const RationalFejerPoleData& pole_data)
{
  // A rational Fejer rule for m finite/infinite poles uses m+1 rational Chebyshev nodes. 
  // The appended infinity pole supplies that extra Chebyshev-like degree.
  const PoleSequence rational_chebyshev_poles = pole_data.canonical_poles_m11.withAppendedInfinity();

  axom::Array<double> nodes;
  axom::Array<double> weights;
  compute_rational_chebyshev_data(rational_chebyshev_poles.view(), nodes, weights);
  return QuadratureRule {std::move(nodes), std::move(weights)};
}

/// \brief Computes a rational Fejer rule on [-1,1] with optional diagnostics.
inline QuadratureRule compute_rational_fejer_rule_m11_impl(
  const RationalFejerPoleData& pole_data,
  RationalFejerDiagnostics* diagnostics = nullptr,
  int allocatorID = axom::getDefaultAllocatorID())
{
  RationalFejerDiagnosticsRecorder diagnostics_recorder {diagnostics,
                                                         allocatorID,
                                                         pole_data.pole_tolerance};
  diagnostics_recorder.recordPoleSetup(pole_data.canonical_poles_m11.view(), pole_data.cayley_poles);

  QuadratureRule companion_chebyshev_rule = computeCompanionRationalChebyshevRule(pole_data);
  diagnostics_recorder.recordCompanionChebyshevRule(companion_chebyshev_rule.view());

  RationalFejerBasis basis(pole_data.cayley_poles, std::move(companion_chebyshev_rule));

  const axom::Array<double> integral_coefficients =
    solveUnweightedIntegralCoefficients(pole_data.canonical_poles_m11,
                                        basis,
                                        pole_data.pole_tolerance,
                                        diagnostics_recorder);

  const axom::Array<double> weight_correction =
    basis.evaluateExpansionAtChebyshevNodes(integral_coefficients);
  QuadratureRule fejer_rule = applyFejerWeightCorrection(basis, weight_correction);

  diagnostics_recorder.recordFinalRule(integral_coefficients,
                                       weight_correction,
                                       fejer_rule.view(),
                                       basis.basisColumns(),
                                       basis.basisColumnCount(),
                                       basis.nodeCount());

  return fejer_rule;
}

/// \brief Computes rational Fejer nodes and weights on [-1,1] with optional diagnostics.
inline void compute_rational_fejer_data_m11_impl(const RationalFejerPoleData& pole_data,
                                                 axom::Array<double>& nodes,
                                                 axom::Array<double>& weights,
                                                 RationalFejerDiagnostics* diagnostics = nullptr,
                                                 int allocatorID = axom::getDefaultAllocatorID())
{
  const QuadratureRule fejer_rule =
    compute_rational_fejer_rule_m11_impl(pole_data, diagnostics, allocatorID);
  copy_array_to_array(fejer_rule.nodes(), nodes, allocatorID);
  copy_array_to_array(fejer_rule.weights(), weights, allocatorID);
}

/// \brief Computes rational Fejer nodes and weights on [-1,1].
void compute_rational_fejer_data_m11(axom::ArrayView<const Complex> poles_m11,
                                     axom::Array<double>& nodes,
                                     axom::Array<double>& weights,
                                     int allocatorID)
{
  const RationalFejerPoleData pole_data = RationalFejerPoleData::fromM11(poles_m11);

  axom::Array<double> nodes_tmp;
  axom::Array<double> weights_tmp;
  compute_rational_fejer_data_m11_impl(pole_data, nodes_tmp, weights_tmp);
  copy_array_to_array(nodes_tmp, nodes, allocatorID);
  copy_array_to_array(weights_tmp, weights, allocatorID);
}

/// \brief Computes rational Fejer diagnostics on [-1,1].
void compute_rational_fejer_diagnostics_m11(axom::ArrayView<const Complex> poles_m11,
                                            RationalFejerDiagnostics& diagnostics,
                                            int allocatorID)
{
  axom::Array<double> nodes_m11;
  axom::Array<double> weights_m11;
  compute_rational_fejer_data_m11_impl(RationalFejerPoleData::fromM11(poles_m11),
                                       nodes_m11,
                                       weights_m11,
                                       &diagnostics,
                                       allocatorID);
}

/// \brief Builds a canonical cache key for an M11 rational Fejer pole sequence.
std::string make_rational_fejer_cache_key_m11(axom::ArrayView<const Complex> poles_m11,
                                              int allocatorID)
{
  return RationalFejerPoleData::fromM11(poles_m11).cacheKey(allocatorID);
}

}  // namespace internal

/// \brief Computes rational Fejer nodes and weights on [0,1].
void compute_rational_fejer_data(axom::ArrayView<const std::complex<double>> poles01,
                                 axom::Array<double>& nodes,
                                 axom::Array<double>& weights,
                                 int allocatorID)
{
  const internal::RationalFejerPoleData pole_data = internal::RationalFejerPoleData::from01(poles01);

  // The implementation below follows the Deckers et al. / van Deun et al.
  // rational Fejer and rational Chebyshev constructions. QuaHOG's MATLAB
  // Spectral_PE implementation served as a secondary reference for how these
  // pieces are assembled in the planar integration workflow.
  const QuadratureRule rule_m11 =
    internal::compute_rational_fejer_rule_m11_impl(pole_data, nullptr, allocatorID);

  internal::map_rule_m11_to_01(rule_m11.nodes(), rule_m11.weights(), nodes, weights, allocatorID);
}

/// \brief Computes rational Fejer diagnostics on [0,1].
void internal::compute_rational_fejer_diagnostics(axom::ArrayView<const std::complex<double>> poles01,
                                                  internal::RationalFejerDiagnostics& diagnostics,
                                                  int allocatorID)
{
  const internal::RationalFejerPoleData pole_data = internal::RationalFejerPoleData::from01(poles01);

  const QuadratureRule rule_m11 =
    internal::compute_rational_fejer_rule_m11_impl(pole_data, &diagnostics, allocatorID);
  internal::map_rule_m11_to_01(rule_m11.nodes(),
                               rule_m11.weights(),
                               diagnostics.nodes_01,
                               diagnostics.weights_01,
                               allocatorID);
}

/// \brief Returns a cached rational Fejer quadrature rule view on [0,1].
QuadratureRuleView get_rational_fejer(axom::ArrayView<const std::complex<double>> poles01,
                                      int allocatorID)
{
  const internal::RationalFejerPoleData pole_data = internal::RationalFejerPoleData::from01(poles01);

  constexpr std::size_t MAX_RATIONAL_FEJER_CACHED_RULES = 1u << 16;
  static axom::LRUCache<std::string, QuadratureRule> rule_library(MAX_RATIONAL_FEJER_CACHED_RULES);
  static std::mutex rule_library_mutex;

  const std::string key = pole_data.cacheKey(allocatorID);

  {
    const std::lock_guard<std::mutex> lock(rule_library_mutex);
    if(auto* rule = rule_library.find(key))
    {
      return rule->view();
    }
  }

  // Build the rule outside the cache mutex so a slow generated pole family
  // does not block unrelated cache hits. Recheck under the lock before
  // insertion in case another thread populated the same key meanwhile.
  axom::Array<double> nodes;
  axom::Array<double> weights;
  const QuadratureRule rule_m11 =
    internal::compute_rational_fejer_rule_m11_impl(pole_data, nullptr, allocatorID);
  internal::map_rule_m11_to_01(rule_m11.nodes(), rule_m11.weights(), nodes, weights, allocatorID);
  QuadratureRule rule {std::move(nodes), std::move(weights)};

  {
    const std::lock_guard<std::mutex> lock(rule_library_mutex);
    if(auto* cached_rule = rule_library.find(key))
    {
      return cached_rule->view();
    }

    auto& cached_rule = rule_library.insert(key, std::move(rule));
    return cached_rule.view();
  }
}

}  // namespace numerics
}  // namespace axom
