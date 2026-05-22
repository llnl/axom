// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/LRUCache.hpp"
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
 * completes missing complex conjugates, builds the companion rational
 * Chebyshev rule, seeds exact pole-moment data, and solves for rational basis
 * coefficients by enforcing orthogonality. The sampled basis expansion is then
 * combined with the rational Chebyshev weights to produce the final rational
 * Fejer nodes and weights. Public entry points map between Axom's `[0,1]`
 * rule domain and the `[-1,1]` reference domain used by the paper formulas.
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

/// \brief Records optional diagnostics for one rational Fejer coefficient solve.
class RationalFejerStepDiagnosticsRecorder
{
public:
  RationalFejerStepDiagnosticsRecorder(RationalFejerDiagnostics::Step* step_diagnostics = nullptr,
                                       int diagnostic_allocator_id = axom::getDefaultAllocatorID())
    : m_step(step_diagnostics)
    , m_allocatorID(diagnostic_allocator_id)
  { }

  bool enabled() const { return m_step != nullptr; }

  void recordWeightedRows(axom::ArrayView<const long double> row0,
                          axom::ArrayView<const long double> row1,
                          bool has_imaginary_component) const
  {
    if(!enabled())
    {
      return;
    }

    copy_array_to_array(row0, m_step->weighted_row0, m_allocatorID);
    if(has_imaginary_component)
    {
      copy_array_to_array(row1, m_step->weighted_row1, m_allocatorID);
    }
  }

  void recordProjectedRow0(axom::ArrayView<const long double> row0,
                           axom::ArrayView<const long double> row0_terms,
                           int node_count,
                           long double rhs0) const
  {
    if(!enabled())
    {
      return;
    }

    copy_array_to_array(row0, m_step->projected_row0, m_allocatorID);
    copy_array_to_array(row0_terms, m_step->projected_row0_terms, m_allocatorID);
    m_step->projected_row_terms_node_count = node_count;
    m_step->rhs0 = static_cast<double>(rhs0);
  }

  void recordProjectedRow1(axom::ArrayView<const long double> row1,
                           axom::ArrayView<const long double> row1_terms,
                           long double rhs1,
                           long double a00,
                           long double a01,
                           long double a10,
                           long double a11) const
  {
    if(!enabled())
    {
      return;
    }

    copy_array_to_array(row1, m_step->projected_row1, m_allocatorID);
    copy_array_to_array(row1_terms, m_step->projected_row1_terms, m_allocatorID);
    m_step->rhs1 = static_cast<double>(rhs1);
    m_step->a00 = static_cast<double>(a00);
    m_step->a01 = static_cast<double>(a01);
    m_step->a10 = static_cast<double>(a10);
    m_step->a11 = static_cast<double>(a11);
  }

  void recordTinySolve(bool used_pivot,
                       long double factor,
                       long double reduced_a11,
                       long double reduced_rhs1,
                       long double x0,
                       long double x1) const
  {
    if(!enabled())
    {
      return;
    }

    m_step->used_pivot = used_pivot;
    m_step->factor = static_cast<double>(factor);
    m_step->reduced_a11 = static_cast<double>(reduced_a11);
    m_step->reduced_rhs1 = static_cast<double>(reduced_rhs1);
    m_step->solved_x0 = static_cast<double>(x0);
    m_step->solved_x1 = static_cast<double>(x1);
  }

private:
  RationalFejerDiagnostics::Step* m_step {nullptr};
  int m_allocatorID {axom::getDefaultAllocatorID()};
};

/// \brief Stores the rational basis sampled at the rational Chebyshev nodes.
///
/// The moment solve and final weight assembly reuse this matrix. It
/// corresponds to the orthonormal rational basis used in the Deckers et al.
/// extended rational Fejer construction; QuaHOG provided a practical reference
/// implementation for this Axom port.
class RationalFejerBasis
{
public:
  RationalFejerBasis(axom::ArrayView<const Complex> cayley_poles,
                     axom::Array<double> nodes,
                     axom::Array<double> lambda)
    : m_nodes(std::move(nodes))
    , m_lambda(std::move(lambda))
    , m_QColumns(buildBasisColumns(cayley_poles, m_nodes))
    , m_basisColumnCount(static_cast<int>(m_nodes.size()))
    , m_nodeCount(static_cast<int>(m_nodes.size()))
  { }

  const axom::Array<double>& nodes() const { return m_nodes; }

  const axom::Array<double>& lambda() const { return m_lambda; }

  axom::ArrayView<const double> basisColumns() const { return m_QColumns.view(); }

  int basisColumnCount() const { return m_basisColumnCount; }

  int nodeCount() const { return m_nodeCount; }

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

    const int n = static_cast<int>(m_nodes.size());
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
      const Complex basis_value = evaluatePoleBasisValue(pole, multiplicity, m_nodes[i]);

      weighted_row0[i] =
        static_cast<long double>(m_lambda[i]) * static_cast<long double>(basis_value.real());
      if(has_imaginary_component)
      {
        weighted_row1[i] =
          static_cast<long double>(m_lambda[i]) * static_cast<long double>(-basis_value.imag());
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

  axom::Array<double> evaluateExpansion(axom::ArrayView<const double> coefficients) const
  {
    axom::Array<double> values(m_nodes.size(), m_nodes.size());
    values.fill(0.0);
    for(int j = 0; j < coefficients.size(); ++j)
    {
      const double coefficient = coefficients[j];
      if(coefficient == 0.0)
      {
        continue;
      }

      for(int i = 0; i < m_nodes.size(); ++i)
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

  axom::Array<double> m_nodes;
  axom::Array<double> m_lambda;
  axom::Array<double> m_QColumns;
  int m_basisColumnCount {0};
  int m_nodeCount {0};
};

/// \brief Copies the sampled rational Fejer basis matrix into diagnostics storage.
///
/// The surrounding construction follows Deckers, Mougaida, and Belhadjsalah,
/// ACM TOMS 43(4), 2017 (Algorithm 973). Axom uses the explicit node/weight
/// construction here, not the semi-automatic adaptive integrator from that paper.
inline void copy_basis_matrix_to_array(const RationalFejerBasis& basis,
                                       axom::Array<double>& array,
                                       int allocatorID)
{
  copy_array_to_array(basis.basisColumns(), array, allocatorID);
}

/// \brief Counts coefficient-solve steps for a canonical rational Fejer pole sequence.
inline int countRationalFejerSteps(axom::ArrayView<const Pole> canonical_poles, double pole_tolerance)
{
  int step_count = 0;
  int coeff_index = 1;
  const int num_poles = static_cast<int>(canonical_poles.size());
  while(coeff_index < num_poles + 1)
  {
    const Pole pole = canonical_poles[coeff_index - 1];
    const int component_count = pole.componentCount(pole_tolerance);
    ++step_count;
    coeff_index += component_count;
  }
  return step_count;
}

/// \brief Records optional diagnostics for the full rational Fejer construction.
class RationalFejerDiagnosticsRecorder
{
public:
  RationalFejerDiagnosticsRecorder(RationalFejerDiagnostics* diagnostics,
                                   int allocatorID,
                                   double pole_tolerance)
    : m_diagnostics(diagnostics)
    , m_allocatorID(allocatorID)
    , m_poleTolerance(pole_tolerance)
  { }

  bool enabled() const { return m_diagnostics != nullptr; }

  void recordPoleSetup(axom::ArrayView<const Pole> canonical_poles,
                       axom::ArrayView<const Complex> cayley_poles)
  {
    if(!enabled())
    {
      return;
    }

    m_diagnostics->canonical_poles_m11 =
      axom::Array<Complex>(canonical_poles.size(), canonical_poles.size(), m_allocatorID);
    for(axom::IndexType i = 0; i < canonical_poles.size(); ++i)
    {
      m_diagnostics->canonical_poles_m11[i] = canonical_poles[i].value();
    }

    m_diagnostics->cayley_poles =
      axom::Array<Complex>(cayley_poles.size(), cayley_poles.size(), m_allocatorID);
    for(axom::IndexType i = 0; i < cayley_poles.size(); ++i)
    {
      m_diagnostics->cayley_poles[i] = cayley_poles[i];
    }

    const int diagnostic_step_count = countRationalFejerSteps(canonical_poles, m_poleTolerance);
    m_diagnostics->steps = axom::Array<RationalFejerDiagnostics::Step>(diagnostic_step_count,
                                                                       diagnostic_step_count,
                                                                       m_allocatorID);
  }

  void recordChebyshevSeedRule(axom::ArrayView<const double> nodes,
                               axom::ArrayView<const double> weights)
  {
    if(!enabled())
    {
      return;
    }

    copy_array_to_array(nodes, m_diagnostics->rational_chebyshev_nodes_m11, m_allocatorID);
    copy_array_to_array(weights, m_diagnostics->rational_chebyshev_weights_m11, m_allocatorID);
  }

  RationalFejerDiagnostics::Step beginStep(int coeff_index,
                                           int component_count,
                                           int pole_multiplicity_so_far,
                                           const Pole& pole,
                                           axom::ArrayView<const double> known_integrals) const
  {
    RationalFejerDiagnostics::Step step;
    if(!enabled())
    {
      return step;
    }

    step.step_index = m_stepIndex;
    step.coeff_index_before = coeff_index;
    step.component_count = component_count;
    step.pole_multiplicity_so_far = pole_multiplicity_so_far;
    step.pole_m11 = pole.value();
    copy_array_to_array(known_integrals, step.basis_coefficients_before, m_allocatorID);
    return step;
  }

  RationalFejerStepDiagnosticsRecorder stepRecorder(RationalFejerDiagnostics::Step& step) const
  {
    return RationalFejerStepDiagnosticsRecorder {enabled() ? &step : nullptr, m_allocatorID};
  }

  void finishStep(RationalFejerDiagnostics::Step&& step,
                  axom::ArrayView<const double> orthogonal_integrals,
                  axom::ArrayView<const double> known_integrals)
  {
    if(!enabled())
    {
      return;
    }

    copy_array_to_array(orthogonal_integrals, step.orthogonal_integrals, m_allocatorID);
    copy_array_to_array(known_integrals, step.basis_coefficients_after, m_allocatorID);
    m_diagnostics->steps[m_stepIndex++] = std::move(step);
  }

  void recordFinalRule(const RationalFejerBasis& basis,
                       axom::ArrayView<const double> basis_coefficients,
                       axom::ArrayView<const double> basis_expansion,
                       axom::ArrayView<const double> weights)
  {
    if(!enabled())
    {
      return;
    }

    copy_array_to_array(basis_coefficients, m_diagnostics->basis_coefficients, m_allocatorID);
    copy_array_to_array(basis_expansion, m_diagnostics->basis_expansion_m11, m_allocatorID);
    copy_array_to_array(weights, m_diagnostics->final_weights_m11, m_allocatorID);
    copy_basis_matrix_to_array(basis, m_diagnostics->basis_matrix_transpose_m11, m_allocatorID);
    m_diagnostics->basis_matrix_row_count = basis.basisColumnCount();
    m_diagnostics->basis_matrix_col_count = basis.nodeCount();
  }

private:
  RationalFejerDiagnostics* m_diagnostics {nullptr};
  int m_allocatorID {axom::getDefaultAllocatorID()};
  double m_poleTolerance {0.0};
  int m_stepIndex {0};
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

/// \brief Initializes rational Fejer basis coefficients with the constant term.
inline axom::Array<double> initializeRationalFejerBasisCoefficients(int num_poles)
{
  axom::Array<double> basis_coefficients(num_poles + 1, num_poles + 1);
  basis_coefficients.fill(0.0);
  basis_coefficients[0] = 2.0 / std::sqrt(M_PI);
  return basis_coefficients;
}

/// \brief Seeds exact non-orthogonal moments for the next newly encountered pole.
inline int seedNewPoleMomentCoefficients(const PoleSequence& canonical_poles,
                                         const Pole& pole,
                                         int coeff_index,
                                         int component_count,
                                         double pole_tolerance,
                                         axom::Array<double>& basis_coefficients)
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
    basis_coefficients[coeff_index + repeated_pole_offsets[i] - 1] = pole_integrals[i];
  }

  return pole_multiplicity_so_far;
}

/// \brief Solves the rational Fejer basis coefficients by enforcing orthogonality.
inline axom::Array<double> solveRationalFejerBasisCoefficients(
  const PoleSequence& canonical_poles,
  const RationalFejerBasis& basis,
  double pole_tolerance,
  RationalFejerDiagnosticsRecorder& diagnostics_recorder)
{
  const int num_poles = static_cast<int>(canonical_poles.size());
  axom::Array<double> basis_coefficients = initializeRationalFejerBasisCoefficients(num_poles);

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
                                                                       basis_coefficients);

    const int num_known_columns = coeff_index + component_count;
    auto known_integrals = basis_coefficients.view().subspan(0, num_known_columns);

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
      basis_coefficients[coeff_index + i] = orthogonal_integrals[i];
    }

    diagnostics_recorder.finishStep(std::move(step_diagnostics), orthogonal_integrals, known_integrals);

    coeff_index += component_count;
  }

  return basis_coefficients;
}

/// \brief Reconstructs quadrature weights from a sampled basis expansion.
inline axom::Array<double> assembleRationalFejerWeights(const RationalFejerBasis& basis,
                                                        axom::ArrayView<const double> basis_expansion)
{
  axom::Array<double> weights(basis_expansion.size(), basis_expansion.size());
  const auto& lambda = basis.lambda();
  // Convert the solved expansion coefficients back into quadrature weights
  // by multiplying the basis expansion by the rational Chebyshev weights.
  for(int i = 0; i < weights.size(); ++i)
  {
    weights[i] = basis_expansion[i] * lambda[i];
  }
  return weights;
}

/// \brief Computes rational Fejer nodes and weights on [-1,1] with optional diagnostics.
inline void compute_rational_fejer_data_m11_impl(const RationalFejerPoleData& pole_data,
                                                 axom::Array<double>& nodes,
                                                 axom::Array<double>& weights,
                                                 RationalFejerDiagnostics* diagnostics = nullptr,
                                                 int allocatorID = axom::getDefaultAllocatorID())
{
  RationalFejerDiagnosticsRecorder diagnostics_recorder {diagnostics,
                                                         allocatorID,
                                                         pole_data.pole_tolerance};
  diagnostics_recorder.recordPoleSetup(pole_data.canonical_poles_m11.view(), pole_data.cayley_poles);

  // A rational Fejer rule for m finite/infinite poles uses m+1 rational Chebyshev nodes.
  // The appended infinity pole supplies that extra Chebyshev-like degree.
  const PoleSequence rational_chebyshev_poles = pole_data.canonical_poles_m11.withAppendedInfinity();

  axom::Array<double> rational_chebyshev_nodes;
  axom::Array<double> rational_chebyshev_weights;
  compute_rational_chebyshev_data(rational_chebyshev_poles.view(),
                                  rational_chebyshev_nodes,
                                  rational_chebyshev_weights);

  diagnostics_recorder.recordChebyshevSeedRule(rational_chebyshev_nodes, rational_chebyshev_weights);

  RationalFejerBasis basis(pole_data.cayley_poles,
                           std::move(rational_chebyshev_nodes),
                           std::move(rational_chebyshev_weights));

  const axom::Array<double> basis_coefficients =
    solveRationalFejerBasisCoefficients(pole_data.canonical_poles_m11,
                                        basis,
                                        pole_data.pole_tolerance,
                                        diagnostics_recorder);

  nodes = basis.nodes();
  const axom::Array<double> basis_expansion = basis.evaluateExpansion(basis_coefficients);
  weights = assembleRationalFejerWeights(basis, basis_expansion);

  diagnostics_recorder.recordFinalRule(basis, basis_coefficients, basis_expansion, weights);
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

  axom::Array<double> nodes_m11;
  axom::Array<double> weights_m11;

  // The implementation below follows the Deckers et al. / van Deun et al.
  // rational Fejer and rational Chebyshev constructions. QuaHOG's MATLAB
  // Spectral_PE implementation served as a secondary reference for how these
  // pieces are assembled in the planar integration workflow.
  internal::compute_rational_fejer_data_m11_impl(pole_data, nodes_m11, weights_m11, nullptr, allocatorID);

  internal::map_rule_m11_to_01(nodes_m11, weights_m11, nodes, weights, allocatorID);
}

/// \brief Computes rational Fejer diagnostics on [0,1].
void internal::compute_rational_fejer_diagnostics(axom::ArrayView<const std::complex<double>> poles01,
                                                  internal::RationalFejerDiagnostics& diagnostics,
                                                  int allocatorID)
{
  const internal::RationalFejerPoleData pole_data = internal::RationalFejerPoleData::from01(poles01);

  axom::Array<double> nodes_m11;
  axom::Array<double> weights_m11;
  internal::compute_rational_fejer_data_m11_impl(pole_data,
                                                 nodes_m11,
                                                 weights_m11,
                                                 &diagnostics,
                                                 allocatorID);
  internal::map_rule_m11_to_01(nodes_m11,
                               weights_m11,
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
  axom::Array<double> nodes_m11;
  axom::Array<double> weights_m11;
  internal::compute_rational_fejer_data_m11_impl(pole_data, nodes_m11, weights_m11, nullptr, allocatorID);
  internal::map_rule_m11_to_01(nodes_m11, weights_m11, nodes, weights, allocatorID);
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
