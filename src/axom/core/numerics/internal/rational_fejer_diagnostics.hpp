// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_NUMERICS_INTERNAL_RATIONAL_FEJER_DIAGNOSTICS_HPP_
#define AXOM_NUMERICS_INTERNAL_RATIONAL_FEJER_DIAGNOSTICS_HPP_

#include "axom/core/numerics/internal/rational_quadrature_common.hpp"

#include <utility>

namespace axom
{
namespace numerics
{
namespace internal
{

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

  void recordCompanionChebyshevRule(const QuadratureRuleView& companion_rule)
  {
    if(!enabled())
    {
      return;
    }

    copy_array_to_array(companion_rule.nodes(),
                        m_diagnostics->rational_chebyshev_nodes_m11,
                        m_allocatorID);
    copy_array_to_array(companion_rule.weights(),
                        m_diagnostics->rational_chebyshev_weights_m11,
                        m_allocatorID);
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

  void recordFinalRule(axom::ArrayView<const double> integral_coefficients,
                       axom::ArrayView<const double> weight_correction,
                       const QuadratureRuleView& fejer_rule,
                       axom::ArrayView<const double> basis_columns,
                       int basis_row_count,
                       int basis_col_count)
  {
    if(!enabled())
    {
      return;
    }

    copy_array_to_array(integral_coefficients, m_diagnostics->basis_coefficients, m_allocatorID);
    copy_array_to_array(weight_correction, m_diagnostics->weight_correction_m11, m_allocatorID);
    copy_array_to_array(fejer_rule.weights(), m_diagnostics->final_weights_m11, m_allocatorID);
    copy_array_to_array(basis_columns, m_diagnostics->basis_matrix_transpose_m11, m_allocatorID);
    m_diagnostics->basis_matrix_row_count = basis_row_count;
    m_diagnostics->basis_matrix_col_count = basis_col_count;
  }

private:
  RationalFejerDiagnostics* m_diagnostics {nullptr};
  int m_allocatorID {axom::getDefaultAllocatorID()};
  double m_poleTolerance {0.0};
  mutable int m_stepIndex {0};
};

}  // namespace internal
}  // namespace numerics
}  // namespace axom

#endif  // AXOM_NUMERICS_INTERNAL_RATIONAL_FEJER_DIAGNOSTICS_HPP_
