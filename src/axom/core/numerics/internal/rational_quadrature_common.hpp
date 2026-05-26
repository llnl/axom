// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_NUMERICS_INTERNAL_RATIONAL_QUADRATURE_COMMON_HPP_
#define AXOM_NUMERICS_INTERNAL_RATIONAL_QUADRATURE_COMMON_HPP_

#include "axom/config.hpp"
#include "axom/core/Array.hpp"
#include "axom/core/NumericLimits.hpp"
#include "axom/core/numerics/internal/rational_quadrature.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/fmt.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <utility>

namespace axom
{
namespace numerics
{
namespace internal
{
/// \brief Reports a rational quadrature precondition failure and aborts the process.
[[noreturn]] inline void failRationalQuadraturePrecondition(const std::string& message)
{
  std::cerr << "ERROR: " << message << "\n";
  axom::utilities::processAbort();
  std::abort();
}

/// \brief Computes an integer power of a complex value.
inline Complex powInteger(Complex base, int exponent)
{
  Complex result {1.0, 0.0};
  while(exponent > 0)
  {
    if((exponent & 1) != 0)
    {
      result *= base;
    }
    exponent >>= 1;
    if(exponent > 0)
    {
      base *= base;
    }
  }
  return result;
}

/// \brief Lightweight wrapper for poles in the reference [-1,1] domain.
///
/// Keeps pole-specific tolerance checks, infinity normalization, and the
/// Cayley map in one place so the construction code can use pole-level
/// operations instead of scattered complex-number checks.
class Pole
{
public:
  Pole() = default;
  Pole(double real, double imag) : m_value(real, imag) { }
  explicit Pole(const Complex& value) : m_value(value) { }

  /// \brief Infinity sentinel used for poles with infinite or very large magnitude.
  static Pole infinity() { return Pole {std::numeric_limits<double>::infinity(), 0.0}; }

  double real() const { return m_value.real(); }
  double imag() const { return m_value.imag(); }
  const Complex& value() const { return m_value; }

  bool hasNaN() const { return std::isnan(m_value.real()) || std::isnan(m_value.imag()); }

  bool hasFiniteCoordinates() const
  {
    return std::isfinite(m_value.real()) && std::isfinite(m_value.imag());
  }

  bool isInfinite() const
  {
    return std::isinf(m_value.real()) || std::isinf(m_value.imag()) || std::isinf(std::abs(m_value));
  }

  bool isEffectivelyReal(double tol) const { return std::abs(m_value.imag()) <= tol; }

  bool liesOnRealInterval(double lower, double upper, double tol) const
  {
    return isEffectivelyReal(tol) && m_value.real() >= lower - tol && m_value.real() <= upper + tol;
  }

  int componentCount(double tol) const { return isEffectivelyReal(tol) ? 1 : 2; }

  /// \brief Compares poles with relative distance where possible.
  ///
  /// Handles zero and infinity explicitly so duplicate poles can be
  /// canonicalized robustly.
  bool closeTo(const Pole& other, double tol) const
  {
    if(isInfinite() || other.isInfinite())
    {
      return isInfinite() && other.isInfinite();
    }

    if(std::abs(m_value) <= tol || std::abs(other.m_value) <= tol)
    {
      return std::abs(m_value - other.m_value) <= tol;
    }

    return std::abs(Complex {1.0, 0.0} - other.m_value / m_value) < tol;
  }

  Pole normalizedInfinite(double threshold) const
  {
    return isInfinite() || std::abs(m_value) >= threshold ? Pole::infinity() : *this;
  }

  /// \brief Returns the upper-half-plane representative of a conjugate pair.
  ///
  /// Rational Chebyshev formulas group conjugate pairs this way before they
  /// are expanded back into real-valued components.
  Pole withPositiveImaginaryMagnitude() const
  {
    return Pole {Complex {m_value.real(), std::abs(m_value.imag())}};
  }

  Pole conjugate() const { return Pole {std::conj(m_value)}; }

  Complex cayleyTransform() const
  {
    // Map poles outside [-1,1] into the unit disk. Infinite poles map to the
    // origin, matching the Chebyshev polynomial limit of the rational formulas.
    if(isInfinite())
    {
      return Complex {0.0, 0.0};
    }

    const double sign = m_value.real() >= 0.0 ? 1.0 : -1.0;
    return Complex {1.0, 0.0} / (m_value + sign * std::sqrt(m_value * m_value - Complex {1.0, 0.0}));
  }

  Complex cayleyTransformPreservingImaginarySign() const
  {
    Complex cayley_value = withPositiveImaginaryMagnitude().cayleyTransform();
    return m_value.imag() < 0.0 ? std::conj(cayley_value) : cayley_value;
  }

  friend bool operator<(const Pole& lhs, const Pole& rhs)
  {
    const bool lhs_inf = lhs.isInfinite();
    const bool rhs_inf = rhs.isInfinite();
    if(lhs_inf != rhs_inf)
    {
      return !lhs_inf && rhs_inf;
    }

    if(lhs.real() != rhs.real())
    {
      return lhs.real() < rhs.real();
    }

    return lhs.imag() < rhs.imag();
  }

private:
  Complex m_value {0.0, 0.0};
};

struct DistinctPoleData
{
  axom::Array<Pole> values;
  axom::Array<int> multiplicities;
};

/// \brief Validates a rational Fejer pole against the construction domain.
inline void validatePole(const Pole& pole,
                         axom::IndexType pole_index,
                         double lower,
                         double upper,
                         double interval_tol,
                         const char* domain_name)
{
  if(pole.hasNaN())
  {
    failRationalQuadraturePrecondition(
      axom::fmt::format("Rational quadrature pole {} in {} contains NaN.", pole_index, domain_name));
  }

  if(pole.isInfinite())
  {
    return;
  }

  if(!pole.hasFiniteCoordinates())
  {
    failRationalQuadraturePrecondition(
      axom::fmt::format("Rational quadrature pole {} in {} is neither finite nor infinite.",
                        pole_index,
                        domain_name));
  }

  if(pole.liesOnRealInterval(lower, upper, interval_tol))
  {
    failRationalQuadraturePrecondition(
      axom::fmt::format("Rational quadrature finite pole {} = ({}, {}) lies on the {} interval.",
                        pole_index,
                        pole.real(),
                        pole.imag(),
                        domain_name));
  }
}

/// \brief Interleaves matching offsets from paired real and imaginary components.
inline axom::Array<int> interleaveOffsets(axom::ArrayView<const int> even_offsets,
                                          axom::ArrayView<const int> odd_offsets)
{
  assert(even_offsets.size() == odd_offsets.size());
  axom::Array<int> interleaved(2 * even_offsets.size(), 2 * even_offsets.size());
  for(axom::IndexType i = 0; i < even_offsets.size(); ++i)
  {
    interleaved[2 * i] = even_offsets[i];
    interleaved[2 * i + 1] = odd_offsets[i];
  }
  return interleaved;
}

/// \brief Returns the pole-matching tolerance used by rational Fejer construction.
inline double rationalFejerPoleTolerance() { return 2.0 * axom::numeric_limits<double>::epsilon(); }

/// \brief Returns the finite-pole magnitude threshold treated as infinity.
inline double rationalFejerInfinitePoleThreshold()
{
  return 2.0 / axom::numeric_limits<double>::epsilon();
}

/// \brief Collapses a pole sequence into distinct representatives plus multiplicities.
///
/// This is the form consumed by the rational Chebyshev and rational Fejer builders.
inline DistinctPoleData collectDistinctPoles(axom::ArrayView<const Pole> poles, double tol)
{
  axom::Array<Pole> sorted(poles);
  std::sort(sorted.begin(), sorted.end());

  DistinctPoleData distinct;
  // Collapse repeated or near-repeated poles into one representative value and
  // record the multiplicity that the rational basis construction should see.
  distinct.values.reserve(sorted.size());
  distinct.multiplicities.reserve(sorted.size());
  for(const auto& pole : sorted)
  {
    if(distinct.values.empty() || !distinct.values[distinct.values.size() - 1].closeTo(pole, tol))
    {
      distinct.values.push_back(pole);
      distinct.multiplicities.push_back(1);
    }
    else
    {
      ++distinct.multiplicities[distinct.multiplicities.size() - 1];
    }
  }

  return distinct;
}

/// \brief Owning sequence of rational quadrature poles.
///
/// This class keeps general pole-sequence transformations close to the pole
/// representation without taking ownership of algorithm-specific steps such as
/// the Fejer coefficient solve.
class PoleSequence
{
public:
  PoleSequence() = default;

  explicit PoleSequence(axom::Array<Pole> poles) : m_poles(std::move(poles)) { }

  explicit PoleSequence(axom::ArrayView<const Pole> poles) : m_poles(poles) { }

  /// \brief Converts complex pole values into a PoleSequence.
  static PoleSequence fromComplex(axom::ArrayView<const Complex> pole_values)
  {
    axom::Array<Pole> poles(pole_values.size(), pole_values.size());
    for(axom::IndexType i = 0; i < pole_values.size(); ++i)
    {
      poles[i] = Pole {pole_values[i]};
    }
    return PoleSequence {std::move(poles)};
  }

  axom::ArrayView<const Pole> view() const { return m_poles.view(); }

  const axom::Array<Pole>& array() const { return m_poles; }

  axom::IndexType size() const { return m_poles.size(); }

  bool empty() const { return m_poles.empty(); }

  const Pole& operator[](axom::IndexType idx) const { return m_poles[idx]; }

  /// \brief Copies this sequence to an owning pole array.
  axom::Array<Pole> toPoleArray() const { return axom::Array<Pole>(m_poles.view()); }

  /// \brief Maps poles from [0,1] coordinates to [-1,1] coordinates.
  PoleSequence mapped01ToM11() const
  {
    axom::Array<Pole> result(m_poles.size(), m_poles.size());
    for(axom::IndexType i = 0; i < m_poles.size(); ++i)
    {
      const Pole& pole = m_poles[i];
      result[i] =
        pole.isInfinite() ? Pole::infinity() : Pole {2.0 * pole.value() - Complex {1.0, 0.0}};
    }
    return PoleSequence {std::move(result)};
  }

  /// \brief Validates the sequence against a construction domain.
  void validate(double lower, double upper, const char* domain_name) const
  {
    const double interval_tol = 16.0 * axom::numeric_limits<double>::epsilon();

    if(m_poles.empty())
    {
      failRationalQuadraturePrecondition(
        axom::fmt::format("Rational quadrature requires at least one pole in {}.", domain_name));
    }

    for(axom::IndexType i = 0; i < m_poles.size(); ++i)
    {
      validatePole(m_poles[i], i, lower, upper, interval_tol, domain_name);
    }
  }

  /// \brief Returns the canonical conjugate-complete form of this sequence.
  ///
  /// In this context, "canonical" means a deterministic pole sequence that the
  /// real-valued rational quadrature construction can consume directly:
  /// duplicate or nearly duplicate poles are represented by the first matching
  /// value seen in caller order, and each non-real pole is immediately followed
  /// by its conjugate. The input order of distinct pole blocks is otherwise
  /// preserved rather than sorted, since pole order is part of the fixed-pole
  /// rational basis and therefore part of the rule/cache key.
  PoleSequence canonicalized(double tol) const
  {
    axom::Array<Pole> poles = toPoleArray();
    const int num_input_poles = static_cast<int>(poles.size());
    axom::Array<int> consumed(num_input_poles, num_input_poles);
    consumed.fill(0);

    axom::Array<Pole> ordered_poles;
    ordered_poles.reserve(2 * num_input_poles);
    for(int i = 0; i < num_input_poles; ++i)
    {
      if(consumed[i])
      {
        continue;
      }

      const Pole pole = poles[i];
      // Coalesce later near-duplicates onto the first representative. This
      // keeps repeated poles deterministic without changing the distinct-block order.
      for(int idx = i + 1; idx < num_input_poles; ++idx)
      {
        if(!consumed[idx] && pole.closeTo(poles[idx], tol))
        {
          poles[idx] = pole;
        }
      }

      ordered_poles.push_back(pole);
      consumed[i] = 1;

      if(!pole.isEffectivelyReal(tol))
      {
        // The real-valued quadrature rule is built from conjugate-complete
        // pole blocks. If the caller supplies only one side of a complex pair,
        // synthesize the missing partner here; repeated poles consume one
        // matching conjugate per occurrence so multiplicities stay balanced.
        const auto conjugate_pole = pole.conjugate();
        for(int idx = i + 1; idx < num_input_poles; ++idx)
        {
          if(!consumed[idx] && conjugate_pole.closeTo(poles[idx], tol))
          {
            consumed[idx] = 1;
            break;
          }
        }

        ordered_poles.push_back(conjugate_pole);
      }
    }

    return PoleSequence {std::move(ordered_poles)};
  }

  /// \brief Returns a sequence with very large finite poles replaced by infinity.
  PoleSequence normalizedInfinity(double threshold) const
  {
    axom::Array<Pole> result = toPoleArray();
    for(auto& pole : result)
    {
      pole = pole.normalizedInfinite(threshold);
    }
    return PoleSequence {std::move(result)};
  }

  /// \brief Returns a canonical sequence with very large finite poles replaced by infinity.
  PoleSequence canonicalizedAndNormalized(double pole_tolerance, double infinite_pole_threshold) const
  {
    return canonicalized(pole_tolerance).normalizedInfinity(infinite_pole_threshold);
  }

  /// \brief Returns the pole form expected by the rational Chebyshev construction.
  ///
  /// Algorithm 882 evaluates complex-conjugate pairs through upper-half-plane
  /// representatives and treats very large finite poles as the Chebyshev
  /// polynomial limit.
  PoleSequence normalizedForRationalChebyshev(double infinite_pole_threshold) const
  {
    axom::Array<Pole> result = toPoleArray();
    for(auto& pole : result)
    {
      pole = pole.withPositiveImaginaryMagnitude().normalizedInfinite(infinite_pole_threshold);
    }
    return PoleSequence {std::move(result)};
  }

  /// \brief Returns this sequence with a trailing infinity pole.
  PoleSequence withAppendedInfinity() const
  {
    axom::Array<Pole> result(m_poles.size() + 1, m_poles.size() + 1);
    for(axom::IndexType i = 0; i < m_poles.size(); ++i)
    {
      result[i] = m_poles[i];
    }
    result[m_poles.size()] = Pole::infinity();
    return PoleSequence {std::move(result)};
  }

  /// \brief Returns the cyclic sequence used by the rational Chebyshev construction.
  PoleSequence cyclic() const
  {
    // Algorithm 882 uses the cyclic sequence (p_1, ..., p_n, p_1, ..., p_{n-1})
    // in its rational Chebyshev phase and weight formulas.
    const axom::IndexType num_cyclic_poles = 2 * m_poles.size() - 1;
    axom::Array<Pole> cyclic_poles(num_cyclic_poles, num_cyclic_poles);
    for(axom::IndexType i = 0; i < m_poles.size(); ++i)
    {
      cyclic_poles[i] = m_poles[i];
    }
    for(axom::IndexType i = 0; i < m_poles.size() - 1; ++i)
    {
      cyclic_poles[m_poles.size() + i] = m_poles[i];
    }
    return PoleSequence {std::move(cyclic_poles)};
  }

  /// \brief Collapses the sequence into distinct representatives plus multiplicities.
  DistinctPoleData distinct(double tol) const { return collectDistinctPoles(m_poles.view(), tol); }

  /// \brief Counts poles in a subrange that match a target pole.
  int countMatching(const Pole& target, int begin, int end, double tol) const
  {
    int count = 0;
    for(int idx = begin; idx < end; ++idx)
    {
      if(target.closeTo(m_poles[idx], tol))
      {
        ++count;
      }
    }
    return count;
  }

  /// \brief Collects one-based offsets of poles in a subrange that match a target pole.
  axom::Array<int> matchingOffsets(const Pole& target, int begin, int end, double tol) const
  {
    axom::Array<int> offsets;
    offsets.reserve(end - begin);

    for(int idx = begin; idx < end; ++idx)
    {
      if(target.closeTo(m_poles[idx], tol))
      {
        // Offsets are one-based relative positions in the coefficient block
        // beginning at `begin`; the Deckers moment data is written into those
        // future coefficient slots.
        offsets.push_back(idx - begin + 1);
      }
    }
    return offsets;
  }

private:
  axom::Array<Pole> m_poles;
};

/// \brief Maps a rule from [-1,1] coordinates to [0,1] coordinates.
inline void map_rule_m11_to_01(axom::ArrayView<const double> nodes_m11,
                               axom::ArrayView<const double> weights_m11,
                               axom::Array<double>& nodes,
                               axom::Array<double>& weights,
                               int allocatorID)
{
  // Rescale a rule from [-1,1] to [0,1]. The public Fejer APIs export
  // unit-interval rules because the downstream curve-parameter integrators
  // operate on Bezier and NURBS parameter spans in that domain.
  const int npts = static_cast<int>(nodes_m11.size());
  nodes = axom::Array<double>(npts, npts, allocatorID);
  weights = axom::Array<double>(npts, npts, allocatorID);

  for(int i = 0; i < npts; ++i)
  {
    nodes[i] = 0.5 * (nodes_m11[i] + 1.0);
    weights[i] = 0.5 * weights_m11[i];
  }
}

/// \brief Copies double values into an axom::Array with the requested allocator.
inline void copy_array_to_array(axom::ArrayView<const double> values,
                                axom::Array<double>& array,
                                int allocatorID)
{
  array = axom::Array<double>(values, allocatorID);
}

/// \brief Copies long double values into a double axom::Array with the requested allocator.
inline void copy_array_to_array(axom::ArrayView<const long double> values,
                                axom::Array<double>& array,
                                int allocatorID)
{
  array = axom::Array<double>(values.size(), values.size(), allocatorID);
  for(int i = 0; i < values.size(); ++i)
  {
    array[i] = static_cast<double>(values[i]);
  }
}

/// \brief Builds the rational Chebyshev rule used by the extended rational Fejer rule.
void compute_rational_chebyshev_data(axom::ArrayView<const Pole> poles_in,
                                     axom::Array<double>& nodes,
                                     axom::Array<double>& weights);

}  // namespace internal
}  // namespace numerics
}  // namespace axom

#endif  // AXOM_NUMERICS_INTERNAL_RATIONAL_QUADRATURE_COMMON_HPP_
