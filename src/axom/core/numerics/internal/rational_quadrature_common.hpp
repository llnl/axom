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

namespace axom
{
namespace numerics
{
namespace internal
{

/// \brief Returns the square of a scalar value.
inline double square(double value) { return value * value; }

/// \brief Reports a rational Fejer precondition failure and aborts the process.
[[noreturn]] inline void failRationalFejerPrecondition(const std::string& message)
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

/// \brief Converts complex pole values to Pole wrappers.
inline axom::Array<Pole> makePoleArray(axom::ArrayView<const Complex> pole_values)
{
  axom::Array<Pole> poles(pole_values.size(), pole_values.size());
  for(axom::IndexType i = 0; i < pole_values.size(); ++i)
  {
    poles[i] = Pole {pole_values[i]};
  }
  return poles;
}

/// \brief Maps poles from [0,1] coordinates to [-1,1] coordinates.
inline axom::Array<Complex> mapUnitIntervalPolesToM11(axom::ArrayView<const Complex> poles01)
{
  // `_m11` denotes the `[-1,1]` interval ("minus-one-to-one"). Public callers
  // provide poles on `[0,1]`, then the implementation maps them into the
  // symmetric interval used by the rational Chebyshev and Fejer formulas.
  axom::Array<Complex> poles_m11(poles01.size(), poles01.size());
  for(axom::IndexType i = 0; i < poles01.size(); ++i)
  {
    const Pole pole {poles01[i]};
    poles_m11[i] =
      pole.isInfinite() ? Pole::infinity().value() : 2.0 * pole.value() - Complex {1.0, 0.0};
  }
  return poles_m11;
}

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
    failRationalFejerPrecondition(
      axom::fmt::format("Rational Fejer pole {} in {} contains NaN.", pole_index, domain_name));
  }

  if(pole.isInfinite())
  {
    return;
  }

  if(!pole.hasFiniteCoordinates())
  {
    failRationalFejerPrecondition(
      axom::fmt::format("Rational Fejer pole {} in {} is neither finite nor infinite.",
                        pole_index,
                        domain_name));
  }

  if(pole.liesOnRealInterval(lower, upper, interval_tol))
  {
    failRationalFejerPrecondition(
      axom::fmt::format("Rational Fejer finite pole {} = ({}, {}) lies on the {} interval.",
                        pole_index,
                        pole.real(),
                        pole.imag(),
                        domain_name));
  }
}

/// \brief Validates a complex-valued pole sequence against the construction domain.
inline void validatePoleSequence(axom::ArrayView<const Complex> poles,
                                 double lower,
                                 double upper,
                                 const char* domain_name)
{
  const double interval_tol = 16.0 * axom::numeric_limits<double>::epsilon();

  if(poles.empty())
  {
    failRationalFejerPrecondition(
      axom::fmt::format("Rational Fejer quadrature requires at least one pole in {}.", domain_name));
  }

  for(axom::IndexType i = 0; i < poles.size(); ++i)
  {
    validatePole(Pole {poles[i]}, i, lower, upper, interval_tol, domain_name);
  }
}

/// \brief Validates a pole sequence against the construction domain.
inline void validatePoleSequence(axom::ArrayView<const Pole> poles,
                                 double lower,
                                 double upper,
                                 const char* domain_name)
{
  const double interval_tol = 16.0 * axom::numeric_limits<double>::epsilon();

  if(poles.empty())
  {
    failRationalFejerPrecondition(
      axom::fmt::format("Rational Fejer quadrature requires at least one pole in {}.", domain_name));
  }

  for(axom::IndexType i = 0; i < poles.size(); ++i)
  {
    validatePole(poles[i], i, lower, upper, interval_tol, domain_name);
  }
}

/// \brief Canonicalizes poles by coalescing near-duplicates and keeping
/// complex-conjugate partners adjacent.
///
/// The result is the stable internal representation used for caching,
/// diagnostics, and the rational recurrences.
inline axom::Array<Pole> canonicalizePoleSequence(axom::Array<Pole> poles, double tol)
{
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

  return ordered_poles;
}

/// \brief Canonicalizes a pole sequence and normalizes very large poles to infinity.
inline axom::Array<Pole> canonicalizeAndNormalizePoleSequence(axom::ArrayView<const Pole> poles,
                                                              double pole_tolerance,
                                                              double infinite_pole_threshold)
{
  axom::Array<Pole> canonical_poles =
    canonicalizePoleSequence(axom::Array<Pole>(poles), pole_tolerance);
  for(auto& pole : canonical_poles)
  {
    pole = pole.normalizedInfinite(infinite_pole_threshold);
  }
  return canonical_poles;
}

/// \brief Appends an infinity pole to a pole sequence.
inline axom::Array<Pole> appendInfinityPole(axom::ArrayView<const Pole> poles)
{
  axom::Array<Pole> result(poles.size() + 1, poles.size() + 1);
  for(axom::IndexType i = 0; i < poles.size(); ++i)
  {
    result[i] = poles[i];
  }
  result[poles.size()] = Pole::infinity();
  return result;
}

/// \brief Forms the cyclic pole sequence used by Algorithm 882.
inline axom::Array<Pole> makeCyclicPoleSequence(axom::ArrayView<const Pole> poles)
{
  // Algorithm 882 uses the cyclic sequence (p_1, ..., p_n, p_1, ..., p_{n-1})
  // in its rational Chebyshev phase and weight formulas.
  const axom::IndexType num_cyclic_poles = 2 * poles.size() - 1;
  axom::Array<Pole> cyclic_poles(num_cyclic_poles, num_cyclic_poles);
  for(axom::IndexType i = 0; i < poles.size(); ++i)
  {
    cyclic_poles[i] = poles[i];
  }
  for(axom::IndexType i = 0; i < poles.size() - 1; ++i)
  {
    cyclic_poles[poles.size() + i] = poles[i];
  }
  return cyclic_poles;
}

/// \brief Counts poles in a subrange that match a target pole.
inline int countMatchingPoles(axom::ArrayView<const Pole> poles,
                              const Pole& target,
                              int begin,
                              int end,
                              double tol)
{
  int count = 0;
  for(int idx = begin; idx < end; ++idx)
  {
    if(target.closeTo(poles[idx], tol))
    {
      ++count;
    }
  }
  return count;
}

/// \brief Collects one-based offsets of poles in a subrange that match a target pole.
inline axom::Array<int> collectMatchingPoleOffsets(axom::ArrayView<const Pole> poles,
                                                   const Pole& target,
                                                   int begin,
                                                   int end,
                                                   double tol)
{
  axom::Array<int> offsets;
  offsets.reserve(end - begin);

  for(int idx = begin; idx < end; ++idx)
  {
    if(target.closeTo(poles[idx], tol))
    {
      // Offsets are one-based relative positions in the coefficient block
      // beginning at `begin`; the Deckers moment data is written into those
      // future coefficient slots.
      offsets.push_back(idx - begin + 1);
    }
  }
  return offsets;
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

/// \brief Maps a rule from [-1,1] coordinates to [0,1] coordinates.
inline void map_rule_m11_to_unit_interval(axom::ArrayView<const double> nodes_m11,
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
