// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core/Array.hpp"
#include "axom/core/LRUCache.hpp"
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
#include <mutex>
#include <string>
#include <utility>

namespace axom
{
namespace numerics
{
namespace internal
{

inline double square(double value) { return value * value; }

[[noreturn]] void failRationalFejerPrecondition(const std::string& message)
{
  std::cerr << "ERROR: " << message << "\n";
  axom::utilities::processAbort();
  std::abort();
}

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

struct RationalFejerRuleStorage
{
  axom::Array<double> nodes;
  axom::Array<double> weights;
};

class Pole
{
  // Lightweight wrapper for poles in the reference [-1,1] domain. Keeping the
  // pole-specific tolerance checks, infinity normalization, and Cayley map here
  // avoids scattering those numerical conventions through the construction.
public:
  Pole() = default;
  Pole(double real, double imag) : m_value(real, imag) { }
  explicit Pole(const Complex& value) : m_value(value) { }

  // Finite poles with very large magnitude are treated as this infinity sentinel.
  static Pole infinity() { return Pole {std::numeric_limits<double>::infinity(), 0.0}; }

  double real() const { return m_value.real(); }
  double imag() const { return m_value.imag(); }
  const Complex& value() const { return m_value; }

  bool isInfinite() const
  {
    return std::isinf(m_value.real()) || std::isinf(m_value.imag()) || std::isinf(std::abs(m_value));
  }

  bool isEffectivelyReal(double tol) const { return std::abs(m_value.imag()) <= tol; }

  // Compare poles using relative distance when possible, with explicit handling
  // for zero and infinity so duplicate poles can be canonicalized robustly.
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

  // Rational Chebyshev formulas group conjugate pairs by upper-half-plane
  // representatives before they are expanded back into real-valued components.
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

struct RationalChebyshevPoleData
{
  axom::Array<Complex> cayley_poles;
  axom::Array<double> multiplicities;
  axom::Array<double> radii;
  axom::Array<double> phases;
  axom::Array<Complex> branch_points;
};

inline axom::Array<Pole> makePoleArray(axom::ArrayView<const Complex> pole_values)
{
  axom::Array<Pole> poles(pole_values.size(), pole_values.size());
  for(axom::IndexType i = 0; i < pole_values.size(); ++i)
  {
    poles[i] = Pole {pole_values[i]};
  }
  return poles;
}

inline axom::Array<Complex> mapUnitIntervalPolesToM11(axom::ArrayView<const Complex> poles01)
{
  // `_m11` denotes the `[-1,1]` interval ("minus-one-to-one"). Public callers
  // provide poles on `[0,1]`, then the implementation maps them into the
  // symmetric interval used by the rational Chebyshev and Fejer formulas.
  axom::Array<Complex> poles_m11(poles01.size(), poles01.size());
  for(axom::IndexType i = 0; i < poles01.size(); ++i)
  {
    const auto& pole = poles01[i];
    const auto pole_magnitude = std::abs(pole);
    poles_m11[i] = std::isinf(pole.real()) || std::isinf(pole.imag()) || std::isinf(pole_magnitude)
      ? Complex {axom::numeric_limits<double>::infinity(), 0.0}
      : 2.0 * pole - Complex {1.0, 0.0};
  }
  return poles_m11;
}

inline bool isInfinitePoleValue(const Complex& pole)
{
  return std::isinf(pole.real()) || std::isinf(pole.imag()) || std::isinf(std::abs(pole));
}

inline bool isInvalidFinitePoleValue(const Complex& pole)
{
  return !std::isfinite(pole.real()) || !std::isfinite(pole.imag());
}

inline bool isPoleOnInterval(const Complex& pole, double lower, double upper)
{
  const double tol = 16.0 * axom::numeric_limits<double>::epsilon();
  return std::abs(pole.imag()) <= tol && pole.real() >= lower - tol && pole.real() <= upper + tol;
}

inline void validatePoleSequence(axom::ArrayView<const Complex> poles,
                                 double lower,
                                 double upper,
                                 const char* domain_name)
{
  if(poles.empty())
  {
    failRationalFejerPrecondition(
      axom::fmt::format("Rational Fejer quadrature requires at least one pole in {}.", domain_name));
  }

  for(axom::IndexType i = 0; i < poles.size(); ++i)
  {
    const Complex pole = poles[i];
    if(std::isnan(pole.real()) || std::isnan(pole.imag()))
    {
      failRationalFejerPrecondition(
        axom::fmt::format("Rational Fejer pole {} in {} contains NaN.", i, domain_name));
    }

    if(isInfinitePoleValue(pole))
    {
      continue;
    }

    if(isInvalidFinitePoleValue(pole))
    {
      failRationalFejerPrecondition(
        axom::fmt::format("Rational Fejer pole {} in {} is neither finite nor infinite.",
                          i,
                          domain_name));
    }

    if(isPoleOnInterval(pole, lower, upper))
    {
      failRationalFejerPrecondition(
        axom::fmt::format("Rational Fejer finite pole {} = ({}, {}) lies on the {} interval.",
                          i,
                          pole.real(),
                          pole.imag(),
                          domain_name));
    }
  }
}

inline void validatePoleSequence(axom::ArrayView<const Pole> poles,
                                 double lower,
                                 double upper,
                                 const char* domain_name)
{
  axom::Array<Complex> pole_values(poles.size(), poles.size());
  for(axom::IndexType i = 0; i < poles.size(); ++i)
  {
    pole_values[i] = poles[i].value();
  }
  validatePoleSequence(pole_values, lower, upper, domain_name);
}

// Normalize the pole sequence by coalescing near-duplicate poles and keeping
// complex-conjugate partners adjacent. This becomes the stable internal
// representation used for caching, diagnostics, and the rational recurrences.
inline axom::Array<Pole> canonicalizePoleSequence(axom::Array<Pole> poles, double tol)
{
  const int num_input_poles = static_cast<int>(poles.size());
  axom::Array<int> consumed(num_input_poles, num_input_poles);
  consumed.fill(0);

  axom::Array<Pole> ordered_poles(2 * num_input_poles, 2 * num_input_poles);
  int ordered_count = 0;
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

    ordered_poles[ordered_count++] = pole;
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

      ordered_poles[ordered_count++] = conjugate_pole;
    }
  }

  axom::Array<Pole> result(ordered_count, ordered_count);
  for(int i = 0; i < ordered_count; ++i)
  {
    result[i] = ordered_poles[i];
  }

  return result;
}

// Collapse a sorted pole sequence into distinct representatives plus multiplicities,
// which is the form consumed by the rational Chebyshev and rational Fejer builders.
DistinctPoleData collectDistinctPoles(axom::ArrayView<const Pole> poles, double tol)
{
  axom::Array<Pole> sorted(poles);
  std::sort(sorted.begin(), sorted.end());

  DistinctPoleData distinct;
  // Collapse repeated or near-repeated poles into one representative value and
  // record the multiplicity that the rational basis construction should see.
  axom::IndexType num_distinct = 0;
  Pole representative;
  bool has_representative = false;
  for(const auto& pole : sorted)
  {
    if(!has_representative || !representative.closeTo(pole, tol))
    {
      representative = pole;
      has_representative = true;
      ++num_distinct;
    }
  }

  distinct.values = axom::Array<Pole>(num_distinct, num_distinct);
  distinct.multiplicities = axom::Array<int>(num_distinct, num_distinct);

  axom::IndexType distinct_index = -1;
  has_representative = false;
  for(const auto& pole : sorted)
  {
    if(!has_representative || !representative.closeTo(pole, tol))
    {
      representative = pole;
      has_representative = true;
      ++distinct_index;
      distinct.values[distinct_index] = pole;
      distinct.multiplicities[distinct_index] = 1;
    }
    else
    {
      ++distinct.multiplicities[distinct_index];
    }
  }

  return distinct;
}

// Standard interior Chebyshev first-kind rule on [-1,1]. In the rational machinery below,
// this is the all-infinite-pole limit of the rational Chebyshev stage.
inline void compute_chebyshev_first_kind_data_m11(int n,
                                                  axom::Array<double>& nodes,
                                                  axom::Array<double>& weights)
{
  nodes = axom::Array<double>(n, n);
  weights = axom::Array<double>(n, n);

  for(int i = 0; i < n; ++i)
  {
    const double theta = M_PI * (n - i - 0.5) / n;
    nodes[i] = std::cos(theta);
    weights[i] = M_PI / n;
  }
}

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

inline Complex expansion(const Complex& pole, int order, int parity, double tol)
{
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

// Exact pole-moment integrals used to seed the rational Fejer basis coefficients.
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
      J[multiplicity - 1] = expansion(pole_value, multiplicity, parity, tol);
      if(multiplicity > 1)
      {
        Complex F = J[multiplicity - 1];
        J[multiplicity - 2] = expansion(pole_value, multiplicity - 1, 1 - parity, tol);
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

class RationalChebyshevHelper
{
  // Constructs the rational Chebyshev rule on [-1,1] that underpins the
  // rational Fejer rule. This corresponds to the fixed-pole rational
  // Gauss-Chebyshev / near-best interpolation construction described in
  // van Deun, Deckers, Bultheel, and Weideman, ACM TOMS 35(2), 2008 (Algorithm 882).
public:
  explicit RationalChebyshevHelper(axom::Array<Pole> poles) : m_poles(std::move(poles))
  {
    for(auto& pole : m_poles)
    {
      pole = pole.withPositiveImaginaryMagnitude().normalizedInfinite(m_infinitePoleThreshold);
    }
  }

  void compute(axom::Array<double>& nodes, axom::Array<double>& weights) const
  {
    axom::Array<Pole> poles = m_poles;
    const int n = static_cast<int>(poles.size());

    const Pole last_pole = poles.back();
    double terminal_cayley_node = 0.0;
    if(last_pole.isInfinite())
    {
      terminal_cayley_node = 0.0;
    }
    else
    {
      terminal_cayley_node = last_pole.cayleyTransform().real();
      if(std::abs(terminal_cayley_node) <= 1.0 / (2.0 * m_infinitePoleThreshold))
      {
        terminal_cayley_node = 0.0;
        poles.back() = Pole::infinity();
      }
      else
      {
        // Algorithm 882 treats the last pole specially: after the Cayley map
        // its real value fixes one endpoint of the phase equation, while the
        // pole list itself is reflected back to the original x-domain.
        poles.back() = Pole {0.5 * (terminal_cayley_node + 1.0 / terminal_cayley_node), 0.0};
      }
    }

    if(n == 1)
    {
      nodes = axom::Array<double> {terminal_cayley_node};
      weights = axom::Array<double> {M_PI};
      return;
    }

    const auto distinct_poles = collectDistinctPoles(poles, m_poleTolerance);
    if(distinct_poles.values.size() == 1 && terminal_cayley_node == 0.0)
    {
      // The all-infinite-pole case reduces to the standard Chebyshev first-kind rule.
      compute_chebyshev_first_kind_data_m11(n, nodes, weights);
      return;
    }

    const auto cayley_pole_data = buildCayleyPoleData(poles);
    const auto residual_grid = buildThetaResidualGrid(n, cayley_pole_data);

    axom::Array<double> node_angles(n, n);
    int bracket_start = 0;
    for(int i = 0; i < n; ++i)
    {
      const double target = M_PI * (n - i - 0.5);
      bool ok = solveTheta(n, cayley_pole_data, residual_grid, target, bracket_start, node_angles[i]);
      if(!ok)
      {
        failRationalFejerPrecondition(
          axom::fmt::format("Failed to construct rational Chebyshev node {} of {}.", i, n));
      }
    }

    axom::Array<std::pair<double, double>> node_angle_pairs(n, n);
    for(int i = 0; i < n; ++i)
    {
      node_angle_pairs[i] = {std::cos(node_angles[i]), node_angles[i]};
    }
    std::sort(node_angle_pairs.begin(), node_angle_pairs.end(), [](const auto& lhs, const auto& rhs) {
      return lhs.first < rhs.first;
    });

    nodes = axom::Array<double>(n, n);
    node_angles = axom::Array<double>(n, n);
    for(int i = 0; i < n; ++i)
    {
      nodes[i] = node_angle_pairs[i].first;
      node_angles[i] = node_angle_pairs[i].second;
    }

    weights = computeWeights(node_angles, cayley_pole_data);
  }

private:
  struct ThetaResidualGrid
  {
    axom::Array<double> theta_samples;
    axom::Array<double> residual_samples;
  };

  RationalChebyshevPoleData buildCayleyPoleData(axom::ArrayView<const Pole> poles) const
  {
    // The rational Chebyshev phase and weight formulas use the cyclic pole
    // sequence (p_1, ..., p_n, p_1, ..., p_{n-1}). We collapse that doubled
    // list into distinct Cayley-domain poles and multiplicities before
    // evaluating the angle contributions.
    const axom::IndexType num_duplicated_poles = 2 * poles.size() - 1;
    axom::Array<Pole> duplicated_poles(num_duplicated_poles, num_duplicated_poles);
    for(axom::IndexType i = 0; i < poles.size(); ++i)
    {
      duplicated_poles[i] = poles[i];
    }
    for(axom::IndexType i = 0; i < poles.size() - 1; ++i)
    {
      duplicated_poles[poles.size() + i] = poles[i];
    }

    const auto distinct_poles = collectDistinctPoles(duplicated_poles, m_poleTolerance);

    const axom::IndexType num_distinct_poles = distinct_poles.values.size();

    RationalChebyshevPoleData data;
    data.cayley_poles = axom::Array<Complex>(num_distinct_poles, num_distinct_poles);
    data.multiplicities = axom::Array<double>(num_distinct_poles, num_distinct_poles);
    data.radii = axom::Array<double>(num_distinct_poles, num_distinct_poles);
    data.phases = axom::Array<double>(num_distinct_poles, num_distinct_poles);
    data.branch_points = axom::Array<Complex>(num_distinct_poles, num_distinct_poles);

    for(axom::IndexType i = 0; i < num_distinct_poles; ++i)
    {
      // The Cayley transform maps poles outside [-1,1] to the unit disk, where
      // the rational Chebyshev phase and weight formulas are expressed.
      const Complex cayley_pole = distinct_poles.values[i].isInfinite()
        ? Complex {0.0, 0.0}
        : distinct_poles.values[i].cayleyTransform();
      const double radius = std::abs(cayley_pole);
      const double phase = std::arg(cayley_pole);
      const Complex branch_point = (1.0 - radius) * std::exp(Complex {0.0, phase});

      data.cayley_poles[i] = cayley_pole;
      data.multiplicities[i] = 0.5 * distinct_poles.multiplicities[i];
      data.radii[i] = radius;
      data.phases[i] = phase;
      data.branch_points[i] = branch_point;
    }

    return data;
  }

  static double evaluatePoleAngleContribution(double theta,
                                              double phi,
                                              const Complex& branch_point,
                                              double multiplicity)
  {
    const double s = branch_point.imag();
    const double c = branch_point.real();
    const double sin1 = std::sin(0.5 * (theta + phi));
    const double sin2 = std::sin(0.5 * (theta - phi));
    const double y1 = 2.0 * sin2 * std::cos(0.5 * (theta + phi)) + s;
    const double y2 = 2.0 * sin1 * std::cos(0.5 * (theta - phi)) - s;
    const double x = -2.0 * sin1 * sin2 + c;
    const double branch = (x < 0.0 && y2 < 0.0) ? 2.0 * M_PI : 0.0;
    return multiplicity * (std::atan2(y1, x) + std::atan2(y2, x) + branch);
  }

  static bool bracketContainsRoot(double fa, double fb)
  {
    return fa == 0.0 || fb == 0.0 || (fa < 0.0 && fb > 0.0) || (fa > 0.0 && fb < 0.0);
  }

  double evaluateThetaResidual(double theta, int num_poles, const RationalChebyshevPoleData& data) const
  {
    // This residual is the phase equation whose roots define the rational Chebyshev node angles.
    // Each pole contributes an angle term after the Cayley transform.
    double value = -(num_poles - 1.0) * theta;
    for(axom::IndexType j = 0; j < data.cayley_poles.size(); ++j)
    {
      value += evaluatePoleAngleContribution(theta,
                                             data.phases[j],
                                             data.branch_points[j],
                                             data.multiplicities[j]);
    }
    return value;
  }

  ThetaResidualGrid buildThetaResidualGrid(int num_poles, const RationalChebyshevPoleData& data) const
  {
    // The residual contains atan2 branch changes, so a dense pre-scan gives the bisection
    // solve reliable brackets without assuming strict monotonicity at machine precision.
    const double eps_theta = 128.0 * axom::numeric_limits<double>::epsilon();
    const double theta_min = eps_theta;
    const double theta_max = M_PI - eps_theta;
    const int sample_count = axom::utilities::max(8192, 256 * num_poles);

    ThetaResidualGrid grid;
    grid.theta_samples = axom::Array<double>(sample_count + 1, sample_count + 1);
    grid.residual_samples = axom::Array<double>(sample_count + 1, sample_count + 1);
    for(int i = 0; i <= sample_count; ++i)
    {
      const double theta =
        theta_min + (theta_max - theta_min) * static_cast<double>(i) / sample_count;
      grid.theta_samples[i] = theta;
      grid.residual_samples[i] = evaluateThetaResidual(theta, num_poles, data);
    }
    return grid;
  }

  // Solve the rational Chebyshev phase condition for one node angle theta.
  // The corresponding quadrature node on [-1,1] is then x = cos(theta).
  bool solveTheta(int num_poles,
                  const RationalChebyshevPoleData& data,
                  const ThetaResidualGrid& grid,
                  double target,
                  int& bracket_start,
                  double& theta) const
  {
    const int sample_count = static_cast<int>(grid.theta_samples.size()) - 1;

    auto findBracket =
      [&](int start_index, double& lower, double& f_lower, double& upper, double& f_upper) {
        for(int i = start_index; i < sample_count; ++i)
        {
          lower = grid.theta_samples[i];
          upper = grid.theta_samples[i + 1];
          f_lower = grid.residual_samples[i] - target;
          f_upper = grid.residual_samples[i + 1] - target;
          if(bracketContainsRoot(f_lower, f_upper))
          {
            bracket_start = i;
            return true;
          }
        }
        return false;
      };

    double lower = 0.0;
    double f_lower = 0.0;
    double upper = 0.0;
    double f_upper = 0.0;
    bool found = findBracket(bracket_start, lower, f_lower, upper, f_upper);
    if(!found && bracket_start > 0)
    {
      found = findBracket(0, lower, f_lower, upper, f_upper);
    }
    if(!found)
    {
      return false;
    }

    for(int iter = 0; iter < 100; ++iter)
    {
      const double mid = 0.5 * (lower + upper);
      const double f_mid = evaluateThetaResidual(mid, num_poles, data) - target;
      if(std::abs(f_mid) <= 1e-14 || std::abs(upper - lower) <= 1e-14)
      {
        theta = mid;
        return true;
      }

      if(bracketContainsRoot(f_lower, f_mid))
      {
        upper = mid;
      }
      else
      {
        lower = mid;
        f_lower = f_mid;
      }
    }

    theta = 0.5 * (lower + upper);
    return true;
  }

  axom::Array<double> computeWeights(axom::ArrayView<const double> theta,
                                     const RationalChebyshevPoleData& data) const
  {
    const int n = static_cast<int>(theta.size());
    axom::Array<double> weights(n, n);
    weights.fill(1.0);
    for(int i = 0; i < n; ++i)
    {
      // Evaluate the closed-form rational Chebyshev weight denominator at the solved node angle.
      double denom = 1.0;
      for(axom::IndexType j = 0; j < data.cayley_poles.size(); ++j)
      {
        const double radius = data.radii[j];
        const double phase = data.phases[j];
        const double multiplicity = data.multiplicities[j];
        const double radius_denom = 1.0 - radius;
        const double left =
          radius_denom + 4.0 * radius * square(std::sin(0.5 * (theta[i] - phase))) / radius_denom;
        const double right =
          radius_denom + 4.0 * radius * square(std::sin(0.5 * (theta[i] + phase))) / radius_denom;
        denom += multiplicity * (1.0 + radius) * (1.0 / left + 1.0 / right);
      }

      weights[i] = 2.0 * M_PI / denom;
    }
    return weights;
  }

  axom::Array<Pole> m_poles;
  double m_infinitePoleThreshold {5.0 / axom::numeric_limits<double>::epsilon()};
  double m_poleTolerance {5.0 * axom::numeric_limits<double>::epsilon()};
};

// Build the rational Chebyshev rule used as the first stage of the extended rational
// Fejer construction. Adding a trailing pole at infinity produces the rational Chebyshev
// node set that the rational Fejer weight solve is built on.
inline void compute_rational_chebyshev_data(axom::ArrayView<const Pole> poles_in,
                                            axom::Array<double>& nodes,
                                            axom::Array<double>& weights)
{
  RationalChebyshevHelper helper {axom::Array<Pole>(poles_in)};
  helper.compute(nodes, weights);
}

inline void copy_array_to_array(axom::ArrayView<const double> values,
                                axom::Array<double>& array,
                                int allocatorID);

inline void copy_array_to_array(axom::ArrayView<const long double> values,
                                axom::Array<double>& array,
                                int allocatorID);

class RationalFejerBasis
{
  // Stores the rational basis sampled at the rational Chebyshev nodes so the
  // moment solve and final weight assembly can reuse the same matrix.
  // This corresponds to the orthonormal rational basis used in the
  // Deckers et al. extended rational Fejer construction; QuaHOG provided
  // a practical reference implementation for this Axom port.

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

  axom::Array<double> solveOrthogonalIntegrals(
    const Pole& pole,
    int multiplicity,
    int num_basis_columns,
    axom::ArrayView<const double> known_integrals,
    int component_count,
    RationalFejerDiagnostics::Step* step_diagnostics = nullptr,
    int diagnostic_allocator_id = axom::getDefaultAllocatorID()) const
  {
    const int n = static_cast<int>(m_nodes.size());
    axom::Array<long double> weighted_row0(n, n);
    weighted_row0.fill(0.0L);
    axom::Array<long double> weighted_row1(component_count == 2 ? n : 0,
                                           component_count == 2 ? n : 0);
    weighted_row1.fill(0.0L);

    const Complex pole_value = pole.value();
    // A real pole contributes one basis component. A non-real pole contributes
    // the real and negative-imaginary parts of the same rational factor, which
    // keeps the eventual quadrature weights real.
    for(int i = 0; i < n; ++i)
    {
      const Complex node_value {m_nodes[i], 0.0};
      const Complex basis_value = pole.isInfinite()
        ? powInteger(node_value, multiplicity)
        : powInteger((Complex {1.0, 0.0} - pole_value * m_nodes[i]) / (node_value - pole_value),
                     multiplicity);

      weighted_row0[i] =
        static_cast<long double>(m_lambda[i]) * static_cast<long double>(basis_value.real());
      if(component_count == 2)
      {
        weighted_row1[i] =
          static_cast<long double>(m_lambda[i]) * static_cast<long double>(-basis_value.imag());
      }
    }

    if(step_diagnostics != nullptr)
    {
      copy_array_to_array(weighted_row0, step_diagnostics->weighted_row0, diagnostic_allocator_id);

      if(component_count == 2)
      {
        copy_array_to_array(weighted_row1, step_diagnostics->weighted_row1, diagnostic_allocator_id);
      }
    }

    axom::Array<long double> bmat_row0(num_basis_columns, num_basis_columns);
    bmat_row0.fill(0.0L);
    axom::Array<long double> bmat_row1(component_count == 2 ? num_basis_columns : 0,
                                       component_count == 2 ? num_basis_columns : 0);
    bmat_row1.fill(0.0L);
    axom::Array<long double> bmat_row0_terms(num_basis_columns * n, num_basis_columns * n);
    bmat_row0_terms.fill(0.0L);
    axom::Array<long double> bmat_row1_terms(component_count == 2 ? num_basis_columns * n : 0,
                                             component_count == 2 ? num_basis_columns * n : 0);
    bmat_row1_terms.fill(0.0L);
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
        bmat_row0_terms[col * n + i] = term0;
        if(component_count == 2)
        {
          const long double term1 = weighted_row1[i] * basis_column_value;
          value1 += term1;
          bmat_row1_terms[col * n + i] = term1;
        }
      }
      bmat_row0[col] = value0;
      if(component_count == 2)
      {
        bmat_row1[col] = value1;
      }
    }

    long double rhs0 = static_cast<long double>(known_integrals[num_basis_columns - component_count]);
    for(int col = 0; col < num_basis_columns - component_count; ++col)
    {
      rhs0 -= bmat_row0[col] * static_cast<long double>(known_integrals[col]);
    }

    if(step_diagnostics != nullptr)
    {
      copy_array_to_array(bmat_row0, step_diagnostics->projected_row0, diagnostic_allocator_id);
      copy_array_to_array(bmat_row0_terms,
                          step_diagnostics->projected_row0_terms,
                          diagnostic_allocator_id);
      step_diagnostics->projected_row_terms_node_count = n;
      step_diagnostics->rhs0 = static_cast<double>(rhs0);
    }

    if(component_count == 1)
    {
      return axom::Array<double> {static_cast<double>(rhs0 / bmat_row0[num_basis_columns - 1])};
    }

    long double rhs1 = static_cast<long double>(known_integrals[num_basis_columns - 1]);
    for(int col = 0; col < num_basis_columns - component_count; ++col)
    {
      rhs1 -= bmat_row1[col] * static_cast<long double>(known_integrals[col]);
    }

    long double a00 = bmat_row0[num_basis_columns - 2];
    long double a01 = bmat_row0[num_basis_columns - 1];
    long double a10 = bmat_row1[num_basis_columns - 2];
    long double a11 = bmat_row1[num_basis_columns - 1];

    if(step_diagnostics != nullptr)
    {
      copy_array_to_array(bmat_row1, step_diagnostics->projected_row1, diagnostic_allocator_id);
      copy_array_to_array(bmat_row1_terms,
                          step_diagnostics->projected_row1_terms,
                          diagnostic_allocator_id);
      step_diagnostics->rhs1 = static_cast<double>(rhs1);
      step_diagnostics->a00 = static_cast<double>(a00);
      step_diagnostics->a01 = static_cast<double>(a01);
      step_diagnostics->a10 = static_cast<double>(a10);
      step_diagnostics->a11 = static_cast<double>(a11);
    }

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
    if(step_diagnostics != nullptr)
    {
      step_diagnostics->used_pivot = used_pivot;
      step_diagnostics->factor = static_cast<double>(factor);
      step_diagnostics->reduced_a11 = static_cast<double>(reduced_a11);
      step_diagnostics->reduced_rhs1 = static_cast<double>(reduced_rhs1);
      step_diagnostics->solved_x0 = static_cast<double>(x0);
      step_diagnostics->solved_x1 = static_cast<double>(x1);
    }
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

  axom::Array<double> assembleWeights(axom::ArrayView<const double> coefficients) const
  {
    // Convert basis coefficients back into point weights by evaluating the
    // sampled basis expansion and multiplying by the rational Chebyshev weights.
    axom::Array<double> weights = evaluateExpansion(coefficients);
    for(int i = 0; i < m_nodes.size(); ++i)
    {
      weights[i] *= m_lambda[i];
    }
    return weights;
  }

private:
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

inline void copy_array_to_array(axom::ArrayView<const double> values,
                                axom::Array<double>& array,
                                int allocatorID)
{
  array = axom::Array<double>(values, allocatorID);
}

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

// Construct the fixed-order rational Fejer nodes and weights on [-1,1]
// following Deckers, Mougaida, and Belhadjsalah, ACM TOMS 43(4), 2017
// (Algorithm 973). Axom uses the explicit node/weight construction here,
// not the semi-automatic adaptive integrator from that paper.
inline void copy_basis_matrix_to_array(const RationalFejerBasis& basis,
                                       axom::Array<double>& array,
                                       int allocatorID)
{
  copy_array_to_array(basis.basisColumns(), array, allocatorID);
}

inline int countRationalFejerSteps(axom::ArrayView<const Pole> canonical_poles, double pole_tolerance)
{
  int step_count = 0;
  int coeff_index = 1;
  const int num_poles = static_cast<int>(canonical_poles.size());
  while(coeff_index < num_poles + 1)
  {
    const Pole pole = canonical_poles[coeff_index - 1];
    const int component_count = pole.isEffectivelyReal(pole_tolerance) ? 1 : 2;
    ++step_count;
    coeff_index += component_count;
  }
  return step_count;
}

inline void compute_rational_fejer_data_m11_impl(axom::ArrayView<const Pole> poles,
                                                 axom::Array<double>& nodes,
                                                 axom::Array<double>& weights,
                                                 RationalFejerDiagnostics* diagnostics = nullptr,
                                                 int allocatorID = axom::getDefaultAllocatorID())
{
  validatePoleSequence(poles, -1.0, 1.0, "[-1,1]");

  const double infinite_pole_threshold = 2.0 / axom::numeric_limits<double>::epsilon();
  const double pole_tolerance = 2.0 * axom::numeric_limits<double>::epsilon();

  axom::Array<Pole> canonical_poles =
    canonicalizePoleSequence(axom::Array<Pole>(poles), pole_tolerance);
  for(auto& pole : canonical_poles)
  {
    pole = pole.normalizedInfinite(infinite_pole_threshold);
  }

  const int num_poles = static_cast<int>(canonical_poles.size());
  axom::Array<Complex> cayley_poles(num_poles, num_poles);
  for(int i = 0; i < num_poles; ++i)
  {
    // The rational basis is parameterized in the Cayley-transformed pole domain,
    // while the exact pole integrals are organized using the original pole ordering.
    const Pole mapped_pole = canonical_poles[i].withPositiveImaginaryMagnitude();
    cayley_poles[i] = mapped_pole.cayleyTransform();
    if(canonical_poles[i].imag() < 0.0)
    {
      cayley_poles[i] = std::conj(cayley_poles[i]);
    }
  }

  int diagnostic_step_index = 0;
  if(diagnostics != nullptr)
  {
    diagnostics->canonical_poles_m11 =
      axom::Array<Complex>(canonical_poles.size(), canonical_poles.size(), allocatorID);
    for(axom::IndexType i = 0; i < canonical_poles.size(); ++i)
    {
      diagnostics->canonical_poles_m11[i] = canonical_poles[i].value();
    }

    diagnostics->cayley_poles =
      axom::Array<Complex>(cayley_poles.size(), cayley_poles.size(), allocatorID);
    for(axom::IndexType i = 0; i < cayley_poles.size(); ++i)
    {
      diagnostics->cayley_poles[i] = cayley_poles[i];
    }

    const int diagnostic_step_count = countRationalFejerSteps(canonical_poles, pole_tolerance);
    diagnostics->steps = axom::Array<RationalFejerDiagnostics::Step>(diagnostic_step_count,
                                                                     diagnostic_step_count,
                                                                     allocatorID);
  }

  axom::Array<Pole> rational_chebyshev_poles(num_poles + 1, num_poles + 1);
  for(int i = 0; i < num_poles; ++i)
  {
    rational_chebyshev_poles[i] = canonical_poles[i];
  }
  // A rational Fejer rule for m finite/infinite poles uses m+1 rational Chebyshev nodes.
  // The appended infinity pole supplies that extra Chebyshev-like degree.
  rational_chebyshev_poles[num_poles] = Pole::infinity();

  axom::Array<double> rational_chebyshev_nodes;
  axom::Array<double> rational_chebyshev_weights;
  compute_rational_chebyshev_data(rational_chebyshev_poles,
                                  rational_chebyshev_nodes,
                                  rational_chebyshev_weights);

  if(diagnostics != nullptr)
  {
    copy_array_to_array(rational_chebyshev_nodes,
                        diagnostics->rational_chebyshev_nodes_m11,
                        allocatorID);
    copy_array_to_array(rational_chebyshev_weights,
                        diagnostics->rational_chebyshev_weights_m11,
                        allocatorID);
  }

  RationalFejerBasis basis(cayley_poles,
                           std::move(rational_chebyshev_nodes),
                           std::move(rational_chebyshev_weights));

  axom::Array<double> basis_coefficients(num_poles + 1, num_poles + 1);
  basis_coefficients.fill(0.0);
  basis_coefficients[0] = 2.0 / std::sqrt(M_PI);

  int coeff_index = 1;
  while(coeff_index < num_poles + 1)
  {
    const Pole pole = canonical_poles[coeff_index - 1];
    const int component_count = pole.isEffectivelyReal(pole_tolerance) ? 1 : 2;

    int pole_multiplicity_so_far = 0;
    for(int idx = 0; idx < coeff_index; ++idx)
    {
      if(pole.closeTo(canonical_poles[idx], pole_tolerance))
      {
        ++pole_multiplicity_so_far;
      }
    }

    if(pole_multiplicity_so_far == 1)
    {
      // Seed all occurrences of a newly encountered pole with the exact
      // non-orthogonal pole moments from the Deckers construction.
      int multiplicity = 0;
      for(int idx = coeff_index - 1; idx < num_poles; ++idx)
      {
        if(pole.closeTo(canonical_poles[idx], pole_tolerance))
        {
          ++multiplicity;
        }
      }

      axom::Array<int> repeated_pole_offsets(multiplicity, multiplicity);
      int repeated_index = 0;
      for(int idx = coeff_index - 1; idx < num_poles; ++idx)
      {
        if(pole.closeTo(canonical_poles[idx], pole_tolerance))
        {
          repeated_pole_offsets[repeated_index++] = idx - (coeff_index - 1) + 1;
        }
      }

      if(component_count == 2)
      {
        // Complex poles occupy adjacent coefficient slots in real/imaginary
        // order. Interleaving the matching conjugate offsets maps the exact
        // complex moments onto those real-valued slots.
        axom::Array<int> conjugate_positions(multiplicity, multiplicity);
        int conjugate_index = 0;
        for(int idx = coeff_index - 1; idx < num_poles; ++idx)
        {
          if(pole.conjugate().closeTo(canonical_poles[idx], pole_tolerance))
          {
            conjugate_positions[conjugate_index++] = idx - (coeff_index - 1) + 1;
          }
        }
        assert(conjugate_index == multiplicity);

        axom::Array<int> interleaved_positions(2 * multiplicity, 2 * multiplicity);
        for(int i = 0; i < multiplicity; ++i)
        {
          interleaved_positions[2 * i] = repeated_pole_offsets[i];
          interleaved_positions[2 * i + 1] = conjugate_positions[i];
        }
        repeated_pole_offsets = std::move(interleaved_positions);
      }

      const auto pole_integrals =
        computePoleMomentIntegrals(pole, multiplicity, component_count, pole_tolerance);
      for(axom::IndexType i = 0; i < repeated_pole_offsets.size(); ++i)
      {
        basis_coefficients[coeff_index + repeated_pole_offsets[i] - 1] = pole_integrals[i];
      }
    }

    const int num_known_columns = coeff_index + component_count;
    auto known_integrals = basis_coefficients.view().subspan(0, num_known_columns);

    RationalFejerDiagnostics::Step step_diagnostics;
    if(diagnostics != nullptr)
    {
      step_diagnostics.step_index = diagnostic_step_index;
      step_diagnostics.coeff_index_before = coeff_index;
      step_diagnostics.component_count = component_count;
      step_diagnostics.pole_multiplicity_so_far = pole_multiplicity_so_far;
      step_diagnostics.pole_m11 = pole.value();
      copy_array_to_array(known_integrals, step_diagnostics.basis_coefficients_before, allocatorID);
    }

    // Enforce orthogonality against the previously determined columns to recover
    // the next one or two rational basis coefficients.
    const auto orthogonal_integrals =
      basis.solveOrthogonalIntegrals(pole,
                                     pole_multiplicity_so_far,
                                     num_known_columns,
                                     known_integrals,
                                     component_count,
                                     diagnostics != nullptr ? &step_diagnostics : nullptr,
                                     allocatorID);
    for(int i = 0; i < component_count; ++i)
    {
      basis_coefficients[coeff_index + i] = orthogonal_integrals[i];
    }

    if(diagnostics != nullptr)
    {
      copy_array_to_array(orthogonal_integrals, step_diagnostics.orthogonal_integrals, allocatorID);
      copy_array_to_array(known_integrals, step_diagnostics.basis_coefficients_after, allocatorID);
      diagnostics->steps[diagnostic_step_index++] = std::move(step_diagnostics);
    }

    coeff_index += component_count;
  }

  nodes = basis.nodes();
  const axom::Array<double> basis_expansion = basis.evaluateExpansion(basis_coefficients);
  weights = basis_expansion;
  const auto& lambda = basis.lambda();
  for(int i = 0; i < weights.size(); ++i)
  {
    weights[i] *= lambda[i];
  }

  if(diagnostics != nullptr)
  {
    copy_array_to_array(basis_coefficients, diagnostics->basis_coefficients, allocatorID);
    copy_array_to_array(basis_expansion, diagnostics->basis_expansion_m11, allocatorID);
    copy_array_to_array(weights, diagnostics->final_weights_m11, allocatorID);
    copy_basis_matrix_to_array(basis, diagnostics->basis_matrix_transpose_m11, allocatorID);
    diagnostics->basis_matrix_row_count = basis.basisColumnCount();
    diagnostics->basis_matrix_col_count = basis.nodeCount();
  }
}

inline void compute_rational_fejer_data_m11(axom::ArrayView<const Pole> poles,
                                            axom::Array<double>& nodes,
                                            axom::Array<double>& weights)
{
  compute_rational_fejer_data_m11_impl(poles, nodes, weights);
}
inline void compute_rational_fejer_diagnostics_m11(axom::ArrayView<const Pole> poles,
                                                   RationalFejerDiagnostics& diagnostics,
                                                   int allocatorID)
{
  axom::Array<double> nodes_m11;
  axom::Array<double> weights_m11;
  compute_rational_fejer_data_m11_impl(poles, nodes_m11, weights_m11, &diagnostics, allocatorID);
}

void compute_rational_chebyshev_data_m11(axom::ArrayView<const Complex> poles_m11,
                                         axom::Array<double>& nodes,
                                         axom::Array<double>& weights,
                                         int allocatorID)
{
  validatePoleSequence(poles_m11, -1.0, 1.0, "[-1,1]");

  axom::Array<double> nodes_tmp;
  axom::Array<double> weights_tmp;
  compute_rational_chebyshev_data(makePoleArray(poles_m11), nodes_tmp, weights_tmp);
  copy_array_to_array(nodes_tmp, nodes, allocatorID);
  copy_array_to_array(weights_tmp, weights, allocatorID);
}

void compute_rational_fejer_data_m11(axom::ArrayView<const Complex> poles_m11,
                                     axom::Array<double>& nodes,
                                     axom::Array<double>& weights,
                                     int allocatorID)
{
  validatePoleSequence(poles_m11, -1.0, 1.0, "[-1,1]");

  axom::Array<double> nodes_tmp;
  axom::Array<double> weights_tmp;
  compute_rational_fejer_data_m11(makePoleArray(poles_m11), nodes_tmp, weights_tmp);
  copy_array_to_array(nodes_tmp, nodes, allocatorID);
  copy_array_to_array(weights_tmp, weights, allocatorID);
}

inline void compute_rational_fejer_data_m11_impl(axom::ArrayView<const Complex> poles_m11,
                                                 axom::Array<double>& nodes,
                                                 axom::Array<double>& weights,
                                                 RationalFejerDiagnostics* diagnostics,
                                                 int allocatorID)
{
  compute_rational_fejer_data_m11_impl(makePoleArray(poles_m11),
                                       nodes,
                                       weights,
                                       diagnostics,
                                       allocatorID);
}

void compute_rational_fejer_diagnostics_m11(axom::ArrayView<const Complex> poles_m11,
                                            RationalFejerDiagnostics& diagnostics,
                                            int allocatorID)
{
  validatePoleSequence(poles_m11, -1.0, 1.0, "[-1,1]");

  compute_rational_fejer_diagnostics_m11(makePoleArray(poles_m11), diagnostics, allocatorID);
}
inline std::string make_rational_fejer_key(axom::ArrayView<const Pole> poles01, int allocatorID)
{
  // Cache by the canonicalized pole sequence that the rule construction
  // actually uses. This preserves correctness while improving reuse when two
  // numerically equivalent pole lists differ only by tiny root-solver noise.
  const double infinite_pole_threshold = 2.0 / axom::numeric_limits<double>::epsilon();
  const double pole_tolerance = 2.0 * axom::numeric_limits<double>::epsilon();

  axom::Array<Pole> canonical_poles =
    canonicalizePoleSequence(axom::Array<Pole>(poles01), pole_tolerance);
  for(auto& pole : canonical_poles)
  {
    pole = pole.normalizedInfinite(infinite_pole_threshold);
  }

  axom::fmt::memory_buffer key;
  axom::fmt::format_to(std::back_inserter(key), "{}|", allocatorID);
  for(const auto& pole : canonical_poles)
  {
    axom::fmt::format_to(std::back_inserter(key), "{:a},{:a};", pole.real(), pole.imag());
  }
  return axom::fmt::to_string(key);
}

std::string make_rational_fejer_key(axom::ArrayView<const Complex> poles01, int allocatorID)
{
  return make_rational_fejer_key(makePoleArray(poles01), allocatorID);
}

}  // namespace internal

void compute_rational_fejer_data(axom::ArrayView<const std::complex<double>> poles01,
                                 axom::Array<double>& nodes,
                                 axom::Array<double>& weights,
                                 int allocatorID)
{
  internal::validatePoleSequence(poles01, 0.0, 1.0, "[0,1]");

  const auto poles_m11 = internal::mapUnitIntervalPolesToM11(poles01);

  axom::Array<double> nodes_m11;
  axom::Array<double> weights_m11;

  // The implementation below follows the Deckers et al. / van Deun et al.
  // rational Fejer and rational Chebyshev constructions. QuaHOG's MATLAB
  // Spectral_PE implementation served as a secondary reference for how these
  // pieces are assembled in the planar integration workflow.
  internal::compute_rational_fejer_data_m11_impl(poles_m11, nodes_m11, weights_m11, nullptr, allocatorID);

  internal::map_rule_m11_to_unit_interval(nodes_m11, weights_m11, nodes, weights, allocatorID);
}

void internal::compute_rational_fejer_diagnostics(axom::ArrayView<const std::complex<double>> poles01,
                                                  internal::RationalFejerDiagnostics& diagnostics,
                                                  int allocatorID)
{
  internal::validatePoleSequence(poles01, 0.0, 1.0, "[0,1]");

  const auto poles_m11 = internal::mapUnitIntervalPolesToM11(poles01);

  axom::Array<double> nodes_m11;
  axom::Array<double> weights_m11;
  internal::compute_rational_fejer_data_m11_impl(poles_m11,
                                                 nodes_m11,
                                                 weights_m11,
                                                 &diagnostics,
                                                 allocatorID);
  internal::map_rule_m11_to_unit_interval(nodes_m11,
                                          weights_m11,
                                          diagnostics.nodes_01,
                                          diagnostics.weights_01,
                                          allocatorID);
}

QuadratureRule get_rational_fejer(axom::ArrayView<const std::complex<double>> poles01, int allocatorID)
{
  internal::validatePoleSequence(poles01, 0.0, 1.0, "[0,1]");

  constexpr std::size_t MAX_RATIONAL_FEJER_CACHED_RULES = 1u << 16;
  static axom::LRUCache<std::string, internal::RationalFejerRuleStorage> rule_library(
    MAX_RATIONAL_FEJER_CACHED_RULES);
  static std::mutex rule_library_mutex;

  const std::string key = internal::make_rational_fejer_key(poles01, allocatorID);

  {
    const std::lock_guard<std::mutex> lock(rule_library_mutex);
    if(auto* storage = rule_library.find(key))
    {
      return QuadratureRule {storage->nodes.view(), storage->weights.view()};
    }
  }

  internal::RationalFejerRuleStorage storage;
  compute_rational_fejer_data(poles01, storage.nodes, storage.weights, allocatorID);

  {
    const std::lock_guard<std::mutex> lock(rule_library_mutex);
    if(auto* cached_storage = rule_library.find(key))
    {
      return QuadratureRule {cached_storage->nodes.view(), cached_storage->weights.view()};
    }

    auto& cached_storage = rule_library.insert(key, std::move(storage));
    return QuadratureRule {cached_storage.nodes.view(), cached_storage.weights.view()};
  }
}

}  // namespace numerics
}  // namespace axom
