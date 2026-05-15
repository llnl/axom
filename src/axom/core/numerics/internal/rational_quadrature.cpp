// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

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
#include <cstdint>
#include <limits>
#include <map>
#include <mutex>
#include <string>
#include <vector>

namespace axom
{
namespace numerics
{
namespace internal
{


inline double square(double value) { return value * value; }

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
public:
  Pole() = default;
  Pole(double real, double imag)
    : m_value(real, imag)
  { }
  explicit Pole(const Complex& value)
    : m_value(value)
  { }

  static Pole infinity() { return Pole {std::numeric_limits<double>::infinity(), 0.0}; }

  double real() const { return m_value.real(); }
  double imag() const { return m_value.imag(); }
  const Complex& value() const { return m_value; }

  bool isInfinite() const
  {
    return std::isinf(m_value.real()) || std::isinf(m_value.imag()) ||
      std::isinf(std::abs(m_value));
  }

  bool isEffectivelyReal(double tol) const { return std::abs(m_value.imag()) <= tol; }

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

  Pole withPositiveImaginaryMagnitude() const
  {
    return Pole {Complex {m_value.real(), std::abs(m_value.imag())}};
  }

  Pole conjugate() const { return Pole {std::conj(m_value)}; }

  Complex cayleyTransform() const
  {
    if(isInfinite())
    {
      return Complex {0.0, 0.0};
    }

    const double sign = m_value.real() >= 0.0 ? 1.0 : -1.0;
    return Complex {1.0, 0.0} /
      (m_value + sign * std::sqrt(m_value * m_value - Complex {1.0, 0.0}));
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
  std::vector<Pole> values;
  std::vector<int> multiplicities;
};

struct RChebPoleData
{
  std::vector<Complex> cayley_poles;
  std::vector<double> multiplicities;
  std::vector<double> radii;
  std::vector<double> phases;
  std::vector<Complex> branch_points;
};

// Normalize the pole sequence by coalescing near-duplicate poles and keeping
// complex-conjugate partners adjacent. This becomes the stable internal
// representation used for caching, diagnostics, and the rational recurrences.
inline std::vector<Pole> canonicalizePoleSequence(std::vector<Pole> poles, double tol)
{
  int j = 0;
  while(j < static_cast<int>(poles.size()))
  {
    for(int idx = j + 1; idx < static_cast<int>(poles.size()); ++idx)
    {
      if(poles[j].closeTo(poles[idx], tol))
      {
        poles[idx] = poles[j];
        break;
      }
    }

    ++j;
    if(j == 0 || j - 1 >= static_cast<int>(poles.size()))
    {
      continue;
    }

    if(!poles[j - 1].isEffectivelyReal(tol))
    {
      const auto conjugate_pole = poles[j - 1].conjugate();
      for(int idx = j; idx < static_cast<int>(poles.size()); ++idx)
      {
        if(conjugate_pole.closeTo(poles[idx], tol))
        {
          poles.erase(poles.begin() + idx);
          break;
        }
      }

      poles.insert(poles.begin() + j, conjugate_pole);
      ++j;
    }
  }

  return poles;
}

// Collapse a sorted pole sequence into distinct representatives plus
// multiplicities, which is the form consumed by the rcheb/rfejer builders.
DistinctPoleData collectDistinctPoles(const std::vector<Pole>& poles, double tol)
{
  std::vector<Pole> sorted = poles;
  std::sort(sorted.begin(), sorted.end());

  DistinctPoleData distinct;
  // Collapse repeated or near-repeated poles into one representative value and
  // record the multiplicity that the rational basis construction should see.
  for(const auto& pole : sorted)
  {
    if(distinct.values.empty() || !distinct.values.back().closeTo(pole, tol))
    {
      distinct.values.push_back(pole);
      distinct.multiplicities.push_back(1);
    }
    else
    {
      ++distinct.multiplicities.back();
    }
  }

  return distinct;
}

// Standard interior Chebyshev first-kind rule on [-1,1]. In the rational
// machinery below, this is the all-infinite-pole limit of the rcheb stage.
inline void compute_chebyshev_first_kind_data_m11(int n,
                                           std::vector<double>& nodes,
                                           std::vector<double>& weights)
{
  nodes.resize(n);
  weights.resize(n);

  for(int i = 0; i < n; ++i)
  {
    const double theta = M_PI * (n - i - 0.5) / n;
    nodes[i] = std::cos(theta);
    weights[i] = M_PI / n;
  }
}

inline void map_rule_m11_to_unit_interval(const std::vector<double>& nodes_m11,
                                   const std::vector<double>& weights_m11,
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
                 (2 * j + order + 1 - parity) * tol) <
        0.0)
  {
    ++j;
  }

  Complex series_sum {0.0, 0.0};
  for(int exponent = 0; exponent <= 2 * j; exponent += 2)
  {
    const double d1 = exponent + order + 3 - parity;
    const double d2 = exponent + order + 1 - parity;
    series_sum += Complex {1.0, 0.0} /
      (d1 * d2 * std::pow(pole, exponent));
  }

  return 2.0 * (1 - parity) / (order + 1.0) + std::pow(-1.0, parity) * 4.0 * order * series_sum /
    std::pow(pole, 2 - parity);
}

// Exact pole-moment integrals used to seed the rational Fejer basis coefficients.
inline std::vector<double> computePoleMomentIntegrals(const Pole& pole,
                                             int multiplicity,
                                             int component_count,
                                             double tol)
{
  // These are the non-orthogonal moments for one repeated pole. The Deckers
  // construction uses them as exact seed data before the basis coefficients
  // are orthogonalized against the sampled rational basis.
  std::vector<Complex> J(multiplicity, Complex {0.0, 0.0});

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
          const Complex next = -((Complex {1.0, 0.0} / (pole_value * pole_value)) - Complex {1.0, 0.0}) *
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
    std::vector<double> result(multiplicity);
    for(int i = 0; i < multiplicity; ++i)
    {
      result[i] = J[i].real();
    }
    return result;
  }

  std::vector<double> result(2 * multiplicity);
  for(int i = 0; i < multiplicity; ++i)
  {
    result[2 * i] = J[i].real();
    result[2 * i + 1] = -J[i].imag();
  }
  return result;
}
class RChebHelper
{
  // Constructs the rational Chebyshev rule on [-1,1] that underpins the
  // rational Fejer rule. This corresponds to the fixed-pole rational
  // Gauss-Chebyshev / near-best interpolation construction described in
  // van Deun, Deckers, Bultheel, and Weideman, ACM TOMS 35(2), 2008
  // (Algorithm 882).

public:
  explicit RChebHelper(std::vector<Pole> poles)
    : m_poles(std::move(poles))
  {
    for(auto& pole : m_poles)
    {
      pole = pole.withPositiveImaginaryMagnitude().normalizedInfinite(m_infinitePoleThreshold);
    }
  }

  void compute(std::vector<double>& nodes, std::vector<double>& weights) const
  {
    std::vector<Pole> poles = m_poles;
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
        poles.back() = Pole {0.5 * (terminal_cayley_node + 1.0 / terminal_cayley_node), 0.0};
      }
    }

    if(n == 1)
    {
      nodes = {terminal_cayley_node};
      weights = {M_PI};
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

    std::vector<double> node_angles(n);
    int bracket_start = 0;
    for(int i = 0; i < n; ++i)
    {
      const double target = M_PI * (n - i - 0.5);
      bool ok =
        solveTheta(n, cayley_pole_data, residual_grid, target, bracket_start, node_angles[i]);
      assert("Failed to construct rational Chebyshev nodes" && ok);
    }

    std::vector<std::pair<double, double>> node_angle_pairs;
    node_angle_pairs.reserve(n);
    for(int i = 0; i < n; ++i)
    {
      node_angle_pairs.push_back({std::cos(node_angles[i]), node_angles[i]});
    }
    std::sort(node_angle_pairs.begin(), node_angle_pairs.end(), [](const auto& lhs, const auto& rhs) {
      return lhs.first < rhs.first;
    });

    nodes.resize(n);
    node_angles.resize(n);
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
    std::vector<double> theta_samples;
    std::vector<double> residual_samples;
  };

  RChebPoleData buildCayleyPoleData(const std::vector<Pole>& poles) const
  {
    std::vector<Pole> duplicated_poles = poles;
    duplicated_poles.insert(duplicated_poles.end(), poles.begin(), poles.end() - 1);

    const auto distinct_poles = collectDistinctPoles(duplicated_poles, m_poleTolerance);

    RChebPoleData data;
    data.cayley_poles.reserve(distinct_poles.values.size());
    data.multiplicities.reserve(distinct_poles.values.size());
    data.radii.reserve(distinct_poles.values.size());
    data.phases.reserve(distinct_poles.values.size());
    data.branch_points.reserve(distinct_poles.values.size());

    for(std::size_t i = 0; i < distinct_poles.values.size(); ++i)
    {
      // The Cayley transform maps poles outside [-1,1] to the unit disk, where
      // the rational Chebyshev phase and weight formulas are expressed.
      const Complex cayley_pole = distinct_poles.values[i].isInfinite()
        ? Complex {0.0, 0.0}
        : distinct_poles.values[i].cayleyTransform();
      const double radius = std::abs(cayley_pole);
      const double phase = std::arg(cayley_pole);
      const Complex branch_point = (1.0 - radius) * std::exp(Complex {0.0, phase});

      data.cayley_poles.push_back(cayley_pole);
      data.multiplicities.push_back(0.5 * distinct_poles.multiplicities[i]);
      data.radii.push_back(radius);
      data.phases.push_back(phase);
      data.branch_points.push_back(branch_point);
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

  double evaluateThetaResidual(double theta, int num_poles, const RChebPoleData& data) const
  {
    // This residual is the phase equation whose roots define the rcheb node
    // angles. Each pole contributes an angle term after the Cayley transform.
    double value = -(num_poles - 1.0) * theta;
    for(std::size_t j = 0; j < data.cayley_poles.size(); ++j)
    {
      value += evaluatePoleAngleContribution(
        theta, data.phases[j], data.branch_points[j], data.multiplicities[j]);
    }
    return value;
  }

  ThetaResidualGrid buildThetaResidualGrid(int num_poles, const RChebPoleData& data) const
  {
    const double eps_theta = 128.0 * axom::numeric_limits<double>::epsilon();
    const double theta_min = eps_theta;
    const double theta_max = M_PI - eps_theta;
    const int sample_count = axom::utilities::max(8192, 256 * num_poles);

    ThetaResidualGrid grid;
    grid.theta_samples.resize(sample_count + 1);
    grid.residual_samples.resize(sample_count + 1);
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
                  const RChebPoleData& data,
                  const ThetaResidualGrid& grid,
                  double target,
                  int& bracket_start,
                  double& theta) const
  {
    const int sample_count = static_cast<int>(grid.theta_samples.size()) - 1;

    auto findBracket = [&](int start_index,
                           double& lower,
                           double& f_lower,
                           double& upper,
                           double& f_upper) {
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

  std::vector<double> computeWeights(const std::vector<double>& theta,
                                     const RChebPoleData& data) const
  {
    const int n = static_cast<int>(theta.size());
    std::vector<double> weights(n, 1.0);
    for(int i = 0; i < n; ++i)
    {
      // Evaluate the closed-form rcheb weight denominator at the solved node angle.
      double denom = 1.0;
      for(std::size_t j = 0; j < data.cayley_poles.size(); ++j)
      {
        const double radius = data.radii[j];
        const double phase = data.phases[j];
        const double multiplicity = data.multiplicities[j];
        const double radius_denom = 1.0 - radius;
        const double left = radius_denom +
          4.0 * radius * square(std::sin(0.5 * (theta[i] - phase))) / radius_denom;
        const double right = radius_denom +
          4.0 * radius * square(std::sin(0.5 * (theta[i] + phase))) / radius_denom;
        denom += multiplicity * (1.0 + radius) * (1.0 / left + 1.0 / right);
      }

      weights[i] = 2.0 * M_PI / denom;
    }
    return weights;
  }

  std::vector<Pole> m_poles;
  double m_infinitePoleThreshold {5.0 / axom::numeric_limits<double>::epsilon()};
  double m_poleTolerance {5.0 * axom::numeric_limits<double>::epsilon()};
};

// Build the rational Chebyshev rule used as the first stage of the
// extended rational Fejer construction. Adding a trailing pole at infinity
// produces the rcheb node set that the rational Fejer weight solve is built on.
inline void compute_rcheb_data(const std::vector<Pole>& poles_in,
                        std::vector<double>& nodes,
                        std::vector<double>& weights)
{
  RChebHelper helper(poles_in);
  helper.compute(nodes, weights);
}

inline void copy_vector_to_array(const std::vector<double>& values,
                          axom::Array<double>& array,
                          int allocatorID);

class RationalFejerBasis
{
  // Stores the rational basis sampled at the rcheb nodes so the
  // moment solve and final weight assembly can reuse the same matrix.
  // This corresponds to the orthonormal rational basis used in the
  // Deckers et al. extended rational Fejer construction; QuaHOG provided
  // a practical reference implementation for this Axom port.

public:
  RationalFejerBasis(const std::vector<Complex>& cayley_poles,
                     std::vector<double> nodes,
                     std::vector<double> lambda)
    : m_nodes(std::move(nodes))
    , m_lambda(std::move(lambda))
    , m_QColumns(buildBasisColumns(cayley_poles, m_nodes))
  { }

  const std::vector<double>& nodes() const { return m_nodes; }

  const std::vector<double>& lambda() const { return m_lambda; }

  const std::vector<std::vector<double>>& basisColumns() const { return m_QColumns; }

  std::vector<double> solveOrthogonalIntegrals(
    const Pole& pole,
    int multiplicity,
    int num_basis_columns,
    const std::vector<double>& known_integrals,
    int component_count,
    RationalFejerDiagnostics::Step* step_diagnostics = nullptr) const
  {
    const int n = static_cast<int>(m_nodes.size());
    std::vector<long double> weighted_row0(n, 0.0L);
    std::vector<long double> weighted_row1(component_count == 2 ? n : 0, 0.0L);

    const Complex pole_value = pole.value();
    for(int i = 0; i < n; ++i)
    {
      const Complex node_value {m_nodes[i], 0.0};
      const Complex basis_value = pole.isInfinite()
        ? powInteger(node_value, multiplicity)
        : powInteger((Complex {1.0, 0.0} - pole_value * m_nodes[i]) /
                       (node_value - pole_value),
                     multiplicity);

      weighted_row0[i] = static_cast<long double>(m_lambda[i]) *
        static_cast<long double>(basis_value.real());
      if(component_count == 2)
      {
        weighted_row1[i] = static_cast<long double>(m_lambda[i]) *
          static_cast<long double>(-basis_value.imag());
      }
    }

    if(step_diagnostics != nullptr)
    {
      std::vector<double> weighted_row0_doubles(n);
      for(int i = 0; i < n; ++i)
      {
        weighted_row0_doubles[i] = static_cast<double>(weighted_row0[i]);
      }
      copy_vector_to_array(weighted_row0_doubles,
                           step_diagnostics->weighted_row0,
                           axom::getDefaultAllocatorID());

      if(component_count == 2)
      {
        std::vector<double> weighted_row1_doubles(n);
        for(int i = 0; i < n; ++i)
        {
          weighted_row1_doubles[i] = static_cast<double>(weighted_row1[i]);
        }
        copy_vector_to_array(weighted_row1_doubles,
                             step_diagnostics->weighted_row1,
                             axom::getDefaultAllocatorID());
      }
    }

    std::vector<long double> bmat_row0(num_basis_columns, 0.0L);
    std::vector<long double> bmat_row1(component_count == 2 ? num_basis_columns : 0, 0.0L);
    std::vector<long double> bmat_row0_terms(num_basis_columns * n, 0.0L);
    std::vector<long double> bmat_row1_terms(component_count == 2 ? num_basis_columns * n : 0, 0.0L);
    // Project the pole-dependent function against the sampled orthonormal basis.
    // The trailing 1x1 or 2x2 block is then solved for the next basis coefficients.
    for(int col = 0; col < num_basis_columns; ++col)
    {
      const auto& basis_column = m_QColumns[col];
      long double value0 = 0.0L;
      long double value1 = 0.0L;
      for(int i = 0; i < n; ++i)
      {
        const long double basis_column_value = static_cast<long double>(basis_column[i]);
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

    long double rhs0 =
      static_cast<long double>(known_integrals[num_basis_columns - component_count]);
    for(int col = 0; col < num_basis_columns - component_count; ++col)
    {
      rhs0 -= bmat_row0[col] * static_cast<long double>(known_integrals[col]);
    }

    if(step_diagnostics != nullptr)
    {
      std::vector<double> projected_row0(num_basis_columns);
      for(int col = 0; col < num_basis_columns; ++col)
      {
        projected_row0[col] = static_cast<double>(bmat_row0[col]);
      }
      copy_vector_to_array(projected_row0, step_diagnostics->projected_row0, axom::getDefaultAllocatorID());

      std::vector<double> projected_row0_terms_doubles(num_basis_columns * n);
      for(int index = 0; index < num_basis_columns * n; ++index)
      {
        projected_row0_terms_doubles[index] = static_cast<double>(bmat_row0_terms[index]);
      }
      copy_vector_to_array(projected_row0_terms_doubles,
                           step_diagnostics->projected_row0_terms,
                           axom::getDefaultAllocatorID());
      step_diagnostics->projected_row_terms_node_count = n;
      step_diagnostics->rhs0 = static_cast<double>(rhs0);
    }

    if(component_count == 1)
    {
      return {static_cast<double>(rhs0 / bmat_row0[num_basis_columns - 1])};
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
      std::vector<double> projected_row1(num_basis_columns);
      for(int col = 0; col < num_basis_columns; ++col)
      {
        projected_row1[col] = static_cast<double>(bmat_row1[col]);
      }
      copy_vector_to_array(projected_row1, step_diagnostics->projected_row1, axom::getDefaultAllocatorID());

      std::vector<double> projected_row1_terms_doubles(num_basis_columns * n);
      for(int index = 0; index < num_basis_columns * n; ++index)
      {
        projected_row1_terms_doubles[index] = static_cast<double>(bmat_row1_terms[index]);
      }
      copy_vector_to_array(projected_row1_terms_doubles,
                           step_diagnostics->projected_row1_terms,
                           axom::getDefaultAllocatorID());
      step_diagnostics->rhs1 = static_cast<double>(rhs1);
      step_diagnostics->a00 = static_cast<double>(a00);
      step_diagnostics->a01 = static_cast<double>(a01);
      step_diagnostics->a10 = static_cast<double>(a10);
      step_diagnostics->a11 = static_cast<double>(a11);
    }

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
    return {static_cast<double>(x0), static_cast<double>(x1)};
  }

  std::vector<double> evaluateExpansion(const std::vector<double>& coefficients) const
  {
    std::vector<double> values(m_nodes.size(), 0.0);
    for(std::size_t j = 0; j < coefficients.size(); ++j)
    {
      const double coefficient = coefficients[j];
      if(coefficient == 0.0)
      {
        continue;
      }

      const auto& basis_column = m_QColumns[j];
      for(std::size_t i = 0; i < m_nodes.size(); ++i)
      {
        values[i] += basis_column[i] * coefficient;
      }
    }

    return values;
  }

  std::vector<double> assembleWeights(const std::vector<double>& coefficients) const
  {
    // Convert basis coefficients back into point weights by evaluating the
    // sampled basis expansion and multiplying by the rcheb weights.
    std::vector<double> weights = evaluateExpansion(coefficients);
    for(std::size_t i = 0; i < m_nodes.size(); ++i)
    {
      weights[i] *= m_lambda[i];
    }
    return weights;
  }

private:
  static std::vector<std::vector<double>> buildBasisColumns(const std::vector<Complex>& cayley_poles,
                                                            const std::vector<double>& nodes)
  {
    // Sample the orthonormal rational basis at every rcheb node. Later stages
    // reuse this matrix both when enforcing orthogonality and when rebuilding
    // the final quadrature weights from the recovered coefficients.
    const int n = static_cast<int>(nodes.size());
    std::vector<std::vector<double>> Q_columns(n, std::vector<double>(n, 0.0));

    for(int i = 0; i < n; ++i)
    {
      Q_columns[0][i] = 1.0 / std::sqrt(M_PI);
    }

    if(n == 1)
    {
      return Q_columns;
    }

    std::vector<Complex> z_values(n);
    for(int i = 0; i < n; ++i)
    {
      z_values[i] = Complex {nodes[i],
                             std::sqrt(axom::utilities::max(0.0, 1.0 - nodes[i] * nodes[i]))};
    }

    std::vector<Complex> blaschke_products(n, Complex {1.0, 0.0});
    std::vector<Complex> conjugate_blaschke_products(n, Complex {1.0, 0.0});
    const double normalizer = std::sqrt(1.0 / (2.0 * M_PI));

    for(int i = 0; i < n; ++i)
    {
      const Complex z = z_values[i];
      const double basis_normalizer =
        normalizer * std::sqrt(axom::utilities::max(0.0, 1.0 - std::norm(cayley_poles[0])));
      const Complex value =
        basis_normalizer *
        (z * conjugate_blaschke_products[i] / (Complex {1.0, 0.0} - cayley_poles[0] * z) +
         Complex {1.0, 0.0} / ((z - cayley_poles[0]) * blaschke_products[i]));
      Q_columns[1][i] = value.real();
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
        blaschke_products[i] *=
          (z - pole) / (Complex {1.0, 0.0} - conjugate_pole * z);
        conjugate_blaschke_products[i] *=
          (z - conjugate_pole) / (Complex {1.0, 0.0} - pole * z);

        const Complex value =
          basis_normalizer *
          (z * conjugate_blaschke_products[i] / (Complex {1.0, 0.0} - cayley_poles[k] * z) +
           Complex {1.0, 0.0} / ((z - cayley_poles[k]) * blaschke_products[i]));
        Q_columns[k + 1][i] = value.real();
      }
    }

    return Q_columns;
  }

  std::vector<double> m_nodes;
  std::vector<double> m_lambda;
  std::vector<std::vector<double>> m_QColumns;
};

inline void copy_vector_to_array(const std::vector<double>& values,
                          axom::Array<double>& array,
                          int allocatorID)
{
  const int size = static_cast<int>(values.size());
  array = axom::Array<double>(size, size, allocatorID);
  for(int i = 0; i < size; ++i)
  {
    array[i] = values[i];
  }
}

// Construct the fixed-order rational Fejer nodes and weights on [-1,1]
// following Deckers, Mougaida, and Belhadjsalah, ACM TOMS 43(4), 2017
// (Algorithm 973). Axom uses the explicit node/weight construction here,
// not the semi-automatic adaptive integrator from that paper.
inline void copy_matrix_to_array(const std::vector<std::vector<double>>& values,
                          axom::Array<double>& array,
                          int allocatorID)
{
  int rows = static_cast<int>(values.size());
  int cols = rows > 0 ? static_cast<int>(values[0].size()) : 0;
  array = axom::Array<double>(rows * cols, rows * cols, allocatorID);

  int index = 0;
  for(int row = 0; row < rows; ++row)
  {
    for(int col = 0; col < cols; ++col)
    {
      array[index++] = values[row][col];
    }
  }
}

inline void compute_rational_fejer_data_m11_impl(const std::vector<Pole>& poles,
                                          std::vector<double>& nodes,
                                          std::vector<double>& weights,
                                          RationalFejerDiagnostics* diagnostics = nullptr,
                                          int allocatorID = axom::getDefaultAllocatorID())
{
  assert("Rational Fejer quadrature requires at least one pole" && !poles.empty());

  const double infinite_pole_threshold = 2.0 / axom::numeric_limits<double>::epsilon();
  const double pole_tolerance = 2.0 * axom::numeric_limits<double>::epsilon();

  std::vector<Pole> canonical_poles = canonicalizePoleSequence(poles, pole_tolerance);
  for(auto& pole : canonical_poles)
  {
    pole = pole.normalizedInfinite(infinite_pole_threshold);
  }

  const int num_poles = static_cast<int>(canonical_poles.size());
  std::vector<Complex> cayley_poles(num_poles);
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

  if(diagnostics != nullptr)
  {
    diagnostics->canonical_poles_m11.clear();
    diagnostics->canonical_poles_m11.reserve(canonical_poles.size());
    for(const auto& pole : canonical_poles)
    {
      diagnostics->canonical_poles_m11.push_back(pole.value());
    }

    diagnostics->cayley_poles.clear();
    diagnostics->cayley_poles.reserve(cayley_poles.size());
    for(const auto& pole : cayley_poles)
    {
      diagnostics->cayley_poles.push_back(pole);
    }

    diagnostics->steps.clear();
  }

  std::vector<Pole> rcheb_poles = canonical_poles;
  rcheb_poles.push_back(Pole::infinity());

  std::vector<double> cheb_nodes;
  std::vector<double> cheb_weights;
  compute_rcheb_data(rcheb_poles, cheb_nodes, cheb_weights);

  if(diagnostics != nullptr)
  {
    copy_vector_to_array(cheb_nodes, diagnostics->rcheb_nodes_m11, allocatorID);
    copy_vector_to_array(cheb_weights, diagnostics->rcheb_weights_m11, allocatorID);
  }

  RationalFejerBasis basis(cayley_poles, std::move(cheb_nodes), std::move(cheb_weights));

  std::vector<double> basis_coefficients(num_poles + 1, 0.0);
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
      std::vector<int> repeated_pole_offsets;
      for(int idx = coeff_index - 1; idx < num_poles; ++idx)
      {
        if(pole.closeTo(canonical_poles[idx], pole_tolerance))
        {
          repeated_pole_offsets.push_back(idx - (coeff_index - 1) + 1);
        }
      }

      const int multiplicity = static_cast<int>(repeated_pole_offsets.size());
      if(component_count == 2)
      {
        std::vector<int> conjugate_positions;
        for(int idx = coeff_index - 1; idx < num_poles; ++idx)
        {
          if(pole.conjugate().closeTo(canonical_poles[idx], pole_tolerance))
          {
            conjugate_positions.push_back(idx - (coeff_index - 1) + 1);
          }
        }

        std::vector<int> interleaved_positions(2 * multiplicity);
        for(int i = 0; i < multiplicity; ++i)
        {
          interleaved_positions[2 * i] = repeated_pole_offsets[i];
          interleaved_positions[2 * i + 1] = conjugate_positions[i];
        }
        repeated_pole_offsets = std::move(interleaved_positions);
      }

      const auto pole_integrals =
        computePoleMomentIntegrals(pole, multiplicity, component_count, pole_tolerance);
      for(std::size_t i = 0; i < repeated_pole_offsets.size(); ++i)
      {
        basis_coefficients[coeff_index + repeated_pole_offsets[i] - 1] = pole_integrals[i];
      }
    }

    const int num_known_columns = coeff_index + component_count;
    std::vector<double> known_integrals(
      basis_coefficients.begin(), basis_coefficients.begin() + num_known_columns);

    RationalFejerDiagnostics::Step step_diagnostics;
    if(diagnostics != nullptr)
    {
      step_diagnostics.step_index = static_cast<int>(diagnostics->steps.size());
      step_diagnostics.coeff_index_before = coeff_index;
      step_diagnostics.component_count = component_count;
      step_diagnostics.pole_multiplicity_so_far = pole_multiplicity_so_far;
      step_diagnostics.pole_m11 = pole.value();
      copy_vector_to_array(known_integrals,
                           step_diagnostics.basis_coefficients_before,
                           allocatorID);
    }

    // Enforce orthogonality against the previously determined columns to recover
    // the next one or two rational basis coefficients.
    const auto orthogonal_integrals = basis.solveOrthogonalIntegrals(pole,
                                                                     pole_multiplicity_so_far,
                                                                     num_known_columns,
                                                                     known_integrals,
                                                                     component_count,
                                                                     diagnostics != nullptr ? &step_diagnostics : nullptr);
    for(int i = 0; i < component_count; ++i)
    {
      basis_coefficients[coeff_index + i] = orthogonal_integrals[i];
    }

    if(diagnostics != nullptr)
    {
      copy_vector_to_array(orthogonal_integrals,
                           step_diagnostics.orthogonal_integrals,
                           allocatorID);
      copy_vector_to_array(
        std::vector<double>(basis_coefficients.begin(), basis_coefficients.begin() + num_known_columns),
        step_diagnostics.basis_coefficients_after,
        allocatorID);
      diagnostics->steps.push_back(std::move(step_diagnostics));
    }

    coeff_index += component_count;
  }

  nodes = basis.nodes();
  const std::vector<double> basis_expansion = basis.evaluateExpansion(basis_coefficients);
  weights = basis_expansion;
  const auto& lambda = basis.lambda();
  for(std::size_t i = 0; i < weights.size(); ++i)
  {
    weights[i] *= lambda[i];
  }

  if(diagnostics != nullptr)
  {
    copy_vector_to_array(basis_coefficients, diagnostics->basis_coefficients, allocatorID);
    copy_vector_to_array(basis_expansion, diagnostics->basis_expansion_m11, allocatorID);
    copy_vector_to_array(weights, diagnostics->final_weights_m11, allocatorID);
    copy_matrix_to_array(basis.basisColumns(), diagnostics->basis_matrix_transpose_m11, allocatorID);
    diagnostics->basis_matrix_row_count = static_cast<int>(basis.basisColumns().size());
    diagnostics->basis_matrix_col_count =
      basis.basisColumns().empty() ? 0 : static_cast<int>(basis.basisColumns()[0].size());
  }
}

inline void compute_rational_fejer_data_m11(const std::vector<Pole>& poles,
                                     std::vector<double>& nodes,
                                     std::vector<double>& weights)
{
  compute_rational_fejer_data_m11_impl(poles, nodes, weights);
}
inline void compute_rational_fejer_diagnostics_m11(const std::vector<Pole>& poles,
                                            RationalFejerDiagnostics& diagnostics,
                                            int allocatorID)
{
  std::vector<double> nodes_m11;
  std::vector<double> weights_m11;
  compute_rational_fejer_data_m11_impl(
    poles, nodes_m11, weights_m11, &diagnostics, allocatorID);
}

void compute_rcheb_data_m11(const std::vector<Complex>& poles_m11,
                            std::vector<double>& nodes,
                            std::vector<double>& weights)
{
  std::vector<Pole> poles;
  poles.reserve(poles_m11.size());
  for(const auto& pole : poles_m11)
  {
    poles.emplace_back(pole);
  }
  compute_rcheb_data(poles, nodes, weights);
}

void compute_rational_fejer_data_m11(const std::vector<Complex>& poles_m11,
                                     std::vector<double>& nodes,
                                     std::vector<double>& weights)
{
  std::vector<Pole> poles;
  poles.reserve(poles_m11.size());
  for(const auto& pole : poles_m11)
  {
    poles.emplace_back(pole);
  }
  compute_rational_fejer_data_m11(poles, nodes, weights);
}

inline void compute_rational_fejer_data_m11_impl(const std::vector<Complex>& poles_m11,
                                          std::vector<double>& nodes,
                                          std::vector<double>& weights,
                                          RationalFejerDiagnostics* diagnostics,
                                          int allocatorID)
{
  std::vector<Pole> poles;
  poles.reserve(poles_m11.size());
  for(const auto& pole : poles_m11)
  {
    poles.emplace_back(pole);
  }
  compute_rational_fejer_data_m11_impl(poles, nodes, weights, diagnostics, allocatorID);
}

void compute_rational_fejer_diagnostics_m11(const std::vector<Complex>& poles_m11,
                                            RationalFejerDiagnostics& diagnostics,
                                            int allocatorID)
{
  std::vector<Pole> poles;
  poles.reserve(poles_m11.size());
  for(const auto& pole : poles_m11)
  {
    poles.emplace_back(pole);
  }
  compute_rational_fejer_diagnostics_m11(poles, diagnostics, allocatorID);
}
inline std::string make_rational_fejer_key(const std::vector<Pole>& poles01, int allocatorID)
{
  // Cache by the canonicalized pole sequence that the rule construction
  // actually uses. This preserves correctness while improving reuse when two
  // numerically equivalent pole lists differ only by tiny root-solver noise.
  const double infinite_pole_threshold = 2.0 / axom::numeric_limits<double>::epsilon();
  const double pole_tolerance = 2.0 * axom::numeric_limits<double>::epsilon();

  std::vector<Pole> canonical_poles = canonicalizePoleSequence(poles01, pole_tolerance);
  for(auto& pole : canonical_poles)
  {
    pole = pole.normalizedInfinite(infinite_pole_threshold);
  }

  axom::fmt::memory_buffer key;
  axom::fmt::format_to(std::back_inserter(key), "{}|", allocatorID);
  for(const auto& pole : canonical_poles)
  {
    axom::fmt::format_to(std::back_inserter(key),
                         "{:a},{:a};",
                         pole.real(),
                         pole.imag());
  }
  return axom::fmt::to_string(key);
}

std::string make_rational_fejer_key(const std::vector<Complex>& poles01, int allocatorID)
{
  std::vector<Pole> poles;
  poles.reserve(poles01.size());
  for(const auto& pole : poles01)
  {
    poles.emplace_back(pole);
  }
  return make_rational_fejer_key(poles, allocatorID);
}

}  // namespace internal

void compute_rational_fejer_data(const std::vector<std::complex<double>>& poles01,
                                 axom::Array<double>& nodes,
                                 axom::Array<double>& weights,
                                 int allocatorID)
{
  assert("Rational Fejer quadrature requires at least one pole" && !poles01.empty());

  std::vector<internal::Complex> poles_m11;
  poles_m11.reserve(poles01.size());
  for(const auto& pole : poles01)
  {
    const auto pole_magnitude = std::abs(pole);
    poles_m11.push_back(std::isinf(pole.real()) || std::isinf(pole.imag()) || std::isinf(pole_magnitude)
      ? internal::Complex {axom::numeric_limits<double>::infinity(), 0.0}
      : 2.0 * pole - internal::Complex {1.0, 0.0});
  }

  std::vector<double> nodes_m11;
  std::vector<double> weights_m11;

  // The implementation below follows the Deckers et al. / van Deun et al.
  // rational Fejer and rational Chebyshev constructions. QuaHOG's MATLAB
  // Spectral_PE implementation served as a secondary reference for how these
  // pieces are assembled in the planar integration workflow.
  internal::compute_rational_fejer_data_m11(poles_m11, nodes_m11, weights_m11);

  internal::map_rule_m11_to_unit_interval(nodes_m11, weights_m11, nodes, weights, allocatorID);
}

void internal::compute_rational_fejer_diagnostics(
  const std::vector<std::complex<double>>& poles01,
  internal::RationalFejerDiagnostics& diagnostics,
  int allocatorID)
{
  assert("Rational Fejer quadrature requires at least one pole" && !poles01.empty());

  std::vector<internal::Complex> poles_m11;
  poles_m11.reserve(poles01.size());
  for(const auto& pole : poles01)
  {
    const auto pole_magnitude = std::abs(pole);
    poles_m11.push_back(std::isinf(pole.real()) || std::isinf(pole.imag()) || std::isinf(pole_magnitude)
      ? internal::Complex {axom::numeric_limits<double>::infinity(), 0.0}
      : 2.0 * pole - internal::Complex {1.0, 0.0});
  }

  std::vector<double> nodes_m11;
  std::vector<double> weights_m11;
  internal::compute_rational_fejer_data_m11_impl(
    poles_m11, nodes_m11, weights_m11, &diagnostics, allocatorID);
  internal::map_rule_m11_to_unit_interval(
    nodes_m11, weights_m11, diagnostics.nodes_01, diagnostics.weights_01, allocatorID);
}

QuadratureRule get_rational_fejer(const std::vector<std::complex<double>>& poles01, int allocatorID)
{
  assert("Rational Fejer quadrature requires at least one pole" && !poles01.empty());

  static std::map<std::string, internal::RationalFejerRuleStorage> rule_library;
  static std::mutex rule_library_mutex;

  std::vector<internal::Complex> poles_key;
  poles_key.reserve(poles01.size());
  for(const auto& pole : poles01)
  {
    poles_key.emplace_back(pole);
  }

  const std::string key = internal::make_rational_fejer_key(poles_key, allocatorID);
  const std::lock_guard<std::mutex> lock(rule_library_mutex);
  auto [it, inserted] = rule_library.emplace(key, internal::RationalFejerRuleStorage {});
  if(inserted)
  {
    compute_rational_fejer_data(poles01, it->second.nodes, it->second.weights, allocatorID);
  }

  return QuadratureRule {it->second.nodes.view(), it->second.weights.view()};
}

}  // namespace numerics
}  // namespace axom
