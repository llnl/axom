// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/numerics/internal/rational_quadrature_common.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace axom
{
namespace numerics
{
namespace internal
{

/*!
 * \file rational_chebyshev.cpp
 *
 * \brief Implements the fixed-pole rational Chebyshev quadrature stage.
 *
 * This file implements the rational Chebyshev rule construction from
 * van Deun, Deckers, Bultheel, and Weideman, "Algorithm 882: Near-best
 * fixed pole rational interpolation with applications in spectral methods",
 * ACM Transactions on Mathematical Software 35(2), 2008.
 *
 * The algorithm maps poles outside `[-1,1]` into the Cayley domain, solves a
 * rational Chebyshev phase equation for the node angles, sorts the resulting
 * nodes on `[-1,1]`, and evaluates the corresponding rational Chebyshev
 * weights. The rational Fejer implementation uses this rule as its node set
 * and seed weighting stage.
 */

/// \brief Cayley-domain pole data used by rational Chebyshev formulas.
struct RationalChebyshevPoleData
{
  axom::Array<Complex> cayley_poles;
  axom::Array<double> multiplicities;
  axom::Array<double> radii;
  axom::Array<double> phases;
  axom::Array<Complex> branch_points;
};

/// \brief Computes the standard interior Chebyshev first-kind rule on [-1,1].
///
/// In the rational machinery below, this is the all-infinite-pole limit of the
/// rational Chebyshev stage.
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

/// \brief Constructs the rational Chebyshev rule on [-1,1].
///
/// This rule underpins the rational Fejer construction and corresponds to the
/// fixed-pole rational Gauss-Chebyshev / near-best interpolation construction
/// described in van Deun, Deckers, Bultheel, and Weideman, ACM TOMS 35(2),
/// 2008 (Algorithm 882).
class RationalChebyshevHelper
{
public:
  explicit RationalChebyshevHelper(PoleSequence poles)
    : m_poles(poles.normalizedForRationalChebyshev(infinitePoleThreshold()).toPoleArray())
  { }

  void compute(axom::Array<double>& nodes, axom::Array<double>& weights) const
  {
    axom::Array<Pole> poles = m_poles;
    const int n = static_cast<int>(poles.size());

    const double terminal_cayley_node = normalizeTerminalPole(poles);

    if(n == 1)
    {
      nodes = axom::Array<double> {terminal_cayley_node};
      weights = axom::Array<double> {M_PI};
      return;
    }

    const auto distinct_poles = collectDistinctPoles(poles, poleTolerance());
    if(distinct_poles.values.size() == 1 && terminal_cayley_node == 0.0)
    {
      // The all-infinite-pole case reduces to the standard Chebyshev first-kind rule.
      compute_chebyshev_first_kind_data_m11(n, nodes, weights);
      return;
    }

    const auto cayley_pole_data = buildCayleyPoleData(poles);
    const auto residual_grid = buildThetaResidualGrid(n, cayley_pole_data);

    axom::Array<double> node_angles = solveNodeAngles(n, cayley_pole_data, residual_grid);
    sortNodesAndAngles(node_angles, nodes);
    weights = computeWeights(node_angles, cayley_pole_data);
  }

private:
  struct ThetaResidualGrid
  {
    axom::Array<double> theta_samples;
    axom::Array<double> residual_samples;
  };

  double normalizeTerminalPole(axom::Array<Pole>& poles) const
  {
    const Pole last_pole = poles.back();
    double terminal_cayley_node = last_pole.cayleyTransform().real();
    if(last_pole.isInfinite())
    {
      return terminal_cayley_node;
    }

    if(std::abs(terminal_cayley_node) <= 1.0 / (2.0 * infinitePoleThreshold()))
    {
      poles.back() = Pole::infinity();
      return 0.0;
    }

    // Algorithm 882 treats the last pole specially: after the Cayley map
    // its real value fixes one endpoint of the phase equation, while the
    // pole list itself is reflected back to the original x-domain.
    poles.back() = Pole {0.5 * (terminal_cayley_node + 1.0 / terminal_cayley_node), 0.0};
    return terminal_cayley_node;
  }

  RationalChebyshevPoleData buildCayleyPoleData(axom::ArrayView<const Pole> poles) const
  {
    // Collapse the cyclic pole list into distinct Cayley-domain poles and
    // multiplicities before evaluating the angle contributions.
    const auto distinct_poles = PoleSequence {poles}.cyclic().distinct(poleTolerance());

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
      const Complex cayley_pole = distinct_poles.values[i].cayleyTransform();
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

  /// \brief Solves the rational Chebyshev phase condition for one node angle.
  ///
  /// The corresponding quadrature node on [-1,1] is then x = cos(theta).
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

  axom::Array<double> solveNodeAngles(int n,
                                      const RationalChebyshevPoleData& data,
                                      const ThetaResidualGrid& grid) const
  {
    axom::Array<double> node_angles(n, n);
    int bracket_start = 0;
    for(int i = 0; i < n; ++i)
    {
      const double target = M_PI * (n - i - 0.5);
      bool ok = solveTheta(n, data, grid, target, bracket_start, node_angles[i]);
      if(!ok)
      {
        failRationalFejerPrecondition(
          axom::fmt::format("Failed to construct rational Chebyshev node {} of {}.", i, n));
      }
    }
    return node_angles;
  }

  static void sortNodesAndAngles(axom::Array<double>& node_angles, axom::Array<double>& nodes)
  {
    const int n = static_cast<int>(node_angles.size());
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

  static double infinitePoleThreshold() { return 5.0 / axom::numeric_limits<double>::epsilon(); }

  static double poleTolerance() { return 5.0 * axom::numeric_limits<double>::epsilon(); }

  axom::Array<Pole> m_poles;
};

/// \brief Builds the rational Chebyshev rule used by the extended rational Fejer rule.
///
/// Adding a trailing pole at infinity produces the rational Chebyshev node set
/// that the rational Fejer weight solve is built on.
void compute_rational_chebyshev_data(axom::ArrayView<const Pole> poles_in,
                                     axom::Array<double>& nodes,
                                     axom::Array<double>& weights)
{
  RationalChebyshevHelper helper {PoleSequence {poles_in}};
  helper.compute(nodes, weights);
}

/// \brief Computes rational Chebyshev nodes and weights on [-1,1].
void compute_rational_chebyshev_data_m11(axom::ArrayView<const Complex> poles_m11,
                                         axom::Array<double>& nodes,
                                         axom::Array<double>& weights,
                                         int allocatorID)
{
  const PoleSequence poles = PoleSequence::fromComplex(poles_m11);
  poles.validate(-1.0, 1.0, "[-1,1]");

  axom::Array<double> nodes_tmp;
  axom::Array<double> weights_tmp;
  compute_rational_chebyshev_data(poles.view(), nodes_tmp, weights_tmp);
  copy_array_to_array(nodes_tmp, nodes, allocatorID);
  copy_array_to_array(weights_tmp, weights, allocatorID);
}

}  // namespace internal
}  // namespace numerics
}  // namespace axom
