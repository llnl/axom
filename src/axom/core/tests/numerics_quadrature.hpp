// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core/numerics/quadrature.hpp"
#include "axom/core/numerics/internal/rational_quadrature.hpp"
#include "axom/core/numerics/internal/rational_quadrature_common.hpp"
#include "axom/core/utilities/Utilities.hpp"

#include <algorithm>
#include <cmath>

#include <complex>
#include <functional>

namespace
{
namespace numerics_internal = axom::numerics::internal;
using Complex = std::complex<double>;
using Pole = numerics_internal::Pole;
using PoleSequence = numerics_internal::PoleSequence;

const double pi = std::acos(-1.0);
axom::Array<Complex> map_interval_poles_m11_to_01(axom::ArrayView<const Complex> poles_m11)
{
  axom::Array<Complex> poles01;
  poles01.reserve(poles_m11.size());
  for(const auto& pole_value : poles_m11)
  {
    const Pole pole {pole_value};
    if(pole.isInfinite())
    {
      poles01.push_back(Pole::infinity().value());
    }
    else
    {
      poles01.push_back(0.5 * (pole.value() + Complex {1.0, 0.0}));
    }
  }
  return poles01;
}

axom::Array<Complex> map_interval_poles_01_to_m11(axom::ArrayView<const Complex> poles01)
{
  return PoleSequence::from01ToM11(poles01).toComplexArray();
}

double integrate_rule(const axom::numerics::QuadratureRuleView& rule,
                      const std::function<double(double)>& integrand)
{
  double value = 0.0;
  for(int i = 0; i < rule.getNumPoints(); ++i)
  {
    value += rule.weight(i) * integrand(rule.node(i));
  }
  return value;
}

double integrate_rule_m11(axom::ArrayView<const double> nodes,
                          axom::ArrayView<const double> weights,
                          const std::function<double(double)>& integrand)
{
  double value = 0.0;
  for(int i = 0; i < nodes.size(); ++i)
  {
    value += weights[i] * integrand(nodes[i]);
  }
  return value;
}

void expect_rules_near(const axom::numerics::QuadratureRuleView& observed,
                       const axom::numerics::QuadratureRuleView& expected,
                       double tol)
{
  ASSERT_EQ(observed.getNumPoints(), expected.getNumPoints());
  for(int i = 0; i < observed.getNumPoints(); ++i)
  {
    EXPECT_NEAR(observed.node(i), expected.node(i), tol);
    EXPECT_NEAR(observed.weight(i), expected.weight(i), tol);
  }
}

double integrate_reference_m11(const std::function<double(double)>& integrand)
{
  constexpr int reference_points = 1024;
  const auto reference_rule = axom::numerics::get_gauss_legendre(reference_points);

  double value = 0.0;
  for(int i = 0; i < reference_rule.getNumPoints(); ++i)
  {
    const double x = 2.0 * reference_rule.node(i) - 1.0;
    value += 2.0 * reference_rule.weight(i) * integrand(x);
  }
  return value;
}

axom::Array<Complex> rotor_first_invalid_unique_poles_m11()
{
  // Regression poles from a two-rotors example dataset that previously exposed
  // unstable rational Fejer weights at moderate padded orders.
  return {Complex {-1.59075549109947922, 0.0},
          Complex {1.20225576665447242, -1.25283898687163120},
          Complex {1.20225576665447242, 1.25283898687163120}};
}

axom::Array<Complex> repeat_pole_sequence(axom::ArrayView<const Complex> poles, int repeat_count)
{
  axom::Array<Complex> repeated;
  repeated.reserve(poles.size() * repeat_count);
  for(int repeat = 0; repeat < repeat_count; ++repeat)
  {
    repeated.insert(repeated.size(), poles);
  }
  return repeated;
}

void pad_poles_with_infinities(axom::Array<Complex>& poles, int total_poles)
{
  const auto inf = axom::numeric_limits<double>::infinity();
  while(static_cast<int>(poles.size()) < total_poles)
  {
    poles.emplace_back(inf, 0.0);
  }
}

struct RationalFejerDiagnosticsSummary
{
  double min_rational_chebyshev_weight {0.0};
  double min_final_weight_m11 {0.0};
  double max_abs_final_weight_m11 {0.0};
  double sum_final_weights_m11 {0.0};
  double max_abs_basis_coefficient {0.0};
};

RationalFejerDiagnosticsSummary summarize_rule_diagnostics(
  const numerics_internal::RationalFejerDiagnostics& diagnostics)
{
  RationalFejerDiagnosticsSummary summary;
  summary.min_rational_chebyshev_weight = axom::numeric_limits<double>::infinity();
  summary.min_final_weight_m11 = axom::numeric_limits<double>::infinity();

  for(int i = 0; i < diagnostics.rational_chebyshev_weights_m11.size(); ++i)
  {
    summary.min_rational_chebyshev_weight =
      axom::utilities::min(summary.min_rational_chebyshev_weight,
                           diagnostics.rational_chebyshev_weights_m11[i]);
  }

  for(int i = 0; i < diagnostics.final_weights_m11.size(); ++i)
  {
    summary.min_final_weight_m11 =
      axom::utilities::min(summary.min_final_weight_m11, diagnostics.final_weights_m11[i]);
    summary.max_abs_final_weight_m11 =
      axom::utilities::max(summary.max_abs_final_weight_m11,
                           std::abs(diagnostics.final_weights_m11[i]));
    summary.sum_final_weights_m11 += diagnostics.final_weights_m11[i];
  }

  for(int i = 0; i < diagnostics.basis_coefficients.size(); ++i)
  {
    summary.max_abs_basis_coefficient =
      axom::utilities::max(summary.max_abs_basis_coefficient,
                           std::abs(diagnostics.basis_coefficients[i]));
  }

  return summary;
}

RationalFejerDiagnosticsSummary compute_rule_diagnostics_summary(axom::ArrayView<const Complex> poles01)
{
  numerics_internal::RationalFejerDiagnostics diagnostics;
  numerics_internal::compute_rational_fejer_diagnostics(poles01, diagnostics);
  return summarize_rule_diagnostics(diagnostics);
}

void expect_healthy_rational_fejer_summary(const RationalFejerDiagnosticsSummary& summary)
{
  EXPECT_GT(summary.min_rational_chebyshev_weight, 0.0);
  EXPECT_GT(summary.min_final_weight_m11, 0.0);
  EXPECT_NEAR(summary.sum_final_weights_m11, 2.0, 1e-10);
  EXPECT_LT(summary.max_abs_basis_coefficient, 1e8);
}

void compute_rational_chebyshev_rule_m11(axom::ArrayView<const Complex> poles_m11,
                                         axom::Array<double>& nodes,
                                         axom::Array<double>& weights)
{
  numerics_internal::compute_rational_chebyshev_data_m11(poles_m11, nodes, weights);
}

void compute_rational_fejer_rule_m11(axom::ArrayView<const Complex> poles_m11,
                                     axom::Array<double>& nodes,
                                     axom::Array<double>& weights)
{
  numerics_internal::compute_rational_fejer_data_m11(poles_m11, nodes, weights);
}
}  // namespace

//------------------------------------------------------------------------------
// Gauss-Legendre tests
//------------------------------------------------------------------------------

TEST(numerics_quadrature, gauss_legendre_math_check)
{
  const int N = 200;

  double coeffs[2 * N - 1];

  // Test that the rules provide exact integration for polynomials of degree 2n - 1
  for(int npts = 1; npts <= N; ++npts)
  {
    // Evaluate using the quadrature rule
    auto rule = axom::numerics::get_gauss_legendre(npts);
    int degree = 2 * npts - 1;

    // Define a collection of random coefficients for a polynomial
    for(int j = 0; j < degree; ++j)
    {
      // Seed the random coefficients for consistency in the test
      coeffs[j] = axom::utilities::random_real(-1.0, 1.0, npts + j);
    }

    // Evaluate the area under the curve from 0 to 1 analytically
    double analytic_result = 0.0;
    for(int j = 0; j < degree; ++j)
    {
      analytic_result += coeffs[j] / (j + 1);
    }

    // Evaluate the polynomial using Horner's rule
    auto eval_polynomial = [&degree, &coeffs](double x) {
      double result = coeffs[degree - 1];
      for(int i = degree - 2; i >= 0; --i)
      {
        result = result * x + coeffs[i];
      }
      return result;
    };

    double quadrature_result = 0.0;
    for(int j = 0; j < npts; ++j)
    {
      quadrature_result += rule.weight(j) * eval_polynomial(rule.node(j));
    }

    EXPECT_NEAR(quadrature_result, analytic_result, 1e-10);

    // Check that nodes aren't duplicated
    for(int j = 1; j < npts; ++j)
    {
      EXPECT_GT(rule.node(j), rule.node(j - 1));
    }

    // Check that the sum of the weights is 1, and that all are positive
    double weight_sum = 0.0;
    for(int j = 0; j < npts; ++j)
    {
      EXPECT_GT(rule.weight(j), 0.0);
      weight_sum += rule.weight(j);
    }

    EXPECT_NEAR(weight_sum, 1.0, 1e-10);

    // Check that each node is the root of the next Legendre polynomial
    for(int j = 0; j < npts; ++j)
    {
      // Rescale the node to [-1, 1]
      double z = 2 * rule.node(j) - 1;

      double P_n = z, P_nm1 = 1.0;
      for(int k = 2; k <= npts; ++k)
      {
        double P_nm2 = P_nm1;
        P_nm1 = P_n;
        P_n = ((2 * k - 1) * z * P_nm1 - (k - 1) * P_nm2) / k;
      }

      // Note that this metric does not directly imply that |z - z*| < tol
      EXPECT_NEAR(P_n, 0.0, 1e-10);
    }
  }
}

template <typename ExecSpace>
struct test_gauss_legendre_device_quadrature
{
  static void test()
  {
    const int npts = 15;

    // Create the rule with the proper allocator
    int allocID;
#if defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_GPU)

    // Use unified memory so the cached rule storage can be accessed from device kernels.
    constexpr bool on_device = axom::execution_space<ExecSpace>::onDevice();

    allocID = on_device ? axom::getUmpireResourceAllocatorID(umpire::resource::Unified)
                        : axom::execution_space<ExecSpace>::allocatorID();
#else
    allocID = axom::execution_space<ExecSpace>::allocatorID();
#endif

    const auto rule = axom::numerics::get_gauss_legendre(npts, allocID);

    // Capture the non-owning rule view in a device lambda.
    axom::ReduceSum<ExecSpace, double> quadrature_sum(0.0);
    axom::for_all<ExecSpace>(
      npts,
      AXOM_LAMBDA(axom::IndexType i) { quadrature_sum += rule.weight(i) * sin(rule.node(i)); });

    EXPECT_NEAR(quadrature_sum.get(), 0.459697694132, 1e-6);
  }
};

TEST(numerics_quadrature, gauss_legendre_device_seq)
{
  test_gauss_legendre_device_quadrature<axom::SEQ_EXEC>::test();
}

#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(numerics_quadrature, gauss_legendre_device_omp)
{
  test_gauss_legendre_device_quadrature<axom::OMP_EXEC>::test();
}
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(numerics_quadrature, gauss_legendre_device_cuda)
{
  test_gauss_legendre_device_quadrature<axom::CUDA_EXEC<256>>::test();
}
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(numerics_quadrature, gauss_legendre_device_hip)
{
  test_gauss_legendre_device_quadrature<axom::HIP_EXEC<256>>::test();
}
#endif

//------------------------------------------------------------------------------
// Algorithm 882 rational Chebyshev tests
//------------------------------------------------------------------------------

TEST(numerics_quadrature,
     rational_chebyshev_internal_matches_algorithm882_pure_imaginary_poles_example)
{
  axom::Array<Complex> poles_m11;
  for(int j = 1; j <= 10; ++j)
  {
    poles_m11.emplace_back(0.0, 0.001 * j);
    poles_m11.emplace_back(0.0, -0.001 * j);
  }

  axom::Array<double> nodes;
  axom::Array<double> weights;
  compute_rational_chebyshev_rule_m11(poles_m11, nodes, weights);

  ASSERT_EQ(nodes.size(), poles_m11.size());
  ASSERT_EQ(weights.size(), poles_m11.size());

  double weight_sum = 0.0;
  for(int i = 0; i < static_cast<int>(nodes.size()); ++i)
  {
    EXPECT_GT(weights[i], 0.0);
    weight_sum += weights[i];

    const int mirrored = static_cast<int>(nodes.size()) - 1 - i;
    EXPECT_NEAR(nodes[i], -nodes[mirrored], 1e-12);
    EXPECT_NEAR(weights[i], weights[mirrored], 1e-12);
  }

  EXPECT_NEAR(weight_sum, pi, 1e-12);
}

TEST(numerics_quadrature,
     rational_chebyshev_internal_matches_algorithm882_all_infinite_poles_chebyshev_limit)
{
  constexpr int pole_count = 6;
  axom::Array<Complex> poles_m11;
  poles_m11.reserve(pole_count);
  for(int i = 0; i < pole_count; ++i)
  {
    poles_m11.emplace_back(axom::numeric_limits<double>::infinity(), 0.0);
  }

  axom::Array<double> nodes;
  axom::Array<double> weights;
  compute_rational_chebyshev_rule_m11(poles_m11, nodes, weights);

  ASSERT_EQ(static_cast<int>(nodes.size()), pole_count);
  ASSERT_EQ(static_cast<int>(weights.size()), pole_count);

  for(int i = 0; i < pole_count; ++i)
  {
    const double theta = pi * (pole_count - i - 0.5) / pole_count;
    EXPECT_NEAR(nodes[i], std::cos(theta), 1e-14);
    EXPECT_NEAR(weights[i], pi / pole_count, 1e-14);
  }
}

TEST(numerics_quadrature,
     rational_chebyshev_internal_matches_algorithm882_boundary_layer_application_poles)
{
  axom::Array<Complex> poles_m11 = {{0.0, 0.0403},
                                    {0.0, -0.0403},
                                    {0.0094, 0.0398},
                                    {0.0094, -0.0398},
                                    {-0.0094, 0.0398},
                                    {-0.0094, -0.0398},
                                    {0.0200, 0.0384},
                                    {0.0200, -0.0384},
                                    {-0.0200, 0.0384},
                                    {-0.0200, -0.0384}};
  for(int i = 0; i < 10; ++i)
  {
    poles_m11.emplace_back(axom::numeric_limits<double>::infinity(), 0.0);
  }

  axom::Array<double> nodes;
  axom::Array<double> weights;
  compute_rational_chebyshev_rule_m11(poles_m11, nodes, weights);

  ASSERT_EQ(nodes.size(), poles_m11.size());
  ASSERT_EQ(weights.size(), poles_m11.size());

  double weight_sum = 0.0;
  double min_abs_node = axom::numeric_limits<double>::infinity();
  int nodes_near_boundary_layer = 0;
  for(int i = 0; i < static_cast<int>(nodes.size()); ++i)
  {
    EXPECT_GT(weights[i], 0.0);
    weight_sum += weights[i];
    min_abs_node = axom::utilities::min(min_abs_node, std::abs(nodes[i]));
    if(std::abs(nodes[i]) < 0.05)
    {
      ++nodes_near_boundary_layer;
    }

    const int mirrored = static_cast<int>(nodes.size()) - 1 - i;
    EXPECT_NEAR(nodes[i], -nodes[mirrored], 1e-10);
    EXPECT_NEAR(weights[i], weights[mirrored], 1e-10);
  }

  EXPECT_NEAR(weight_sum, pi, 1e-10);
  EXPECT_LT(min_abs_node, 0.02);
  EXPECT_GE(nodes_near_boundary_layer, 4);
}

TEST(numerics_quadrature, rational_chebyshev_internal_matches_algorithm882_near_interval_stress)
{
  constexpr double eps = axom::numeric_limits<double>::epsilon();
  axom::Array<Complex> poles_m11;
  poles_m11.reserve(70);
  for(int repeat = 0; repeat < 10; ++repeat)
  {
    for(int j = 0; j < 7; ++j)
    {
      poles_m11.emplace_back(-0.6 + 0.2 * j, 100.0 * eps);
    }
  }

  axom::Array<double> nodes;
  axom::Array<double> weights;
  compute_rational_chebyshev_rule_m11(poles_m11, nodes, weights);

  ASSERT_EQ(nodes.size(), poles_m11.size());
  ASSERT_EQ(weights.size(), poles_m11.size());

  double weight_sum = 0.0;
  for(int i = 0; i < static_cast<int>(nodes.size()); ++i)
  {
    EXPECT_GE(nodes[i], -1.0);
    EXPECT_LE(nodes[i], 1.0);
    EXPECT_TRUE(std::isfinite(weights[i]));
    weight_sum += weights[i];
  }

  EXPECT_NEAR(weight_sum / pi, 1.0, 1e-7);
}

TEST(numerics_quadrature, rational_chebyshev_internal_matches_algorithm882_repeated_pole_scaling_guard)
{
  axom::Array<Complex> poles_m11;
  constexpr int repeat_count = 100;
  poles_m11.reserve(3 * repeat_count);
  for(int repeat = 0; repeat < repeat_count; ++repeat)
  {
    poles_m11.emplace_back(-1.1, 0.0);
    poles_m11.emplace_back(0.0, 0.1);
    poles_m11.emplace_back(1.1, 0.0);
  }

  axom::Array<double> nodes;
  axom::Array<double> weights;
  compute_rational_chebyshev_rule_m11(poles_m11, nodes, weights);

  ASSERT_EQ(nodes.size(), poles_m11.size());
  ASSERT_EQ(weights.size(), poles_m11.size());

  double weight_sum = 0.0;
  for(int i = 0; i < static_cast<int>(nodes.size()); ++i)
  {
    EXPECT_GE(nodes[i], -1.0);
    EXPECT_LE(nodes[i], 1.0);
    EXPECT_GT(weights[i], 0.0);
    weight_sum += weights[i];
  }

  EXPECT_NEAR(weight_sum, pi, 1e-11);
}

//------------------------------------------------------------------------------
// Rational Fejer API and cache tests
//------------------------------------------------------------------------------

TEST(numerics_quadrature_DeathTest, rational_fejer_rejects_invalid_public_poles)
{
  axom::Array<double> nodes;
  axom::Array<double> weights;

  axom::Array<Complex> empty_poles;
  EXPECT_DEATH_IF_SUPPORTED(axom::numerics::compute_rational_fejer_data(empty_poles, nodes, weights),
                            "");

  const axom::Array<Complex> pole_on_unit_interval {Complex {0.5, 0.0}};
  EXPECT_DEATH_IF_SUPPORTED(
    axom::numerics::compute_rational_fejer_data(pole_on_unit_interval, nodes, weights),
    "");

  const axom::Array<Complex> nan_pole {Complex {axom::numeric_limits<double>::quiet_NaN(), 0.0}};
  EXPECT_DEATH_IF_SUPPORTED(axom::numerics::get_rational_fejer(nan_pole), "");
}

TEST(numerics_quadrature_DeathTest, rational_fejer_rejects_invalid_internal_poles)
{
  axom::Array<double> nodes;
  axom::Array<double> weights;

  axom::Array<Complex> empty_poles;
  EXPECT_DEATH_IF_SUPPORTED(compute_rational_fejer_rule_m11(empty_poles, nodes, weights), "");

  const axom::Array<Complex> pole_on_reference_interval {Complex {0.0, 0.0}};
  EXPECT_DEATH_IF_SUPPORTED(
    compute_rational_chebyshev_rule_m11(pole_on_reference_interval, nodes, weights),
    "");
}

TEST(numerics_quadrature, rational_fejer_infinite_poles_are_polynomial_exact)
{
  constexpr double inf = axom::numeric_limits<double>::infinity();

  for(int npoles = 1; npoles <= 8; ++npoles)
  {
    axom::Array<Complex> poles;
    poles.reserve(npoles);
    for(int i = 0; i < npoles; ++i)
    {
      poles.emplace_back(inf, 0.0);
    }
    auto rule = axom::numerics::get_rational_fejer(poles);

    const int degree = npoles;
    axom::Array<double> coeffs(degree + 1, degree + 1);
    for(int i = 0; i <= degree; ++i)
    {
      coeffs[i] = axom::utilities::random_real(-1.0, 1.0, 100 + 13 * npoles + i);
    }

    auto polynomial = [&coeffs](double x) {
      double value = coeffs.back();
      for(int i = static_cast<int>(coeffs.size()) - 2; i >= 0; --i)
      {
        value = value * x + coeffs[i];
      }
      return value;
    };

    double exact = 0.0;
    for(int i = 0; i <= degree; ++i)
    {
      exact += coeffs[i] / (i + 1);
    }

    double observed = 0.0;
    double weight_sum = 0.0;
    for(int i = 0; i < rule.getNumPoints(); ++i)
    {
      EXPECT_GE(rule.node(i), 0.0);
      EXPECT_LE(rule.node(i), 1.0);
      EXPECT_GT(rule.weight(i), 0.0);
      observed += rule.weight(i) * polynomial(rule.node(i));
      weight_sum += rule.weight(i);
    }

    EXPECT_NEAR(weight_sum, 1.0, 1e-10);
    EXPECT_NEAR(observed, exact, 1e-10);
  }
}

TEST(numerics_quadrature, rational_fejer_matches_simple_pole_moments)
{
  axom::Array<Complex> poles;
  for(int a = 2; a <= 10; ++a)
  {
    poles.emplace_back(0.5 * (a + 1.0), 0.0);
  }

  auto rule = axom::numerics::get_rational_fejer(poles);

  for(double pole_value = 1.5; pole_value <= 5.5; pole_value += 0.5)
  {
    double exact = std::log(std::abs(1.0 - pole_value)) - std::log(std::abs(pole_value));
    double observed = 0.0;
    for(int i = 0; i < rule.getNumPoints(); ++i)
    {
      observed += rule.weight(i) / (rule.node(i) - pole_value);
    }

    EXPECT_NEAR(observed, exact, 1e-9);
  }
}

TEST(numerics_quadrature, rational_fejer_auto_adds_missing_complex_conjugate_pole)
{
  const double a = 1.4;
  const double b = 0.8;
  const Complex pole {a, b};

  const axom::Array<Complex> implicit_poles {pole};
  const axom::Array<Complex> explicit_poles {pole, std::conj(pole)};
  const auto implicit_rule = axom::numerics::get_rational_fejer(implicit_poles);
  const auto explicit_rule = axom::numerics::get_rational_fejer(explicit_poles);

  expect_rules_near(implicit_rule, explicit_rule, 1e-14);

  const auto integrand = [a, b](double x) {
    const double dx = x - a;
    return 1.0 / (dx * dx + b * b);
  };
  const double exact = (std::atan((1.0 - a) / b) - std::atan(-a / b)) / b;

  EXPECT_NEAR(integrate_rule(implicit_rule, integrand), exact, 5e-13);
  EXPECT_NEAR(integrate_rule(explicit_rule, integrand), exact, 5e-13);
}

TEST(numerics_quadrature, rational_fejer_matches_repeated_complex_pair_exactness)
{
  const double a = 1.25;
  const double b = 0.6;
  const Complex pole {a, b};
  const axom::Array<Complex> poles {pole, pole};
  const auto rule = axom::numerics::get_rational_fejer(poles);

  const auto integrand = [a, b](double x) {
    const double dx = x - a;
    const double denom = dx * dx + b * b;
    return 1.0 / (denom * denom);
  };

  const auto antiderivative = [a, b](double x) {
    const double dx = x - a;
    const double denom = dx * dx + b * b;
    return dx / (2.0 * b * b * denom) + std::atan(dx / b) / (2.0 * b * b * b);
  };

  const double exact = antiderivative(1.0) - antiderivative(0.0);
  EXPECT_NEAR(integrate_rule(rule, integrand), exact, 5e-12);
}

TEST(numerics_quadrature, rational_fejer_point_count_follows_canonical_poles)
{
  const axom::Array<Complex> real_poles {Complex {1.25, 0.0}, Complex {1.5, 0.0}, Complex {2.0, 0.0}};
  axom::Array<double> nodes;
  axom::Array<double> weights;
  axom::numerics::compute_rational_fejer_data(real_poles, nodes, weights);
  EXPECT_EQ(nodes.size(), real_poles.size() + 1);
  EXPECT_EQ(weights.size(), real_poles.size() + 1);
  EXPECT_EQ(axom::numerics::get_rational_fejer(real_poles).getNumPoints(), real_poles.size() + 1);

  const axom::Array<Complex> explicit_pair {Complex {1.4, 0.8}, Complex {1.4, -0.8}};
  axom::numerics::compute_rational_fejer_data(explicit_pair, nodes, weights);
  EXPECT_EQ(nodes.size(), explicit_pair.size() + 1);
  EXPECT_EQ(weights.size(), explicit_pair.size() + 1);
  EXPECT_EQ(axom::numerics::get_rational_fejer(explicit_pair).getNumPoints(),
            explicit_pair.size() + 1);

  const axom::Array<Complex> implicit_pair {Complex {1.4, 0.8}};
  axom::numerics::compute_rational_fejer_data(implicit_pair, nodes, weights);
  EXPECT_EQ(nodes.size(), 3);
  EXPECT_EQ(weights.size(), 3);
  EXPECT_EQ(axom::numerics::get_rational_fejer(implicit_pair).getNumPoints(), 3);
}

TEST(numerics_quadrature, rational_fejer_cached_rule_matches_uncached_compute)
{
  const axom::Array<Complex> poles {Complex {1.4, 0.8}, Complex {1.4, -0.8}, Complex {2.5, 0.0}};

  axom::Array<double> computed_nodes;
  axom::Array<double> computed_weights;
  axom::numerics::compute_rational_fejer_data(poles, computed_nodes, computed_weights);

  const auto cached_rule = axom::numerics::get_rational_fejer(poles);
  const axom::numerics::QuadratureRuleView computed_rule {computed_nodes.view(),
                                                          computed_weights.view()};

  expect_rules_near(cached_rule, computed_rule, 1e-14);
}

TEST(numerics_quadrature, rational_fejer_cache_key_uses_canonical_m11_poles)
{
  const axom::Array<Complex> explicit_poles_m11 {Complex {1.8, 0.6},
                                                 Complex {1.8, -0.6},
                                                 Complex {2.5, 0.0}};
  const axom::Array<Complex> reordered_poles_m11 {Complex {1.8, -0.6},
                                                  Complex {1.8, 0.6},
                                                  Complex {2.5, 0.0}};
  const axom::Array<Complex> implicit_conjugate_m11 {Complex {1.8, 0.6}, Complex {2.5, 0.0}};

  // Conjugates are canonicalized, but pole order is part of the rule definition.
  const auto explicit_key = numerics_internal::make_rational_fejer_cache_key_m11(explicit_poles_m11);
  EXPECT_EQ(explicit_key,
            numerics_internal::make_rational_fejer_cache_key_m11(implicit_conjugate_m11));
  EXPECT_NE(explicit_key, numerics_internal::make_rational_fejer_cache_key_m11(reordered_poles_m11));
}

TEST(numerics_quadrature, rational_fejer_cache_key_canonicalizes_m11_near_duplicates)
{
  constexpr double eps = axom::numeric_limits<double>::epsilon();
  const Complex base_pole {0.0, 0.2};
  const axom::Array<Complex> base_poles_m11 {base_pole, base_pole};
  const axom::Array<Complex> close_poles_m11 {base_pole, base_pole * (1.0 + eps)};
  const axom::Array<Complex> distinct_poles_m11 {base_pole, base_pole * (1.0 + 64.0 * eps)};

  // The cache key follows PoleSequence's tolerance for coalescing repeated poles.
  const auto base_key = numerics_internal::make_rational_fejer_cache_key_m11(base_poles_m11);
  EXPECT_EQ(base_key, numerics_internal::make_rational_fejer_cache_key_m11(close_poles_m11));
  EXPECT_NE(base_key, numerics_internal::make_rational_fejer_cache_key_m11(distinct_poles_m11));
}

TEST(numerics_quadrature, rational_fejer_cached_rule_uses_m11_canonical_key)
{
  const axom::Array<Complex> explicit_poles01 {Complex {1.4, 0.3},
                                               Complex {1.4, -0.3},
                                               Complex {1.75, 0.0}};
  const axom::Array<Complex> implicit_conjugate_poles01 {Complex {1.4, 0.3}, Complex {1.75, 0.0}};

  const auto explicit_key = numerics_internal::make_rational_fejer_cache_key_m11(
    map_interval_poles_01_to_m11(explicit_poles01));
  const auto implicit_key = numerics_internal::make_rational_fejer_cache_key_m11(
    map_interval_poles_01_to_m11(implicit_conjugate_poles01));
  ASSERT_EQ(explicit_key, implicit_key);

  axom::Array<double> computed_nodes;
  axom::Array<double> computed_weights;
  axom::numerics::compute_rational_fejer_data(explicit_poles01, computed_nodes, computed_weights);

  const auto explicit_cached_rule = axom::numerics::get_rational_fejer(explicit_poles01);
  const auto cached_rule = axom::numerics::get_rational_fejer(implicit_conjugate_poles01);
  const axom::numerics::QuadratureRuleView computed_rule {computed_nodes.view(),
                                                          computed_weights.view()};

  expect_rules_near(explicit_cached_rule, cached_rule, 1e-14);
  expect_rules_near(cached_rule, computed_rule, 1e-14);
}

TEST(numerics_quadrature, rational_fejer_cached_rule_views_can_be_copied_to_owned_rule)
{
  const axom::Array<Complex> poles {Complex {1.4, 0.8}, Complex {1.4, -0.8}, Complex {2.5, 0.0}};
  const auto cached_rule = axom::numerics::get_rational_fejer(poles);

  const auto owned_rule = cached_rule.copy();
  const auto copied_rule = owned_rule.view();

  expect_rules_near(owned_rule.view(), cached_rule, 1e-14);
  expect_rules_near(copied_rule, cached_rule, 1e-14);
}

#if defined(AXOM_USE_UMPIRE)
TEST(numerics_quadrature, rational_fejer_diagnostics_respect_allocator)
{
  const int allocatorID = axom::getUmpireResourceAllocatorID(umpire::resource::Pinned);
  const axom::Array<Complex> poles {Complex {1.4, 0.8}, Complex {1.4, -0.8}, Complex {2.5, 0.0}};

  numerics_internal::RationalFejerDiagnostics diagnostics;
  numerics_internal::compute_rational_fejer_diagnostics(poles, diagnostics, allocatorID);

  EXPECT_EQ(diagnostics.canonical_poles_m11.getAllocatorID(), allocatorID);
  EXPECT_EQ(diagnostics.cayley_poles.getAllocatorID(), allocatorID);
  EXPECT_EQ(diagnostics.rational_chebyshev_nodes_m11.getAllocatorID(), allocatorID);
  EXPECT_EQ(diagnostics.rational_chebyshev_weights_m11.getAllocatorID(), allocatorID);
  EXPECT_EQ(diagnostics.basis_coefficients.getAllocatorID(), allocatorID);
  EXPECT_EQ(diagnostics.final_weights_m11.getAllocatorID(), allocatorID);
  EXPECT_EQ(diagnostics.nodes_01.getAllocatorID(), allocatorID);
  EXPECT_EQ(diagnostics.weights_01.getAllocatorID(), allocatorID);
  EXPECT_EQ(diagnostics.steps.getAllocatorID(), allocatorID);

  ASSERT_GT(diagnostics.steps.size(), 0);
  const auto& step = diagnostics.steps[0];
  EXPECT_EQ(step.basis_coefficients_before.getAllocatorID(), allocatorID);
  EXPECT_EQ(step.weighted_row0.getAllocatorID(), allocatorID);
  EXPECT_EQ(step.projected_row0.getAllocatorID(), allocatorID);
  EXPECT_EQ(step.projected_row0_terms.getAllocatorID(), allocatorID);
  EXPECT_EQ(step.orthogonal_integrals.getAllocatorID(), allocatorID);
  EXPECT_EQ(step.basis_coefficients_after.getAllocatorID(), allocatorID);
}
#endif

//------------------------------------------------------------------------------
// Algorithm 973 examples and regression tests
//------------------------------------------------------------------------------

TEST(numerics_quadrature, rotor_pole_family_stays_healthy_for_moderate_orders)
{
  const auto base_poles01 = map_interval_poles_m11_to_01(rotor_first_invalid_unique_poles_m11());

  // Exercise the original poles, repeated poles, and infinity-padded variants from the regression.
  const auto unique_summary = compute_rule_diagnostics_summary(base_poles01);
  expect_healthy_rational_fejer_summary(unique_summary);

  auto repeated_poles01 = repeat_pole_sequence(base_poles01, 5);
  const auto repeated_summary = compute_rule_diagnostics_summary(repeated_poles01);
  expect_healthy_rational_fejer_summary(repeated_summary);

  auto summary_for_total_poles = [&base_poles01](int total_poles) {
    auto poles01 = base_poles01;
    pad_poles_with_infinities(poles01, total_poles);
    return compute_rule_diagnostics_summary(poles01);
  };

  const auto summary24 = summary_for_total_poles(24);
  expect_healthy_rational_fejer_summary(summary24);
  EXPECT_LT(summary24.max_abs_basis_coefficient, 10.0);
}

TEST(numerics_quadrature, rational_fejer_internal_matches_algorithm973_example1_exactness)
{
  axom::Array<Complex> poles_m11;
  for(int a = 2; a <= 10; ++a)
  {
    poles_m11.emplace_back(static_cast<double>(a), 0.0);
  }

  axom::Array<double> forward_nodes;
  axom::Array<double> forward_weights;
  compute_rational_fejer_rule_m11(poles_m11, forward_nodes, forward_weights);

  auto reversed_poles_m11 = poles_m11;
  std::reverse(reversed_poles_m11.begin(), reversed_poles_m11.end());
  axom::Array<double> reversed_nodes;
  axom::Array<double> reversed_weights;
  compute_rational_fejer_rule_m11(reversed_poles_m11, reversed_nodes, reversed_weights);

  ASSERT_EQ(forward_nodes.size(), reversed_nodes.size());
  ASSERT_EQ(forward_weights.size(), reversed_weights.size());

  for(int a = 2; a <= 10; ++a)
  {
    const double exact = std::log(std::abs(1.0 - a)) - std::log(std::abs(-1.0 - a));
    const auto integrand = [a](double x) { return 1.0 / (x - a); };
    const double forward_value = integrate_rule_m11(forward_nodes, forward_weights, integrand);
    const double reversed_value = integrate_rule_m11(reversed_nodes, reversed_weights, integrand);

    EXPECT_NEAR(forward_value, exact, 5e-13);
    EXPECT_NEAR(reversed_value, exact, 5e-13);
  }
}

TEST(numerics_quadrature, rational_fejer_internal_matches_algorithm973_example2_repeated_real_pole)
{
  constexpr double A = 2.0;
  axom::Array<Complex> poles_m11;
  poles_m11.reserve(32);
  for(int i = 0; i < 32; ++i)
  {
    poles_m11.emplace_back(A, 0.0);
  }
  axom::Array<double> nodes;
  axom::Array<double> weights;
  compute_rational_fejer_rule_m11(poles_m11, nodes, weights);

  const auto integrand = [A](double x) { return 1.0 + 2.0 / std::pow(x - A, 100); };

  const double exact = 2.0 - (2.0 / 99.0) * (std::pow(1.0 - A, -99) + std::pow(1.0 + A, -99));
  const double value = integrate_rule_m11(nodes, weights, integrand);

  EXPECT_NEAR(value, exact, std::abs(exact) * 1e-11);
}

TEST(numerics_quadrature, rational_fejer_internal_matches_algorithm973_example7_1_table_scale)
{
  struct Case
  {
    double omega;
    int node_count;
    double paper_relative_error;
  };

  const Case cases[] = {{1.1, 2, 1.59e-01},
                        {1.1, 4, 3.77e-04},
                        {1.1, 8, 1.39e-08},
                        {1.1, 12, 9.42e-14},
                        {1.1, 16, 1.99e-16},
                        {1.001, 2, 4.06e-01},
                        {1.001, 4, 4.20e-03},
                        {1.001, 8, 4.98e-08},
                        {1.001, 12, 3.21e-13},
                        {1.001, 16, 1.29e-13}};

  for(const auto& test_case : cases)
  {
    axom::Array<Complex> poles_m11;
    poles_m11.reserve(test_case.node_count - 1);
    for(int k = 1; k < test_case.node_count; ++k)
    {
      const double sign = (k % 2 == 1) ? 1.0 : -1.0;
      const int multiple = (k + 1) / 2;
      poles_m11.emplace_back(sign * multiple * test_case.omega, 0.0);
    }

    axom::Array<double> nodes;
    axom::Array<double> weights;
    compute_rational_fejer_rule_m11(poles_m11, nodes, weights);
    ASSERT_EQ(nodes.size(), test_case.node_count);
    ASSERT_EQ(weights.size(), test_case.node_count);

    const double omega = test_case.omega;
    const auto integrand = [omega](double x) {
      const double scaled_x = pi * x / omega;
      return scaled_x / std::sin(scaled_x);
    };

    const double reference = integrate_reference_m11(integrand);
    const double observed = integrate_rule_m11(nodes, weights, integrand);
    const double relative_error = std::abs(observed - reference) / std::abs(reference);
    const double tolerance = axom::utilities::max(1e-12, 100.0 * test_case.paper_relative_error);
    EXPECT_LE(relative_error, tolerance)
      << "omega=" << test_case.omega << ", node_count=" << test_case.node_count;
  }
}

TEST(numerics_quadrature, rational_fejer_internal_matches_algorithm973_example7_4_table_scale)
{
  struct Case
  {
    int node_count;
    double paper_relative_error;
  };

  constexpr double omega = 1.1;
  const Case cases[] = {{10, 2.79e-05}, {20, 9.07e-14}, {30, 9.50e-15}};

  const auto integrand = [omega](double x) { return std::sin(1.0 / (x - omega)); };
  const double reference = integrate_reference_m11(integrand);

  for(const auto& test_case : cases)
  {
    axom::Array<Complex> poles_m11;
    poles_m11.reserve(test_case.node_count - 1);
    for(int i = 0; i < test_case.node_count - 1; ++i)
    {
      poles_m11.emplace_back(omega, 0.0);
    }

    axom::Array<double> nodes;
    axom::Array<double> weights;
    compute_rational_fejer_rule_m11(poles_m11, nodes, weights);
    ASSERT_EQ(nodes.size(), test_case.node_count);
    ASSERT_EQ(weights.size(), test_case.node_count);

    const double observed = integrate_rule_m11(nodes, weights, integrand);
    const double relative_error = std::abs(observed - reference) / std::abs(reference);
    const double tolerance = axom::utilities::max(1e-12, 100.0 * test_case.paper_relative_error);
    EXPECT_LE(relative_error, tolerance) << "node_count=" << test_case.node_count;
  }
}

TEST(numerics_quadrature, rational_fejer_internal_matches_algorithm973_example7_5_symmetry)
{
  constexpr double omega = 0.1;
  constexpr int half_pole_count = 8;
  constexpr int node_count = 2 * half_pole_count + 1;

  axom::Array<Complex> poles_m11;
  poles_m11.reserve(2 * half_pole_count);
  for(int k = -(half_pole_count - 1); k <= half_pole_count; ++k)
  {
    poles_m11.emplace_back(0.0, (2 * k - 1) * omega);
  }

  axom::Array<double> nodes;
  axom::Array<double> weights;
  compute_rational_fejer_rule_m11(poles_m11, nodes, weights);
  ASSERT_EQ(nodes.size(), node_count);
  ASSERT_EQ(weights.size(), node_count);

  for(int i = 0; i < node_count; ++i)
  {
    const int mirrored = node_count - 1 - i;
    EXPECT_NEAR(nodes[i], -nodes[mirrored], 1e-12);
    EXPECT_NEAR(weights[i], weights[mirrored], 1e-10);
  }

  const auto integrand = [omega](double x) { return 1.0 / (std::exp(pi * x / omega) + 1.0); };
  const double observed = integrate_rule_m11(nodes, weights, integrand);
  EXPECT_NEAR(observed, 1.0, 1e-12);
}

TEST(numerics_quadrature, rational_fejer_internal_matches_algorithm973_example7_6_table_scale)
{
  struct Case
  {
    int node_count;
    double paper_relative_error;
  };

  constexpr double pole = -2.5;
  const Case cases[] = {{2, 2.49e-03}, {4, 1.85e-06}, {8, 5.36e-12}, {12, 2.22e-17}};

  const auto integrand = [](double x) { return 1.0 / std::sqrt((x + 3.0) * (x + 2.0)); };
  const double reference = integrate_reference_m11(integrand);

  for(const auto& test_case : cases)
  {
    axom::Array<Complex> poles_m11;
    poles_m11.reserve(test_case.node_count - 1);
    for(int i = 0; i < test_case.node_count - 1; ++i)
    {
      poles_m11.emplace_back(pole, 0.0);
    }

    axom::Array<double> nodes;
    axom::Array<double> weights;
    compute_rational_fejer_rule_m11(poles_m11, nodes, weights);
    ASSERT_EQ(nodes.size(), test_case.node_count);
    ASSERT_EQ(weights.size(), test_case.node_count);

    const double observed = integrate_rule_m11(nodes, weights, integrand);
    const double relative_error = std::abs(observed - reference) / std::abs(reference);
    const double tolerance = axom::utilities::max(1e-12, 100.0 * test_case.paper_relative_error);
    EXPECT_LE(relative_error, tolerance) << "node_count=" << test_case.node_count;
  }
}

TEST(numerics_quadrature,
     rational_fejer_internal_matches_algorithm973_example7_7_outperforms_polynomial)
{
  constexpr int node_count = 10;
  const double inf = axom::numeric_limits<double>::infinity();

  axom::Array<Complex> rational_poles_m11 {{0.0, 0.2}, {0.0, -0.2}};
  while(static_cast<int>(rational_poles_m11.size()) < node_count - 1)
  {
    rational_poles_m11.emplace_back(inf, 0.0);
  }

  axom::Array<Complex> polynomial_poles_m11;
  polynomial_poles_m11.reserve(node_count - 1);
  for(int i = 0; i < node_count - 1; ++i)
  {
    polynomial_poles_m11.emplace_back(inf, 0.0);
  }

  const auto integrand = [](double x) { return std::exp(x) / (25.0 * x * x + 1.0); };
  const double reference = integrate_reference_m11(integrand);

  axom::Array<double> rational_nodes;
  axom::Array<double> rational_weights;
  compute_rational_fejer_rule_m11(rational_poles_m11, rational_nodes, rational_weights);

  axom::Array<double> polynomial_nodes;
  axom::Array<double> polynomial_weights;
  compute_rational_fejer_rule_m11(polynomial_poles_m11, polynomial_nodes, polynomial_weights);

  const double rational_error =
    std::abs(integrate_rule_m11(rational_nodes, rational_weights, integrand) - reference);
  const double polynomial_error =
    std::abs(integrate_rule_m11(polynomial_nodes, polynomial_weights, integrand) - reference);

  EXPECT_LT(rational_error, 1e-10);
  EXPECT_LT(rational_error, polynomial_error * 1e-4);
}
