// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core/numerics/quadrature.hpp"
#include "axom/core/numerics/internal/rational_quadrature.hpp"
#include "axom/core/utilities/Utilities.hpp"

#include <algorithm>
#include <cmath>

#include <complex>
#include <functional>

namespace
{
namespace numerics_internal = axom::numerics::internal;
using Complex = std::complex<double>;

const double pi = std::acos(-1.0);
axom::Array<Complex> map_interval_poles_m11_to_01(axom::ArrayView<const Complex> poles_m11)
{
  axom::Array<Complex> poles01;
  poles01.reserve(poles_m11.size());
  for(const auto& pole : poles_m11)
  {
    if(std::isinf(pole.real()))
    {
      poles01.emplace_back(axom::numeric_limits<double>::infinity(), 0.0);
    }
    else
    {
      poles01.push_back(0.5 * (pole + Complex {1.0, 0.0}));
    }
  }
  return poles01;
}

double integrate_rule(const axom::numerics::QuadratureRule& rule,
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

axom::Array<Complex> rotor_first_invalid_unique_poles_m11()
{
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

template <typename ExecSpace>
struct test_device_quadrature
{
  static void test()
  {
    const int npts = 15;

    // Create the rule with the proper allocator
    int allocID;
#if defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_GPU)

    // TODO QuadratureRule class needs to be ported for CUDA
    constexpr bool on_device = axom::execution_space<ExecSpace>::onDevice();

    allocID = on_device ? axom::getUmpireResourceAllocatorID(umpire::resource::Unified)
                        : axom::execution_space<ExecSpace>::allocatorID();
#else
    allocID = axom::execution_space<ExecSpace>::allocatorID();
#endif

    const auto rule = axom::numerics::get_gauss_legendre(npts, allocID);

    // Use the rule in a lambda to integrate the volume under std::sin(pi * x) on [0, 1]
    axom::ReduceSum<ExecSpace, double> quadrature_sum(0.0);
    axom::for_all<ExecSpace>(
      npts,
      AXOM_LAMBDA(axom::IndexType i) { quadrature_sum += rule.weight(i) * sin(rule.node(i)); });

    EXPECT_NEAR(quadrature_sum.get(), 0.459697694132, 1e-6);
  }
};

TEST(numerics_quadrature, get_nodes_seq) { test_device_quadrature<axom::SEQ_EXEC>::test(); }

#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
TEST(numerics_quadrature, get_nodes_omp) { test_device_quadrature<axom::OMP_EXEC>::test(); }
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
TEST(numerics_quadrature, get_nodes_cuda) { test_device_quadrature<axom::CUDA_EXEC<256>>::test(); }
#endif

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
TEST(numerics_quadrature, get_nodes_hip) { test_device_quadrature<axom::HIP_EXEC<256>>::test(); }
#endif

TEST(numerics_quadrature, rational_fejer_auto_adds_missing_complex_conjugate_pole)
{
  const double a = 1.4;
  const double b = 0.8;
  const Complex pole {a, b};

  const axom::Array<Complex> implicit_poles {pole};
  const axom::Array<Complex> explicit_poles {pole, std::conj(pole)};
  const auto implicit_rule = axom::numerics::get_rational_fejer(implicit_poles);
  const auto explicit_rule = axom::numerics::get_rational_fejer(explicit_poles);

  ASSERT_EQ(implicit_rule.getNumPoints(), explicit_rule.getNumPoints());
  for(int i = 0; i < implicit_rule.getNumPoints(); ++i)
  {
    EXPECT_NEAR(implicit_rule.node(i), explicit_rule.node(i), 1e-14);
    EXPECT_NEAR(implicit_rule.weight(i), explicit_rule.weight(i), 1e-14);
  }

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

TEST(numerics_quadrature, rotor_pole_family_stays_healthy_until_infinite_padding_reaches_32)
{
  const auto base_poles01 = map_interval_poles_m11_to_01(rotor_first_invalid_unique_poles_m11());

  auto summary_for_total_poles = [&base_poles01](int total_poles) {
    auto poles01 = base_poles01;
    pad_poles_with_infinities(poles01, total_poles);
    return compute_rule_diagnostics_summary(poles01);
  };

  const auto summary24 = summary_for_total_poles(24);
  const auto summary32 = summary_for_total_poles(32);
  const auto summary48 = summary_for_total_poles(48);
  const auto summary64 = summary_for_total_poles(64);

  EXPECT_GT(summary24.min_rational_chebyshev_weight, 0.0);
  EXPECT_GT(summary24.min_final_weight_m11, 0.0);
  EXPECT_NEAR(summary24.sum_final_weights_m11, 2.0, 1e-10);
  EXPECT_LT(summary24.max_abs_basis_coefficient, 10.0);

  // The first observed breakdown for this rotor pole family occurs when the
  // same three finite poles are padded to 32 total poles with infinities.
  EXPECT_GT(summary32.min_rational_chebyshev_weight, 0.0);
  EXPECT_LT(summary32.min_final_weight_m11, 0.0);
  EXPECT_NEAR(summary32.sum_final_weights_m11, 2.0, 1e-10);
  EXPECT_GT(summary32.max_abs_basis_coefficient, 10.0);

  EXPECT_GT(summary48.min_rational_chebyshev_weight, 0.0);
  EXPECT_LT(summary48.min_final_weight_m11, -1e6);
  EXPECT_GT(summary48.max_abs_basis_coefficient, 1e8);

  EXPECT_GT(summary64.min_rational_chebyshev_weight, 0.0);
  EXPECT_LT(summary64.min_final_weight_m11, -1e10);
  EXPECT_GT(std::abs(summary64.sum_final_weights_m11 - 2.0), 1e-2);
  EXPECT_GT(summary64.max_abs_basis_coefficient, 1e12);
}

TEST(numerics_quadrature, rotor_pole_family_shows_infinite_padding_drives_instability)
{
  const auto base_poles01 = map_interval_poles_m11_to_01(rotor_first_invalid_unique_poles_m11());

  auto make_summary = [](axom::Array<Complex> poles01) {
    return compute_rule_diagnostics_summary(poles01);
  };

  const auto unique_summary = make_summary(base_poles01);

  auto padded_poles01 = base_poles01;
  pad_poles_with_infinities(padded_poles01, 64);
  const auto padded_summary = make_summary(padded_poles01);

  auto repeated_poles01 = repeat_pole_sequence(base_poles01, 5);
  const auto repeated_summary = make_summary(repeated_poles01);

  auto repeated_padded_poles01 = repeated_poles01;
  pad_poles_with_infinities(repeated_padded_poles01, 64);
  const auto repeated_padded_summary = make_summary(repeated_padded_poles01);

  EXPECT_GT(unique_summary.min_rational_chebyshev_weight, 0.0);
  EXPECT_GT(unique_summary.min_final_weight_m11, 0.0);
  EXPECT_NEAR(unique_summary.sum_final_weights_m11, 2.0, 1e-10);

  EXPECT_GT(repeated_summary.min_rational_chebyshev_weight, 0.0);
  EXPECT_GT(repeated_summary.min_final_weight_m11, 0.0);
  EXPECT_NEAR(repeated_summary.sum_final_weights_m11, 2.0, 1e-10);

  EXPECT_GT(padded_summary.min_rational_chebyshev_weight, 0.0);
  EXPECT_LT(padded_summary.min_final_weight_m11, -1e6);
  EXPECT_GT(std::abs(padded_summary.sum_final_weights_m11 - 2.0), 1e-2);
  EXPECT_GT(padded_summary.max_abs_basis_coefficient, 1e12);

  EXPECT_GT(repeated_padded_summary.min_rational_chebyshev_weight, 0.0);
  EXPECT_LT(repeated_padded_summary.min_final_weight_m11, -1e6);
  EXPECT_GT(std::abs(repeated_padded_summary.sum_final_weights_m11 - 2.0), 1e-2);
  EXPECT_GT(repeated_padded_summary.max_abs_basis_coefficient, 1e12);

  EXPECT_LT(padded_summary.min_final_weight_m11, repeated_summary.min_final_weight_m11);
  EXPECT_LT(repeated_padded_summary.min_final_weight_m11, repeated_summary.min_final_weight_m11);
}

TEST(numerics_quadrature, rational_chebyshev_internal_matches_vandeun_pure_imaginary_poles_example)
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

  EXPECT_NEAR(weight_sum, pi, 1e-10);
}

TEST(numerics_quadrature,
     rational_chebyshev_internal_all_infinite_poles_reduces_to_chebyshev_first_kind)
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

TEST(numerics_quadrature, rational_chebyshev_internal_matches_vandeun_boundary_layer_application_poles)
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
  for(int i = 0; i < static_cast<int>(nodes.size()); ++i)
  {
    EXPECT_GT(weights[i], 0.0);
    weight_sum += weights[i];

    const int mirrored = static_cast<int>(nodes.size()) - 1 - i;
    EXPECT_NEAR(nodes[i], -nodes[mirrored], 1e-10);
    EXPECT_NEAR(weights[i], weights[mirrored], 1e-10);
  }

  EXPECT_NEAR(weight_sum, pi, 1e-10);
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
