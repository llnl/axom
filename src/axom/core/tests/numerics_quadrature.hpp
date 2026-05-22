// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core/Array.hpp"
#include "axom/core/numerics/quadrature.hpp"
#include "axom/core/utilities/Utilities.hpp"

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

namespace
{
template <typename RuleGetter>
void check_polynomial_exactness(RuleGetter&& getRule, int maxNpts)
{
  for(int npts = 1; npts <= maxNpts; ++npts)
  {
    const auto rule = getRule(npts);
    const int exactDegree = npts - 1 + npts % 2;
    axom::Array<double> coeffs(exactDegree + 1, exactDegree + 1);

    for(int j = 0; j <= exactDegree; ++j)
    {
      coeffs[j] = axom::utilities::random_real(-1.0, 1.0, 1000 * npts + j);
    }

    double analyticResult = 0.0;
    for(int j = 0; j <= exactDegree; ++j)
    {
      analyticResult += coeffs[j] / (j + 1);
    }

    auto evalPolynomial = [&coeffs, exactDegree](double x) {
      double result = coeffs[exactDegree];
      for(int i = exactDegree - 1; i >= 0; --i)
      {
        result = result * x + coeffs[i];
      }
      return result;
    };

    double quadratureResult = 0.0;
    double weightSum = 0.0;
    for(int j = 0; j < npts; ++j)
    {
      quadratureResult += rule.weight(j) * evalPolynomial(rule.node(j));
      weightSum += rule.weight(j);
      if(j > 0)
      {
        EXPECT_GT(rule.node(j), rule.node(j - 1));
      }
    }

    EXPECT_NEAR(quadratureResult, analyticResult, 1e-12);
    EXPECT_NEAR(weightSum, 1.0, 1e-12);
  }
}
}  // namespace

TEST(numerics_quadrature, open_uniform_small_rules)
{
  auto rule = axom::numerics::get_open_uniform(1);
  ASSERT_EQ(rule.getNumPoints(), 1);
  EXPECT_DOUBLE_EQ(rule.node(0), 0.5);
  EXPECT_DOUBLE_EQ(rule.weight(0), 1.0);

  rule = axom::numerics::get_open_uniform(3);
  ASSERT_EQ(rule.getNumPoints(), 3);
  EXPECT_DOUBLE_EQ(rule.node(0), 0.25);
  EXPECT_DOUBLE_EQ(rule.node(1), 0.5);
  EXPECT_DOUBLE_EQ(rule.node(2), 0.75);
  EXPECT_DOUBLE_EQ(rule.weight(0), 2.0 / 3.0);
  EXPECT_DOUBLE_EQ(rule.weight(1), -1.0 / 3.0);
  EXPECT_DOUBLE_EQ(rule.weight(2), 2.0 / 3.0);
}

TEST(numerics_quadrature, closed_uniform_small_rules)
{
  auto rule = axom::numerics::get_closed_uniform(1);
  ASSERT_EQ(rule.getNumPoints(), 1);
  EXPECT_DOUBLE_EQ(rule.node(0), 0.5);
  EXPECT_DOUBLE_EQ(rule.weight(0), 1.0);

  rule = axom::numerics::get_closed_uniform(3);
  ASSERT_EQ(rule.getNumPoints(), 3);
  EXPECT_DOUBLE_EQ(rule.node(0), 0.0);
  EXPECT_DOUBLE_EQ(rule.node(1), 0.5);
  EXPECT_DOUBLE_EQ(rule.node(2), 1.0);
  EXPECT_DOUBLE_EQ(rule.weight(0), 1.0 / 6.0);
  EXPECT_DOUBLE_EQ(rule.weight(1), 4.0 / 6.0);
  EXPECT_DOUBLE_EQ(rule.weight(2), 1.0 / 6.0);
}

TEST(numerics_quadrature, quadrature_type_dispatch)
{
  using axom::numerics::QuadratureType;

  auto rule = axom::numerics::get_quadrature_rule(QuadratureType::Invalid, 2);
  EXPECT_DOUBLE_EQ(rule.node(0), 0.21132486540518711775);
  EXPECT_DOUBLE_EQ(rule.node(1), 0.78867513459481288225);

  rule = axom::numerics::get_quadrature_rule(QuadratureType::GaussLegendre, 2);
  EXPECT_DOUBLE_EQ(rule.node(0), 0.21132486540518711775);
  EXPECT_DOUBLE_EQ(rule.node(1), 0.78867513459481288225);

  rule = axom::numerics::get_quadrature_rule(QuadratureType::OpenUniform, 3);
  EXPECT_DOUBLE_EQ(rule.node(0), 0.25);
  EXPECT_DOUBLE_EQ(rule.node(1), 0.5);
  EXPECT_DOUBLE_EQ(rule.node(2), 0.75);

  rule = axom::numerics::get_quadrature_rule(QuadratureType::ClosedUniform, 3);
  EXPECT_DOUBLE_EQ(rule.node(0), 0.0);
  EXPECT_DOUBLE_EQ(rule.node(1), 0.5);
  EXPECT_DOUBLE_EQ(rule.node(2), 1.0);
}

TEST(numerics_quadrature, open_uniform_exactness)
{
  check_polynomial_exactness([](int npts) { return axom::numerics::get_open_uniform(npts); }, 10);
}

TEST(numerics_quadrature, closed_uniform_exactness)
{
  check_polynomial_exactness([](int npts) { return axom::numerics::get_closed_uniform(npts); }, 10);
}
