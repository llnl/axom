// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
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

TEST(numerics_quadrature, gauss_jacobi_math_check)
{
  const int N = 10;
  // Test with a few alpha, beta pairs
  std::vector<std::pair<double, double>> params = {{0.0, 0.0},
                                                   {0.5, 0.5},
                                                   {-0.5, -0.5},
                                                   {1.0, 2.0},
                                                   {0.5, -0.2}};

  for(auto p : params)
  {
    double alpha = p.first;
    double beta = p.second;

    for(int npts = 1; npts <= N; ++npts)
    {
      auto rule = axom::numerics::get_gauss_jacobi(npts, alpha, beta);
      int max_degree = 2 * npts - 1;

      for(int k = 0; k <= max_degree; ++k)
      {
        // Compute analytic integral of x^k * (1-x)^alpha * (1+x)^beta on [-1, 1]
        // = 2^(alpha+beta+1) * Integral_0^1 (1-t)^alpha * t^beta * (2t-1)^k dt

        double analytic = 0.0;
        // Expand (2t-1)^k using binomial theorem
        // (2t-1)^k = sum_{j=0}^k C(k, j) * (2t)^j * (-1)^{k-j}
        // Integral term: 2^j * (-1)^{k-j} * Integral_0^1 (1-t)^alpha * t^{beta+j} dt
        // Integral is Beta(alpha+1, beta+j+1)
        // Beta(x,y) = Gamma(x)Gamma(y)/Gamma(x+y)

        double prefactor = std::pow(2.0, alpha + beta + 1.0);

        for(int j = 0; j <= k; ++j)
        {
          double binom = 1.0;  // C(k,j)
          for(int m = 1; m <= j; ++m)
          {
            binom = binom * (k - m + 1) / m;
          }

          double term_integral = std::exp(std::lgamma(alpha + 1.0) + std::lgamma(beta + j + 1.0) -
                                          std::lgamma(alpha + beta + j + 2.0));

          analytic += binom * std::pow(2.0, j) * std::pow(-1.0, k - j) * term_integral;
        }
        analytic *= prefactor;

        // Quadrature sum
        double quad_sum = 0.0;
        for(int i = 0; i < npts; ++i)
        {
          double x_node = rule.node(i);   // in [0, 1]
          double z = 2.0 * x_node - 1.0;  // in [-1, 1]
          quad_sum += rule.weight(i) * std::pow(z, k);
        }

        // increase tolerance as the number of points (N) increases
        double tol = std::pow(10.0, -12.0 + static_cast<double>(npts) / 2.0);
        EXPECT_NEAR(quad_sum, analytic, tol) << "Failed for npts=" << npts << ", alpha=" << alpha
                                             << ", beta=" << beta << ", degree=" << k;
      }
    }
  }
}

TEST(numerics_quadrature, jacobi_legendre_consistency)
{
  int npts = 5;
  auto leg_rule = axom::numerics::get_gauss_legendre(npts);
  auto jac_rule = axom::numerics::get_gauss_jacobi(npts, 0.0, 0.0);

  for(int i = 0; i < npts; ++i)
  {
    EXPECT_NEAR(leg_rule.node(i), jac_rule.node(i), 1e-12);
    // Legendre rule weights sum to 1 (normalized for [0,1])
    // Jacobi rule weights sum to 2 (integral of 1 on [-1,1])
    EXPECT_NEAR(leg_rule.weight(i) * 2.0, jac_rule.weight(i), 1e-12);
  }
}

template <typename ExecSpace>
struct test_device_quadrature
{
  static void test()
  {
    const int npts = 15;

    // Create the rule with the proper allocator
    const auto rule =
      axom::numerics::get_gauss_legendre(npts, axom::execution_space<ExecSpace>::allocatorID());

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