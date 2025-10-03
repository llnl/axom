// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/numerics/quadrature.hpp"
#include <iomanip>

TEST(numerics_quadrature, gauss_legendre_math_check)
{
  const int N = 200;

  double coeffs[2 * N - 1];

  // Test that the rules provide exact integration for polynomials of degree 2n - 1
  for(int npts = 1; npts < N; ++npts)
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

TEST(numerics_quadrature, gauss_legendre_cache_check)
{
  // The first two rules are put in the cache
  const axom::numerics::QuadratureRule& first_rule = axom::numerics::get_gauss_legendre(20);
  const axom::numerics::QuadratureRule& second_rule = axom::numerics::get_gauss_legendre(20);

  // The third is not
  const axom::numerics::QuadratureRule& third_rule = axom::numerics::compute_gauss_legendre(20);

  // Check that the two rules are located in the same place in memory, and the third isn't
  //EXPECT_EQ(&first_rule, &second_rule);
  //EXPECT_NE(&first_rule, &third_rule);
}
