// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/numerics/quadrature.hpp"
#include <iomanip>

TEST(numerics_quadrature, generate_rules)
{
  const int N = 50;

  // Define a collection of random coefficients for a polynomial
  double coeffs[2 * N - 1];
  for(int j = 0; j < 2 * N - 1; ++j)
  {
    coeffs[j] = axom::utilities::random_real(-5.0, 5.0);
  }

  // Test that the rules provide exact integration for polynomials of degree 2n - 1
  for(int npts = 1; npts < N; ++npts)
  {
    int degree = 2 * npts - 1;

    // Evaluate the area under the curve from 0 to 1 analytically
    double analytic_result = 0.0;
    for(int j = 0; j < degree; ++j)
    {
      analytic_result += coeffs[j] / (j + 1);
    }

    // Evaluate using the quadrature rule
    auto rule = axom::numerics::compute_gauss_legendre(npts);

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

    std::cout << std::setprecision(20);
    std::cout << "=====" << std::endl;
    std::cout << analytic_result << std::endl;
    std::cout << quadrature_result << std::endl;
    std::cout << "=====" << std::endl;
  }
}
