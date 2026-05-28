// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! \file core_quadrature.cpp
 *  \brief This example demonstrates Axom's numerical quadrature capabilities.
 */

/* This example code contains snippets used in the Core Sphinx documentation.
 * They begin and end with comments such as
 *
 * gauss_legendre_start
 * gauss_legendre_end
 *
 * each prepended with an underscore.
 */

// Axom includes
#include "axom/core/numerics/quadrature.hpp"
#include "axom/core/Array.hpp"
#include "axom/fmt.hpp"

// C/C++ includes
#include <iostream>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>

//------------------------------------------------------------------------------
// Helper function to demonstrate quadrature usage
//------------------------------------------------------------------------------
// _quadrature_rule_start
double integrate(const axom::numerics::QuadratureRuleView& rule, std::function<double(double)> f)
{
  double result = 0.0;
  for(int i = 0; i < rule.getNumPoints(); ++i)
  {
    result += rule.weight(i) * f(rule.node(i));
  }
  return result;
}
// _quadrature_rule_end

//------------------------------------------------------------------------------
// Gauss-Legendre quadrature examples
//------------------------------------------------------------------------------
void demoGaussLegendre()
{
  std::cout << axom::fmt::format("\n{:=^40}\n", " Gauss-Legendre Quadrature Examples ");
  std::cout << "\n";
  {
    // _gauss_legendre_basic_start
    // Get a Gauss-Legendre quadrature rule with 5 points
    auto rule = axom::numerics::get_gauss_legendre(5);

    std::cout << axom::fmt::format("Gauss-Legendre rule with {} points:\n", rule.getNumPoints());

    // Access nodes and weights
    for(int i = 0; i < rule.getNumPoints(); ++i)
    {
      std::cout << axom::fmt::format("  node[{}] = {:.6f}, weight[{}] = {:.6f}\n",
                                     i,
                                     rule.node(i),
                                     i,
                                     rule.weight(i));
    }
    // _gauss_legendre_basic_end
  }
  std::cout << "\n";

  {
    // _gauss_legendre_integrate_start
    // Get a Gauss-Legendre quadrature rule with 2 points
    auto rule = axom::numerics::get_gauss_legendre(2);

    // Integrate f(x) = x^3 from 0 to 1
    // Exact answer: 1/4 = 0.25
    auto cubic = [](double x) { return x * x * x; };

    double integral = 0.0;
    for(int i = 0; i < rule.getNumPoints(); ++i)
    {
      integral += rule.weight(i) * cubic(rule.node(i));
    }

    std::cout << axom::fmt::format("Integrate f(x) = x^3 from 0 to 1:\n");
    std::cout << axom::fmt::format("  Result: {:.6f}\n", integral);
    std::cout << axom::fmt::format("  Exact:  0.25\n");
    std::cout << axom::fmt::format("  Error:  {:.2e}\n", std::abs(integral - 0.25));
    // _gauss_legendre_integrate_end
  }
  std::cout << "\n";

  {
    // _gauss_legendre_exactness_start
    // Demonstrate polynomial exactness
    // A 5-point rule can exactly integrate polynomials up to degree 2*5-1 = 9
    auto rule = axom::numerics::get_gauss_legendre(5);

    auto polynomial_degree_7 = [](double x) {
      // x^7 + 2x^5 - 3x^3 + x
      return std::pow(x, 7) + 2 * std::pow(x, 5) - 3 * std::pow(x, 3) + x;
    };

    // Analytic integral from 0 to 1
    constexpr double exact = 1. / 8. + 2. / 6. - 3. / 4. + 1. / 2.;

    const double result = integrate(rule, polynomial_degree_7);

    std::cout << axom::fmt::format("Integrate polynomial of degree 7:\n");
    std::cout << axom::fmt::format("  Result: {:.15f}\n", result);
    std::cout << axom::fmt::format("  Exact:  {:.15f}\n", exact);
    std::cout << axom::fmt::format("  Error:  {:.2e} (should be ~machine epsilon)\n",
                                   std::abs(result - exact));
    // _gauss_legendre_exactness_end
  }
  std::cout << "\n";

  {
    // _gauss_legendre_smooth_start
    // Integrate a smooth function: sin(pi*x) from 0 to 1
    // Exact: 2/pi ≈ 0.636620
    auto sine_func = [](double x) { return std::sin(M_PI * x); };
    double exact_sine = 2.0 / M_PI;

    std::cout << axom::fmt::format("Integrate sin(pi*x) from 0 to 1 with varying orders:\n");
    for(int npts : {3, 5, 10, 20})
    {
      auto r = axom::numerics::get_gauss_legendre(npts);
      double res = integrate(r, sine_func);
      double error = std::abs(res - exact_sine);
      std::cout << axom::fmt::format("  {:2} points: error = {:.3e}\n", npts, error);
    }
    // _gauss_legendre_smooth_end
  }
}

//------------------------------------------------------------------------------
// Rational Fejer quadrature examples
//------------------------------------------------------------------------------
void demoRationalFejer()
{
  std::cout << axom::fmt::format("\n{:=^40}\n", " Rational Fejer Quadrature Examples ");
  std::cout << "\n";

  {
    // _rational_fejer_singularity_start
    // Integrate f(x) = 1/(x - 1.5) from 0 to 1
    // This has a singularity at x = 1.5 (outside the domain)
    // Exact integral: log(|1-1.5|/|0-1.5|) = log(0.5/1.5) ≈ -1.0986

    // Define poles at the singularity location
    axom::Array<std::complex<double>> poles = {std::complex<double>(1.5, 0.0)};

    // Get rational Fejer rule
    auto rule = axom::numerics::get_rational_fejer(poles);

    std::cout << axom::fmt::format("Rational Fejer rule with pole at x = 1.5:\n");
    std::cout << axom::fmt::format("  Number of points: {}\n", rule.getNumPoints());

    // Integrate
    auto singular_func = [](double x) { return 1.0 / (x - 1.5); };
    double result = integrate(rule, singular_func);
    double exact = std::log(0.5 / 1.5);

    std::cout << axom::fmt::format("  Result: {:.6f}\n", result);
    std::cout << axom::fmt::format("  Exact:  {:.6f}\n", exact);
    std::cout << axom::fmt::format("  Error:  {:.2e}\n", std::abs(result - exact));
    // _rational_fejer_singularity_end
  }
  std::cout << "\n";

  {
    // _rational_fejer_vs_gauss_start
    // Compare Gauss-Legendre vs Rational Fejer for near-singularity
    std::cout << axom::fmt::format("Comparison: Gauss-Legendre vs Rational Fejer\n");
    std::cout << axom::fmt::format("Function: f(x) = 1/(x - 1.2) on [0, 1]\n");
    std::cout << axom::fmt::format("Singularity at x = 1.2 (close to domain boundary)\n\n");

    auto near_singular = [](double x) { return 1.0 / (x - 1.2); };
    double exact_near = std::log(0.2 / 1.2);

    // Try Gauss-Legendre with increasing points
    std::cout << axom::fmt::format("Gauss-Legendre:\n");
    for(int npts : {5, 10, 20, 50})
    {
      auto gl_rule = axom::numerics::get_gauss_legendre(npts);
      double res = integrate(gl_rule, near_singular);
      double error = std::abs(res - exact_near);
      std::cout << axom::fmt::format("  {:2} points: error = {:.3e}\n", npts, error);
    }

    // Rational Fejer with pole at singularity
    std::cout << axom::fmt::format("\nRational Fejer (pole at x = 1.2):\n");
    axom::Array<std::complex<double>> poles_near = {std::complex<double>(1.2, 0.0)};
    auto rf_rule = axom::numerics::get_rational_fejer(poles_near);
    double rf_result = integrate(rf_rule, near_singular);
    double rf_error = std::abs(rf_result - exact_near);
    std::cout << axom::fmt::format("  {:2} points: error = {:.3e}\n", rf_rule.getNumPoints(), rf_error);
    std::cout << axom::fmt::format("  (Achieves machine precision with only {} points!)\n",
                                   rf_rule.getNumPoints());
    // _rational_fejer_vs_gauss_end
  }
  std::cout << "\n";

  {
    // _rational_fejer_corner_start
    // Integrate sqrt(x) from 0 to 1
    // Has infinite derivative at x = 0 (corner singularity)
    // Exact: 2/3 ≈ 0.666667

    std::cout << axom::fmt::format("Corner singularity: integrate sqrt(x) from 0 to 1\n");
    std::cout << axom::fmt::format("Exact integral: 2/3 = 0.666667\n\n");

    auto sqrt_func = [](double x) { return std::sqrt(x); };
    double exact_sqrt = 2.0 / 3.0;

    // For algebraic endpoint singularities, nearby poles alone are often too crude.
    // Mixing a couple of nearby poles with infinity poles keeps some polynomial
    // Fejer character and works much better for sqrt(x).
    constexpr double inf = std::numeric_limits<double>::infinity();
    axom::Array<std::complex<double>> poles_corner = {std::complex<double>(-0.1, 0.0),
                                                      std::complex<double>(-0.2, 0.0),
                                                      std::complex<double>(inf, 0.0),
                                                      std::complex<double>(inf, 0.0)};

    auto corner_rule = axom::numerics::get_rational_fejer(poles_corner);
    double corner_result = integrate(corner_rule, sqrt_func);
    double corner_error = std::abs(corner_result - exact_sqrt);

    std::cout << axom::fmt::format("Rational Fejer with poles at x = {{-0.1, -0.2, inf, inf}}:\n");
    std::cout << axom::fmt::format("  {} points\n", corner_rule.getNumPoints());
    std::cout << axom::fmt::format("  Result: {:.6f}\n", corner_result);
    std::cout << axom::fmt::format("  Error:  {:.3e}\n", corner_error);
    // _rational_fejer_corner_end
  }
  std::cout << "\n";

  {
    // _rational_fejer_complex_poles_start
    // Integrate with complex poles
    // f(x) = 1/((x - 1.5)^2 + 0.04)
    // Has complex poles at x = 1.5 + .2i and x = 1.5 - .2i

    std::cout << axom::fmt::format("Complex poles: f(x) = 1/((x - 1.5)^2 + 0.04)\n");
    std::cout << axom::fmt::format("Poles at x = {{1.5 + 0.2i, 1.5 - 0.2i}}\n\n");

    // Only need to specify one pole; conjugate is auto-added
    axom::Array<std::complex<double>> poles_complex = {std::complex<double>(1.5, 0.2)};

    auto complex_rule = axom::numerics::get_rational_fejer(poles_complex);

    std::cout << axom::fmt::format("Specified 1 pole, rule has {} points\n",
                                   complex_rule.getNumPoints());
    std::cout << axom::fmt::format(
      "(Complex poles come in conjugate pairs: 1 input --> 2 canonical --> 3 points)\n\n");

    auto complex_func = [](double x) {
      double denom = (x - 1.5) * (x - 1.5) + 0.04;
      return 1.0 / denom;
    };

    double complex_result = integrate(complex_rule, complex_func);
    std::cout << axom::fmt::format("  Integral: {:.6f}\n", complex_result);
    // _rational_fejer_complex_poles_end
  }
  std::cout << "\n";

  {
    // _rational_fejer_multiple_poles_start
    // Multiple singularities
    // f(x) = 1/((x - 1.3)(x - 1.7))
    std::cout << axom::fmt::format("Multiple singularities: f(x) = 1/((x - 1.3)(x - 1.7))\n");
    std::cout << axom::fmt::format("Poles at x = 1.3 and x = 1.7\n\n");

    axom::Array<std::complex<double>> poles_multiple = {std::complex<double>(1.3, 0.0),
                                                        std::complex<double>(1.7, 0.0)};

    auto multi_rule = axom::numerics::get_rational_fejer(poles_multiple);

    auto multi_func = [](double x) { return 1.0 / ((x - 1.3) * (x - 1.7)); };

    double multi_result = integrate(multi_rule, multi_func);

    // Analytic integral: (1/(1.7-1.3)) * log(|(0-1.3)(1-1.7)|/|(1-1.3)(0-1.7)|)
    double exact_multi = (1.0 / 0.4) * std::log((1.3 * 0.7) / (0.3 * 1.7));

    std::cout << axom::fmt::format("  {} points\n", multi_rule.getNumPoints());
    std::cout << axom::fmt::format("  Result: {:.6f}\n", multi_result);
    std::cout << axom::fmt::format("  Exact:  {:.6f}\n", exact_multi);
    std::cout << axom::fmt::format("  Error:  {:.3e}\n", std::abs(multi_result - exact_multi));
    // _rational_fejer_multiple_poles_end
  }
}

//------------------------------------------------------------------------------
// Advanced usage examples
//------------------------------------------------------------------------------
void demoAdvancedUsage()
{
  std::cout << axom::fmt::format("\n{:=^40}\n", " Advanced Quadrature Usage ");
  std::cout << "\n";

  // _quadrature_rule_copy_start
  // Demonstrating QuadratureRule vs QuadratureRuleView

  // Get a cached view (lightweight, non-owning)
  axom::Array<std::complex<double>> poles = {std::complex<double>(1.5, 0.0)};
  auto view = axom::numerics::get_rational_fejer(poles);

  std::cout << axom::fmt::format("QuadratureRuleView from cache:\n");
  std::cout << axom::fmt::format("  Points: {}\n", view.getNumPoints());

  // Create an owned copy if needed beyond immediate use
  axom::numerics::QuadratureRule owned = view.copy();

  std::cout << axom::fmt::format("\nOwned QuadratureRule from copy:\n");
  std::cout << axom::fmt::format("  Points: {}\n", owned.getNumPoints());
  std::cout << axom::fmt::format("  (Data is now independent of cache)\n");

  // Get view back from owned rule
  auto view_from_owned = owned.view();
  std::cout << axom::fmt::format("\nView from owned rule:\n");
  std::cout << axom::fmt::format("  Points: {}\n", view_from_owned.getNumPoints());
  // _quadrature_rule_copy_end

  std::cout << "\n";

  // _quadrature_compute_vs_cached_start
  // compute_rational_fejer_data vs get_rational_fejer

  std::cout << axom::fmt::format("Two ways to get quadrature rules:\n\n");

  // Method 1: Compute directly (no caching)
  axom::Array<double> nodes, weights;
  axom::numerics::compute_rational_fejer_data(poles, nodes, weights);

  std::cout << axom::fmt::format("1. compute_rational_fejer_data() - no caching:\n");
  std::cout << axom::fmt::format("   Returns nodes/weights arrays directly\n");
  std::cout << axom::fmt::format("   {} nodes, {} weights\n\n", nodes.size(), weights.size());

  // Method 2: Get cached (reuses previous computation)
  auto cached_rule = axom::numerics::get_rational_fejer(poles);

  std::cout << axom::fmt::format("2. get_rational_fejer() - uses LRU cache:\n");
  std::cout << axom::fmt::format("   Returns QuadratureRuleView over cached data\n");
  std::cout << axom::fmt::format("   {} points\n", cached_rule.getNumPoints());
  std::cout << axom::fmt::format("   Subsequent calls with same poles are fast (cache hit)\n");
  // _quadrature_compute_vs_cached_end

  std::cout << "\n";

  // _quadrature_infinite_poles_start
  // Using infinite poles (polynomial limit)

  std::cout << axom::fmt::format("Infinite poles (polynomial Fejer limit):\n\n");

  // All poles at infinity → standard Fejer (Chebyshev-based) rule
  constexpr double inf = std::numeric_limits<double>::infinity();
  axom::Array<std::complex<double>> infinite_poles = {std::complex<double>(inf, 0.0),
                                                      std::complex<double>(inf, 0.0),
                                                      std::complex<double>(inf, 0.0)};

  auto infinite_rule = axom::numerics::get_rational_fejer(infinite_poles);

  std::cout << axom::fmt::format("  3 infinite poles → {} point rule\n",
                                 infinite_rule.getNumPoints());
  std::cout << axom::fmt::format("  This is equivalent to a polynomial Fejer rule\n");
  std::cout << axom::fmt::format("  Good for smooth functions without specific pole structure\n");
  // _quadrature_infinite_poles_end
}

//------------------------------------------------------------------------------
// Main function
//------------------------------------------------------------------------------
int main()
{
  demoGaussLegendre();
  demoRationalFejer();
  demoAdvancedUsage();

  std::cout << axom::fmt::format("\n{:=^40}\n", " All quadrature examples completed! ");
  std::cout << "\n";

  return 0;
}
