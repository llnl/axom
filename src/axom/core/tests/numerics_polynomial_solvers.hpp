// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/core/Array.hpp"
#include "axom/core/numerics/polynomial_solvers.hpp"
#include "axom/core/utilities/Utilities.hpp"  // for isNearlyEqual()

// Google Test include
#include "gtest/gtest.h"

int count_mismatches(const axom::Array<double>& standard,
                     const axom::Array<double>& test,
                     int n,
                     double thresh = 1e-8);

int count_mismatches(const axom::Array<double>& standard,
                     const axom::Array<double>& test,
                     int n,
                     double thresh)
{
  int mcount = 0;

  for(int i = 0; i < n; ++i)
  {
    if(!axom::utilities::isNearlyEqual(standard[i], test[i], thresh))
    {
      ++mcount;
    }
  }

  return mcount;
}

TEST(numerics_polynomial_solvers, solve_linear)
{
  int n = 1;
  int rc = -1;

  // In these tests, we are solving ax + b = 0, so
  // coeff[1] = a, coeff[0] = b.
  {
    SCOPED_TRACE("Line 1 through origin.");
    axom::Array<double> coeff {0., 1.};
    axom::Array<double> roots {0.};
    axom::Array<double> expected {0.};
    rc = axom::numerics::solve_linear(coeff.view(), roots.view(), n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(count_mismatches(expected, roots, 1), 0);
  }

  {
    SCOPED_TRACE("Line 2 through origin.");
    axom::Array<double> coeff {0., 18.};
    axom::Array<double> roots {0.};
    axom::Array<double> expected {0.};
    rc = axom::numerics::solve_linear(coeff.view(), roots.view(), n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(count_mismatches(expected, roots, 1), 0);
  }

  {
    SCOPED_TRACE("Off origin 1");
    axom::Array<double> coeff {-1., 0.5};
    axom::Array<double> roots {0.};
    axom::Array<double> expected {2.};
    rc = axom::numerics::solve_linear(coeff.view(), roots.view(), n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(count_mismatches(expected, roots, 1), 0);
  }

  {
    SCOPED_TRACE("Off origin 2");
    axom::Array<double> coeff {0.5, -1};
    axom::Array<double> roots {0.};
    axom::Array<double> expected {.5};
    rc = axom::numerics::solve_linear(coeff.view(), roots.view(), n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(count_mismatches(expected, roots, 1), 0);
  }

  {
    SCOPED_TRACE("X-axis");
    axom::Array<double> coeff {0., 0.};
    axom::Array<double> roots {0.};
    axom::Array<double> expected {0.};
    rc = axom::numerics::solve_linear(coeff.view(), roots.view(), n);
    // rc == 0 because there are real solutions
    EXPECT_EQ(rc, 0);
    // n == -1 because there are infinitely many solutions
    EXPECT_EQ(n, -1);
    EXPECT_EQ(count_mismatches(expected, roots, 1), 0);
  }
}

TEST(numerics_polynomial_solvers, solve_quadratic)
{
  int n = 2;
  int rc = -1;

  // In these tests, we are solving ax^2 + bx + c = 0, so
  // coeff[2] = a, coeff[1] = b, coeff[0] = c.
  {
    // y = (x + 2.3)(x + 2.3)
    SCOPED_TRACE("Double root at x = -2.3");
    axom::Array<double> coeff {5.29, 4.6, 1.};
    axom::Array<double> roots {0., 0.};
    axom::Array<double> expected {-2.3, -2.3};
    rc = axom::numerics::solve_quadratic(coeff.view(), roots.view(), n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(count_mismatches(expected, roots, 2), 0);
  }

  {
    // y = (-x + 1.5)(x - 1.5)
    SCOPED_TRACE("Double root at x = 1.5 (opening down)");
    axom::Array<double> coeff {-2.25, 3., -1.};
    axom::Array<double> roots {0., 0.};
    axom::Array<double> expected {1.5, 1.5};
    rc = axom::numerics::solve_quadratic(coeff.view(), roots.view(), n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(count_mismatches(expected, roots, 2), 0);
  }

  {
    // y = 3.2(x + 0.7)(x - 2)
    SCOPED_TRACE("Roots at -0.7 and 2");
    axom::Array<double> coeff {-4.48, -4.16, 3.2};
    axom::Array<double> roots {0., 0.};
    axom::Array<double> expected {2, -0.7};
    rc = axom::numerics::solve_quadratic(coeff.view(), roots.view(), n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 2);
    EXPECT_EQ(count_mismatches(expected, roots, 2), 0);
  }

  {
    // y = 0.1x^2 + 0.2x + 6
    SCOPED_TRACE("No real roots (opening up)");
    axom::Array<double> coeff {6, 0.2, 0.1};
    axom::Array<double> roots {0., 0.};
    axom::Array<double> expected {0., 0.};
    rc = axom::numerics::solve_quadratic(coeff.view(), roots.view(), n);
    EXPECT_EQ(rc, -1);
    EXPECT_EQ(n, 0);
    EXPECT_EQ(count_mismatches(expected, roots, 2), 0);
  }

  {
    // y = -5x^2 + 0.2x -20
    SCOPED_TRACE("No real roots (opening up)");
    axom::Array<double> coeff {6, 0.2, 0.1};
    axom::Array<double> roots {0., 0.};
    axom::Array<double> expected {0., 0.};
    rc = axom::numerics::solve_quadratic(coeff.view(), roots.view(), n);
    EXPECT_EQ(rc, -1);
    EXPECT_EQ(n, 0);
    EXPECT_EQ(count_mismatches(expected, roots, 2), 0);
  }
}

TEST(numerics_polynomial_solvers, solve_cubic)
{
  int n = 3;
  int rc = -1;

  // In these tests, we are solving ax^3 + bx^2 + cx + d = 0, so
  // coeff[3] = a, coeff[2] = b, coeff[1] = c, coeff[0] = d.
  {
    // y = (x - 1.2)^3 = 1.0 x^3 - 3.6 x^2 + 4.32 x - 1.728
    // Repeated roots have inherently poor numerical conditioning in polynomial solvers.
    // We use a looser tolerance here (5e-5 vs default 1e-8) because the discriminant
    // is nearly zero and floating-point error accumulates during the triple root calculation.
    SCOPED_TRACE("Triple root at x = 1.2 (loose tolerance due to repeated root conditioning)");
    axom::Array<double> coeff {-1.728, 4.32, -3.6, 1};
    axom::Array<double> roots {0., 0., 0.};
    axom::Array<double> expected {1.2, 1.2, 1.2};
    rc = axom::numerics::solve_cubic(coeff.view(), roots.view(), n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(count_mismatches(expected, roots, 3, 5.0e-5), 0);
  }

  {
    // y = -x^3 - x + 10
    SCOPED_TRACE("Single real root at x = 2");
    axom::Array<double> coeff {10, -1, 0, -1};
    axom::Array<double> roots {0., 0., 0.};
    axom::Array<double> expected {2., 0., 0.};
    rc = axom::numerics::solve_cubic(coeff.view(), roots.view(), n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(count_mismatches(expected, roots, 3), 0);
  }

  {
    // y = x^3 + x^2 - 8*x - 12
    SCOPED_TRACE("Two real roots, at x = 3, twice at x = -2");
    axom::Array<double> coeff {-12, -8, 1, 1};
    axom::Array<double> roots {0., 0., 0.};
    axom::Array<double> expected {3, -2, -2};
    rc = axom::numerics::solve_cubic(coeff.view(), roots.view(), n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 2);
    EXPECT_EQ(count_mismatches(expected, roots, 3), 0);
  }

  {
    // y = (x + 0.8)(-x + 1)(x - 8) = -3x^3 + 24.6x^2 - 2.4x - 19.2
    SCOPED_TRACE("Three real roots, at x = -0.8, 1, 8");
    axom::Array<double> coeff {-19.2, -2.4, 24.6, -3};
    axom::Array<double> roots {0., 0., 0.};
    axom::Array<double> expected {8, 1, -0.8};
    rc = axom::numerics::solve_cubic(coeff.view(), roots.view(), n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 3);
    EXPECT_EQ(count_mismatches(expected, roots, 3), 0);
  }

  {
    // y = 4.3(x + 38)(x + 1)(x - 0.001)
    //   = 4.3x^3 + 167.6957x^2 + 163.2323x - 0.1634
    SCOPED_TRACE("Three real roots, at x = -38, -1, 0.001");
    axom::Array<double> coeff {-0.1634, 163.2323, 167.6957, 4.3};
    axom::Array<double> roots {0., 0., 0.};
    axom::Array<double> expected {0.001, -1, -38};
    rc = axom::numerics::solve_cubic(coeff.view(), roots.view(), n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 3);
    EXPECT_EQ(count_mismatches(expected, roots, 3), 0);
  }
}

TEST(numerics_polynomial_solvers, bernstein_to_monomial)
{
  {
    SCOPED_TRACE("Quadratic Bernstein conversion");
    const axom::Array<double> bernstein_coeffs {1., 2., 4.};
    const axom::Array<double> expected_monomial {1., 2., 1.};
    const auto monomial_coeffs = axom::numerics::bernstein_to_monomial(bernstein_coeffs.view());

    ASSERT_EQ(expected_monomial.size(), monomial_coeffs.size());
    for(axom::IndexType i = 0; i < expected_monomial.size(); ++i)
    {
      EXPECT_DOUBLE_EQ(expected_monomial[i], monomial_coeffs[i]);
    }
  }

  {
    SCOPED_TRACE("Cubic Bernstein conversion");
    const axom::Array<double> bernstein_coeffs {2., 3., 5., 8.};
    const axom::Array<double> expected_monomial {2., 3., 3., 0.};
    const auto monomial_coeffs = axom::numerics::bernstein_to_monomial(bernstein_coeffs.view());

    ASSERT_EQ(expected_monomial.size(), monomial_coeffs.size());
    for(axom::IndexType i = 0; i < expected_monomial.size(); ++i)
    {
      EXPECT_DOUBLE_EQ(expected_monomial[i], monomial_coeffs[i]);
    }
  }
}

TEST(numerics_polynomial_solvers, evaluate_polynomial)
{
  const axom::Array<double> coeffs_descending {1., -3., 2.};
  const std::complex<double> x {1., 1.};
  const std::complex<double> value = axom::numerics::evaluate_polynomial(coeffs_descending.view(), x);

  EXPECT_DOUBLE_EQ(-1., value.real());
  EXPECT_DOUBLE_EQ(-1., value.imag());
}

TEST(numerics_polynomial_solvers, solve_polynomial_durand_kerner)
{
  auto expect_complex_eq = [](const axom::Array<std::complex<double>>& expected,
                              const axom::Array<std::complex<double>>& actual,
                              double tol = 1e-10) {
    ASSERT_EQ(expected.size(), actual.size());
    for(axom::IndexType i = 0; i < expected.size(); ++i)
    {
      EXPECT_NEAR(expected[i].real(), actual[i].real(), tol);
      EXPECT_NEAR(expected[i].imag(), actual[i].imag(), tol);
    }
  };

  const axom::Array<std::complex<double>> seeds {std::complex<double> {0.4, 0.9},
                                                 std::complex<double> {0.8, 0.6},
                                                 std::complex<double> {0.99, 0.1}};

  auto seed_trace = [](std::complex<double> seed) {
    return ::testing::Message() << "Seed=(" << seed.real() << ", " << seed.imag() << ")";
  };

  constexpr int max_iters = 200;
  constexpr double tol = 1e-12;

  {
    SCOPED_TRACE("Two distinct real roots");
    const axom::Array<double> coeffs_ascending {6., -5., 1.};
    axom::Array<std::complex<double>> expected_roots {std::complex<double> {2., 0.},
                                                      std::complex<double> {3., 0.}};
    for(const auto& seed : seeds)
    {
      SCOPED_TRACE(seed_trace(seed));
      const auto result =
        axom::numerics::solve_polynomial_durand_kerner_checked(coeffs_ascending.view(),
                                                               tol,
                                                               max_iters,
                                                               seed);
      EXPECT_TRUE(result.converged);
      EXPECT_EQ(result.effective_degree, 2);
      EXPECT_LE(result.max_residual, 1e-10);
      expect_complex_eq(expected_roots, result.roots);
    }
  }

  {
    SCOPED_TRACE("Purely imaginary roots");
    const axom::Array<double> coeffs_ascending {1., 0., 1.};
    axom::Array<std::complex<double>> expected_roots {std::complex<double> {0., -1.},
                                                      std::complex<double> {0., 1.}};
    for(const auto& seed : seeds)
    {
      SCOPED_TRACE(seed_trace(seed));
      const auto result =
        axom::numerics::solve_polynomial_durand_kerner_checked(coeffs_ascending.view(),
                                                               tol,
                                                               max_iters,
                                                               seed);
      EXPECT_TRUE(result.converged);
      EXPECT_EQ(result.effective_degree, 2);
      EXPECT_LE(result.max_residual, 1e-10);
      expect_complex_eq(expected_roots, result.roots);
    }
  }

  {
    SCOPED_TRACE("Trailing zero high-order coefficient is ignored");
    const axom::Array<double> coeffs_ascending {2., -3., 1., 0.};
    axom::Array<std::complex<double>> expected_roots {std::complex<double> {1., 0.},
                                                      std::complex<double> {2., 0.}};
    for(const auto& seed : seeds)
    {
      SCOPED_TRACE(seed_trace(seed));
      const auto roots =
        axom::numerics::solve_polynomial_durand_kerner(coeffs_ascending.view(), tol, seed);
      expect_complex_eq(expected_roots, roots);
    }
  }

  {
    // Durand-Kerner converges slowly near repeated roots due to the denominator
    // becoming small (product of differences approaches zero). We use 1e-4 tolerance
    // here because the roots cluster rather than separate to machine precision.
    SCOPED_TRACE("Repeated-root cubic is retained when residual is small");
    const axom::Array<double> coeffs_ascending {-1., 3., -3., 1.};
    for(const auto& seed : seeds)
    {
      SCOPED_TRACE(seed_trace(seed));
      const auto result =
        axom::numerics::solve_polynomial_durand_kerner_checked(coeffs_ascending.view(),
                                                               tol,
                                                               max_iters,
                                                               seed);
      const auto roots =
        axom::numerics::solve_polynomial_durand_kerner(coeffs_ascending.view(), tol, seed);

      EXPECT_TRUE(result.converged);
      EXPECT_EQ(result.effective_degree, 3);
      EXPECT_LE(result.max_residual, 1e-10);
      ASSERT_EQ(result.roots.size(), 3);
      ASSERT_EQ(roots.size(), 3);
      for(axom::IndexType i = 0; i < roots.size(); ++i)
      {
        EXPECT_NEAR(roots[i].real(), 1., 1e-4);
        EXPECT_NEAR(roots[i].imag(), 0., 1e-4);
      }
    }
  }

  {
    // Despite tiny coefficients (1e-20 scale), the polynomial (t-1)^2 has the same
    // root structure. The repeated root at t=1 still exhibits clustering behavior,
    // so we use 1e-4 tolerance to account for the slow Durand-Kerner convergence.
    SCOPED_TRACE("Low-scale polynomial keeps its effective degree");
    const axom::Array<double> coeffs_ascending {1e-20, -2e-20, 1e-20};
    for(const auto& seed : seeds)
    {
      SCOPED_TRACE(seed_trace(seed));
      const auto result =
        axom::numerics::solve_polynomial_durand_kerner_checked(coeffs_ascending.view(),
                                                               tol,
                                                               max_iters,
                                                               seed);
      const auto roots =
        axom::numerics::solve_polynomial_durand_kerner(coeffs_ascending.view(), tol, seed);

      EXPECT_TRUE(result.converged);
      EXPECT_EQ(result.effective_degree, 2);
      EXPECT_LE(result.max_residual, 1e-10);
      ASSERT_EQ(result.roots.size(), 2);
      ASSERT_EQ(roots.size(), 2);
      for(axom::IndexType i = 0; i < roots.size(); ++i)
      {
        EXPECT_NEAR(roots[i].real(), 1., 1e-4);
        EXPECT_NEAR(roots[i].imag(), 0., 1e-4);
      }
    }
  }
}
