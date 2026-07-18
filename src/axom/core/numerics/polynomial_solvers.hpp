// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"

#include <complex>
#include <type_traits>

/*!
 * \file polynomial_solve.hpp
 * The functions declared in this header file find real roots of polynomials
 * of the form
 *
 *   \f[ \sum a_i * x^i = 0. \f]
 *
 * The functions take an input array of values for coefficients, an
 * output array for roots found, and an output int indicating the number
 * of distinct real roots found.
 *
 * Note that coeff[i] = \f$ a_i \f$. The constant term goes in coeff[0],
 * the linear in coeff[1], quadratic in coeff[2], and so forth.
 */

namespace axom
{
namespace numerics
{

/*!
 * \brief Evaluate a polynomial with coefficients stored in descending power order.
 *
 * \tparam ScalarType The scalar type used for evaluation, e.g. `double` or `std::complex<double>`.
 * \param [in] coeffs_descending Polynomial coefficients in descending power order.
 * \param [in] x The evaluation point.
 *
 * \return The polynomial value at `x`.
 */
template <typename ScalarType>
ScalarType evaluate_polynomial(ArrayView<const double> coeffs_descending, const ScalarType& x)
{
  static_assert(
    std::is_same_v<ScalarType, double> || std::is_same_v<ScalarType, std::complex<double>>,
    "evaluate_polynomial requires ScalarType to be double or std::complex<double>.");

  ScalarType value {0.0};
  for(double coeff : coeffs_descending)
  {
    value = value * x + coeff;
  }
  return value;
}

/*!
 * \brief Convert Bernstein coefficients to monomial coefficients on `[0, 1]`.
 *
 * Given Bernstein coefficients \f$ b_i \f$ of degree \f$ n \f$, this returns
 * monomial coefficients \f$ a_i \f$ such that
 * \f$ \sum_i b_i B_i^n(t) = \sum_i a_i t^i \f$.
 *
 * \param [in] bernstein_coeffs Coefficients in the Bernstein basis.
 *
 * \return Coefficients in ascending monomial order.
 */
axom::Array<double> bernstein_to_monomial(ArrayView<const double> bernstein_coeffs);

/*!
 * \brief Return the effective degree of a polynomial after trimming near-zero
 * leading coefficients.
 *
 * \param [in] coeffs_ascending Polynomial coefficients in ascending power order.
 * \param [in] tol Relative trimming tolerance. Coefficients with magnitude at
 * or below `tol * max(abs(coeffs_ascending))` are treated as zero when
 * trimming the highest-degree terms.
 *
 * \return The largest exponent whose coefficient exceeds the scaled trimming
 * threshold, or `0` if the polynomial is constant to within the requested
 * tolerance.
 */
int effective_polynomial_degree(ArrayView<const double> coeffs_ascending, double tol = 1e-12);

/// \brief Result metadata for an all-roots polynomial solve.
struct PolynomialRootResult
{
  axom::Array<std::complex<double>> roots;
  bool converged {false};
  int iterations {0};
  int effective_degree {0};
  double max_update {0.0};
  double max_residual {0.0};
};

/*!
 * \brief Approximate all roots of a polynomial using the Durand-Kerner method.
 *
 * \param [in] coeffs_ascending Polynomial coefficients in ascending power order.
 * \param [in] tol Iteration tolerance and relative effective-zero threshold
 * used when trimming leading coefficients.
 * \param [in] seed Seed used to generate deterministic, distinct initial guesses.
 * This should be a complex number on or near the unit circle (and not purely real),
 * so that `seed^(i+1)` spreads the starting points around a reasonable radius.
 *
 * \return The polynomial roots, sorted by real part and then imaginary part,
 * or an empty array if the iteration does not converge.
 */
axom::Array<std::complex<double>> solve_polynomial_durand_kerner(
  ArrayView<const double> coeffs_ascending,
  double tol = 1e-12,
  std::complex<double> seed = std::complex<double> {0.4, 0.9});

/*!
 * \brief Approximate all roots of a polynomial using the Durand-Kerner method,
 * including convergence diagnostics.
 *
 * \param [in] coeffs_ascending Polynomial coefficients in ascending power order.
 * \param [in] tol Iteration tolerance and relative effective-zero threshold
 * used when trimming leading coefficients.
 * \param [in] max_iters Maximum number of Durand-Kerner iterations. Default is 200,
 * which provides reliable convergence for well-conditioned polynomials up to degree ~15-20. 
 * Higher-degree or ill-conditioned polynomials (repeated roots, clustered roots) 
 * may require more iterations or specialized methods.
 * \param [in] seed Seed used to generate deterministic, distinct initial guesses.
 * This should be a complex number on or near the unit circle (and not purely real),
 * so that `seed^(i+1)` spreads the starting points around a reasonable radius.
 *
 * \return The roots and convergence metadata for the solve.
 */
PolynomialRootResult solve_polynomial_durand_kerner_checked(
  ArrayView<const double> coeffs_ascending,
  double tol = 1e-12,
  int max_iters = 200,
  std::complex<double> seed = std::complex<double> {0.4, 0.9});

/*!
 * \brief Find the real root for a linear equation of form \f$ ax + b = 0 \f$.
 *
 * \param [in] coeff Equation coefficients: coeff[i] multiplies \f$ x^i \f$.
 *   So, coeff[0] = b and coeff[1] = a.
 * \param [out] roots The real roots of the equation.
 * \param [out] numRoots The number of distinct, real roots found (max: 1).
 *   If the line lies on the X-axis, numRoots is assigned -1 to indicate
 *   infinitely many solutions. If the line does not intersect the X-axis,
 *   numRoots is assigned 0.
 *
 * \return 0 for success, or -1 to indicate an inconsistent equation
 *   (all coefficients 0 except for coeff[0]).
 *
 * \pre coeff.size() >= 2
 * \pre roots.size() >= 1
 */
int solve_linear(ArrayView<const double> coeff, ArrayView<double> roots, int& numRoots);

/*!
 * \brief Deprecated pointer-based overload for \ref solve_linear(ArrayView<const
 * double>, ArrayView<double>, int&).
 */
[[deprecated(
  "Use solve_linear(axom::ArrayView<const double>, axom::ArrayView<double>, int&) instead.")]]
int solve_linear(const double* coeff, double* roots, int& numRoots);

/*!
 * \brief For a quadratic equation of the form \f$ ax^2 + bx + c = 0 \f$,
 * find real roots using the quadratic formula.
 *
 * \param [in] coeff Equation coefficients: coeff[i] multiplies \f$ x^i \f$.
 *   So, coeff[0] = c, coeff[1] = b, coeff[2] = a.
 * \param [out] roots The real roots of the equation.
 * \param [out] numRoots The number of distinct, real roots found (max: 2).
 *   If the equation degenerates to a line lying on the X-axis, numRoots is
 *   assigned -1 to indicate infinitely many solutions. If the equation does
 *   not intersect the X-axis, numRoots is assigned 0.
 *
 * \return 0 for success, or -1 to indicate an inconsistent equation
 *   (all coefficients 0 except for coeff[0]).
 *
 * \pre coeff.size() >= 3
 * \pre roots.size() >= 2
 */
int solve_quadratic(ArrayView<const double> coeff, ArrayView<double> roots, int& numRoots);

/*!
 * \brief Deprecated pointer-based overload for \ref
 * solve_quadratic(ArrayView<const double>, ArrayView<double>, int&).
 */
[[deprecated(
  "Use solve_quadratic(axom::ArrayView<const double>, axom::ArrayView<double>, int&) instead.")]]
int solve_quadratic(const double* coeff, double* roots, int& numRoots);

/*!
 * \brief For a cubic equation of the form \f$ ax^3 + bx^2 + cx + d = 0 \f$,
 * find real roots.
 *
 * A closed-form solution for cubic equations was published in Cardano's
 * *Ars Magna* of 1545. This can be summarized as follows:
 *
 * Start with the cubic equation
 *
 *   \f[ ax^3 + bx^2 + cx + d = 0. \f]
 *
 * Divide through by a to normalize, then define the following:
 *
 *   \f[ p = b^2 - 3c \f]
 *   \f[ q = -\frac{27}{2}d - b^3 + \frac{9}{2}cb \f]
 *   \f[ t = 2p^{-3/2}q \f]
 *   \f[ y = \frac{3}{\sqrt{p}}(x + \frac{1}{3}b) \f]
 *
 * Then the original cubic equation can be written as \f$ y^3 - 3y = t \f$
 * with a solution \f$ y = \frac{1}{u} + u \f$, with
 *
 *   \f[ u = \sqrt[3]{\frac{t}{2} \pm \sqrt{\frac{t^2}{4} - 1}}. \f]
 *
 * Because the cubic is an odd function, there can be either one, two, or
 * three distinct real roots. If there is one distinct real root, the root can
 * either have multiplicity 3, or have multiplicity 1 with two complex roots.
 * In the case of two real roots, one root will have multiplicity 2. If the
 * discriminant \f$ d = -27t^2 - 4(-3^3) \f$ is zero, the equation has at
 * least one root with multiplicity greater than 1. If \f$ d > 0, \f$ there
 * are three distinct real roots, and if \f$ d < 0, \f$ there is one real
 * root.
 *
 * See J. Kopp, Efficient numerical diagonalization of hermitian 3x3 matrices,
 * Int.J.Mod.Phys. C19:523-548, 2008 (https://arxiv.org/abs/physics/0610206)
 * and G. A. Korn and T. M. Korn, "Mathematical Handbook for Scientists and
 * Engineers," QA37 K84 1968 in the library; section 1.8-3 "Cubic Equations",
 * p. 23.
 *
 * \param [in] coeff Equation coefficients: coeff[i] multiplies \f$ x^i \f$.
 *   So, coeff[0] = d, coeff[1] = c, coeff[2] = b, coeff[3] = a.
 * \param [out] roots The real roots of the equation.
 * \param [out] numRoots The number of distinct, real roots found (max: 3).
 *   If the equation degenerates to a line lying on the X-axis, numRoots is
 *   assigned -1 to indicate infinitely many solutions. If the equation does
 *   not intersect the X-axis, numRoots is assigned 0.
 *
 * \return 0 for success, or -1 to indicate an inconsistent equation
 *   (all coefficients 0 except for coeff[0]).
 *
 * \pre coeff.size() >= 4
 * \pre roots.size() >= 3
 */
int solve_cubic(ArrayView<const double> coeff, ArrayView<double> roots, int& numRoots);

/*!
 * \brief Deprecated pointer-based overload for \ref
 * solve_cubic(ArrayView<const double>, ArrayView<double>, int&).
 */
[[deprecated(
  "Use solve_cubic(axom::ArrayView<const double>, axom::ArrayView<double>, int&) instead.")]]
int solve_cubic(const double* coeff, double* roots, int& numRoots);

}  // namespace numerics
}  // namespace axom
