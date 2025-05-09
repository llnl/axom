// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_NUMERICS_POLY_SOLVE_HPP_
#define AXOM_NUMERICS_POLY_SOLVE_HPP_

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
 * Note that coeff[i] = \f$ a_i \f$.  The constant term goes in coeff[0],
 * the linear in coeff[1], quadratic in coeff[2], and so forth.
 */

namespace axom
{
namespace numerics
{
/*!
 * \brief Find the real root for a linear equation of form \f$ ax + b = 0 \f$.
 *
 * \param [in] coeff Equation coefficients: coeff[i] multiplies \f$ x^i \f$.
 *   So, coeff[0] = b and coeff[1] = a.
 * \param [out] roots The real roots of the equation.
 * \param [out] numRoots The number of distinct, real roots found (max: 1).
 *   If the line lies on the X-axis, numRoots is assigned -1 to indicate
 *   infinitely many solutions.  If the line doesn't intersect the X-axis,
 *   numRoots is assigned 0.
 * \return 0 for success, or -1 to indicate an inconsistent equation
 *   (all coefficients 0 except for coeff[0]).
 *
 * \pre coeff points to an array of length at least 2.
 * \pre roots points to an array of length at least 1.
 */
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
 *   assigned -1 to indicate infinitely many solutions.  If the equation
 *   doesn't intersect the X-axis, numRoots is assigned 0.
 * \return 0 for success, or -1 to indicate an inconsistent equation
 *   (all coefficients 0 except for coeff[0]).
 *
 * \pre coeff points to an array of length at least 3.
 * \pre roots points to an array of length at least 2.
 */
int solve_quadratic(const double* coeff, double* roots, int& numRoots);

/*!
 * \brief For a cubic equation of the form \f$ ax^3 + bx^2 + cx + d = 0 \f$,
 * find real roots.
 *
 * A closed-form solution for cubic equations was published in Cardano's
 * *Ars Magna* of 1545.  This can be summarized as follows:
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
 * three distinct real roots.  If there is one distinct real root, the root
 * can either have multiplicity 3, or have multiplicity 1 with two complex
 * roots.  In the case of two real roots, one root will have multiplicity 2.
 * If the discriminant \f$ d = -27t^2 - 4(-3^3) \f$ is zero, the equation has
 * at least one root with multiplicity greater than 1.  If \f$ d > 0, \f$
 * there are three distinct real roots, and if \f$ d > 0, \f$ there is one
 * real root.
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
 *   assigned -1 to indicate infinitely many solutions.  If the equation
 *   doesn't intersect the X-axis, numRoots is assigned 0.
 * \return 0 for success, or -1 to indicate an inconsistent equation
 *   (all coefficients 0 except for coeff[0]).
 *
 * \pre coeff points to an array of length at least 4.
 * \pre roots points to an array of length at least 3.
 */
int solve_cubic(const double* coeff, double* roots, int& numRoots);

} /* end namespace numerics */
} /* end namespace axom */

#endif  // AXOM_NUMERICS_POLY_SOLVE_HPP_
