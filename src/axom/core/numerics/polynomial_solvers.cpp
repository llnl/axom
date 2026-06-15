// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/core/numerics/polynomial_solvers.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>

namespace axom
{
namespace numerics
{
namespace
{
using Complex = std::complex<double>;

double coefficient_scale(ArrayView<const double> coeffs)
{
  double scale = 0.0;
  for(double coeff : coeffs)
  {
    scale = axom::utilities::max(scale, std::abs(coeff));
  }
  return scale;
}
}  // namespace

//------------------------------------------------------------------------------
axom::Array<double> bernstein_to_monomial(ArrayView<const double> bernstein_coeffs)
{
  axom::Array<double> monomial_coeffs(bernstein_coeffs);

  const int degree = static_cast<int>(bernstein_coeffs.size()) - 1;

  for(int order = degree; order >= 1; --order)
  {
    for(int j = order; j <= degree; ++j)
    {
      monomial_coeffs[j] -= monomial_coeffs[j - 1];
    }
  }

  for(int k = 0; k <= degree; ++k)
  {
    monomial_coeffs[k] *= axom::utilities::binomialCoefficient(degree, k);
  }

  return monomial_coeffs;
}

//------------------------------------------------------------------------------
int effective_polynomial_degree(ArrayView<const double> coeffs_ascending, double tol)
{
  const double trim_threshold = tol * coefficient_scale(coeffs_ascending);
  int degree = static_cast<int>(coeffs_ascending.size()) - 1;
  while(degree > 0 && std::abs(coeffs_ascending[degree]) <= trim_threshold)
  {
    --degree;
  }
  return degree;
}

//------------------------------------------------------------------------------
PolynomialRootResult solve_polynomial_durand_kerner_checked(ArrayView<const double> coeffs_ascending,
                                                            double tol,
                                                            int max_iters,
                                                            Complex seed)
{
  PolynomialRootResult result;
  const int degree = static_cast<int>(coeffs_ascending.size()) - 1;

  // Degree 0 or negative indicates a constant polynomial with no finite roots.
  // Empty root array with converged=true signals successful handling of the degenerate case.
  if(degree <= 0)
  {
    result.converged = true;
    return result;
  }

  const int effective_degree = effective_polynomial_degree(coeffs_ascending, tol);
  result.effective_degree = effective_degree;

  // Effective degree 0 means the polynomial has been trimmed to a near-constant
  // after removing near-zero high-order terms. Same handling as degree 0 above.
  if(effective_degree <= 0)
  {
    result.converged = true;
    return result;
  }

  // Normalize so the Durand-Kerner iteration works with a monic polynomial (highest coefficient 1)
  // in descending power order, regardless of the input scale.
  axom::Array<double> coeffs_descending(effective_degree + 1, effective_degree + 1);
  for(int i = 0; i <= effective_degree; ++i)
  {
    coeffs_descending[i] =
      coeffs_ascending[effective_degree - i] / coeffs_ascending[effective_degree];
  }

  result.roots.resize(effective_degree);
  // Powers of a non-real base with magnitude around 1 give deterministic, distinct initial guesses.
  // The seed should be a complex number on or near the unit circle so that the initial guesses
  // are neither too large nor too small.
  for(int i = 0; i < effective_degree; ++i)
  {
    result.roots[i] = std::pow(seed, i + 1);
  }

  // Each sweep updates every root estimate by subtracting p(z_i) divided by
  // the product of differences from the current guesses for all other roots.
  for(int iter = 0; iter < max_iters; ++iter)
  {
    double max_update = 0.0;
    for(int i = 0; i < effective_degree; ++i)
    {
      Complex denom {1.0, 0.0};
      for(int j = 0; j < effective_degree; ++j)
      {
        if(i != j)
        {
          denom *= result.roots[i] - result.roots[j];
        }
      }

      // Repeated or tightly clustered roots can make the product nearly zero,
      // so guard the update against division by an unstable denominator.
      if(std::abs(denom) <= tol)
      {
        denom = Complex {tol, tol};
      }

      const Complex update = evaluate_polynomial(coeffs_descending.view(), result.roots[i]) / denom;
      result.roots[i] -= update;
      max_update = axom::utilities::max(max_update, std::abs(update));
    }

    result.iterations = iter + 1;
    result.max_update = max_update;
    if(max_update <= tol)
    {
      result.converged = true;
      break;
    }
  }

  for(const auto& root : result.roots)
  {
    result.max_residual =
      axom::utilities::max(result.max_residual,
                           std::abs(evaluate_polynomial(coeffs_descending.view(), root)));
  }

  // Residual-based acceptance lets repeated-root cases succeed even when the update size stalls
  // above tol near a multiple root. The 100x multiplier accounts for accumulated roundoff error.
  constexpr double residual_acceptance_factor = 100.0;
  result.converged = result.converged || result.max_residual <= residual_acceptance_factor * tol;

  std::sort(result.roots.begin(), result.roots.end(), [](const Complex& lhs, const Complex& rhs) {
    if(lhs.real() != rhs.real())
    {
      return lhs.real() < rhs.real();
    }
    return lhs.imag() < rhs.imag();
  });

  return result;
}

//------------------------------------------------------------------------------
axom::Array<Complex> solve_polynomial_durand_kerner(ArrayView<const double> coeffs_ascending,
                                                    double tol,
                                                    Complex seed)
{
  const auto result = solve_polynomial_durand_kerner_checked(coeffs_ascending, tol, 200, seed);
  return result.converged ? result.roots : axom::Array<Complex> {};
}

//------------------------------------------------------------------------------
int solve_linear(ArrayView<const double> coeff, ArrayView<double> roots, int& numRoots)
{
  assert(coeff.size() >= 2);
  assert(roots.size() >= 1);

  int status = -1;

  // solve ax + b = 0
  const double a = coeff[1];
  const double b = coeff[0];

  if(utilities::isNearlyEqual(a, 0.))
  {
    if(utilities::isNearlyEqual(b, 0.))
    {
      // Infinite solutions: a horizontal line on the X-axis.
      status = 0;
      numRoots = -1;
    }
    else
    {
      // No solutions: a horizontal line not on the X-axis.
      numRoots = 0;
    }
  }
  else
  {
    // One solution, where the line crosses the X-axis.
    status = 0;
    numRoots = 1;
    roots[0] = -b / a;
  }

  return status;
}

//------------------------------------------------------------------------------
#ifdef _MSC_VER
  #pragma warning(push)
  #pragma warning(disable : 4244)
#endif
int solve_linear(const double* coeff, double* roots, int& numRoots)
{
  return solve_linear(ArrayView<const double>(coeff,
                                              axom::StackArray<axom::IndexType, 1> {2},
                                              axom::StackArray<axom::IndexType, 1> {1}),
                      ArrayView<double>(roots,
                                        axom::StackArray<axom::IndexType, 1> {1},
                                        axom::StackArray<axom::IndexType, 1> {1}),
                      numRoots);
}

//------------------------------------------------------------------------------
int solve_quadratic(ArrayView<const double> coeff, ArrayView<double> roots, int& numRoots)
{
  assert(coeff.size() >= 3);
  assert(roots.size() >= 2);

  int status = -1;

  // solve ax^2 + bx + c = 0
  const double a = coeff[2];
  const double b = coeff[1];
  const double c = coeff[0];

  if(utilities::isNearlyEqual(a, 0.))
  {
    // If this system is nearly linear, solve it as such.
    return solve_linear(coeff, roots, numRoots);
  }

  const double discriminant = b * b - 4 * a * c;
  const double overtwoa = 1. / (2 * a);

  if(utilities::isNearlyEqual(discriminant, 0.))
  {
    // One unique real root
    status = 0;
    numRoots = 1;
    roots[0] = roots[1] = -b * overtwoa;
  }
  else if(discriminant < 0)
  {
    // No real roots
    numRoots = 0;
  }
  else
  {
    // Two real roots
    status = 0;
    numRoots = 2;
    const double sqrtdisc = std::sqrt(discriminant);
    roots[0] = (-b + sqrtdisc) * overtwoa;
    roots[1] = (-b - sqrtdisc) * overtwoa;
  }

  return status;
}

//------------------------------------------------------------------------------
int solve_quadratic(const double* coeff, double* roots, int& numRoots)
{
  return solve_quadratic(ArrayView<const double>(coeff,
                                                 axom::StackArray<axom::IndexType, 1> {3},
                                                 axom::StackArray<axom::IndexType, 1> {1}),
                         ArrayView<double>(roots,
                                           axom::StackArray<axom::IndexType, 1> {2},
                                           axom::StackArray<axom::IndexType, 1> {1}),
                         numRoots);
}

inline double cuberoot(double x)
{
  // pow(x, y) returns NaN for negative finite x and noninteger y.
  if(x < 0)
  {
    return -pow(-x, 1. / 3.);
  }
  else
  {
    return pow(x, 1. / 3.);
  }
}

//------------------------------------------------------------------------------
int solve_cubic(ArrayView<const double> coeff, ArrayView<double> roots, int& numRoots)
{
  assert(coeff.size() >= 4);
  assert(roots.size() >= 3);

  int status = -1;

  // Here I use variable names as presented in Korn & Korn:
  // x^3 + ax^2 + bx + c = 0
  double cubecoeff = coeff[3];
  double a = coeff[2];
  double b = coeff[1];
  double c = coeff[0];

  if(utilities::isNearlyEqual(cubecoeff, 0.))
  {
    // If this system is nearly quadratic, solve it as such.
    return solve_quadratic(coeff, roots, numRoots);
  }

  // We normalize by dividing all by the cubic coefficient.
  const double invcubecoeff = 1. / cubecoeff;
  a *= invcubecoeff;
  b *= invcubecoeff;
  c *= invcubecoeff;

  // Note that p and q differ by a multiplicative constant from Korn,
  // because they're always used with that multiplication.
  const double p = (-a * a + 3 * b) / 9;                      //  1/3 Korn's p
  const double q = (a * (-2 * a * a + 9 * b) - 27 * c) / 54;  // -1/2 Korn's q

  const double Q = p * p * p + q * q;  // actual discriminant == -108Q
  // the term chVar occurs because we've changed variables
  // (x = y - a/3) and we need to change back to x.
  const double onethird = 1. / 3.;
  const double chVar = -a * onethird;

  if(utilities::isNearlyEqual(Q, 0.))
  {
    // We have three real roots, and at least two are equal.
    if(utilities::isNearlyEqual(q, 0.))
    {
      numRoots = 1;
    }
    else
    {
      numRoots = 2;
    }
    status = 0;

    const double cuberootq = cuberoot(q);
    roots[0] = chVar + 2 * cuberootq;
    roots[1] = chVar - cuberootq;
    roots[2] = roots[1];
  }
  else if(Q > 0)
  {
    // We have one real root, and two complex roots.
    // Right now we're calculating the real root only, but the complex roots
    // can be easily added by un-commenting the calculation.
    numRoots = 1;
    status = 0;

    const double sqrtQ = sqrt(Q);
    const double A = cuberoot(q + sqrtQ);
    const double B = cuberoot(q - sqrtQ);

    roots[0] = chVar + A + B;
    roots[1] = 0;
    roots[2] = 0;
    // double imagpart = sqrt(3) * (A - B)/2;
    // double stavg = (A + B)/2;
    // Here note use of imaginary i:
    // roots[1] = chVar - stavg - i*imagpart;
    // roots[2] = chVar - stavg + i*imagpart;
  }
  else
  {
    // Q < 0, and we have three distinct real roots.
    // The case for Q == 0 is a special case of the case for Q > 0, and
    // both correspond to Cardano's method proper, reported in section
    // 1.8-3 of Korn, where the quadratic term ax^2 is eliminated by a
    // change of variable x = y - a/3.  For "irreducible" cubics (this
    // case, where the variable substitution didn't eliminate the quadratic),
    // we use the trigonometric solution to the cubic equation, reported
    // in section 1.8-4a of Korn.
    numRoots = 3;
    status = 0;

    const double alpha = acos(q / sqrt(-p * p * p));
    const double m = 2 * sqrt(-p);

    roots[0] = chVar + m * cos(alpha * onethird);
    roots[1] = chVar - m * cos((alpha + M_PI) * onethird);
    roots[2] = chVar - m * cos((alpha - M_PI) * onethird);
  }

  return status;
}

//------------------------------------------------------------------------------
int solve_cubic(const double* coeff, double* roots, int& numRoots)
{
  return solve_cubic(ArrayView<const double>(coeff,
                                             axom::StackArray<axom::IndexType, 1> {4},
                                             axom::StackArray<axom::IndexType, 1> {1}),
                     ArrayView<double>(roots,
                                       axom::StackArray<axom::IndexType, 1> {3},
                                       axom::StackArray<axom::IndexType, 1> {1}),
                     numRoots);
}

#ifdef _MSC_VER
  #pragma warning(pop)
#endif

} /* end namespace numerics */
} /* end namespace axom */
