// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/Array.hpp"
#include "axom/core/FlatMap.hpp"
#include "axom/core/NumericLimits.hpp"
#include "axom/core/numerics/quadrature.hpp"

// For math constants and includes
#include "axom/config.hpp"
#include <cmath>

namespace axom
{
namespace numerics
{

/*!
 * \brief Computes a 1D quadrature rule of Gauss-Legendre points 
 *
 * \param [in] npts The number of points in the rule
 * 
 * A Gauss-Legendre rule with \a npts points can exactly integrate
 *  polynomials of order 2 * npts - 1
 *
 * Algorithm adapted from the MFEM implementation in `mfem/fem/intrules.cpp`
 * 
 * \note This method constructs the points by scratch each time, without caching
 * \sa get_gauss_legendre(int)
 *
 * \return The `QuadratureRule` object which contains axom::Array<double>'s of nodes and weights
 */
QuadratureRule compute_gauss_legendre(int npts)
{
  SLIC_ASSERT(npts >= 1, "Quadrature rules must have >= 1 point");

  QuadratureRule rule(npts);

  if(npts == 1)
  {
    rule.m_nodes[0] = 0.5;
    rule.m_weights[0] = 1.0;
    return rule;
  }
  if(npts == 2)
  {
    rule.m_nodes[0] = 0.21132486540518711775;
    rule.m_nodes[1] = 0.78867513459481288225;

    rule.m_weights[0] = rule.m_weights[1] = 0.5;
    return rule;
  }
  if(npts == 3)
  {
    rule.m_nodes[0] = 0.11270166537925831148207345;
    rule.m_nodes[1] = 0.5;
    rule.m_nodes[2] = 0.88729833462074168851792655;

    rule.m_weights[0] = 0.2777777777777777777777778;
    rule.m_weights[1] = 0.4444444444444444444444444;
    rule.m_weights[2] = 0.2777777777777777777777778;
    return rule;
  }

  const int n = npts;
  const int m = (npts + 1) / 2;

  // Nodes are mirrored across x = 0.5, so only need to evaluate half
  for(int i = 1; i <= m; ++i)
  {
    // Each node is the root of a Legendre polynomial,
    //  which are approximately uniformly distributed in arccos(xi).
    // This makes cos a good initial guess for subsequent Newton iterations
    double z = std::cos(M_PI * (i - 0.25) / (n + 0.5));
    double Pp_n, P_n, dz, xi = 0.0;

    bool done = false;
    while(true)
    {
      // Recursively evaluate P_n(z) through the recurrence relation
      //  n * P_n(z) = (2n - 1) * P_{n-1}(z) - (n - 1) * P_{n - 2}(z)
      double P_nm1 = 1.0;  // P_0(z) = 1
      P_n = z;             // P_1(z) = z
      for(int j = 2; j <= n; ++j)
      {
        double P_nm2 = P_nm1;
        P_nm1 = P_n;
        P_n = ((2 * j - 1) * z * P_nm1 - (j - 1) * P_nm2) / j;
      }

      // Evaluate P'_n(z) using another recurrence relation
      //  (z^2 - 1) * P'_n(z) = n * z * P_n(z) - n * P_{n-1}(z)
      Pp_n = n * (z * P_n - P_nm1) / (z * z - 1);

      if(done)
      {
        break;
      }

      // Compute the Newton method step size
      dz = P_n / Pp_n;

      if(std::fabs(dz) < axom::numeric_limits<double>::epsilon())
      {
        done = true;

        // Scale back to [0, 1]
        xi = ((1 - z) + dz) / 2;
      }

      // Take the Newton step, repeat the process
      z -= dz;
    }

    rule.m_nodes[i - 1] = xi;
    rule.m_nodes[n - i] = 1.0 - xi;

    // Evaluate the weight as w_i = 2 / (1 - z^2) / P'_n(z)^2, with z \in [-1, 1]
    rule.m_weights[i - 1] = rule.m_weights[n - i] = 1.0 / (4.0 * xi * (1.0 - xi) * Pp_n * Pp_n);
  }

  return rule;
}

/*!
 * \brief Computes or accesses a precomputed 1D quadrature rule of Gauss-Legendre points 
 *
 * \param [in] npts The number of points in the rule
 * 
 * A Gauss-Legendre rule with \a npts points can exactly integrate
 *  polynomials of order 2 * npts - 1
 *
 * \note If this method has already been called for a given order, it will reuse the same quadrature points
 *  without needing to recompute them
 *
 * \warning The use of a static variable to store cached nodes makes this method not threadsafe.
 * 
 * \return The `QuadratureRule` object which contains axom::Array<double>'s of nodes and weights
 */
const QuadratureRule& get_gauss_legendre(int npts)
{
  SLIC_ASSERT(npts >= 1, "Quadrature rules must have >= 1 point");

  // Define a static map that stores the GL quadrature rule for a given order
  static axom::FlatMap<int, QuadratureRule> rule_library;
  if(rule_library.find(npts) == rule_library.end())
  {
    rule_library[npts] = compute_gauss_legendre(npts);
  }

  return rule_library[npts];
}

} /* end namespace numerics */
} /* end namespace axom */
