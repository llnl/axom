// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/Utilities.hpp"
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
 * \param [out] nodes The array of 1D nodes
 * \param [out] weights The array of weights
 * 
 * A Gauss-Legendre rule with \a npts points can exactly integrate
 *  polynomials of order 2 * npts - 1
 *
 * Algorithm adapted from the MFEM implementation in `mfem/fem/intrules.cpp`
 * 
 * \note This method constructs the points by scratch each time, without caching
 * \sa get_gauss_legendre(int)
 */
void compute_gauss_legendre_data(int npts,
                                 axom::Array<double>& nodes,
                                 axom::Array<double>& weights,
                                 int allocatorID)
{
  assert("Quadrature rules must have >= 1 point" && (npts >= 1));

  nodes = axom::Array<double>(npts, npts, allocatorID);
  weights = axom::Array<double>(npts, npts, allocatorID);

  if(npts == 1)
  {
    nodes[0] = 0.5;
    weights[0] = 1.0;
    return;
  }
  if(npts == 2)
  {
    nodes[0] = 0.21132486540518711775;
    nodes[1] = 0.78867513459481288225;

    weights[0] = weights[1] = 0.5;
    return;
  }
  if(npts == 3)
  {
    nodes[0] = 0.11270166537925831148207345;
    nodes[1] = 0.5;
    nodes[2] = 0.88729833462074168851792655;

    weights[0] = 0.2777777777777777777777778;
    weights[1] = 0.4444444444444444444444444;
    weights[2] = 0.2777777777777777777777778;
    return;
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

    nodes[i - 1] = xi;
    nodes[n - i] = 1.0 - xi;

    // For z \in [-1, 1], w_i = 2 / (1 - z^2) / P'_n(z)^2.
    // For nodes[i] = xi = (1 - z)/2 \in [0, 1], weights[i] = 0.5 * w_i
    weights[i - 1] = weights[n - i] = 1.0 / (4.0 * xi * (1.0 - xi) * Pp_n * Pp_n);
  }
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
 * \return The `QuadratureRule` object which contains axom::ArrayView<double>'s of stored nodes and weights
 */
QuadratureRule get_gauss_legendre(int npts, int allocatorID)
{
  assert("Quadrature rules must have >= 1 point" && (npts >= 1));

  // Define a static map that stores the GL quadrature rule for a given order
  static std::map<std::pair<int, int>, std::pair<axom::Array<double>, axom::Array<double>>> rule_library;

  const std::pair<int, int> key = std::make_pair(npts, allocatorID);

  auto value_it = rule_library.find(key);
  if(value_it == rule_library.end())
  {
    auto& vals = rule_library[key];
    compute_gauss_legendre_data(npts, vals.first, vals.second, allocatorID);
    value_it = rule_library.find(key);
  }

  return QuadratureRule {value_it->second.first.view(), value_it->second.second.view()};
}

} /* end namespace numerics */
} /* end namespace axom */
