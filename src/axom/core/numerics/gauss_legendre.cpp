// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core/numerics/gauss_legendre.hpp"

#include <cmath>

namespace axom
{
namespace numerics
{

void compute_rule(int np, double* nodes, double* weights)
{
  if(np == 1)
  {
    nodes[0] = 0.5;
    weights[0] = 1.0;
    return;
  }
  if(np == 2)
  {
    nodes[0] = 0.21132486540518711775;
    nodes[1] = 0.78867513459481288225;

    weights[0] = weights[1] = 0.5;
    return;
  }

  const int n = np;
  const int m = (n + 1) / 2;

  // Nodes are mirrored across x = 0.5, so only need to evaluate half
  for(int i = 1; i <= m; ++i)
  {
    double z = std::cos(M_PI * (i - 0.25) / (n + 0.5));
    double pp, p1, dz, xi = 0.0;

    bool done = false;
    while(true)
    {
      double p2 = 1.0;
      p1 = z;
      for(int j = 2; j <= n; ++j)
      {
        double p3 = p2;
        p2 = p1;
        p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;
      }
      // p1 is Legendre polynomial

      pp = n * (z * p1 - p2) / (z * z - 1);

      if(done)
      {
        break;
      }

      dz = p1 / pp;

      if(std::fabs(dz) < 1e-16)
      {
        done = true;
        xi = ((1 - z) + dz) / 2;
      }

      z -= dz;
    }

    nodes[i - 1] = xi;
    nodes[n - i] = 1.0 - xi;
    weights[i - 1] = weights[n - i] = 1.0 / (4.0 * xi * (1.0 - xi) * pp * pp);
  }
}
} /* end namespace numerics */
} /* end namespace axom */