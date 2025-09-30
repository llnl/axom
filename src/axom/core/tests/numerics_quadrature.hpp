// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/numerics/gauss_legendre.hpp"

TEST(numerics_quadrature, generate_rules)
{
  const int N = 3;

  double nodes[N];
  double weights[N];

  axom::numerics::compute_rule(N, nodes, weights);

  for(int i = 0; i < N; ++i)
  {
    std::cout << nodes[i] << " " << weights[i] << std::endl;
  }
}
