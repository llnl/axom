// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core.hpp"

#include "axom/primal/geometry/Octahedron.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"

#include "axom/primal/operators/split.hpp"

#include <cmath>

using OctType = axom::primal::Octahedron<double, 3>;
using TetType = axom::primal::Tetrahedron<double, 3>;
using PointType = OctType::PointType;
using ArrayType = axom::Array<TetType>;

TEST(primal_split, split_basics)
{
  // Here is a triangular prism stored as an octahedron.  It has two parallel end-caps:
  //
  //   - PRT, at z=0: (0, 0, 0), (1, 0, 0), (0, 1, 0)
  //   - UQS, at z=0: (0, 0, 1), (1, 0, 1), (0, 1, 1)
  //
  // This prism-stored-as-octahedron and the arrangement of vertices is typical of
  // the discretization algorithm in quest.  It also is easy to compute the volume, 0.5.

  // clang-format off
  PointType P = PointType { 0., 0., 0. };
  PointType Q = PointType { 1., 0., 1. };
  PointType R = PointType { 1., 0., 0. };
  PointType S = PointType { 0., 1., 1. };
  PointType T = PointType { 0., 1., 0. };
  PointType U = PointType { 0., 0., 1. };
  // clang-format on

  OctType oct(P, Q, R, S, T, U);

  ArrayType tets;
  int old_tetlength = tets.size();

  axom::primal::split(oct, tets);

  constexpr int tets_in_oct = 8;
  EXPECT_EQ(tets.size() - old_tetlength, tets_in_oct);

  double octvol = 0.;

  for(int i = 0; i < tets_in_oct; ++i)
  {
    EXPECT_FALSE(tets[old_tetlength + i].degenerate());
    octvol += tets[old_tetlength + i].signedVolume();
  }

  double sign_adjustment = -1.;
  octvol *= sign_adjustment;

  constexpr double expected_octvol = 0.5;
  EXPECT_DOUBLE_EQ(octvol, expected_octvol);
}
