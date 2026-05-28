// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/detail/clipping/TetrahedronClipUtils.hpp"

#include "gtest/gtest.h"

#include <algorithm>
#include <array>

namespace
{

using axom::quest::experimental::detail::clipTetByVertexValues;

TEST(quest_tetrahedron_clip_utils, affine_volume_cases)
{
  constexpr double volume = 2.5;

  const double outside[4] = {-1.0, -2.0, -3.0, -4.0};
  EXPECT_DOUBLE_EQ(0.0, clipTetByVertexValues(outside, volume));

  const double inside[4] = {1.0, 2.0, 3.0, 4.0};
  EXPECT_DOUBLE_EQ(volume, clipTetByVertexValues(inside, volume));

  const double oneInside[4] = {1.0, -1.0, -1.0, -1.0};
  EXPECT_DOUBLE_EQ(volume / 8.0, clipTetByVertexValues(oneInside, volume));

  const double threeInside[4] = {-1.0, 1.0, 1.0, 1.0};
  EXPECT_DOUBLE_EQ(7.0 * volume / 8.0, clipTetByVertexValues(threeInside, volume));

  const double twoInside[4] = {1.0, 1.0, -1.0, -1.0};
  EXPECT_DOUBLE_EQ(volume / 2.0, clipTetByVertexValues(twoInside, volume));
}

TEST(quest_tetrahedron_clip_utils, two_inside_is_permutation_invariant)
{
  constexpr double volume = 3.0;
  constexpr double expectedFraction = 37.0 / 105.0;
  std::array<double, 4> values {{-4.0, -2.0, 1.0, 3.0}};
  do
  {
    EXPECT_NEAR(volume * expectedFraction, clipTetByVertexValues(values.data(), volume), 1e-14);
  } while(std::next_permutation(values.begin(), values.end()));
}

}  // namespace
