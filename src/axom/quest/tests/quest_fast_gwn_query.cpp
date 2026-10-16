// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/core.hpp"
#include "axom/mint.hpp"
#include "axom/primal.hpp"
#include "axom/quest/GWNMethods.hpp"

#include "gtest/gtest.h"

namespace primal = axom::primal;
namespace quest = axom::quest;

TEST(quest_fast_gwn_query, direct_and_accelerated_segment_evaluation)
{
  using Point2D = primal::Point<double, 2>;
  using Segment = primal::Segment<double, 2>;
  using Query = quest::FastGWNQuery<axom::SEQ_EXEC, 2>;

  axom::Array<Segment> segments;
  segments.emplace_back(Point2D {-1.0, 0.0}, Point2D {0.0, 1.0});
  segments.emplace_back(Point2D {0.0, 1.0}, Point2D {1.0, 0.0});
  const axom::Array<Segment>& const_segments = segments;
  const auto segment_view = const_segments.view();

  const Point2D query {20.0, -15.0};
  const primal::WindingTolerances tolerances;

  double direct_gwn = 0.0;
  for(const auto& segment : segment_view)
  {
    direct_gwn += primal::winding_number(query, segment, tolerances);
  }

  Query fast_query;
  EXPECT_FALSE(fast_query.isInitialized());
  EXPECT_DOUBLE_EQ(fast_query.evaluate(query, segment_view, tolerances), direct_gwn);

  fast_query.initialize(segment_view);
  EXPECT_TRUE(fast_query.isInitialized());
  EXPECT_NEAR(fast_query.evaluate(query, segment_view, tolerances), direct_gwn, 1.e-4);
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
