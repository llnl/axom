// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/primal/operators/detail/intersect_bounding_box_impl.hpp"

#include "axom/slic.hpp"
#include "gtest/gtest.h"

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(primal_bounding_box_intersect, aabb_aabb_adjacent)
{
  constexpr double LO = 0.0;
  constexpr double HI = 1.0;
  constexpr int NUM_OFFSETS = 3;
  const double OFFSETS[] = {-1.0, 0.0, 1.0};

  using namespace axom::primal::detail;

  // Check 2D adjacent and self-intersecting bounding boxes for intersection
  for(int i = 0; i < NUM_OFFSETS; i++)
  {
    for(int j = 0; j < NUM_OFFSETS; j++)
    {
      EXPECT_TRUE(intersect_bounding_box(LO + OFFSETS[i],
                                         HI + OFFSETS[i],
                                         LO + OFFSETS[j],
                                         HI + OFFSETS[j],
                                         LO,
                                         HI,
                                         LO,
                                         HI));
    }
  }

  // Check 3D adjacent and self-intersecting bounding boxes for intersection
  for(int i = 0; i < NUM_OFFSETS; i++)
  {
    for(int j = 0; j < NUM_OFFSETS; j++)
    {
      for(int k = 0; k < NUM_OFFSETS; k++)
      {
        EXPECT_TRUE(intersect_bounding_box(LO + OFFSETS[i],
                                           HI + OFFSETS[i],
                                           LO + OFFSETS[j],
                                           HI + OFFSETS[j],
                                           LO + OFFSETS[k],
                                           HI + OFFSETS[k],
                                           LO,
                                           HI,
                                           LO,
                                           HI,
                                           LO,
                                           HI));
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_bounding_box_intersect, aabb_aabb_non_intersecting)
{
  constexpr double LO = 0.0;
  constexpr double HI = 1.0;
  constexpr int NUM_OFFSETS = 3;
  const double OFFSETS[] = {-1.01, 0.0, 1.01};

  using namespace axom::primal::detail;

  // Check 2D bounding boxes for non-intersection
  for(int i = 0; i < NUM_OFFSETS; i++)
  {
    for(int j = 0; j < NUM_OFFSETS; j++)
    {
      // Ignore identity box
      if(i != 1 && j != 1)
      {
        EXPECT_FALSE(intersect_bounding_box(LO + OFFSETS[i],
                                            HI + OFFSETS[i],
                                            LO + OFFSETS[j],
                                            HI + OFFSETS[j],
                                            LO,
                                            HI,
                                            LO,
                                            HI));
      }
    }
  }

  // Check 3D bounding boxes for non-intersection
  for(int i = 0; i < NUM_OFFSETS; i++)
  {
    for(int j = 0; j < NUM_OFFSETS; j++)
    {
      for(int k = 0; k < NUM_OFFSETS; k++)
      {
        // Ignore identity box
        if(i != 1 && j != 1)
        {
          EXPECT_FALSE(intersect_bounding_box(LO + OFFSETS[i],
                                              HI + OFFSETS[i],
                                              LO + OFFSETS[j],
                                              HI + OFFSETS[j],
                                              LO + OFFSETS[j],
                                              HI + OFFSETS[j],
                                              LO,
                                              HI,
                                              LO,
                                              HI,
                                              LO,
                                              HI));
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
