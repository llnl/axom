// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/bump.hpp"
#include "axom/bump/views/StridedStructuredIndexing.hpp"

//----------------------------------------------------------------------
TEST(bump_views_indexing, strided_structured_indexing_2d)
{
  /*

   x---x---x---x---x---x---x
   |   |   |   |   |   |   |
   x---x---*---*---*---*---x   *=real node, x=ignored node, O=origin node 16
   |   |   |   |   |   |   |
   x---x---*---*---*---*---x
   |   |   |   |   |   |   |
   x---x---O---*---*---*---x
   |   |   |   |   |   |   |
   x---x---x---x---x---x---x
   |   |   |   |   |   |   |
   x---x---x---x---x---x---x

   */
  using Indexing = axom::bump::views::StridedStructuredIndexing<int, 2>;
  using LogicalIndex = typename Indexing::LogicalIndex;
  LogicalIndex dims {4, 3};  // window size in 4*3 elements in 7,6 overall
  LogicalIndex origin {2, 2};
  LogicalIndex stride {1, 7};
  Indexing indexing(dims, origin, stride);

  EXPECT_EQ(indexing.dimension(), 2);
  EXPECT_EQ(indexing.size(), dims[0] * dims[1]);

  // Iterate over local
  for(int j = 0; j < dims[1]; j++)
  {
    for(int i = 0; i < dims[0]; i++)
    {
      LogicalIndex logical {i, j};
      const auto flat = indexing.LogicalIndexToIndex(logical);
      const auto logical2 = indexing.IndexToLogicalIndex(flat);
      EXPECT_EQ(logical, logical2);

      EXPECT_EQ(logical, indexing.GlobalToLocal(indexing.LocalToGlobal(logical)));
      EXPECT_EQ(flat, indexing.GlobalToLocal(indexing.LocalToGlobal(flat)));
    }
  }

  // Iterate over global
  int index = 0;
  for(int j = 0; j < dims[1] + origin[1]; j++)
  {
    for(int i = 0; i < stride[1]; i++, index++)
    {
      LogicalIndex logical {i, j};
      const auto flat = indexing.GlobalToGlobal(logical);

      // flat should start at 0 and increase
      EXPECT_EQ(flat, index);

      // Global flat back to logical.
      LogicalIndex logical2 = indexing.GlobalToGlobal(flat);
      EXPECT_EQ(logical, logical2);

      // If we're in a valid region for the local window, try some other things.
      if(i >= origin[0] && i < origin[0] + dims[0] && j >= origin[1] && j < origin[1] + dims[1])
      {
        const auto logicalLocal = indexing.GlobalToLocal(logical);
        const auto flatLocal = indexing.GlobalToLocal(flat);

        // Flat local back to flat global
        EXPECT_EQ(flat, indexing.LocalToGlobal(flatLocal));

        // Logical local back to logical global
        EXPECT_EQ(logical, indexing.LocalToGlobal(logicalLocal));
      }
    }
  }

  // Check whether these local points exist.
  EXPECT_TRUE(indexing.contains(LogicalIndex {0, 0}));
  EXPECT_TRUE(indexing.contains(LogicalIndex {dims[0] - 1, dims[1] - 1}));
  EXPECT_FALSE(indexing.contains(LogicalIndex {4, 0}));
  EXPECT_FALSE(indexing.contains(LogicalIndex {4, 3}));
}

//----------------------------------------------------------------------
TEST(bump_views_indexing, strided_structured_indexing_3d)
{
  using Indexing3D = axom::bump::views::StridedStructuredIndexing<int, 3>;
  using LogicalIndex = typename axom::bump::views::StridedStructuredIndexing<int, 3>::LogicalIndex;
  LogicalIndex dims {4, 3, 3};  // window size in 4*3*3 elements in 6*5*5 overall
  LogicalIndex origin {2, 2, 2};
  LogicalIndex stride {1, 6, 30};
  Indexing3D indexing(dims, origin, stride);

  EXPECT_EQ(indexing.dimension(), 3);
  EXPECT_EQ(indexing.size(), dims[0] * dims[1] * dims[2]);

  const LogicalIndex logical0_0_0 {0, 0, 0};
  const auto index0_0_0 = indexing.LogicalIndexToIndex(logical0_0_0);
  EXPECT_EQ(index0_0_0, 0);

  const LogicalIndex logical2_2_2 {2, 2, 2};
  const auto index2_2_2 = indexing.LogicalIndexToIndex(logical2_2_2);
  EXPECT_EQ(index2_2_2, 2 + 2 * dims[0] + 2 * dims[0] * dims[1]);

  LogicalIndex logical = indexing.IndexToLogicalIndex(index2_2_2);
  EXPECT_TRUE(logical == logical2_2_2);

  for(int k = 0; k < dims[2]; k++)
  {
    for(int j = 0; j < dims[1]; j++)
    {
      for(int i = 0; i < dims[0]; i++)
      {
        LogicalIndex logical {i, j, k};
        const auto flat = indexing.LogicalIndexToIndex(logical);
        const auto logical2 = indexing.IndexToLogicalIndex(flat);
        EXPECT_EQ(logical, logical2);
      }
    }
  }

  EXPECT_TRUE(indexing.contains(logical0_0_0));
  EXPECT_TRUE(indexing.contains(LogicalIndex {dims[0] - 1, dims[1] - 1, dims[2] - 1}));
  EXPECT_FALSE(indexing.contains(LogicalIndex {4, 0, 0}));
  EXPECT_FALSE(indexing.contains(LogicalIndex {4, 3, 0}));
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();
  return result;
}
