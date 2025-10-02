// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/spin/SpatialOctree.hpp"

#include "axom/slic.hpp"

TEST(spin_spatial_octree, spatial_octree_point_location)
{
  SLIC_INFO("*** This test verifies that a query point falls into " << " a child block.");

  static const int DIM = 3;
  using LeafNodeType = axom::spin::BlockData;

  using OctreeType = axom::spin::SpatialOctree<DIM, LeafNodeType>;
  using BlockIndex = OctreeType::BlockIndex;
  using SpacePt = OctreeType::SpacePt;
  using GeometricBoundingBox = OctreeType::GeometricBoundingBox;

  GeometricBoundingBox bb(SpacePt(10), SpacePt(20));

  // Generate a point within the bounding box
  double alpha = 2. / 3.;
  SpacePt queryPt = SpacePt::lerp(bb.getMin(), bb.getMax(), alpha);
  EXPECT_TRUE(bb.contains(queryPt));

  OctreeType octree(bb);

  // Check that the point lies in a leaf of the tree
  // and that this is the root of the tree
  BlockIndex leafBlock = octree.findLeafBlock(queryPt);
  EXPECT_TRUE(octree.isLeaf(leafBlock));
  EXPECT_EQ(octree.root(), leafBlock);

  GeometricBoundingBox leafBB = octree.blockBoundingBox(leafBlock);
  EXPECT_TRUE(leafBB.contains(queryPt));
  EXPECT_TRUE(bb.contains(leafBB));

  SLIC_INFO("Query pt: " << queryPt << "\n\t" << (leafBB.contains(queryPt) ? " was" : " was NOT")
                         << " contained in bounding box " << leafBB << "\n\t of octree root "
                         << leafBlock);

  for(int i = 0; i < octree.maxInternalLevel(); ++i)
  {
    EXPECT_EQ(0, octree.getOctreeLevel(i + 1).numLeafBlocks());
    octree.refineLeaf(leafBlock);
    EXPECT_EQ(1 << DIM, octree.getOctreeLevel(i + 1).numLeafBlocks());

    leafBlock = octree.findLeafBlock(queryPt);
    EXPECT_TRUE(octree.isLeaf(leafBlock));

    leafBB = octree.blockBoundingBox(leafBlock);
    EXPECT_TRUE(leafBB.contains(queryPt));
    EXPECT_TRUE(bb.contains(leafBB));

    SLIC_INFO("Level " << i << " -- Query pt: " << queryPt << "\n\t"
                       << (leafBB.contains(queryPt) ? " was" : " was not")
                       << " contained in bounding box " << leafBB << "\n\t of leaf " << leafBlock
                       << " in the octree. ");
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
