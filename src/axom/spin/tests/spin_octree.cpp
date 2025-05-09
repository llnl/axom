// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/spin/OctreeBase.hpp"

#include "axom/slic.hpp"

//------------------------------------------------------------------------------
TEST(spin_octree, topological_octree_parent_child)
{
  SLIC_INFO("*** This test exercises the parent/child relation in spin::OctreeBase");

  static const int DIM = 2;
  using CoordType = int;
  using LeafNodeType = axom::spin::BlockData;

  using OctreeType = axom::spin::OctreeBase<DIM, LeafNodeType>;
  using GridPt = OctreeType::GridPt;
  using BlockIndex = OctreeType::BlockIndex;
  using OctreeChildIndexSet = BlockIndex::ChildIndexSet;

  OctreeType octree;

  BlockIndex rootBlock = octree.root();
  EXPECT_EQ(rootBlock.pt(), GridPt());
  EXPECT_EQ(rootBlock.level(), 0);

  OctreeChildIndexSet childIndexSet;
  int numChildren = childIndexSet.size();
  for(int i = 0; i < numChildren; ++i)
  {
    // Root's children are at level one and have coordinates that are either
    // zero or one
    BlockIndex childBlock = octree.child(rootBlock, childIndexSet[i]);
    EXPECT_EQ(childBlock.level(), 1);
    EXPECT_EQ(childBlock.level(), rootBlock.childLevel());

    int recombineIndex = 0;
    for(int dim = 0; dim < DIM; ++dim)
    {
      CoordType coordVal = childBlock.pt()[dim];
      EXPECT_TRUE(coordVal == 0 || coordVal == 1);

      bool expBit = (childIndexSet[i] & 1 << dim) > 0;
      EXPECT_EQ(coordVal, expBit ? 1 : 0);

      recombineIndex += (coordVal << dim);
    }
    EXPECT_EQ(childIndexSet[i], recombineIndex);

    // Parent of child is self
    EXPECT_EQ(octree.parent(childBlock), rootBlock);
  }
}

//------------------------------------------------------------------------------
TEST(spin_octree, topological_octree_refine)
{
  SLIC_INFO("*** This test exercises the block refinement in spin::OctreeBase"
            << "\nSpecifically, that refining the root block adds "
            << "all its children to the octree.");

  static const int DIM = 3;
  using LeafNodeType = axom::spin::BlockData;

  using OctreeType = axom::spin::OctreeBase<DIM, LeafNodeType>;
  using BlockIndex = OctreeType::BlockIndex;

  OctreeType octree;

  BlockIndex rootBlock = octree.root();

  EXPECT_TRUE(octree.hasBlock(rootBlock));
  EXPECT_TRUE(octree.isLeaf(rootBlock));
  EXPECT_FALSE(octree.isInternal(rootBlock));

  int numChildren = BlockIndex::numChildren();
  for(int i = 0; i < numChildren; ++i)
  {
    EXPECT_FALSE(octree.hasBlock(octree.child(rootBlock, i)));
  }

  octree.refineLeaf(rootBlock);

  EXPECT_TRUE(octree.hasBlock(rootBlock));
  EXPECT_FALSE(octree.isLeaf(rootBlock));
  EXPECT_TRUE(octree.isInternal(rootBlock));

  for(int i = 0; i < numChildren; ++i)
  {
    EXPECT_TRUE(octree.hasBlock(octree.child(rootBlock, i)));
    EXPECT_TRUE(octree.isLeaf(octree.child(rootBlock, i)));
  }
}

//------------------------------------------------------------------------------
TEST(spin_octree, octree_coveringLeafBlocks)
{
  SLIC_INFO("*** This test exercises the coveringLeafBlock function of OctreeBase");

  static const int DIM = 2;
  using LeafNodeType = axom::spin::BlockData;
  using OctreeType = axom::spin::OctreeBase<DIM, LeafNodeType>;
  using BlockIndex = OctreeType::BlockIndex;
  using GridPt = OctreeType::GridPt;

  OctreeType octree;

  BlockIndex rootBlock = octree.root();
  SLIC_INFO("Root block of octree is " << rootBlock);
  EXPECT_EQ(rootBlock, octree.coveringLeafBlock(rootBlock));

  EXPECT_EQ(2 * DIM, rootBlock.numFaceNeighbors());

  // All children blocks of a leaf block have the leaf as their covering block
  int numChildren = BlockIndex::numChildren();
  for(int i = 0; i < numChildren; ++i)
  {
    BlockIndex blk = octree.child(rootBlock, i);

    EXPECT_EQ(rootBlock, octree.coveringLeafBlock(blk));
  }

  // All neighbors of the root block are invalid
  for(int i = 0; i < rootBlock.numFaceNeighbors(); ++i)
  {
    BlockIndex neighborBlk = rootBlock.faceNeighbor(i);

    SLIC_INFO(" Face neighbor " << i << " is " << neighborBlk);

    // The root has no valid neighbors at the same level
    EXPECT_EQ(BlockIndex::invalid_index(), octree.coveringLeafBlock(neighborBlk));
  }

  octree.refineLeaf(rootBlock);

  // No leaf covers a block after it is refined
  EXPECT_EQ(BlockIndex::invalid_index(), octree.coveringLeafBlock(rootBlock));

  // Check neighbors of the root's children
  for(int i = 0; i < numChildren; ++i)
  {
    BlockIndex blk = octree.child(rootBlock, i);

    EXPECT_EQ(blk, octree.coveringLeafBlock(blk));

    // Each child or the root has two valid face neighbors
    int validNeighborCount = 0;
    for(int j = 0; j < blk.numFaceNeighbors(); ++j)
    {
      BlockIndex neighborBlk = blk.faceNeighbor(j);

      SLIC_INFO(" Face neighbor " << j << " is " << neighborBlk);

      if(octree.coveringLeafBlock(neighborBlk) != BlockIndex::invalid_index())
      {
        validNeighborCount++;
      }
    }
    EXPECT_EQ(2, validNeighborCount);
  }

  octree.refineLeaf(octree.child(rootBlock, 0));

  // Iterate through level 2 blocks and find face-neighbors
  int lev = 2;
  for(int i = 0; i < (1 << lev); ++i)
  {
    for(int j = 0; j < (1 << lev); ++j)
    {
      GridPt pt = GridPt::make_point(i, j);
      BlockIndex blk(pt, lev);

      BlockIndex coveringBlock = octree.coveringLeafBlock(blk);

      SLIC_INFO("Covering block of " << blk << " is " << coveringBlock);

      // Every blk has a valid covering block
      EXPECT_NE(BlockIndex::invalid_index(), coveringBlock);

      // .. at a coarser level
      EXPECT_GE(blk.level(), coveringBlock.level());

      for(int k = 0; k < blk.numFaceNeighbors(); ++k)
      {
        BlockIndex neighborBlk = blk.faceNeighbor(k);
        BlockIndex coveringBlk = octree.coveringLeafBlock(neighborBlk);

        SLIC_INFO("\tFace neighbor "
                  << k << " is " << neighborBlk << " -- Covering block is " << coveringBlk
                  << (coveringBlk == BlockIndex::invalid_index() ? " -- invalid_index" : ""));

        if(octree.coveringLeafBlock(neighborBlk) != BlockIndex::invalid_index())
        {
          EXPECT_GE(blk.level(), neighborBlk.level());
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(spin_octree, octree_block_is_descendant)
{
  SLIC_INFO(
    "*** This test exercises the isDescendantOf() function of "
    "OctreeBase::BlockIndex");

  static const int DIM = 2;
  using LeafNodeType = axom::spin::BlockData;
  using OctreeType = axom::spin::OctreeBase<DIM, LeafNodeType>;
  using BlockIndex = OctreeType::BlockIndex;

  OctreeType octree;

  BlockIndex rootBlock = octree.root();
  SLIC_INFO("Root block of octree is " << rootBlock);

  EXPECT_TRUE(rootBlock.isDescendantOf(rootBlock));
  EXPECT_FALSE(rootBlock.isDescendantOf(BlockIndex::invalid_index()));
  EXPECT_FALSE(BlockIndex::invalid_index().isDescendantOf(rootBlock));

  SLIC_INFO("-- checking children of root");
  octree.refineLeaf(rootBlock);
  int numChildren = BlockIndex::numChildren();
  for(int i = 0; i < numChildren; ++i)
  {
    BlockIndex blk = octree.child(rootBlock, i);

    EXPECT_TRUE(blk.isDescendantOf(rootBlock));
    EXPECT_FALSE(rootBlock.isDescendantOf(blk));

    EXPECT_FALSE(blk.isDescendantOf(octree.child(rootBlock, (i + 1) % numChildren)));
  }

  BlockIndex child3 = octree.child(rootBlock, 3);
  octree.refineLeaf(child3);

  for(int i = 0; i < numChildren; ++i)
  {
    BlockIndex grandchild = child3.child(i);
    EXPECT_TRUE(grandchild.isDescendantOf(grandchild));
    EXPECT_TRUE(grandchild.isDescendantOf(child3));
    EXPECT_TRUE(grandchild.isDescendantOf(rootBlock));

    if(i != 3)
    {
      EXPECT_FALSE(grandchild.isDescendantOf(rootBlock.child(i)));
    }

    EXPECT_FALSE(grandchild.isDescendantOf(BlockIndex::invalid_index()));
    EXPECT_FALSE(BlockIndex::invalid_index().isDescendantOf(grandchild));
  }
}

TEST(spin_octree, count_octree_blocks)
{
  static const int DIM = 2;
  using LeafNodeType = axom::spin::BlockData;
  using OctreeType = axom::spin::OctreeBase<DIM, LeafNodeType>;

  OctreeType octree;

  // A default initialized octree should have a single leaf block -- the root
  SLIC_INFO("Check number of blocks in empty octree");
  EXPECT_EQ(1, octree.getOctreeLevel(0).numBlocks());
  EXPECT_EQ(0, octree.getOctreeLevel(0).numInternalBlocks());
  EXPECT_EQ(1, octree.getOctreeLevel(0).numLeafBlocks());

  for(int lev = 1; lev < octree.maxLeafLevel(); ++lev)
  {
    EXPECT_EQ(0, octree.getOctreeLevel(lev).numBlocks());
    EXPECT_EQ(0, octree.getOctreeLevel(lev).numInternalBlocks());
    EXPECT_EQ(0, octree.getOctreeLevel(lev).numLeafBlocks());
  }

  SLIC_INFO("Refine the root and check number of blocks in first three levels");
  octree.refineLeaf(octree.root());
  EXPECT_EQ(1, octree.getOctreeLevel(0).numBlocks());
  EXPECT_EQ(1, octree.getOctreeLevel(0).numInternalBlocks());
  EXPECT_EQ(0, octree.getOctreeLevel(0).numLeafBlocks());

  EXPECT_EQ(4, octree.getOctreeLevel(1).numBlocks());
  EXPECT_EQ(0, octree.getOctreeLevel(1).numInternalBlocks());
  EXPECT_EQ(4, octree.getOctreeLevel(1).numLeafBlocks());

  EXPECT_EQ(0, octree.getOctreeLevel(2).numBlocks());
  EXPECT_EQ(0, octree.getOctreeLevel(2).numInternalBlocks());
  EXPECT_EQ(0, octree.getOctreeLevel(2).numLeafBlocks());

  SLIC_INFO("Refine a child of the root and check blocks in first few levels");
  octree.refineLeaf(octree.root().child(3));
  EXPECT_EQ(1, octree.getOctreeLevel(0).numBlocks());
  EXPECT_EQ(1, octree.getOctreeLevel(0).numInternalBlocks());
  EXPECT_EQ(0, octree.getOctreeLevel(0).numLeafBlocks());

  EXPECT_EQ(4, octree.getOctreeLevel(1).numBlocks());
  EXPECT_EQ(1, octree.getOctreeLevel(1).numInternalBlocks());
  EXPECT_EQ(3, octree.getOctreeLevel(1).numLeafBlocks());

  EXPECT_EQ(4, octree.getOctreeLevel(2).numBlocks());
  EXPECT_EQ(0, octree.getOctreeLevel(2).numInternalBlocks());
  EXPECT_EQ(4, octree.getOctreeLevel(2).numLeafBlocks());

  SLIC_INFO("Refine another child of the root and check blocks in first few levels");
  octree.refineLeaf(octree.root().child(2));
  EXPECT_EQ(4, octree.getOctreeLevel(1).numBlocks());
  EXPECT_EQ(2, octree.getOctreeLevel(1).numInternalBlocks());
  EXPECT_EQ(2, octree.getOctreeLevel(1).numLeafBlocks());

  EXPECT_EQ(8, octree.getOctreeLevel(2).numBlocks());
  EXPECT_EQ(0, octree.getOctreeLevel(2).numInternalBlocks());
  EXPECT_EQ(8, octree.getOctreeLevel(2).numLeafBlocks());

  SLIC_INFO("Refine a grandchild of the root and check blocks in first few levels");
  octree.refineLeaf(octree.root().child(2).child(1));
  EXPECT_EQ(4, octree.getOctreeLevel(1).numBlocks());
  EXPECT_EQ(2, octree.getOctreeLevel(1).numInternalBlocks());
  EXPECT_EQ(2, octree.getOctreeLevel(1).numLeafBlocks());

  EXPECT_EQ(8, octree.getOctreeLevel(2).numBlocks());
  EXPECT_EQ(1, octree.getOctreeLevel(2).numInternalBlocks());
  EXPECT_EQ(7, octree.getOctreeLevel(2).numLeafBlocks());

  EXPECT_EQ(4, octree.getOctreeLevel(3).numBlocks());
  EXPECT_EQ(0, octree.getOctreeLevel(3).numInternalBlocks());
  EXPECT_EQ(4, octree.getOctreeLevel(3).numLeafBlocks());
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
