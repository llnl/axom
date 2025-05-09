// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file InOutOctreeValidator.hpp
 *
 * \brief Defines helper class to validate an InOutOctree instance
 */

#ifndef AXOM_QUEST_INOUT_OCTREE_VALIDATOR__HPP_
#define AXOM_QUEST_INOUT_OCTREE_VALIDATOR__HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"

#include "BlockData.hpp"
#include "MeshWrapper.hpp"

#include "axom/fmt.hpp"

namespace axom
{
namespace quest
{
// Predeclare InOutOctree class
template <int DIM>
class InOutOctree;

namespace detail
{
template <int DIM>
class InOutOctreeValidator
{
public:
  using InOutOctreeType = InOutOctree<DIM>;
  using CellIndexSet = typename InOutOctreeType::CellIndexSet;

  using OctreeBaseType = typename InOutOctreeType::OctreeBaseType;
  using OctreeLevels = typename OctreeBaseType::OctreeLevels;
  using BlockIndex = typename OctreeBaseType::BlockIndex;

  using SpacePt = typename InOutOctreeType::SpacePt;
  using VertexIndex = typename InOutOctreeType::VertexIndex;
  using CellIndex = typename InOutOctreeType::CellIndex;
  using CellVertIndices = typename MeshWrapper<DIM>::CellVertIndices;
  using GeometricBoundingBox = typename InOutOctreeType::GeometricBoundingBox;

public:
  InOutOctreeValidator(const InOutOctreeType& octree)
    : m_octree(octree)
    , m_generationState(m_octree.m_generationState)
  { }

  /// This function checks that all leaf blocks in the tree were assigned a color
  void checkAllLeavesColored() const
  {
    SLIC_DEBUG("--Checking that all leaves have a color -- black, white and gray");

    for(int lev = 0; lev < m_octree.m_levels.size(); ++lev)
    {
      checkAllLeavesColoredAtLevel(lev);
    }
  }

  /// This function checks that each leaf block at level \a level has a color
  void checkAllLeavesColoredAtLevel(int level) const
  {
    const auto& levelLeafMap = m_octree.getOctreeLevel(level);
    auto itEnd = levelLeafMap.end();
    for(auto it = levelLeafMap.begin(); it != itEnd; ++it)
    {
      if(!it->isLeaf())
      {
        continue;
      }

      SLIC_ASSERT_MSG(it->isColored(),
                      fmt::format("Error after coloring level {} leaf block {} was not colored.",
                                  level,
                                  fmt::streamed(BlockIndex(it.pt(), level))));
    }
  }

  /// This function checks that each vertex in the mesh is indexed by a leaf block of the octree
  void checkEachVertexIsIndexed() const
  {
    SLIC_DEBUG("--Checking that each vertex is in a leaf block of the tree.");

    const VertexIndex numVertices = m_octree.m_meshWrapper.numMeshVertices();
    for(VertexIndex i = 0; i < numVertices; ++i)
    {
      const SpacePt& pos = m_octree.m_meshWrapper.vertexPosition(i);
      BlockIndex vertBlock = m_octree.findLeafBlock(pos);
      const InOutBlockData& leafData = m_octree[vertBlock];

      VertexIndex vertInBlock = m_octree.leafVertex(vertBlock, leafData);
      AXOM_UNUSED_VAR(vertInBlock);

      // Check that we can find the leaf block indexing each vertex from its position
      SLIC_ASSERT_MSG(
        leafData.hasData() && vertInBlock == i,
        fmt::format(" Vertex {} at position {} was not indexed in block {} with bounding "
                    "box {} (point {} contained in block)."
                    " \n\n *** \n Leaf data: {} \n ***",
                    i,
                    pos,
                    fmt::streamed(vertBlock),
                    m_octree.blockBoundingBox(vertBlock),
                    (m_octree.blockBoundingBox(vertBlock).contains(pos) ? "is" : "is NOT"),
                    leafData));

      // Check that our cached value of the vertex's block is accurate
      SLIC_ASSERT_MSG(vertBlock == m_octree.m_vertexToBlockMap[i],
                      fmt::format("Cached block for vertex {} differs from its found block. "
                                  "\n\t -- cached block {} -- is leaf? {}"
                                  "\n\t -- actual block {} -- is leaf? {}"
                                  "\n\t -- vertex set size: {}",
                                  i,
                                  fmt::streamed(m_octree.m_vertexToBlockMap[i]),
                                  (m_octree[m_octree.m_vertexToBlockMap[i]].isLeaf() ? "yes" : "no"),
                                  fmt::streamed(vertBlock),
                                  (m_octree[vertBlock].isLeaf() ? "yes" : "no"),
                                  numVertices));
    }
  }

  // Check that for each a cell's vertices, that cell is indexed by the block containing that vertex
  void checkCellsReferencedInBoundaryVertexBlocks() const
  {
    SLIC_DEBUG(
      "--Checking that each cell is referenced by the leaf blocks "
      "containing its vertices.");

    const axom::IndexType numCells = m_octree.m_meshWrapper.numMeshCells();
    for(axom::IndexType cIdx = 0; cIdx < numCells; ++cIdx)
    {
      CellVertIndices cvRel = m_octree.m_meshWrapper.cellVertexIndices(cIdx);
      for(int j = 0; j < cvRel.size(); ++j)
      {
        VertexIndex vIdx = cvRel[j];
        BlockIndex vertBlock = m_octree.m_vertexToBlockMap[vIdx];
        const InOutBlockData& leafData = m_octree[vertBlock];

        // Check that this cell is referenced in this block
        bool foundCell = false;
        CellIndexSet leafCells = m_octree.leafCells(vertBlock, leafData);
        for(int k = 0; !foundCell && k < leafCells.size(); ++k)
        {
          if(leafCells[k] == cIdx)
          {
            foundCell = true;
          }
        }

        SLIC_ASSERT_MSG(foundCell,
                        fmt::format("Did not find cell with index {} and vertices [{}] "
                                    "in block {} containing vertex {}"
                                    " \n\n *** "
                                    "\n Leaf data: {} "
                                    "\n\t Cells in block? [{}] "
                                    "\n ***",
                                    cIdx,
                                    fmt::join(cvRel, ", "),
                                    fmt::streamed(vertBlock),
                                    vIdx,
                                    leafData,
                                    fmt::join(leafCells.begin(), leafCells.end(), ", ")));
      }
    }
  }

  /// Checks some indexing contstaints on internal and leaf blocks
  void checkBlockIndexingConsistency() const
  {
    SLIC_DEBUG(
      "--Checking that internal blocks have no data, and that leaves satisfy "
      "all PM conditions");

    const double bb_scale_factor = m_octree.m_boundingBoxScaleFactor;

    for(int lev = 0; lev < m_octree.m_levels.size(); ++lev)
    {
      const auto& levelLeafMap = m_octree.getOctreeLevel(lev);
      auto itEnd = levelLeafMap.end();
      for(auto it = levelLeafMap.begin(); it != itEnd; ++it)
      {
        const BlockIndex block(it.pt(), lev);
        const InOutBlockData& data = *it;

        if(!data.isLeaf())
        {
          SLIC_ASSERT(!data.hasData());  // internal blocks should not have data
        }
        else  // leaf block
        {
          if(data.hasData())
          {
            VertexIndex vIdx = m_octree.leafVertex(block, data);
            CellIndexSet blockCells = m_octree.leafCells(block, data);
            for(int i = 0; i < blockCells.size(); ++i)
            {
              CellIndex cIdx = blockCells[i];

              // Check that vIdx is one of this cell's vertices
              CellVertIndices cvRel = m_octree.m_meshWrapper.cellVertexIndices(cIdx);

              SLIC_ASSERT_MSG(m_octree.m_meshWrapper.incidentInVertex(cvRel, vIdx),
                              fmt::format("All cells in a gray block must be incident "
                                          " in a common vertex, but cell {} with "
                                          "vertices [{}] in block {}"
                                          " is not incident in vertex {}",
                                          cIdx,
                                          fmt::join(cvRel, ", "),
                                          fmt::streamed(block),
                                          vIdx));

              // Check that this cell intersects the bounding box of the block
              GeometricBoundingBox blockBB = m_octree.blockBoundingBox(block);
              blockBB.expand(bb_scale_factor);
              SLIC_ASSERT_MSG(m_octree.blockIndexesElementVertex(cIdx, block) ||
                                intersect(m_octree.m_meshWrapper.cellPositions(cIdx), blockBB),
                              fmt::format("Cell {} was indexed in block {}"
                                          " but it does not intersect the block."
                                          "\n\tBlock bounding box: {}"
                                          "\n\tCell positions: {}"
                                          "\n\tCell vertex indices: [{}]"
                                          "\n\tLeaf vertex is: {}"
                                          "\n\tLeaf cells: {} ({})",
                                          cIdx,
                                          fmt::streamed(block),
                                          blockBB,
                                          m_octree.m_meshWrapper.cellPositions(cIdx),
                                          fmt::join(cvRel, ", "),
                                          vIdx,
                                          fmt::join(blockCells, ", "),
                                          blockCells.size()));
            }
          }
        }
      }
    }
  }

  /// Checks conditions about block colors
  void checkNeighboringBlockColors() const
  {
    SLIC_DEBUG("--Checking that inside blocks do not neighbor outside blocks");
    for(int lev = m_octree.maxLeafLevel() - 1; lev >= 0; --lev)
    {
      const auto& levelLeafMap = m_octree.getOctreeLevel(lev);
      auto itEnd = levelLeafMap.end();
      for(auto it = levelLeafMap.begin(); it != itEnd; ++it)
      {
        const BlockIndex block(it.pt(), lev);
        const InOutBlockData& data = *it;

        if(data.isLeaf() && data.color() != InOutBlockData::Gray)
        {
          // When leaf is inside or outside -- face neighbors need same label
          for(int i = 0; i < block.numFaceNeighbors(); ++i)
          {
            BlockIndex neighborBlk = m_octree.coveringLeafBlock(block.faceNeighbor(i));
            if(neighborBlk != BlockIndex::invalid_index())
            {
              const InOutBlockData& neighborData = m_octree[neighborBlk];
              switch(neighborData.color())
              {
              case InOutBlockData::Black:  // intentional fallthrough
              case InOutBlockData::White:
                SLIC_CHECK_MSG(data.color() == neighborData.color(),
                               fmt::format("Problem at block {}  with data {} "
                                           "--- neighbor {} with data {}."
                                           " Neighboring leaves that are not "
                                           "gray must have the same color.",
                                           fmt::streamed(block),
                                           data,
                                           fmt::streamed(neighborBlk),
                                           neighborData));
                break;
              case InOutBlockData::Gray:  // intentional fallthrough
              case InOutBlockData::Undetermined:
                break;
              }
            }
          }
        }
      }
    }
  }

  void checkValid() const
  {
    // We are assumed to be valid before we insert the vertices
    if(m_generationState < InOutOctreeType::INOUTOCTREE_VERTICES_INSERTED)
    {
      return;
    }

    // Iterate through the tree
    // Internal blocks should not have associated vertex data
    // Leaf block consistency depends on 'color'
    //   Black or White blocks should not have vertex data
    //   Gray blocks should have a vertex reference; not necessarily located within the block
    //   Gray blocks should have one or more cells.
    // All cells should be incident in a common vertex -- the indexed vertex -- if it exists.
    // The total sum of vertices located within a block should equal the number of mesh vertices.

    if(m_generationState > InOutOctreeType::INOUTOCTREE_VERTICES_INSERTED)
    {
      checkEachVertexIsIndexed();
    }

    if(m_generationState >= InOutOctreeType::INOUTOCTREE_ELEMENTS_INSERTED)
    {
      checkCellsReferencedInBoundaryVertexBlocks();
      checkBlockIndexingConsistency();
    }

    if(m_generationState >= InOutOctreeType::INOUTOCTREE_LEAVES_COLORED)
    {
      checkAllLeavesColored();
      checkNeighboringBlockColors();
    }
  }

private:
  const InOutOctreeType& m_octree;
  typename InOutOctreeType::GenerationState m_generationState;
};

}  // namespace detail
}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_INOUT_OCTREE_VALIDATOR__HPP_
