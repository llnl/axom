// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file InOutOctreeStats.hpp
 *
 * \brief Defines helper class to generate statistics about an InOutOctree.
 */

#ifndef AXOM_QUEST_INOUT_OCTREE_STATS__HPP_
#define AXOM_QUEST_INOUT_OCTREE_STATS__HPP_

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
class InOutOctreeStats
{
public:
  using InOutOctreeType = InOutOctree<DIM>;
  using CellIndexSet = typename InOutOctreeType::CellIndexSet;

  using OctreeBaseType = typename InOutOctreeType::OctreeBaseType;
  using OctreeLevels = typename OctreeBaseType::OctreeLevels;
  using BlockIndex = typename OctreeBaseType::BlockIndex;

  using LeafCountMap = slam::Map<int>;
  using CellCountMap = slam::Map<int>;
  using CardinalityVCMap = slam::Map<int>;

  using LogHistogram = std::map<int, int>;
  using MinMaxRange = primal::BoundingBox<double, 1>;
  using LengthType = MinMaxRange::PointType;
  using LogRangeMap = std::map<int, MinMaxRange>;

  /** A simple struct to track totals within the octree levels */
  struct Totals
  {
    int blocks {0};
    int leaves {0};
    int leavesWithVert {0};
    int cellRefCount {0};
    int whiteBlocks {0};
    int blackBlocks {0};
    int grayBlocks {0};
  };

public:
  InOutOctreeStats(const InOutOctreeType& octree)
    : m_octree(octree)
    , m_generationState(m_octree.m_generationState)
    , m_levelBlocks(&m_octree.m_levels)
    , m_levelLeaves(&m_octree.m_levels)
    , m_levelLeavesWithVert(&m_octree.m_levels)
    , m_levelCellRefCount(&m_octree.m_levels)
    , m_levelWhiteBlockCount(&m_octree.m_levels)
    , m_levelBlackBlockCount(&m_octree.m_levels)
    , m_levelGrayBlockCount(&m_octree.m_levels)
  {
    if(m_generationState >= InOutOctreeType::INOUTOCTREE_ELEMENTS_INSERTED)
    {
      m_cellCount = CellCountMap(&m_octree.m_meshWrapper.elementSet());
    }

    // Iterate through blocks -- count the numbers of internal and leaf blocks
    for(int lev = 0; lev < m_octree.m_levels.size(); ++lev)
    {
      const auto& levelLeafMap = m_octree.getOctreeLevel(lev);

      m_levelBlocks[lev] = levelLeafMap.numBlocks();
      m_levelLeaves[lev] = levelLeafMap.numLeafBlocks();
      m_levelLeavesWithVert[lev] = 0;
      m_levelCellRefCount[lev] = 0;
      m_levelWhiteBlockCount[lev] = 0;
      m_levelBlackBlockCount[lev] = 0;
      m_levelGrayBlockCount[lev] = 0;

      auto itEnd = levelLeafMap.end();
      for(auto it = levelLeafMap.begin(); it != itEnd; ++it)
      {
        const InOutBlockData& blockData = *it;
        BlockIndex block(it.pt(), lev);

        if(blockData.isLeaf())
        {
          if(blockData.hasData())
          {
            ++m_levelLeavesWithVert[lev];

            if(m_generationState >= InOutOctreeType::INOUTOCTREE_ELEMENTS_INSERTED)
            {
              m_levelCellRefCount[lev] += m_octree.leafCells(block, blockData).size();

              BlockIndex blk(it.pt(), lev);
              CellIndexSet cells = m_octree.leafCells(blk, blockData);
              for(int i = 0; i < cells.size(); ++i)
              {
                ++m_cellCount[cells[i]];
              }
            }
          }

          if(m_generationState >= InOutOctreeType::INOUTOCTREE_LEAVES_COLORED)
          {
            switch(blockData.color())
            {
            case InOutBlockData::Black:
              ++m_levelBlackBlockCount[lev];
              break;
            case InOutBlockData::White:
              ++m_levelWhiteBlockCount[lev];
              break;
            case InOutBlockData::Gray:
              ++m_levelGrayBlockCount[lev];
              break;
            case InOutBlockData::Undetermined:
              break;
            }
          }
        }
      }

      m_totals.blocks += m_levelBlocks[lev];
      m_totals.leaves += m_levelLeaves[lev];
      m_totals.leavesWithVert += m_levelLeavesWithVert[lev];
      m_totals.cellRefCount += m_levelCellRefCount[lev];
      m_totals.whiteBlocks += m_levelWhiteBlockCount[lev];
      m_totals.blackBlocks += m_levelBlackBlockCount[lev];
      m_totals.grayBlocks += m_levelGrayBlockCount[lev];
    }
  }

  /// Generates a string summarizing information about the leaves and blocks of the octree
  std::string blockDataStats() const
  {
    std::stringstream sstr;

    for(int lev = 0; lev < m_octree.m_levels.size(); ++lev)
    {
      if(m_levelBlocks[lev] > 0)
      {
        sstr << fmt::format("\t Level {} has {} blocks -- {} internal; {} leaves ({}% w/ vert);",
                            lev,
                            m_levelBlocks[lev],
                            m_levelBlocks[lev] - m_levelLeaves[lev],
                            m_levelLeaves[lev],
                            integerPercentage(m_levelLeavesWithVert[lev], m_levelLeaves[lev]));

        if(m_generationState >= InOutOctreeType::INOUTOCTREE_LEAVES_COLORED)
        {
          sstr << fmt::format(" Leaf counts: {} Black, {} White, {} Gray w/ {} cell refs.",
                              m_levelBlackBlockCount[lev],
                              m_levelWhiteBlockCount[lev],
                              m_levelGrayBlockCount[lev],
                              m_levelCellRefCount[lev]);
        }
        //sstr <<"Hash load factor: "
        //     << this->m_leavesLevelMap[ lev ].load_factor()
        //     << " -- max lf: " << this->m_leavesLevelMap[ lev
        // ].max_load_factor();
        sstr << "\n";
      }
    }

    return sstr.str();
  }

  /// Generates a string summarizing information about the mesh elements indexed by the octree
  std::string meshDataStats() const
  {
    std::stringstream sstr;

    double meshNumCells = m_octree.m_meshWrapper.numMeshCells();

    sstr << fmt::format(
      "  Mesh has {} vertices."
      "\n  Octree has {} blocks; {} internal; {} leaves ({}% w/ vert); ",
      meshNumCells,
      m_totals.blocks,
      m_totals.blocks - m_totals.leaves,
      m_totals.leaves,
      integerPercentage(m_totals.leavesWithVert, m_totals.leaves));

    if(m_generationState >= InOutOctreeType::INOUTOCTREE_ELEMENTS_INSERTED)
    {
      sstr << fmt::format(
        " \n\t There were {} cell references "
        " (avg. {} refs per cell).",
        m_totals.cellRefCount,
        (m_totals.cellRefCount / meshNumCells));
    }

    return sstr.str();
  }

  std::string cellCountHistogram() const
  {
    std::stringstream sstr;

    // Generate and output a histogram of the bucket counts on a lg-scale
    LogHistogram cellCountHist;  // Create histogram of edge lengths (log scale)
    LogRangeMap cellCountRange;

    int numCells = m_octree.m_meshWrapper.numMeshCells();

    for(int i = 0; i < numCells; ++i)
    {
      LengthType count(m_cellCount[i]);
      int expBase2;
      std::frexp(m_cellCount[i], &expBase2);
      cellCountHist[expBase2]++;
      cellCountRange[expBase2].addPoint(count);
    }

    std::stringstream cellCountStr;
    cellCountStr << "\tCell index count "
                 << "(lg-arithmic bins for number of references per cell):";
    for(auto it = cellCountHist.begin(); it != cellCountHist.end(); ++it)
    {
      cellCountStr << fmt::format("\n\t exp: {}\t count: {}\tRange: {}",
                                  it->first,
                                  it->second,
                                  cellCountRange[it->first]);
    }

    return cellCountStr.str();
  }

  std::string vertexCardinalityHistogram() const
  {
    std::stringstream sstr;

    using CellVertIndices = typename MeshWrapper<DIM>::CellVertIndices;

    // Generate and output histogram of VT relation
    CardinalityVCMap cardVC(&m_octree.m_meshWrapper.vertexSet());

    int numCells = m_octree.m_meshWrapper.numMeshCells();
    for(int i = 0; i < numCells; ++i)
    {
      CellVertIndices cvRel = m_octree.m_meshWrapper.cellVertexIndices(i);
      for(int j = 0; j < cvRel.size(); ++j)
      {
        cardVC[cvRel[j]]++;
      }
    }

    using LinHistogram = std::map<int, int>;
    LinHistogram vtCardHist;
    int numVerts = m_octree.m_meshWrapper.numMeshVertices();
    for(int i = 0; i < numVerts; ++i)
    {
      LengthType count(cardVC[i]);
      vtCardHist[cardVC[i]]++;
    }

    sstr << "\tCardinality VT relation histogram (linear): ";
    for(auto it = vtCardHist.begin(); it != vtCardHist.end(); ++it)
    {
      sstr << fmt::format("\n\t exp: {}\t count: {}", it->first, it->second);
    }

    return sstr.str();
  }

  std::string summaryStats() const
  {
    std::stringstream octreeStatsStr;

    octreeStatsStr << fmt::format(
      "*** {} octree summary *** \n",
      (m_generationState >= InOutOctreeType::INOUTOCTREE_ELEMENTS_INSERTED ? "PM" : "PR"));

    octreeStatsStr << blockDataStats() << "\n" << meshDataStats();

    if(m_generationState >= InOutOctreeType::INOUTOCTREE_LEAVES_COLORED)
    {
      octreeStatsStr << fmt::format("\n\tColors: {} Black, {} White, {} Gray",
                                    m_totals.blackBlocks,
                                    m_totals.whiteBlocks,
                                    m_totals.grayBlocks);
    }

    if(m_generationState >= InOutOctreeType::INOUTOCTREE_ELEMENTS_INSERTED)
    {
      octreeStatsStr << "\n" << cellCountHistogram() << "\n" << vertexCardinalityHistogram();
    }

    return octreeStatsStr.str();
  }

private:
  int integerPercentage(double val, double size) const
  {
    return (size > 0) ? static_cast<int>((100. * val) / size) : 0;
  }

private:
  const InOutOctreeType& m_octree;
  typename InOutOctreeType::GenerationState m_generationState;

  LeafCountMap m_levelBlocks;
  LeafCountMap m_levelLeaves;
  LeafCountMap m_levelLeavesWithVert;
  LeafCountMap m_levelCellRefCount;

  LeafCountMap m_levelWhiteBlockCount;
  LeafCountMap m_levelBlackBlockCount;
  LeafCountMap m_levelGrayBlockCount;

  CellCountMap m_cellCount;

  Totals m_totals;
};

}  // namespace detail
}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_INOUT_OCTREE_STATS__HPP_
