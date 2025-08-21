// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_DENSE_OCTREE_LEVEL__HPP_
#define AXOM_SPIN_DENSE_OCTREE_LEVEL__HPP_

#include "axom/config.hpp"
#include "axom/core/Types.hpp"

#include "axom/spin/Brood.hpp"
#include "axom/spin/OctreeLevel.hpp"

namespace axom
{
namespace spin
{
/**
 * \class
 * \brief A representation of a dense OctreeLevel.
 *  that allocates space for all possible blocks at the given level.
 *
 *  It uses a Brood-based organization of the data, where the data for all
 *  octree siblings are stored contiguously, and uses a Morton-based order on
 *  the Broods in the level.  Associated with each block is a BlockDataType
 *  For efficiency, the data is associated with an entire brood, a collection of
 *  siblings that are created simultaneously. In dimension DIM, there are 2^DIM
 *  siblings in a brood.
 *
 *  \see OctreeLevel
 *  \see Brood
 */
template <int DIM, typename BlockDataType, typename MortonIndexType>
class DenseOctreeLevel : public OctreeLevel<DIM, BlockDataType>
{
public:
  using Base = OctreeLevel<DIM, BlockDataType>;
  using GridPt = typename Base::GridPt;
  using BroodData = typename Base::BroodData;
  using BaseBlockIteratorHelper = typename Base::BlockIteratorHelper;
  using ConstBaseBlockIteratorHelper = typename Base::ConstBlockIteratorHelper;

  template <typename OctreeLevelType, typename ParentType>
  class IteratorHelper;

  using IterHelper = IteratorHelper<DenseOctreeLevel, BaseBlockIteratorHelper>;
  using ConstIterHelper = IteratorHelper<const DenseOctreeLevel, ConstBaseBlockIteratorHelper>;

  using BroodType = Brood<GridPt, MortonIndexType>;

public:
  /**
   * \brief Concrete instance of the BlockIteratorHelper class defined in the
   * OctreeLevel base class.
   *
   * \note ParenType must be either BlockIteratorHelper or ConstBlockIteratorHelper,
   *       both are defined in the OctreeLevel base class
   */
  template <typename OctreeLevelType, typename ParentType>
  class IteratorHelper : public ParentType
  {
  public:
    using self = IteratorHelper<OctreeLevelType, ParentType>;
    using BaseBlockItType = ParentType;

    IteratorHelper(OctreeLevelType* octLevel, bool begin)
      : m_octreeLevel(octLevel)
      , m_endIdx(static_cast<MortonIndexType>(octLevel->m_broodCapacity))
      , m_offset(0)
      , m_isLevelZero(octLevel->level() == 0)
    {
      m_currentIdx = begin ? 0 : m_endIdx;

      // Advance the iterator to point to a valid Block
      if(begin && !data()->isBlock())
      {
        increment();
      }
    }

    /// Increment to next block of the level
    void increment()
    {
      // Note, must skip blocks that are not in the tree
      do
      {
        ++m_offset;
        if(m_offset == Base::BROOD_SIZE || m_isLevelZero)
        {
          ++m_currentIdx;
          m_offset = 0;
        }
      } while(m_currentIdx < m_endIdx && !data()->isBlock());
    }

    /// Access to point associated with the block pointed to by the iterator
    GridPt pt() const
    {
      // Reconstruct the grid point from its brood representation
      return BroodType::reconstructGridPt(m_currentIdx, m_offset);
    }

    /// Accessor for data associated with the iterator's block
    BlockDataType* data() { return &m_octreeLevel->m_data[m_currentIdx][m_offset]; }
    /// Const accessor for data associated with the iterator's block
    const BlockDataType* data() const { return &m_octreeLevel->m_data[m_currentIdx][m_offset]; }

    /// \brief Predicate to determine if two block iterators are the same
    bool equal(const BaseBlockItType* other)
    {
      const self* pother = dynamic_cast<const self*>(other);

      return (pother != nullptr) && (m_currentIdx == pother->m_currentIdx)  // iterators are the same
        && (m_offset == pother->m_offset);  // brood indices are the same
    }

  private:
    OctreeLevelType* m_octreeLevel {nullptr};
    MortonIndexType m_currentIdx, m_endIdx;
    int m_offset;
    bool m_isLevelZero;
  };

public:
  /// \brief Default constructor for an octree level
  DenseOctreeLevel(int level = -1) : Base(level), m_blockCount(0)
  {
    if(level < 0)
    {
      m_broodCapacity = 0;
      m_data = nullptr;
    }
    else
    {
      const int rowsize = (level > 0) ? 1 << (level - 1) : 1;
      m_broodCapacity = 1;
      for(int i = 0; i < DIM; ++i)
      {
        m_broodCapacity *= rowsize;
      }

      m_data = new BroodData[m_broodCapacity];
    }

    // Mark all blocks as not in tree
    for(int i = 0; i < m_broodCapacity; ++i)
    {
      BroodData& bd = m_data[i];
      for(int j = 0; j < Base::BROOD_SIZE; ++j)
      {
        bd[j].setNonBlock();
      }
    }
  }

  ~DenseOctreeLevel()
  {
    if(m_data != nullptr)
    {
      delete[] m_data;
      m_data = nullptr;
    }

    m_broodCapacity = 0;
    m_blockCount = 0;
  }

  /**
   * \brief Factory function to return a GridBlockIterHelper for this level
   *
   *  \param begin A boolean to determine if this is to be
   * a begin (true) or end (false) iterator
   */
  BaseBlockIteratorHelper* getIteratorHelper(bool begin) { return new IterHelper(this, begin); }

  /**
   * \brief Factory function to return a ConstGridBlockIterHelper for this level
   *
   * \param begin A boolean to determine if this is to be
   *   a begin (true) or end (false) iterator
   */
  ConstBaseBlockIteratorHelper* getIteratorHelper(bool begin) const
  {
    return new ConstIterHelper(this, begin);
  }

  /**
   * \brief Predicate to check whether the block associated with
   * the given GridPt pt is in the current level
   */
  bool hasBlock(const GridPt& pt) const
  {
    const BroodType brood(pt);
    return m_data[brood.base()][brood.offset()].isBlock();
  }

  /**
   * \brief Adds all children of the given grid point to the octree level
   *
   * \param [in] pt The gridPoint associated with the parent
   * of the children that are being added
   * \pre pt must be in bounds for the level
   * \sa inBounds()
   */
  void addAllChildren(const GridPt& pt)
  {
    SLIC_ASSERT_MSG(this->inBounds(pt),
                    "Problem while inserting children of point "
                      << pt << " into octree level " << this->m_level << ". Point was out of bounds -- "
                      << "each coordinate must be between 0 and " << this->maxCoord() << ".");

    getBroodData(pt) = BroodData();

    // Handle level 0 -- only add the root, mark its 'siblings' as non-blocks
    if(this->level() == 0)
    {
      for(int j = 1; j < Base::BROOD_SIZE; ++j)
      {
        m_data[0][j].setNonBlock();
      }
      ++m_blockCount;
    }
    else
    {
      m_blockCount += Base::BROOD_SIZE;
    }
  }

  /// \brief Accessor for the data associated with pt
  BlockDataType& operator[](const GridPt& pt)
  {
    const BroodType brood(pt);
    return m_data[brood.base()][brood.offset()];
  }

  /// \brief Const accessor for the data associated with pt
  const BlockDataType& operator[](const GridPt& pt) const
  {
    SLIC_ASSERT_MSG(hasBlock(pt),
                    "(" << pt << ", " << this->m_level << ") was not a block in the tree at level.");

    const BroodType brood(pt);
    return m_data[brood.base()][brood.offset()];
  }

  /// \brief Access the data associated with the entire brood
  BroodData& getBroodData(const GridPt& pt)
  {
    return m_data[BroodType::MortonizerType::mortonize(pt)];
  }

  /// \brief Const access to data associated with the entire brood
  const BroodData& getBroodData(const GridPt& pt) const
  {
    return m_data[BroodType::MortonizerType::mortonize(pt)];
  }

  /// \brief Predicate to check if there are any blocks in this octree level
  bool empty() const { return m_blockCount == 0; }

  /// \brief Returns the number of blocks (internal and leaf) in the level
  int numBlocks() const { return m_blockCount; }

  /// \brief Returns the number of internal blocks in the level
  int numInternalBlocks() const { return numBlocks() - numLeafBlocks(); }

  /// \brief Returns the number of leaf blocks in the level
  int numLeafBlocks() const
  {
    if(empty())
    {
      return 0;
    }

    int count = 0;
    for(int i = 0; i < m_broodCapacity; ++i)
    {
      const BroodData& bd = m_data[i];
      for(int j = 0; j < Base::BROOD_SIZE; ++j)
      {
        if(bd[j].isLeaf())
        {
          ++count;
        }
      }
    }
    return count;
  }

  /**
   * \brief Helper function to determine the status of an octree block within this octree level
   *
   * \param pt The grid point of the block index that we are testing
   * \return The status of the grid point pt (e.g. LeafBlock, InternalBlock, ...)
   */
  TreeBlockStatus blockStatus(const GridPt& pt) const
  {
    if(!this->inBounds(pt))
    {
      return BlockNotInTree;
    }

    const BroodType brood(pt);
    const BlockDataType& blockData = m_data[brood.base()][brood.offset()];

    return blockData.isBlock() ? (blockData.isLeaf() ? LeafBlock : InternalBlock) : BlockNotInTree;
  }

private:
  DISABLE_COPY_AND_ASSIGNMENT(DenseOctreeLevel);
  DISABLE_MOVE_AND_ASSIGNMENT(DenseOctreeLevel);

private:
  BroodData* m_data {nullptr};
  int m_broodCapacity;
  int m_blockCount;
};

}  // end namespace spin
}  // end namespace axom

#endif  // AXOM_SPIN_DENSE_OCTREE_LEVEL__HPP_
