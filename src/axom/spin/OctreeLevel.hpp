// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file OctreeLevel.hpp
 * \brief Defines templated OctreeLevel class
 *
 * An OctreeLevel associates data with the integer points on a sparse grid.
 * OctreeLevel is an abstract base class.
 *
 * This file also defines two concrete instantiations:
 * GridPointOctreeLevel uses a GridPoint as a hash table key for its octree
 * blocks  MortonOctreeLevel uses a Morton index (of the given bit width) as a
 * hash key for its octree blocks.
 */

#ifndef AXOM_SPIN_OCTREE_LEVEL__HPP_
#define AXOM_SPIN_OCTREE_LEVEL__HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/spin/Brood.hpp"

#include <iterator>

namespace axom
{
namespace spin
{
/// \brief Helper enumeration for status of a BlockIndex within an OctreeLevel instance
enum TreeBlockStatus
{
  BlockNotInTree, /** Status of blocks that are not in the tree */
  LeafBlock,      /** Status of blocks that are leaves in the tree */
  InternalBlock   /** Status of blocks that are internal to the tree */
};

/**
 * \class
 * \brief An abstract base class to represent a sparse level of blocks within an octree.
 *
 * Each block is associated with an integer grid point whose coordinates
 * have values between 0 and 2^L (where L = this->level() is the encoded level).
 * The OctreeLevel associates data of (templated) type BlockDataType with each such block.
 * \note For efficiency, the data is stored within a brood 
 * (collection of siblings that are created simultaneously).
 * \note BlockDataType must define a predicate function with the signature: bool isLeaf() const;
 */
template <int DIM, typename BlockDataType>
class OctreeLevel
{
public:
  /** The coordinate type of a block in the octree */
  using CoordType = axom::IndexType;

  /**
   * \brief A type for the grid points of the octree.
   * \note CoordType must be an integral type
   */
  using GridPt = primal::Point<CoordType, DIM>;

  enum
  {
    BROOD_SIZE = 1 << DIM
  };

  /// A brood is a collection of sibling blocks that are generated simultaneously
  using BroodData = NumericArray<BlockDataType, BROOD_SIZE>;

  /// Predeclare the BlockIterator type
  template <typename OctreeLevel, typename IterHelper, typename DataType>
  class BlockIterator;

protected:
  /// \brief A virtual base class to help with iteration of an OctreeLevel's blocks
  class BlockIteratorHelper
  {
  public:
    /// Virtual destructor
    virtual ~BlockIteratorHelper() { }

    /// \brief A function to increment to the next Block in the level
    virtual void increment() = 0;

    /// \brief Predicate to determine if two BlockIteratorHelpers are equivalent
    virtual bool equal(const BlockIteratorHelper* other) = 0;

    /// Accessor for the point associated with the current octree block
    virtual GridPt pt() const = 0;

    /// Accessor for the data associated with the current octree block
    virtual BlockDataType* data() = 0;

    /// Const accessor for the data associated with the current octree block
    virtual const BlockDataType* data() const = 0;
  };

  /// \brief A virtual base class to help with constant iteration of an OctreeLevel's blocks
  class ConstBlockIteratorHelper
  {
  public:
    /// Virtual destructor
    virtual ~ConstBlockIteratorHelper() { }

    /// \brief A function to increment to the next Block in the level
    virtual void increment() = 0;

    /// \brief Predicate to determine if two BlockIteratorHelpers are equivalent
    virtual bool equal(const ConstBlockIteratorHelper* other) = 0;

    /// Accessor for the point associated with the current octree block
    virtual GridPt pt() const = 0;

    /// Const accessor for the data associated with the current octree block
    virtual const BlockDataType* data() const = 0;
  };

public:
  using BlockIter = BlockIterator<OctreeLevel, BlockIteratorHelper, BlockDataType>;
  using ConstBlockIter =
    BlockIterator<const OctreeLevel, ConstBlockIteratorHelper, const BlockDataType>;

public:
  /// \brief Constructor of an OctreeLevel at level \a lev
  OctreeLevel(int level = -1) : m_level(level) { }

  /// \brief Virtual destructor of an OctreeLevel
  virtual ~OctreeLevel() { };

  /**
   * \brief Returns the maximum coordinate value in the level
   *
   * \note This is (2^l -1), where L is the current level
   */
  CoordType maxCoord() const { return (CoordType(1) << m_level) - CoordType(1); }

  /**
   * \brief Returns a GridPt whose coordinates are set to maxCoord
   *
   * \sa maxCoord()
   */
  GridPt maxGridCell() const { return GridPt(maxCoord()); }

  /// Accessor for the instance's level
  int level() const { return m_level; }

  /**
   * \brief Predicate to check whether the block associated with the
   * given GridPt pt is an allowed block in the level
   *
   * \param [in] pt The gridpoint of the block to check
   * \note pt is inBounds if each of its coordinates is a non-negative
   * integer less than maxCoord()
   * \sa maxCoord()
   */
  bool inBounds(const GridPt& pt) const
  {
    const CoordType maxVal = maxCoord();
    for(int i = 0; i < DIM; ++i)
    {
      if(pt[i] < 0 || pt[i] > maxVal)
      {
        return false;
      }
    }
    return true;
  }

public:
  /**
   * \class
   * \brief An iterator type for the blocks of an octree level
   *
   * \note Uses a helper class to manage the polymorphic OctreeLevel's
   * iteration. The helper defines the following functions:
   *   \a increment(), \a equal() \a update(), \a pt() and \a data()
   */
  template <typename OctreeLevel, typename IterHelper, typename DataType>
  class BlockIterator
  {
  public:
    using GridPt = typename OctreeLevel::GridPt;
    using iter = BlockIterator<OctreeLevel, IterHelper, DataType>;

    // In C++17, inheriting from std::iterator was deprecated.
    // We provide these typedefs for class BlockIterator to avoid inheriting
    // from std::iterator and causing warnings for those compiling to C++17 or newer.
    using value_type = DataType;
    using difference_type = std::ptrdiff_t;
    using pointer = DataType*;
    using reference = DataType&;
    using iterator_category = std::forward_iterator_tag;

    BlockIterator(OctreeLevel* octLevel, bool begin = false) : m_octLevel(octLevel)
    {
      SLIC_ASSERT(octLevel != nullptr);
      m_iterHelper = octLevel->getIteratorHelper(begin);  // factory function
    }

    ~BlockIterator()
    {
      if(m_iterHelper != nullptr)
      {
        delete m_iterHelper;
        m_iterHelper = nullptr;
      }
    }

    /// \brief A const dereference function to access the data
    DataType& operator*() const { return *m_iterHelper->data(); }

    /// \brief A const pointer dereference function to access the data
    DataType* operator->() const { return m_iterHelper->data(); }

    /// \brief Const accessor for the iterator's current grid point
    GridPt pt() const { return m_iterHelper->pt(); }

    /// \brief Equality test against another iterator
    bool operator==(const iter& other) const
    {
      return (m_octLevel == other.m_octLevel)  // point to same object
        && m_iterHelper->equal(other.m_iterHelper);
    }

    /// \brief Inequality test against another iterator
    bool operator!=(const iter& other) const { return !operator==(other); }

    /// \brief Increment the iterator to the next point
    iter& operator++()
    {
      m_iterHelper->increment();
      return *this;
    }

    /// \brief Increment the iterator to the next point
    iter operator++(int)
    {
      iter res = *this;
      m_iterHelper->increment();
      return res;
    }

  private:
    OctreeLevel* m_octLevel {nullptr};   /// Pointer to the iterator's container class
    IterHelper* m_iterHelper {nullptr};  /// Instance of iterator helper class
  };

public:
  /// \brief Predicate to check whether the block associated with the given GridPt pt is a leaf block
  bool isLeaf(const GridPt& pt) const { return this->blockStatus(pt) == LeafBlock; }

  /**
   * \brief Predicate to check whether the block associated with the
   *  given GridPt pt is an internal block
   */
  bool isInternal(const GridPt& pt) const { return this->blockStatus(pt) == InternalBlock; }

  /// \brief Begin iterator to points and data in tree level
  BlockIter begin() { return BlockIter(this, true); }

  /// \brief Const begin iterator to points and data in tree level
  ConstBlockIter begin() const { return ConstBlockIter(this, true); }

  /// \brief End iterator to points and data in tree level
  BlockIter end() { return BlockIter(this, false); }

  /// \brief Const end iterator to points and data in tree level
  ConstBlockIter end() const { return ConstBlockIter(this, false); }

  /** 
   * \brief Virtual function to check the status of a block
   *  (e.g. Leaf, Internal, NotInTree)
   */
  virtual TreeBlockStatus blockStatus(const GridPt& pt) const = 0;

  /// \brief Virtual predicate to determine if the OctreeLevel is empty
  virtual bool empty() const = 0;

  /** 
   * \brief Virtual predicate to determine if the OctreeLevel has
   * a block with the given grid point pt
   */
  virtual bool hasBlock(const GridPt& pt) const = 0;

  /** 
   * \brief Virtual function to add all children of the given
   *   grid point pt to the OctreeLevel
   */
  virtual void addAllChildren(const GridPt& pt) = 0;

  /// \brief Virtual const accessor for the data associated with grid point pt
  virtual const BlockDataType& operator[](const GridPt& pt) const = 0;

  /// \brief Virtual accessor for the data associated with grid point pt
  virtual BlockDataType& operator[](const GridPt& pt) = 0;

  /** 
   * \brief Virtual accessor for the data associated with
   * all children of the given grid point (i.e. the brood)
   */
  virtual BroodData& getBroodData(const GridPt& pt) = 0;

  /** 
   * \brief Virtual const accessor for the data associated
   * with all children of the given grid point (i.e. the brood)
   */
  virtual const BroodData& getBroodData(const GridPt& pt) const = 0;

  /**
   *  \brief Virtual factory function to create an iterator helper
   *  \param A boolean to determine if the iterator should be
   *  a begin iterator (true) or an end iterator (false)
   */
  virtual BlockIteratorHelper* getIteratorHelper(bool) = 0;

  /**
   * \brief Virtual factory function to create a const iterator helper
   *
   * \param A boolean to determine if the iterator should be
   * a begin iterator (true) or an end iterator (false)
   */
  virtual ConstBlockIteratorHelper* getIteratorHelper(bool) const = 0;

  /// \brief Virtual function to compute the number of blocks (internal and leaf) in the level
  virtual int numBlocks() const = 0;

  /** \brief Virtual function to compute the number of
   * internal blocks in the level
   */
  virtual int numInternalBlocks() const = 0;

  /// \brief Virtual function to compute the number of leaf blocks in the level
  virtual int numLeafBlocks() const = 0;

protected:
  int m_level;
};

}  // end namespace spin
}  // end namespace axom

#endif  // AXOM_SPIN_OCTREE_LEVEL__HPP_
