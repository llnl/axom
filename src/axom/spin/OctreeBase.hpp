// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file OctreeBase.hpp
 * \brief Defines templated OctreeBase class and its inner class BlockIndex
 */

#ifndef AXOM_SPIN_OCTREE_BASE__HPP_
#define AXOM_SPIN_OCTREE_BASE__HPP_

#include "axom/config.hpp"
#include "axom/core/NumericLimits.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"

#include "axom/spin/DenseOctreeLevel.hpp"
#include "axom/spin/OctreeLevel.hpp"
#include "axom/spin/SparseOctreeLevel.hpp"

#include <ostream>

namespace axom
{
namespace spin
{
/**
 * \brief Minimal implementation of a BlockDataType for an OctreeBase.
 *
 * BlockData is default constructible and provides the following functions:
 *  isLeaf(), setInternal(), setNonBlock() and isBlock().
 *
 * \note This implementation uses ones-complement to differentiate between leaf
 * and internal blocks. This has the nice properties that
 *    (a) all internal LeafData have an id with a negative number,
 *    (b) all BlockDatas -- representing leaf and internal blocks
 *        can live in the same index space.
 */
class BlockData
{
  static const int NON_BLOCK = ~0;

public:
  BlockData()
  {
    static int idGenerator = 1;
    m_id = idGenerator++;
  }

  BlockData(const BlockData& other) : m_id(other.m_id) { }
  BlockData& operator=(const BlockData& other)
  {
    m_id = other.m_id;
    return *this;
  }

  /// \brief Predicate to determine if this BlockData represents a leaf block in the tree
  bool isLeaf() const { return m_id >= 0; }

  /// \brief Sets the data associated with this block
  void setData(int blockID) { m_id = blockID; }

  /// Marks the block as not in the octree
  void setNonBlock() { m_id = NON_BLOCK; }

  /// Predicate to check if the associated block is in the octree
  bool isBlock() const { return m_id != NON_BLOCK; }

  /**
   * Returns the normalized form of the id for this BlockData instance
   * \note The normalized form is a non-negative integer.
   */
  int getID() const { return isLeaf() ? m_id : ~m_id; }

  const int& dataIndex() const { return m_id; }

  /// \brief Sets the block type to internal
  void setInternal()
  {
    if(isLeaf())
    {
      m_id = ~m_id;
    }
  }

  /// Equality operator for comparing two BlockData instances
  friend bool operator==(const BlockData& lhs, const BlockData& rhs)
  {
    return lhs.m_id == rhs.m_id;
  }

protected:
  int m_id;
};

/**
 * \class OctreeBase
 *
 * \brief Handles the non-geometric operations for our octree such as
 * refinement, finding the parents and children of a node and determining
 * whether a leaf node exists
 *
 * There are two concepts here:
 * A set of nested dyadic integer grids -- A grid at level n has 2^n cells
 * in each of the DIM dimensions and an adaptive octree defined by a subset
 * of the grid cells from these nested grids. The former is a conceptual aid
 * in the sense that it is implicitly encoded.
 *
 * An octree block at level l is refined by adding its 2^DIM children
 * blocks at level (l+1). The root of the octree covers the entire domains
 * and has no parent. The leaf blocks of the octree have no children.
 * The interior of the children do not overlap, and their union covers that
 * of their parent block. Non-leaf blocks are referred to as 'internal'
 *
 * Requirements for BlockDataType: it must be default constructible and
 * provide an isLeaf() predicate as well as a setInternal() function that
 * changes its state from representing a leaf block to an internal block.
 */
template <int DIM, typename BlockDataType>
class OctreeBase
{
public:
  using CoordType = axom::IndexType;
  using GridPt = primal::Point<CoordType, DIM>;
  using GridVec = primal::Vector<CoordType, DIM>;

  using MAX_LEVEL_SIZE =
    slam::policies::CompileTimeSize<CoordType, axom::numeric_limits<CoordType>::digits>;
  using OctreeLevels = slam::OrderedSet<CoordType, CoordType, MAX_LEVEL_SIZE>;

  using OctreeLevelType = OctreeLevel<DIM, BlockDataType>;
  using LeafIndicesLevelMap = slam::Map<OctreeLevelType*>;

  /**
   * \brief Inner class encapsulating the index of an octree <em>block</em>.
   *
   * Each block index is represented as a point on an integer grid (the minimum
   * point of the block's extent) at a given level of resolution.
   *
   * Each level of resolution is a regular grid with \f$ 2^{level} \f$
   * grid points along each dimension.  The <em>root</em> block (at level 0)
   * covers the entire domain.
   * An octree block at level \f$ \ell \f$ has \f$ 2^{DIM} \f$ <em>children</em>
   * at level \f$ \ell + 1 \f$ covering its domain.
   */
  class BlockIndex
  {
  public:
    enum
    {
      /// The number of children of an octree block (2^D in dimension D)
      NUM_CHILDREN = 1 << DIM,
      /// The number of face neighbors of an octree block (2*D in dimension D)
      NUM_FACE_NEIGHBORS = 2 * DIM
    };

  private:
    using OCTREE_CHILDREN_SIZE = slam::policies::CompileTimeSize<int, NUM_CHILDREN>;
    using OCTREE_FACE_NEIGHBORS_SIZE = slam::policies::CompileTimeSize<int, NUM_FACE_NEIGHBORS>;

  public:
    using ChildIndexSet = slam::OrderedSet<int, int, OCTREE_CHILDREN_SIZE>;
    using FaceNeighborIndexSet = slam::OrderedSet<int, int, OCTREE_FACE_NEIGHBORS_SIZE>;

  public:
    /// \brief Default constructor
    BlockIndex() : m_pt(GridPt()), m_lev(0) { }

    /// \brief Constructor from a point and a level
    BlockIndex(const GridPt& pt, int level) : m_pt(pt), m_lev(level) { }

    /**
     * \brief Accessor for the BlockIndex instance's point
     *
     * \returns const reference to the instance's point
     */
    const GridPt& pt() const { return m_pt; }

    /**
     * \brief Accessor for the BlockIndex instance's point
     *
     * \returns reference to the instance's point
     */
    GridPt& pt() { return m_pt; }

    /**
     * \brief Accessor for the BlockIndex instance's level
     *
     * \returns const reference to the instance's level
     */
    const int& level() const { return m_lev; }

    /**
     * \brief Accessor for the BlockIndex instance's level
     *
     * \returns reference to the instance's level
     */
    int& level() { return m_lev; }

    /// \brief The level of the block index's parent
    int parentLevel() const { return m_lev - 1; }

    /// \brief The level of the block index's child
    int childLevel() const { return m_lev + 1; }

    /// \brief Returns the grid point of the block's parent
    GridPt parentPt() const { return GridPt(m_pt.array() / 2); }

    /**
     * \brief Returns the grid point of the block's child at index childIndex
     *
     * \param [in] childIndex The index of the child whose grid point we are finding
     * \pre \f$ 0 \le childIndex < \f$ Octree::NUM_CHILDREN
     */
    GridPt childPt(int childIndex) const
    {
      SLIC_ASSERT(ChildIndexSet().isValidIndex(childIndex));

      GridPt cPoint;

      // Child is at next level of resolution (multiply by two)
      // and offset according to whether corresponding bit
      // in child index is set or not
      for(int dim = 0; dim < DIM; ++dim)
      {
        cPoint[dim] = (m_pt[dim] << 1) + (childIndex & (CoordType(1) << dim) ? 1 : 0);
      }

      return cPoint;
    }

    /// \brief Returns a grid point at the specified offset from the current block index's point
    GridPt neighborPt(const GridPt& offset) const
    {
      GridPt nPoint(m_pt);
      for(int i = 0; i < DIM; ++i)
      {
        nPoint[i] += offset[i];
      }

      return nPoint;
    }

    /**
     * \brief Returns the parent BlockIndex of this block
     *
     * \note Returns an invalid BlockIndex if we attempt to find
     *       the parent of the root block
     */
    BlockIndex parent() const { return BlockIndex(parentPt(), parentLevel()); }

    /**
     * \brief Returns the child BlockIndex of this block
     *
     * \param [in] childIndex The index of the child whose grid point we are finding
     * \pre \f$ 0 \le childIndex < \f$ Octree::NUM_CHILDREN
     */
    BlockIndex child(int childIndex) const { return BlockIndex(childPt(childIndex), childLevel()); }

    /**
     * \brief Returns the face neighbor grid point of this block
     *
     * \pre 0 <= neighborIndex < 2 * DIM
     * \note The face neighbors indices cycle through the dimensions, two per
     * dimension, e.g. Neighbor 0 is at offset (-1, 0,0,...,0), 
     * neighbor 1 is at offset (1,0,0,..,0) and neighbor 2 is at offset (0,-1, 0, 0, ..., 0) etc...
     */
    BlockIndex faceNeighbor(int neighborIndex) const
    {
      SLIC_ASSERT(FaceNeighborIndexSet().isValidIndex(neighborIndex));

      GridPt offset;
      offset[neighborIndex / 2] = (neighborIndex % 2 == 0) ? -1 : 1;

      return BlockIndex(neighborPt(offset), m_lev);
    }

    bool operator==(const BlockIndex& other) const
    {
      return (m_lev == other.m_lev) && (m_pt == other.m_pt);
    }

    bool operator!=(const BlockIndex& other) const { return !(*this == other); }

    bool operator<(const BlockIndex& other) const
    {
      if(m_lev < other.m_lev)
      {
        return true;
      }
      if(m_lev > other.m_lev)
      {
        return false;
      }

      for(int i = 0; i < DIM; ++i)
      {
        if(m_pt[i] < other.m_pt[i])
        {
          return true;
        }
        if(m_pt[i] > other.m_pt[i])
        {
          return false;
        }
      }

      return false;
    }

    /**
     * \brief Checks the validity of the index.
     *
     * A block index is valid when its level is \f$ \ge 0 \f$ and  it is inBounds
     * 
     * \returns true if the block index is valid, else false
     * \see inBounds()
     */
    bool isValid() const { return (m_lev >= 0) && inBounds(); }

    /**
     * \brief Checks the if the block is in bounds for the level.
     *
     * A block index is in bounds when each coordinate p[i]
     * of its grid point is \f$ 0 \le p[i] < 2^{level()} \f$.
     *
     * \returns true if the block index is inBounds, else false
     */
    bool inBounds() const
    {
      const CoordType maxVal = (CoordType(1) << m_lev) - CoordType(1);
      for(int i = 0; i < DIM; ++i)
      {
        if((m_pt[i] < 0) || (m_pt[i] > maxVal))
        {
          return false;
        }
      }
      return true;
    }

    /**
     * \brief Predicate to determine if the block instance is a descendant of ancestor block
     *
     * \param ancestor The potential ancestor of the block
     * \note A block is an ancestor of another block if neither block is an invalid_index() 
     *   and the block's are equivalent after 0 or more calls to BlockIndex::parent()
     * \return True, if the block instance is a descendant of the ancestor block
     */
    bool isDescendantOf(const BlockIndex& ancestor) const
    {
      const int ancestorLevel = ancestor.level();
      const int levelDiff = m_lev - ancestorLevel;
      if(levelDiff < 0 || m_lev < 0 || ancestorLevel < 0)
      {
        return false;
      }

      BlockIndex blk(*this);
      for(int i = 0; i < levelDiff; ++i)
      {
        blk = blk.parent();
      }

      SLIC_ASSERT(blk.level() == ancestorLevel);
      return blk.pt() == ancestor.pt();
    }

    std::ostream& print(std::ostream& os) const
    {
      os << "{grid pt: " << m_pt << "; level: " << m_lev << "}";

      return os;
    }

    friend std::ostream& operator<<(std::ostream& os, const BlockIndex& block)
    {
      block.print(os);
      return os;
    }

    /**
     * \brief Helper function to generate an invalid block index.
     *
     * \return A new BlockIndex instance blk
     * \post  blk.isValid() will return false
     */
    static BlockIndex invalid_index() { return BlockIndex(GridPt(), -1); }

    /// \brief The number of children that an octree block can have
    static int numChildren() { return ChildIndexSet().size(); }

    /// \brief The number of face neighbors that an octree block can have (ignoring boundaries)
    static int numFaceNeighbors() { return FaceNeighborIndexSet().size(); }

  private:
    GridPt m_pt;
    int m_lev;
  };

private:
  enum
  {
    MAX_DENSE_LEV = 4,
    MAX_SPARSE16_LEV = 16 / DIM,
    MAX_SPARSE32_LEV = 32 / DIM,
    MAX_SPARSE64_LEV = 64 / DIM
  };

  using DenseOctLevType = DenseOctreeLevel<DIM, BlockDataType, std::uint16_t>;
  using Sparse16OctLevType = SparseOctreeLevel<DIM, BlockDataType, std::uint16_t>;
  using Sparse32OctLevType = SparseOctreeLevel<DIM, BlockDataType, std::uint32_t>;
  using Sparse64OctLevType = SparseOctreeLevel<DIM, BlockDataType, std::uint64_t>;
  using SparsePtOctLevType = SparseOctreeLevel<DIM, BlockDataType, GridPt>;

  using DenseOctLevPtr = DenseOctLevType*;
  using Sparse16OctLevPtr = Sparse16OctLevType*;
  using Sparse32OctLevPtr = Sparse32OctLevType*;
  using Sparse64OctLevPtr = Sparse64OctLevType*;
  using SparsePtOctLevPtr = SparsePtOctLevType*;

  /// \brief Simple utility to check if a pointer of type BasePtrType can be cast to a pointer of type DerivedPtrType
  template <typename DerivedPtrType, typename BasePtrType>
  bool checkCast(BasePtrType base) const
  {
    return dynamic_cast<DerivedPtrType>(base) != nullptr;
  }

public:
  /// Sets up an octree containing only the root block
  OctreeBase() : m_leavesLevelMap(&m_levels)
  {
    for(int i = 0; i < maxLeafLevel(); ++i)
    {
      // Use DenseOctreeLevel on first few levels to reduce allocations
      // and fragmentation Use Morton-based SparseOctreeLevel
      // (key is smallest possible integer) on next few levels.
      // Use point bases SparseOctreeLevel (key is Point<int, DIM>,
      // hashed using a MortonIndex)  when MortonIndex requires more than 64
      if(i <= MAX_DENSE_LEV)
      {
        m_leavesLevelMap[i] = new DenseOctLevType(i);
      }
      else if(i <= MAX_SPARSE16_LEV)
      {
        m_leavesLevelMap[i] = new Sparse16OctLevType(i);
      }
      else if(i <= MAX_SPARSE32_LEV)
      {
        m_leavesLevelMap[i] = new Sparse32OctLevType(i);
      }
      else if(i <= MAX_SPARSE64_LEV)
      {
        m_leavesLevelMap[i] = new Sparse64OctLevType(i);
      }
      else
      {
        m_leavesLevelMap[i] = new SparsePtOctLevType(i);
      }
    }

    // Add the root block to the octree
    BlockIndex rootBlock = root();
    (*m_leavesLevelMap[rootBlock.level()]).addAllChildren(rootBlock.pt());
  }

  /// \brief OctreeBase destructor
  ~OctreeBase()
  {
    for(int i = 0; i < maxLeafLevel(); ++i)
    {
      delete m_leavesLevelMap[i];
      m_leavesLevelMap[i] = nullptr;
    }
  }

  /// \brief The max level for leaf blocks of the octree
  int maxLeafLevel() const { return m_levels.size(); }

  /// \brief The max level for internal blocks of the octree
  int maxInternalLevel() const { return m_levels.size() - 1; }

public:
  ///@{

  /**
   * \brief Utility function to find the number of (possible) grid cells at a given level or resolution
   *
   * \param [in] level The level or resolution.
   * \pre \f$ 0 \le lev
   */
  static GridPt maxGridCellAtLevel(int level) { return GridPt(maxCoordAtLevel(level)); }

  /**
   * \brief Finds the highest coordinate value at a given level or resolution
   *
   * \param [in] level The level or resolution.
   * \pre \f$ 0 \le lev
   */
  static CoordType maxCoordAtLevel(int level) { return (CoordType(1) << level) - CoordType(1); }

  /**
   * Auxiliary function to return the root of the octree
   * \note The root block has no parent. Its parent is an invalid BlockIndex.
   *       I.e. octree.parent( octree.root()).isValid() = false.
   */
  static BlockIndex root() { return BlockIndex(); }

  ///@}

public:
  ///@{
  // KW: The following four functions are probably not necessary any more
  //     Since their functionality is in the BlockIndex inner class.

  /**
   * \brief Finds the grid index and level of the current octree block's parent.
   *
   * \note The root level is 0 and its children are at level 1
   * \note The root node has no parent.  The returned level will be '-1'
   * \param [in] pt The grid index of the block whose parent we want to find.
   * \param [in] level The level of the block whose parent we want to find.
   * \return The parent of the provided octree leaf.
   */
  BlockIndex parent(const GridPt& pt, int level) const { return BlockIndex(pt, level).parent(); }

  /**
   * \brief Finds the BlockIndex of the given block's parent.
   *
   * \param [in] block The block whose parent we want to find
   * \return The BlockIndex of the parent of the provided octree leaf.
   */
  BlockIndex parent(const BlockIndex& block) const { return block.parent(); }

  /**
   * \brief Finds the BlockIndex of the given block's child
   *
   * \param [in] pt The grid index of the block whose child we want to find.
   * \param [in] level The level of the block whose child we want to find.
   * \param [in] childIndex The index of the child to find
   * \pre \f$ 0 \le childIndex < 2^{DIM} \f$
   * \return The BlockIndex of the child of the provided octree leaf.
   */
  BlockIndex child(const GridPt& pt, int level, int childIndex) const
  {
    return BlockIndex(pt, level).child(childIndex);
  }

  /**
   * \brief Finds the BlockIndex of the given block's child.
   *
   * \param [in] block The block whose child we want to find
   * \param [in] childIndex The index of the child to find
   * \pre \f$ 0 \le childIndex < 2^{DIM} \f$
   * \return The BlockIndex of the child of the provided octree leaf.
   */
  BlockIndex child(const BlockIndex& block, int childIndex) const
  {
    return block.child(childIndex);
  }

  ///@}

  /// \brief Accessor for a reference to the octree level instance at level lev
  OctreeLevelType& getOctreeLevel(int lev) { return *m_leavesLevelMap[lev]; }

  /// \brief Const accessor for a reference to the octree level instance at level lev
  const OctreeLevelType& getOctreeLevel(int lev) const { return *m_leavesLevelMap[lev]; }

public:
  /**
   * \brief Predicate to determine if level lev is in the range
   *
   * \note lev is in range if 0 <= lev <= maxLeafLevel()
   */
  bool isLevelValid(int lev) const { return lev >= 0 && lev <= maxLeafLevel(); }

  /**
   * \brief Determine whether the octree contains a leaf block associated with
   * grid point pt at level lev
   *
   * \param [in] pt The grid point to check
   * \param [in] lev The level of the grid point
   * \returns true if the associated block is a leaf in the octree, false otherwise
   */
  bool isLeaf(const GridPt& pt, int lev) const
  {
    return isLevelValid(lev) && getOctreeLevel(lev).isLeaf(pt);
  }

  /**
   * \brief Determine whether the octree contains a leaf block associated with this BlockIndex
   *
   * \param [in] block The BlockIndex of the tree to check
   * \returns true if the associated block is a leaf in the octree, false otherwise
   */
  bool isLeaf(const BlockIndex& block) const
  {
    return isLevelValid(block.level()) && getOctreeLevel(block.level()).isLeaf(block.pt());
  }

  /**
   * \brief Determine whether the octree contains an internal block associated
   * with grid point pt at level lev
   *
   * \param [in] pt The grid point to check
   * \param [in] lev The level of the grid point
   * \returns true if the associated block is an internal block of the octree, false otherwise
   */
  bool isInternal(const GridPt& pt, int lev) const
  {
    return isLevelValid(lev) && getOctreeLevel(lev).isInternal(pt);
  }

  /**
   * \brief Determine whether the octree contains an internal block associated with this BlockIndex
   *
   * \param [in] block The BlockIndex of the tree to check
   * \returns true if the associated block is an internal block of the octree,
   * false otherwise
   */
  bool isInternal(const BlockIndex& block) const
  {
    return isLevelValid(block.level()) && getOctreeLevel(block.level()).isInternal(block.pt());
  }

  /**
   * \brief Determine whether the octree contains a block (internal or leaf)
   * associated with grid point pt at level lev
   *
   * \param [in] pt The grid point to check
   * \param [in] lev The level of the grid point
   * \returns true if the associated block is in the octree, false otherwise
   */
  bool hasBlock(const GridPt& pt, int lev) const
  {
    bool ret = false;
    if(!isLevelValid(lev))
    {
      // No-op
    }
    else if(lev <= MAX_DENSE_LEV)
    {
      SLIC_ASSERT(checkCast<DenseOctLevPtr>(m_leavesLevelMap[lev]));
      ret = static_cast<DenseOctLevPtr>(m_leavesLevelMap[lev])->hasBlock(pt);
    }
    else if(lev <= MAX_SPARSE16_LEV)
    {
      SLIC_ASSERT(checkCast<Sparse16OctLevPtr>(m_leavesLevelMap[lev]));
      ret = static_cast<Sparse16OctLevPtr>(m_leavesLevelMap[lev])->hasBlock(pt);
    }
    else if(lev <= MAX_SPARSE32_LEV)
    {
      SLIC_ASSERT(checkCast<Sparse32OctLevPtr>(m_leavesLevelMap[lev]));
      ret = static_cast<Sparse32OctLevPtr>(m_leavesLevelMap[lev])->hasBlock(pt);
    }
    else if(lev <= MAX_SPARSE64_LEV)
    {
      SLIC_ASSERT(checkCast<Sparse64OctLevPtr>(m_leavesLevelMap[lev]));
      ret = static_cast<Sparse64OctLevPtr>(m_leavesLevelMap[lev])->hasBlock(pt);
    }
    else
    {
      SLIC_ASSERT(checkCast<SparsePtOctLevPtr>(m_leavesLevelMap[lev]));
      ret = static_cast<SparsePtOctLevPtr>(m_leavesLevelMap[lev])->hasBlock(pt);
    }
    return ret;
  }

  /**
   * \brief Determine whether the octree contains a block (internal or leaf)
   * associated with this BlockIndex
   *
   * \param [in] block The BlockIndex of the tree to check
   * \returns true if the associated block is a block of the octree, false otherwise
   */
  bool hasBlock(const BlockIndex& block) const { return this->hasBlock(block.pt(), block.level()); }

  /**
   * \brief Determine whether the octree block associated with grid point pt and
   * level lev is a possible block in this octree
   *
   * \note A block index is out of bounds if its level is not in the tree, or
   * its grid point is out of the range of possible grid points for its level
   */
  bool inBounds(const GridPt& pt, int lev) const
  {
    return isLevelValid(lev) && BlockIndex(pt, lev).inBounds();
  }

  /**
   * \brief Determine whether the octree block associated with BlockIndex is a
   * possible block in this octree
   *
   * \note A block index is out of bounds if its level is not in the tree, or
   * its grid point is out of the range of possible grid points for its level
   */
  bool inBounds(const BlockIndex& block) const
  {
    return isLevelValid(block.level()) && block.inBounds();
  }

  /**
   * \brief Refines the given leaf block in the octree
   *
   * Marks leafBlock as internal (non-leaf) and adds its children to the tree
   * \pre leafBlock is a valid leaf block in the octree.
   */
  void refineLeaf(const BlockIndex& leafBlock)
  {
    SLIC_ASSERT(isLeaf(leafBlock));

    // Find the leaf node and set as internal
    OctreeLevelType& currentNodeLevelMap = getOctreeLevel(leafBlock.level());
    currentNodeLevelMap[leafBlock.pt()].setInternal();

    // Add its children to the tree
    OctreeLevelType& childLevelMap = getOctreeLevel(leafBlock.childLevel());
    childLevelMap.addAllChildren(leafBlock.pt());
  }

  /**
   * \brief Accessor to the data associated with block
   *
   * \param block A block (internal or leaf) in the tree
   * \pre block is a leaf in the tree
   */
  BlockDataType& operator[](const BlockIndex& block)
  {
    SLIC_ASSERT_MSG(hasBlock(block), "Block " << block << " was not a block in the tree.");

    const int& lev = block.level();
    const GridPt& pt = block.pt();

    if(lev <= MAX_DENSE_LEV)
    {
      SLIC_ASSERT(checkCast<DenseOctLevPtr>(m_leavesLevelMap[lev]));
      return (*static_cast<DenseOctLevPtr>(m_leavesLevelMap[lev]))[pt];
    }
    else if(lev <= MAX_SPARSE16_LEV)
    {
      SLIC_ASSERT(checkCast<Sparse16OctLevPtr>(m_leavesLevelMap[lev]));
      return (*static_cast<Sparse16OctLevPtr>(m_leavesLevelMap[lev]))[pt];
    }
    else if(lev <= MAX_SPARSE32_LEV)
    {
      SLIC_ASSERT(checkCast<Sparse32OctLevPtr>(m_leavesLevelMap[lev]));
      return (*static_cast<Sparse32OctLevPtr>(m_leavesLevelMap[lev]))[pt];
    }
    else if(lev <= MAX_SPARSE64_LEV)
    {
      SLIC_ASSERT(checkCast<Sparse64OctLevPtr>(m_leavesLevelMap[lev]));
      return (*static_cast<Sparse64OctLevPtr>(m_leavesLevelMap[lev]))[pt];
    }
    else
    {
      SLIC_ASSERT(checkCast<SparsePtOctLevPtr>(m_leavesLevelMap[lev]));
      return (*static_cast<SparsePtOctLevPtr>(m_leavesLevelMap[lev]))[pt];
    }
  }

  /**
   * \brief Const accessor to the data associated with block
   *
   * \param block A block (internal or leaf) in the tree
   * \pre block is a leaf in the tree
   */
  const BlockDataType& operator[](const BlockIndex& block) const
  {
    SLIC_ASSERT_MSG(hasBlock(block), "Block " << block << " was not a block in the tree.");

    const int& lev = block.level();
    const GridPt& pt = block.pt();

    if(lev <= MAX_DENSE_LEV)
    {
      SLIC_ASSERT(checkCast<DenseOctLevPtr>(m_leavesLevelMap[lev]));
      return (*static_cast<DenseOctLevPtr>(m_leavesLevelMap[lev]))[pt];
    }
    else if(lev <= MAX_SPARSE16_LEV)
    {
      SLIC_ASSERT(checkCast<Sparse16OctLevPtr>(m_leavesLevelMap[lev]));
      return (*static_cast<Sparse16OctLevPtr>(m_leavesLevelMap[lev]))[pt];
    }
    else if(lev <= MAX_SPARSE32_LEV)
    {
      SLIC_ASSERT(checkCast<Sparse32OctLevPtr>(m_leavesLevelMap[lev]));
      return (*static_cast<Sparse32OctLevPtr>(m_leavesLevelMap[lev]))[pt];
    }
    else if(lev <= MAX_SPARSE64_LEV)
    {
      SLIC_ASSERT(checkCast<Sparse64OctLevPtr>(m_leavesLevelMap[lev]));
      return (*static_cast<Sparse64OctLevPtr>(m_leavesLevelMap[lev]))[pt];
    }
    else
    {
      SLIC_ASSERT(checkCast<SparsePtOctLevPtr>(m_leavesLevelMap[lev]));
      return (*static_cast<SparsePtOctLevPtr>(m_leavesLevelMap[lev]))[pt];
    }
  }

  /**
   * \brief Finds the finest octree leaf covering BlockIndex blk
   *
   * \param blk A BlockIndex, not necessarily in the octree
   * \param checkInBounds A flag to determine if we should check that
   *        the block lies within the octree bounds (default=true)
   * \post The returned block, if valid, is blk or one of its ancestor blocks.
   * \return The blockIndex of the finest octree leaf covering blk, if it exists, 
   *    BlockIndex::invalid_index otherwise (e.g. blk is an internal block of
   *   the tree or is out of bounds)
   */
  BlockIndex coveringLeafBlock(const BlockIndex& blk, bool checkInBounds = true) const
  {
    // Check that point is in bounds
    if(checkInBounds && !this->inBounds(blk))
    {
      return BlockIndex::invalid_index();
    }

    switch(blockStatus(blk))
    {
    case BlockNotInTree:  // Find its nearest ancestor in the tree (it will be a leaf)
    {
      BlockIndex ancBlk = blk.parent();
      while(!this->hasBlock(ancBlk))
      {
        ancBlk = ancBlk.parent();
      }

      SLIC_ASSERT(this->isLeaf(ancBlk));
      return ancBlk;
    }
    case LeafBlock:  // Already a leaf -- nothing to do
      return blk;
    case InternalBlock:  // An internal block -- no tree leaf can contain it
      return BlockIndex::invalid_index();
    }

    SLIC_ASSERT_MSG(false,
                    "OctreeBase::coveringLeafBlock -- Should never get past the switch statement.  "
                    " Perhaps a new case was added to the TreeBlock enum");
    return BlockIndex::invalid_index();
  }

protected:
  /**
   * \brief Helper function to determine the status of a BlockIndex within an octree instance
   *
   * \note This function is meant to help with implementing basic octree
   * functionality and is not meant to be exposed in the public API
   * \param pt The grid point of the block index that we are testing
   * \param lev The level of the block index that we are testing
   */
  TreeBlockStatus blockStatus(const GridPt& pt, int lev) const
  {
    TreeBlockStatus bStat = BlockNotInTree;
    if(!isLevelValid(lev))
    {
      // No-op
    }
    else if(lev <= MAX_DENSE_LEV)
    {
      SLIC_ASSERT(checkCast<DenseOctLevPtr>(m_leavesLevelMap[lev]));
      bStat = static_cast<DenseOctLevPtr>(m_leavesLevelMap[lev])->blockStatus(pt);
    }
    else if(lev <= MAX_SPARSE16_LEV)
    {
      SLIC_ASSERT(checkCast<Sparse16OctLevPtr>(m_leavesLevelMap[lev]));
      bStat = static_cast<Sparse16OctLevPtr>(m_leavesLevelMap[lev])->blockStatus(pt);
    }
    else if(lev <= MAX_SPARSE32_LEV)
    {
      SLIC_ASSERT(checkCast<Sparse32OctLevPtr>(m_leavesLevelMap[lev]));
      bStat = static_cast<Sparse32OctLevPtr>(m_leavesLevelMap[lev])->blockStatus(pt);
    }
    else if(lev <= MAX_SPARSE64_LEV)
    {
      SLIC_ASSERT(checkCast<Sparse64OctLevPtr>(m_leavesLevelMap[lev]));
      bStat = static_cast<Sparse64OctLevPtr>(m_leavesLevelMap[lev])->blockStatus(pt);
    }
    else
    {
      SLIC_ASSERT(checkCast<SparsePtOctLevPtr>(m_leavesLevelMap[lev]));
      bStat = static_cast<SparsePtOctLevPtr>(m_leavesLevelMap[lev])->blockStatus(pt);
    }
    return bStat;
  }

  /**
   * \brief Helper function to determine the status of a BlockIndex within an octree instance
   *
   * \note This function is meant to help with implementing basic octree functionality
   *       and is not meant to be exposed in the public API
   * \param blk The block index we are testing
   */
  TreeBlockStatus blockStatus(const BlockIndex& blk) const
  {
    return this->blockStatus(blk.pt(), blk.level());
  }

private:
  DISABLE_COPY_AND_ASSIGNMENT(OctreeBase);
  DISABLE_MOVE_AND_ASSIGNMENT(OctreeBase);

protected:
  OctreeLevels m_levels;
  LeafIndicesLevelMap m_leavesLevelMap;
};

}  // end namespace spin
}  // end namespace axom

#endif  // AXOM_SPIN_OCTREE_BASE__HPP_
