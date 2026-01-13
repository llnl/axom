// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file BlockData.hpp
 *
 * \brief Defines helper classes for data associated with InOutOctree blocks.
 */

#ifndef AXOM_QUEST_INOUT_OCTREE_BLOCKDATA__HPP_
#define AXOM_QUEST_INOUT_OCTREE_BLOCKDATA__HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include <iostream>
#include <vector>

namespace axom
{
namespace quest
{
/**
 * \brief Compact BlockDataType for an InOutOctree
 *
 * Storage requirement is one integer per block to hold the color of a block
 * and for gray block, the index of the associated triangles
 */
class InOutBlockData
{
  // Some internal constants for keeping tracking of the associated block
  // A block is a leaf block when its m_idx is not INTERNAL_BLOCK
  // Leaf blocks can be uncolored or colored (without additional data)
  //      or m_idx be the index of the data associated with a gray block
  enum
  {
    LEAF_BLOCK_UNCOLORED = -1,
    LEAF_BLOCK_WHITE = -2,
    LEAF_BLOCK_BLACK = -3,
    INTERNAL_BLOCK = -4,
    NON_BLOCK = -5
  };

public:
  enum LeafColor
  {
    Undetermined = -2,
    White = -1,
    Gray = 0,
    Black = 1
  };

public:
  /**
   * \brief Default constructor for an InOutBlockData
   *
   * \note Default constructed instances are assumed to be leaf blocks
   */
  InOutBlockData() : m_idx(LEAF_BLOCK_UNCOLORED) { }

  /// \brief Constructor from a given index
  explicit InOutBlockData(int dataIdx) : m_idx(dataIdx) { }

  /// \brief Copy constructor for an InOutBlockData instance
  InOutBlockData(const InOutBlockData& other) : m_idx(other.m_idx) { }

  /// \brief Assignment operator for an InOutBlockData instance
  InOutBlockData& operator=(const InOutBlockData& other)
  {
    this->m_idx = other.m_idx;
    return *this;
  }

public:  // API for a BlockData
  /**
   * \brief Predicate to determine if the associated block is a leaf
   *
   * \return True, if the block is a leaf, False otherwise
   */
  bool isLeaf() const { return m_idx > INTERNAL_BLOCK; }

  /// \brief Marks the associated block as internal
  void setInternal() { m_idx = INTERNAL_BLOCK; }

  /// \brief Marks the associated block as a non-block (i.e. not in the tree)
  void setNonBlock() { m_idx = NON_BLOCK; }

  /**
   * \brief Predicate to determine if the associated block is in the tree
   *
   * \return True, if the block is in the tree (internal or leaf), False otherwise
   */
  bool isBlock() const { return m_idx != NON_BLOCK; }

public:  // Other functions
  /**
   * Clears the data associated with the block
   * \note This function is currently a no-op
   * */
  void clear()
  {
    // No-op for now -- eventually, will need to do something about the index
  }

  /**
   * Predicate to determine if the associated block has data (i.e. it is a gray block)
   * \return True, if the block has data, False otherwise
   * */
  bool hasData() const { return m_idx >= 0; }

  /// Returns the index of the data associated with the block
  const int& dataIndex() const
  {
    //SLIC_ASSERT(hasData());
    return m_idx;
  }

  /**
   * \brief Sets the block as gray, and provides index of its associated data
   *
   * \param idx The index of the data associated with the gray leaf block
   * \pre The block must be a leaf block
   * \pre The passed in index, idx, must be a non-negative integer
   */
  void setGray(int idx)
  {
    SLIC_ASSERT(isLeaf());
    SLIC_ASSERT(idx >= 0);
    m_idx = idx;
  }

  /// Marks the block as Black (the entire domain is inside the surface)
  void setBlack()
  {
    SLIC_ASSERT(isLeaf());
    m_idx = LEAF_BLOCK_BLACK;
  }

  /// Marks the block as Black (the entire domain is outside the surface)
  void setWhite()
  {
    SLIC_ASSERT(isLeaf());
    m_idx = LEAF_BLOCK_WHITE;
  }

  /// Sets the data associated with the block to the given index idx
  void setData(int idx) { m_idx = idx; }

  /// Marks the block as uncolored
  void setUncoloredLeaf()
  {
    SLIC_ASSERT(isLeaf());
    m_idx = LEAF_BLOCK_UNCOLORED;
  }

  /**
   * \brief Find the 'color' of this LeafBlock
   *
   * 'Black' indicates that the entire block is within the surface
   * 'White' indicates that the entire block is outside the surface
   * 'Gray' indicates that the block intersects the surface geometry
   * Leaves that haven't been colored yet are 'Undetermined'
   */
  LeafColor color() const
  {
    if(hasData())
    {
      return Gray;
    }

    switch(m_idx)
    {
    case LEAF_BLOCK_BLACK:
      return Black;
    case LEAF_BLOCK_WHITE:
      return White;
    case LEAF_BLOCK_UNCOLORED:
      return Undetermined;
    }

    SLIC_ASSERT_MSG(false, "Invalid state in InOuLeafData::color()");
    return Undetermined;
  }

  /** 
   * \brief Predicate to determine if the associated block has a color
   * \return True if the block has a color, false otherwise
   * \sa color()
   */
  bool isColored() const { return color() != Undetermined; }

  /// Friend function to compare equality of two InOutBlockData instances
  friend bool operator==(const InOutBlockData& lhs, const InOutBlockData& rhs)
  {
    return lhs.m_idx == rhs.m_idx;
  }

private:
  int m_idx;
};

/**
 * \brief Free function to print an InOutBlockData to an output stream
 * \param os The output stream to write to
 * \param iob The InOUtBlockData instance that we are writing
 */
inline std::ostream& operator<<(std::ostream& os, const InOutBlockData& iob)
{
  os << "InOutBlockData{"
     << "isLeaf: " << (iob.isLeaf() ? "yes" : "no");

  bool showData = true;

  if(iob.isLeaf())
  {
    os << ", color: ";
    switch(iob.color())
    {
    case InOutBlockData::Gray:
      os << "Gray";
      break;
    case InOutBlockData::White:
      os << "White";
      showData = false;
      break;
    case InOutBlockData::Black:
      os << "Black";
      showData = false;
      break;
    default:
      os << "Undetermined";
      break;
    }
  }

  if(showData)
  {
    os << ", dataIndex: ";
    if(!iob.hasData())
    {
      os << "<no data>";
    }
    else
    {
      os << iob.dataIndex();
    }
  }

  os << "}";

  return os;
}

/**
 * \brief Verbose BlockDataType for an InOutOctree
 *
 * \note Used when generating the octree.
 */
class DynamicGrayBlockData
{
public:
  enum
  {
    NO_VERTEX = -1
  };

  using VertexIndex = axom::IndexType;
  using CellIndex = axom::IndexType;

  using CellList = std::vector<CellIndex>;

public:
  /// \brief Default constructor for an InOutLeafData
  DynamicGrayBlockData() : m_vertIndex(NO_VERTEX), m_isLeaf(true) { }

  /**
   * \brief Constructor for an InOutLeafData
   *
   * \param vInd The index of a vertex
   * (optional; default is to not set a vertex)
   */
  DynamicGrayBlockData(VertexIndex vInd, bool isLeaf) : m_vertIndex(vInd), m_isLeaf(isLeaf) { }

  /// \brief Copy constructor for an DynamicGrayBlockData instance
  DynamicGrayBlockData(const DynamicGrayBlockData& other)
    : m_vertIndex(other.m_vertIndex)
    , m_cells(other.m_cells)
    , m_isLeaf(other.m_isLeaf)
  { }

  /// \brief Assignment operator for an InOutLeafData instance
  DynamicGrayBlockData& operator=(const DynamicGrayBlockData& other)
  {
    this->m_vertIndex = other.m_vertIndex;

    this->m_cells.reserve(other.m_cells.size());
    std::copy(other.m_cells.begin(), other.m_cells.end(), std::back_inserter(this->m_cells));

    this->m_isLeaf = other.m_isLeaf;

    return *this;
  }

  //        /**
  //         * \brief Removes all indexed data from this leaf
  //         */
  //        void clear()
  //        {
  //            m_isLeaf = false;
  //            m_cellIndex = NO_VERTEX;
  //            m_cells.clear();
  //            m_cells = CellList(0);    // reconstruct to deallocate memory
  //        }

  /**
   * \brief Equality operator to determine if two
   * DynamicGrayBlockData instances are equivalent
   */
  friend bool operator==(const DynamicGrayBlockData& lhs, const DynamicGrayBlockData& rhs)
  {
    // Note: We are not checking the contents of the cells array, only the size
    return (lhs.m_vertIndex == rhs.m_vertIndex) && (lhs.m_cells.size() == rhs.m_cells.size()) &&
      lhs.m_isLeaf == rhs.m_isLeaf;
  }

public:  // Functions related to whether this is a leaf
  /// Predicate to determine if the associated block is a leaf in the octree
  bool isLeaf() const { return m_isLeaf; }

  /// Sets a flag to determine whether the associated block is a leaf or internal
  void setLeafFlag(bool isLeaf) { m_isLeaf = isLeaf; }

public:  // Functions related to the associated vertex
  /// \brief Checks whether there is a vertex associated with this leaf
  bool hasVertex() const { return m_vertIndex >= 0; }

  /// Sets the vertex associated with this leaf
  void setVertex(VertexIndex vInd) { m_vertIndex = vInd; }

  /// Clears the associated vertex index
  void clearVertex() { m_vertIndex = NO_VERTEX; }

  /// Accessor for the index of the vertex associated with this leaf
  VertexIndex& vertexIndex() { return m_vertIndex; }

  /// Constant accessor for the index of the vertex associated with this leaf
  const VertexIndex& vertexIndex() const { return m_vertIndex; }

public:  // Functions related to the associated triangles
  /// Check whether this Leaf has any associated triangles
  bool hasCells() const { return !m_cells.empty(); }

  /**
   * \brief Reserves space for a given number of triangles
   * \param count The number of triangles for which to reserve space
   */
  void reserveCells(int count) { m_cells.reserve(count); }

  /// Find the number of triangles associated with this leaf
  int numCells() const { return static_cast<int>(m_cells.size()); }

  /// Associates the surface triangle with the given index with this block
  void addCell(CellIndex tInd) { m_cells.push_back(tInd); }

  /// Returns a const reference to the list of triangle indexes associated with the block
  const CellList& cells() const { return m_cells; }

  /// Returns a reference to the list of triangle indexes associated with the block
  CellList& cells() { return m_cells; }

private:
  VertexIndex m_vertIndex;
  CellList m_cells;
  bool m_isLeaf;
};

/// Free function to print a DynamicGrayBlockData instance to an output stream
inline std::ostream& operator<<(std::ostream& os, const DynamicGrayBlockData& bData)
{
  os << "DynamicGrayBlockData{";

  os << "isLeaf: " << (bData.isLeaf() ? "yes" : "no");

  os << ", vertex: ";
  if(bData.hasVertex())
  {
    os << bData.vertexIndex();
  }
  else
  {
    os << "<none>";
  }

  os << ", cells: ";
  if(bData.hasCells())
  {
    int numCell = bData.numCells();
    os << "(" << numCell << ") {";
    for(int i = 0; i < numCell; ++i)
    {
      os << bData.cells()[i] << ((i == numCell - 1) ? "} " : ",");
    }
  }

  os << "}";

  return os;
}

}  // namespace quest
}  // namespace axom

/// Overload to format a quest::InOutBlockData using fmt
template <>
struct axom::fmt::formatter<axom::quest::InOutBlockData> : ostream_formatter
{ };

/// Overload to format a quest::DynamicGrayBlockData using fmt
template <>
struct axom::fmt::formatter<axom::quest::DynamicGrayBlockData> : ostream_formatter
{ };

#endif  // AXOM_QUEST_INOUT_OCTREE_BLOCKDATA__HPP_
