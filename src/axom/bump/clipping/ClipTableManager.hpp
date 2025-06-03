// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_CLIPPING_CLIP_TABLE_MANAGER_HPP_
#define AXOM_MIR_CLIPPING_CLIP_TABLE_MANAGER_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/mir/clipping/ClipCases.h"

namespace axom
{
namespace mir
{
namespace clipping
{
/*!
 * \accelerated
 * \brief This class contains a view of table data and it provides an
 *        iterator for traversing shapes in a case.
 */
class TableView
{
public:
  using IndexData = int;
  using TableData = unsigned char;
  using IndexView = axom::ArrayView<IndexData>;
  using TableDataView = axom::ArrayView<TableData>;

  /*!
   * \brief An iterator for shapes within a table case.
   */
  class iterator
  {
  public:
    /*!
     * \brief Return the index of the iterator's current shape.
     * \return The index of the iterator's current shape.
     */
    AXOM_HOST_DEVICE
    inline int index() const { return m_currentShape; }

    /*!
     * \brief Return the number of shapes in the iterator's case.
     * \return The number of shapes in the iterator's case.
     */
    AXOM_HOST_DEVICE
    inline int size() const { return m_numShapes; }

    /*!
     * \brief Increment the iterator, moving it to the next shape.
     */
    AXOM_HOST_DEVICE
    inline void operator++()
    {
      if(m_currentShape < m_numShapes)
      {
        const TableData *ptr = m_shapeStart + m_offset;
        m_offset += shapeLength(ptr);
        m_currentShape++;
      }
    }

    /*!
     * \brief Increment the iterator, moving it to the next shape.
     */
    AXOM_HOST_DEVICE
    inline void operator++(int)
    {
      if(m_currentShape < m_numShapes)
      {
        const TableData *ptr = m_shapeStart + m_offset;
        m_offset += shapeLength(ptr);
        m_currentShape++;
      }
    }

    /*!
     * \brief Compare 2 iterators for equality.
     * \param it The iterator to be compared to this.
     * \return true if the iterators are equal; false otherwise.
     */
    AXOM_HOST_DEVICE
    inline bool operator==(const iterator &it) const
    {
      // Do not worry about m_offset
      return m_shapeStart == it.m_shapeStart && m_currentShape == it.m_currentShape &&
        m_numShapes == it.m_numShapes;
    }

    /*!
     * \brief Compare 2 iterators to see if not equal.
     * \param it The iterator to be compared to this.
     * \return true if the iterators are different; false otherwise.
     */
    AXOM_HOST_DEVICE
    inline bool operator!=(const iterator &it) const
    {
      // Do not worry about m_offset
      return m_shapeStart != it.m_shapeStart || m_currentShape != it.m_currentShape ||
        m_numShapes != it.m_numShapes;
    }

    /*!
     * \brief Dereference operator that wraps the current shape data in an array
     *        view so the caller can use the shape data.
     */
    AXOM_HOST_DEVICE
    inline TableDataView operator*() const
    {
      TableData *ptr = m_shapeStart + m_offset;
      const auto len = shapeLength(ptr);
      return TableDataView(ptr, len);
    }
#if !defined(AXOM_DEVICE_CODE)
  private:
    void printShape(std::ostream &os, TableData shape) const
    {
      switch(shape)
      {
      case ST_PNT:
        os << "ST_PNT";
        break;
      case ST_TRI:
        os << "ST_TRI";
        break;
      case ST_QUA:
        os << "ST_QUA";
        break;
      case ST_TET:
        os << "ST_TET";
        break;
      case ST_PYR:
        os << "ST_PYR";
        break;
      case ST_WDG:
        os << "ST_WDG";
        break;
      case ST_HEX:
        os << "ST_HEX";
        break;
      }
    }
    void printColor(std::ostream &os, TableData color) const
    {
      switch(color)
      {
      case COLOR0:
        os << "COLOR0";
        break;
      case COLOR1:
        os << "COLOR1";
        break;
      case NOCOLOR:
        os << "NOCOLOR";
        break;
      }
    }
    void printIds(std::ostream &os, const TableData *ids, int n) const
    {
      for(int i = 0; i < n; i++)
      {
        if(/*ids[i] >= P0 &&*/ ids[i] <= P7)
          os << "P" << static_cast<int>(ids[i]);
        else if(ids[i] >= EA && ids[i] <= EL)
        {
          char c = 'A' + (ids[i] - EA);
          os << "E" << c;
        }
        else if(ids[i] >= N0 && ids[i] <= N3)
        {
          os << "N" << static_cast<int>(ids[i] - N0);
        }
        os << " ";
      }
    }

  public:
    void print(std::ostream &os) const
    {
      TableData *ptr = m_shapeStart + m_offset;
      printShape(os, ptr[0]);
      os << " ";
      int offset = 2;
      if(ptr[0] == ST_PNT)
      {
        os << static_cast<int>(ptr[1]);  // point number.
        os << " ";

        printColor(os, ptr[2]);
        os << " ";

        os << static_cast<int>(ptr[3]);  // npts
        os << " ";
        offset = 4;
      }
      else
      {
        printColor(os, ptr[1]);
        os << " ";
      }

      const auto n = shapeLength(ptr) - offset;
      printIds(os, ptr + offset, n);
    }
#endif
  private:
    friend class TableView;

    /*!
     * \brief Given the input shape, return how many values to advance to get to the next shape.
     *
     * \param shape The shape type.
     *
     * \return The number of values to advance.
     */
    AXOM_HOST_DEVICE
    size_t shapeLength(const TableData *caseData) const
    {
      size_t retval = 0;
      const auto shape = caseData[0];
      switch(shape)
      {
      case ST_PNT:
        retval = 4 + caseData[3];
        break;
      case ST_TRI:
        retval = 2 + 3;
        break;
      case ST_QUA:
        retval = 2 + 4;
        break;
      case ST_TET:
        retval = 2 + 4;
        break;
      case ST_PYR:
        retval = 2 + 5;
        break;
      case ST_WDG:
        retval = 2 + 6;
        break;
      case ST_HEX:
        retval = 2 + 8;
        break;
      }
      return retval;
    }

    TableData *m_shapeStart {nullptr};
    int m_offset {0};
    int m_currentShape {0};
    int m_numShapes {0};
  };

  /*!
   * \brief Constructor
   */
  AXOM_HOST_DEVICE
  TableView() : m_shapes(), m_offsets(), m_table() { }

  /*!
   * \brief Constructor
   *
   * \param shapes  The number of shapes in each table case.
   * \param offsets The offsets to each shape case in the \a table.
   * \param table   The table data that contains all cases.
   */
  AXOM_HOST_DEVICE
  TableView(const IndexView &shapes, const IndexView &offsets, const TableDataView &table)
    : m_shapes(shapes)
    , m_offsets(offsets)
    , m_table(table)
  { }

  /*!
   * \brief Return the number of cases for the clipping table.
   *
   * \return The number of cases for the clipping table.
   */
  AXOM_HOST_DEVICE
  size_t size() const { return m_shapes.size(); }

  /*!
   * \brief Return the iterator for the beginning of a case.
   *
   * \param caseId The case whose begin iterator we want.
   * \return The iterator at the begin of the case.
   */
  AXOM_HOST_DEVICE
  iterator begin(size_t caseId) const
  {
    SLIC_ASSERT(static_cast<IndexType>(caseId) < m_shapes.size());
    iterator it;
    it.m_shapeStart = const_cast<TableData *>(m_table.data() + m_offsets[caseId]);
    it.m_offset = 0;
    it.m_currentShape = 0;
    it.m_numShapes = m_shapes[caseId];
    return it;
  }

  /*!
   * \brief Return the iterator for the end of a case.
   *
   * \param caseId The case whose end iterator we want.
   * \return The iterator at the end of the case.
   */
  AXOM_HOST_DEVICE
  iterator end(size_t caseId) const
  {
    SLIC_ASSERT(static_cast<IndexType>(caseId) < m_shapes.size());
    iterator it;
    it.m_shapeStart = const_cast<TableData *>(m_table.data() + m_offsets[caseId]);
    it.m_offset = 0;  // not checked in iterator::operator==
    it.m_currentShape = m_shapes[caseId];
    it.m_numShapes = m_shapes[caseId];
    return it;
  }

private:
  IndexView m_shapes;     // The number of shapes in each case.
  IndexView m_offsets;    // The offset to the case in the table.
  TableDataView m_table;  // The table data that contains the shapes.
};

/*!
 * \brief This class manages data table arrays and can produce a view for the data.
 */
template <typename ExecSpace>
class Table
{
public:
  using IndexData = int;
  using TableData = unsigned char;
  using IndexDataArray = axom::Array<IndexData>;
  using TableDataArray = axom::Array<TableData>;

  /*!
   * \brief Returns whether the table data have been loaded.
   * \return True if the data have been loaded; false otherwise.
   */
  bool isLoaded() const { return m_shapes.size() > 0; }

  /*!
   * \brief Load table data into the arrays, moving data as needed.
   *
   * \param n The number of cases in the clip table.
   * \param shapes The number of shapes produced by clipping cases.
   * \param offsets The offset into the table for each clipping case.
   * \param table The clipping table data.
   * \param tableLen The size of the clipping table data.
   */
  void load(size_t n,
            const IndexData *shapes,
            const IndexData *offsets,
            const TableData *table,
            size_t tableLen)
  {
    const int allocatorID = execution_space<ExecSpace>::allocatorID();

    // Allocate space.
    m_shapes = IndexDataArray(n, n, allocatorID);
    m_offsets = IndexDataArray(n, n, allocatorID);
    m_table = TableDataArray(tableLen, tableLen, allocatorID);

    // Copy data to the arrays.
    axom::copy(m_shapes.data(), shapes, n * sizeof(int));
    axom::copy(m_offsets.data(), offsets, n * sizeof(int));
    axom::copy(m_table.data(), table, tableLen * sizeof(unsigned char));
  }

  /*!
   * \brief Create a view to access the table data.
   *
   * \return A view of the table data.
   */
  TableView view() { return TableView(m_shapes.view(), m_offsets.view(), m_table.view()); }

private:
  IndexDataArray m_shapes;
  IndexDataArray m_offsets;
  TableDataArray m_table;
};

/*!
 * \brief Manage several clipping tables.
 */
template <typename ExecSpace>
class ClipTableManager
{
public:
  static constexpr int NumberOfTables = ST_MAX - ST_MIN;

  /*!
   * \brief Return a reference to the clipping table, which is loaded on demand.
   *
   * \param shape The shape type to be retrieved.
   *
   * \return A reference to the clipping table. 
   */
  Table<ExecSpace> &operator[](size_t shape)
  {
    const size_t index = shapeToIndex(shape);
    SLIC_ASSERT(shape < ST_MAX);
    loadShape(shape);
    return m_tables[index];
  }

  /*!
   * \brief Load tables based on dimension.
   * \param dim The dimension of shapes to load.
   */
  void load(int dim)
  {
    for(const auto shape : shapes(dim))
    {
      loadShape(shape);
    }
  }

  /*!
   * \brief Return a vector of clipping shape ids for the given dimension.
   *
   * \param The spatial dimension.
   *
   * \return A vector of clipping shape ids.
   */
  std::vector<size_t> shapes(int dim) const
  {
    std::vector<size_t> s;
    if(dim == -1 || dim == 2)
    {
      for(const auto value : std::vector<size_t> {ST_TRI, ST_QUA})
      {
        s.push_back(value);
      }
    }
    if(dim == -1 || dim == 3)
    {
      for(const auto value : std::vector<size_t> {ST_TET, ST_PYR, ST_WDG, ST_HEX})
      {
        s.push_back(value);
      }
    }
    return s;
  }

private:
  /*!
   * \brief Turn a shape into an table index.
   *
   * \param shape The shape type ST_XXX.
   *
   * \return An index into the m_tables array.
   */
  constexpr static size_t shapeToIndex(size_t shape) { return shape - ST_MIN; }

  /*!
   * \brief Load the clipping table for a shape.
   *
   * \param shape The shape whose table will be loaded.
   */
  void loadShape(size_t shape)
  {
    const auto index = shapeToIndex(shape);
    if(!m_tables[index].isLoaded())
    {
      if(shape == ST_TRI)
      {
        m_tables[index].load(axom::mir::clipping::visit::numClipCasesTri,
                             axom::mir::clipping::visit::numClipShapesTri,
                             axom::mir::clipping::visit::startClipShapesTri,
                             axom::mir::clipping::visit::clipShapesTri,
                             axom::mir::clipping::visit::clipShapesTriSize);
      }
      else if(shape == ST_QUA)
      {
        m_tables[index].load(axom::mir::clipping::visit::numClipCasesQua,
                             axom::mir::clipping::visit::numClipShapesQua,
                             axom::mir::clipping::visit::startClipShapesQua,
                             axom::mir::clipping::visit::clipShapesQua,
                             axom::mir::clipping::visit::clipShapesQuaSize);
      }
      else if(shape == ST_TET)
      {
        m_tables[index].load(axom::mir::clipping::visit::numClipCasesTet,
                             axom::mir::clipping::visit::numClipShapesTet,
                             axom::mir::clipping::visit::startClipShapesTet,
                             axom::mir::clipping::visit::clipShapesTet,
                             axom::mir::clipping::visit::clipShapesTetSize);
      }
      else if(shape == ST_PYR)
      {
        m_tables[index].load(axom::mir::clipping::visit::numClipCasesPyr,
                             axom::mir::clipping::visit::numClipShapesPyr,
                             axom::mir::clipping::visit::startClipShapesPyr,
                             axom::mir::clipping::visit::clipShapesPyr,
                             axom::mir::clipping::visit::clipShapesPyrSize);
      }
      else if(shape == ST_WDG)
      {
        m_tables[index].load(axom::mir::clipping::visit::numClipCasesWdg,
                             axom::mir::clipping::visit::numClipShapesWdg,
                             axom::mir::clipping::visit::startClipShapesWdg,
                             axom::mir::clipping::visit::clipShapesWdg,
                             axom::mir::clipping::visit::clipShapesWdgSize);
      }
      else if(shape == ST_HEX)
      {
        m_tables[index].load(axom::mir::clipping::visit::numClipCasesHex,
                             axom::mir::clipping::visit::numClipShapesHex,
                             axom::mir::clipping::visit::startClipShapesHex,
                             axom::mir::clipping::visit::clipShapesHex,
                             axom::mir::clipping::visit::clipShapesHexSize);
      }
    }
  }

  axom::StackArray<Table<ExecSpace>, NumberOfTables> m_tables {};
};

}  // end namespace clipping
}  // end namespace mir
}  // end namespace axom

#endif
