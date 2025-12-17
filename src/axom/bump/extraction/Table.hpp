// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_EXTRACTION_TABLE_HPP_
#define AXOM_BUMP_EXTRACTION_TABLE_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/bump/extraction/ExtractionConstants.hpp"

namespace axom
{
namespace bump
{
namespace extraction
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
      case ST_LIN:
        os << "ST_LIN";
        break;
      case ST_TRI:
        os << "ST_TRI";
        break;
      case ST_QUA:
        os << "ST_QUA";
        break;
      case ST_POLY5:
        os << "ST_POLY5";
        break;
      case ST_POLY6:
        os << "ST_POLY6";
        break;
      case ST_POLY7:
        os << "ST_POLY7";
        break;
      case ST_POLY8:
        os << "ST_POLY8";
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
      case ST_LIN:
        retval = 2 + 2;
        break;
      case ST_TRI:
        retval = 2 + 3;
        break;
      case ST_QUA:
        retval = 2 + 4;
        break;
      case ST_POLY5:
        retval = 2 + 5;
        break;
      case ST_POLY6:
        retval = 2 + 6;
        break;
      case ST_POLY7:
        retval = 2 + 7;
        break;
      case ST_POLY8:
        retval = 2 + 8;
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
   * \brief Return the number of cases for the table.
   *
   * \return The number of cases for the table.
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
   * \param shapes The number of shapes produced by cases.
   * \param offsets The offset into the table for each case.
   * \param table The table data.
   * \param tableLen The size of the table data.
   * \param allocatorID The allocator ID to use when allocating memory.
   */
  void load(size_t n,
            const IndexData *shapes,
            const IndexData *offsets,
            const TableData *table,
            size_t tableLen,
            int allocatorID);
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

}  // end namespace extraction
}  // end namespace bump
}  // end namespace axom

#endif
