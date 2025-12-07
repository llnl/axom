// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_EXTRACTION_TABLE_MANAGER_HPP_
#define AXOM_BUMP_EXTRACTION_TABLE_MANAGER_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/bump/extraction/Table.hpp"
#include <vector>

namespace axom
{
namespace bump
{
namespace extraction
{

/*!
 * \brief Manage several tables.
 */
template <typename ExecSpace>
class TableManager
{
public:
  static constexpr int NumberOfTables = ST_MAX - ST_MIN;

  /*!
   * \brief Return a reference to the table, which is loaded on demand.
   *
   * \param shape The shape type to be retrieved.
   *
   * \return A reference to the table. 
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
   * \brief Return a vector of shape ids for the given dimension.
   *
   * \param The spatial dimension.
   *
   * \return A vector of shape ids.
   */
  std::vector<size_t> shapes(int dim) const
  {
    std::vector<size_t> s;
    if(dim == -1 || dim == 2)
    {
      for(const auto value :
          std::vector<size_t> {ST_TRI, ST_QUA, ST_POLY5, ST_POLY6, ST_POLY7, ST_POLY8})
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

protected:
  /*!
   * \brief Turn a shape into an table index.
   *
   * \param shape The shape type ST_XXX.
   *
   * \return An index into the m_tables array.
   */
  constexpr static size_t shapeToIndex(size_t shape)
  {
    return shape - ST_MIN;
  }

  /*!
   * \brief Load the table for a shape.
   *
   * \param shape The shape whose table will be loaded.
   */
  virtual void loadShape(size_t shape) = 0;

protected:
  axom::StackArray<Table<ExecSpace>, NumberOfTables> m_tables {};
};

}  // end namespace extraction
}  // end namespace bump
}  // end namespace axom

#endif
