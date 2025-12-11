// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/bump/extraction/TableManager.hpp"

namespace axom
{
namespace bump
{
namespace extraction
{

void Table::load(size_t n,
                 const IndexData *shapes,
                 const IndexData *offsets,
                 const TableData *table,
                 size_t tableLen,
                 int allocatorID)
{
  // Allocate space.
  m_shapes = IndexDataArray(n, n, allocatorID);
  m_offsets = IndexDataArray(n, n, allocatorID);
  m_table = TableDataArray(tableLen, tableLen, allocatorID);

  // Copy data to the arrays.
  axom::copy(m_shapes.data(), shapes, n * sizeof(IndexData));
  axom::copy(m_offsets.data(), offsets, n * sizeof(IndexData));
  axom::copy(m_table.data(), table, tableLen * sizeof(TableData));
}

}  // end namespace extraction
}  // end namespace bump
}  // end namespace axom
