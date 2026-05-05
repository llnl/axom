// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/bump/extraction/TableManager.hpp"

namespace axom
{
namespace bump
{
namespace extraction
{

TableManager::TableManager() { m_allocator_id = axom::getDefaultAllocatorID(); }

void TableManager::setAllocatorID(int allocatorID) { m_allocator_id = allocatorID; }

Table &TableManager::operator[](size_t shape)
{
  const size_t index = shapeToIndex(shape);
  SLIC_ASSERT(shape < ST_MAX);
  loadShape(shape);
  return m_tables[index];
}

void TableManager::load(int dim)
{
  for(const auto shape : shapes(dim))
  {
    loadShape(shape);
  }
}

std::vector<size_t> TableManager::shapes(int dim) const
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

}  // end namespace extraction
}  // end namespace bump
}  // end namespace axom
