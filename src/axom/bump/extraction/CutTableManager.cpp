// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/core.hpp"
#include "axom/bump/extraction/CutTableManager.hpp"
#include "axom/bump/extraction/tables/cutting/CutCases.hpp"

namespace axom
{
namespace bump
{
namespace extraction
{

void CutTableManager::loadShape(size_t shape)
{
  using namespace axom::bump::extraction::tables::cutting;
  const auto index = shapeToIndex(shape);
  if(!m_tables[index].isLoaded())
  {
    if(shape == ST_TRI)
    {
      m_tables[index].load(numCutCasesTri,
                           numCutShapesTri,
                           startCutShapesTri,
                           cutShapesTri,
                           cutShapesTriSize,
                           m_allocator_id);
    }
    else if(shape == ST_QUA)
    {
      m_tables[index].load(numCutCasesQua,
                           numCutShapesQua,
                           startCutShapesQua,
                           cutShapesQua,
                           cutShapesQuaSize,
                           m_allocator_id);
    }
    else if(shape == ST_POLY5)
    {
      m_tables[index].load(numCutCasesPoly5,
                           numCutShapesPoly5,
                           startCutShapesPoly5,
                           cutShapesPoly5,
                           cutShapesPoly5Size,
                           m_allocator_id);
    }
    else if(shape == ST_POLY6)
    {
      m_tables[index].load(numCutCasesPoly6,
                           numCutShapesPoly6,
                           startCutShapesPoly6,
                           cutShapesPoly6,
                           cutShapesPoly6Size,
                           m_allocator_id);
    }
    else if(shape == ST_POLY7)
    {
      m_tables[index].load(numCutCasesPoly7,
                           numCutShapesPoly7,
                           startCutShapesPoly7,
                           cutShapesPoly7,
                           cutShapesPoly7Size,
                           m_allocator_id);
    }
    else if(shape == ST_POLY8)
    {
      m_tables[index].load(numCutCasesPoly8,
                           numCutShapesPoly8,
                           startCutShapesPoly8,
                           cutShapesPoly8,
                           cutShapesPoly8Size,
                           m_allocator_id);
    }
    else if(shape == ST_TET)
    {
      m_tables[index].load(numCutCasesTet,
                           numCutShapesTet,
                           startCutShapesTet,
                           cutShapesTet,
                           cutShapesTetSize,
                           m_allocator_id);
    }
    else if(shape == ST_PYR)
    {
      m_tables[index].load(numCutCasesPyr,
                           numCutShapesPyr,
                           startCutShapesPyr,
                           cutShapesPyr,
                           cutShapesPyrSize,
                           m_allocator_id);
    }
    else if(shape == ST_WDG)
    {
      m_tables[index].load(numCutCasesWdg,
                           numCutShapesWdg,
                           startCutShapesWdg,
                           cutShapesWdg,
                           cutShapesWdgSize,
                           m_allocator_id);
    }
    else if(shape == ST_HEX)
    {
      m_tables[index].load(numCutCasesHex,
                           numCutShapesHex,
                           startCutShapesHex,
                           cutShapesHex,
                           cutShapesHexSize,
                           m_allocator_id);
    }
  }
}

}  // end namespace extraction
}  // end namespace bump
}  // end namespace axom
