// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/core.hpp"
#include "axom/bump/extraction/ClipTableManager.hpp"
#include "axom/bump/extraction/tables/clipping/ClipCases.h"

namespace axom
{
namespace bump
{
namespace extraction
{

void ClipTableManager::loadShape(size_t shape)
{
  using namespace axom::bump::extraction::tables::clipping;
  const auto index = shapeToIndex(shape);
  if(!m_tables[index].isLoaded())
  {
    if(shape == ST_TRI)
    {
      m_tables[index].load(numClipCasesTri,
                           numClipShapesTri,
                           startClipShapesTri,
                           clipShapesTri,
                           clipShapesTriSize,
                           m_allocatorID);
    }
    else if(shape == ST_QUA)
    {
      m_tables[index].load(numClipCasesQua,
                           numClipShapesQua,
                           startClipShapesQua,
                           clipShapesQua,
                           clipShapesQuaSize,
                           m_allocatorID);
    }
    else if(shape == ST_POLY5)
    {
      m_tables[index].load(numClipCasesPoly5,
                           numClipShapesPoly5,
                           startClipShapesPoly5,
                           clipShapesPoly5,
                           clipShapesPoly5Size,
                           m_allocatorID);
    }
    else if(shape == ST_POLY6)
    {
      m_tables[index].load(numClipCasesPoly6,
                           numClipShapesPoly6,
                           startClipShapesPoly6,
                           clipShapesPoly6,
                           clipShapesPoly6Size,
                           m_allocatorID);
    }
    else if(shape == ST_POLY7)
    {
      m_tables[index].load(numClipCasesPoly7,
                           numClipShapesPoly7,
                           startClipShapesPoly7,
                           clipShapesPoly7,
                           clipShapesPoly7Size,
                           m_allocatorID);
    }
    else if(shape == ST_POLY8)
    {
      m_tables[index].load(numClipCasesPoly8,
                           numClipShapesPoly8,
                           startClipShapesPoly8,
                           clipShapesPoly8,
                           clipShapesPoly8Size,
                           m_allocatorID);
    }
    else if(shape == ST_TET)
    {
      m_tables[index].load(numClipCasesTet,
                           numClipShapesTet,
                           startClipShapesTet,
                           clipShapesTet,
                           clipShapesTetSize,
                           m_allocatorID);
    }
    else if(shape == ST_PYR)
    {
      m_tables[index].load(numClipCasesPyr,
                           numClipShapesPyr,
                           startClipShapesPyr,
                           clipShapesPyr,
                           clipShapesPyrSize,
                           m_allocatorID);
    }
    else if(shape == ST_WDG)
    {
      m_tables[index].load(numClipCasesWdg,
                           numClipShapesWdg,
                           startClipShapesWdg,
                           clipShapesWdg,
                           clipShapesWdgSize,
                           m_allocatorID);
    }
    else if(shape == ST_HEX)
    {
      m_tables[index].load(numClipCasesHex,
                           numClipShapesHex,
                           startClipShapesHex,
                           clipShapesHex,
                           clipShapesHexSize,
                           m_allocatorID);
    }
  }
}

}  // end namespace extraction
}  // end namespace bump
}  // end namespace axom
