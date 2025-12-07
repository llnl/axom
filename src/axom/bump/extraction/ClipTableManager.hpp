// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_EXTRACTION_CLIP_TABLE_MANAGER_HPP_
#define AXOM_BUMP_EXTRACTION_CLIP_TABLE_MANAGER_HPP_

#include "axom/core.hpp"
#include "axom/bump/extraction/TableManager.hpp"
#include "axom/bump/extraction/tables/clipping/ClipCases.h"

namespace axom
{
namespace bump
{
namespace extraction
{

/*!
 * \brief Manage several clipping tables.
 */
template <typename ExecSpace>
class ClipTableManager : public TableManager<ExecSpace>
{
protected:
  /*!
   * \brief Load the clipping table for a shape.
   *
   * \param shape The shape whose table will be loaded.
   */
  virtual void loadShape(size_t shape) override
  {
    using namespace axom::bump::extraction::tables::clipping;
    auto &tables = TableManager<ExecSpace>::m_tables;
    const auto index = TableManager<ExecSpace>::shapeToIndex(shape);
    if(!tables[index].isLoaded())
    {
      if(shape == ST_TRI)
      {
        tables[index].load(numClipCasesTri,
                             numClipShapesTri,
                             startClipShapesTri,
                             clipShapesTri,
                             clipShapesTriSize);
      }
      else if(shape == ST_QUA)
      {
        tables[index].load(numClipCasesQua,
                             numClipShapesQua,
                             startClipShapesQua,
                             clipShapesQua,
                             clipShapesQuaSize);
      }
      else if(shape == ST_POLY5)
      {
        tables[index].load(numClipCasesPoly5,
                             numClipShapesPoly5,
                             startClipShapesPoly5,
                             clipShapesPoly5,
                             clipShapesPoly5Size);
      }
      else if(shape == ST_POLY6)
      {
        tables[index].load(numClipCasesPoly6,
                             numClipShapesPoly6,
                             startClipShapesPoly6,
                             clipShapesPoly6,
                             clipShapesPoly6Size);
      }
      else if(shape == ST_POLY7)
      {
        tables[index].load(numClipCasesPoly7,
                             numClipShapesPoly7,
                             startClipShapesPoly7,
                             clipShapesPoly7,
                             clipShapesPoly7Size);
      }
      else if(shape == ST_POLY8)
      {
        tables[index].load(numClipCasesPoly8,
                             numClipShapesPoly8,
                             startClipShapesPoly8,
                             clipShapesPoly8,
                             clipShapesPoly8Size);
      }
      else if(shape == ST_TET)
      {
        tables[index].load(numClipCasesTet,
                             numClipShapesTet,
                             startClipShapesTet,
                             clipShapesTet,
                             clipShapesTetSize);
      }
      else if(shape == ST_PYR)
      {
        tables[index].load(numClipCasesPyr,
                             numClipShapesPyr,
                             startClipShapesPyr,
                             clipShapesPyr,
                             clipShapesPyrSize);
      }
      else if(shape == ST_WDG)
      {
        tables[index].load(numClipCasesWdg,
                             numClipShapesWdg,
                             startClipShapesWdg,
                             clipShapesWdg,
                             clipShapesWdgSize);
      }
      else if(shape == ST_HEX)
      {
        tables[index].load(numClipCasesHex,
                             numClipShapesHex,
                             startClipShapesHex,
                             clipShapesHex,
                             clipShapesHexSize);
      }
    }
  }
};

}  // end namespace extraction
}  // end namespace bump
}  // end namespace axom

#endif
