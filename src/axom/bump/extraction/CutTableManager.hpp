// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_EXTRACTION_CUT_TABLE_MANAGER_HPP_
#define AXOM_BUMP_EXTRACTION_CUT_TABLE_MANAGER_HPP_

#include "axom/core.hpp"
#include "axom/bump/extraction/TableManager.hpp"

namespace axom
{
namespace bump
{
namespace extraction
{

/*!
 * \brief Manage several cutting tables.
 */
class CutTableManager : public TableManager
{
public:
  /// Return false since the tables cannot generate ST_PNT points.
  AXOM_HOST_DEVICE static constexpr bool generates_points() { return false; }

protected:
  /*!
   * \brief Load the cutting table for a shape.
   *
   * \param shape The shape whose table will be loaded.
   * \param allocatorID The allocatorID to use for allocating memory.
   */
  virtual void loadShape(size_t shape) override;
};

}  // end namespace extraction
}  // end namespace bump
}  // end namespace axom

#endif
