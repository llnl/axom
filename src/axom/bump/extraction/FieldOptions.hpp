// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BUMP_FIELD_OPTIONS_HPP_
#define AXOM_BUMP_FIELD_OPTIONS_HPP_

#include "axom/bump/extraction/ExtractorOptions.hpp"

namespace axom
{
namespace bump
{
namespace extraction
{
/**
 * \brief This class provides a kind of schema over the clipping options, as well
 *        as default values, and some utilities functions.
 */
class FieldOptions : public axom::bump::extraction::ExtractorOptions
{
public:
  /**
   * \brief Constructor
   *
   * \param options The node that contains the clipping options.
   */
  FieldOptions(const conduit::Node &options) : axom::bump::extraction::ExtractorOptions(options) { }

  /**
   * \brief Return the name of the field used for clipping.
   * \return The name of the field used for clipping.
   */
  std::string field() const { return options().fetch_existing("field").as_string(); }

  /**
   * \brief Return the clip value.
   * \return The clip value.
   */
  float value() const
  {
    return options().has_child("value") ? options().fetch_existing("value").to_float() : 0.f;
  }
};

}  // end namespace extraction
}  // end namespace bump
}  // end namespace axom

#endif
