// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BUMP_EXTRACTOR_OPTIONS_HPP_
#define AXOM_BUMP_EXTRACTOR_OPTIONS_HPP_

#include "axom/bump/Options.hpp"

namespace axom
{
namespace bump
{
namespace extraction
{
/**
 * \brief This class provides a kind of schema over the extractor options, as well
 *        as default values, and some utilities functions.
 */
class ExtractorOptions : public axom::bump::Options
{
public:
  /**
   * \brief Constructor
   *
   * \param options The node that contains the clipping options.
   */
  ExtractorOptions(const conduit::Node &options) : axom::bump::Options(options) { }

  /**
   * \brief Return the name of the field used for clipping.
   * \return The name of the field used for clipping.
   */
  std::string clipField() const { return options().fetch_existing("clipField").as_string(); }

  /**
   * \brief Return the clip value.
   * \return The clip value.
   */
  float clipValue() const
  {
    return options().has_child("clipValue") ? options().fetch_existing("clipValue").to_float() : 0.f;
  }

  /**
   * \brief Return the name of the new color field to be created.
   * \return The name of the new color field to be created.
   */
  std::string colorField() const
  {
    std::string name("color");
    if(options().has_child("colorField")) name = options().fetch_existing("colorField").as_string();
    return name;
  }

  /**
   * \brief Return the name of the new original elements field to be created.
   * \return The name of the new original elements to be created.
   */
  std::string originalElementsField() const
  {
    std::string name("originalElements");
    if(options().has_child("originalElementsField"))
      name = options().fetch_existing("originalElementsField").as_string();
    return name;
  }

  /**
   * \brief Return the name of the new nodes field to be created. If the name
   *        is not set then it will not be created.
   * \return The name of the new field to be created.
   */
  std::string newNodesField() const
  {
    std::string name;
    if(options().has_child("newNodesField"))
      name = options().fetch_existing("newNodesField").as_string();
    return name;
  }

  /**
   * \brief Whether the "inside" of the clipping field is selected.
   * \return 1 of the inside clipping is selected, false otherwise.
   */
  bool inside() const
  {
    return options().has_path("inside") ? (options().fetch_existing("inside").to_int() > 0) : true;
  }

  /**
   * \brief Whether the "outside" of the clipping field is selected.
   * \return 1 of the outside clipping is selected, false otherwise.
   */
  bool outside() const
  {
    return options().has_path("outside") ? (options().fetch_existing("outside").to_int() > 0) : false;
  }

protected:
  /// Access the base class' options.
  const conduit::Node &options() const { return this->m_options; }
};

}  // end namespace extraction
}  // end namespace bump
}  // end namespace axom

#endif
