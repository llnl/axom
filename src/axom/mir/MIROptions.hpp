// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_MIROPTIONS_HPP_
#define AXOM_MIR_MIROPTIONS_HPP_

#include "axom/mir/Options.hpp"

namespace axom
{
namespace mir
{
/**
 * \brief This class provides a kind of schema over the MIR options, as well
 *        as default values, and some utilities functions.
 */
class MIROptions : public axom::mir::Options
{
public:
  /**
   * \brief Constructor
   *
   * \param options The node that contains the clipping options.
   */
  MIROptions(const conduit::Node &options) : axom::mir::Options(options) { }

  /**
   * \brief Get the name of the matset on which we'll operate.
   * \return The name of the matset.
   */
  std::string matset() const { return options().fetch_existing("matset").as_string(); }

  /**
   * \brief Return the name of the matset to make in the output.
   * \param default_value The name to use if the option is not defined.
   * \return The name of the matset to make in the output.
   */
  std::string matsetName(const std::string &default_value = std::string()) const
  {
    std::string name(default_value.empty() ? matset() : default_value);
    if(options().has_child("matsetName")) name = options().fetch_existing("matsetName").as_string();
    return name;
  }

private:
  /// Access the base class' options.
  const conduit::Node &options() const { return this->m_options; }
};

}  // end namespace mir
}  // end namespace axom

#endif
