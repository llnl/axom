// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file CommandLineUtilities.hpp
 *
 * \brief Defines utilities in support of validating command line input
 */

#ifndef AXOM_CORE_COMMANDLINE_UTILITIES_HPP_
#define AXOM_CORE_COMMANDLINE_UTILITIES_HPP_

#include "axom/config.hpp"
#include "axom/core/utilities/Annotations.hpp"

#ifdef AXOM_USE_CLI11
  #include "axom/CLI11.hpp"
#endif
#include "axom/fmt.hpp"

#include <string>

namespace axom
{
namespace utilities
{
#if defined(AXOM_USE_CLI11) && defined(AXOM_USE_CALIPER)
/// Helper class for CLI11 to validate a caliper \a mode string passed into an axom app
struct CaliperModeValidator : public axom::CLI::Validator
{
  CaliperModeValidator()
  {
    name_ = "MODE";
    func_ = [](const std::string &str) {
      if(str == "help")
      {
        return axom::fmt::format("Valid caliper modes are:\n{}\n",
                                 axom::utilities::annotations::detail::mode_help_string());
      }
      return axom::utilities::annotations::detail::is_mode_valid(str)
        ? std::string("")
        : axom::fmt::format(
            "'{}' invalid caliper mode. "
            "Run with '--caliper help' to see all valid options",
            str);
    };
  }
};

const static CaliperModeValidator ValidCaliperMode;
#endif  // defined(AXOM_USE_CLI11) && defined(AXOM_USE_CALIPER)

}  // namespace utilities
}  // namespace axom

#endif  // AXOM_CORE_COMMANDLINE_UTILITIES_HPP_