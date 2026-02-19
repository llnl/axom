// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef CORE_SYSTEM_UTILITIES_H_
#define CORE_SYSTEM_UTILITIES_H_

#include <string>
#include <locale>

namespace axom
{
namespace utilities
{
/**
 * @brief Returns the name of the machine
 *
 * @return The name of the current machine, empty string on failure
 */
std::string getHostName();

/**
 * @brief Returns the name of the current user
 *
 * @return The name of the current user, empty string on failure
 */
std::string getUserName();

/**
 * @brief Returns a valid locale for the current system
 * 
 * @param name The name of the desired locale
 */
std::locale locale(const std::string& name = "en_US.UTF-8");

}  // end namespace utilities
}  // end namespace axom

#endif  //  CORE_SYSTEM_UTILITIES_H_
