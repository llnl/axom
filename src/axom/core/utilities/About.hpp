// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_UTILITIES_ABOUT_H_
#define AXOM_UTILITIES_ABOUT_H_

#include <ostream>
#include <string>

namespace axom
{
/*!
 * \brief Returns the Git SHA if Axom was built in Git repository, empty if not
 *
 * Note: This will not update unless you re-run CMake between commits.
 */
std::string gitSHA();

/*!
 * \brief Prints info about how Axom was configured and built to stdout
 */
void about();

/*!
 * \brief Prints info about how Axom was configured and built to a stream
 *
 * \param [in,out] oss the target stream where to append the Axom info
 */
void about(std::ostream &oss);

/*!
 * \brief Returns a string consisting of the Axom version.
 *
 * \return string corresponding to the Axom version
 */
std::string getVersion();

}  // end namespace axom

#endif  //  AXOM_UTILITIES_ABOUT_H_
