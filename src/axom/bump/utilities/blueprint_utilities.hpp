// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_BLUEPRINT_UTILITIES_HPP_
#define AXOM_BUMP_BLUEPRINT_UTILITIES_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace bump
{
namespace utilities
{

//------------------------------------------------------------------------------

/*!
 * \brief Return the names of the axes for a coordset.
 *
 * \param n_input A Conduit node containing a coordset.
 *
 * \return A vector containing the names of the coordset's axes.
 */
std::vector<std::string> coordsetAxes(const conduit::Node &n_input);

}  // end namespace utilities
}  // end namespace bump
}  // end namespace axom

#endif
