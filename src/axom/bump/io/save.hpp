// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_IO_SAVE_HPP_
#define AXOM_BUMP_IO_SAVE_HPP_

#include <conduit/conduit.hpp>
#include <string>

namespace axom
{
namespace bump
{
namespace io
{
//------------------------------------------------------------------------------
/*!
 * \brief Save a Blueprint mesh to a legacy ASCII VTK file.
 *
 * \param node The node that contains the mesh data.
 * \param path The file path to save.
 *
 * \note This function currently handles only unstructured topos with explicit coordsets.
 */
void save_vtk(const conduit::Node &node, const std::string &path);

}  // end namespace io
}  // end namespace bump
}  // end namespace axom

#endif
