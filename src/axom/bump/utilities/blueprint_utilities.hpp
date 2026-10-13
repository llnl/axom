// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

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
std::vector<std::string> coordsetAxes(const conduit::Node& n_input);

/*!
 * \brief Validate a vertex field's indexing layout against its topology.
 *
 * \param[in] n_topology The Conduit node containing the topology.
 * \param[in] n_field    The Conduit node containing the vertex-associated field.
 * \param[in] fieldName  The field's name, used only in diagnostics.
 *
 * Bump uses the node ids returned by \c TopologyView::zone() to index flat vertex fields.
 * Strided field layouts require strides, while offsets default to zero if not provided.
 * The field layout must match the topology's node layout, and the topology's i-stride must be 1.
 *
 * \note Without field layout metadata, topology node ids directly index the field values.
 *       The function rejects layout metadata on a topology without \c elements/dims.
 * \note Executes on the host. The metadata can be in device memory.
 */
void validateVertexFieldIndexing(const conduit::Node& n_topology,
                                 const conduit::Node& n_field,
                                 const std::string& fieldName);

}  // end namespace utilities
}  // end namespace bump
}  // end namespace axom
