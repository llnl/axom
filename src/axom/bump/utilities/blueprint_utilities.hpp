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
 * \brief Verify that topology node ids can index a vertex-associated field.
 *
 * \param[in] n_topology The Conduit node containing the topology.
 * \param[in] n_field    The Conduit node containing the vertex-associated field.
 * \param[in] fieldName  The field's name, used only in diagnostics.
 *
 * Bump kernels read a vertex field as a flat array indexed by the node ids from \c TopologyView::zone().
 * When a field supplies Blueprint \c offsets or \c strides, each array must match
 * the corresponding topology array under \c elements/dims.
 *
 * Bump does not support fields with independent layout metadata.
 * This function reports mismatches before a kernel reads the wrong values.
 *
 * \c StructuredTopologyView::zone() assumes an i-stride of 1 when it constructs adjacent corner ids,
 * so this function rejects other values.
 *
 * \note A field with neither offsets nor strides passes this check to preserve existing behavior.
 */
void validateVertexFieldIndexing(const conduit::Node& n_topology,
                                 const conduit::Node& n_field,
                                 const std::string& fieldName);

}  // end namespace utilities
}  // end namespace bump
}  // end namespace axom
