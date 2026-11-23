// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/config.hpp"
#include "axom/core.hpp"

#include <conduit/conduit.hpp>

namespace axom::bump::utilities
{

//------------------------------------------------------------------------------

/*!
 * \brief Find a referenced Blueprint sibling node using a device-safe read of the
 *        reference string.
 *
 * \param[in] node    The node containing the Blueprint reference child.
 * \param[in] refKey  The name of the reference child, such as `"coordset"` or
 *                    `"topology"`.
 *
 * \return A pointer to the referenced node, or `nullptr` if the reference child
 *         or referenced node does not exist.
 *
 * This mirrors `conduit::blueprint::mesh::utils::find_reference_node()`, but it
 * reads the reference string with `hostReadString()` so it remains safe when the
 * referenced leaf storage is device-resident.
 */
const conduit::Node* findReferenceNode(const conduit::Node& node,
                                       const std::string& refKey);

/*!
 * \brief Find a referenced Blueprint sibling node using a device-safe read of the
 *        reference string.
 *
 * \param[in,out] node  The node containing the Blueprint reference child.
 * \param[in] refKey    The name of the reference child, such as `"coordset"` or
 *                      `"topology"`.
 *
 * \return A pointer to the referenced node, or `nullptr` if the reference child
 *         or referenced node does not exist.
 */
conduit::Node* findReferenceNode(conduit::Node& node, const std::string& refKey);

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
 * Strided field layouts require strides. Missing offsets default to zero.
 * The field layout must match the topology's node layout, and the topology's i-stride must be 1.
 *
 * \note Without field layout metadata, topology node ids directly index the field values.
 *       The function rejects layout metadata on a topology without \c elements/dims.
 * \note Runs on the host and reads metadata from device memory when needed.
 */
void validateVertexFieldIndexing(const conduit::Node& n_topology,
                                 const conduit::Node& n_field,
                                 const std::string& fieldName);

/*!
 * \brief Get the coordset type string.
 *
 * \param n_coordset The coordset node.
 * 
 * \return The string for "type".
 *
 * \note We use hostReadString in case the string is in device memory.
 */
std::string coordsetType(const conduit::Node &n_coordset);

/*!
 * \brief Get the topology type string.
 *
 * \param n_topo The topology node.
 * 
 * \return The string for "type".
 *
 * \note We use hostReadString in case the string is in device memory.
 */
std::string topologyType(const conduit::Node &n_topo);

/*!
 * \brief Get the unstructured topology shape string.
 *
 * \param n_topo The topology node.
 * 
 * \return The string for elements/shape.
 *
 * \note We use hostReadString in case the string is in device memory.
 */
std::string unstructuredTopologyShape(const conduit::Node &n_topo);

}  // end namespace axom::bump::utilities
