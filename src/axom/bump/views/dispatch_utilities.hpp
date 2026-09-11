// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include <conduit/conduit_node.hpp>

namespace axom
{
namespace bump
{
namespace views
{
/// \name Dimension selection utilities
/// \brief \ref select_dimensions converts dimension indices to flags and uses
/// \ref encode_dimensions to combine them into a mask; \ref dimension_selected
/// queries the resulting mask.
/// @{

/*!
 * \brief Combines encoded dimension flags into a single bit mask.
 *
 * \tparam Dimensions The types of the encoded dimension flags.
 *
 * \param[in] dims The encoded dimension flags to combine.
 *
 * \return The bitwise union of the supplied flags.
 */
template <typename... Dimensions>
constexpr int encode_dimensions(Dimensions... dims)
{
  return (... | dims);
}

/*!
 * \brief Encodes a list of selected spatial dimensions as a bit mask.
 *
 * The bit at each supplied dimension index is set in the returned mask. The
 * mask can be passed as a template argument to limit which dimensions a view
 * dispatcher instantiates.
 *
 * \tparam Dimensions The types of the dimension indices.
 *
 * \param[in] dims The dimension indices to select.
 *
 * \return A bit mask encoding the selected dimensions.
 */
template <typename... Dimensions>
constexpr int select_dimensions(Dimensions... dims)
{
  return encode_dimensions((1 << dims)...);
}

/*!
 * \brief Determines whether a dimension is present in an encoded dimension mask.
 *
 * \param[in] encoded_dims A bit mask returned by \ref select_dimensions.
 * \param[in] dim The dimension index to query.
 *
 * \return True when \a dim is selected; otherwise, false.
 */
constexpr bool dimension_selected(int encoded_dims, int dim) { return encoded_dims & (1 << dim); }

/// @}

/*!
 * \brief Call Blueprint mesh verify functions and convert the output to SLIC_ERROR
 *        if the verify method failed.
 *
 * \param obj The node that contains the object being checked.
 * \param protocol The name of the item to check in the mesh. If the string is empty,
 *                 \a obj node is treated as a mesh and it all gets checked.
 */
void verify(const conduit::Node& obj, const std::string& protocol = std::string());

}  // end namespace views
}  // end namespace bump
}  // end namespace axom
