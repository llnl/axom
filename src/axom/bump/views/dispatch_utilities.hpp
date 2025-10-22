// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_DISPATCH_UTILITIES_HPP_
#define AXOM_BUMP_DISPATCH_UTILITIES_HPP_

#include <conduit/conduit_node.hpp>

namespace axom
{
namespace bump
{
namespace views
{
template <typename... Dimensions>
constexpr int encode_dimensions(Dimensions... dims)
{
  return (... | dims);
}

template <typename... Dimensions>
constexpr int select_dimensions(Dimensions... dims)
{
  return encode_dimensions((1 << dims)...);
}

constexpr bool dimension_selected(int encoded_dims, int dim) { return encoded_dims & (1 << dim); }

/*!
 * \brief Call Blueprint mesh verify functions and convert the output to SLIC_ERROR
 *        if the verify method failed.
 *
 * \param obj The node that contains the object being checked.
 * \param protocol The name of the item to check in the mesh. If the string is empty,
 *                 \a obj node is treated as a mesh and it all gets checked.
 */
void verify(const conduit::Node &obj, const std::string &protocol = std::string());

}  // end namespace views
}  // end namespace bump
}  // end namespace axom

#endif
