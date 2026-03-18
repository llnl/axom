// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_CONDUIT_ARRAY_VIEW_HPP_
#define AXOM_BUMP_CONDUIT_ARRAY_VIEW_HPP_

#include "axom/bump/utilities/conduit_traits.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/slic.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace bump
{
namespace utilities
{
namespace detail
{

/*!
 * \brief Make a 1D axom::ArrayView from a typed Conduit node.
 *
 * \note This helper honors Conduit offset and stride metadata, so it supports
 *       regular interleaved layouts in addition to dense arrays.
 */
template <typename T>
inline axom::ArrayView<T> make_conduit_array_view(conduit::Node &n)
{
  SLIC_ASSERT_MSG(cpp2conduit<T>::id == n.dtype().id(),
                  "Cannot create ArrayView with a type that does not match the Conduit node.");

  const auto stride_bytes = n.dtype().stride();
  const auto element_bytes = n.dtype().element_bytes();
  SLIC_ERROR_IF(element_bytes != static_cast<conduit::index_t>(sizeof(T)),
                "Conduit element size does not match the selected node type.");
  SLIC_ERROR_IF(stride_bytes % static_cast<conduit::index_t>(sizeof(T)) != 0,
                "Conduit stride is not compatible with the selected node type.");

  auto *data = static_cast<T *>(n.element_ptr(0));
  const auto stride = stride_bytes / static_cast<conduit::index_t>(sizeof(T));
  return axom::ArrayView<T>(
    data,
    axom::StackArray<axom::IndexType, 1> {{n.dtype().number_of_elements()}},
    axom::StackArray<axom::IndexType, 1> {{stride}});
}

template <typename T>
inline axom::ArrayView<T> make_conduit_array_view(const conduit::Node &n)
{
  SLIC_ASSERT_MSG(cpp2conduit<T>::id == n.dtype().id(),
                  "Cannot create ArrayView with a type that does not match the Conduit node.");

  const auto stride_bytes = n.dtype().stride();
  const auto element_bytes = n.dtype().element_bytes();
  SLIC_ERROR_IF(element_bytes != static_cast<conduit::index_t>(sizeof(T)),
                "Conduit element size does not match the selected node type.");
  SLIC_ERROR_IF(stride_bytes % static_cast<conduit::index_t>(sizeof(T)) != 0,
                "Conduit stride is not compatible with the selected node type.");

  auto *data = const_cast<T *>(static_cast<const T *>(n.element_ptr(0)));
  const auto stride = stride_bytes / static_cast<conduit::index_t>(sizeof(T));
  return axom::ArrayView<T>(
    data,
    axom::StackArray<axom::IndexType, 1> {{n.dtype().number_of_elements()}},
    axom::StackArray<axom::IndexType, 1> {{stride}});
}

}  // namespace detail
}  // namespace utilities
}  // namespace bump
}  // namespace axom

#endif
