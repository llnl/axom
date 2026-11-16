// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/core/Array.hpp"
#include "axom/core/AnnotationMacros.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"

namespace axom::spin::internal::linear_bvh
{
/*!
 * \brief A binary radix tree of axis-aligned boxes sorted by spatial key.
 *
 * \note LinearBVH uses this intermediate tree to build its traversal nodes in parallel.
 */
template <typename FloatType, int NDIMS>
struct RadixTree
{
  using BoxType = primal::BoundingBox<FloatType, NDIMS>;

  std::int32_t m_size;
  std::int32_t m_inner_size;

  axom::Array<std::int32_t> m_left_children;
  axom::Array<std::int32_t> m_right_children;
  axom::Array<std::int32_t> m_parents;
  axom::Array<BoxType> m_inner_aabbs;

  axom::Array<std::int32_t> m_leafs;
  axom::Array<std::uint32_t> m_keys;
  axom::Array<BoxType> m_leaf_aabbs;

  void allocate(std::int32_t size, int allocID)
  {
    AXOM_ANNOTATE_SCOPE("RadixTree::allocate");

    m_size = size;
    m_inner_size = m_size - 1;
    std::int32_t parent_size = m_size + m_inner_size;

    m_left_children = axom::Array<std::int32_t>(m_inner_size, m_inner_size, allocID);
    m_right_children = axom::Array<std::int32_t>(m_inner_size, m_inner_size, allocID);
    m_parents = axom::Array<std::int32_t>(parent_size, parent_size, allocID);
    m_inner_aabbs =
      axom::Array<BoxType>(ArrayOptions::Uninitialized {}, m_inner_size, m_inner_size, allocID);

    m_leafs = axom::Array<std::int32_t>(m_size, m_size, allocID);
    m_keys = axom::Array<std::uint32_t>(m_size, m_size, allocID);
    m_leaf_aabbs = axom::Array<BoxType>(ArrayOptions::Uninitialized {}, m_size, m_size, allocID);
  }
};

}  // namespace axom::spin::internal::linear_bvh
