// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_RADIXTREE_HPP_
#define AXOM_SPIN_RADIXTREE_HPP_

#include "axom/core/Array.hpp"
#include "axom/core/AnnotationMacros.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"

namespace axom
{
namespace spin
{
namespace internal
{
namespace linear_bvh
{
/*!
 * \brief RadixTree provides a binary radix tree representation that stores a
 *  list of axis-aligned bounding boxes sorted according to their Morton code.
 *
 * \note This data-structure provides an intermediate representation that serves
 *  as the building-block to construct a BVH in parallel.
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
  axom::Array<std::uint32_t> m_mcodes;
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
    m_mcodes = axom::Array<std::uint32_t>(m_size, m_size, allocID);
    m_leaf_aabbs = axom::Array<BoxType>(ArrayOptions::Uninitialized {}, m_size, m_size, allocID);
  }
};

} /* namespace linear_bvh */
} /* namespace internal */
} /* namespace spin */
} /* namespace axom */

#endif /* AXOM_RADIXTREE_HPP_ */
