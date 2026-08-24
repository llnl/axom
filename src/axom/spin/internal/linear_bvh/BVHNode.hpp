// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef Axom_Spin_BVHNode_HH
#define Axom_Spin_BVHNode_HH

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
 * \brief Node structure for a 2-wide BVH tree.
 */
template <typename FloatType, int NDIMS>
struct BVH2Node
{
  using BoxType = primal::BoundingBox<FloatType, NDIMS>;

  BoxType left;
  BoxType right;
  std::int32_t left_child;
  std::int32_t right_child;
};

}  // namespace linear_bvh
}  // namespace internal
}  // namespace spin
}  // namespace axom

#endif
