// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_KLEE_DIMENSIONS_HPP_
#define AXOM_KLEE_DIMENSIONS_HPP_

namespace axom
{
namespace klee
{
/// The dimensions that are supported for specifying operations in Klee.
enum class Dimensions : int
{
  Two = 2,
  Three = 3,
  Unspecified = -1
};

}  // namespace klee
}  // namespace axom

#endif  // AXOM_KLEE_DIMENSIONS_HPP_
