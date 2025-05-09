// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_DIMENSIONS_HPP
#define AXOM_DIMENSIONS_HPP

namespace axom
{
namespace klee
{
/**
 * The dimensions that are supported for specifying operations in Klee.
 */
enum class Dimensions : int
{
  Two = 2,
  Three = 3,
};

}  // namespace klee
}  // namespace axom

#endif  //AXOM_DIMENSIONS_HPP
