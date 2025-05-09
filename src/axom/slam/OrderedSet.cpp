// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file OrderedSet.cpp
 */

#include "OrderedSet.hpp"

namespace axom
{
namespace slam
{
namespace policies
{
const NullSet<> NoSubset::s_nullSet;
NullSet<> VirtualParentSubset::s_nullSet;

}  // namespace policies
}  // namespace slam
}  // namespace axom
