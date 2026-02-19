// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
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
