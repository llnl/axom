// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *
 * \file MemoryTesting.hpp
 *
 * \brief Header file containing utility functions for memory spaces in Axom tests
 *
 */

#ifndef AXOM_MEMORY_TESTING_HPP_
#define AXOM_MEMORY_TESTING_HPP_

#include "axom/config.hpp"
#include "axom/core/execution/execution_space.hpp"

namespace axom
{
namespace utilities
{

template <typename ExecSpace>
int globalDefaultAllocatorForExecSpace()
{
#if defined(AXOM_USE_UMPIRE)
  return axom::execution_space<ExecSpace>::onDevice()
    ? axom::execution_space<ExecSpace>::allocatorID()
    : axom::getUmpireResourceAllocatorID(umpire::resource::Host);
#else
  return axom::execution_space<ExecSpace>::allocatorID();
#endif
}

}  // namespace utilities
}  // namespace axom

#endif  // AXOM_MEMORY_TESTING_HPP_
