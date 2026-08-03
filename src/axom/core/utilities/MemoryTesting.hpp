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

#pragma once

#include "axom/config.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/memory_management.hpp"

#include <exception>

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

inline bool runtimeMemorySpaceAvailable(axom::MemorySpace space)
{
  if(!axom::isMemorySpaceAvailable(space))
  {
    return false;
  }

#if defined(AXOM_USE_UMPIRE)
  try
  {
    switch(space)
    {
    case axom::MemorySpace::Host:
      axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Host);
      break;
    case axom::MemorySpace::Device:
  #if defined(UMPIRE_ENABLE_DEVICE)
      axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Device);
      break;
  #else
      return false;
  #endif
    case axom::MemorySpace::Unified:
  #if defined(UMPIRE_ENABLE_UM)
      axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Unified);
      break;
  #else
      return false;
  #endif
    case axom::MemorySpace::Pinned:
  #if defined(UMPIRE_ENABLE_PINNED)
      axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Pinned);
      break;
  #else
      return false;
  #endif
    case axom::MemorySpace::Constant:
  #if defined(UMPIRE_ENABLE_CONST)
      axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Constant);
      break;
  #else
      return false;
  #endif
    case axom::MemorySpace::Malloc:
    case axom::MemorySpace::Dynamic:
      break;
    }
  }
  catch(const std::exception&)
  {
    return false;
  }
#endif

  return true;
}

}  // namespace utilities
}  // namespace axom
