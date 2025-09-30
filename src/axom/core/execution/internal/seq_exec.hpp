// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SEQ_EXEC_HPP_
#define AXOM_SEQ_EXEC_HPP_

#include "axom/config.hpp"
#include "axom/core/memory_management.hpp"

// RAJA includes
#ifdef AXOM_USE_RAJA
  #include "RAJA/RAJA.hpp"
#endif

// Umpire includes
#ifdef AXOM_USE_UMPIRE
  #include "umpire/Umpire.hpp"
#endif

namespace axom
{
/*!
 * \brief Indicates sequential execution on the CPU.
 */
struct SEQ_EXEC
{ };

/*!
 * \brief execution_space traits specialization for SEQ_EXEC
 */
template <>
struct execution_space<SEQ_EXEC>
{
#ifdef AXOM_USE_RAJA
  #if RAJA_VERSION_MAJOR > 2022
  using loop_policy = RAJA::seq_exec;
  using reduce_policy = RAJA::seq_reduce;
  using atomic_policy = RAJA::seq_atomic;
  #else
  using loop_policy = RAJA::loop_exec;
  using reduce_policy = RAJA::loop_reduce;
  using atomic_policy = RAJA::loop_atomic;
  #endif
#else
  using loop_policy = void;
  using reduce_policy = void;
  using atomic_policy = void;
#endif

  using sync_policy = void;

#ifdef AXOM_USE_UMPIRE
  static constexpr MemorySpace memory_space = MemorySpace::Host;
#else
  static constexpr MemorySpace memory_space = MemorySpace::Dynamic;
#endif

  static constexpr bool async() noexcept { return false; }
  static constexpr bool valid() noexcept { return true; }
  static constexpr bool onDevice() noexcept { return false; }
  static constexpr char* name() noexcept { return (char*)"[SEQ_EXEC]"; }
  static int allocatorID() noexcept
  {
#ifdef AXOM_USE_UMPIRE
    return axom::getUmpireResourceAllocatorID(umpire::resource::Host);
#else
    return axom::getDefaultAllocatorID();
#endif
  }
  static constexpr runtime_policy::Policy runtimePolicy() noexcept
  {
    return runtime_policy::Policy::seq;
  }
  static bool usesMemorySpace(axom::MemorySpace m) noexcept
  {
    return m == MemorySpace::Dynamic || m == MemorySpace::Malloc
#ifdef AXOM_USE_UMPIRE
      || m == MemorySpace::Host || m == MemorySpace::Unified
#endif
      ;
  }
  static bool usesAllocId(int allocId) noexcept
  {
    return allocId == axom::INVALID_ALLOCATOR_ID
      ? false
      : usesMemorySpace(axom::detail::getAllocatorSpace(allocId));
  }
};

}  // namespace axom

#endif  // AXOM_SEQ_EXEC_HPP_
