// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/config.hpp"
#include "axom/core/memory_management.hpp"

// RAJA includes
#ifdef AXOM_USE_RAJA
  #include "RAJA/RAJA.hpp"
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

  static constexpr MemorySpace memory_space = MemorySpace::Host;

  static constexpr bool async() noexcept { return false; }
  static constexpr bool valid() noexcept { return true; }
  static constexpr bool onDevice() noexcept { return false; }
  static constexpr const char* name() noexcept { return "[SEQ_EXEC]"; }
  static int allocatorID() noexcept
  {
    return axom::getAllocatorIDFromMemorySpace(memory_space);
  }
  static constexpr runtime_policy::Policy runtimePolicy() noexcept
  {
    return runtime_policy::Policy::seq;
  }
  static bool usesMemorySpace(axom::MemorySpace m) noexcept
  {
#if defined(AXOM_USE_UMPIRE)
    return m == MemorySpace::Dynamic || m == MemorySpace::Malloc || m == MemorySpace::Host ||
      (m == MemorySpace::Unified && axom::isMemorySpaceAvailable(MemorySpace::Unified));
#else
    return m == MemorySpace::Dynamic || m == MemorySpace::Malloc || m == MemorySpace::Host;
#endif
  }
  static bool usesAllocId(int allocId) noexcept
  {
    return allocId == axom::INVALID_ALLOCATOR_ID
      ? false
      : usesMemorySpace(axom::detail::getAllocatorSpace(allocId));
  }
};

}  // namespace axom
