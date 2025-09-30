// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_OMP_EXEC_HPP_
#define AXOM_OMP_EXEC_HPP_

#include "axom/config.hpp"
#include "axom/core/memory_management.hpp"

// RAJA includes
#include "RAJA/RAJA.hpp"

#ifndef RAJA_ENABLE_OPENMP
  #error OMP_EXEC requires an OpenMP enabled RAJA
#endif

// Umpire includes
#ifdef AXOM_USE_UMPIRE
  #include "umpire/Umpire.hpp"
#endif

namespace axom
{
/*!
 * \brief Indicates parallel execution on the CPU using OpenMP.
 */
struct OMP_EXEC
{ };

/*!
 * \brief execution_space traits specialization for OMP_EXEC
 */
template <>
struct execution_space<OMP_EXEC>
{
  using loop_policy = RAJA::omp_parallel_for_exec;

  using reduce_policy = RAJA::omp_reduce;
  using atomic_policy = RAJA::omp_atomic;
  using sync_policy = RAJA::omp_synchronize;

#ifdef AXOM_USE_UMPIRE
  static constexpr MemorySpace memory_space = MemorySpace::Host;
#else
  static constexpr MemorySpace memory_space = MemorySpace::Dynamic;
#endif

  static constexpr bool async() noexcept { return false; }
  static constexpr bool valid() noexcept { return true; }
  static constexpr bool onDevice() noexcept { return false; }
  static constexpr char* name() noexcept { return (char*)"[OMP_EXEC]"; }

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
    return runtime_policy::Policy::omp;
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

#endif  // AXOM_OMP_EXEC_HPP_
