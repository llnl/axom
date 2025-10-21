// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_CUDA_EXEC_HPP_
#define AXOM_CUDA_EXEC_HPP_

#include "axom/config.hpp"
#include "axom/core/memory_management.hpp"

#include "RAJA/RAJA.hpp"
#include "umpire/Umpire.hpp"

#ifndef RAJA_ENABLE_CUDA
  #error CUDA_EXEC requires a CUDA enabled RAJA
#endif

#if !defined(UMPIRE_ENABLE_CUDA) && !defined(UMPIRE_ENABLE_UM)
  #error CUDA_EXEC requires a CUDA enabled UMPIRE with UM support
#endif

namespace axom
{
enum ExecutionMode
{
  SYNCHRONOUS,
  ASYNC
};

// _cuda_exec_start
/*!
 * \brief Indicates parallel execution on the GPU with CUDA.
 *
 * \tparam BLOCK_SIZE the number of CUDA threads in a block.
 * \tparam ExecutionMode indicates synchronous or asynchronous execution.
 */
template <int BLOCK_SIZE, ExecutionMode EXEC_MODE = SYNCHRONOUS>
struct CUDA_EXEC
{ };
// _cuda_exec_end

/*!
 * \brief execution_space traits specialization for CUDA_EXEC.
 *
 * \tparam BLOCK_SIZE the number of CUDA threads to launch
 *
 */
template <int BLOCK_SIZE>
struct execution_space<CUDA_EXEC<BLOCK_SIZE, SYNCHRONOUS>>
{
  using loop_policy = RAJA::cuda_exec<BLOCK_SIZE>;

  using reduce_policy = RAJA::cuda_reduce;
  using atomic_policy = RAJA::cuda_atomic;
  using sync_policy = RAJA::cuda_synchronize;

  static constexpr MemorySpace memory_space = MemorySpace::Device;

  static constexpr bool async() noexcept { return false; }
  static constexpr bool valid() noexcept { return true; }
  static constexpr bool onDevice() noexcept { return true; }
  static constexpr char* name() noexcept { return (char*)"[CUDA_EXEC]"; }
  static int allocatorID() noexcept
  {
    return axom::getUmpireResourceAllocatorID(umpire::resource::Device);
  }
  static constexpr runtime_policy::Policy runtimePolicy() noexcept
  {
    return runtime_policy::Policy::cuda;
  }
  static bool usesMemorySpace(axom::MemorySpace m) noexcept
  {
    return m == memory_space || m == MemorySpace::Unified;
  }
  static bool usesAllocId(int allocId) noexcept
  {
    return allocId == axom::INVALID_ALLOCATOR_ID
      ? false
      : usesMemorySpace(axom::detail::getAllocatorSpace(allocId));
  }
};

/*!
 * \brief execution_space traits specialization for CUDA_EXEC.
 *
 * \tparam BLOCK_SIZE the number of CUDA threads to launch
 *
 */
template <int BLOCK_SIZE>
struct execution_space<CUDA_EXEC<BLOCK_SIZE, ASYNC>>
{
  using loop_policy = RAJA::cuda_exec_async<BLOCK_SIZE>;

  using reduce_policy = RAJA::cuda_reduce;
  using atomic_policy = RAJA::cuda_atomic;
  using sync_policy = RAJA::cuda_synchronize;

  static constexpr MemorySpace memory_space = MemorySpace::Device;

  static constexpr bool async() noexcept { return true; }
  static constexpr bool valid() noexcept { return true; }
  static constexpr bool onDevice() noexcept { return true; }
  static constexpr char* name() noexcept { return (char*)"[CUDA_EXEC] (async)"; }
  static int allocatorID() noexcept
  {
    return axom::getUmpireResourceAllocatorID(umpire::resource::Device);
  }
  static constexpr runtime_policy::Policy runtimePolicy() noexcept
  {
    return runtime_policy::Policy::cuda;
  }
  static bool usesMemorySpace(axom::MemorySpace m) noexcept
  {
    return m == memory_space || m == MemorySpace::Unified;
  }
  static bool usesAllocId(int allocId) noexcept
  {
    return allocId == axom::INVALID_ALLOCATOR_ID
      ? false
      : usesMemorySpace(axom::detail::getAllocatorSpace(allocId));
  }
};
}  // namespace axom

#endif  // AXOM_CUDA_EXEC_HPP_
