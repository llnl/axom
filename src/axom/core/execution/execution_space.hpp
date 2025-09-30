// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_EXECUTIONSPACE_HPP_
#define AXOM_EXECUTIONSPACE_HPP_

#include "axom/config.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/core/execution/runtime_policy.hpp"

/*!
 * \file execution_space.hpp
 *
 * \brief Defines the list of available execution spaces for axom.
 *
 * The list of defined execution spaces are the following:
 *
 *  * <b>SEQ_EXEC<b> <br />
 *
 *    Indicates sequential execution on the CPU. Always defined.
 *
 *    When using this execution space, the data must reside on CPU/host memory.
 *
 *  * <b>OMP_EXEC<b> <br />
 *
 *    Indicates parallel execution with OpenMP on the CPU.
 *
 *    Defined when AXOM_USE_OPENMP and AXOM_USE_RAJA are defined. In addition,
 *    using this execution space requires linking to a RAJA that is configured
 *    with OpenMp, i.e., RAJA_ENABLE_OPENMP must be defined in the generated
 *    RAJA/config.hpp.
 *
 *    The default memory allocator when using this execution space is HOST.
 *
 *    When using this execution space, the data must reside on CPU/host memory.
 *
 *  * <b>CUDA_EXEC<BLOCKSIZE></b> <br />
 *
 *    Indicates parallel execution with CUDA on the GPU.
 *
 *    Defined when AXOM_USE_CUDA and AXOM_USE_RAJA are defined. In  addition,
 *    using this execution space requires linking to a RAJA that is configured
 *    with CUDA, i.e., RAJA_ENABLE_CUDA must be defined in the generated
 *    RAJA/config.hpp.
 *
 *    The CUDA_EXEC execution space strictly requires Umpire for memory
 *    management. The data must reside either on device memory which, is only
 *    accessible on the GPU, or in Unified memory which, is accessible on both
 *    the CPU or GPU. Note, when using Unified memory, the hardware transfers
 *    the data accordingly as needed. Consequently, using unified memory must
 *    be handled with extreme care, in order to avoid the latencies associated
 *    with data transfers between the CPU and GPU.
 *
 *    The default memory allocator when using this execution space is UNIFIED.
 *
 * The execution spaces bind the corresponding RAJA execution policies and
 * default memory space.
 *
 * @see runtime_policy.hpp
 */

namespace axom
{
/*!
 * \brief The execution_space is a traits class that binds the execution
 *  space to a corresponding RAJA execution policies and default memory
 *  allocator.
 *
 * \tparam ExecSpace the execution space
 *
 * \note This class is specialized for each execution space.
 *
 */
template <typename ExecSpace>
struct execution_space
{
  using loop_policy = void;

  using reduce_policy = void;
  using atomic_policy = void;
  using sync_policy = void;

  static constexpr MemorySpace memory_space = MemorySpace::Dynamic;

  static constexpr bool async() noexcept { return false; }
  static constexpr bool valid() noexcept { return false; }
  static constexpr bool onDevice() noexcept { return false; }
  static constexpr char* name() noexcept { return (char*)"[UNDEFINED]"; }
  static int allocatorID() noexcept { return axom::INVALID_ALLOCATOR_ID; }
  static constexpr runtime_policy::Policy runtimePolicy() noexcept
  {
    return runtime_policy::Policy::seq;
  }
  //!@brief Returns whether @c ExecSpace can use the given @c MemorySpace.
  static bool usesMemorySpace(axom::MemorySpace m) noexcept { return m == memory_space; }
  //!@brief Returns whether @c ExecSpace can use the given allocator id.
  static bool usesAllocId(int allocId) noexcept
  {
    return allocId == axom::INVALID_ALLOCATOR_ID
      ? false
      : usesMemorySpace(axom::detail::getAllocatorSpace(allocId));
  }
};

}  // namespace axom

// execution_space traits specialization
#include "axom/core/execution/internal/seq_exec.hpp"

#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  #include "axom/core/execution/internal/omp_exec.hpp"
#endif

#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && \
  defined(__CUDACC__)
  #include "axom/core/execution/internal/cuda_exec.hpp"
#endif

#if defined(AXOM_USE_HIP) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && \
  defined(__HIPCC__)
  #include "axom/core/execution/internal/hip_exec.hpp"
#endif

namespace axom
{

/// \brief Return default allocator id for a runtime policy.
inline int policyToDefaultAllocatorID(axom::runtime_policy::Policy policy)
{
  return policy == axom::runtime_policy::Policy::seq
    ? axom::execution_space<axom::SEQ_EXEC>::allocatorID()
    :
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
    policy == axom::runtime_policy::Policy::omp
    ? axom::execution_space<axom::OMP_EXEC>::allocatorID()
    :
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
    policy == axom::runtime_policy::Policy::cuda
    ? axom::execution_space<axom::CUDA_EXEC<256>>::allocatorID()
    :
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
    policy == axom::runtime_policy::Policy::hip
    ? axom::execution_space<axom::HIP_EXEC<256>>::allocatorID()
    :
#endif
    INVALID_ALLOCATOR_ID;
}

}  // namespace axom

#endif  // AXOM_EXECUTIONSPACE_HPP_
