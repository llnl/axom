// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! \file core_containers.cpp
 *  \brief This example code is a demonstration of the Axom Core containers.
 */

/* This example code contains snippets used in the Core Sphinx documentation.
  * They begin and end with comments such as
  *
  * timer_start
  * timer_end
  *
  * each prepended with an underscore.
  */

// Axom includes
#include "axom/core/StackArray.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/atomics.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/execution/reductions.hpp"
#include "axom/core/execution/synchronize.hpp"
#include "axom/core/memory_management.hpp"

// C/C++ includes
#include <iostream>
#include <vector>

constexpr int N = 1000;

void demoMemoryManageBasic()
{
  // _membasic_start
  int *dynamic_memory_array;
  int *dyn_array_dst;
  int len = 20;

  //Allocation looks similar to use of malloc() in C -- just template
  //return type instead of casting.
  dynamic_memory_array = axom::allocate<int>(len);

  for(int i = 0; i < len; i++)
  {
    dynamic_memory_array[i] = i;
  }

  //Print array values after initialization
  for(int i = 0; i < len; i++)
  {
    std::cout << i << " Current value: " << dynamic_memory_array[i] << std::endl;
  }

  dyn_array_dst = axom::allocate<int>(len);

  //Now, a copy operation. It's used exactly like memcpy --
  //destination, source, number of bytes.
  axom::copy(dyn_array_dst, dynamic_memory_array, sizeof(int) * len);

  //Print array values and compare to copy
  for(int i = 0; i < len; i++)
  {
    std::cout << i << " Current value: " << dyn_array_dst[i] << std::endl;
    std::cout << "Matches old value? " << std::boolalpha
              << (dynamic_memory_array[i] == dyn_array_dst[i]) << std::endl;
  }

  //Deallocate is exactly like free. Of course, you cannot access the
  //now-deallocated memory after this:
  axom::deallocate(dyn_array_dst);

  //Reallocate is like realloc -- copies existing contents into new
  //memory allocation.
  //Slight deviation from realloc() in that second arg is item count,
  //rather than bytes.
  dynamic_memory_array = axom::reallocate(dynamic_memory_array, len * 2);
  for(int i = 20; i < len * 2; i++)
  {
    dynamic_memory_array[i] = i;
  }

  for(int i = 0; i < len * 2; i++)
  {
    std::cout << i << " Current value: " << dynamic_memory_array[i] << std::endl;
  }

  axom::deallocate(dynamic_memory_array);
  // _membasic_end
}

void demoAxomExecution()
{
  // _exebasic_start
  //This part of the code works regardless of Umpire's presence, allowing for generic
  //use of axom::allocate in C++ code.
  int *A = axom::allocate<int>(N);
  int *B = axom::allocate<int>(N);
  int *C = axom::allocate<int>(N);

  for(int i = 0; i < N; i++)
  {
    A[i] = i * 5;
    B[i] = i * 2;
    C[i] = 0;
  }

  //Axom provides an API for the most basic usage of RAJA, the for_all loop.
  axom::for_all<axom::SEQ_EXEC>(0, N, AXOM_LAMBDA(axom::IndexType i) { C[i] = A[i] + B[i]; });

  std::cout << "Sums: " << std::endl;
  for(int i = 0; i < N; i++)
  {
    std::cout << C[i] << " ";
    C[i] = 0;
  }

  axom::deallocate(A);
  axom::deallocate(B);
  axom::deallocate(C);

// _exebasic_end

//Now, let's say we want to try out use of CUDA or HIP. We just change that execution space.
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_GPU) && \
  defined(AXOM_GPUCC)
  // _deviceexebasic_start
  //This example requires Umpire to be in use, and Unified memory available.
  const int allocator_id =
    axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Unified);
  A = axom::allocate<int>(N, allocator_id);
  B = axom::allocate<int>(N, allocator_id);
  C = axom::allocate<int>(N, allocator_id);

  for(int i = 0; i < N; i++)
  {
    A[i] = i * 5;
    B[i] = i * 2;
    C[i] = 0;
  }

  #if defined(__CUDACC__)
  using ExecSpace = axom::CUDA_EXEC<256>;
  #elif defined(__HIPCC__)
  using ExecSpace = axom::HIP_EXEC<256>;
  #else
  using ExecSpace = axom::SEQ_EXEC;
  #endif

  axom::for_all<ExecSpace>(0, N, AXOM_LAMBDA(axom::IndexType i) { C[i] = A[i] + B[i]; });

  std::cout << "\nSums (" << axom::execution_space<ExecSpace>::name() << ") :" << std::endl;
  for(int i = 0; i < N; i++)
  {
    std::cout << C[i] << " ";
  }
  std::cout << std::endl;

  axom::deallocate(A);
  axom::deallocate(B);
  axom::deallocate(C);
  // _deviceexebasic_end

  //_gpu_reduce_start
  axom::ReduceSum<ExecSpace, axom::IndexType> totalSum(0);

  // Sum integers [0,99]
  axom::for_all<ExecSpace>(100, AXOM_LAMBDA(axom::IndexType i) { totalSum += i; });

  std::cout << "\nTotal Reduction Sum (" << axom::execution_space<ExecSpace>::name()
            << ") :" << totalSum.get() << std::endl;
  //_gpu_reduce_end

  //_gpu_atomic_start

  int *sum = axom::allocate<int>(1, allocator_id);
  *sum = 0;

  // Increment sum 100 times
  axom::for_all<ExecSpace>(100, AXOM_LAMBDA(axom::IndexType) { axom::atomicAdd<ExecSpace>(sum, 1); });

  std::cout << "\nTotal Atomic Sum (" << axom::execution_space<ExecSpace>::name() << ") :" << sum[0]
            << std::endl;

  axom::deallocate(sum);
  //_gpu_atomic_end

#endif
}

int main(int AXOM_UNUSED_PARAM(argc), char **AXOM_UNUSED_PARAM(argv))
{
  demoMemoryManageBasic();
  demoAxomExecution();
  return 0;
}
