// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/config.hpp"                  // for compile-time definitions
#include "axom/core/memory_management.hpp"  // alloc() / free() methods
#include "axom/core/Macros.hpp"             // for AXOM_LAMBDA

// RAJA includes
#include "RAJA/RAJA.hpp"  // for RAJA

// Google Test
#include "gtest/gtest.h"  // for google test functions

// C/C++ includes
#include <iostream>  // for std::cout

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
template <typename execution_policy>
void raja_basic_usage_test()
{
  constexpr int N = 100;
  int* a = axom::allocate<int>(N);
  int* b = axom::allocate<int>(N);
  int* c = axom::allocate<int>(N);

  // initialize
  RAJA::forall<execution_policy>(
    RAJA::RangeSegment(0, N),
    AXOM_LAMBDA(int i) {
      a[i] = b[i] = 1;
      c[i] = 0;
    });

  // add vectors
  RAJA::forall<execution_policy>(RAJA::RangeSegment(0, N), AXOM_LAMBDA(int i) { c[i] = a[i] + b[i]; });

  // check result in serial
  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(c[i], 2);
  }

  axom::deallocate(a);
  axom::deallocate(b);
  axom::deallocate(c);
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
AXOM_CUDA_TEST(raja_smoke, basic_use)
{
  std::cout << "Testing RAJA Sequential execution" << std::endl;
  raja_basic_usage_test<RAJA::seq_exec>();

#if defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)
  std::cout << "Testing RAJA OpenMP(CPU) execution" << std::endl;
  raja_basic_usage_test<RAJA::omp_parallel_for_exec>();
#endif

#if defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA) && defined(AXOM_USE_UMPIRE)
  const int prev_allocator = axom::getDefaultAllocatorID();
  const int UnifiedAllocatorID = axom::getUmpireResourceAllocatorID(umpire::resource::Unified);
  axom::setDefaultAllocator(UnifiedAllocatorID);

  std::cout << "Testing RAJA CUDA execution" << std::endl;
  constexpr int BLOCKSIZE = 256;
  raja_basic_usage_test<RAJA::cuda_exec<BLOCKSIZE>>();

  axom::setDefaultAllocator(prev_allocator);
#endif

#if defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP) && defined(AXOM_USE_UMPIRE)
  const int prev_allocator = axom::getDefaultAllocatorID();
  const int UnifiedAllocatorID = axom::getUmpireResourceAllocatorID(umpire::resource::Unified);
  axom::setDefaultAllocator(UnifiedAllocatorID);

  std::cout << "Testing RAJA HIP execution" << std::endl;
  constexpr int BLOCKSIZE = 256;
  raja_basic_usage_test<RAJA::hip_exec<BLOCKSIZE>>();

  axom::setDefaultAllocator(prev_allocator);
#endif
}

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();
  return (result);
}
