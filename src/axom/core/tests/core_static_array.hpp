// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"  // for compile time definitions

#include "axom/core/StaticArray.hpp"

// for gtest macros
#include "gtest/gtest.h"

//------------------------------------------------------------------------------
//  HELPER METHOD
//------------------------------------------------------------------------------
namespace
{
template <typename ExecSpace>
void check_static_array_policy()
{
  const int MAX_SIZE = 10;
  using StaticArrayType = axom::StaticArray<int, MAX_SIZE>;

  // Get ids of necessary allocators
  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int kernel_allocator = axom::execution_space<ExecSpace>::allocatorID();

  axom::Array<StaticArrayType> s_arrays_device(MAX_SIZE, MAX_SIZE, kernel_allocator);
  auto s_arrays_view = s_arrays_device.view();

  axom::Array<StaticArrayType> sizes_device(1, 1, kernel_allocator);
  auto sizes_view = sizes_device.view();

  axom::for_all<ExecSpace>(
    MAX_SIZE,
    AXOM_LAMBDA(int i) {
      // Sanity check - function is callable on device
      s_arrays_view[i].clear();

      for(int idx = 0; idx <= i; idx++)
      {
        s_arrays_view[i].push_back(idx);
      }

      sizes_view[0][i] = s_arrays_view[i].size();
    });

  if(axom::execution_space<ExecSpace>::async())
  {
    axom::synchronize<ExecSpace>();
  }

  // Copy static arrays and their sizes back to host
  axom::Array<StaticArrayType> s_arrays_host =
    axom::Array<StaticArrayType>(s_arrays_device, host_allocator);
  axom::Array<StaticArrayType> sizes_host =
    axom::Array<StaticArrayType>(sizes_device, host_allocator);

  // Verify values
  for(int i = 0; i < MAX_SIZE; i++)
  {
    EXPECT_EQ(sizes_host[0][i], i + 1);
    for(int j = 0; j <= i; j++)
    {
      EXPECT_EQ(s_arrays_host[i][j], j);
    }
  }
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
TEST(core_static_array, check_static_array_seq) { check_static_array_policy<axom::SEQ_EXEC>(); }

#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
TEST(core_static_array, check_static_array_omp) { check_static_array_policy<axom::OMP_EXEC>(); }
#endif

#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
TEST(core_static_array, check_static_array_cuda)
{
  check_static_array_policy<axom::CUDA_EXEC<256>>();
  check_static_array_policy<axom::CUDA_EXEC<256, axom::ASYNC>>();
}
#endif

#if defined(AXOM_USE_HIP) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
TEST(core_static_array, check_static_array_hip)
{
  check_static_array_policy<axom::HIP_EXEC<256>>();
  check_static_array_policy<axom::HIP_EXEC<256, axom::ASYNC>>();
}
#endif
