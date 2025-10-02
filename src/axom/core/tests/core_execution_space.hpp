// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"  // for compile time definitions

// spin includes
#include "axom/core/execution/execution_space.hpp"

#ifdef AXOM_USE_UMPIRE
  #include "umpire/Umpire.hpp"  // for Umpire
#endif

#ifdef AXOM_USE_RAJA
  #include "RAJA/RAJA.hpp"  // for RAJA
#endif

// for gtest macros
#include "gtest/gtest.h"

// C/C++ includes
#include <cstring>      // for strlen(),strcmp()
#include <type_traits>  // for std::is_same

//------------------------------------------------------------------------------
//  HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
template <typename ExecSpace>
void check_valid()
{
  std::cout << "checking execution space:" << axom::execution_space<ExecSpace>::name() << std::endl;

  EXPECT_TRUE(axom::execution_space<ExecSpace>::valid());
  EXPECT_TRUE(strlen(axom::execution_space<ExecSpace>::name()) > 0);
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
void check_invalid()
{
  std::cout << "checking execution space:" << axom::execution_space<ExecSpace>::name() << std::endl;

  EXPECT_FALSE(axom::execution_space<ExecSpace>::valid());

  EXPECT_EQ(axom::execution_space<ExecSpace>::allocatorID(), axom::INVALID_ALLOCATOR_ID);

  EXPECT_EQ(strcmp(axom::execution_space<ExecSpace>::name(), "[UNDEFINED]"), 0);
}

//------------------------------------------------------------------------------
template <typename ExecSpace, typename LoopPolicy, typename ReducePolicy, typename AtomicPolicy, typename SyncPolicy>
void check_execution_mappings(int expectedAllocatorID, bool is_async, bool is_onDevice)
{
  std::cout << "checking execution space: " << axom::execution_space<ExecSpace>::name() << std::endl;

  using loop_pol = typename axom::execution_space<ExecSpace>::loop_policy;
  using reduce_pol = typename axom::execution_space<ExecSpace>::reduce_policy;
  using atomic_pol = typename axom::execution_space<ExecSpace>::atomic_policy;
  using sync_pol = typename axom::execution_space<ExecSpace>::sync_policy;

  bool valid_loop_policy = std::is_same<loop_pol, LoopPolicy>::value;
  bool valid_reduce_policy = std::is_same<reduce_pol, ReducePolicy>::value;
  bool valid_atomic_policy = std::is_same<atomic_pol, AtomicPolicy>::value;
  bool valid_sync_policy = std::is_same<sync_pol, SyncPolicy>::value;

  EXPECT_TRUE(valid_loop_policy);
  EXPECT_TRUE(valid_reduce_policy);
  EXPECT_TRUE(valid_atomic_policy);
  EXPECT_TRUE(valid_sync_policy);

  const bool async = axom::execution_space<ExecSpace>::async();
  EXPECT_EQ(async, is_async);

  const bool onDevice = axom::execution_space<ExecSpace>::onDevice();
  EXPECT_EQ(onDevice, is_onDevice);

  int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
  EXPECT_EQ(expectedAllocatorID, allocatorID);
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
TEST(core_execution_space, check_valid)
{
  check_valid<axom::SEQ_EXEC>();

#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  check_valid<axom::OMP_EXEC>();
#endif

#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  check_valid<axom::CUDA_EXEC<256>>();
  check_valid<axom::CUDA_EXEC<256, axom::ASYNC>>();
#endif

#if defined(AXOM_USE_HIP) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  check_valid<axom::HIP_EXEC<256>>();
  check_valid<axom::HIP_EXEC<256, axom::ASYNC>>();
#endif
}

//------------------------------------------------------------------------------
TEST(core_execution_space, check_invalid)
{
  struct InvalidSpace
  { };
  check_invalid<InvalidSpace>();
}

//==============================================================================
// The following tests require RAJA and UMPIRE
//==============================================================================
#if defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_RAJA)

TEST(core_execution_space, check_seq_exec)
{
  check_valid<axom::SEQ_EXEC>();

  constexpr bool IS_ASYNC = false;
  constexpr bool ON_DEVICE = false;

  int allocator_id = axom::getUmpireResourceAllocatorID(umpire::resource::Host);
  check_execution_mappings<axom::SEQ_EXEC,
  #if RAJA_VERSION_MAJOR > 2022
                           RAJA::seq_exec,
                           RAJA::seq_reduce,
                           RAJA::seq_atomic,
  #else
                           RAJA::loop_exec,
                           RAJA::loop_reduce,
                           RAJA::loop_atomic,
  #endif
                           void>(allocator_id, IS_ASYNC, ON_DEVICE);
}

  //------------------------------------------------------------------------------
  #if defined(AXOM_USE_OPENMP)
TEST(core_execution_space, check_omp_exec)
{
  check_valid<axom::OMP_EXEC>();

  constexpr bool IS_ASYNC = false;
  constexpr bool ON_DEVICE = false;

  int allocator_id = axom::getUmpireResourceAllocatorID(umpire::resource::Host);
  check_execution_mappings<axom::OMP_EXEC,
                           RAJA::omp_parallel_for_exec,
                           RAJA::omp_reduce,
                           RAJA::omp_atomic,
                           RAJA::omp_synchronize>(allocator_id, IS_ASYNC, ON_DEVICE);
}
  #endif  // defined(AXOM_USE_OPENMP)

  //------------------------------------------------------------------------------
  #if defined(AXOM_USE_CUDA)

TEST(core_execution_space, check_cuda_exec)
{
  constexpr int BLOCK_SIZE = 256;

  check_valid<axom::CUDA_EXEC<BLOCK_SIZE>>();

  constexpr bool IS_ASYNC = false;
  constexpr bool ON_DEVICE = true;

  int allocator_id = axom::getUmpireResourceAllocatorID(umpire::resource::Device);
  check_execution_mappings<axom::CUDA_EXEC<BLOCK_SIZE>,
                           RAJA::cuda_exec<BLOCK_SIZE>,
                           RAJA::cuda_reduce,
                           RAJA::cuda_atomic,
                           RAJA::cuda_synchronize>(allocator_id, IS_ASYNC, ON_DEVICE);
}

//------------------------------------------------------------------------------
TEST(core_execution_space, check_cuda_exec_async)
{
  constexpr int BLOCK_SIZE = 256;

  check_valid<axom::CUDA_EXEC<BLOCK_SIZE, axom::ASYNC>>();

  constexpr bool IS_ASYNC = true;
  constexpr bool ON_DEVICE = true;

  int allocator_id = axom::getUmpireResourceAllocatorID(umpire::resource::Device);
  check_execution_mappings<axom::CUDA_EXEC<BLOCK_SIZE, axom::ASYNC>,
                           RAJA::cuda_exec_async<BLOCK_SIZE>,
                           RAJA::cuda_reduce,
                           RAJA::cuda_atomic,
                           RAJA::cuda_synchronize>(allocator_id, IS_ASYNC, ON_DEVICE);
}
//------------------------------------------------------------------------------
void build(axom::Array<axom::IndexType> &values, axom::Array<axom::IndexType> &ids, int allocatorID)
{
  const std::vector<axom::IndexType> data {{0, 1, 2, 3}};
  const axom::IndexType n = static_cast<axom::IndexType>(data.size());

  values = axom::Array<axom::IndexType>(n, n, allocatorID);
  ids = axom::Array<axom::IndexType>(n, n, allocatorID);

  axom::copy(values.data(), data.data(), sizeof(axom::IndexType) * n);
  axom::copy(ids.data(), data.data(), sizeof(axom::IndexType) * n);
}

template <typename ExecSpace>
void test_check_cuda_array()
{
  int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

  axom::Array<axom::IndexType> values, ids;
  build(values, ids, allocatorID);

  EXPECT_EQ(values.size(), ids.size());
  EXPECT_TRUE(axom::isDeviceAllocator(values.getAllocatorID()));
  EXPECT_TRUE(axom::isDeviceAllocator(ids.getAllocatorID()));
}

TEST(core_execution_space, check_cuda_array)
{
  constexpr int BLOCK_SIZE = 256;
  using ExecSpace = axom::CUDA_EXEC<BLOCK_SIZE>;
  test_check_cuda_array<ExecSpace>();
}
  #endif  // defined(AXOM_USE_CUDA)

  //------------------------------------------------------------------------------
  #if defined(AXOM_USE_HIP)

TEST(core_execution_space, check_hip_exec)
{
  constexpr int BLOCK_SIZE = 256;

  check_valid<axom::HIP_EXEC<BLOCK_SIZE>>();

  constexpr bool IS_ASYNC = false;
  constexpr bool ON_DEVICE = true;

  int allocator_id = axom::getUmpireResourceAllocatorID(umpire::resource::Device);
  check_execution_mappings<axom::HIP_EXEC<BLOCK_SIZE>,
                           RAJA::hip_exec<BLOCK_SIZE>,
                           RAJA::hip_reduce,
                           RAJA::hip_atomic,
                           RAJA::hip_synchronize>(allocator_id, IS_ASYNC, ON_DEVICE);
}

//------------------------------------------------------------------------------
TEST(core_execution_space, check_hip_exec_async)
{
  constexpr int BLOCK_SIZE = 256;

  check_valid<axom::HIP_EXEC<BLOCK_SIZE, axom::ASYNC>>();

  constexpr bool IS_ASYNC = true;
  constexpr bool ON_DEVICE = true;

  int allocator_id = axom::getUmpireResourceAllocatorID(umpire::resource::Device);
  check_execution_mappings<axom::HIP_EXEC<BLOCK_SIZE, axom::ASYNC>,
                           RAJA::hip_exec_async<BLOCK_SIZE>,
                           RAJA::hip_reduce,
                           RAJA::hip_atomic,
                           RAJA::hip_synchronize>(allocator_id, IS_ASYNC, ON_DEVICE);
}
  #endif  // defined(AXOM_USE_HIP)

#endif  // defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_RAJA)
