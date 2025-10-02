// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"                         /* for compile time defs */
#include "axom/core/Macros.hpp"                    /* for axom macros */
#include "axom/core/execution/execution_space.hpp" /* execution_space traits */
#include "axom/core/execution/for_all.hpp"         /* for_all() traversals */
#include "axom/core/execution/synchronize.hpp"     /* synchronize() */
#include "axom/core/memory_management.hpp"         /* allocate/deallocate */

// gtest includes
#include "gtest/gtest.h"

//------------------------------------------------------------------------------
//  HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

/// Base template
template <typename ExecSpace, int NDIMS>
struct utility
{ };

//------------------------------------------------------------------------------
/// 1D version
template <typename ExecSpace>
struct utility<ExecSpace, 1>
{
  static axom::IndexType numValues(axom::IndexType N) { return N; }
  static void initialize(int *array, axom::IndexType N, int value)
  {
    axom::for_all<ExecSpace>(
      N,
      AXOM_LAMBDA(axom::IndexType index) { array[index] = static_cast<int>(index + value); });
  }
  static void modify(int *array, axom::IndexType N, int value)
  {
    axom::for_all<ExecSpace>(
      N,
      AXOM_LAMBDA(axom::IndexType index) { array[index] -= static_cast<int>(index + value); });
  }
};

//------------------------------------------------------------------------------
/// 2D version
template <typename ExecSpace>
struct utility<ExecSpace, 2>
{
  static axom::IndexType numValues(axom::IndexType N) { return N * N; }
  static void initialize(int *array, axom::IndexType N, int value)
  {
    axom::StackArray<axom::IndexType, 2> shape {{N, N}};
    axom::for_all<ExecSpace>(
      shape,
      AXOM_LAMBDA(axom::IndexType i, axom::IndexType j) {
        const auto index = j * N + i;
        array[index] = static_cast<int>(index + value);
      });
  }
  static void modify(int *array, axom::IndexType N, int value)
  {
    axom::StackArray<axom::IndexType, 2> shape {{N, N}};
    axom::for_all<ExecSpace>(
      shape,
      AXOM_LAMBDA(axom::IndexType i, axom::IndexType j) {
        const auto index = j * N + i;
        array[index] -= static_cast<int>(index + value);
      });
  }
};

//------------------------------------------------------------------------------
/// 3D version
template <typename ExecSpace>
struct utility<ExecSpace, 3>
{
  static axom::IndexType numValues(axom::IndexType N) { return N * N * N; }
  static void initialize(int *array, axom::IndexType N, int value)
  {
    axom::StackArray<axom::IndexType, 3> shape {{N, N, N}};
    axom::for_all<ExecSpace>(
      shape,
      AXOM_LAMBDA(axom::IndexType i, axom::IndexType j, axom::IndexType k) {
        const auto index = (k * N * N) + (j * N) + i;
        array[index] = static_cast<int>(index + value);
      });
  }
  static void modify(int *array, axom::IndexType N, int value)
  {
    axom::StackArray<axom::IndexType, 3> shape {{N, N, N}};
    axom::for_all<ExecSpace>(
      shape,
      AXOM_LAMBDA(axom::IndexType i, axom::IndexType j, axom::IndexType k) {
        const auto index = (k * N * N) + (j * N) + i;
        array[index] -= static_cast<int>(index + value);
      });
  }
};

//------------------------------------------------------------------------------
template <typename ExecSpace, int NDIMS>
void check_for_all(axom::IndexType N)
{
  using utils = utility<ExecSpace, NDIMS>;

  EXPECT_TRUE(axom::execution_space<ExecSpace>::valid());
  std::cout << "checking axom::for_all() with [" << axom::execution_space<ExecSpace>::name() << "]\n";

  // STEP 1: define a constant
  constexpr int VALUE = 42;

  // STEP 2: set allocators for the execution spaces
  const int hostID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int allocID = axom::execution_space<ExecSpace>::allocatorID();

  // STEP 3: allocate buffer
  const auto arraySize = utils::numValues(N);
  int *a = axom::allocate<int>(arraySize, allocID);

  // STEP 4: initialize to (index + VALUE)
  utils::initialize(a, N, VALUE);

  if(axom::execution_space<ExecSpace>::async())
  {
    axom::synchronize<ExecSpace>();
  }

  // STEP 5: check array
  int *a_host = axom::allocate<int>(arraySize, hostID);
  axom::copy(a_host, a, arraySize * sizeof(int));

  for(int i = 0; i < arraySize; ++i)
  {
    EXPECT_EQ(a_host[i], i + VALUE);
  }

  // STEP 6: Subtract (index + VALUE) from all entries resulting in zero
  utils::modify(a, N, VALUE);

  if(axom::execution_space<ExecSpace>::async())
  {
    axom::synchronize<ExecSpace>();
  }

  // STEP 7: check result
  axom::copy(a_host, a, arraySize * sizeof(int));

  for(int i = 0; i < arraySize; ++i)
  {
    EXPECT_EQ(a_host[i], 0);
  }

  // STEP 8: cleanup
  axom::deallocate(a);
  axom::deallocate(a_host);
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------
TEST(core_execution_for_all, seq_exec)
{
  check_for_all<axom::SEQ_EXEC, 1>(256);
  check_for_all<axom::SEQ_EXEC, 2>(64);
  check_for_all<axom::SEQ_EXEC, 3>(32);
}

//------------------------------------------------------------------------------

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP)

TEST(core_execution_for_all, omp_exec)
{
  check_for_all<axom::OMP_EXEC, 1>(256);
  check_for_all<axom::OMP_EXEC, 2>(64);
  check_for_all<axom::OMP_EXEC, 3>(32);
}

#endif

//------------------------------------------------------------------------------

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)

TEST(core_execution_for_all, cuda_exec)
{
  constexpr int BLOCK_SIZE = 256;
  check_for_all<axom::CUDA_EXEC<BLOCK_SIZE>, 1>(256);
  check_for_all<axom::CUDA_EXEC<BLOCK_SIZE>, 2>(64);
  check_for_all<axom::CUDA_EXEC<BLOCK_SIZE>, 3>(32);
}

//------------------------------------------------------------------------------
TEST(core_execution_for_all, cuda_exec_async)
{
  constexpr int BLOCK_SIZE = 256;
  check_for_all<axom::CUDA_EXEC<BLOCK_SIZE, axom::ASYNC>, 1>(256);
  check_for_all<axom::CUDA_EXEC<BLOCK_SIZE, axom::ASYNC>, 2>(64);
  check_for_all<axom::CUDA_EXEC<BLOCK_SIZE, axom::ASYNC>, 3>(32);
}

#endif

//------------------------------------------------------------------------------

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(AXOM_USE_UMPIRE)

TEST(core_execution_for_all, hip_exec)
{
  constexpr int BLOCK_SIZE = 256;
  check_for_all<axom::HIP_EXEC<BLOCK_SIZE>, 1>(256);
  check_for_all<axom::HIP_EXEC<BLOCK_SIZE>, 2>(64);
  check_for_all<axom::HIP_EXEC<BLOCK_SIZE>, 3>(32);
}

//------------------------------------------------------------------------------
TEST(core_execution_for_all, hip_exec_async)
{
  constexpr int BLOCK_SIZE = 256;
  check_for_all<axom::HIP_EXEC<BLOCK_SIZE, axom::ASYNC>, 1>(256);
  check_for_all<axom::HIP_EXEC<BLOCK_SIZE, axom::ASYNC>, 2>(64);
  check_for_all<axom::HIP_EXEC<BLOCK_SIZE, axom::ASYNC>, 3>(32);
}

#endif
