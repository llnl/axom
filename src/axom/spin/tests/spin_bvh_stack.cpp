// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/spin/internal/linear_bvh/bvh_traverse.hpp"

#include "gtest/gtest.h"

namespace
{
namespace lbvh = axom::spin::internal::linear_bvh;

TEST(spin_bvh_stack, initializes_with_barrier_and_is_lifo)
{
  lbvh::BVHStack stack;
  lbvh::BVHStack::LocalStack local_stack;
  stack.setLocalStack(local_stack);

  EXPECT_EQ(stack.pop(), lbvh::BVHStack::BARRIER);

  stack.setLocalStack(local_stack);
  stack.push(17);
  stack.push(23);
  stack.push(42);

  EXPECT_EQ(stack.pop(), 42);
  EXPECT_EQ(stack.pop(), 23);
  EXPECT_EQ(stack.pop(), 17);
  EXPECT_EQ(stack.pop(), lbvh::BVHStack::BARRIER);
}

TEST(spin_bvh_stack, reset_discards_entries_and_supports_maximum_depth)
{
  lbvh::BVHStack stack;
  lbvh::BVHStack::LocalStack local_stack;
  stack.setLocalStack(local_stack);

  stack.push(1);
  stack.push(2);
  stack.setLocalStack(local_stack);
  EXPECT_EQ(stack.pop(), lbvh::BVHStack::BARRIER);

  stack.setLocalStack(local_stack);
  constexpr int max_depth = lbvh::BVHStack::STACK_SIZE - 1;
  for(int i = 0; i < max_depth; ++i)
  {
    stack.push(i);
  }

  for(int i = max_depth - 1; i >= 0; --i)
  {
    EXPECT_EQ(stack.pop(), i);
  }
  EXPECT_EQ(stack.pop(), lbvh::BVHStack::BARRIER);
}

#if defined(AXOM_USE_GPU) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && \
  (defined(__CUDACC__) || defined(__HIPCC__))

template <typename ExecSpace>
void check_shared_stack_chunk_save_and_restore()
{
  constexpr int block_size = axom::execution_space<ExecSpace>::BlockSize;
  using SharedStack = lbvh::SharedBVHStack<block_size>;
  constexpr int chunk_size = SharedStack::CHUNK_SIZE;
  constexpr int shared_capacity = SharedStack::SHMEM_SIZE_PER_THREAD;
  constexpr int num_chunk_spills = 8;

  // After the first shared_capacity pushes, push() saves a chunk to g_stack.
  // Every additional chunk_size pushes saves another chunk. This depth causes
  // exactly num_chunk_spills saves from shared memory to g_stack. It leaves
  // five values in shared memory, so popping first drains shared memory; each
  // subsequent pop restores a saved chunk from g_stack to shared memory.
  constexpr int stack_depth = shared_capacity + 1 + (num_chunk_spills - 1) * chunk_size;
  static_assert(shared_capacity == 2 * chunk_size, "SharedBVHStack capacity changed");
  static_assert(stack_depth == 37, "Update the spill/restore test depth");

  constexpr int num_threads = 2 * block_size;
  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int device_allocator = axom::execution_space<ExecSpace>::allocatorID();
  axom::Array<int> device_results(num_threads, num_threads, device_allocator);

  auto results = device_results.view();
  axom::for_all<ExecSpace>(
    num_threads,
    AXOM_LAMBDA(axom::IndexType thread_idx) {
      SharedStack stack;
      typename SharedStack::LocalStack local_stack;
      stack.setLocalStack(local_stack);

      bool passed = stack.pop() == SharedStack::BARRIER;
      const int seed = static_cast<int>(thread_idx) * stack_depth;
      for(int i = 0; i < stack_depth; ++i)
      {
        // At full shared capacity, push() saves its oldest Chunk to g_stack.
        stack.push(seed + i);
      }
      for(int i = stack_depth - 1; i >= 0; --i)
      {
        // Once shared storage is empty, pop() restores a Chunk from g_stack.
        passed = passed && stack.pop() == seed + i;
      }
      passed = passed && stack.pop() == SharedStack::BARRIER;
      results[thread_idx] = passed ? 1 : 0;
    });

  axom::Array<int> host_results(device_results, host_allocator);
  for(axom::IndexType i = 0; i < host_results.size(); ++i)
  {
    EXPECT_EQ(host_results[i], 1) << "thread " << i;
  }
}

TEST(spin_bvh_stack, shared_stack_saves_and_restores_global_chunks)
{
  #if defined(__CUDACC__)
  check_shared_stack_chunk_save_and_restore<axom::CUDA_EXEC<256>>();
  #elif defined(__HIPCC__)
  check_shared_stack_chunk_save_and_restore<axom::HIP_EXEC<256>>();
  #endif
}

#endif

}  // namespace

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;
  return RUN_ALL_TESTS();
}
