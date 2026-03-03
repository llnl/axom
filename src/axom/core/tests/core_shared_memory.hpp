// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"
#include "axom/core/memory_management.hpp"

// gtest includes
#include "gtest/gtest.h"

#include <mpi.h>
#include <string>

//------------------------------------------------------------------------------
#if defined(AXOM_USE_UMPIRE_SHARED_MEMORY)
TEST(core_shared_memory, shared_memory_allocator_min_segment_size)
{
  constexpr std::size_t min_segment_size = 32u * 1024u * 1024u;  // 32 MiB

  const int shared_alloc_id = axom::getSharedMemoryAllocatorID(min_segment_size);
  EXPECT_NE(shared_alloc_id, axom::INVALID_ALLOCATOR_ID);
  EXPECT_TRUE(axom::isSharedMemoryAllocator(shared_alloc_id));

  auto& rm = umpire::ResourceManager::getInstance();
  auto shared_alloc = rm.getAllocator(shared_alloc_id);

  const auto created_segment_size = shared_alloc.getAllocationStrategy()->getTraits().size;
  EXPECT_GE(created_segment_size, min_segment_size);

  // Subsequent calls should return the same allocator id when the allocator already exists.
  EXPECT_EQ(axom::getSharedMemoryAllocatorID(), shared_alloc_id);
  EXPECT_EQ(axom::getSharedMemoryAllocatorID(min_segment_size), shared_alloc_id);
}

TEST(core_memory_management, shared_memory_allocator_resize_aborts)
{
  // Run in a child process so the created shared-memory allocator does not affect other tests.
  constexpr std::size_t initial_min_segment_size = 0;
  constexpr std::size_t larger_min_segment_size = 64u * 1024u * 1024u;  // 64 MiB

  EXPECT_DEATH(
    {
      axom::getSharedMemoryAllocatorID(initial_min_segment_size);
      axom::getSharedMemoryAllocatorID(larger_min_segment_size);
    },
    "cannot be increased");
}

TEST(core_shared_memory, shared_memory_allocator)
{
  EXPECT_FALSE(axom::isSharedMemoryAllocator(axom::getDefaultAllocatorID()));
  EXPECT_TRUE(axom::isSharedMemoryAllocator(axom::getSharedMemoryAllocatorID()));

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const int shared_alloc_id = axom::getSharedMemoryAllocatorID();

  const int sizes[] = {1024, 4096, 1000000};
  for(int si = 0; si < 3; si++)
  {
    // Allocate shared memory
    const int N = sizes[si];
    const std::string allocation_name =
      std::string {"axom::core_shared_memory::buffer_"} + std::to_string(si);
    int* buffer = axom::allocate<int>(N, allocation_name, shared_alloc_id);

    // Populate the shared memory on rank 0
    if(rank == 0)
    {
      for(int i = 0; i < N; i++)
      {
        buffer[i] = i;
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Test the shared memory on all ranks.
    for(int i = 0; i < N; i++)
    {
      EXPECT_EQ(buffer[i], i);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    axom::deallocate(buffer);
  }
}
#endif
