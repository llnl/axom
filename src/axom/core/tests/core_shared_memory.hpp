// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"
#include "axom/core/memory_management.hpp"

// gtest includes
#include "gtest/gtest.h"

#include <mpi.h>

//------------------------------------------------------------------------------
#if defined(AXOM_USE_UMPIRE_SHARED_MEMORY)
TEST(core_shared_memory, shared_memory_allocator)
{
  EXPECT_FALSE(axom::isSharedMemoryAllocator(axom::getDefaultAllocatorID()));
  EXPECT_TRUE(axom::isSharedMemoryAllocator(axom::getSharedMemoryAllocatorID()));

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const int sizes[] = {1024, 4096, 1000000};
  for(int si = 0; si < 3; si++)
  {
    // Allocate shared memory
    const int N = sizes[si];
    int *buffer = axom::allocate<int>(N, axom::getSharedMemoryAllocatorID());

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
