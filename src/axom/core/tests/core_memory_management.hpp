// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/memory_management.hpp"

#ifdef AXOM_USE_UMPIRE
  #include "umpire/config.hpp"
  #include "umpire/Allocator.hpp"
  #include "umpire/ResourceManager.hpp"
#endif

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------

// This value is such that the 64Kb limit on device constant memory is not hit
// in check_alloc_realloc_free when reallocating to 3 * ARRAY_SIZE.
constexpr int ARRAY_SIZE = 5345;

class CopyTest : public ::testing::TestWithParam<::testing::tuple<std::string, std::string>>
{
public:
  void SetUp() override
  {
#ifdef AXOM_USE_UMPIRE
    umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
#endif

    src_string = ::testing::get<0>(GetParam());
    dst_string = ::testing::get<1>(GetParam());

    if(src_string == "NEW")
    {
      src_array = new int[ARRAY_SIZE];
    }
    else if(src_string == "MALLOC")
    {
      src_array = static_cast<int*>(std::malloc(ARRAY_SIZE * sizeof(int)));
    }
    else if(src_string == "STATIC")
    {
      src_array = m_static_src_array;
    }
#ifdef AXOM_USE_UMPIRE
    else
    {
      auto source_allocator = rm.getAllocator(src_string);
      src_array = static_cast<int*>(source_allocator.allocate(ARRAY_SIZE * sizeof(int)));
    }
#endif

    if(dst_string == "NEW")
    {
      dst_array = new int[ARRAY_SIZE];
    }
    else if(dst_string == "MALLOC")
    {
      dst_array = static_cast<int*>(std::malloc(ARRAY_SIZE * sizeof(int)));
    }
    else if(dst_string == "STATIC")
    {
      dst_array = m_static_dst_array;
    }
#ifdef AXOM_USE_UMPIRE
    else
    {
      auto source_allocator = rm.getAllocator(dst_string);
      dst_array = static_cast<int*>(source_allocator.allocate(ARRAY_SIZE * sizeof(int)));
    }
#endif
  }

  void TearDown() override
  {
#ifdef AXOM_USE_UMPIRE
    auto& rm = umpire::ResourceManager::getInstance();
#endif

    if(src_string == "NEW")
    {
      delete[] src_array;
    }
    else if(src_string == "MALLOC")
    {
      std::free(src_array);
    }
#ifdef AXOM_USE_UMPIRE
    else if(src_string != "STATIC")
    {
      rm.deallocate(src_array);
    }
#endif

    if(dst_string == "NEW")
    {
      delete[] dst_array;
    }
    else if(dst_string == "MALLOC")
    {
      std::free(dst_array);
    }
#ifdef AXOM_USE_UMPIRE
    else if(dst_string != "STATIC")
    {
      rm.deallocate(dst_array);
    }
#endif
  }

  std::string src_string;
  std::string dst_string;

  int* src_array = nullptr;
  int* dst_array = nullptr;
  int host_array[ARRAY_SIZE];

private:
  int m_static_src_array[ARRAY_SIZE];
  int m_static_dst_array[ARRAY_SIZE];
};

#ifdef AXOM_USE_UMPIRE
void check_alloc_and_free(int allocatorID = axom::getDefaultAllocatorID(), bool hostAccessible = true)
#else
void check_alloc_and_free(bool hostAccessible = true)
#endif
{
  for(int size = 0; size <= ARRAY_SIZE; size = size * 2 + 1)
  {
#ifdef AXOM_USE_UMPIRE
    int* buffer = axom::allocate<int>(size, allocatorID);

    if(size > 0)
    {
      umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
      EXPECT_EQ(allocatorID, rm.getAllocator(buffer).getId());
    }
#else
    int* buffer = axom::allocate<int>(size);
#endif

    if(hostAccessible)
    {
      for(int i = 0; i < size; ++i)
      {
        buffer[i] = i;
      }

      for(int i = 0; i < size; ++i)
      {
        EXPECT_EQ(buffer[i], i);
      }
    }

    axom::deallocate(buffer);
    EXPECT_TRUE(buffer == nullptr);
  }
}

#ifdef AXOM_USE_UMPIRE
void check_alloc_realloc_free(int allocatorID = axom::getDefaultAllocatorID(),
                              bool hostAccessible = true)
#else
void check_alloc_realloc_free(bool hostAccessible = true)
#endif
{
  for(int size = 0; size <= ARRAY_SIZE; size = size * 2 + 1)
  {
    int buffer_size = size;

#ifdef AXOM_USE_UMPIRE
    int* buffer = axom::allocate<int>(buffer_size, allocatorID);

    umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
    if(buffer_size > 0)
    {
      ASSERT_EQ(allocatorID, rm.getAllocator(buffer).getId());
    }
#else
    int* buffer = axom::allocate<int>(buffer_size);
#endif

    if(hostAccessible)
    {
      // Populate the buffer.
      for(int i = 0; i < buffer_size; ++i)
      {
        buffer[i] = i;
      }

      // Check the values.
      for(int i = 0; i < buffer_size; ++i)
      {
        EXPECT_EQ(buffer[i], i);
      }
    }

    // Reallocate to a larger size.
    buffer_size *= 3;
    buffer = axom::reallocate<int>(buffer, buffer_size);
#ifdef AXOM_USE_UMPIRE
    if(buffer_size > 0)
    {
      ASSERT_EQ(allocatorID, rm.getAllocator(buffer).getId());
    }
#endif

    if(hostAccessible)
    {
      // Populate the new values.
      for(int i = size; i < buffer_size; ++i)
      {
        buffer[i] = i;
      }

      // Check all the values.
      for(int i = 0; i < buffer_size; ++i)
      {
        EXPECT_EQ(buffer[i], i);
      }
    }

    // Reallocate to a smaller size.
    buffer_size /= 5;
    buffer = axom::reallocate<int>(buffer, buffer_size);
#ifdef AXOM_USE_UMPIRE
    if(buffer_size > 0)
    {
      ASSERT_EQ(allocatorID, rm.getAllocator(buffer).getId());
    }
#endif

    if(hostAccessible)
    {
      // Check all the values.
      for(int i = 0; i < buffer_size; ++i)
      {
        EXPECT_EQ(buffer[i], i);
      }
    }

    // Free
    axom::deallocate(buffer);
    EXPECT_TRUE(buffer == nullptr);
  }
}

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

#ifdef AXOM_USE_UMPIRE

TEST(core_memory_management, set_get_default_memory_space)
{
  const int HostAllocatorID = axom::getUmpireResourceAllocatorID(umpire::resource::Host);
  EXPECT_EQ(HostAllocatorID, axom::getDefaultAllocatorID());

  #if defined(AXOM_USE_GPU)

    #ifdef UMPIRE_ENABLE_PINNED
  const int PinnedAllocatorID = axom::getUmpireResourceAllocatorID(umpire::resource::Pinned);

  axom::setDefaultAllocator(PinnedAllocatorID);
  EXPECT_EQ(PinnedAllocatorID, axom::getDefaultAllocatorID());
    #endif

    #ifdef UMPIRE_ENABLE_DEVICE
  const int DeviceAllocatorID = axom::getUmpireResourceAllocatorID(umpire::resource::Device);
  axom::setDefaultAllocator(DeviceAllocatorID);
  EXPECT_EQ(DeviceAllocatorID, axom::getDefaultAllocatorID());
    #endif

    #ifdef UMPIRE_ENABLE_CONST
  const int ConstantAllocatorID = axom::getUmpireResourceAllocatorID(umpire::resource::Constant);
  axom::setDefaultAllocator(ConstantAllocatorID);
  EXPECT_EQ(ConstantAllocatorID, axom::getDefaultAllocatorID());
    #endif

    #ifdef UMPIRE_ENABLE_UM
  const int UnifiedAllocatorID = axom::getUmpireResourceAllocatorID(umpire::resource::Unified);
  axom::setDefaultAllocator(UnifiedAllocatorID);
  EXPECT_EQ(UnifiedAllocatorID, axom::getDefaultAllocatorID());
    #endif

  #endif  // AXOM_USE_GPU

  axom::setDefaultAllocator(HostAllocatorID);
  EXPECT_EQ(HostAllocatorID, axom::getDefaultAllocatorID());
}
#endif /* AXOM_USE_UMPIRE */

//------------------------------------------------------------------------------
TEST(core_memory_management, alloc_free)
{
#ifdef AXOM_USE_UMPIRE

  constexpr bool HOST_ACCESSIBLE = true;

  const int HostAllocatorID = axom::getUmpireResourceAllocatorID(umpire::resource::Host);
  check_alloc_and_free(HostAllocatorID, HOST_ACCESSIBLE);

  #if defined(AXOM_USE_GPU)

  constexpr bool NOT_HOST_ACCESSIBLE = false;

    #ifdef UMPIRE_ENABLE_PINNED
  const int PinnedAllocatorID = axom::getUmpireResourceAllocatorID(umpire::resource::Pinned);
  check_alloc_and_free(PinnedAllocatorID, HOST_ACCESSIBLE);
    #endif

    #ifdef UMPIRE_ENABLE_DEVICE
  const int DeviceAllocatorID = axom::getUmpireResourceAllocatorID(umpire::resource::Device);
  check_alloc_and_free(DeviceAllocatorID, NOT_HOST_ACCESSIBLE);
    #endif

    #ifdef UMPIRE_ENABLE_CONST
  const int ConstantAllocatorID = axom::getUmpireResourceAllocatorID(umpire::resource::Constant);
  check_alloc_and_free(ConstantAllocatorID, NOT_HOST_ACCESSIBLE);
    #endif

    #ifdef UMPIRE_ENABLE_UM
  const int UnifiedAllocatorID = axom::getUmpireResourceAllocatorID(umpire::resource::Unified);
  check_alloc_and_free(UnifiedAllocatorID, HOST_ACCESSIBLE);
    #endif

  #endif  // AXOM_USE_GPU

#endif  // AXOM_USE_UMPIRE

  check_alloc_and_free();
}

//------------------------------------------------------------------------------
TEST(core_memory_management, alloc_realloc_free)
{
#ifdef AXOM_USE_UMPIRE

  constexpr bool HOST_ACCESSIBLE = true;

  const int HostAllocatorID = axom::getUmpireResourceAllocatorID(umpire::resource::Host);
  check_alloc_realloc_free(HostAllocatorID, HOST_ACCESSIBLE);

  #if defined(AXOM_USE_GPU)

  constexpr bool NOT_HOST_ACCESSIBLE = false;

    #ifdef UMPIRE_ENABLE_PINNED
  const int PinnedAllocatorID = axom::getUmpireResourceAllocatorID(umpire::resource::Pinned);
  check_alloc_realloc_free(PinnedAllocatorID, HOST_ACCESSIBLE);
    #endif

    #ifdef UMPIRE_ENABLE_DEVICE
  const int DeviceAllocatorID = axom::getUmpireResourceAllocatorID(umpire::resource::Device);
  check_alloc_realloc_free(DeviceAllocatorID, NOT_HOST_ACCESSIBLE);
    #endif

    // Umpire doesn't allow reallocation of Constant memory.
    // check_alloc_realloc_free( axom::getAllocator( umpire::resource::Constant ),
    // false );

    #ifdef UMPIRE_ENABLE_UM

  const int UnifiedAllocatorID = axom::getUmpireResourceAllocatorID(umpire::resource::Unified);
  check_alloc_realloc_free(UnifiedAllocatorID, HOST_ACCESSIBLE);

    #endif

  #endif /* AXOM_USE_GPU */

#endif /* AXOM_USE_UMPIRE */

  check_alloc_realloc_free();
}

TEST_P(CopyTest, Copy)
{
  std::cout << "SRC = " << src_string << ", DST = " << dst_string << std::endl;
  for(int i = 0; i < ARRAY_SIZE; ++i)
  {
    host_array[i] = i;
  }

  axom::copy(src_array, host_array, ARRAY_SIZE * sizeof(int));
  axom::copy(dst_array, src_array, ARRAY_SIZE * sizeof(int));

  for(int i = 0; i < ARRAY_SIZE; ++i)
  {
    host_array[i] = -i;
  }

  axom::copy(host_array, src_array, ARRAY_SIZE * sizeof(int));

  for(int i = 0; i < ARRAY_SIZE; ++i)
  {
    ASSERT_EQ(host_array[i], i);
  }
}

const std::string copy_locations[] = {"NEW",
                                      "MALLOC",
                                      "STATIC"
#if defined(AXOM_USE_UMPIRE)
                                      ,
                                      "HOST"
  #if defined(UMPIRE_ENABLE_DEVICE)
                                      ,
                                      "DEVICE"
  #endif
  #if defined(UMPIRE_ENABLE_UM)
                                      ,
                                      "UM"
  #endif
  #if defined(UMPIRE_ENABLE_PINNED)
                                      ,
                                      "PINNED"
  #endif
#endif
};

INSTANTIATE_TEST_SUITE_P(core_memory_management,
                         CopyTest,
                         ::testing::Combine(::testing::ValuesIn(copy_locations),
                                            ::testing::ValuesIn(copy_locations)));

//------------------------------------------------------------------------------
TEST(core_memory_management, basic_alloc_realloc_dealloc)
{
  // A basic test for axom's allocate, reallocate and deallocate functionality
  // for the default memory allocator

  constexpr std::size_t N = 5;

  int* buf = nullptr;

  // allocate an array of size N
  buf = axom::allocate<int>(N);
  EXPECT_NE(buf, nullptr);

  // free the array
  axom::deallocate<int>(buf);
  EXPECT_EQ(buf, nullptr);

  // reallocate array to size 0
  buf = axom::reallocate<int>(buf, 0);
  EXPECT_NE(buf, nullptr);

  // reallocate array to size N
  buf = axom::reallocate<int>(buf, N);
  EXPECT_NE(buf, nullptr);

  // reallocate array to size 0
  buf = axom::reallocate<int>(buf, 0);
  EXPECT_NE(buf, nullptr);

  // reallocate array to size 0, again
  buf = axom::reallocate<int>(buf, 0);
  EXPECT_NE(buf, nullptr);

  // Finally, free the array
  axom::deallocate<int>(buf);
  EXPECT_EQ(buf, nullptr);
}

//------------------------------------------------------------------------------
#ifdef AXOM_USE_UMPIRE
TEST(core_memory_management, allocator_id_from_pointer)
{
  constexpr std::size_t N = 5;

  int* buf = nullptr;

  // Allocate through allocator.
  buf = axom::allocate<int>(N);
  EXPECT_NE(buf, nullptr);
  int id = axom::getAllocatorIDFromPointer(buf);
  EXPECT_EQ(id, axom::getDefaultAllocatorID());
  axom::deallocate<int>(buf);

  // Allocate directly (not through allocator).
  buf = (int*)std::malloc(N * sizeof(int));
  id = axom::getAllocatorIDFromPointer(buf);
  EXPECT_EQ(id, axom::MALLOC_ALLOCATOR_ID);
  std::free(buf);
}
#endif

//------------------------------------------------------------------------------
TEST(core_memory_management, interspace_reallocation)
{
  // Allocator ids to test.
  std::vector<int> allocIds(1, axom::MALLOC_ALLOCATOR_ID);
#ifdef AXOM_USE_UMPIRE
  allocIds.push_back(axom::detail::getAllocatorID<axom::MemorySpace::Host>());
  #ifdef AXOM_USE_GPU
  allocIds.push_back(axom::detail::getAllocatorID<axom::MemorySpace::Device>());
  allocIds.push_back(axom::detail::getAllocatorID<axom::MemorySpace::Unified>());
    // Does it make sense to check Pinned and Constant memory spaces?
  #endif
#endif

  // We'll allocate N items, reallocate to K items, reallocate back to N.
  constexpr std::size_t N = 5;
  constexpr std::size_t K = 8;
  constexpr std::size_t maxNK = std::max(N, K);
  constexpr std::size_t minNK = std::min(N, K);

  // origOnHost and tempOnHost are for initialization and results-checking on host.
  int* origOnHost = axom::allocate<int>(maxNK, axom::MALLOC_ALLOCATOR_ID);
  for(std::size_t i = 0; i < maxNK; ++i)
  {
    origOnHost[i] = 100 + i;
  }
  int* tempOnHost = axom::allocate<int>(maxNK, axom::MALLOC_ALLOCATOR_ID);

  // Count differences between origOnHost and tempOnHost.
  auto countDiffs = [=]() {
    std::size_t diffCount = 0;
    for(std::size_t j = 0; j < minNK; ++j)
    {
      diffCount += tempOnHost[j] != origOnHost[j];
    }
    return diffCount;
  };

  std::size_t diffCount = 0;
  for(auto srcAllocId : allocIds)
  {
    for(auto dstAllocId : allocIds)
    {
      std::cout << "Testing allocator ids " << srcAllocId << " and " << dstAllocId << std::endl;
      // For each combination of srcAllocId and dstAllocId,
      // allocate src, reallocate to dst, reallocate back to src.

      int* src = axom::allocate<int>(N, srcAllocId);
      axom::copy(src, origOnHost, N * sizeof(int));

      int* dst = axom::reallocate(src, K, dstAllocId);
      axom::copy(tempOnHost, dst, N * sizeof(int));
      diffCount = countDiffs();
      EXPECT_EQ(diffCount, 0);

      src = axom::reallocate(dst, N, srcAllocId);
      axom::copy(tempOnHost, src, N * sizeof(int));
      diffCount = countDiffs();
      EXPECT_EQ(diffCount, 0);

      axom::deallocate(src);
    }
  }

  axom::deallocate(origOnHost);
  axom::deallocate(tempOnHost);
}

//------------------------------------------------------------------------------
TEST(core_memory_management, test_fill)
{
  // Allocator ids to test.
  std::vector<int> allocIds(1, axom::MALLOC_ALLOCATOR_ID);
#ifdef AXOM_USE_UMPIRE
  allocIds.push_back(axom::detail::getAllocatorID<axom::MemorySpace::Host>());
  #ifdef AXOM_USE_GPU
  allocIds.push_back(axom::detail::getAllocatorID<axom::MemorySpace::Device>());
  allocIds.push_back(axom::detail::getAllocatorID<axom::MemorySpace::Unified>());
    // Does it make sense to check Pinned and Constant memory spaces?
  #endif
#endif

  constexpr std::size_t N = 500;
  int* hostData = axom::allocate<int>(N, axom::MALLOC_ALLOCATOR_ID);
  int iteration = 0;
  for(auto allocId : allocIds)
  {
    std::cout << "Testing allocator id " << allocId << std::endl;

    // Allocate buffer using allocId and fill it.
    int* buffer = axom::allocate<int>(N, allocId);
    const int fillValue = 12345 + iteration;
    axom::fill(buffer, N, fillValue);

    // Copy back to host
    axom::copy(hostData, buffer, N * sizeof(int));

    // Make sure elements have the right fill value.
    for(std::size_t i = 0; i < N; i++)
    {
      EXPECT_EQ(hostData[i], fillValue);
    }

    axom::deallocate(buffer);
    iteration++;
  }

  axom::deallocate(hostData);
}
