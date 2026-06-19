// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/memory_management.hpp"

#if defined(AXOM_USE_UMPIRE)
  #include "umpire/config.hpp"
  #include "umpire/Allocator.hpp"
  #include "umpire/ResourceManager.hpp"
#endif

#include <cstdlib>

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------

// This value is such that the 64Kb limit on device constant memory is not hit
// in check_alloc_realloc_free when reallocating to 3 * ARRAY_SIZE.
constexpr int ARRAY_SIZE = 5345;

struct ScopedDefaultAllocatorState
{
  ScopedDefaultAllocatorState()
    : m_defaultAllocator(axom::getDefaultAllocatorID())
    , m_defaultHostAllocator(axom::getDefaultHostAllocatorID())
  { }

  ~ScopedDefaultAllocatorState()
  {
    axom::setDefaultAllocator(m_defaultAllocator);
    axom::setDefaultHostAllocator(m_defaultHostAllocator);
  }

  int m_defaultAllocator;
  int m_defaultHostAllocator;
};

class CopyTest : public ::testing::TestWithParam<::testing::tuple<std::string, std::string>>
{
public:
  void SetUp() override
  {
#if defined(AXOM_USE_UMPIRE)
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
#if defined(AXOM_USE_UMPIRE)
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
#if defined(AXOM_USE_UMPIRE)
    else
    {
      auto source_allocator = rm.getAllocator(dst_string);
      dst_array = static_cast<int*>(source_allocator.allocate(ARRAY_SIZE * sizeof(int)));
    }
#endif
  }

  void TearDown() override
  {
#if defined(AXOM_USE_UMPIRE)
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
#if defined(AXOM_USE_UMPIRE)
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
#if defined(AXOM_USE_UMPIRE)
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

#if defined(AXOM_USE_UMPIRE)
void check_alloc_and_free(int allocatorID = axom::getDefaultAllocatorID(), bool hostAccessible = true)
#else
void check_alloc_and_free(bool hostAccessible = true)
#endif
{
  for(int size = 0; size <= ARRAY_SIZE; size = size * 2 + 1)
  {
#if defined(AXOM_USE_UMPIRE)
    int* buffer = axom::allocate<int>(size, allocatorID);

    if(size > 0)
    {
      EXPECT_EQ(allocatorID, axom::getAllocatorIDFromPointer(buffer));
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

#if defined(AXOM_USE_UMPIRE)
void check_alloc_realloc_free(int allocatorID = axom::getDefaultAllocatorID(),
                              bool hostAccessible = true)
#else
void check_alloc_realloc_free(bool hostAccessible = true)
#endif
{
  for(int size = 0; size <= ARRAY_SIZE; size = size * 2 + 1)
  {
    int buffer_size = size;

#if defined(AXOM_USE_UMPIRE)
    int* buffer = axom::allocate<int>(buffer_size, allocatorID);

    if(buffer_size > 0)
    {
      ASSERT_EQ(allocatorID, axom::getAllocatorIDFromPointer(buffer));
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
#if defined(AXOM_USE_UMPIRE)
    buffer = axom::reallocate<int>(buffer, buffer_size, allocatorID);
    if(buffer_size > 0)
    {
      ASSERT_EQ(allocatorID, axom::getAllocatorIDFromPointer(buffer));
    }
#else
    buffer = axom::reallocate<int>(buffer, buffer_size);
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
#if defined(AXOM_USE_UMPIRE)
    buffer = axom::reallocate<int>(buffer, buffer_size, allocatorID);
    if(buffer_size > 0)
    {
      ASSERT_EQ(allocatorID, axom::getAllocatorIDFromPointer(buffer));
    }
#else
    buffer = axom::reallocate<int>(buffer, buffer_size);
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

TEST(core_memory_management, memory_space_availability)
{
  EXPECT_TRUE(axom::isMemorySpaceAvailable(axom::MemorySpace::Malloc));
  EXPECT_TRUE(axom::isMemorySpaceAvailable(axom::MemorySpace::Dynamic));
  EXPECT_TRUE(axom::isMemorySpaceAvailable(axom::MemorySpace::Host));

#if defined(AXOM_USE_UMPIRE)

  #if defined(UMPIRE_ENABLE_DEVICE)
  EXPECT_TRUE(axom::isMemorySpaceAvailable(axom::MemorySpace::Device));
  #endif

  #if defined(UMPIRE_ENABLE_UM)
  EXPECT_TRUE(axom::isMemorySpaceAvailable(axom::MemorySpace::Unified));
  #endif

  #if defined(UMPIRE_ENABLE_PINNED)
  EXPECT_TRUE(axom::isMemorySpaceAvailable(axom::MemorySpace::Pinned));
  #endif

  #if defined(UMPIRE_ENABLE_CONST)
  EXPECT_TRUE(axom::isMemorySpaceAvailable(axom::MemorySpace::Constant));
  #endif

#endif
}

TEST(core_memory_management, get_allocator_id_from_memory_space)
{
  EXPECT_EQ(axom::MALLOC_ALLOCATOR_ID,
            axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Malloc));
  EXPECT_EQ(axom::getDefaultAllocatorID(),
            axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Dynamic));
  EXPECT_EQ(axom::getDefaultHostAllocatorID(),
            axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Host));
}

TEST(core_memory_management, allocator_memory_space_compatibility)
{
  EXPECT_TRUE(axom::isAllocatorCompatibleWithMemorySpace(axom::MALLOC_ALLOCATOR_ID,
                                                         axom::MemorySpace::Malloc));
  EXPECT_TRUE(
    axom::isAllocatorCompatibleWithMemorySpace(axom::MALLOC_ALLOCATOR_ID, axom::MemorySpace::Host));
  EXPECT_TRUE(axom::isAllocatorCompatibleWithMemorySpace(axom::MALLOC_ALLOCATOR_ID,
                                                         axom::MemorySpace::Dynamic));

#if defined(AXOM_USE_UMPIRE)
  const int platformHostAllocatorID = axom::getUmpireResourceAllocatorID(umpire::resource::Host);
  EXPECT_TRUE(
    axom::isAllocatorCompatibleWithMemorySpace(platformHostAllocatorID, axom::MemorySpace::Host));
  EXPECT_FALSE(
    axom::isAllocatorCompatibleWithMemorySpace(platformHostAllocatorID, axom::MemorySpace::Malloc));
#endif
}

TEST(core_memory_management, set_get_default_host_allocator)
{
  ScopedDefaultAllocatorState scopedState;
  const int defaultAllocatorID = axom::getDefaultAllocatorID();

#if defined(AXOM_USE_UMPIRE)
  const int platformHostAllocatorID = axom::getUmpireResourceAllocatorID(umpire::resource::Host);
  EXPECT_EQ(axom::MALLOC_ALLOCATOR_ID, axom::getDefaultHostAllocatorID());
#else
  const int platformHostAllocatorID = axom::MALLOC_ALLOCATOR_ID;
  EXPECT_EQ(platformHostAllocatorID, axom::getDefaultHostAllocatorID());
#endif

  axom::setDefaultHostAllocator(axom::MemorySpace::Malloc);
  EXPECT_EQ(axom::MALLOC_ALLOCATOR_ID, axom::getDefaultHostAllocatorID());
  EXPECT_EQ(axom::MALLOC_ALLOCATOR_ID, axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Host));
  EXPECT_EQ(defaultAllocatorID, axom::getDefaultAllocatorID());
  EXPECT_EQ(defaultAllocatorID, axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Dynamic));

  axom::setDefaultHostAllocator(axom::MemorySpace::Host);
  EXPECT_EQ(platformHostAllocatorID, axom::getDefaultHostAllocatorID());
  EXPECT_EQ(platformHostAllocatorID, axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Host));
  EXPECT_EQ(defaultAllocatorID, axom::getDefaultAllocatorID());

#if defined(AXOM_USE_UMPIRE)
  axom::setDefaultHostAllocator(axom::MemorySpace::Malloc);
  axom::setDefaultAllocator(axom::MemorySpace::Host);
  EXPECT_EQ(platformHostAllocatorID, axom::getDefaultAllocatorID());
  EXPECT_EQ(axom::MALLOC_ALLOCATOR_ID, axom::getDefaultHostAllocatorID());
#endif
}

#if defined(AXOM_USE_UMPIRE)

bool hostAllocationUsesExpectedAllocator(int selectedHostAllocId,
                                         int expectedAllocId,
                                         bool setGlobalDefaultToHost = false)
{
  axom::setDefaultHostAllocator(selectedHostAllocId);

  if(setGlobalDefaultToHost)
  {
    axom::setDefaultAllocator(axom::MemorySpace::Host);
  }

  const int resolvedHostAllocatorID = axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Host);
  if(resolvedHostAllocatorID != selectedHostAllocId)
  {
    return false;
  }

  int* buffer = axom::allocate<int>(ARRAY_SIZE, resolvedHostAllocatorID);
  if(buffer == nullptr)
  {
    return false;
  }

  const bool allocatorMatches = axom::getAllocatorIDFromPointer(buffer) == expectedAllocId;
  axom::deallocate(buffer);
  return allocatorMatches;
}

TEST(core_memory_management, host_space_allocation_uses_umpire_host_default)
{
  EXPECT_EXIT(([]() {
                if(!hostAllocationUsesExpectedAllocator(
                     axom::getUmpireResourceAllocatorID(umpire::resource::Host),
                     axom::getUmpireResourceAllocatorID(umpire::resource::Host)))
                {
                  std::exit(1);
                }

                std::exit(0);
              })(),
              ::testing::ExitedWithCode(0),
              "");
}

TEST(core_memory_management, host_space_allocation_uses_malloc_host_default)
{
  EXPECT_EXIT(
    ([]() {
      if(!hostAllocationUsesExpectedAllocator(axom::MALLOC_ALLOCATOR_ID, axom::MALLOC_ALLOCATOR_ID))
      {
        std::exit(1);
      }

      std::exit(0);
    })(),
    ::testing::ExitedWithCode(0),
    "");
}

TEST(core_memory_management, host_space_allocation_ignores_global_default_allocator)
{
  EXPECT_EXIT(
    ([]() {
      const int hostAllocatorID = axom::getUmpireResourceAllocatorID(umpire::resource::Host);
      axom::setDefaultHostAllocator(axom::MemorySpace::Malloc);
      axom::setDefaultAllocator(axom::MemorySpace::Host);

      if(axom::getDefaultAllocatorID() != hostAllocatorID)
      {
        std::exit(1);
      }

      if(!hostAllocationUsesExpectedAllocator(axom::MALLOC_ALLOCATOR_ID, axom::MALLOC_ALLOCATOR_ID))
      {
        std::exit(1);
      }

      std::exit(0);
    })(),
    ::testing::ExitedWithCode(0),
    "");
}

TEST(core_memory_management, changing_default_host_allocator_after_host_allocation_fails)
{
  EXPECT_DEATH_IF_SUPPORTED(
    []() {
      const int hostAllocatorID =
        axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Host);
      axom::setDefaultHostAllocator(hostAllocatorID);
      int* buffer =
        axom::allocate<int>(ARRAY_SIZE, axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Host));
      AXOM_UNUSED_VAR(buffer);
      axom::setDefaultHostAllocator(axom::MALLOC_ALLOCATOR_ID);
    }(),
    "Default host allocator cannot be changed");
}

TEST(core_memory_management, set_get_default_memory_space)
{
  ScopedDefaultAllocatorState scopedState;
  const int HostAllocatorID = axom::getDefaultHostAllocatorID();
  const int platformHostAllocatorID = axom::getUmpireResourceAllocatorID(umpire::resource::Host);
  EXPECT_EQ(HostAllocatorID, axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Host));

  #if defined(AXOM_USE_GPU)

    #if defined(UMPIRE_ENABLE_PINNED)
  const int PinnedAllocatorID = axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Pinned);
  axom::setDefaultAllocator(axom::MemorySpace::Pinned);
  EXPECT_EQ(PinnedAllocatorID, axom::getDefaultAllocatorID());
    #endif

    #if defined(UMPIRE_ENABLE_DEVICE)
  const int DeviceAllocatorID = axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Device);
  axom::setDefaultAllocator(axom::MemorySpace::Device);
  EXPECT_EQ(DeviceAllocatorID, axom::getDefaultAllocatorID());
    #endif

    #if defined(UMPIRE_ENABLE_CONST)
  const int ConstantAllocatorID = axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Constant);
  axom::setDefaultAllocator(axom::MemorySpace::Constant);
  EXPECT_EQ(ConstantAllocatorID, axom::getDefaultAllocatorID());
    #endif

    #if defined(UMPIRE_ENABLE_UM)
  const int UnifiedAllocatorID = axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Unified);
  axom::setDefaultAllocator(axom::MemorySpace::Unified);
  EXPECT_EQ(UnifiedAllocatorID, axom::getDefaultAllocatorID());
    #endif

  #endif  // AXOM_USE_GPU

  axom::setDefaultAllocator(axom::MemorySpace::Host);
  EXPECT_EQ(platformHostAllocatorID, axom::getDefaultAllocatorID());
}
#endif /* AXOM_USE_UMPIRE */

//------------------------------------------------------------------------------
TEST(core_memory_management, alloc_free)
{
#if defined(AXOM_USE_UMPIRE)

  constexpr bool HOST_ACCESSIBLE = true;

  const int HostAllocatorID = axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Host);
  check_alloc_and_free(HostAllocatorID, HOST_ACCESSIBLE);

  #if defined(AXOM_USE_GPU)

  constexpr bool NOT_HOST_ACCESSIBLE = false;

    #if defined(UMPIRE_ENABLE_PINNED)
  const int PinnedAllocatorID = axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Pinned);
  check_alloc_and_free(PinnedAllocatorID, HOST_ACCESSIBLE);
    #endif

    #if defined(UMPIRE_ENABLE_DEVICE)
  const int DeviceAllocatorID = axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Device);
  check_alloc_and_free(DeviceAllocatorID, NOT_HOST_ACCESSIBLE);
    #endif

    #if defined(UMPIRE_ENABLE_CONST)
  const int ConstantAllocatorID = axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Constant);
  check_alloc_and_free(ConstantAllocatorID, NOT_HOST_ACCESSIBLE);
    #endif

    #if defined(UMPIRE_ENABLE_UM)
  const int UnifiedAllocatorID = axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Unified);
  check_alloc_and_free(UnifiedAllocatorID, HOST_ACCESSIBLE);
    #endif

  #endif  // AXOM_USE_GPU

#endif  // AXOM_USE_UMPIRE

  check_alloc_and_free();
}

//------------------------------------------------------------------------------
TEST(core_memory_management, alloc_realloc_free)
{
#if defined(AXOM_USE_UMPIRE)

  constexpr bool HOST_ACCESSIBLE = true;

  const int HostAllocatorID = axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Host);
  check_alloc_realloc_free(HostAllocatorID, HOST_ACCESSIBLE);

  #if defined(AXOM_USE_GPU)

  constexpr bool NOT_HOST_ACCESSIBLE = false;

    #if defined(UMPIRE_ENABLE_PINNED)
  const int PinnedAllocatorID = axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Pinned);
  check_alloc_realloc_free(PinnedAllocatorID, HOST_ACCESSIBLE);
    #endif

    #if defined(UMPIRE_ENABLE_DEVICE)
  const int DeviceAllocatorID = axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Device);
  check_alloc_realloc_free(DeviceAllocatorID, NOT_HOST_ACCESSIBLE);
    #endif

    // Umpire doesn't allow reallocation of Constant memory.
    // check_alloc_realloc_free( axom::getAllocator( umpire::resource::Constant ),
    // false );

    #if defined(UMPIRE_ENABLE_UM)
  const int UnifiedAllocatorID = axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Unified);
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
#if defined(AXOM_USE_UMPIRE)
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

TEST(core_memory_management, foreign_malloc_to_umpire_reallocate_fails)
{
  EXPECT_DEATH_IF_SUPPORTED(
    []() {
      constexpr std::size_t localN = 5;
      const int localUmpireHostAllocId =
        axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Host);
      int* foreignBuffer = static_cast<int*>(std::malloc(localN * sizeof(int)));
      axom::reallocate(foreignBuffer, localN + 1, localUmpireHostAllocId);
    }(),
    "Cannot reallocate across allocator backends");
}

TEST(core_memory_management, axom_malloc_to_umpire_reallocate_fails)
{
  EXPECT_DEATH_IF_SUPPORTED(
    []() {
      constexpr std::size_t localN = 5;
      const int localUmpireHostAllocId =
        axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Host);
      int* buffer = axom::allocate<int>(localN, axom::MALLOC_ALLOCATOR_ID);
      axom::reallocate(buffer, localN + 1, localUmpireHostAllocId);
    }(),
    "Cannot reallocate across allocator backends");
}

TEST(core_memory_management, umpire_to_axom_malloc_reallocate_fails)
{
  EXPECT_DEATH_IF_SUPPORTED(
    []() {
      constexpr std::size_t localN = 5;
      const int localUmpireHostAllocId =
        axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Host);
      int* buffer = axom::allocate<int>(localN, localUmpireHostAllocId);
      axom::reallocate(buffer, localN + 1, axom::MALLOC_ALLOCATOR_ID);
    }(),
    "Cannot reallocate across allocator backends");
}
#endif

//------------------------------------------------------------------------------
TEST(core_memory_management, test_fill)
{
  // Allocator ids to test.
  std::vector<int> allocIds(1, axom::MALLOC_ALLOCATOR_ID);
#if defined(AXOM_USE_UMPIRE)
  allocIds.push_back(axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Host));
  #ifdef AXOM_USE_GPU
    #if defined(UMPIRE_ENABLE_DEVICE)
  allocIds.push_back(axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Device));
    #endif
    #if defined(UMPIRE_ENABLE_UM)
  allocIds.push_back(axom::getAllocatorIDFromMemorySpace(axom::MemorySpace::Unified));
    #endif
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
