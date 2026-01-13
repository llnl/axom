// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/sidre/core/ConduitMemory.hpp"

namespace axom
{
namespace sidre
{

std::map<int, std::shared_ptr<ConduitMemory>> ConduitMemory::s_axomToInstance;
std::map<conduit::index_t, std::shared_ptr<ConduitMemory>> ConduitMemory::s_conduitToInstance;
const conduit::index_t ConduitMemory::s_defaultConduitId = conduit::Node().allocator();

void ConduitMemory::privateRegisterAllocator()
{
  using conduit::utils::register_allocator;
  auto deallocator = [](void* ptr) -> void {
    char* cPtr = (char*)(ptr);
    axom::deallocate<char>(cPtr);
  };
  m_deallocCallback = deallocator;
#if defined(AXOM_CONDUIT_USES_STD_FUNCTION)
  m_allocCallback = [=](size_t itemCount, size_t itemByteSize) -> void* {
    void* ptr = axom::allocate<char>(itemCount * itemByteSize, m_axomId);
    return ptr;
  };
  m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
#else
  /*
   * Note: conduit-0.9.4 allows the callbacks as std::function types.
   * Once we are there, we can use a single allocator, eliminating
   * the need for these if-else blocks.
   */
  if(m_axomId == MALLOC_ALLOCATOR_ID)
  {
    m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, MALLOC_ALLOCATOR_ID);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
  }
  else if(m_axomId == axom::INVALID_ALLOCATOR_ID)
  {
    m_allocCallback = nullptr;
    m_conduitId = -1;
  }
  else if(m_axomId == 0)
  {
    m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, 0);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
  }
  else if(m_axomId == 1)
  {
    m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, 1);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
  }
  else if(m_axomId == 2)
  {
    m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, 2);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
  }
  else if(m_axomId == 3)
  {
    m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, 3);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
  }
  else if(m_axomId == 4)
  {
    m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, 4);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
  }
  else if(m_axomId == 5)
  {
    m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, 5);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
  }
  else if(m_axomId == 6)
  {
    m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, 6);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
  }
  else if(m_axomId == 7)
  {
    m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, 7);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
  }
  else if(m_axomId == 8)
  {
    m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, 8);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
  }
  else if(m_axomId == 9)
  {
    m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, 9);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
  }
  else if(m_axomId == 10)
  {
    m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, 10);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
  }
  else if(m_axomId == 11)
  {
    m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, 11);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
  }
  else if(m_axomId == 12)
  {
    m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, 12);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
  }
  else if(m_axomId == 13)
  {
    m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, 13);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
  }
  else if(m_axomId == 14)
  {
    m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, 14);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
  }
  else if(m_axomId == 15)
  {
    m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, 15);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
  }
  else if(m_axomId == 15)
  {
    m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, 15);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
  }
  else if(m_axomId == 17)
  {
    m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, 17);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
  }
  else if(m_axomId == 18)
  {
    m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, 18);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
  }
  else if(m_axomId == 19)
  {
    m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, 19);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
  }
  else if(m_axomId == 20)
  {
    m_allocCallback = [](size_t itemCount, size_t itemByteSize) {
      void* ptr = axom::allocate<char>(itemCount * itemByteSize, 20);
      return ptr;
    };
    m_conduitId = register_allocator(m_allocCallback, m_deallocCallback);
  }
  else
  {
    std::cerr << "*** Work-around for conduit::utils::register_allocator "
                 "needs case for "
                 "m_axomId = "
              << std::to_string(m_axomId) << ".  Please add it to ConduitMemory.hpp.";
    axom::utilities::processAbort();
  }
#endif
}

const ConduitMemory& ConduitMemory::instanceForAxomId(int axomAllocId)
{
  if(s_axomToInstance.empty())
  {
    // Required one-time actions
    static auto axomMemcopy = [](void* dst, const void* src, size_t byteCount) {
      axom::copy(dst, src, byteCount);
    };
    static auto axomMemset = [](void* ptr, int value, size_t count) {
      if(axom::getAllocatorIDFromPointer(ptr) == axom::MALLOC_ALLOCATOR_ID)
      {
        std::memset(ptr, value, count);
      }
      else
      {
#if defined(AXOM_USE_UMPIRE)
        umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
        rm.memset(ptr, value, count);
#else
        std::cerr << "*** Error: Unrecognized axom allocator id" << std::endl;
        axom::utilities::processAbort();
#endif
      }
    };
    conduit::utils::set_memcpy_handler(axomMemcopy);
    conduit::utils::set_memset_handler(axomMemset);
  }

  auto it = s_axomToInstance.find(axomAllocId);
  if(it == s_axomToInstance.end())
  {
    std::shared_ptr<ConduitMemory> newInstance(new ConduitMemory(axomAllocId));
    s_axomToInstance[axomAllocId] = newInstance;
    it = s_axomToInstance.insert({axomAllocId, newInstance}).first;

    auto conduitAllocId = newInstance->m_conduitId;
    assert(s_conduitToInstance.find(conduitAllocId) == s_conduitToInstance.end());
    s_conduitToInstance[conduitAllocId] = newInstance;
  }
  assert(it->first == axomAllocId);

  return *it->second;
}

const ConduitMemory& ConduitMemory::instanceForConduitId(conduit::index_t conduitAllocId)
{
  auto it = s_conduitToInstance.find(conduitAllocId);
  if(it == s_conduitToInstance.end())
  {
    // conduitAllocId doesn't have a corresponding axom allocator.
    return instanceForAxomId(axom::INVALID_ALLOCATOR_ID);
  }
  assert(it->first == conduitAllocId);

  return *it->second;
}

}  // end namespace sidre
}  // end namespace axom
