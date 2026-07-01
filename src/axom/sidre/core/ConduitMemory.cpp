// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
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
  const auto axomId = m_axomId;
  m_allocCallback = [=](size_t itemCount, size_t itemByteSize) -> void* {
    void* ptr = axom::allocate<char>(itemCount * itemByteSize, axomId);
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

namespace
{
axom::utilities::CheckSum checksumNodeMetadata(const conduit::Node& n)
{
  auto cs = axom::utilities::CheckSum {0};
  // Give each metadata component different coefficients so they contribute differently
  cs += axom::utilities::checksum(static_cast<conduit::index_t>(n.dtype().id()), 2.0);
  cs += axom::utilities::checksum(n.number_of_children(), 3.0);
  cs += axom::utilities::checksum(n.dtype().number_of_elements(), 5.0);
  return cs;
}

/// Operate on conduit::DataArray so we can handle strided data.
template <typename T>
axom::utilities::CheckSum checksumArray(const conduit::DataArray<T>& arr)
{
  return axom::utilities::calculateChecksum(
           [=](axom::IndexType i) { return static_cast<axom::utilities::CheckSum>(arr[i]); },
           arr.number_of_elements());
}

/// Compute a checksum on the conduit tree \a n.
axom::utilities::CheckSum checksumImpl(const conduit::Node& n, bool include_name)
{
  auto cs = axom::utilities::CheckSum {0};
  if(include_name)
  {
    std::string name(n.name());
    axom::ArrayView<const char> view(name.data(), name.size());
    cs = axom::utilities::checksum(view);
  }

  cs += checksumNodeMetadata(n);

  if(n.number_of_children() > 0)
  {
    if(n.dtype().is_list())
    {
      cs += axom::utilities::calculateChecksum(
        [&](axom::IndexType i) -> axom::utilities::CheckSum {
          return checksumImpl(n[static_cast<conduit::index_t>(i)], true);
        },
        n.number_of_children());
    }
    else
    {
      for(conduit::index_t i = 0; i < n.number_of_children(); i++)
      {
        cs += checksumImpl(n[i], true);
      }
    }
  }
  else
  {
    if(n.dtype().is_string())
    {
      axom::ArrayView<const char> view(static_cast<const char*>(n.data_ptr()),
                                       n.dtype().number_of_elements());
      cs += axom::utilities::checksum(view);
    }
    else if(n.dtype().is_int8())
    {
      cs += checksumArray(n.as_int8_array());
    }
    else if(n.dtype().is_int16())
    {
      cs += checksumArray(n.as_int16_array());
    }
    else if(n.dtype().is_int32())
    {
      cs += checksumArray(n.as_int32_array());
    }
    else if(n.dtype().is_int64())
    {
      cs += checksumArray(n.as_int64_array());
    }
    else if(n.dtype().is_uint8())
    {
      cs += checksumArray(n.as_uint8_array());
    }
    else if(n.dtype().is_uint16())
    {
      cs += checksumArray(n.as_uint16_array());
    }
    else if(n.dtype().is_uint32())
    {
      cs += checksumArray(n.as_uint32_array());
    }
    else if(n.dtype().is_uint64())
    {
      cs += checksumArray(n.as_uint64_array());
    }
    else if(n.dtype().is_index_t())
    {
      cs += checksumArray(n.as_index_t_array());
    }
    else if(n.dtype().is_float32())
    {
      cs += checksumArray(n.as_float32_array());
    }
    else if(n.dtype().is_float64())
    {
      cs += checksumArray(n.as_float64_array());
    }
    else if(n.dtype().is_empty() || n.dtype().is_object() || n.dtype().is_list())
    {
      // no-op
    }
  }
  return cs;
}
} // end namespace

axom::utilities::CheckSum checksum(const conduit::Node& n, axom::utilities::ScaleFactor scaleFactor, bool include_name)
{
  return checksumImpl(n, include_name) * static_cast<axom::utilities::CheckSum>(scaleFactor);
}

axom::utilities::CheckSum checksum(const conduit::Node& n, bool include_name)
{
  return checksumImpl(n, include_name);
}

}  // end namespace sidre
}  // end namespace axom
