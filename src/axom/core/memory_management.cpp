// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/memory_management.hpp"

#if defined(AXOM_USE_UMPIRE)
  #include "axom/fmt.hpp"
  #include "umpire/Umpire.hpp"
  #include "umpire/util/MemoryResourceTraits.hpp"
  #if defined(AXOM_USE_UMPIRE_SHARED_MEMORY)
    #include "umpire/strategy/NamedAllocationStrategy.hpp"
  #endif
#endif

namespace axom
{

namespace
{
const char* memorySpaceName(MemorySpace space) noexcept
{
  switch(space)
  {
  case MemorySpace::Malloc:
    return "Malloc";
  case MemorySpace::Dynamic:
    return "Dynamic";
  case MemorySpace::Host:
    return "Host";
  case MemorySpace::Device:
    return "Device";
  case MemorySpace::Unified:
    return "Unified";
  case MemorySpace::Pinned:
    return "Pinned";
  case MemorySpace::Constant:
    return "Constant";
  }

  return "Unknown";
}

int platformHostAllocatorID() noexcept
{
#if defined(AXOM_USE_UMPIRE)
  return getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Host);
#else
  return MALLOC_ALLOCATOR_ID;
#endif
}

int& defaultHostAllocatorIDStorage() noexcept
{
  static int allocId = platformHostAllocatorID();
  return allocId;
}
}  // namespace

bool isMemorySpaceAvailable(MemorySpace space) noexcept
{
  switch(space)
  {
  case MemorySpace::Malloc:
  case MemorySpace::Dynamic:
  case MemorySpace::Host:
    return true;
  case MemorySpace::Device:
#if defined(AXOM_USE_UMPIRE) && defined(UMPIRE_ENABLE_DEVICE)
    return true;
#else
    return false;
#endif
  case MemorySpace::Unified:
#if defined(AXOM_USE_UMPIRE) && defined(UMPIRE_ENABLE_UM)
    return true;
#else
    return false;
#endif
  case MemorySpace::Pinned:
#if defined(AXOM_USE_UMPIRE) && defined(UMPIRE_ENABLE_PINNED)
    return true;
#else
    return false;
#endif
  case MemorySpace::Constant:
#if defined(AXOM_USE_UMPIRE) && defined(UMPIRE_ENABLE_CONST)
    return true;
#else
    return false;
#endif
  }

  return false;
}

int getAllocatorIDFromMemorySpace(MemorySpace space)
{
  switch(space)
  {
  case MemorySpace::Dynamic:
    return getDefaultAllocatorID();
  case MemorySpace::Malloc:
    return MALLOC_ALLOCATOR_ID;
  case MemorySpace::Host:
    return getDefaultHostAllocatorID();
  case MemorySpace::Device:
#if defined(AXOM_USE_UMPIRE) && defined(UMPIRE_ENABLE_DEVICE)
    return getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Device);
#else
    break;
#endif
  case MemorySpace::Unified:
#if defined(AXOM_USE_UMPIRE) && defined(UMPIRE_ENABLE_UM)
    return getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Unified);
#else
    break;
#endif
  case MemorySpace::Pinned:
#if defined(AXOM_USE_UMPIRE) && defined(UMPIRE_ENABLE_PINNED)
    return getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Pinned);
#else
    break;
#endif
  case MemorySpace::Constant:
#if defined(AXOM_USE_UMPIRE) && defined(UMPIRE_ENABLE_CONST)
    return getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Constant);
#else
    break;
#endif
  }

  std::cerr << "Axom memory space \"" << memorySpaceName(space)
            << "\" is not available in this build." << std::endl;
  axom::utilities::processAbort();

  return INVALID_ALLOCATOR_ID;  // Silence warning.
}

void setDefaultAllocator(MemorySpace space)
{
#if defined(AXOM_USE_UMPIRE)
  if(space == MemorySpace::Host)
  {
    setDefaultAllocator(platformHostAllocatorID());
    return;
  }
#endif
  setDefaultAllocator(getAllocatorIDFromMemorySpace(space));
}

void setDefaultHostAllocator(MemorySpace space)
{
  switch(space)
  {
  case MemorySpace::Malloc:
    setDefaultHostAllocator(MALLOC_ALLOCATOR_ID);
    return;
  case MemorySpace::Host:
    setDefaultHostAllocator(platformHostAllocatorID());
    return;
  default:
    break;
  }

  std::cerr << "Axom memory space \"" << memorySpaceName(space)
            << "\" is not a valid default host allocator." << std::endl;
  axom::utilities::processAbort();
}

void setDefaultHostAllocator(int allocId)
{
  if(!isAllocatorCompatibleWithMemorySpace(allocId, MemorySpace::Host))
  {
    std::cerr << "Allocator id " << allocId
              << " is not compatible with Axom's host memory space." << std::endl;
    axom::utilities::processAbort();
  }

  defaultHostAllocatorIDStorage() = allocId;
}

int getDefaultHostAllocatorID()
{
  return defaultHostAllocatorIDStorage();
}

bool isSharedMemoryAllocator(int allocID)
{
  bool isShared = false;
#if defined(AXOM_USE_UMPIRE)
  if(umpire::ResourceManager& rm = umpire::ResourceManager::getInstance(); rm.isAllocator(allocID))
  {
    umpire::Allocator allocator = rm.getAllocator(allocID);

    isShared = allocator.getAllocationStrategy()->getTraits().resource ==
      umpire::MemoryResourceTraits::resource_type::shared;
  }
#else
  AXOM_UNUSED_VAR(allocID);
#endif

  return isShared;
}

int getSharedMemoryAllocatorID(std::size_t minSegmentSize)
{
  int allocator_id = INVALID_ALLOCATOR_ID;
#if defined(AXOM_USE_UMPIRE_SHARED_MEMORY)
  const std::string allocatorName("axom_shared_allocator");
  auto& rm = umpire::ResourceManager::getInstance();
  if(!rm.isAllocator(allocatorName))
  {
    // Create the allocator
    const auto name = axom::fmt::format("SHARED::{}", UMPIRE_DEFAULT_SHARED_MEMORY_RESOURCE);
    auto traits {umpire::get_default_resource_traits(name)};
    traits.scope = umpire::MemoryResourceTraits::shared_scope::node;
    if(minSegmentSize > traits.size)
    {
      traits.size = minSegmentSize;
    }
    auto axom_node_allocator {
      rm.makeResource(axom::fmt::format("{}::axom_node_allocator", name), traits)};
    auto axom_shared_allocator {
      rm.makeAllocator<umpire::strategy::NamedAllocationStrategy>(allocatorName, axom_node_allocator)};
    allocator_id = axom_shared_allocator.getId();
  }
  else
  {
    auto alloc = rm.getAllocator(allocatorName);
    if(minSegmentSize > 0)
    {
      const auto existing_size = alloc.getAllocationStrategy()->getTraits().size;
      if(existing_size > 0 && minSegmentSize > existing_size)
      {
        std::cerr << "Axom shared-memory allocator \"" << allocatorName
                  << "\" already exists with size " << existing_size
                  << " bytes; requested minimum segment size is " << minSegmentSize
                  << " bytes. The shared-memory segment size cannot be increased after creation."
                  << std::endl;
        axom::utilities::processAbort();
      }
    }
    allocator_id = alloc.getId();
  }
#else
  AXOM_UNUSED_VAR(minSegmentSize);
#endif
  return allocator_id;
}

}  // end namespace axom
