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
