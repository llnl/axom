// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/memory_management.hpp"

#if defined(AXOM_USE_UMPIRE_SHARED_MEMORY)
  #include "umpire/Umpire.hpp"
  #include "umpire/strategy/NamedAllocationStrategy.hpp"
  #include "umpire/util/MemoryResourceTraits.hpp"
#endif

namespace axom
{

bool isSharedMemoryAllocator(int allocID)
{
  bool isShared = false;
#if defined(AXOM_USE_UMPIRE)
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  if(rm.isAllocator(allocID))
  {
    umpire::Allocator allocator = rm.getAllocator(allocID);

    isShared = allocator.getAllocationStrategy()->getTraits().resource ==
                 umpire::MemoryResourceTraits::resource_type::shared;
  }
#endif
  return isShared;
}

int getSharedMemoryAllocatorID()
{
  int allocator_id = INVALID_ALLOCATOR_ID;
#if defined(AXOM_USE_UMPIRE_SHARED_MEMORY)
  const std::string allocatorName("axom_named_allocator");
  auto& rm = umpire::ResourceManager::getInstance();
  if(!rm.isAllocator(allocatorName))
  {
    // Create the allocator
    auto traits {umpire::get_default_resource_traits("SHARED")};
    traits.scope = umpire::MemoryResourceTraits::shared_scope::node;
    auto axom_node_allocator {rm.makeResource("SHARED::axom_node_allocator", traits)};
    auto axom_named_allocator {
          rm.makeAllocator<umpire::strategy::NamedAllocationStrategy>(allocatorName,
                                                                      axom_node_allocator)};
    allocator_id = axom_named_allocator.getId();
  }
  else
  {
    allocator_id = rm.getAllocator(allocatorName).getId();
  }
#endif
  return allocator_id;
}

} // end namespace axom
