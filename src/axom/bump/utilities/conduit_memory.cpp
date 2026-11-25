// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/core.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"

namespace axom::bump::utilities
{

bool isDeviceAllocated(const conduit::Node& n)
{
#if defined(AXOM_USE_UMPIRE)
  return isDeviceAllocator(getAllocatorIDFromNode(n));
#else
  AXOM_UNUSED_VAR(n);
  return false;
#endif
}

int getAllocatorIDFromNode(const conduit::Node &n)
{
  // We use Node::data_ptr() here on purpose. The data_ptr() may be a base pointer
  // for multiple nodes that were allocated at the same time.
  return axom::getAllocatorIDFromPointer(n.data_ptr());
}

namespace internal
{

void copyImpl(conduit::Node& dest,
              const conduit::Node& src,
              int destConduitAllocatorID)
{
  if(src.number_of_children() > 0)
  {
    for(conduit::index_t i = 0; i < src.number_of_children(); i++)
    {
      copyImpl(dest[src[i].name()], src[i], destConduitAllocatorID);
    }
  }
  else
  {
    dest.set_allocator(destConduitAllocatorID);
    dest.set(conduit::DataType(src.dtype().id(), src.dtype().number_of_elements()));
    // Copy the data to the destination node. Axom uses Umpire to manage that.
    if(src.is_compact())
      axom::copy(dest.data_ptr(), src.data_ptr(), src.dtype().bytes_compact());
    else
    {
      // NOTE: This assumes that src is on the host.
      conduit::Node tmp;
      src.compact_to(tmp);
      axom::copy(dest.data_ptr(), tmp.data_ptr(), tmp.dtype().bytes_compact());
    }
  }
}

} // end namespace internal
} // end namespace axom::bump::utilities
