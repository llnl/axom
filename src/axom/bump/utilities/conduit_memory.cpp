// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/core.hpp"
#include "axom/bump/utilities/conduit_memory.hpp"

namespace axom::bump::utilities
{

bool isDeviceAllocated(const conduit::Node &n)
{
#if defined(AXOM_USE_UMPIRE)
  return isDeviceAllocator(axom::getAllocatorIDFromPointer(n.data_ptr()));
#else
  AXOM_UNUSED_VAR(n);
  return false;
#endif
}

}  // end namespace axom::bump::utilities
