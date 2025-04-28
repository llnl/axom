// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_RADIUSS_TUTORIAL_HIP_PATCH
#define AXOM_RADIUSS_TUTORIAL_HIP_PATCH

#include "axom/config.hpp"

#ifdef AXOM_USE_RAJA
  #include "RAJA/RAJA.hpp"
#endif

/// Need to patch rocprim with some versions of hip when using C++14
#if defined(AXOM_USE_HIP) && __cplusplus == 201402L && HIP_VERSION >= 50600000 && \
  HIP_VERSION < 60100000

// Error is due to missing an out-of-class definition for `block_radix_sort<...>::radix_bits_per_pass`.
namespace rocprim
{
template <class Key,
          unsigned int BlockSizeX,
          unsigned int ItemsPerThread,
          class Value,
          unsigned int BlockSizeY,
          unsigned int BlockSizeZ>
constexpr unsigned int
  block_radix_sort<Key, BlockSizeX, ItemsPerThread, Value, BlockSizeY, BlockSizeZ>::radix_bits_per_pass;
}  // namespace rocprim

#endif

#endif  // AXOM_RADIUSS_TUTORIAL_HIP_PATCH