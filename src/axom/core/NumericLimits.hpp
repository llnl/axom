// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *
 * \file NumericLimits.hpp
 *
 * \brief Header file containing portability layer for std::numeric_limits
 *        capabilities
 *
 */

#ifndef AXOM_NUMERICLIMITS_HPP_
#define AXOM_NUMERICLIMITS_HPP_

#include "axom/config.hpp"  // for compile-time definitions

#include <limits>

#if defined(AXOM_USE_CUDA)
  #include <cuda/std/limits>
#endif

namespace axom
{
#if defined(AXOM_USE_CUDA)
// Note: cuda::std types work in host and device code as long as Axom is
//       configured with CUDA enabled. No need to rely on two different
//       header files in that case.
template <typename T>
using numeric_limits = cuda::std::numeric_limits<T>;
#else
template <typename T>
using numeric_limits = std::numeric_limits<T>;
#endif

}  // namespace axom

#endif  // AXOM_NUMERICLIMITS_HPP_
