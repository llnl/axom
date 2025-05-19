// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_CORE_EXECUTION_ATOMICS_HPP_
#define AXOM_CORE_EXECUTION_ATOMICS_HPP_

#include "axom/config.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"
#include "axom/core/utilities/Utilities.hpp"

#if defined(AXOM_USE_RAJA)
  //------------------------------------------------------------------------------
  #include "RAJA/RAJA.hpp"

namespace axom
{
// Some RAJA atomic operations.

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicAdd(T* address, T value)
{
  using atomic_policy = typename axom::execution_space<ExecSpace>::atomic_policy;
  return RAJA::atomicAdd<atomic_policy>(address, value);
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicSub(T* address, T value)
{
  using atomic_policy = typename axom::execution_space<ExecSpace>::atomic_policy;
  return RAJA::atomicSub<atomic_policy>(address, value);
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicMin(T* address, T value)
{
  using atomic_policy = typename axom::execution_space<ExecSpace>::atomic_policy;
  return RAJA::atomicMin<atomic_policy>(address, value);
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicMax(T* address, T value)
{
  using atomic_policy = typename axom::execution_space<ExecSpace>::atomic_policy;
  return RAJA::atomicMax<atomic_policy>(address, value);
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicAnd(T* address, T value)
{
  using atomic_policy = typename axom::execution_space<ExecSpace>::atomic_policy;
  return RAJA::atomicAnd<atomic_policy>(address, value);
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicOr(T* address, T value)
{
  using atomic_policy = typename axom::execution_space<ExecSpace>::atomic_policy;
  return RAJA::atomicOr<atomic_policy>(address, value);
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicXor(T* address, T value)
{
  using atomic_policy = typename axom::execution_space<ExecSpace>::atomic_policy;
  return RAJA::atomicXor<atomic_policy>(address, value);
}

}  // namespace axom

#else
//------------------------------------------------------------------------------
namespace axom
{
// NOTE: There is nothing atomic about these functions but that is okay because
//       they are strictly for serial.

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicAdd(T* address, T value)
{
  AXOM_STATIC_ASSERT(std::is_same<ExecSpace, SEQ_EXEC>::value);
  *address += value;
  return *address;
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicSub(T* address, T value)
{
  AXOM_STATIC_ASSERT(std::is_same<ExecSpace, SEQ_EXEC>::value);
  *address -= value;
  return *address;
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicMin(T* address, T value)
{
  AXOM_STATIC_ASSERT(std::is_same<ExecSpace, SEQ_EXEC>::value);
  *address = axom::utilities::min(*address, value);
  return *address;
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicMax(T* address, T value)
{
  AXOM_STATIC_ASSERT(std::is_same<ExecSpace, SEQ_EXEC>::value);
  *address = axom::utilities::max(*address, value);
  return *address;
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicAnd(T* address, T value)
{
  AXOM_STATIC_ASSERT(std::is_same<ExecSpace, SEQ_EXEC>::value);
  *address &= value;
  return *address;
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicOr(T* address, T value)
{
  AXOM_STATIC_ASSERT(std::is_same<ExecSpace, SEQ_EXEC>::value);
  *address |= value;
  return *address;
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicXor(T* address, T value)
{
  AXOM_STATIC_ASSERT(std::is_same<ExecSpace, SEQ_EXEC>::value);
  *address ^= value;
  return *address;
}

}  // namespace axom
#endif  // AXOM_HAVE_RAJA

#endif
