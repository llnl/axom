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

// Type that indicates that the atomic policy should be selected to match the
// current loop policy where the atomic function is being called.
struct auto_atomic
{ };

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicAdd(T* address, T value)
{
  using atomic_policy =
    typename std::conditional<std::is_same<ExecSpace, auto_atomic>::value,
                              RAJA::auto_atomic,
                              typename axom::execution_space<ExecSpace>::atomic_policy>::type;
  return RAJA::atomicAdd<atomic_policy>(address, value);
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicSub(T* address, T value)
{
  using atomic_policy =
    typename std::conditional<std::is_same<ExecSpace, auto_atomic>::value,
                              RAJA::auto_atomic,
                              typename axom::execution_space<ExecSpace>::atomic_policy>::type;
  return RAJA::atomicSub<atomic_policy>(address, value);
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicMin(T* address, T value)
{
  using atomic_policy =
    typename std::conditional<std::is_same<ExecSpace, auto_atomic>::value,
                              RAJA::auto_atomic,
                              typename axom::execution_space<ExecSpace>::atomic_policy>::type;
  return RAJA::atomicMin<atomic_policy>(address, value);
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicMax(T* address, T value)
{
  using atomic_policy =
    typename std::conditional<std::is_same<ExecSpace, auto_atomic>::value,
                              RAJA::auto_atomic,
                              typename axom::execution_space<ExecSpace>::atomic_policy>::type;
  return RAJA::atomicMax<atomic_policy>(address, value);
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicAnd(T* address, T value)
{
  using atomic_policy =
    typename std::conditional<std::is_same<ExecSpace, auto_atomic>::value,
                              RAJA::auto_atomic,
                              typename axom::execution_space<ExecSpace>::atomic_policy>::type;
  return RAJA::atomicAnd<atomic_policy>(address, value);
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicOr(T* address, T value)
{
  using atomic_policy =
    typename std::conditional<std::is_same<ExecSpace, auto_atomic>::value,
                              RAJA::auto_atomic,
                              typename axom::execution_space<ExecSpace>::atomic_policy>::type;
  return RAJA::atomicOr<atomic_policy>(address, value);
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicXor(T* address, T value)
{
  using atomic_policy =
    typename std::conditional<std::is_same<ExecSpace, auto_atomic>::value,
                              RAJA::auto_atomic,
                              typename axom::execution_space<ExecSpace>::atomic_policy>::type;
  return RAJA::atomicXor<atomic_policy>(address, value);
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicExchange(T* address, T value)
{
  using atomic_policy =
    typename std::conditional<std::is_same<ExecSpace, auto_atomic>::value,
                              RAJA::auto_atomic,
                              typename axom::execution_space<ExecSpace>::atomic_policy>::type;
  return RAJA::atomicExchange<atomic_policy>(address, value);
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicLoad(T* address)
{
  using atomic_policy =
    typename std::conditional<std::is_same<ExecSpace, auto_atomic>::value,
                              RAJA::auto_atomic,
                              typename axom::execution_space<ExecSpace>::atomic_policy>::type;
  return RAJA::atomicLoad<atomic_policy>(address);
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE void atomicStore(T* address, T value)
{
  using atomic_policy =
    typename std::conditional<std::is_same<ExecSpace, auto_atomic>::value,
                              RAJA::auto_atomic,
                              typename axom::execution_space<ExecSpace>::atomic_policy>::type;
  RAJA::atomicStore<atomic_policy>(address, value);
}

}  // namespace axom

#else   // AXOM_HAVE_RAJA
//------------------------------------------------------------------------------
namespace axom
{
// NOTE: There is nothing atomic about these functions but that is okay because
//       they are strictly for serial.

// Type that indicates that the atomic policy should be selected to match the
// current loop policy where the atomic function is being called.
struct auto_atomic
{ };

template <typename ExecSpace>
struct is_serial_atomic_exec
{
  static constexpr bool value =
    std::is_same<ExecSpace, SEQ_EXEC>::value || std::is_same<ExecSpace, auto_atomic>::value;
};

template <typename ExecSpace>
constexpr bool is_serial_atomic_exec_v = is_serial_atomic_exec<ExecSpace>::value;

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicAdd(T* address, T value)
{
  static_assert(is_serial_atomic_exec_v<ExecSpace>);
  const T retval = *address;
  *address += value;
  return retval;
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicSub(T* address, T value)
{
  static_assert(is_serial_atomic_exec_v<ExecSpace>);
  const T retval = *address;
  *address -= value;
  return retval;
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicMin(T* address, T value)
{
  static_assert(is_serial_atomic_exec_v<ExecSpace>);
  const T retval = *address;
  *address = axom::utilities::min(*address, value);
  return retval;
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicMax(T* address, T value)
{
  static_assert(is_serial_atomic_exec_v<ExecSpace>);
  const T retval = *address;
  *address = axom::utilities::max(*address, value);
  return retval;
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicAnd(T* address, T value)
{
  static_assert(is_serial_atomic_exec_v<ExecSpace>);
  const T retval = *address;
  *address &= value;
  return retval;
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicOr(T* address, T value)
{
  static_assert(is_serial_atomic_exec_v<ExecSpace>);
  const T retval = *address;
  *address |= value;
  return retval;
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicXor(T* address, T value)
{
  static_assert(is_serial_atomic_exec_v<ExecSpace>);
  const T retval = *address;
  *address ^= value;
  return retval;
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicExchange(T* address, T value)
{
  static_assert(is_serial_atomic_exec_v<ExecSpace>);
  const T retval = *address;
  *address = value;
  return retval;
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE T atomicLoad(T* address)
{
  static_assert(is_serial_atomic_exec_v<ExecSpace>);
  const T retval = *address;
  return retval;
}

template <typename ExecSpace, typename T>
inline AXOM_HOST_DEVICE void atomicStore(T* address, T value)
{
  static_assert(is_serial_atomic_exec_v<ExecSpace>);
  *address = value;
}

}  // namespace axom
#endif  // AXOM_HAVE_RAJA

#endif
