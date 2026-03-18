// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_VIEWS_NODE_ARRAY_VIEW_HPP_
#define AXOM_BUMP_VIEWS_NODE_ARRAY_VIEW_HPP_

#include "axom/slic/interface/slic.hpp"

#include <conduit/conduit.hpp>
#include <iostream>

namespace axom
{
namespace bump
{
namespace views
{
namespace detail
{
struct Delimiter
{ };

/// Used to separate arguments.
constexpr Delimiter ArgumentDelimiter;

#if __cplusplus >= 201703L
// C++17 and later.
template <typename... Args>
constexpr int encode_types(Args... args)
{
  return (... | args);
}
#else
template <typename T>
constexpr int encode_types_impl(T arg)
{
  return arg;
}

template <typename T, typename... Args>
constexpr int encode_types_impl(T arg, Args... args)
{
  return (arg | encode_types_impl(args...));
}

template <typename... Args>
constexpr int encode_types(Args... args)
{
  return encode_types_impl(args...);
}
#endif

template <typename... Args>
constexpr int select_types(Args... args)
{
  return encode_types((1 << args)...);
}

constexpr bool type_selected(int flag, int bit) { return flag & (1 << bit); }

constexpr int select_all_types() { return -1; }

constexpr int select_index_types()
{
  return select_types(conduit::DataType::INT32_ID,
                      conduit::DataType::INT64_ID,
                      conduit::DataType::UINT32_ID,
                      conduit::DataType::UINT64_ID);
}

constexpr int select_float_types()
{
  return select_types(conduit::DataType::FLOAT32_ID, conduit::DataType::FLOAT64_ID);
}

//------------------------------------------------------------------------------
// General Node to ArrayView. Handle all types.
//------------------------------------------------------------------------------
/// NOTE: Some of these functions use const_cast to get data into the ArrayView.
///       Is there a better way that does not let const bleed all over?
///
/// TODO: Handle strided data from the Conduit node.

template <bool Enabled, typename FuncType>
std::enable_if_t<Enabled, void> nodeToArrayViewSingleInt8(const conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::int8> view(const_cast<conduit::int8 *>(n.as_int8_ptr()), size);
  func(view);
}

template <bool Enabled, typename FuncType>
std::enable_if_t<!Enabled, void> nodeToArrayViewSingleInt8(const conduit::Node &AXOM_UNUSED_PARAM(n),
                                                               FuncType &&AXOM_UNUSED_PARAM(func))
{
  SLIC_WARNING("Unsupported int8 node.");
}

template <bool Enabled, typename FuncType>
std::enable_if_t<Enabled, void> nodeToArrayViewSingleInt8(conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::int8> view(n.as_int8_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType>
std::enable_if_t<!Enabled, void> nodeToArrayViewSingleInt8(conduit::Node &AXOM_UNUSED_PARAM(n),
                                                               FuncType &&AXOM_UNUSED_PARAM(func))
{
  SLIC_WARNING("Unsupported int8 node.");
}

template <bool Enabled, typename FuncType>
std::enable_if_t<Enabled, void> nodeToArrayViewSingleInt16(const conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::int16> view(const_cast<conduit::int16 *>(n.as_int16_ptr()), size);
  func(view);
}

template <bool Enabled, typename FuncType>
std::enable_if_t<!Enabled, void> nodeToArrayViewSingleInt16(
  const conduit::Node &AXOM_UNUSED_PARAM(n),
  FuncType &&AXOM_UNUSED_PARAM(func))
{
  SLIC_WARNING("Unsupported int16 node.");
}

template <bool Enabled, typename FuncType>
std::enable_if_t<Enabled, void> nodeToArrayViewSingleInt16(conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::int16> view(n.as_int16_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType>
std::enable_if_t<!Enabled, void> nodeToArrayViewSingleInt16(conduit::Node &AXOM_UNUSED_PARAM(n),
                                                                FuncType &&AXOM_UNUSED_PARAM(func))
{
  SLIC_WARNING("Unsupported int16 node.");
}

template <bool Enabled, typename FuncType>
std::enable_if_t<Enabled, void> nodeToArrayViewSingleInt32(const conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::int32> view(const_cast<conduit::int32 *>(n.as_int32_ptr()), size);
  func(view);
}

template <bool Enabled, typename FuncType>
std::enable_if_t<!Enabled, void> nodeToArrayViewSingleInt32(
  const conduit::Node &AXOM_UNUSED_PARAM(n),
  FuncType &&AXOM_UNUSED_PARAM(func))
{
  SLIC_WARNING("Unsupported int32 node.");
}

template <bool Enabled, typename FuncType>
std::enable_if_t<Enabled, void> nodeToArrayViewSingleInt32(conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::int32> view(n.as_int32_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType>
std::enable_if_t<!Enabled, void> nodeToArrayViewSingleInt32(conduit::Node &AXOM_UNUSED_PARAM(n),
                                                                FuncType &&AXOM_UNUSED_PARAM(func))
{
  SLIC_WARNING("Unsupported int32 node.");
}

template <bool Enabled, typename FuncType>
std::enable_if_t<Enabled, void> nodeToArrayViewSingleInt64(const conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::int64> view(const_cast<conduit::int64 *>(n.as_int64_ptr()), size);
  func(view);
}

template <bool Enabled, typename FuncType>
std::enable_if_t<!Enabled, void> nodeToArrayViewSingleInt64(
  const conduit::Node &AXOM_UNUSED_PARAM(n),
  FuncType &&AXOM_UNUSED_PARAM(func))
{
  SLIC_WARNING("Unsupported int64 node.");
}

template <bool Enabled, typename FuncType>
std::enable_if_t<Enabled, void> nodeToArrayViewSingleInt64(conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::int64> view(n.as_int64_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType>
std::enable_if_t<!Enabled, void> nodeToArrayViewSingleInt64(conduit::Node &AXOM_UNUSED_PARAM(n),
                                                                FuncType &&AXOM_UNUSED_PARAM(func))
{
  SLIC_WARNING("Unsupported int64 node.");
}

template <bool Enabled, typename FuncType>
std::enable_if_t<Enabled, void> nodeToArrayViewSingleUint8(const conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::uint8> view(const_cast<conduit::uint8 *>(n.as_uint8_ptr()), size);
  func(view);
}

template <bool Enabled, typename FuncType>
std::enable_if_t<!Enabled, void> nodeToArrayViewSingleUint8(
  const conduit::Node &AXOM_UNUSED_PARAM(n),
  FuncType &&AXOM_UNUSED_PARAM(func))
{
  SLIC_WARNING("Unsupported uint8 node.");
}

template <bool Enabled, typename FuncType>
std::enable_if_t<Enabled, void> nodeToArrayViewSingleUint8(conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::uint8> view(n.as_uint8_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType>
std::enable_if_t<!Enabled, void> nodeToArrayViewSingleUint8(conduit::Node &AXOM_UNUSED_PARAM(n),
                                                                FuncType &&AXOM_UNUSED_PARAM(func))
{
  SLIC_WARNING("Unsupported uint8 node.");
}

template <bool Enabled, typename FuncType>
std::enable_if_t<Enabled, void> nodeToArrayViewSingleUint16(const conduit::Node &n,
                                                                FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::uint16> view(const_cast<conduit::uint16 *>(n.as_uint16_ptr()), size);
  func(view);
}

template <bool Enabled, typename FuncType>
std::enable_if_t<!Enabled, void> nodeToArrayViewSingleUint16(
  const conduit::Node &AXOM_UNUSED_PARAM(n),
  FuncType &&AXOM_UNUSED_PARAM(func))
{
  SLIC_WARNING("Unsupported uint16 node.");
}

template <bool Enabled, typename FuncType>
std::enable_if_t<Enabled, void> nodeToArrayViewSingleUint16(conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::uint16> view(n.as_uint16_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType>
std::enable_if_t<!Enabled, void> nodeToArrayViewSingleUint16(conduit::Node &AXOM_UNUSED_PARAM(n),
                                                                 FuncType &&AXOM_UNUSED_PARAM(func))
{
  SLIC_WARNING("Unsupported uint16 node.");
}

template <bool Enabled, typename FuncType>
std::enable_if_t<Enabled, void> nodeToArrayViewSingleUint32(const conduit::Node &n,
                                                                FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::uint32> view(const_cast<conduit::uint32 *>(n.as_uint32_ptr()), size);
  func(view);
}

template <bool Enabled, typename FuncType>
std::enable_if_t<!Enabled, void> nodeToArrayViewSingleUint32(
  const conduit::Node &AXOM_UNUSED_PARAM(n),
  FuncType &&AXOM_UNUSED_PARAM(func))
{
  SLIC_WARNING("Unsupported uint32 node.");
}

template <bool Enabled, typename FuncType>
std::enable_if_t<Enabled, void> nodeToArrayViewSingleUint32(conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::uint32> view(n.as_uint32_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType>
std::enable_if_t<!Enabled, void> nodeToArrayViewSingleUint32(conduit::Node &AXOM_UNUSED_PARAM(n),
                                                                 FuncType &&AXOM_UNUSED_PARAM(func))
{
  SLIC_WARNING("Unsupported uint32 node.");
}

template <bool Enabled, typename FuncType>
std::enable_if_t<Enabled, void> nodeToArrayViewSingleUint64(const conduit::Node &n,
                                                                FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::uint64> view(const_cast<conduit::uint64 *>(n.as_uint64_ptr()), size);
  func(view);
}

template <bool Enabled, typename FuncType>
std::enable_if_t<!Enabled, void> nodeToArrayViewSingleUint64(
  const conduit::Node &AXOM_UNUSED_PARAM(n),
  FuncType &&AXOM_UNUSED_PARAM(func))
{
  SLIC_WARNING("Unsupported uint64 node.");
}

template <bool Enabled, typename FuncType>
std::enable_if_t<Enabled, void> nodeToArrayViewSingleUint64(conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::uint64> view(n.as_uint64_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType>
std::enable_if_t<!Enabled, void> nodeToArrayViewSingleUint64(conduit::Node &AXOM_UNUSED_PARAM(n),
                                                                 FuncType &&AXOM_UNUSED_PARAM(func))
{
  SLIC_WARNING("Unsupported uint64 node.");
}

template <bool Enabled, typename FuncType>
std::enable_if_t<Enabled, void> nodeToArrayViewSingleFloat32(const conduit::Node &n,
                                                                 FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::float32> view(const_cast<conduit::float32 *>(n.as_float32_ptr()), size);
  func(view);
}

template <bool Enabled, typename FuncType>
std::enable_if_t<!Enabled, void> nodeToArrayViewSingleFloat32(
  const conduit::Node &AXOM_UNUSED_PARAM(n),
  FuncType &&AXOM_UNUSED_PARAM(func))
{
  SLIC_WARNING("Unsupported float32 node.");
}

template <bool Enabled, typename FuncType>
std::enable_if_t<Enabled, void> nodeToArrayViewSingleFloat32(conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::float32> view(n.as_float32_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType>
std::enable_if_t<!Enabled, void> nodeToArrayViewSingleFloat32(conduit::Node &AXOM_UNUSED_PARAM(n),
                                                                  FuncType &&AXOM_UNUSED_PARAM(func))
{
  SLIC_WARNING("Unsupported float32 node.");
}

template <bool Enabled, typename FuncType>
std::enable_if_t<Enabled, void> nodeToArrayViewSingleFloat64(const conduit::Node &n,
                                                                 FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::float64> view(const_cast<conduit::float64 *>(n.as_float64_ptr()), size);
  func(view);
}

template <bool Enabled, typename FuncType>
std::enable_if_t<!Enabled, void> nodeToArrayViewSingleFloat64(
  const conduit::Node &AXOM_UNUSED_PARAM(n),
  FuncType &&AXOM_UNUSED_PARAM(func))
{
  SLIC_WARNING("Unsupported float64 node.");
}

template <bool Enabled, typename FuncType>
std::enable_if_t<Enabled, void> nodeToArrayViewSingleFloat64(conduit::Node &n, FuncType &&func)
{
  const auto size = n.dtype().number_of_elements();
  axom::ArrayView<conduit::float64> view(n.as_float64_ptr(), size);
  func(view);
}

template <bool Enabled, typename FuncType>
std::enable_if_t<!Enabled, void> nodeToArrayViewSingleFloat64(conduit::Node &AXOM_UNUSED_PARAM(n),
                                                                  FuncType &&AXOM_UNUSED_PARAM(func))
{
  SLIC_WARNING("Unsupported float64 node.");
}

template <int Types = select_all_types(), typename FuncType>
void Node_to_ArrayView_single(const conduit::Node &n, FuncType &&func)
{
  /* Later, with C++17, we can do this instead of using all of the SFINAE functions above:
    if constexpr (type_selected(Types, conduit::DataType::INT8_ID))
    {
      if(n.dtype().is_int8())
      {
        axom::ArrayView<conduit::int8> view(n.as_int8_ptr(), size);
        func(view);
      }
    }
  */

  if(n.dtype().is_int8())
  {
    nodeToArrayViewSingleInt8<type_selected(Types, conduit::DataType::INT8_ID)>(n, func);
  }
  else if(n.dtype().is_int16())
  {
    nodeToArrayViewSingleInt16<type_selected(Types, conduit::DataType::INT16_ID)>(n, func);
  }
  else if(n.dtype().is_int32())
  {
    nodeToArrayViewSingleInt32<type_selected(Types, conduit::DataType::INT32_ID)>(n, func);
  }
  else if(n.dtype().is_int64())
  {
    nodeToArrayViewSingleInt64<type_selected(Types, conduit::DataType::INT64_ID)>(n, func);
  }
  else if(n.dtype().is_uint8())
  {
    nodeToArrayViewSingleUint8<type_selected(Types, conduit::DataType::UINT8_ID)>(n, func);
  }
  else if(n.dtype().is_uint16())
  {
    nodeToArrayViewSingleUint16<type_selected(Types, conduit::DataType::UINT16_ID)>(n, func);
  }
  else if(n.dtype().is_uint32())
  {
    nodeToArrayViewSingleUint32<type_selected(Types, conduit::DataType::UINT32_ID)>(n, func);
  }
  else if(n.dtype().is_uint64())
  {
    nodeToArrayViewSingleUint64<type_selected(Types, conduit::DataType::UINT64_ID)>(n, func);
  }
  else if(n.dtype().is_float32())
  {
    nodeToArrayViewSingleFloat32<type_selected(Types, conduit::DataType::FLOAT32_ID)>(n, func);
  }
  else if(n.dtype().is_float64())
  {
    nodeToArrayViewSingleFloat64<type_selected(Types, conduit::DataType::FLOAT64_ID)>(n, func);
  }
  else
  {
    SLIC_ERROR("Unsupported data type " << n.dtype().name() << " on node " << n.path());
  }
}

template <int Types = select_all_types(), typename FuncType>
void Node_to_ArrayView_single(conduit::Node &n, FuncType &&func)
{
  if(n.dtype().is_int8())
  {
    nodeToArrayViewSingleInt8<type_selected(Types, conduit::DataType::INT8_ID)>(n, func);
  }
  else if(n.dtype().is_int16())
  {
    nodeToArrayViewSingleInt16<type_selected(Types, conduit::DataType::INT16_ID)>(n, func);
  }
  else if(n.dtype().is_int32())
  {
    nodeToArrayViewSingleInt32<type_selected(Types, conduit::DataType::INT32_ID)>(n, func);
  }
  else if(n.dtype().is_int64())
  {
    nodeToArrayViewSingleInt64<type_selected(Types, conduit::DataType::INT64_ID)>(n, func);
  }
  else if(n.dtype().is_uint8())
  {
    nodeToArrayViewSingleUint8<type_selected(Types, conduit::DataType::UINT8_ID)>(n, func);
  }
  else if(n.dtype().is_uint16())
  {
    nodeToArrayViewSingleUint16<type_selected(Types, conduit::DataType::UINT16_ID)>(n, func);
  }
  else if(n.dtype().is_uint32())
  {
    nodeToArrayViewSingleUint32<type_selected(Types, conduit::DataType::UINT32_ID)>(n, func);
  }
  else if(n.dtype().is_uint64())
  {
    nodeToArrayViewSingleUint64<type_selected(Types, conduit::DataType::UINT64_ID)>(n, func);
  }
  else if(n.dtype().is_float32())
  {
    nodeToArrayViewSingleFloat32<type_selected(Types, conduit::DataType::FLOAT32_ID)>(n, func);
  }
  else if(n.dtype().is_float64())
  {
    nodeToArrayViewSingleFloat64<type_selected(Types, conduit::DataType::FLOAT64_ID)>(n, func);
  }
  else
  {
    SLIC_ERROR("Unsupported data type " << n.dtype().name() << " on node " << n.path());
  }
}

template <int Types, typename FuncType, typename... View>
void nodeToArrayViewInternal(FuncType &&func, Delimiter, View &...views)
{
  func(views...);
}

template <int Types = select_all_types(), typename... Args>
void nodeToArrayViewInternal(const conduit::Node &first, Args &&...args)
{
  Node_to_ArrayView_single<Types>(first, [&](auto view) {
    nodeToArrayViewInternal<Types>(args..., view);
  });
}

template <int Types = select_all_types(), typename... Args>
void nodeToArrayViewInternal(conduit::Node &first, Args &&...args)
{
  Node_to_ArrayView_single<Types>(first, [&](auto view) {
    nodeToArrayViewInternal<Types>(args..., view);
  });
}

//------------------------------------------------------------------------------
/// NOTE: handle const conduit::Node& better. For now, const_cast.

template <bool Enabled, typename FuncType, typename... Args>
std::enable_if_t<Enabled, void> nodeToArrayViewSameInternalInt8(FuncType &&func, Args &&...args)
{
  func(axom::ArrayView<conduit::int8>(const_cast<conduit::int8 *>(args.as_int8_ptr()),
                                      args.dtype().number_of_elements())...);
}

template <bool Enabled, typename FuncType, typename... Args>
std::enable_if_t<!Enabled, void> nodeToArrayViewSameInternalInt8(
  FuncType &&AXOM_UNUSED_PARAM(func),
  Args &&...AXOM_UNUSED_PARAM(args))
{ }

template <bool Enabled, typename FuncType, typename... Args>
std::enable_if_t<Enabled, void> nodeToArrayViewSameInternalInt16(FuncType &&func, Args &&...args)
{
  func(axom::ArrayView<conduit::int16>(const_cast<conduit::int16 *>(args.as_int16_ptr()),
                                       args.dtype().number_of_elements())...);
}

template <bool Enabled, typename FuncType, typename... Args>
std::enable_if_t<!Enabled, void> nodeToArrayViewSameInternalInt16(
  FuncType &&AXOM_UNUSED_PARAM(func),
  Args &&...AXOM_UNUSED_PARAM(args))
{ }

template <bool Enabled, typename FuncType, typename... Args>
std::enable_if_t<Enabled, void> nodeToArrayViewSameInternalInt32(FuncType &&func, Args &&...args)
{
  func(axom::ArrayView<conduit::int32>(const_cast<conduit::int32 *>(args.as_int32_ptr()),
                                       args.dtype().number_of_elements())...);
}

template <bool Enabled, typename FuncType, typename... Args>
std::enable_if_t<!Enabled, void> nodeToArrayViewSameInternalInt32(
  FuncType &&AXOM_UNUSED_PARAM(func),
  Args &&...AXOM_UNUSED_PARAM(args))
{ }

template <bool Enabled, typename FuncType, typename... Args>
std::enable_if_t<Enabled, void> nodeToArrayViewSameInternalInt64(FuncType &&func, Args &&...args)
{
  func(axom::ArrayView<conduit::int64>(const_cast<conduit::int64 *>(args.as_int64_ptr()),
                                       args.dtype().number_of_elements())...);
}

template <bool Enabled, typename FuncType, typename... Args>
std::enable_if_t<!Enabled, void> nodeToArrayViewSameInternalInt64(
  FuncType &&AXOM_UNUSED_PARAM(func),
  Args &&...AXOM_UNUSED_PARAM(args))
{ }

template <bool Enabled, typename FuncType, typename... Args>
std::enable_if_t<Enabled, void> nodeToArrayViewSameInternalUint8(FuncType &&func, Args &&...args)
{
  func(axom::ArrayView<conduit::uint8>(const_cast<conduit::uint8 *>(args.as_uint8_ptr()),
                                       args.dtype().number_of_elements())...);
}

template <bool Enabled, typename FuncType, typename... Args>
std::enable_if_t<!Enabled, void> nodeToArrayViewSameInternalUint8(
  FuncType &&AXOM_UNUSED_PARAM(func),
  Args &&...AXOM_UNUSED_PARAM(args))
{ }

template <bool Enabled, typename FuncType, typename... Args>
std::enable_if_t<Enabled, void> nodeToArrayViewSameInternalUint16(FuncType &&func, Args &&...args)
{
  func(axom::ArrayView<conduit::uint16>(const_cast<conduit::uint16 *>(args.as_uint16_ptr()),
                                        args.dtype().number_of_elements())...);
}

template <bool Enabled, typename FuncType, typename... Args>
std::enable_if_t<!Enabled, void> nodeToArrayViewSameInternalUint16(
  FuncType &&AXOM_UNUSED_PARAM(func),
  Args &&...AXOM_UNUSED_PARAM(args))
{ }

template <bool Enabled, typename FuncType, typename... Args>
std::enable_if_t<Enabled, void> nodeToArrayViewSameInternalUint32(FuncType &&func, Args &&...args)
{
  func(axom::ArrayView<conduit::uint32>(const_cast<conduit::uint32 *>(args.as_uint32_ptr()),
                                        args.dtype().number_of_elements())...);
}

template <bool Enabled, typename FuncType, typename... Args>
std::enable_if_t<!Enabled, void> nodeToArrayViewSameInternalUint32(
  FuncType &&AXOM_UNUSED_PARAM(func),
  Args &&...AXOM_UNUSED_PARAM(args))
{ }

template <bool Enabled, typename FuncType, typename... Args>
std::enable_if_t<Enabled, void> nodeToArrayViewSameInternalUint64(FuncType &&func, Args &&...args)
{
  func(axom::ArrayView<conduit::uint64>(const_cast<conduit::uint64 *>(args.as_uint64_ptr()),
                                        args.dtype().number_of_elements())...);
}

template <bool Enabled, typename FuncType, typename... Args>
std::enable_if_t<!Enabled, void> nodeToArrayViewSameInternalUint64(
  FuncType &&AXOM_UNUSED_PARAM(func),
  Args &&...AXOM_UNUSED_PARAM(args))
{ }

template <bool Enabled, typename FuncType, typename... Args>
std::enable_if_t<Enabled, void> nodeToArrayViewSameInternalFloat32(FuncType &&func,
                                                                        Args &&...args)
{
  func(axom::ArrayView<conduit::float32>(const_cast<conduit::float32 *>(args.as_float32_ptr()),
                                         args.dtype().number_of_elements())...);
}

template <bool Enabled, typename FuncType, typename... Args>
std::enable_if_t<!Enabled, void> nodeToArrayViewSameInternalFloat32(
  FuncType &&AXOM_UNUSED_PARAM(func),
  Args &&...AXOM_UNUSED_PARAM(args))
{ }

template <bool Enabled, typename FuncType, typename... Args>
std::enable_if_t<Enabled, void> nodeToArrayViewSameInternalFloat64(FuncType &&func,
                                                                        Args &&...args)
{
  func(axom::ArrayView<conduit::float64>(const_cast<conduit::float64 *>(args.as_float64_ptr()),
                                         args.dtype().number_of_elements())...);
}

template <bool Enabled, typename FuncType, typename... Args>
std::enable_if_t<!Enabled, void> nodeToArrayViewSameInternalFloat64(
  FuncType &&AXOM_UNUSED_PARAM(func),
  Args &&...AXOM_UNUSED_PARAM(args))
{ }

template <int Types = select_all_types(), typename FuncType, typename... Args>
void nodeToArrayViewSameInternal(FuncType &&func,
                                     Delimiter,
                                     const conduit::Node &first,
                                     Args &&...args)
{
  if(first.dtype().is_int8())
  {
    nodeToArrayViewSameInternalInt8<type_selected(Types, conduit::DataType::INT8_ID)>(func,
                                                                                           first,
                                                                                           args...);
  }
  else if(first.dtype().is_int16())
  {
    nodeToArrayViewSameInternalInt16<type_selected(Types, conduit::DataType::INT16_ID)>(func,
                                                                                             first,
                                                                                             args...);
  }
  else if(first.dtype().is_int32())
  {
    nodeToArrayViewSameInternalInt32<type_selected(Types, conduit::DataType::INT32_ID)>(func,
                                                                                             first,
                                                                                             args...);
  }
  else if(first.dtype().is_int64())
  {
    nodeToArrayViewSameInternalInt64<type_selected(Types, conduit::DataType::INT64_ID)>(func,
                                                                                             first,
                                                                                             args...);
  }
  else if(first.dtype().is_uint8())
  {
    nodeToArrayViewSameInternalUint8<type_selected(Types, conduit::DataType::UINT8_ID)>(func,
                                                                                             first,
                                                                                             args...);
  }
  else if(first.dtype().is_uint16())
  {
    nodeToArrayViewSameInternalUint16<type_selected(Types, conduit::DataType::UINT16_ID)>(
      func,
      first,
      args...);
  }
  else if(first.dtype().is_uint32())
  {
    nodeToArrayViewSameInternalUint32<type_selected(Types, conduit::DataType::UINT32_ID)>(
      func,
      first,
      args...);
  }
  else if(first.dtype().is_uint64())
  {
    nodeToArrayViewSameInternalUint64<type_selected(Types, conduit::DataType::UINT64_ID)>(
      func,
      first,
      args...);
  }
  else if(first.dtype().is_float32())
  {
    nodeToArrayViewSameInternalFloat32<type_selected(Types, conduit::DataType::FLOAT32_ID)>(
      func,
      first,
      args...);
  }
  else if(first.dtype().is_float64())
  {
    nodeToArrayViewSameInternalFloat64<type_selected(Types, conduit::DataType::FLOAT64_ID)>(
      func,
      first,
      args...);
  }
  else
  {
    SLIC_ERROR("Unsupported data type " << first.dtype().name() << " on node " << first.path());
  }
}

template <int Types = select_all_types(), typename FuncType, typename... Args>
void nodeToArrayViewSameInternal(FuncType &&func, Delimiter, conduit::Node &first, Args &&...args)
{
  if(first.dtype().is_int8())
  {
    nodeToArrayViewSameInternalInt8<type_selected(Types, conduit::DataType::INT8_ID)>(func,
                                                                                           first,
                                                                                           args...);
  }
  else if(first.dtype().is_int16())
  {
    nodeToArrayViewSameInternalInt16<type_selected(Types, conduit::DataType::INT16_ID)>(func,
                                                                                             first,
                                                                                             args...);
  }
  else if(first.dtype().is_int32())
  {
    nodeToArrayViewSameInternalInt32<type_selected(Types, conduit::DataType::INT32_ID)>(func,
                                                                                             first,
                                                                                             args...);
  }
  else if(first.dtype().is_int64())
  {
    nodeToArrayViewSameInternalInt64<type_selected(Types, conduit::DataType::INT64_ID)>(func,
                                                                                             first,
                                                                                             args...);
  }
  else if(first.dtype().is_uint8())
  {
    nodeToArrayViewSameInternalUint8<type_selected(Types, conduit::DataType::UINT8_ID)>(func,
                                                                                             first,
                                                                                             args...);
  }
  else if(first.dtype().is_uint16())
  {
    nodeToArrayViewSameInternalUint16<type_selected(Types, conduit::DataType::UINT16_ID)>(
      func,
      first,
      args...);
  }
  else if(first.dtype().is_uint32())
  {
    nodeToArrayViewSameInternalUint32<type_selected(Types, conduit::DataType::UINT32_ID)>(
      func,
      first,
      args...);
  }
  else if(first.dtype().is_uint64())
  {
    nodeToArrayViewSameInternalUint64<type_selected(Types, conduit::DataType::UINT64_ID)>(
      func,
      first,
      args...);
  }
  else if(first.dtype().is_float32())
  {
    nodeToArrayViewSameInternalFloat32<type_selected(Types, conduit::DataType::FLOAT32_ID)>(
      func,
      first,
      args...);
  }
  else if(first.dtype().is_float64())
  {
    nodeToArrayViewSameInternalFloat64<type_selected(Types, conduit::DataType::FLOAT64_ID)>(
      func,
      first,
      args...);
  }
  else
  {
    SLIC_ERROR("Unsupported data type " << first.dtype().name() << " on node " << first.path());
  }
}

/// Reorder args
template <int Types = select_all_types(), typename... Args>
void nodeToArrayViewSameInternal(const conduit::Node &first, Args &&...args)
{
  nodeToArrayViewSameInternal<Types>(args..., first);
}

template <int Types = select_all_types(), typename... Args>
void nodeToArrayViewSameInternal(conduit::Node &first, Args &&...args)
{
  nodeToArrayViewSameInternal<Types>(args..., first);
}

}  // namespace detail

//------------------------------------------------------------------------------
// Node to ArrayView. Handle all types.
//------------------------------------------------------------------------------

/*!
 * \brief Convert a series of Conduit nodes to axom::ArrayView and pass the concrete
 *        views to a lambda function passed as the last argument.
 *
 * \note  This method handles all Conduit array types and will instantiate views
 *        of any type. In other words, mixed node types can be used.
 *
 * \param first A Conduit node to be convered to a view.
 * \param args  A sequence of Conduit nodes, followed by a lambda that can accept
 *              the same number of views of any type.
 *
 * nodeToArrayView(node1, node2, [](auto &view1, auto &view2) { });
 * 
 */
template <typename... Args>
void nodeToArrayView(const conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewInternal(first, args..., detail::ArgumentDelimiter);
}

template <typename... Args>
void nodeToArrayView(conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewInternal(first, args..., detail::ArgumentDelimiter);
}

/*!
 * \brief Convert a series of Conduit nodes to axom::ArrayView and pass the concrete
 *        views to a lambda function passed as the last argument.
 *
 * \note  This method handles all Conduit array types. All nodes will be treated
 *        as the same type as the first Conduit node. Use this when all nodes
 *        are assumed to contain the same type since it results in less code.
 *
 * \param first A Conduit node to be convered to a view.
 * \param args  A sequence of Conduit nodes, followed by a lambda that can accept
 *              the same number of views of any type.
 *
 * nodeToArrayViewSame(node1, node2, [](auto &view1, auto &view2) { });
 * 
 */
template <typename... Args>
void nodeToArrayViewSame(const conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewSameInternal(first, args..., detail::ArgumentDelimiter);
}

template <typename... Args>
void nodeToArrayViewSame(conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewSameInternal(first, args..., detail::ArgumentDelimiter);
}

//------------------------------------------------------------------------------
// Index Node to ArrayView. Handle types used for indexing.
//------------------------------------------------------------------------------

template <typename... Args>
void indexNodeToArrayView(const conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewInternal<detail::select_index_types()>(first,
                                                                   args...,
                                                                   detail::ArgumentDelimiter);
}

template <typename... Args>
void indexNodeToArrayView(conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewInternal<detail::select_index_types()>(first,
                                                                   args...,
                                                                   detail::ArgumentDelimiter);
}

template <typename... Args>
void indexNodeToArrayViewSame(const conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewSameInternal<detail::select_index_types()>(first,
                                                                        args...,
                                                                        detail::ArgumentDelimiter);
}

template <typename... Args>
void indexNodeToArrayViewSame(conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewSameInternal<detail::select_index_types()>(first,
                                                                        args...,
                                                                        detail::ArgumentDelimiter);
}

//------------------------------------------------------------------------------
// Float Node to ArrayView. Handle float types.
//------------------------------------------------------------------------------
template <typename... Args>
void floatNodeToArrayView(const conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewInternal<detail::select_float_types()>(first,
                                                                   args...,
                                                                   detail::ArgumentDelimiter);
}

template <typename... Args>
void floatNodeToArrayView(conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewInternal<detail::select_float_types()>(first,
                                                                   args...,
                                                                   detail::ArgumentDelimiter);
}

template <typename... Args>
void floatNodeToArrayViewSame(const conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewSameInternal<detail::select_float_types()>(first,
                                                                        args...,
                                                                        detail::ArgumentDelimiter);
}

template <typename... Args>
void floatNodeToArrayViewSame(conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewSameInternal<detail::select_float_types()>(first,
                                                                        args...,
                                                                        detail::ArgumentDelimiter);
}

}  // namespace views
}  // namespace bump
}  // namespace axom

#endif
