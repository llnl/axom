// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_VIEWS_NODE_ARRAY_VIEW_HPP_
#define AXOM_BUMP_VIEWS_NODE_ARRAY_VIEW_HPP_

#include "axom/bump/utilities/conduit_array_view.hpp"
#include "axom/slic/interface/slic.hpp"

#include <conduit/conduit.hpp>

#include <type_traits>
#include <utility>

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

// NOTE: Const Conduit nodes still expose mutable ArrayViews for compatibility.

#define AXOM_BUMP_NODE_ARRAY_VIEW_TYPES(MACRO)                               \
  MACRO(conduit::int8, INT8_ID, is_int8, as_int8_ptr, "int8")                \
  MACRO(conduit::int16, INT16_ID, is_int16, as_int16_ptr, "int16")           \
  MACRO(conduit::int32, INT32_ID, is_int32, as_int32_ptr, "int32")           \
  MACRO(conduit::int64, INT64_ID, is_int64, as_int64_ptr, "int64")           \
  MACRO(conduit::uint8, UINT8_ID, is_uint8, as_uint8_ptr, "uint8")           \
  MACRO(conduit::uint16, UINT16_ID, is_uint16, as_uint16_ptr, "uint16")      \
  MACRO(conduit::uint32, UINT32_ID, is_uint32, as_uint32_ptr, "uint32")      \
  MACRO(conduit::uint64, UINT64_ID, is_uint64, as_uint64_ptr, "uint64")      \
  MACRO(conduit::float32, FLOAT32_ID, is_float32, as_float32_ptr, "float32") \
  MACRO(conduit::float64, FLOAT64_ID, is_float64, as_float64_ptr, "float64")

template <typename T>
struct NodeTypeTraits;

#define AXOM_BUMP_DECLARE_NODE_TYPE_TRAITS(CppType, DTypeID, IsMethod, PtrMethod, Label)          \
  template <>                                                                                     \
  struct NodeTypeTraits<CppType>                                                                  \
  {                                                                                               \
    static bool matches(const conduit::Node &n) { return n.dtype().IsMethod(); }                  \
    static CppType *data(conduit::Node &n) { return n.PtrMethod(); }                              \
    static CppType *data(const conduit::Node &n) { return const_cast<CppType *>(n.PtrMethod()); } \
    static const char *label() { return Label; }                                                  \
    static int dtypeId() { return conduit::DataType::DTypeID; }                                   \
  };

AXOM_BUMP_NODE_ARRAY_VIEW_TYPES(AXOM_BUMP_DECLARE_NODE_TYPE_TRAITS)

#undef AXOM_BUMP_DECLARE_NODE_TYPE_TRAITS

template <bool Enabled, typename T, typename FuncType, typename NodeType>
std::enable_if_t<Enabled, void> invoke_single_array_view(NodeType &n, FuncType &&func)
{
  func(axom::bump::utilities::detail::make_conduit_array_view<T>(n));
}

template <bool Enabled, typename T, typename FuncType, typename NodeType>
std::enable_if_t<!Enabled, void> invoke_single_array_view(NodeType &AXOM_UNUSED_PARAM(n),
                                                          FuncType &&AXOM_UNUSED_PARAM(func))
{
  SLIC_WARNING("Unsupported " << NodeTypeTraits<T>::label() << " node.");
}

template <bool Enabled, typename T, typename FuncType, typename... NodeTypes>
std::enable_if_t<Enabled, void> invoke_same_array_views(FuncType &&func, NodeTypes &&...nodes)
{
  func(axom::bump::utilities::detail::make_conduit_array_view<T>(nodes)...);
}

template <bool Enabled, typename T, typename FuncType, typename... NodeTypes>
std::enable_if_t<!Enabled, void> invoke_same_array_views(FuncType &&AXOM_UNUSED_PARAM(func),
                                                         NodeTypes &&...AXOM_UNUSED_PARAM(nodes))
{ }

template <int Types = select_all_types(), typename NodeType, typename FuncType>
void dispatch_single_array_view(NodeType &n, FuncType &&func)
{
#define AXOM_BUMP_DISPATCH_SINGLE(CppType, DTypeID, IsMethod, PtrMethod, Label)          \
  if(NodeTypeTraits<CppType>::matches(n))                                                \
  {                                                                                      \
    invoke_single_array_view<type_selected(Types, conduit::DataType::DTypeID), CppType>( \
      n,                                                                                 \
      std::forward<FuncType>(func));                                                     \
  }                                                                                      \
  else

  AXOM_BUMP_NODE_ARRAY_VIEW_TYPES(AXOM_BUMP_DISPATCH_SINGLE)
  {
    SLIC_ERROR("Unsupported data type " << n.dtype().name() << " on node " << n.path());
  }

#undef AXOM_BUMP_DISPATCH_SINGLE
}

template <int Types = select_all_types(), typename FirstNodeType, typename FuncType, typename... NodeTypes>
void dispatch_same_array_views(FirstNodeType &first, FuncType &&func, NodeTypes &&...nodes)
{
#define AXOM_BUMP_DISPATCH_SAME(CppType, DTypeID, IsMethod, PtrMethod, Label)           \
  if(NodeTypeTraits<CppType>::matches(first))                                           \
  {                                                                                     \
    invoke_same_array_views<type_selected(Types, conduit::DataType::DTypeID), CppType>( \
      std::forward<FuncType>(func),                                                     \
      std::forward<NodeTypes>(nodes)...);                                               \
  }                                                                                     \
  else

  AXOM_BUMP_NODE_ARRAY_VIEW_TYPES(AXOM_BUMP_DISPATCH_SAME)
  {
    SLIC_ERROR("Unsupported data type " << first.dtype().name() << " on node " << first.path());
  }

#undef AXOM_BUMP_DISPATCH_SAME
}

template <int Types, typename FuncType, typename... ViewTypes>
void nodeToArrayViewInternal(FuncType &&func, Delimiter, ViewTypes... views)
{
  func(views...);
}

template <int Types = select_all_types(), typename NodeType, typename... Args>
void nodeToArrayViewInternal(NodeType &first, Args &&...args)
{
  dispatch_single_array_view<Types>(first, [&](auto view) {
    nodeToArrayViewInternal<Types>(std::forward<Args>(args)..., view);
  });
}

template <int Types = select_all_types(), typename FuncType, typename FirstNode, typename... Args>
void nodeToArrayViewSameInternal(FuncType &&func, Delimiter, FirstNode &first, Args &&...args)
{
  dispatch_same_array_views<Types>(first,
                                   std::forward<FuncType>(func),
                                   first,
                                   std::forward<Args>(args)...);
}

template <int Types = select_all_types(), typename FirstNode, typename... Args>
void nodeToArrayViewSameInternal(FirstNode &first, Args &&...args)
{
  nodeToArrayViewSameInternal<Types>(std::forward<Args>(args)..., first);
}

#undef AXOM_BUMP_NODE_ARRAY_VIEW_TYPES

}  // namespace detail

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
 */
template <typename... Args>
void nodeToArrayView(const conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewInternal(first, std::forward<Args>(args)..., detail::ArgumentDelimiter);
}

template <typename... Args>
void nodeToArrayView(conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewInternal(first, std::forward<Args>(args)..., detail::ArgumentDelimiter);
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
 */
template <typename... Args>
void nodeToArrayViewSame(const conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewSameInternal(first, std::forward<Args>(args)..., detail::ArgumentDelimiter);
}

template <typename... Args>
void nodeToArrayViewSame(conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewSameInternal(first, std::forward<Args>(args)..., detail::ArgumentDelimiter);
}

template <typename... Args>
void indexNodeToArrayView(const conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewInternal<detail::select_index_types()>(first,
                                                                std::forward<Args>(args)...,
                                                                detail::ArgumentDelimiter);
}

template <typename... Args>
void indexNodeToArrayView(conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewInternal<detail::select_index_types()>(first,
                                                                std::forward<Args>(args)...,
                                                                detail::ArgumentDelimiter);
}

template <typename... Args>
void indexNodeToArrayViewSame(const conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewSameInternal<detail::select_index_types()>(first,
                                                                    std::forward<Args>(args)...,
                                                                    detail::ArgumentDelimiter);
}

template <typename... Args>
void indexNodeToArrayViewSame(conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewSameInternal<detail::select_index_types()>(first,
                                                                    std::forward<Args>(args)...,
                                                                    detail::ArgumentDelimiter);
}

template <typename... Args>
void floatNodeToArrayView(const conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewInternal<detail::select_float_types()>(first,
                                                                std::forward<Args>(args)...,
                                                                detail::ArgumentDelimiter);
}

template <typename... Args>
void floatNodeToArrayView(conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewInternal<detail::select_float_types()>(first,
                                                                std::forward<Args>(args)...,
                                                                detail::ArgumentDelimiter);
}

template <typename... Args>
void floatNodeToArrayViewSame(const conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewSameInternal<detail::select_float_types()>(first,
                                                                    std::forward<Args>(args)...,
                                                                    detail::ArgumentDelimiter);
}

template <typename... Args>
void floatNodeToArrayViewSame(conduit::Node &first, Args &&...args)
{
  detail::nodeToArrayViewSameInternal<detail::select_float_types()>(first,
                                                                    std::forward<Args>(args)...,
                                                                    detail::ArgumentDelimiter);
}

}  // namespace views
}  // namespace bump
}  // namespace axom

#endif
