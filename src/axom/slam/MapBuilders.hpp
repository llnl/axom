// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file MapBuilders.hpp
 *
 * \brief Construction helpers for SLAM maps that deduce a full policy stack
 *        from a set, stride, and backing buffer.
 */

#pragma once

#include "axom/core/ArrayView.hpp"
#include "axom/slic.hpp"

#include "axom/slam/Concepts.hpp"
#include "axom/slam/Map.hpp"
#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/policies/StridePolicies.hpp"

namespace axom::slam
{
namespace detail
{
template <typename SetType, typename PosType>
axom::IndexType map_storage_size(const SetType* set, PosType stride)
{
  const axom::IndexType sz = set ? static_cast<axom::IndexType>(set->size()) : axom::IndexType {0};
  return sz * static_cast<axom::IndexType>(stride);
}

/*!
 * \brief Debug-only check that an ArrayView backing a map is correctly sized.
 *
 * An ArrayView-backed map indexes through `pos * stride + offset`, so the backing
 * storage must contain exactly `set->size() * stride` elements. Unlike the raw-pointer
 * make_map overloads (which size the view themselves), the ArrayView overloads trust the
 * caller's view length; an undersized view would index out of bounds. This asserts the
 * invariant in debug builds and is a no-op in release builds.
 */
template <typename SetType, typename T, typename PosType>
inline void check_map_view_size(const SetType* set,
                                PosType stride,
                                const axom::ArrayView<T>& AXOM_DEBUG_PARAM(data))
{
#ifdef AXOM_DEBUG
  const axom::IndexType expected = map_storage_size(set, stride);
  SLIC_ASSERT_MSG(static_cast<axom::IndexType>(data.size()) == expected,
                  "slam::make_map -- ArrayView backing storage has "
                    << data.size() << " elements, but the set (size "
                    << (set ? static_cast<axom::IndexType>(set->size()) : axom::IndexType {0})
                    << ") with stride " << stride << " requires exactly " << expected << ".");
#else
  AXOM_UNUSED_VAR(set);
  AXOM_UNUSED_VAR(stride);
#endif
}
}  // namespace detail

/// \name Map construction helpers
/// \brief Construct a SLAM map while deducing its policy stack from the set and
/// backing buffer. Runtime strides must model PositionLike and be convertible
/// to the set's position type; the returned map always uses that position type.
/// Compile-time strides must be positive.
/// \{

/*!
 * \brief Make a strided SLAM map backed by ArrayView storage.
 *
 * \param set   pointer to the map's set (must outlive the map)
 * \param stride runtime stride (#values per set element)
 * \param data  backing storage as an ArrayView (must outlive the map)
 *
 * \pre `data.size() == set->size() * stride`. The view must be sized to back every
 *  element of the set at the given stride; this is checked in debug builds.
 */
template <typename SetType, typename T, typename StrideType>
  requires detail::SetPositionConvertible<SetType, StrideType>
auto make_map(const SetType* set, StrideType stride, axom::ArrayView<T> data)
{
  using PosType = typename SetType::PositionType;
  using Indirection = policies::ArrayViewIndirection<PosType, T>;
  using Stride = policies::RuntimeStride<PosType>;
  using MapType = Map<T, SetType, Indirection, Stride>;
  const auto canonicalStride = static_cast<PosType>(stride);
  detail::check_map_view_size(set, canonicalStride, data);
  return MapType(set, data, canonicalStride);
}

/*!
 * \brief Make a stride-one SLAM map backed by ArrayView storage.
 *
 * \pre `data.size() == set->size()`. The view must be sized to back every element of the
 *  set; this is checked in debug builds.
 */
template <typename SetType, typename T, typename ExplicitPosType = void>
  requires detail::OptionalSetPositionSame<SetType, ExplicitPosType>
auto make_map(const SetType* set, axom::ArrayView<T> data)
{
  using PosType = typename SetType::PositionType;
  using Indirection = policies::ArrayViewIndirection<PosType, T>;
  using Stride = policies::StrideOne<PosType>;
  using MapType = Map<T, SetType, Indirection, Stride>;
  detail::check_map_view_size(set, PosType {1}, data);
  return MapType(set, data);
}

/*!
 * \brief Make a strided SLAM map backed by a raw pointer buffer.
 *
 * This overload wraps the buffer as an ArrayView with length `set->size() * stride`
 * and returns an ArrayView-backed map.
 */
template <typename SetType, typename T, typename StrideType>
  requires detail::SetPositionConvertible<SetType, StrideType>
auto make_map(const SetType* set, StrideType stride, T* data)
{
  using PosType = typename SetType::PositionType;
  const auto canonicalStride = static_cast<PosType>(stride);
  const auto n = detail::map_storage_size(set, canonicalStride);
  return make_map(set, canonicalStride, axom::ArrayView<T>(data, n));
}

/*!
 * \brief Make a stride-one SLAM map backed by a raw pointer buffer.
 */
template <typename SetType, typename T, typename ExplicitPosType = void>
  requires detail::OptionalSetPositionSame<SetType, ExplicitPosType>
auto make_map(const SetType* set, T* data)
{
  using PosType = typename SetType::PositionType;
  const auto n = detail::map_storage_size(set, PosType {1});
  return make_map(set, axom::ArrayView<T>(data, n));
}

/*!
 * \brief Make a compile-time strided SLAM map backed by ArrayView storage.
 *
 * \tparam STRIDE number of values per set element
 *
 * \pre `data.size() == set->size() * STRIDE`. The view must be sized to back every element
 *  of the set at the compile-time stride; this is checked in debug builds.
 */
template <int STRIDE, typename SetType, typename T, typename ExplicitPosType = void>
  requires detail::PositiveStaticStrideFor<STRIDE, SetType> &&
  detail::OptionalSetPositionSame<SetType, ExplicitPosType>
auto make_map_ct(const SetType* set, axom::ArrayView<T> data)
{
  using PosType = typename SetType::PositionType;
  using Indirection = policies::ArrayViewIndirection<PosType, T>;
  using Stride = policies::CompileTimeStride<PosType, static_cast<PosType>(STRIDE)>;
  using MapType = Map<T, SetType, Indirection, Stride>;
  detail::check_map_view_size(set, static_cast<PosType>(STRIDE), data);
  return MapType(set, data);
}

/*!
 * \brief Make a compile-time strided SLAM map backed by a raw pointer buffer.
 *
 * This overload wraps the buffer as an ArrayView with length `set->size() * STRIDE`
 * and returns an ArrayView-backed map.
 */
template <int STRIDE, typename SetType, typename T, typename ExplicitPosType = void>
  requires detail::PositiveStaticStrideFor<STRIDE, SetType> &&
  detail::OptionalSetPositionSame<SetType, ExplicitPosType>
auto make_map_ct(const SetType* set, T* data)
{
  using PosType = typename SetType::PositionType;
  const PosType stride = static_cast<PosType>(STRIDE);
  const auto n = detail::map_storage_size(set, stride);
  return make_map_ct<STRIDE>(set, axom::ArrayView<T>(data, n));
}

/// \}

}  // namespace axom::slam
