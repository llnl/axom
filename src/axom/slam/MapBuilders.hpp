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

#ifndef SLAM_MAP_BUILDERS_HPP_
#define SLAM_MAP_BUILDERS_HPP_

#include "axom/core/ArrayView.hpp"

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
}  // namespace detail

/*!
 * \brief Make a strided SLAM map backed by ArrayView storage.
 *
 * \param set   pointer to the map's set (must outlive the map)
 * \param stride runtime stride (#values per set element)
 * \param data  backing storage as an ArrayView (must outlive the map)
 */
template <typename SetType, typename T, typename PosType = typename SetType::PositionType>
auto make_map(SetType* set, PosType stride, axom::ArrayView<T> data)
{
  using Indirection = policies::ArrayViewIndirection<PosType, T>;
  using Stride = policies::RuntimeStride<PosType>;
  using MapType = Map<T, SetType, Indirection, Stride>;
  return MapType(set, data, stride);
}

/*!
 * \brief Make a stride-one SLAM map backed by ArrayView storage.
 */
template <typename SetType, typename T, typename PosType = typename SetType::PositionType>
auto make_map(SetType* set, axom::ArrayView<T> data)
{
  using Indirection = policies::ArrayViewIndirection<PosType, T>;
  using Stride = policies::StrideOne<PosType>;
  using MapType = Map<T, SetType, Indirection, Stride>;
  return MapType(set, data);
}

/*!
 * \brief Make a strided SLAM map backed by a raw pointer buffer.
 *
 * This overload wraps the buffer as an ArrayView with length `set->size() * stride`
 * and returns an ArrayView-backed map.
 */
template <typename SetType, typename T, typename PosType = typename SetType::PositionType>
auto make_map(SetType* set, PosType stride, T* data)
{
  const auto n = detail::map_storage_size(set, stride);
  return make_map(set, stride, axom::ArrayView<T>(data, n));
}

/*!
 * \brief Make a stride-one SLAM map backed by a raw pointer buffer.
 */
template <typename SetType, typename T, typename PosType = typename SetType::PositionType>
auto make_map(SetType* set, T* data)
{
  const auto n = detail::map_storage_size(set, PosType {1});
  return make_map(set, axom::ArrayView<T>(data, n));
}

/*!
 * \brief Make a compile-time strided SLAM map backed by ArrayView storage.
 *
 * \tparam STRIDE number of values per set element
 */
template <int STRIDE, typename SetType, typename T, typename PosType = typename SetType::PositionType>
auto make_map_ct(SetType* set, axom::ArrayView<T> data)
{
  using Indirection = policies::ArrayViewIndirection<PosType, T>;
  using Stride = policies::CompileTimeStride<PosType, static_cast<PosType>(STRIDE)>;
  using MapType = Map<T, SetType, Indirection, Stride>;
  return MapType(set, data);
}

/*!
 * \brief Make a compile-time strided SLAM map backed by a raw pointer buffer.
 *
 * This overload wraps the buffer as an ArrayView with length `set->size() * STRIDE`
 * and returns an ArrayView-backed map.
 */
template <int STRIDE, typename SetType, typename T, typename PosType = typename SetType::PositionType>
auto make_map_ct(SetType* set, T* data)
{
  const PosType stride = static_cast<PosType>(STRIDE);
  const auto n = detail::map_storage_size(set, stride);
  return make_map_ct<STRIDE>(set, axom::ArrayView<T>(data, n));
}

}  // namespace axom::slam

#endif  // SLAM_MAP_BUILDERS_HPP_
