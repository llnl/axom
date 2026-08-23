// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file Traits.hpp
 *
 * \brief Compatibility trait predicates backed by SLAM's C++20 concepts.
 *
 * New constrained code should use the concepts in Concepts.hpp directly.
 * The variable templates in this header are for callers that need a Boolean constant.
 */

#pragma once

#include "axom/slam/Concepts.hpp"

namespace axom::slam
{
/// \brief Compatibility wrapper for BivariateSetLike.
template <typename T>
inline constexpr bool is_bivariate_set_like_v = BivariateSetLike<T>;

/// \brief Compatibility wrapper for SetLike.
template <typename T>
inline constexpr bool is_set_like_v = SetLike<T>;

/// \brief Compatibility wrapper for OrderedSetLike.
template <typename T>
inline constexpr bool is_ordered_set_like_v = OrderedSetLike<T>;

/// \brief Compatibility wrapper for RelationLike.
template <typename T>
inline constexpr bool is_relation_like_v = RelationLike<T>;

/// \brief Compatibility wrapper for MapLike.
template <typename T>
inline constexpr bool is_map_like_v = MapLike<T>;

/// \brief Compatibility wrapper for MapOver.
template <typename M, typename S>
inline constexpr bool is_map_over_v = MapOver<M, S>;

/// \brief Compatibility wrapper for ValuePolicy.
template <typename T>
inline constexpr bool is_value_policy_v = ValuePolicy<T>;

/// \brief Compatibility wrapper for SizePolicy.
template <typename T>
inline constexpr bool is_size_policy_v = SizePolicy<T>;

/// \brief Compatibility wrapper for StridePolicy.
template <typename T>
inline constexpr bool is_stride_policy_v = StridePolicy<T>;

/// \brief Compatibility wrapper for OffsetPolicy.
template <typename T>
inline constexpr bool is_offset_policy_v = OffsetPolicy<T>;

/// \brief Compatibility wrapper for IndirectionPolicy.
template <typename T>
inline constexpr bool is_indirection_policy_v = IndirectionPolicy<T>;

/// \brief Compatibility wrapper for PositionLike.
template <typename T>
inline constexpr bool is_position_like_v = PositionLike<T>;

/*!
 * \brief Compatibility wrapper for the former handle/device-capture predicate.
 *
 * Trivial copyability is a deployment property, not the definition of an element handle.
 * New code should use DeviceCapturable when that is the requirement.
 */
template <typename T>
inline constexpr bool is_handle_like_v = DeviceCapturable<T>;

}  // namespace axom::slam
