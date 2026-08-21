// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file Concepts.hpp
 *
 * \brief C++20 concepts for SLAM containers, policies, and index properties.
 *
 * The concepts in this header describe public operations and relationships
 * between associated types. They intentionally do not require concrete SLAM
 * implementations or include the standard ranges library.
 */

#pragma once

#include <concepts>
#include <type_traits>
#include <utility>

namespace axom::slam
{
namespace detail
{
template <typename T>
using model_t = std::remove_cvref_t<T>;

template <typename T>
concept HasSetAssociatedTypes = requires {
  typename model_t<T>::PositionType;
  typename model_t<T>::ElementType;
};

template <typename T>
concept HasBivariateSetAssociatedTypes = requires {
  typename model_t<T>::FirstSetType;
  typename model_t<T>::SecondSetType;
};

template <typename T>
concept HasRelationAssociatedTypes = requires {
  typename model_t<T>::FromSetType;
  typename model_t<T>::ToSetType;
  typename model_t<T>::SetPosition;
  typename model_t<T>::SetElement;
};

template <typename T>
concept HasMapAssociatedTypes = requires {
  typename model_t<T>::DataType;
  typename model_t<T>::SetType;
  typename model_t<T>::SetPosition;
  typename model_t<T>::SetElement;
  typename model_t<T>::ValueType;
  typename model_t<T>::ConstValueType;
};

template <typename T>
concept HasBivariateMapAssociatedTypes = requires { typename model_t<T>::BivariateSetType; };

template <typename T>
concept HasIndirectionAssociatedTypes = requires {
  typename model_t<T>::IndirectionResult;
  typename model_t<T>::ConstIndirectionResult;
  typename model_t<T>::IndirectionBufferType;
  typename model_t<T>::IndirectionPtrType;
};

template <typename T>
using policy_default_t = std::remove_cv_t<decltype(model_t<T>::DEFAULT_VALUE)>;
}  // namespace detail

/*!
 * \brief A univariate set with position-based size and element access.
 *
 * Bivariate sets are excluded because their positions describe pairs from two
 * different sets and have a distinct access contract.
 */
template <typename T>
concept SetLike = detail::HasSetAssociatedTypes<T> &&
  !detail::HasBivariateSetAssociatedTypes<T> &&
  requires(const detail::model_t<T>& set, typename detail::model_t<T>::PositionType pos) {
    { set.size() } -> std::same_as<typename detail::model_t<T>::PositionType>;
    { set.empty() } -> std::convertible_to<bool>;
    { set.at(pos) } -> std::convertible_to<typename detail::model_t<T>::ElementType>;
  };

/*!
 * \brief A SetLike type with const iteration over its elements.
 *
 * This concept states SLAM's ordered-set surface without importing the standard ranges taxonomy.
 * Standard iterator and range categories are tested separately.
 */
template <typename T>
concept OrderedSetLike = SetLike<T> && requires(const detail::model_t<T>& set) {
  set.begin();
  { set.end() } -> std::same_as<decltype(set.begin())>;
  { *set.begin() } -> std::convertible_to<typename detail::model_t<T>::ElementType>;
};

/// \brief A set whose elements are indexed by positions from two component sets.
template <typename T>
concept BivariateSetLike = detail::HasSetAssociatedTypes<T> &&
  detail::HasBivariateSetAssociatedTypes<T> &&
  SetLike<typename detail::model_t<T>::FirstSetType> &&
  SetLike<typename detail::model_t<T>::SecondSetType> &&
  requires(const detail::model_t<T>& set, typename detail::model_t<T>::PositionType pos) {
    { set.size() } -> std::same_as<typename detail::model_t<T>::PositionType>;
    { set.at(pos) } -> std::convertible_to<typename detail::model_t<T>::ElementType>;
    { set.getFirstSet() } ->
      std::same_as<const typename detail::model_t<T>::FirstSetType*>;
    { set.getSecondSet() } ->
      std::same_as<const typename detail::model_t<T>::SecondSetType*>;
    { set.findElementFlatIndex(pos, pos) } ->
      std::same_as<typename detail::model_t<T>::PositionType>;
    { set.flatToFirstIndex(pos) } ->
      std::same_as<typename detail::model_t<T>::PositionType>;
    { set.flatToSecondIndex(pos) } ->
      std::same_as<typename detail::model_t<T>::PositionType>;
  };

/// \brief A relation that exposes its two sets and a const iterable row for a position in the from-set.
template <typename T>
concept RelationLike = detail::HasRelationAssociatedTypes<T> &&
  SetLike<typename detail::model_t<T>::FromSetType> &&
  SetLike<typename detail::model_t<T>::ToSetType> &&
  std::same_as<typename detail::model_t<T>::SetPosition,
               typename detail::model_t<T>::FromSetType::PositionType> &&
  requires(const detail::model_t<T>& relation,
           typename detail::model_t<T>::SetPosition fromPosition) {
    { relation.fromSet() } ->
      std::same_as<const typename detail::model_t<T>::FromSetType*>;
    { relation.toSet() } ->
      std::same_as<const typename detail::model_t<T>::ToSetType*>;
    relation[fromPosition];
    relation.begin(fromPosition);
    { relation.end(fromPosition) } -> std::same_as<decltype(relation.begin(fromPosition))>;
    { *relation.begin(fromPosition) } ->
      std::convertible_to<typename detail::model_t<T>::SetElement>;
  };

namespace detail
{
template <typename T>
concept CommonMapModel = HasMapAssociatedTypes<T> &&
  requires(const model_t<T>& map, typename model_t<T>::SetPosition pos) {
    { map.size() } -> std::same_as<typename model_t<T>::SetPosition>;
    { map[pos] } -> std::convertible_to<typename model_t<T>::ConstValueType>;
  };
}  // namespace detail

/// \brief A map whose domain is a univariate SetType.
template <typename T>
concept UnivariateMapLike = detail::CommonMapModel<T> &&
  !detail::HasBivariateMapAssociatedTypes<T> &&
  SetLike<typename detail::model_t<T>::SetType> &&
  std::same_as<typename detail::model_t<T>::SetPosition,
               typename detail::model_t<T>::SetType::PositionType> &&
  std::same_as<typename detail::model_t<T>::SetElement,
               typename detail::model_t<T>::SetType::ElementType> &&
  requires(const detail::model_t<T>& map) {
    { map.set() } -> std::same_as<const typename detail::model_t<T>::SetType*>;
  };

/// \brief A map whose domain is a BivariateSetType.
template <typename T>
concept BivariateMapLike = detail::CommonMapModel<T> &&
  detail::HasBivariateMapAssociatedTypes<T> &&
  BivariateSetLike<typename detail::model_t<T>::BivariateSetType> &&
  std::same_as<typename detail::model_t<T>::SetPosition,
               typename detail::model_t<T>::BivariateSetType::PositionType> &&
  std::same_as<typename detail::model_t<T>::SetElement,
               typename detail::model_t<T>::BivariateSetType::ElementType> &&
  requires(const detail::model_t<T>& map) {
    { map.set() } ->
      std::same_as<const typename detail::model_t<T>::BivariateSetType*>;
  };

/// \brief A univariate or bivariate SLAM map.
template <typename T>
concept MapLike = UnivariateMapLike<T> || BivariateMapLike<T>;

/// \brief A map whose semantic domain is exactly S.
template <typename M, typename S>
concept MapOver =
  (UnivariateMapLike<M> &&
   std::same_as<typename detail::model_t<M>::SetType, detail::model_t<S>>) ||
  (BivariateMapLike<M> &&
   std::same_as<typename detail::model_t<M>::BivariateSetType, detail::model_t<S>>);

/// \brief A scalar runtime or compile-time value policy.
template <typename T>
concept ValuePolicy = requires(const detail::model_t<T>& policy) {
  typename detail::model_t<T>::TagType;
  typename detail::model_t<T>::IntType;
  { policy.value() } -> std::same_as<typename detail::model_t<T>::IntType>;
  { policy.isValid(false) } -> std::convertible_to<bool>;
};

/// \brief A policy that reports a size and whether that size is empty.
template <typename T>
concept SizePolicy = requires(const detail::model_t<T>& policy) {
  detail::model_t<T>::DEFAULT_VALUE;
  { policy.size() } -> std::same_as<detail::policy_default_t<T>>;
  { policy.empty() } -> std::convertible_to<bool>;
  { policy.isValid(false) } -> std::convertible_to<bool>;
};

/// \brief A one- or multi-dimensional map/set stride policy.
template <typename T>
concept StridePolicy = requires(const detail::model_t<T>& policy) {
  typename detail::model_t<T>::IndexType;
  typename detail::model_t<T>::ShapeType;
  detail::model_t<T>::NumDims;
  { detail::model_t<T>::DefaultSize() } ->
    std::same_as<typename detail::model_t<T>::ShapeType>;
  { policy.stride() } -> std::same_as<typename detail::model_t<T>::IndexType>;
  { policy.shape() } -> std::same_as<typename detail::model_t<T>::ShapeType>;
};

/// \brief A scalar value policy that reports an offset.
template <typename T>
concept OffsetPolicy = ValuePolicy<T> && requires(const detail::model_t<T>& policy) {
  detail::model_t<T>::DEFAULT_VALUE;
  { policy.offset() } -> std::same_as<typename detail::model_t<T>::IntType>;
};

/*!
 * \brief A storage/indirection policy independent of a particular position type.
 *
 * Use IndirectionPolicyFor when the calling position type is available and the
 * indirection operation itself should also be checked.
 */
template <typename T>
concept IndirectionPolicy = detail::HasIndirectionAssociatedTypes<T> &&
  requires(const detail::model_t<T>& policy) {
    detail::model_t<T>::DeviceAccessible;
    { policy.hasIndirection() } -> std::convertible_to<bool>;
  };

/// \brief An IndirectionPolicy callable with Position.
template <typename T, typename Position>
concept IndirectionPolicyFor = IndirectionPolicy<T> &&
  requires(detail::model_t<T>& policy,
           const detail::model_t<T>& constPolicy,
           detail::model_t<Position> pos) {
    { policy.indirection(pos) } ->
      std::convertible_to<typename detail::model_t<T>::IndirectionResult>;
    { constPolicy.indirection(pos) } ->
      std::convertible_to<typename detail::model_t<T>::ConstIndirectionResult>;
  };

/*!
 * \brief Opt-in customization for a non-integral SLAM position type.
 *
 * A specialized type is still responsible for providing the arithmetic and
 * ordering required by the SLAM APIs in which it is used.
 */
template <typename T>
inline constexpr bool enable_position_like = false;

/// \brief A built-in integral or explicitly opted-in SLAM position type.
template <typename T>
concept PositionLike = std::integral<detail::model_t<T>> ||
  enable_position_like<detail::model_t<T>>;

/// \brief A non-reference type that can be copied byte-for-byte into device code.
template <typename T>
concept DeviceCapturable = !std::is_reference_v<T> &&
  std::is_trivially_copyable_v<std::remove_cv_t<T>>;

}  // namespace axom::slam
