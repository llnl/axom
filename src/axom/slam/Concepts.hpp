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
 * The concepts in this header describe public operations and relationships between associated types.
 * They do not require concrete SLAM implementations or include the standard ranges library.
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
};

template <typename T>
concept HasMapAssociatedTypes = requires {
  typename model_t<T>::DataType;
  typename model_t<T>::SetPosition;
  typename model_t<T>::SetElement;
  typename model_t<T>::ValueType;
  typename model_t<T>::ConstValueType;
};

template <typename T>
concept HasUnivariateMapAssociatedTypes = requires { typename model_t<T>::SetType; };

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
 * \brief Opt-in customization for a non-integral SLAM position type.
 *
 * A specialized type is still responsible for providing the arithmetic
 * and ordering required by the SLAM APIs in which it is used.
 */
template <typename T>
inline constexpr bool enable_position_like = false;

/// \brief A built-in integral or explicitly opted-in SLAM position type.
template <typename T>
concept PositionLike = std::integral<detail::model_t<T>> || enable_position_like<detail::model_t<T>>;

/*!
 * \brief A univariate set with position-based size and element access.
 *
 * \note SetLike and BivariateSetLike are intentionally disjoint.
 */
template <typename T>
concept SetLike = detail::HasSetAssociatedTypes<T> && !detail::HasBivariateSetAssociatedTypes<T> &&
  PositionLike<typename detail::model_t<T>::PositionType> &&
  requires(const detail::model_t<T>& set, typename detail::model_t<T>::PositionType pos) {
    { set.size() } -> std::same_as<typename detail::model_t<T>::PositionType>;
    { set.empty() } -> std::convertible_to<bool>;
    { set.at(pos) } -> std::convertible_to<typename detail::model_t<T>::ElementType>;
  };

namespace detail
{
// Shared glue for constrained construction helpers.
template <typename Set, typename Value>
concept SetPositionConvertible = SetLike<Set> && PositionLike<Value> &&
  std::convertible_to<model_t<Value>, typename model_t<Set>::PositionType>;

template <typename Set, typename Position>
concept SetPositionSame =
  SetLike<Set> && std::same_as<model_t<Position>, typename model_t<Set>::PositionType>;

template <typename Set, typename Position>
concept OptionalSetPositionSame =
  SetLike<Set> && (std::same_as<model_t<Position>, void> || SetPositionSame<Set, Position>);

template <int Stride, typename Position>
consteval bool positiveStaticStrideRepresentable()
{
  using PositionType = model_t<Position>;
  if constexpr(std::integral<PositionType>)
  {
    return std::in_range<PositionType>(Stride);
  }
  else
  {
    // Opted-in position types state their own construction contract.
    return std::constructible_from<PositionType, int>;
  }
}

template <int Stride, typename Position>
concept PositiveStaticStrideForPosition =
  PositionLike<Position> && (Stride > 0) && positiveStaticStrideRepresentable<Stride, Position>();

template <int Stride, typename Set>
concept PositiveStaticStrideFor =
  SetLike<Set> && PositiveStaticStrideForPosition<Stride, typename model_t<Set>::PositionType>;
}  // namespace detail

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

/*!
 * \brief A set whose elements are indexed by positions from two component sets.
 *
 * PositionType indexes the flattened sequence. ElementType is a coordinate
 * whose `first` and `second` members have the exact endpoint position types.
 * Row-local positions are inferred from row operations rather than imposed as
 * another required associated type.
 */
template <typename T>
concept BivariateSetLike =
  detail::HasSetAssociatedTypes<T> && detail::HasBivariateSetAssociatedTypes<T> &&
  PositionLike<typename detail::model_t<T>::PositionType> &&
  SetLike<typename detail::model_t<T>::FirstSetType> &&
  SetLike<typename detail::model_t<T>::SecondSetType> &&
  requires(typename detail::model_t<T>::ElementType coordinate) {
    requires std::same_as<std::remove_cvref_t<decltype(coordinate.first)>,
                          typename detail::model_t<T>::FirstSetType::PositionType>;
    requires std::same_as<std::remove_cvref_t<decltype(coordinate.second)>,
                          typename detail::model_t<T>::SecondSetType::PositionType>;
  } &&
  requires(const detail::model_t<T>& set,
           typename detail::model_t<T>::PositionType flatPosition,
           typename detail::model_t<T>::FirstSetType::PositionType firstPosition,
           typename detail::model_t<T>::SecondSetType::PositionType secondPosition) {
    { set.size() } -> std::same_as<typename detail::model_t<T>::PositionType>;
    { set.at(flatPosition) } -> std::convertible_to<typename detail::model_t<T>::ElementType>;
    { set.getFirstSet() } -> std::same_as<const typename detail::model_t<T>::FirstSetType*>;
    { set.getSecondSet() } -> std::same_as<const typename detail::model_t<T>::SecondSetType*>;
    { set.getElements(firstPosition).size() } -> PositionLike;
    set.getElements(firstPosition).begin();
    {
      set.getElements(firstPosition).end()
    } -> std::same_as<decltype(set.getElements(firstPosition).begin())>;
    {
      *set.getElements(firstPosition).begin()
    } -> std::convertible_to<typename detail::model_t<T>::SecondSetType::PositionType>;
    {
      set.findElementIndex(firstPosition, secondPosition)
    } -> std::same_as<decltype(set.getElements(firstPosition).size())>;
    {
      set.findElementFlatIndex(firstPosition, secondPosition)
    } -> std::same_as<typename detail::model_t<T>::PositionType>;
    {
      set.flatToFirstIndex(flatPosition)
    } -> std::same_as<typename detail::model_t<T>::FirstSetType::PositionType>;
    {
      set.flatToSecondIndex(flatPosition)
    } -> std::same_as<typename detail::model_t<T>::SecondSetType::PositionType>;
  };

/*!
 * \brief A relation that exposes its two sets and a const iterable row for a
 * position in the from-set.
 *
 * A relation row is selected by the FromSetType's PositionType and contains
 * positions in the ToSetType. Endpoint elements are obtained by projecting
 * those positions through the endpoint sets. A flattened representation is an
 * implementation capability, not a requirement of the relation abstraction.
 */
template <typename T>
concept RelationLike =
  detail::HasRelationAssociatedTypes<T> && SetLike<typename detail::model_t<T>::FromSetType> &&
  SetLike<typename detail::model_t<T>::ToSetType> &&
  requires(const detail::model_t<T>& relation,
           typename detail::model_t<T>::FromSetType::PositionType fromPosition) {
    { relation.fromSet() } -> std::same_as<const typename detail::model_t<T>::FromSetType*>;
    { relation.toSet() } -> std::same_as<const typename detail::model_t<T>::ToSetType*>;
    { relation[fromPosition].size() } -> PositionLike;
    relation[fromPosition].begin();
    { relation[fromPosition].end() } -> std::same_as<decltype(relation[fromPosition].begin())>;
    {
      *relation[fromPosition].begin()
    } -> std::convertible_to<typename detail::model_t<T>::ToSetType::PositionType>;
  };

namespace detail
{
template <typename Value, typename Data>
concept MapValueFor = std::same_as<std::remove_cvref_t<Value>, std::remove_cvref_t<Data>>;

template <typename T>
concept CommonMapModel = HasMapAssociatedTypes<T> &&
  MapValueFor<typename model_t<T>::ValueType, typename model_t<T>::DataType> &&
  MapValueFor<typename model_t<T>::ConstValueType, typename model_t<T>::DataType> &&
  requires(model_t<T>& map, const model_t<T>& constMap, typename model_t<T>::SetPosition pos) {
    { constMap.size() } -> std::same_as<typename model_t<T>::SetPosition>;
    { map[pos] } -> std::same_as<typename model_t<T>::ValueType>;
    { constMap[pos] } -> std::same_as<typename model_t<T>::ConstValueType>;
  };
}  // namespace detail

/// \brief A map whose domain is a univariate SetType.
template <typename T>
concept UnivariateMapLike = detail::CommonMapModel<T> &&
  detail::HasUnivariateMapAssociatedTypes<T> && SetLike<typename detail::model_t<T>::SetType> &&
  std::same_as<typename detail::model_t<T>::SetPosition, typename detail::model_t<T>::SetType::PositionType> &&
  std::same_as<typename detail::model_t<T>::SetElement, typename detail::model_t<T>::SetType::ElementType> &&
  requires(const detail::model_t<T>& map) {
    { map.set() } -> std::same_as<const typename detail::model_t<T>::SetType*>;
  };

/// \brief A map whose domain is a BivariateSetType.
template <typename T>
concept BivariateMapLike = detail::CommonMapModel<T> && detail::HasBivariateMapAssociatedTypes<T> &&
  BivariateSetLike<typename detail::model_t<T>::BivariateSetType> &&
  std::same_as<typename detail::model_t<T>::SetPosition,
               typename detail::model_t<T>::BivariateSetType::PositionType> &&
  std::same_as<typename detail::model_t<T>::SetElement,
               typename detail::model_t<T>::BivariateSetType::ElementType> &&
  requires(const detail::model_t<T>& map) {
    { map.set() } -> std::same_as<const typename detail::model_t<T>::BivariateSetType*>;
  };

/// \brief A univariate or bivariate SLAM map.
template <typename T>
concept MapLike = UnivariateMapLike<T> || BivariateMapLike<T>;

/// \brief A map whose semantic domain is exactly S.
template <typename M, typename S>
concept MapOver =
  (UnivariateMapLike<M> && std::same_as<typename detail::model_t<M>::SetType, detail::model_t<S>>) ||
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

/// \brief A SizePolicy usable by a set whose position type is Position.
template <typename T, typename Position>
concept SetSizePolicyFor = SizePolicy<T> && PositionLike<Position> &&
  std::same_as<detail::policy_default_t<T>, detail::model_t<Position>> &&
  std::constructible_from<detail::model_t<T>, detail::model_t<Position>>;

/*!
 * \brief The common capability shared by scalar and multi-dimensional stride policies.
 *
 * Use OrderedSetStridePolicyFor or MapStridePolicyFor when checking whether a stride
 * can actually be substituted into one of those owners.
 */
template <typename T>
concept StridePolicy = requires(const detail::model_t<T>& policy) {
  typename detail::model_t<T>::IndexType;
  typename detail::model_t<T>::ShapeType;
  detail::model_t<T>::NumDims;
  { detail::model_t<T>::DefaultSize() } -> std::same_as<typename detail::model_t<T>::ShapeType>;
  { policy.stride() } -> std::same_as<typename detail::model_t<T>::IndexType>;
  { policy.shape() } -> std::same_as<typename detail::model_t<T>::ShapeType>;
};

/*!
 * \brief A scalar stride policy usable by OrderedSet with Position.
 *
 * OrderedSet constructs its stride from a position and validates it as a scalar value policy.
 * Multi-dimensional map strides do not satisfy this refinement.
 */
template <typename T, typename Position>
concept OrderedSetStridePolicyFor = StridePolicy<T> && ValuePolicy<T> && PositionLike<Position> &&
  std::same_as<typename detail::model_t<T>::IndexType, detail::model_t<Position>> &&
  std::same_as<typename detail::model_t<T>::IntType, detail::model_t<Position>> &&
  std::same_as<typename detail::model_t<T>::ShapeType, detail::model_t<Position>> &&
  (detail::model_t<T>::NumDims == 1) &&
  std::constructible_from<detail::model_t<T>, detail::model_t<Position>> && requires {
    detail::model_t<T>::DEFAULT_VALUE;
    requires std::same_as<detail::policy_default_t<T>, detail::model_t<Position>>;
  };

/// \brief A scalar or multi-dimensional stride policy usable by Map with Position.
template <typename T, typename Position>
concept MapStridePolicyFor = StridePolicy<T> && PositionLike<Position> &&
  PositionLike<typename detail::model_t<T>::IndexType> &&
  std::convertible_to<typename detail::model_t<T>::IndexType, detail::model_t<Position>> &&
  (detail::model_t<T>::NumDims > 0) &&
  std::constructible_from<detail::model_t<T>, typename detail::model_t<T>::ShapeType> &&
  ((detail::model_t<T>::NumDims == 1) || requires(const detail::model_t<T>& policy) {
                               {
                                 policy.strides()
                               } -> std::same_as<typename detail::model_t<T>::ShapeType>;
                             });

/// \brief A scalar value policy that reports an offset.
template <typename T>
concept OffsetPolicy = ValuePolicy<T> && requires(const detail::model_t<T>& policy) {
  detail::model_t<T>::DEFAULT_VALUE;
  { policy.offset() } -> std::same_as<typename detail::model_t<T>::IntType>;
};

/// \brief An OffsetPolicy usable by an OrderedSet whose position type is Position.
template <typename T, typename Position>
concept OrderedSetOffsetPolicyFor = OffsetPolicy<T> && PositionLike<Position> &&
  std::same_as<typename detail::model_t<T>::IntType, detail::model_t<Position>> &&
  std::same_as<detail::policy_default_t<T>, detail::model_t<Position>> &&
  std::constructible_from<detail::model_t<T>, detail::model_t<Position>>;

/*!
 * \brief The common storage/indirection-policy capability.
 *
 * Use IndirectionPolicyFor when the calling position type is available and the
 * indirection operation itself should also be checked.
 * Use OrderedSetIndirectionPolicyFor or MapIndirectionPolicyFor
 * when checking substitutability into those owners.
 */
template <typename T>
concept IndirectionPolicy =
  detail::HasIndirectionAssociatedTypes<T> && requires(const detail::model_t<T>& policy) {
    detail::model_t<T>::DeviceAccessible;
    { policy.hasIndirection() } -> std::convertible_to<bool>;
  };

/// \brief An IndirectionPolicy callable with Position.
template <typename T, typename Position>
concept IndirectionPolicyFor = IndirectionPolicy<T> &&
  requires(detail::model_t<T>& policy,
           const detail::model_t<T>& constPolicy,
           detail::model_t<Position> pos) {
    {
      policy.indirection(pos)
    } -> std::convertible_to<typename detail::model_t<T>::IndirectionResult>;
    {
      constPolicy.indirection(pos)
    } -> std::convertible_to<typename detail::model_t<T>::ConstIndirectionResult>;
  };

namespace detail
{
template <typename T>
concept HasTypedIndirectionAssociatedTypes = requires {
  typename model_t<T>::PositionType;
  typename model_t<T>::ElementType;
};

template <typename T>
concept HasMapIndirectionAssociatedTypes = requires {
  typename model_t<T>::IndirectionRefType;
  typename model_t<T>::IndirectionConstRefType;
  typename model_t<T>::ResultPtr;
  typename model_t<T>::ConstResultPtr;
  model_t<T>::IsMutableBuffer;
  std::integral_constant<bool, model_t<T>::IsMutableBuffer> {};
};

template <typename T, typename Position>
concept MapBufferFor = HasMapIndirectionAssociatedTypes<T> &&
  requires(const typename model_t<T>::IndirectionBufferType& buffer) {
    { buffer.size() } -> std::convertible_to<model_t<Position>>;
    { buffer.empty() } -> std::convertible_to<bool>;
  } &&
  (!model_t<T>::IsMutableBuffer ||
   requires(typename model_t<T>::IndirectionBufferType& buffer, model_t<Position> size) {
     buffer.resize(size);
   });
}  // namespace detail

/// \brief An indirection policy usable by OrderedSet over Position and Element.
template <typename T, typename Position, typename Element>
concept OrderedSetIndirectionPolicyFor = IndirectionPolicyFor<T, Position> &&
  PositionLike<Position> && detail::HasTypedIndirectionAssociatedTypes<T> &&
  std::same_as<typename detail::model_t<T>::PositionType, detail::model_t<Position>> &&
  std::same_as<typename detail::model_t<T>::ElementType, std::remove_reference_t<Element>> &&
  std::same_as<std::remove_cvref_t<typename detail::model_t<T>::IndirectionResult>,
               detail::model_t<Element>> &&
  std::same_as<std::remove_cvref_t<typename detail::model_t<T>::ConstIndirectionResult>,
               detail::model_t<Element>> &&
  std::default_initializable<detail::model_t<T>> &&
  std::constructible_from<detail::model_t<T>, typename detail::model_t<T>::IndirectionPtrType> &&
  requires(const detail::model_t<T>& policy,
           detail::model_t<Position> size,
           detail::model_t<Position> offset,
           detail::model_t<Position> stride) {
    { policy.isValid(size, offset, stride, false) } -> std::convertible_to<bool>;
  };

/// \brief An indirection policy providing Map's buffer and static access API.
template <typename T, typename Position, typename Data>
concept MapIndirectionPolicyFor =
  IndirectionPolicy<T> && PositionLike<Position> && detail::HasTypedIndirectionAssociatedTypes<T> &&
  detail::HasMapIndirectionAssociatedTypes<T> && detail::MapBufferFor<T, Position> &&
  std::same_as<typename detail::model_t<T>::PositionType, detail::model_t<Position>> &&
  std::same_as<typename detail::model_t<T>::ElementType, std::remove_reference_t<Data>> &&
  detail::MapValueFor<typename detail::model_t<T>::IndirectionResult, Data> &&
  detail::MapValueFor<typename detail::model_t<T>::ConstIndirectionResult, Data> &&
  requires(typename detail::model_t<T>::IndirectionBufferType& buffer,
           const typename detail::model_t<T>::IndirectionBufferType& constBuffer,
           detail::model_t<Position> pos) {
    {
      detail::model_t<T>::getIndirection(buffer, pos)
    } -> std::same_as<typename detail::model_t<T>::ResultPtr>;
    {
      detail::model_t<T>::getConstIndirection(constBuffer, pos)
    } -> std::same_as<typename detail::model_t<T>::ConstResultPtr>;
  };

/// \brief A MapIndirectionPolicyFor that can allocate and initialize its buffer.
template <typename T, typename Position, typename Data>
concept AllocatingMapIndirectionPolicyFor = MapIndirectionPolicyFor<T, Position, Data> &&
  requires(detail::model_t<Position> size, const detail::model_t<Data>& value, int allocatorId) {
    {
      detail::model_t<T>::create(size, value, allocatorId)
    } -> std::same_as<typename detail::model_t<T>::IndirectionBufferType>;
  };

/// \brief A non-reference type that can be copied byte-for-byte into device code.
template <typename T>
concept DeviceCapturable =
  !std::is_reference_v<T> && std::is_trivially_copyable_v<std::remove_cv_t<T>>;

}  // namespace axom::slam
