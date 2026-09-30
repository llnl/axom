// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file Concepts.hpp
 * \brief Semantic concepts for Slam sets, relations, and maps.
 *
 * Object concepts describe positional access, related elements and mapped values.
 * They ignore top-level const and reference qualification, but preserve the
 * constness of elements and values.
 */

#pragma once

#include <concepts>
#include <limits>
#include <type_traits>
#include <utility>

namespace axom::slam
{

// ----------------------------------------------------------------------------
// Helper type aliases for public concepts
// ----------------------------------------------------------------------------
namespace detail
{
template <typename T>
using model_t = std::remove_cvref_t<T>;

template <typename T>
using position_t = typename model_t<T>::PositionType;

template <typename T>
using element_t = typename model_t<T>::ElementType;

/// Integer arguments may be unsigned, but cannot be Boolean.
template <typename T>
concept PositionValueLike = std::integral<model_t<T>> && !std::same_as<model_t<T>, bool>;

/// The operations needed to traverse values, without prescribing iterator aliases.
template <typename Iterator, typename Sentinel, typename Value>
concept IteratesAs = requires(Iterator it, Sentinel end) {
  { *it } -> std::convertible_to<Value>;
  ++it;
  { it != end } -> std::convertible_to<bool>;
};

template <typename R, typename Value>
concept IterableValues = requires(const model_t<R>& range) {
  requires IteratesAs<decltype(range.begin()), decltype(range.end()), Value>;
};

template <typename R, typename Value>
concept SizedValues = IterableValues<R, Value> && requires(const model_t<R>& range) {
  { range.size() } -> PositionValueLike;
};

/// Scalar references may add constness, but cannot discard element qualification.
template <typename Reference, typename Data>
concept MapReferenceFor = std::is_lvalue_reference_v<Reference> &&
  (std::same_as<std::remove_reference_t<Reference>, Data> ||
   std::same_as<std::remove_reference_t<Reference>, const Data>);
}  // namespace detail

// ----------------------------------------------------------------------------
// Public concepts
// ----------------------------------------------------------------------------

/// \brief A signed integral type for Slam positions and indexing arithmetic.
template <typename T>
concept PositionLike = std::signed_integral<detail::model_t<T>>;

/// \brief A sized collection with positional access to its elements.
/// size() is nonnegative, and empty() agrees with size() == 0. Access requires a valid position.
template <typename T>
concept SetLike = PositionLike<detail::position_t<T>> &&
  requires(const detail::model_t<T>& set, detail::position_t<T> pos) {
    typename detail::element_t<T>;
    { set.size() } -> std::same_as<detail::position_t<T>>;
    { set.empty() } -> std::convertible_to<bool>;
    { set.at(pos) } -> std::convertible_to<detail::element_t<T>>;
  };

/// \brief A SetLike type whose const traversal visits the elements in positional order.
template <typename T>
concept IterableSetLike = SetLike<T> && detail::IterableValues<T, detail::element_t<T>>;

/// \brief A container that can check its own internal consistency.
/// Validation is a separate capability, not a requirement of every set or relation.
template <typename T>
concept Validatable = requires(const detail::model_t<T>& container) {
  { container.isValid(false) } -> std::convertible_to<bool>;
};

/// \brief A set of coordinate pairs with access to the subset for each first-set position.
/// getElements(first) provides a sized, iterable collection of second-set positions.
/// Visiting these subsets in first-set order agrees with at(flat), and their sizes
/// sum to size(). Search and flat-index conversion are not required.
template <typename T>
concept BivariateSetLike = SetLike<T> && SetLike<typename detail::model_t<T>::FirstSetType> &&
  SetLike<typename detail::model_t<T>::SecondSetType> &&
  requires(detail::element_t<T> coordinate) {
    requires std::same_as<std::remove_cvref_t<decltype(coordinate.first)>,
                          typename detail::model_t<T>::FirstSetType::PositionType>;
    requires std::same_as<std::remove_cvref_t<decltype(coordinate.second)>,
                          typename detail::model_t<T>::SecondSetType::PositionType>;
  } &&
  requires(const detail::model_t<T>& set,
           typename detail::model_t<T>::FirstSetType::PositionType first) {
    { set.getFirstSet() } -> std::same_as<const typename detail::model_t<T>::FirstSetType*>;
    { set.getSecondSet() } -> std::same_as<const typename detail::model_t<T>::SecondSetType*>;
    requires detail::SizedValues<decltype(set.getElements(first)),
                                 typename detail::model_t<T>::SecondSetType::PositionType>;
  };

/// \brief A relation that provides to-set positions for each from-set position.
/// relation[from] returns a sized, iterable collection of valid to-set positions.
/// The sets may contain coordinates or other element types. Flat storage,
/// validation methods and a named subset type are not required.
template <typename T>
concept RelationLike = SetLike<typename detail::model_t<T>::FromSetType> &&
  SetLike<typename detail::model_t<T>::ToSetType> &&
  requires(const detail::model_t<T>& relation,
           typename detail::model_t<T>::FromSetType::PositionType from) {
    { relation.fromSet() } -> std::same_as<const typename detail::model_t<T>::FromSetType*>;
    { relation.toSet() } -> std::same_as<const typename detail::model_t<T>::ToSetType*>;
    requires detail::SizedValues<decltype(relation[from]),
                                 typename detail::model_t<T>::ToSetType::PositionType>;
  };

/// \brief Values associated with entries, addressed by entry position and local component.
/// index(pos) identifies the associated set element. flatValue(pos, component)
/// accesses one scalar component, including for tensor-valued maps. "flat" here
/// selects an entry of a bivariate set, not a global component-storage position.
/// Bound maps have positive numComp(). A default unbound SubMap may be empty with zero components.
template <typename T>
concept MapLike = PositionLike<detail::position_t<T>> &&
  requires(detail::model_t<T>& map,
           const detail::model_t<T>& constMap,
           detail::position_t<T> pos,
           detail::position_t<T> component) {
    typename detail::model_t<T>::DataType;
    { constMap.size() } -> std::same_as<detail::position_t<T>>;
    { constMap.numComp() } -> std::same_as<detail::position_t<T>>;
    constMap.index(pos);
    requires std::is_object_v<std::remove_cvref_t<decltype(constMap.index(pos))>>;
    {
      map.flatValue(pos, component)
    } -> detail::MapReferenceFor<typename detail::model_t<T>::DataType>;
    {
      constMap.flatValue(pos, component)
    } -> detail::MapReferenceFor<typename detail::model_t<T>::DataType>;
  };

/// \brief A MapLike type explicitly bound to all positions of exactly S.
/// MappedSetType names that binding and set() returns it. size() agrees with
/// set()->size(), and index(pos) agrees with set()->at(pos). SubMap has no such
/// binding. Its set() selects parent positions, and index() returns their associated elements.
template <typename M, typename S>
concept MapOver = MapLike<M> && SetLike<S> &&
  std::same_as<typename detail::model_t<M>::MappedSetType, detail::model_t<S>> &&
  std::same_as<detail::position_t<M>, detail::position_t<S>> &&
  requires(const detail::model_t<M>& map, detail::position_t<M> pos) {
    { map.set() } -> std::same_as<const detail::model_t<S>*>;
    { map.index(pos) } -> std::convertible_to<detail::element_t<S>>;
  };

/// \brief A non-reference type with a trivially copyable C++ representation.
/// This does not certify device-callable operations or the accessibility/lifetime
/// of referenced objects. Check those requirements separately before device use.
template <typename T>
concept TriviallyCopyableRepresentation =
  !std::is_reference_v<T> && std::is_trivially_copyable_v<std::remove_cv_t<T>>;

// ----------------------------------------------------------------------------
// Helper concepts built from the public concepts
// ----------------------------------------------------------------------------
namespace detail
{
/// Type-level representability of nonnegative positions and sizes.
template <typename Position, typename RepresentedPosition>
concept PositionCanRepresent = PositionLike<Position> && PositionLike<RepresentedPosition> &&
  (std::numeric_limits<model_t<Position>>::digits >=
   std::numeric_limits<model_t<RepresentedPosition>>::digits);

template <int Stride, typename Position>
concept PositiveStaticStrideForPosition = PositionLike<Position> && (Stride > 0) &&
  (Stride <= std::numeric_limits<model_t<Position>>::max());

template <int Stride, typename Set>
concept PositiveStaticStrideFor =
  SetLike<Set> && PositiveStaticStrideForPosition<Stride, position_t<Set>>;

template <typename Set, typename Value>
concept SetPositionConvertible =
  SetLike<Set> && PositionValueLike<Value> && std::convertible_to<model_t<Value>, position_t<Set>>;

template <typename Set, typename Position>
concept SetPositionSame = SetLike<Set> && std::same_as<model_t<Position>, position_t<Set>>;

template <typename Set, typename Position>
concept OptionalSetPositionSame =
  SetLike<Set> && (std::same_as<model_t<Position>, void> || SetPositionSame<Set, Position>);

/// Flat storage consumed by RelationSet, not a public relation requirement.
template <typename T>
using relation_row_t =
  std::remove_cvref_t<decltype(std::declval<const model_t<T>&>()
                                 [std::declval<typename model_t<T>::FromSetType::PositionType>()])>;

template <typename T>
concept RelationSetSource = RelationLike<model_t<T>> && Validatable<model_t<T>> &&
  requires { typename model_t<T>::FlatPositionType; } &&
  PositionLike<typename model_t<T>::FlatPositionType> &&
  requires(const model_t<T>& relation,
           typename model_t<T>::FromSetType::PositionType fromPosition,
           typename model_t<T>::FlatPositionType flatPosition) {
    { relation.offset(fromPosition) } -> std::convertible_to<typename model_t<T>::FlatPositionType>;
    {
      static_cast<typename model_t<T>::FromSetType::PositionType>(relation.firstIndex(flatPosition))
    } -> std::same_as<typename model_t<T>::FromSetType::PositionType>;
    {
      relation.relationData().size()
    } -> std::convertible_to<typename model_t<T>::FlatPositionType>;
    requires PositionValueLike<decltype(relation.relationData().size())>;
    {
      relation.relationData()[flatPosition]
    } -> std::convertible_to<typename model_t<T>::ToSetType::PositionType>;
  };

/// Exact policy arguments used as base classes must be ordinary inheritable types.
template <typename T>
concept InheritablePolicy = std::is_class_v<T> && std::same_as<T, std::remove_cvref_t<T>> &&
  !std::is_final_v<T> && !std::is_abstract_v<T>;

template <typename T>
using policy_default_t = std::remove_cv_t<decltype(T::DEFAULT_VALUE)>;

/// An indirection result may add constness, but cannot discard element qualification.
template <typename Result, typename Element>
concept IndirectionResultFor = std::same_as<std::remove_reference_t<Result>, Element> ||
  std::same_as<std::remove_reference_t<Result>, const Element>;

}  // namespace detail

/// Reports size, emptiness, and validity, with a default size value.
template <typename T>
concept SizePolicy = requires(const detail::model_t<T>& policy) {
  detail::model_t<T>::DEFAULT_VALUE;
  { policy.size() } -> std::same_as<detail::policy_default_t<detail::model_t<T>>>;
  { policy.empty() } -> std::convertible_to<bool>;
  { policy.isValid(false) } -> std::convertible_to<bool>;
};

/// Reports the signed stride. Shape is an additional requirement of map owners.
template <typename T>
concept StridePolicy = requires(const detail::model_t<T>& policy) {
  { policy.stride() } -> PositionLike;
};

/// Reports a scalar offset, its default value, and validity.
template <typename T>
concept OffsetPolicy = requires(const detail::model_t<T>& policy) {
  detail::model_t<T>::DEFAULT_VALUE;
  { policy.offset() } -> std::same_as<detail::policy_default_t<detail::model_t<T>>>;
  { policy.isValid(false) } -> std::convertible_to<bool>;
};

/// Reports whether a set has a parent and provides its pointer.
/// OrderedSet additionally checks construction and validation with its actual iterators.
template <typename T>
concept SubsetPolicy = requires(const detail::model_t<T>& policy) {
  typename detail::model_t<T>::ParentSetType;
  { policy.isSubset() } -> std::convertible_to<bool>;
  { policy.parentSet() } -> std::convertible_to<const typename detail::model_t<T>::ParentSetType*>;
};

/// Indirection that an OrderedSet can construct, copy, bind and use to access elements.
/// No buffer-container aliases, map accessors, or device flags are required.
template <typename T, typename Position, typename Element>
concept OrderedSetIndirectionPolicyFor = detail::InheritablePolicy<T> && PositionLike<Position> &&
  std::default_initializable<T> && std::copyable<T> &&
  requires {
    typename T::IndirectionPtrType;
    typename T::IndirectionResult;
    typename T::ConstIndirectionResult;
  } && std::constructible_from<T, typename T::IndirectionPtrType> &&
  detail::IndirectionResultFor<typename T::IndirectionResult, Element> &&
  detail::IndirectionResultFor<typename T::ConstIndirectionResult, Element> &&
  requires(T& policy, const T& constPolicy, Position pos) {
    { policy.indirection(pos) } -> std::convertible_to<typename T::IndirectionResult>;
    { constPolicy.indirection(pos) } -> std::convertible_to<typename T::ConstIndirectionResult>;
    requires(!std::is_reference_v<typename T::IndirectionResult> ||
             std::same_as<decltype(policy.indirection(pos)), typename T::IndirectionResult>);
    requires(!std::is_reference_v<typename T::ConstIndirectionResult> ||
             std::same_as<decltype(constPolicy.indirection(pos)), typename T::ConstIndirectionResult>);
    { constPolicy.isValid(pos, pos, pos, false) } -> std::convertible_to<bool>;
  };

/// Static buffer access for Map, which does not inherit the descriptor.
/// Both access paths return stable scalar references. Buffer ownership and
/// referenced allocation accessibility are separate from this type check.
template <typename T, typename Position, typename Data>
concept MapIndirectionPolicyFor =
  std::is_class_v<T> && std::same_as<T, std::remove_cvref_t<T>> && PositionLike<Position> &&
  requires {
    typename T::IndirectionBufferType;
    typename T::IndirectionResult;
    typename T::ConstIndirectionResult;
    typename T::ResultPtr;
    typename T::ConstResultPtr;
    std::integral_constant<bool, T::IsMutableBuffer> {};
  } && detail::MapReferenceFor<typename T::IndirectionResult, Data> &&
  detail::MapReferenceFor<typename T::ConstIndirectionResult, Data> &&
  std::same_as<typename T::ResultPtr,
               std::add_pointer_t<std::remove_reference_t<typename T::IndirectionResult>>> &&
  std::same_as<typename T::ConstResultPtr,
               std::add_pointer_t<std::remove_reference_t<typename T::ConstIndirectionResult>>> &&
  requires(typename T::IndirectionBufferType& buffer,
           const typename T::IndirectionBufferType& constBuffer,
           Position pos) {
    { constBuffer.size() } -> detail::PositionValueLike;
    { constBuffer.empty() } -> std::convertible_to<bool>;
    { T::getIndirection(buffer, pos) } -> std::same_as<typename T::ResultPtr>;
    { T::getConstIndirection(constBuffer, pos) } -> std::same_as<typename T::ConstResultPtr>;
    { T::getIndirection(buffer) } -> std::same_as<typename T::ResultPtr>;
    { T::getConstIndirection(constBuffer) } -> std::same_as<typename T::ConstResultPtr>;
  } &&
  (!T::IsMutableBuffer ||
   requires(typename T::IndirectionBufferType& buffer, Position size) { buffer.resize(size); });

/// Map storage that can allocate and initialize its buffer.
template <typename T, typename Position, typename Data>
concept AllocatingMapIndirectionPolicyFor = MapIndirectionPolicyFor<T, Position, Data> &&
  requires(Position size, const std::remove_cvref_t<Data>& value, int allocatorId) {
    { T::create(size, value, allocatorId) } -> std::same_as<typename T::IndirectionBufferType>;
  };

namespace detail
{
/// Scalar policies are copied into builders, assigned there, and inherited by sets.
template <typename T, typename Position>
concept PolicyDefaultedOver = InheritablePolicy<T> && PositionLike<Position> &&
  std::default_initializable<T> && std::copyable<T> && requires { T::DEFAULT_VALUE; } &&
  std::same_as<policy_default_t<T>, Position> && std::constructible_from<T, Position>;

template <typename T, typename Position>
concept SetSizePolicyFor = SizePolicy<T> && PolicyDefaultedOver<T, Position>;

template <typename T, typename Position>
concept DynamicSetSizePolicyFor = SetSizePolicyFor<T, Position> && requires(T& policy) {
  { policy.size() } -> std::same_as<Position&>;
};

template <typename T, typename Position>
concept OrderedSetOffsetPolicyFor = OffsetPolicy<T> && PolicyDefaultedOver<T, Position>;

template <typename T, typename Position>
concept OrderedSetStridePolicyFor =
  StridePolicy<T> && PolicyDefaultedOver<T, Position> && requires(const T& policy) {
    { policy.stride() } -> std::same_as<Position>;
    { policy.isValid(false) } -> std::convertible_to<bool>;
  };

/// Map constructs this base from a shape. It need not be default-constructible.
template <typename T, typename Position>
concept MapStridePolicyFor = InheritablePolicy<T> && PositionLike<Position> && StridePolicy<T> &&
  requires(const T& policy) {
    typename T::IndexType;
    typename T::ShapeType;
    std::integral_constant<int, T::NumDims> {};
    requires(T::NumDims > 0);
    requires PositionLike<typename T::IndexType>;
    requires std::constructible_from<T, typename T::ShapeType>;
    { T::DefaultSize() } -> std::same_as<typename T::ShapeType>;
    { policy.stride() } -> std::same_as<typename T::IndexType>;
    { policy.shape() } -> std::same_as<typename T::ShapeType>;
  } &&
  ((T::NumDims == 1 && std::convertible_to<typename T::ShapeType, typename T::IndexType>) ||
   (T::NumDims > 1 && requires(const T& policy, const typename T::ShapeType& shape, int dim) {
     { policy.strides() } -> std::same_as<typename T::ShapeType>;
     { shape[dim] } -> std::convertible_to<typename T::IndexType>;
   }));

template <typename T>
concept OrderedSetSubsetPolicy =
  InheritablePolicy<T> && SubsetPolicy<T> && std::default_initializable<T> && std::copyable<T> &&
  std::constructible_from<T, typename T::ParentSetType*>;

template <typename T, typename Iterator>
concept ValidatesSubset = requires(const T& policy, Iterator begin, Iterator end) {
  { policy.isValid(begin, end, false) } -> std::convertible_to<bool>;
};

template <typename Position, typename Element, typename Size, typename Offset, typename Stride, typename Indirection, typename Subset>
concept OrderedSetPoliciesFor = SetSizePolicyFor<Size, Position> &&
  OrderedSetOffsetPolicyFor<Offset, Position> && OrderedSetStridePolicyFor<Stride, Position> &&
  OrderedSetIndirectionPolicyFor<Indirection, Position, Element> && OrderedSetSubsetPolicy<Subset>;

template <typename Data, typename Set, typename Indirection, typename Stride>
concept MapParameters = SetLike<Set> && MapStridePolicyFor<Stride, typename Set::PositionType> &&
  MapIndirectionPolicyFor<Indirection, typename Set::PositionType, Data>;
}  // namespace detail

}  // namespace axom::slam
