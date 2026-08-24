// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file slam_concepts.cpp
 *
 * \brief Positive and negative compile-time tests for SLAM's C++20 concepts.
 */

#include "gtest/gtest.h"

#include "axom/slam/Aliases.hpp"
#include "axom/slam/BivariateMap.hpp"
#include "axom/slam/Concepts.hpp"
#include "axom/slam/DynamicConstantRelation.hpp"
#include "axom/slam/DynamicMap.hpp"
#include "axom/slam/DynamicVariableRelation.hpp"
#include "axom/slam/Map.hpp"
#include "axom/slam/ProductSet.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/slam/RelationSet.hpp"
#include "axom/slam/Traits.hpp"
#include "axom/slam/policies/CardinalityPolicies.hpp"
#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/policies/OffsetPolicies.hpp"
#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/policies/StridePolicies.hpp"

#include <cstddef>
#include <cstdint>

namespace slam_concept_test
{
namespace slam = axom::slam;
namespace policies = axom::slam::policies;

using Position = slam::DefaultPositionType;
using Element = slam::DefaultElementType;
using Range = slam::RangeSet<Position, Element>;
using ConcreteRange = typename Range::ConcreteSet;
using Product = typename slam::ProductSet<ConcreteRange, ConcreteRange>::ConcreteSet;
using ViewIndirection = policies::ArrayViewIndirection<Position, double>;
using UnaryMap = slam::Map<double, ConcreteRange, ViewIndirection>;
using BinaryMap = slam::BivariateMap<double, Product, ViewIndirection>;
using DynamicSet = slam::DynamicSet<Position, Element>;
using DynamicMap = slam::DynamicMap<DynamicSet, double>;
using WrongDataIndirection = policies::ArrayViewIndirection<Position, int>;
using VariableRelation = slam::VariableRelationView<ConcreteRange, ConcreteRange>;
using DynamicVariableRelation = slam::DynamicVariableRelation<ConcreteRange, ConcreteRange>;
using NarrowRange = slam::RangeSet<std::int32_t, std::int32_t>;
using WideRange = slam::RangeSet<std::int64_t, std::int64_t>;
using FirstHandleRange = slam::RangeSet<std::int32_t, double>;
using SecondHandleRange = slam::RangeSet<std::int32_t, float>;
using DistinctHandleProduct =
  typename slam::ProductSet<FirstHandleRange, SecondHandleRange>::ConcreteSet;
using HeterogeneousProduct = typename slam::ProductSet<NarrowRange, WideRange>::ConcreteSet;
using HeterogeneousNullBivariateSet = slam::NullBivariateSet<NarrowRange, WideRange>;
using HeterogeneousMapIndirection = policies::ArrayViewIndirection<std::int64_t, double>;
using HeterogeneousBinaryMap =
  slam::BivariateMap<double, HeterogeneousProduct, HeterogeneousMapIndirection>;
using HeterogeneousVariableRelation = slam::VariableRelationView<WideRange, NarrowRange>;
using HeterogeneousDynamicVariableRelation = slam::DynamicVariableRelation<NarrowRange, WideRange>;
using DynamicConstantCardinality =
  policies::ConstantCardinality<Position, policies::CompileTimeStride<Position, 3>>;
using DynamicConstantRelation =
  slam::DynamicConstantRelation<Position, Element, DynamicConstantCardinality>;
using DistinctElementDynamicConstantRelation =
  slam::DynamicConstantRelation<Position, std::int64_t, DynamicConstantCardinality>;

template <typename FirstSet, typename SecondSet>
concept CanFormProductSet = requires { sizeof(slam::ProductSet<FirstSet, SecondSet>); };

template <typename RelationType>
concept CanFormRelationSet = requires {
  typename std::type_identity_t<
    slam::RelationSet<RelationType, typename RelationType::FromSetType, typename RelationType::ToSetType>>;
};

struct StrongPosition
{
  std::int64_t value;
};

struct ExplicitFirstPosition
{
  explicit ExplicitFirstPosition(std::int32_t value = 0) : value(value) { }
  std::int32_t value;
};

struct ExplicitSecondPosition
{
  explicit ExplicitSecondPosition(std::int64_t value = 0) : value(value) { }
  std::int64_t value;
};

struct ExplicitFlatPosition
{
  explicit ExplicitFlatPosition(std::int64_t value = 0) : value(value) { }
  explicit ExplicitFlatPosition(ExplicitFirstPosition position) : value(position.value) { }
  explicit ExplicitFlatPosition(ExplicitSecondPosition position) : value(position.value) { }
  std::int64_t value;
};

struct ExplicitRowPosition
{
  explicit ExplicitRowPosition(std::int32_t value = 0) : value(value) { }
  std::int32_t value;
};

template <typename PositionType_, typename ElementType_>
struct MinimalSet
{
  using PositionType = PositionType_;
  using ElementType = ElementType_;

  PositionType size() const;
  bool empty() const;
  ElementType at(PositionType) const;
};

using ExplicitFirstSet = MinimalSet<ExplicitFirstPosition, double>;
using ExplicitSecondSet = MinimalSet<ExplicitSecondPosition, float>;

struct ExplicitPositionRow
{
  ExplicitRowPosition size() const;
  const ExplicitSecondPosition* begin() const;
  const ExplicitSecondPosition* end() const;
};

struct ExplicitPositionBivariateSet
{
  using FirstSetType = ExplicitFirstSet;
  using SecondSetType = ExplicitSecondSet;
  using PositionType = ExplicitFlatPosition;
  using ElementType = std::pair<ExplicitFirstPosition, ExplicitSecondPosition>;

  PositionType size() const;
  ElementType at(PositionType) const;
  const FirstSetType* getFirstSet() const;
  const SecondSetType* getSecondSet() const;
  ExplicitPositionRow getElements(ExplicitFirstPosition) const;
  ExplicitRowPosition findElementIndex(ExplicitFirstPosition, ExplicitSecondPosition) const;
  PositionType findElementFlatIndex(ExplicitFirstPosition, ExplicitSecondPosition) const;
  ExplicitFirstPosition flatToFirstIndex(PositionType) const;
  ExplicitSecondPosition flatToSecondIndex(PositionType) const;
};

// Associated types alone must not satisfy a semantic concept.
struct TypedefOnlySet
{
  using PositionType = Position;
  using ElementType = Element;
};

struct WrongSizeSet
{
  using PositionType = Position;
  using ElementType = Element;

  double size() const;
  bool empty() const;
  Element at(Position) const;
};

struct MutableAccessOnlySet
{
  using PositionType = Position;
  using ElementType = Element;

  Position size() const;
  bool empty() const;
  Element at(Position);
};

struct FloatingPositionSet
{
  using PositionType = double;
  using ElementType = Element;

  double size() const;
  bool empty() const;
  Element at(double) const;
};

struct TypedefOnlyBivariateSet
{
  using FirstSetType = ConcreteRange;
  using SecondSetType = ConcreteRange;
  using PositionType = Position;
  using ElementType =
    std::pair<typename FirstSetType::PositionType, typename SecondSetType::PositionType>;
};

// A valid model needs no aliases that duplicate endpoint or set types.
struct MinimalCoordinate
{
  NarrowRange::PositionType first;
  WideRange::PositionType second;
};

struct MinimalBivariateSet
{
  using FirstSetType = NarrowRange;
  using SecondSetType = WideRange;
  using PositionType = std::int64_t;
  using ElementType = MinimalCoordinate;

  PositionType size() const;
  ElementType at(PositionType) const;
  const FirstSetType* getFirstSet() const;
  const SecondSetType* getSecondSet() const;
  SecondSetType getElements(FirstSetType::PositionType) const;
  PositionType findElementIndex(FirstSetType::PositionType, SecondSetType::PositionType) const;
  PositionType findElementFlatIndex(FirstSetType::PositionType, SecondSetType::PositionType) const;
  FirstSetType::PositionType flatToFirstIndex(PositionType) const;
  SecondSetType::PositionType flatToSecondIndex(PositionType) const;
};

struct HeterogeneousPositionBivariateSet
{
  using PositionType = typename NarrowRange::PositionType;
  using ElementType = typename WideRange::PositionType;
  using FirstSetType = NarrowRange;
  using SecondSetType = WideRange;

  PositionType size() const;
  ElementType at(PositionType) const;
  const FirstSetType* getFirstSet() const;
  const SecondSetType* getSecondSet() const;
  PositionType findElementFlatIndex(PositionType, PositionType) const;
  PositionType flatToFirstIndex(PositionType) const;
  PositionType flatToSecondIndex(PositionType) const;
};

struct WrongElementBivariateSet
{
  using FirstSetType = NarrowRange;
  using SecondSetType = NarrowRange;
  using FirstPositionType = typename FirstSetType::PositionType;
  using SecondPositionType = typename SecondSetType::PositionType;
  using PositionType = std::int32_t;
  using ElementType = double;

  PositionType size() const;
  PositionType size(FirstPositionType) const;
  ElementType at(PositionType) const;
  const FirstSetType* getFirstSet() const;
  const SecondSetType* getSecondSet() const;
  FirstPositionType firstSetSize() const;
  SecondPositionType secondSetSize() const;
  SecondSetType getElements(FirstPositionType) const;
  PositionType findElementIndex(FirstPositionType, SecondPositionType) const;
  PositionType findElementFlatIndex(FirstPositionType, SecondPositionType) const;
  FirstPositionType flatToFirstIndex(PositionType) const;
  SecondPositionType flatToSecondIndex(PositionType) const;
};

struct WrongCoordinateBivariateSet
{
  using FirstSetType = NarrowRange;
  using SecondSetType = WideRange;
  using FirstPositionType = typename FirstSetType::PositionType;
  using SecondPositionType = typename SecondSetType::PositionType;
  using PositionType = std::int64_t;
  using ElementType = std::pair<SecondPositionType, SecondPositionType>;

  PositionType size() const;
  PositionType size(FirstPositionType) const;
  ElementType at(PositionType) const;
  const FirstSetType* getFirstSet() const;
  const SecondSetType* getSecondSet() const;
  FirstPositionType firstSetSize() const;
  SecondPositionType secondSetSize() const;
  SecondSetType getElements(FirstPositionType) const;
  PositionType findElementIndex(FirstPositionType, SecondPositionType) const;
  PositionType findElementFlatIndex(FirstPositionType, SecondPositionType) const;
  FirstPositionType flatToFirstIndex(PositionType) const;
  SecondPositionType flatToSecondIndex(PositionType) const;
};

struct TypedefOnlyRelation
{
  using FromSetType = ConcreteRange;
  using ToSetType = ConcreteRange;
};

struct MinimalRelationRow
{
  std::size_t size() const;
  WideRange::PositionType* begin() const;
  WideRange::PositionType* end() const;
};

struct MinimalRelation
{
  using FromSetType = NarrowRange;
  using ToSetType = WideRange;

  const FromSetType* fromSet() const;
  const ToSetType* toSet() const;
  MinimalRelationRow operator[](FromSetType::PositionType) const;
};

struct NotAPosition
{ };

struct WrongRelationRow
{
  std::size_t size() const;
  NotAPosition* begin() const;
  NotAPosition* end() const;
};

struct WrongRelationEntry
{
  using FromSetType = ConcreteRange;
  using ToSetType = WideRange;

  const FromSetType* fromSet() const;
  const ToSetType* toSet() const;
  WrongRelationRow operator[](FromSetType::PositionType) const;
};

struct TypedefOnlyMap
{
  using DataType = double;
  using SetType = ConcreteRange;
  using SetPosition = Position;
  using SetElement = Element;
  using ValueType = double&;
  using ConstValueType = const double&;
};

struct MinimalBivariateMap
{
  using DataType = double;
  using BivariateSetType = MinimalBivariateSet;
  using SetPosition = BivariateSetType::PositionType;
  using SetElement = BivariateSetType::ElementType;
  using ValueType = double&;
  using ConstValueType = const double&;

  SetPosition size() const;
  ValueType operator[](SetPosition);
  ConstValueType operator[](SetPosition) const;
  const BivariateSetType* set() const;
};

struct WrongDomainMap
{
  using DataType = double;
  using SetType = ConcreteRange;
  using SetPosition = Position;
  using SetElement = Element;
  using ValueType = double&;
  using ConstValueType = const double&;

  SetPosition size() const;
  ConstValueType operator[](SetPosition) const;
  const TypedefOnlySet* set() const;
};

struct WrongPositionMap
{
  using DataType = double;
  using SetType = ConcreteRange;
  using SetPosition = short;
  using SetElement = Element;
  using ValueType = double&;
  using ConstValueType = const double&;

  SetPosition size() const;
  ConstValueType operator[](SetPosition) const;
  const SetType* set() const;
};

struct WrongValueTypeMap
{
  using DataType = double;
  using SetType = ConcreteRange;
  using SetPosition = Position;
  using SetElement = Element;
  using ValueType = int&;
  using ConstValueType = const int&;

  SetPosition size() const;
  ValueType operator[](SetPosition);
  ConstValueType operator[](SetPosition) const;
  const SetType* set() const;
};

struct WrongMutableAccessMap
{
  using DataType = double;
  using SetType = ConcreteRange;
  using SetPosition = Position;
  using SetElement = Element;
  using ValueType = double&;
  using ConstValueType = const double&;

  SetPosition size() const;
  ConstValueType operator[](SetPosition);
  ConstValueType operator[](SetPosition) const;
  const SetType* set() const;
};

struct TypedefOnlyValuePolicy
{
  struct TagType;
  using IntType = int;
};

struct WrongValuePolicy
{
  struct TagType;
  using IntType = int;

  double value() const;
  bool isValid(bool) const;
};

struct TypedefOnlyIndirection
{
  using IndirectionResult = double&;
  using ConstIndirectionResult = const double&;
  using IndirectionBufferType = double*;
  using IndirectionPtrType = double*;

  static constexpr bool DeviceAccessible = true;
};

struct TrivialCapture
{
  int value;
};

struct NonTrivialCapture
{
  NonTrivialCapture(const NonTrivialCapture&) { }
};

}  // namespace slam_concept_test

template <>
inline constexpr bool axom::slam::enable_position_like<slam_concept_test::StrongPosition> = true;

template <>
inline constexpr bool axom::slam::enable_position_like<slam_concept_test::ExplicitFirstPosition> =
  true;

template <>
inline constexpr bool axom::slam::enable_position_like<slam_concept_test::ExplicitSecondPosition> =
  true;

template <>
inline constexpr bool axom::slam::enable_position_like<slam_concept_test::ExplicitFlatPosition> =
  true;

template <>
inline constexpr bool axom::slam::enable_position_like<slam_concept_test::ExplicitRowPosition> =
  true;

namespace slam_concept_test
{
// Sets and bivariate sets
static_assert(slam::SetLike<slam::Set<>>);
static_assert(slam::SetLike<Range>);
static_assert(slam::SetLike<const Range&>);
static_assert(slam::OrderedSetLike<Range>);
static_assert(!slam::BivariateSetLike<Range>);
static_assert(slam::BivariateSetLike<Product>);
static_assert(std::same_as<typename Product::ElementType, std::pair<Position, Position>>);
static_assert(slam::BivariateSetLike<DistinctHandleProduct>);
static_assert(
  std::same_as<typename DistinctHandleProduct::ElementType, std::pair<std::int32_t, std::int32_t>>);
static_assert(slam::BivariateSetLike<HeterogeneousProduct>);
static_assert(slam::BivariateSetLike<HeterogeneousNullBivariateSet>);
static_assert(slam::BivariateSetLike<MinimalBivariateSet>);
static_assert(
  !std::same_as<typename MinimalBivariateSet::ElementType, std::pair<std::int32_t, std::int64_t>>);
static_assert(slam::BivariateSetLike<ExplicitPositionBivariateSet>);
static_assert(!std::convertible_to<ExplicitFirstPosition, ExplicitFlatPosition>);
static_assert(!std::convertible_to<ExplicitSecondPosition, ExplicitFlatPosition>);
static_assert(!std::convertible_to<ExplicitRowPosition, ExplicitFlatPosition>);
static_assert(std::same_as<typename HeterogeneousProduct::FirstPositionType, std::int32_t>);
static_assert(std::same_as<typename HeterogeneousProduct::SecondPositionType, std::int64_t>);
static_assert(std::same_as<typename HeterogeneousProduct::PositionType, std::int64_t>);
static_assert(
  std::same_as<typename HeterogeneousProduct::ElementType, std::pair<std::int32_t, std::int64_t>>);
static_assert(!slam::SetLike<Product>);
static_assert(!slam::SetLike<TypedefOnlySet>);
static_assert(!slam::SetLike<WrongSizeSet>);
static_assert(!slam::SetLike<MutableAccessOnlySet>);
static_assert(!slam::SetLike<FloatingPositionSet>);
static_assert(!slam::BivariateSetLike<TypedefOnlyBivariateSet>);
static_assert(!slam::BivariateSetLike<HeterogeneousPositionBivariateSet>);
static_assert(!slam::BivariateSetLike<WrongElementBivariateSet>);
static_assert(!slam::BivariateSetLike<WrongCoordinateBivariateSet>);
static_assert(CanFormProductSet<NarrowRange, NarrowRange>);
static_assert(CanFormProductSet<NarrowRange, WideRange>);
static_assert(!slam::SetLike<int>);

// Relations
static_assert(slam::RelationLike<VariableRelation>);
static_assert(slam::RelationLike<HeterogeneousVariableRelation>);
static_assert(slam::RelationLike<DynamicVariableRelation>);
static_assert(slam::RelationLike<HeterogeneousDynamicVariableRelation>);
static_assert(std::same_as<typename HeterogeneousDynamicVariableRelation::SetPosition,
                           typename NarrowRange::PositionType>);
static_assert(std::same_as<typename HeterogeneousDynamicVariableRelation::SetElement,
                           typename WideRange::PositionType>);
static_assert(CanFormRelationSet<VariableRelation>);
static_assert(CanFormRelationSet<HeterogeneousVariableRelation>);
static_assert(slam::RelationLike<DynamicConstantRelation>);
static_assert(slam::RelationLike<DistinctElementDynamicConstantRelation>);
static_assert(slam::RelationLike<const DynamicConstantRelation&>);
static_assert(slam::is_relation_like_v<DynamicConstantRelation>);
static_assert(slam::RelationLike<MinimalRelation>);
static_assert(!slam::RelationLike<TypedefOnlyRelation>);
static_assert(!slam::RelationLike<WrongRelationEntry>);

// Maps
// BivariateMap::SetType is its flat backing set.
// Its semantic domain is BivariateSetType, which is the type MapOver checks.
static_assert(slam::UnivariateMapLike<UnaryMap>);
static_assert(!slam::BivariateMapLike<UnaryMap>);
static_assert(slam::BivariateMapLike<BinaryMap>);
static_assert(slam::BivariateMapLike<HeterogeneousBinaryMap>);
static_assert(slam::BivariateMapLike<MinimalBivariateMap>);
static_assert(std::same_as<typename HeterogeneousBinaryMap::FirstPositionType, std::int32_t>);
static_assert(std::same_as<typename HeterogeneousBinaryMap::SecondPositionType, std::int64_t>);
static_assert(std::same_as<decltype(HeterogeneousBinaryMap::iterator {}.firstIndex()), std::int32_t>);
static_assert(std::same_as<decltype(HeterogeneousBinaryMap::iterator {}.secondIndex()), std::int64_t>);
static_assert(
  std::same_as<decltype(HeterogeneousBinaryMap::range_iterator {}.firstIndex()), std::int32_t>);
static_assert(
  std::same_as<decltype(HeterogeneousBinaryMap::range_iterator {}.secondIndex()), std::int64_t>);
static_assert(
  std::same_as<typename HeterogeneousBinaryMap::SetElement, std::pair<std::int32_t, std::int64_t>>);
static_assert(slam::MapLike<UnaryMap>);
static_assert(slam::MapLike<const UnaryMap&>);
static_assert(slam::MapLike<BinaryMap>);
static_assert(slam::UnivariateMapLike<DynamicMap>);
static_assert(slam::MapLike<DynamicMap>);
static_assert(slam::MapLike<const DynamicMap&>);
static_assert(slam::MapOver<UnaryMap, ConcreteRange>);
static_assert(slam::MapOver<BinaryMap, Product>);
static_assert(slam::MapOver<DynamicMap, DynamicSet>);
static_assert(!slam::MapOver<BinaryMap, typename BinaryMap::SetType>);
static_assert(!slam::MapLike<TypedefOnlyMap>);
static_assert(!slam::MapLike<WrongDomainMap>);
static_assert(!slam::MapLike<WrongPositionMap>);
static_assert(!slam::MapLike<WrongValueTypeMap>);
static_assert(!slam::MapLike<WrongMutableAccessMap>);
static_assert(!slam::MapLike<int>);
static_assert(!slam::MapOver<int, ConcreteRange>);

// Policies
using Size = policies::CompileTimeSize<int, 5>;
using EmptySize = policies::ZeroSize<int>;
using ScalarStride = policies::CompileTimeStride<int, 3>;
using MatrixStride = policies::MultiDimStride<int, 2>;
using Offset = policies::CompileTimeOffset<int, 4>;
using NoIndirection = policies::NoIndirection<Position, Element>;
using OwningIndirection = policies::ArrayIndirection<Position, double>;
using OwningMap = slam::Map<double, ConcreteRange, OwningIndirection>;

static_assert(slam::ValuePolicy<Size>);
static_assert(slam::SizePolicy<Size>);
static_assert(slam::SizePolicy<EmptySize>);
static_assert(!slam::ValuePolicy<EmptySize>);
static_assert(slam::StridePolicy<ScalarStride>);
static_assert(slam::StridePolicy<MatrixStride>);
static_assert(slam::OrderedSetStridePolicyFor<ScalarStride, int>);
static_assert(!slam::OrderedSetStridePolicyFor<MatrixStride, int>);
static_assert(slam::MapStridePolicyFor<ScalarStride, int>);
static_assert(slam::MapStridePolicyFor<MatrixStride, int>);
// Maps support a stride index that converts to their (possibly wider) position type
// OrderedSet's scalar value policy must use that exact position type.
static_assert(slam::MapStridePolicyFor<ScalarStride, Position>);
static_assert(!slam::OrderedSetStridePolicyFor<ScalarStride, Position>);
static_assert(slam::OffsetPolicy<Offset>);
static_assert(!slam::ValuePolicy<TypedefOnlyValuePolicy>);
static_assert(!slam::ValuePolicy<WrongValuePolicy>);
static_assert(!slam::SizePolicy<int>);
static_assert(!slam::StridePolicy<int>);
static_assert(!slam::OffsetPolicy<int>);
static_assert(slam::IndirectionPolicy<ViewIndirection>);
static_assert(slam::IndirectionPolicyFor<ViewIndirection, Position>);
static_assert(slam::OrderedSetIndirectionPolicyFor<ViewIndirection, Position, double>);
static_assert(slam::MapIndirectionPolicyFor<ViewIndirection, Position, double>);
static_assert(!slam::AllocatingMapIndirectionPolicyFor<ViewIndirection, Position, double>);
static_assert(slam::IndirectionPolicy<NoIndirection>);
static_assert(slam::IndirectionPolicyFor<NoIndirection, Position>);
static_assert(slam::OrderedSetIndirectionPolicyFor<NoIndirection, Position, Element>);
static_assert(!slam::MapIndirectionPolicyFor<NoIndirection, Position, Element>);
static_assert(slam::MapIndirectionPolicyFor<OwningIndirection, Position, double>);
static_assert(slam::AllocatingMapIndirectionPolicyFor<OwningIndirection, Position, double>);
static_assert(std::default_initializable<OwningMap>);
static_assert(!std::default_initializable<UnaryMap>);
static_assert(std::constructible_from<UnaryMap, const ConcreteRange*, typename UnaryMap::OrderedMap>);
static_assert(!slam::MapIndirectionPolicyFor<WrongDataIndirection, Position, double>);
static_assert(!slam::IndirectionPolicy<TypedefOnlyIndirection>);
static_assert(!slam::IndirectionPolicy<int>);

// Index and deployment properties
static_assert(slam::PositionLike<int>);
static_assert(!slam::PositionLike<double>);
static_assert(slam::PositionLike<StrongPosition>);
static_assert(slam::DeviceCapturable<TrivialCapture>);
static_assert(!slam::DeviceCapturable<NonTrivialCapture>);
static_assert(!slam::DeviceCapturable<TrivialCapture&>);

// Compatibility trait spellings remain exact Boolean wrappers around the concepts.
static_assert(slam::is_set_like_v<Range> == slam::SetLike<Range>);
static_assert(slam::is_bivariate_set_like_v<Product> == slam::BivariateSetLike<Product>);
static_assert(slam::is_relation_like_v<VariableRelation> == slam::RelationLike<VariableRelation>);
static_assert(slam::is_map_like_v<BinaryMap> == slam::MapLike<BinaryMap>);
static_assert(slam::is_map_over_v<BinaryMap, Product> == slam::MapOver<BinaryMap, Product>);
static_assert(slam::is_size_policy_v<EmptySize> == slam::SizePolicy<EmptySize>);
static_assert(slam::is_stride_policy_v<MatrixStride> == slam::StridePolicy<MatrixStride>);
static_assert(slam::is_indirection_policy_v<ViewIndirection> ==
              slam::IndirectionPolicy<ViewIndirection>);
static_assert(slam::is_position_like_v<StrongPosition> == slam::PositionLike<StrongPosition>);
static_assert(slam::is_handle_like_v<TrivialCapture> == slam::DeviceCapturable<TrivialCapture>);
static_assert(!slam::is_set_like_v<int>);
static_assert(!slam::is_map_over_v<int, ConcreteRange>);

TEST(slam_concepts, compile_time_contracts) { SUCCEED(); }
}  // namespace slam_concept_test
