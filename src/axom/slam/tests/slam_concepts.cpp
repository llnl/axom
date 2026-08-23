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
#include "axom/slam/Traits.hpp"
#include "axom/slam/policies/CardinalityPolicies.hpp"
#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/policies/OffsetPolicies.hpp"
#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/policies/StridePolicies.hpp"

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
using MismatchedPolicyMap = slam::Map<double, ConcreteRange, WrongDataIndirection>;
using VariableRelation = slam::VariableRelationView<ConcreteRange, ConcreteRange>;
using DynamicVariableRelation =
  slam::DynamicVariableRelation<ConcreteRange, ConcreteRange>;
using NarrowRange = slam::RangeSet<std::int32_t, std::int32_t>;
using WideRange = slam::RangeSet<std::int64_t, std::int64_t>;
using HeterogeneousVariableRelation =
  slam::VariableRelationView<WideRange, NarrowRange>;
using HeterogeneousDynamicVariableRelation =
  slam::DynamicVariableRelation<NarrowRange, WideRange>;
using DynamicConstantCardinality =
  policies::ConstantCardinality<Position,
                                policies::CompileTimeStride<Position, 3>>;
using DynamicConstantRelation = slam::DynamicConstantRelation<Position,
                                                              Element,
                                                              DynamicConstantCardinality>;
using DistinctElementDynamicConstantRelation =
  slam::DynamicConstantRelation<Position, std::int64_t, DynamicConstantCardinality>;

struct StrongPosition
{
  std::int64_t value;
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
  using PositionType = Position;
  using ElementType = Element;
  using FirstSetType = ConcreteRange;
  using SecondSetType = ConcreteRange;
};

struct TypedefOnlyRelation
{
  using FromSetType = ConcreteRange;
  using ToSetType = ConcreteRange;
  using SetPosition = Position;
  using SetElement = Element;
};

struct WrongSetPositionRelation
{
  using FromSetType = ConcreteRange;
  using ToSetType = ConcreteRange;
  using SetPosition = short;
  using SetElement = Element;

  const FromSetType* fromSet() const;
  const ToSetType* toSet() const;
  int operator[](SetPosition) const;
  int* begin(SetPosition) const;
  int* end(SetPosition) const;
};

struct WrongSetElementRelation
{
  using FromSetType = ConcreteRange;
  using ToSetType = WideRange;
  using SetPosition = typename FromSetType::PositionType;
  using SetElement = short;

  const FromSetType* fromSet() const;
  const ToSetType* toSet() const;
  int operator[](SetPosition) const;
  SetElement* begin(SetPosition) const;
  SetElement* end(SetPosition) const;
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
inline constexpr bool
  axom::slam::enable_position_like<slam_concept_test::StrongPosition> = true;

namespace slam_concept_test
{
// Sets and bivariate sets
static_assert(slam::SetLike<slam::Set<>>);
static_assert(slam::SetLike<Range>);
static_assert(slam::SetLike<const Range&>);
static_assert(slam::OrderedSetLike<Range>);
static_assert(!slam::BivariateSetLike<Range>);
static_assert(slam::BivariateSetLike<Product>);
static_assert(!slam::SetLike<Product>);
static_assert(!slam::SetLike<TypedefOnlySet>);
static_assert(!slam::SetLike<WrongSizeSet>);
static_assert(!slam::SetLike<MutableAccessOnlySet>);
static_assert(!slam::SetLike<FloatingPositionSet>);
static_assert(!slam::BivariateSetLike<TypedefOnlyBivariateSet>);
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
static_assert(slam::RelationLike<DynamicConstantRelation>);
static_assert(slam::RelationLike<DistinctElementDynamicConstantRelation>);
static_assert(slam::RelationLike<const DynamicConstantRelation&>);
static_assert(slam::is_relation_like_v<DynamicConstantRelation>);
static_assert(!slam::RelationLike<TypedefOnlyRelation>);
static_assert(!slam::RelationLike<WrongSetPositionRelation>);
static_assert(!slam::RelationLike<WrongSetElementRelation>);

// Maps
// BivariateMap::SetType is its flat backing set.
// Its semantic domain is BivariateSetType, which is the type MapOver checks.
static_assert(slam::UnivariateMapLike<UnaryMap>);
static_assert(!slam::BivariateMapLike<UnaryMap>);
static_assert(slam::BivariateMapLike<BinaryMap>);
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
static_assert(!slam::MapLike<MismatchedPolicyMap>);
static_assert(!slam::MapLike<int>);
static_assert(!slam::MapOver<int, ConcreteRange>);

// Policies
using Size = policies::CompileTimeSize<int, 5>;
using EmptySize = policies::ZeroSize<int>;
using ScalarStride = policies::CompileTimeStride<int, 3>;
using MatrixStride = policies::MultiDimStride<int, 2>;
using Offset = policies::CompileTimeOffset<int, 4>;
using NoIndirection = policies::NoIndirection<Position, Element>;

static_assert(slam::ValuePolicy<Size>);
static_assert(slam::SizePolicy<Size>);
static_assert(slam::SizePolicy<EmptySize>);
static_assert(!slam::ValuePolicy<EmptySize>);
static_assert(slam::StridePolicy<ScalarStride>);
static_assert(slam::StridePolicy<MatrixStride>);
static_assert(slam::OffsetPolicy<Offset>);
static_assert(!slam::ValuePolicy<TypedefOnlyValuePolicy>);
static_assert(!slam::ValuePolicy<WrongValuePolicy>);
static_assert(!slam::SizePolicy<int>);
static_assert(!slam::StridePolicy<int>);
static_assert(!slam::OffsetPolicy<int>);
static_assert(slam::IndirectionPolicy<ViewIndirection>);
static_assert(slam::IndirectionPolicyFor<ViewIndirection, Position>);
static_assert(slam::IndirectionPolicy<NoIndirection>);
static_assert(slam::IndirectionPolicyFor<NoIndirection, Position>);
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
