// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file slam_concepts.cpp
 *
 * \brief Positive and negative compile-time tests for Slam's C++20 concepts.
 */

#include "axom/slam/Concepts.hpp"

#include "axom/config.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/policies/OffsetPolicies.hpp"
#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/policies/StridePolicies.hpp"
#include "axom/slam/policies/SubsettingPolicies.hpp"

#include "gtest/gtest.h"

#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <vector>

namespace slam_concept_test
{
namespace slam = axom::slam;
namespace policies = axom::slam::policies;

using Position = slam::DefaultPositionType;
using Element = slam::DefaultElementType;
using ConcreteRange = slam::RangeSet<Position, Element>::ConcreteSet;
using NarrowRange = slam::RangeSet<std::int32_t, std::int32_t>;
using WideRange = slam::RangeSet<std::int64_t, std::int64_t>;
using ViewIndirection = policies::ArrayViewIndirection<Position, double>;
using ConstViewIndirection = policies::ArrayViewIndirection<Position, const double>;
using WrongDataIndirection = policies::ArrayViewIndirection<Position, int>;
using RawPointerMapIndirection = policies::CArrayIndirection<Position, double>;
using VectorMapIndirection = policies::STLVectorIndirection<Position, double>;

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

  PositionType size() const { return 3; }
  bool empty() const { return false; }
  ElementType at(PositionType pos) const { return {pos == 0 ? 0 : 1, pos}; }
  const FirstSetType* getFirstSet() const { return &first; }
  const SecondSetType* getSecondSet() const { return &second; }
  SecondSetType getElements(FirstSetType::PositionType pos) const
  {
    return pos == 0 ? SecondSetType(0, 1) : SecondSetType(1, 3);
  }

  FirstSetType first {2};
  SecondSetType second {3};
};

struct BivariateSetMissingEmpty
{
  using FirstSetType = NarrowRange;
  using SecondSetType = WideRange;
  using PositionType = std::int64_t;
  using ElementType = MinimalCoordinate;

  PositionType size() const;
  typename FirstSetType::PositionType firstSetSize() const;
  typename SecondSetType::PositionType secondSetSize() const;
  PositionType size(typename FirstSetType::PositionType) const;
  ElementType at(PositionType) const;
  const FirstSetType* getFirstSet() const;
  const SecondSetType* getSecondSet() const;
  SecondSetType getElements(FirstSetType::PositionType) const;
  PositionType findElementIndex(FirstSetType::PositionType, SecondSetType::PositionType) const;
  PositionType findElementFlatIndex(FirstSetType::PositionType, SecondSetType::PositionType) const;
  FirstSetType::PositionType flatToFirstIndex(PositionType) const;
  SecondSetType::PositionType flatToSecondIndex(PositionType) const;
};

struct TypedefOnlyRelation
{
  using FromSetType = ConcreteRange;
  using ToSetType = ConcreteRange;
};

struct MinimalRelationRow
{
  std::size_t size() const { return count; }
  const WideRange::PositionType* begin() const { return data; }
  const WideRange::PositionType* end() const { return data + count; }

  const WideRange::PositionType* data;
  std::size_t count;
};

struct MinimalRelation
{
  using FromSetType = NarrowRange;
  using ToSetType = WideRange;

  const FromSetType* fromSet() const { return &from; }
  const ToSetType* toSet() const { return &to; }
  MinimalRelationRow operator[](FromSetType::PositionType pos) const
  {
    return pos == 0 ? MinimalRelationRow {indices, 1} : MinimalRelationRow {indices + 1, 2};
  }

  FromSetType from {2};
  ToSetType to {3};
  WideRange::PositionType indices[3] {2, 0, 1};
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
  using PositionType = Position;
  using SetElement = Element;
  using ValueType = double&;
  using ConstValueType = const double&;
};

struct MinimalBivariateMap
{
  using DataType = double;
  using PositionType = std::int64_t;

  PositionType size() const { return 2; }
  PositionType numComp() const { return 3; }
  MinimalCoordinate index(PositionType pos) const { return {1, 2 * pos}; }
  double& flatValue(PositionType pos, PositionType component)
  {
    return values[3 * pos + component];
  }
  const double& flatValue(PositionType pos, PositionType component) const
  {
    return values[3 * pos + component];
  }

  double values[6] {};
};

struct WrongDomainMap : MinimalBivariateMap
{
  using MappedSetType = MinimalBivariateSet;
  const TypedefOnlySet* set() const;
};

struct WrongPositionMap : MinimalBivariateMap
{
  using PositionType = short;
};

struct WrongValueTypeMap : MinimalBivariateMap
{
  int& flatValue(PositionType, PositionType);
  const int& flatValue(PositionType, PositionType) const;
};

struct WrongMutableAccessMap : MinimalBivariateMap
{
  double flatValue(PositionType, PositionType);
  const double& flatValue(PositionType, PositionType) const;
};

struct DropsDataConstnessMap : MinimalBivariateMap
{
  using DataType = const double;
};

struct MissingComponentCountMap : MinimalBivariateMap
{
  void numComp() const = delete;
};

struct MissingIndexMap : MinimalBivariateMap
{
  void index(PositionType) const;
};

struct MissingConstAccessMap : MinimalBivariateMap
{
  using MinimalBivariateMap::flatValue;
  const double& flatValue(PositionType, PositionType) const = delete;
};

struct ExtraAliasesMap : MinimalBivariateMap
{
  using SetType = TypedefOnlySet;
  using BivariateSetType = void;
  using ValueType = int;
};

struct TypedefOnlyIndirection
{
  using IndirectionResult = double&;
  using ConstIndirectionResult = const double&;
  using IndirectionBufferType = double*;
  using IndirectionPtrType = double*;

  static constexpr bool DeviceAccessible = true;
};

struct PrvalueMapIndirection : ViewIndirection
{
  using IndirectionResult = double;
  using ConstIndirectionResult = const double&;
  using ResultPtr = double*;
  using ConstResultPtr = const double*;
};

struct MismatchedConstPointerIndirection : ViewIndirection
{
  using IndirectionResult = double&;
  using ConstIndirectionResult = double&;
  using ResultPtr = double*;
  using ConstResultPtr = const double*;
};

struct DropsReferentConst : ViewIndirection
{
  using ElementType = const double;
};

struct AccessDropsReferentConst : ConstViewIndirection
{
  double& indirection(Position);
  double& indirection(Position) const;
};

struct LeanMapIndirection
{
  using IndirectionResult = double&;
  using ConstIndirectionResult = const double&;
  using IndirectionBufferType = std::vector<double>;
  using ResultPtr = double*;
  using ConstResultPtr = const double*;

  static constexpr bool IsMutableBuffer = true;

  static ResultPtr getIndirection(IndirectionBufferType& buffer, Position pos = 0)
  {
    return buffer.data() + pos;
  }
  static ConstResultPtr getConstIndirection(const IndirectionBufferType& buffer, Position pos = 0)
  {
    return buffer.data() + pos;
  }
};

struct PositionedAccessOnlyIndirection : LeanMapIndirection
{
  static ResultPtr getIndirection(IndirectionBufferType&, Position) { return nullptr; }
  static ConstResultPtr getConstIndirection(const IndirectionBufferType&, Position)
  {
    return nullptr;
  }
};

struct StrongPosition
{
  std::int64_t value;
};

struct TrivialCapture
{
  int value;
};

struct NonTrivialCapture
{
  NonTrivialCapture(const NonTrivialCapture&) { }
};

struct ExternalSet
{
  using PositionType = Position;
  using ElementType = int;

  PositionType size() const { return 2; }
  bool empty() const { return false; }
  int at(PositionType pos) const { return pos == 0 ? 10 : 40; }
};

struct NonAdvancingSet : ExternalSet
{
  struct Iterator
  {
    int operator*() const;
    bool operator!=(Iterator) const;
  };
  Iterator begin() const;
  Iterator end() const;
};

struct ExtraSetAliases : ExternalSet
{
  using FirstSetType = void;
  using SecondSetType = void;
};

struct NonComparingSet : ExternalSet
{
  struct Iterator
  {
    int operator*() const;
    Iterator& operator++();
  };
  Iterator begin() const;
  Iterator end() const;
};

struct WrongCoordinateBivariateSet : MinimalBivariateSet
{
  using ElementType = std::pair<std::int64_t, std::int64_t>;
  ElementType at(PositionType) const;
};

struct BoundMap : MinimalBivariateMap
{
  struct Domain
  {
    using PositionType = MinimalBivariateMap::PositionType;
    using ElementType = MinimalCoordinate;

    PositionType size() const { return 2; }
    bool empty() const { return false; }
    ElementType at(PositionType pos) const { return {1, 2 * pos}; }
  };
  using MappedSetType = Domain;
  const Domain* set() const { return &domain; }

  Domain domain;
};

struct ValidatingSet : ExternalSet
{
  bool isValid(bool = false) const { return true; }
};

struct IterableExternalSet : ExternalSet
{
  const int* begin() const { return values; }
  const int* end() const { return values + 2; }

  int values[2] {10, 40};
};

struct FinalIndirection final : ViewIndirection
{ };

struct AbstractIndirection : ViewIndirection
{
  virtual void bind() = 0;
};

// Object queries normalize cv/ref qualification while preserving element qualification.
static_assert(slam::SetLike<ExternalSet>);
static_assert(slam::SetLike<const ExternalSet&>);
static_assert(slam::SetLike<ExtraSetAliases>);
static_assert(slam::SetLike<MinimalBivariateSet>);
static_assert(!slam::SetLike<TypedefOnlySet>);
static_assert(!slam::SetLike<WrongSizeSet>);
static_assert(!slam::SetLike<MutableAccessOnlySet>);
static_assert(!slam::SetLike<FloatingPositionSet>);
static_assert(!slam::SetLike<int>);
static_assert(!slam::SetLike<BivariateSetMissingEmpty>);

static_assert(slam::IterableSetLike<IterableExternalSet>);
static_assert(!slam::IterableSetLike<ExternalSet>);
static_assert(!slam::IterableSetLike<NonAdvancingSet>);
static_assert(!slam::IterableSetLike<NonComparingSet>);
static_assert(!slam::IterableSetLike<MinimalBivariateSet>);

static_assert(slam::Validatable<ValidatingSet>);
static_assert(!slam::Validatable<ExternalSet>);
static_assert(!slam::Validatable<MinimalRelation>);

static_assert(slam::BivariateSetLike<const MinimalBivariateSet&>);
static_assert(!slam::BivariateSetLike<ExtraSetAliases>);
static_assert(!slam::BivariateSetLike<BivariateSetMissingEmpty>);
static_assert(!slam::BivariateSetLike<WrongCoordinateBivariateSet>);

static_assert(slam::RelationLike<MinimalRelation>);
static_assert(slam::RelationLike<const MinimalRelation&>);
static_assert(!slam::RelationLike<TypedefOnlyRelation>);
static_assert(!slam::RelationLike<WrongRelationEntry>);
static_assert(!slam::RelationLike<int>);

static_assert(slam::MapLike<MinimalBivariateMap>);
static_assert(slam::MapLike<const MinimalBivariateMap&>);
static_assert(slam::MapLike<ExtraAliasesMap>);
static_assert(slam::MapLike<WrongDomainMap>);
static_assert(!slam::MapLike<TypedefOnlyMap>);
static_assert(!slam::MapLike<DropsDataConstnessMap>);
static_assert(!slam::MapLike<MissingComponentCountMap>);
static_assert(!slam::MapLike<MissingIndexMap>);
static_assert(!slam::MapLike<MissingConstAccessMap>);
static_assert(!slam::MapLike<WrongPositionMap>);
static_assert(!slam::MapLike<WrongValueTypeMap>);
static_assert(!slam::MapLike<WrongMutableAccessMap>);
static_assert(!slam::MapLike<int>);

static_assert(slam::MapOver<BoundMap, BoundMap::Domain>);
static_assert(slam::MapOver<const BoundMap&, const BoundMap::Domain&>);
static_assert(!slam::MapOver<WrongDomainMap, MinimalBivariateSet>);
static_assert(!slam::MapOver<BoundMap, ExternalSet>);
static_assert(!slam::MapOver<MinimalBivariateMap, BoundMap::Domain>);
static_assert(!slam::MapOver<int, ExternalSet>);

// Refinement must participate in overload ordering, not merely yield true.
template <slam::SetLike S>
std::integral_constant<int, 1> selectByConstraint();
template <slam::BivariateSetLike S>
std::integral_constant<int, 2> selectByConstraint();

static_assert(decltype(selectByConstraint<ExternalSet>())::value == 1);
static_assert(decltype(selectByConstraint<MinimalBivariateSet>())::value == 2);

// Policies
using Size = policies::CompileTimeSize<int, 5>;
using EmptySize = policies::ZeroSize<int>;
using RuntimeSize = policies::RuntimeSize<int>;
using WrongRuntimeSize = policies::RuntimeSize<double>;
using WrongCompileTimeSize = policies::CompileTimeSize<std::int64_t, 5>;
using WrongEmptySize = policies::ZeroSize<std::int64_t>;
using ScalarStride = policies::CompileTimeStride<int, 3>;
using MatrixStride = policies::MultiDimStride<int, 2>;
using Offset = policies::CompileTimeOffset<int, 4>;
using RuntimeOffset = policies::RuntimeOffset<int>;
using EmptyOffset = policies::ZeroOffset<int>;
using WrongRuntimeOffset = policies::RuntimeOffset<std::int64_t>;
using WrongCompileTimeOffset = policies::CompileTimeOffset<std::int64_t, 4>;
using WrongEmptyOffset = policies::ZeroOffset<std::int64_t>;
using NoIndirection = policies::NoIndirection<Position, Element>;
using OwningIndirection = policies::ArrayIndirection<Position, double>;

static_assert(slam::SizePolicy<Size>);
static_assert(!slam::StridePolicy<Size>);
static_assert(!slam::SizePolicy<ScalarStride>);
static_assert(!slam::SetLike<Size>);
static_assert(slam::SizePolicy<EmptySize>);
static_assert(slam::detail::SetSizePolicyFor<RuntimeSize, int>);
static_assert(slam::detail::SetSizePolicyFor<Size, int>);
static_assert(slam::detail::SetSizePolicyFor<EmptySize, int>);
static_assert(slam::SizePolicy<WrongRuntimeSize>);
static_assert(!slam::detail::SetSizePolicyFor<WrongRuntimeSize, int>);
static_assert(!slam::detail::SetSizePolicyFor<WrongCompileTimeSize, int>);
static_assert(!slam::detail::SetSizePolicyFor<WrongEmptySize, int>);
static_assert(slam::StridePolicy<ScalarStride>);
static_assert(slam::StridePolicy<MatrixStride>);
static_assert(slam::detail::OrderedSetStridePolicyFor<ScalarStride, int>);
static_assert(!slam::detail::OrderedSetStridePolicyFor<MatrixStride, int>);
static_assert(slam::detail::MapStridePolicyFor<ScalarStride, int>);
static_assert(slam::detail::MapStridePolicyFor<MatrixStride, int>);
// Maps support a stride index that converts to their (possibly wider) position type
// OrderedSet's scalar value policy must use that exact position type.
#if !defined(AXOM_NO_INT64_T)
// Use explicit types for wide/narrow pairs so they do not depend on AXOM_USE_64BIT_INDEXTYPE.
using WidePosition = std::int64_t;
using NarrowStride = policies::CompileTimeStride<std::int32_t, 3>;
static_assert(!std::is_same_v<NarrowStride::IndexType, WidePosition>,
              "the wide and narrow types must actually differ for this to test anything");
static_assert(slam::detail::MapStridePolicyFor<NarrowStride, WidePosition>);
static_assert(!slam::detail::OrderedSetStridePolicyFor<NarrowStride, WidePosition>);
#endif

static_assert(slam::OffsetPolicy<Offset>);
static_assert(slam::detail::OrderedSetOffsetPolicyFor<RuntimeOffset, int>);
static_assert(slam::detail::OrderedSetOffsetPolicyFor<Offset, int>);
static_assert(slam::detail::OrderedSetOffsetPolicyFor<EmptyOffset, int>);
static_assert(!slam::detail::OrderedSetOffsetPolicyFor<WrongRuntimeOffset, int>);
static_assert(!slam::detail::OrderedSetOffsetPolicyFor<WrongCompileTimeOffset, int>);
static_assert(!slam::detail::OrderedSetOffsetPolicyFor<WrongEmptyOffset, int>);
static_assert(!slam::SizePolicy<int>);
static_assert(!slam::StridePolicy<int>);
static_assert(!slam::OffsetPolicy<int>);
static_assert(!slam::OrderedSetIndirectionPolicyFor<TypedefOnlyIndirection, Position, double>);
static_assert(!slam::MapIndirectionPolicyFor<TypedefOnlyIndirection, Position, double>);
static_assert(slam::OrderedSetIndirectionPolicyFor<ViewIndirection, Position, double>);
// Access may add constness to the requested element type, but cannot remove it.
static_assert(slam::OrderedSetIndirectionPolicyFor<ConstViewIndirection, Position, const double>);
static_assert(!slam::OrderedSetIndirectionPolicyFor<ViewIndirection, Position, const double>);
static_assert(slam::OrderedSetIndirectionPolicyFor<ConstViewIndirection, Position, double>);
static_assert(slam::MapIndirectionPolicyFor<ViewIndirection, Position, double>);
static_assert(slam::MapIndirectionPolicyFor<ConstViewIndirection, Position, const double>);
static_assert(!slam::MapIndirectionPolicyFor<ViewIndirection, Position, const double>);
static_assert(slam::MapIndirectionPolicyFor<ConstViewIndirection, Position, double>);
static_assert(!slam::MapIndirectionPolicyFor<DropsReferentConst, Position, const double>);
static_assert(!slam::OrderedSetIndirectionPolicyFor<DropsReferentConst, Position, const double>);
static_assert(!slam::OrderedSetIndirectionPolicyFor<AccessDropsReferentConst, Position, const double>);
static_assert(!slam::AllocatingMapIndirectionPolicyFor<ViewIndirection, Position, double>);
static_assert(!slam::MapIndirectionPolicyFor<RawPointerMapIndirection, Position, double>);
static_assert(!slam::MapIndirectionPolicyFor<PrvalueMapIndirection, Position, double>);
static_assert(!slam::MapIndirectionPolicyFor<MismatchedConstPointerIndirection, Position, double>);
static_assert(slam::MapIndirectionPolicyFor<VectorMapIndirection, Position, double>);
static_assert(slam::AllocatingMapIndirectionPolicyFor<VectorMapIndirection, Position, double>);
static_assert(slam::OrderedSetIndirectionPolicyFor<NoIndirection, Position, Element>);
static_assert(!slam::MapIndirectionPolicyFor<NoIndirection, Position, Element>);
static_assert(slam::MapIndirectionPolicyFor<OwningIndirection, Position, double>);
static_assert(slam::AllocatingMapIndirectionPolicyFor<OwningIndirection, Position, double>);
static_assert(!slam::MapIndirectionPolicyFor<WrongDataIndirection, Position, double>);

// The Map indirection contract that Map family reads.
// Indirection[Const]RefType are StaticRelation's aliases
static_assert(slam::MapIndirectionPolicyFor<LeanMapIndirection, Position, double>);
// ... but the whole-buffer accessor behind Map::data_ptr() is required.
static_assert(!slam::MapIndirectionPolicyFor<PositionedAccessOnlyIndirection, Position, double>);
// Map's data_ptr() (private) is declared in terms of the policy's ResultPtr alias
static_assert(std::same_as<typename ViewIndirection::ResultPtr, double*>);
static_assert(std::same_as<typename ViewIndirection::ConstResultPtr, double*>,
              "ArrayView indirection has shallow constness");

static_assert(slam::SizePolicy<const RuntimeSize&>);
static_assert(slam::StridePolicy<const ScalarStride&>);
static_assert(slam::OffsetPolicy<const Offset&>);
static_assert(slam::SubsetPolicy<policies::NoSubset>);
static_assert(!slam::SubsetPolicy<int>);
static_assert(!slam::detail::SetSizePolicyFor<const RuntimeSize, int>);
static_assert(!slam::OrderedSetIndirectionPolicyFor<const ViewIndirection, Position, double>);
static_assert(!slam::OrderedSetIndirectionPolicyFor<ViewIndirection&, Position, double>);
static_assert(!slam::OrderedSetIndirectionPolicyFor<FinalIndirection, Position, double>);
static_assert(!slam::OrderedSetIndirectionPolicyFor<AbstractIndirection, Position, double>);
static_assert(!slam::MapIndirectionPolicyFor<const ViewIndirection, Position, double>);
static_assert(!slam::AllocatingMapIndirectionPolicyFor<LeanMapIndirection, Position, double>);

// Position types and representation properties
static_assert(slam::PositionLike<int>);
static_assert(slam::PositionLike<const std::int64_t&>);
static_assert(!slam::PositionLike<unsigned>);
static_assert(!slam::PositionLike<bool>);
static_assert(!slam::PositionLike<double>);
static_assert(!slam::PositionLike<StrongPosition>);
static_assert(slam::TriviallyCopyableRepresentation<TrivialCapture>);
static_assert(slam::TriviallyCopyableRepresentation<const TrivialCapture>);
static_assert(!slam::TriviallyCopyableRepresentation<NonTrivialCapture>);
static_assert(!slam::TriviallyCopyableRepresentation<TrivialCapture&>);

// An algorithm over the public core needs only row traversal and coordinate access.
template <slam::BivariateSetLike S>
void checkRowCoordinates(const S& set)
{
  typename S::PositionType flat = 0;
  for(typename S::FirstSetType::PositionType first = 0; first < set.getFirstSet()->size(); ++first)
  {
    const auto row = set.getElements(first);
    decltype(row.size()) count = 0;
    for(const auto second : row)
    {
      ASSERT_LT(flat, set.size());
      const auto coordinate = set.at(flat++);
      EXPECT_EQ(first, coordinate.first);
      EXPECT_EQ(second, coordinate.second);
      EXPECT_GE(second, 0);
      EXPECT_LT(second, set.getSecondSet()->size());
      ++count;
    }
    EXPECT_EQ(count, row.size());
  }
  EXPECT_EQ(flat, set.size());
}

template <slam::RelationLike R>
auto relationPositions(const R& relation)
{
  using Pair = std::pair<typename R::FromSetType::PositionType, typename R::ToSetType::PositionType>;
  std::vector<Pair> pairs;
  for(typename R::FromSetType::PositionType from = 0; from < relation.fromSet()->size(); ++from)
  {
    for(const auto to : relation[from])
    {
      pairs.emplace_back(from, to);
    }
  }
  return pairs;
}

// These tests use maps with no invalid entries within their extent.
// Writable access is constrained separately from MapLike.
template <slam::MapLike M>
  requires requires(M& map) { map.flatValue(0, 0) = 0.; }
void fillComponents(M& map)
{
  using P = typename M::PositionType;
  for(P pos = 0; pos < map.size(); ++pos)
  {
    for(P component = 0; component < map.numComp(); ++component)
    {
      map.flatValue(pos, component) = 100. * pos + component;
    }
  }
}

TEST(slam_concepts, external_models_use_only_core_operations)
{
  MinimalBivariateSet set;
  checkRowCoordinates(set);
  EXPECT_EQ(0, set.at(0).first);
  EXPECT_EQ(1, set.at(2).first);
  EXPECT_EQ(2, set.at(2).second);

  MinimalRelation relation;
  using Pair = std::pair<std::int32_t, std::int64_t>;
  const std::vector<Pair> expected {{0, 2}, {1, 0}, {1, 1}};
  EXPECT_EQ(expected, relationPositions(relation));

  MinimalBivariateMap map;
  fillComponents(map);
  const auto& const_map = map;
  EXPECT_DOUBLE_EQ(102., const_map.flatValue(1, 2));
  EXPECT_EQ(1, const_map.index(1).first);
  EXPECT_EQ(2, const_map.index(1).second);
  EXPECT_EQ(&map.values[5], &const_map.flatValue(1, 2));
}

TEST(slam_concepts, whole_set_binding_agrees_with_element_lookup)
{
  BoundMap map;
  fillComponents(map);
  EXPECT_EQ(map.size(), map.set()->size());
  for(BoundMap::PositionType pos = 0; pos < map.size(); ++pos)
  {
    EXPECT_EQ(map.index(pos).first, map.set()->at(pos).first);
    EXPECT_EQ(map.index(pos).second, map.set()->at(pos).second);
  }
  EXPECT_DOUBLE_EQ(102., map.flatValue(1, 2));
}
}  // namespace slam_concept_test
