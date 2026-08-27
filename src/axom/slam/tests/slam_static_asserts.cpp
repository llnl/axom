// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file slam_static_asserts.cpp
 *
 * \brief Compile-time conformance harness for slam's constexpr arithmetic
 *
 * The tests in this file exercise the compile time set/value asserts
 * for sizes, offsets, strides, modular wraparound, and policy composition, so the checks cost
 * nothing at run time and run on every backend on every build. 
 * 
 * Many of the checks in this file are static_asserts. If the file compiles, its invariants hold.
 */

#include "gtest/gtest.h"

#include "axom/slam/ModularInt.hpp"
#include "axom/slam/Traits.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/slam/ProductSet.hpp"
#include "axom/slam/StaticRelation.hpp"
#include "axom/slam/BivariateMap.hpp"
#include "axom/slam/Map.hpp"
#include "axom/slam/Aliases.hpp"
#include "axom/slam/MapBuilders.hpp"
#include "axom/slam/RelationBuilders.hpp"
#include "axom/slam/SetBuilders.hpp"
#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/policies/StridePolicies.hpp"
#include "axom/slam/policies/OffsetPolicies.hpp"
#include "axom/slam/policies/ValuePolicies.hpp"
#include "axom/slam/policies/CardinalityPolicies.hpp"
#include "axom/slam/policies/IndirectionPolicies.hpp"

#include <iterator>
#include <optional>
#include <utility>

namespace
{
namespace slam = axom::slam;
namespace policies = axom::slam::policies;

//------------------------------------------------------------------------------
// Value policies: compile-time-known values report their value, validity,
// and (for size) emptiness without any runtime state.
//------------------------------------------------------------------------------

// --- Size ---
using Size0 = policies::CompileTimeSize<int, 0>;
using Size5 = policies::CompileTimeSize<int, 5>;

static_assert(Size5().size() == 5, "CompileTimeSize reports its size");
static_assert(Size0().size() == 0, "CompileTimeSize of zero reports zero");
static_assert(!Size5().empty(), "non-zero CompileTimeSize is not empty");
static_assert(Size0().empty(), "zero CompileTimeSize is empty");
static_assert(Size5().isValid(false), "non-negative size is valid");
static_assert(Size5().DEFAULT_VALUE == 5, "DEFAULT_VALUE matches the NTTP");

// --- Offset ---
using Off0 = policies::CompileTimeOffset<int, 0>;
using Off3 = policies::CompileTimeOffset<int, 3>;

static_assert(Off3().offset() == 3, "CompileTimeOffset reports its offset");
static_assert(Off0().offset() == 0, "ZeroOffset-equivalent reports zero");
static_assert(policies::ZeroOffset<int>().offset() == 0, "ZeroOffset is zero");
static_assert(Off3().isValid(false), "all offsets are valid");

// --- Stride ---
using Stride1 = policies::CompileTimeStride<int, 1>;
using Stride4 = policies::CompileTimeStride<int, 4>;

static_assert(Stride4().stride() == 4, "CompileTimeStride reports its stride");
static_assert(Stride4().shape() == 4, "CompileTimeStride shape equals stride");
static_assert(policies::StrideOne<int>().stride() == 1, "StrideOne is one");
static_assert(Stride4().isValid(false), "non-zero stride is valid");
static_assert(Stride1::IS_COMPILE_TIME, "CompileTimeStride is compile-time");
static_assert(Stride1::NumDims == 1, "scalar stride is one-dimensional");

//------------------------------------------------------------------------------
// The unified core: RuntimeValue / CompileTimeValue behave through the tags.
//------------------------------------------------------------------------------
static_assert(policies::CompileTimeValue<7, policies::SizeTag<int>>().value() == 7,
              "CompileTimeValue forwards its NTTP");
static_assert(policies::SizeTag<int>::isValidValue(0), "size 0 is valid");
static_assert(!policies::SizeTag<int>::isValidValue(-1), "negative size is invalid");
static_assert(!policies::StrideTag<int>::isValidValue(0), "stride 0 is invalid");
static_assert(policies::StrideTag<int>::isValidValue(-2), "negative stride is valid");
static_assert(policies::OffsetTag<int>::isValidValue(-5), "any offset is valid");
static_assert(policies::StrideTag<int>::defaultValue() == 1, "default stride is one");
static_assert(policies::SizeTag<int>::defaultValue() == 0, "default size is zero");

//------------------------------------------------------------------------------
// Offset / stride composition: the flat index of element i in a strided,
// offset range is  offset + i * stride.  Verify the policy values compose.
//------------------------------------------------------------------------------
constexpr int flatIndex(int offset, int stride, int i) { return offset + i * stride; }

static_assert(flatIndex(Off3().offset(), Stride4().stride(), 0) == 3,
              "element 0 lands at the offset");
static_assert(flatIndex(Off3().offset(), Stride4().stride(), 2) == 11,
              "offset 3 + 2*stride 4 == 11");
static_assert(flatIndex(policies::ZeroOffset<int>().offset(), policies::StrideOne<int>().stride(), 9) ==
                9,
              "unit stride, zero offset is the identity index");

//------------------------------------------------------------------------------
// ModularInt: cyclic arithmetic is fully constexpr-evaluable.
//------------------------------------------------------------------------------
using Mod5 = slam::ModularInt<policies::CompileTimeSize<int, 5>>;

static_assert(int(Mod5(0)) == 0, "0 mod 5 == 0");
static_assert(int(Mod5(5)) == 0, "5 mod 5 == 0 (wraparound)");
static_assert(int(Mod5(7)) == 2, "7 mod 5 == 2");
static_assert(int(Mod5(-1)) == 4, "-1 mod 5 normalizes to 4");
static_assert(int(Mod5(-5)) == 0, "-5 mod 5 == 0");
static_assert(int(Mod5(13)) == 3, "13 mod 5 == 3");

// arithmetic operators
static_assert(int(Mod5(3) + 4) == 2, "(3+4) mod 5 == 2");
static_assert(int(Mod5(4) + 1) == 0, "(4+1) mod 5 wraps to 0");
static_assert(int(4 + Mod5(3)) == 2, "int + ModularInt commutes");
static_assert(int(Mod5(1) - 3) == 3, "(1-3) mod 5 == 3");
static_assert(int(Mod5(2) * 4) == 3, "(2*4) mod 5 == 3");

// pre/post increment in a constexpr lambda
constexpr int incTwice(int start)
{
  Mod5 m(start);
  ++m;
  ++m;
  return int(m);
}
static_assert(incTwice(4) == 1, "++ twice from 4 mod 5 == 1");

// equality across the modulus
static_assert(Mod5(2) == Mod5(7), "2 and 7 are equal mod 5");
static_assert(Mod5(2) != Mod5(3), "2 and 3 differ mod 5");

//------------------------------------------------------------------------------
// Check trait predicates for policies and sets
//------------------------------------------------------------------------------

// Value policies satisfy their policy concept, and are distinguished from each
// other and from sets (a set has size() but no value()).
static_assert(slam::is_value_policy_v<Size5>, "a size policy is a value policy");
static_assert(slam::is_size_policy_v<Size5>, "CompileTimeSize is a size policy");
static_assert(slam::is_stride_policy_v<Stride4>, "CompileTimeStride is a stride policy");
static_assert(slam::is_offset_policy_v<Off3>, "CompileTimeOffset is an offset policy");
static_assert(!slam::is_stride_policy_v<Size5>, "a size policy is not a stride policy");
static_assert(!slam::is_size_policy_v<Stride4>, "a stride policy is not a size policy");
static_assert(!slam::is_set_like_v<Size5>, "a policy is not a set");

// RangeSet models the (ordered) set concept and nothing else.
static_assert(slam::is_set_like_v<slam::RangeSet<>>, "RangeSet is set-like");
static_assert(slam::is_ordered_set_like_v<slam::RangeSet<>>, "RangeSet is an ordered set");
static_assert(!slam::is_relation_like_v<slam::RangeSet<>>, "RangeSet is not a relation");
static_assert(!slam::is_map_like_v<slam::RangeSet<>>, "RangeSet is not a map");
static_assert(!slam::is_bivariate_set_like_v<slam::RangeSet<>>, "RangeSet is not bivariate");

// Index-space concepts: raw integrals are positions
static_assert(slam::is_position_like_v<int>, "int is a position");
static_assert(slam::is_position_like_v<axom::slam::DefaultPositionType>,
              "slam's default position type is a position");
static_assert(!slam::is_position_like_v<double>, "double is not a position");

// An element handle must be trivially copyable so it survives capture-by-value into a device kernel
struct TrivialHandle
{
  int id;
};
static_assert(slam::is_handle_like_v<TrivialHandle>, "a trivially-copyable struct is handle-like");

struct UserCopyHandle
{
  UserCopyHandle(const UserCopyHandle&) { }
};
static_assert(!slam::is_handle_like_v<UserCopyHandle>,
              "a user-declared copy ctor breaks the handle contract");

//------------------------------------------------------------------------------
// Device-capture contract: a slam container survives capture-by-value into a
// kernel only if it is trivially copyable, which requires the concrete interface.
//------------------------------------------------------------------------------
using SetPos = slam::DefaultPositionType;
using SetElem = slam::DefaultElementType;
using ViewInd = policies::ArrayViewIndirection<SetPos, SetElem>;
using ViewPosInd = policies::ArrayViewIndirection<SetPos, SetPos>;
using ConcreteRangeSet = slam::RangeSet<SetPos, SetElem>::ConcreteSet;
using ConcreteProductSet = slam::ProductSet<ConcreteRangeSet, ConcreteRangeSet>::ConcreteSet;
using ViewMap = slam::Map<SetElem, ConcreteRangeSet, ViewInd>;
using ViewBivariateMap = slam::BivariateMap<SetElem, ConcreteProductSet, ViewInd>;
using ViewVariableRelation =  //
  slam::StaticRelation<SetPos,
                       SetElem,
                       policies::VariableCardinality<SetPos, ViewPosInd>,
                       ViewInd,
                       ConcreteRangeSet,
                       ConcreteRangeSet>;
using ViewConstantRelation =
  slam::StaticRelation<SetPos,
                       SetElem,
                       policies::ConstantCardinality<SetPos, policies::CompileTimeStride<SetPos, 3>>,
                       ViewInd,
                       ConcreteRangeSet,
                       ConcreteRangeSet>;

static_assert(std::is_trivially_copyable_v<ViewInd>,
              "ArrayViewIndirection holds a view by value and must be trivially copyable");
static_assert(std::is_trivially_copyable_v<ConcreteRangeSet>,
              "a concrete-interface RangeSet is trivially copyable");
static_assert(std::is_trivially_copyable_v<ConcreteProductSet>,
              "a concrete-interface ProductSet is trivially copyable");
static_assert(std::is_trivially_copyable_v<ViewMap>, "a view-backed Map is trivially copyable");
static_assert(std::is_trivially_copyable_v<ViewBivariateMap>,
              "a view-backed BivariateMap is trivially copyable");
static_assert(std::is_trivially_copyable_v<ViewVariableRelation>,
              "a view-backed variable relation is trivially copyable");
static_assert(std::is_trivially_copyable_v<ViewConstantRelation>,
              "a view-backed constant relation is trivially copyable");

// The default set uses the virtual interface so it is not trivially copyable
static_assert(!std::is_trivially_copyable_v<slam::RangeSet<>>,
              "the virtual-interface set is not device-capturable by design");

//------------------------------------------------------------------------------
// OrderedSetIterator models C++17 RandomAccessIterator
// std::iterator_traits is complete for both the mutable and const iterator
// and reports the random-access category. The single const-qualified
// accessor path also makes the iterator const-dereferenceable.
//------------------------------------------------------------------------------
using RangeIter = slam::RangeSet<>::iterator;
using RangeConstIter = slam::RangeSet<>::const_iterator;

static_assert(
  std::is_same_v<std::iterator_traits<RangeIter>::iterator_category, std::random_access_iterator_tag>,
  "the mutable ordered-set iterator is random-access");
static_assert(std::is_same_v<std::iterator_traits<RangeConstIter>::iterator_category,
                             std::random_access_iterator_tag>,
              "the const ordered-set iterator is random-access");
static_assert(
  std::is_same_v<std::iterator_traits<RangeIter>::difference_type, slam::RangeSet<>::PositionType>,
  "iterator difference_type is the set's position type");
static_assert(std::is_integral_v<std::iterator_traits<RangeIter>::difference_type> &&
                std::is_signed_v<std::iterator_traits<RangeIter>::difference_type>,
              "a random-access difference_type is a signed integral");
static_assert(!std::is_void_v<std::iterator_traits<RangeConstIter>::value_type>,
              "iterator_traits exposes a value_type");

// maybe_const_t adds const exactly when requested.
static_assert(std::is_same_v<slam::maybe_const_t<true, int>, const int>,
              "maybe_const_t<true, T> is const T");
static_assert(std::is_same_v<slam::maybe_const_t<false, int>, int>, "maybe_const_t<false, T> is T");

//------------------------------------------------------------------------------
// Checks alias layer (Aliases.hpp): each alias instantiates with Slam's
// canonical Array/ArrayView indirection and satisfies the matching predicate.
//------------------------------------------------------------------------------
using AliasArraySet = slam::ArraySet<>;
using AliasCustomArraySet = slam::ArraySet<SetPos, double>;
using AliasArrayViewSet = slam::ArrayViewSet<>;

using AliasVarRelation = slam::VariableRelation<slam::RangeSet<>, slam::RangeSet<>>;
using AliasVarRelationView = slam::VariableRelationView<slam::RangeSet<>, slam::RangeSet<>>;
using AliasConstRelation = slam::ConstantRelation<slam::RangeSet<>, slam::RangeSet<>, 3>;
using AliasConstRelationView = slam::ConstantRelationView<slam::RangeSet<>, slam::RangeSet<>, 3>;
using AliasRuntimeConstRelation = slam::RuntimeConstantRelation<slam::RangeSet<>, slam::RangeSet<>>;
using AliasRuntimeConstRelationView =
  slam::RuntimeConstantRelationView<slam::RangeSet<>, slam::RangeSet<>>;

using MapArrayCT = slam::Map<double,
                             slam::RangeSet<>,
                             policies::ArrayIndirection<SetPos, double>,
                             policies::CompileTimeStride<SetPos, 3>>;
using MapArrayViewCT = slam::Map<double,
                                 slam::RangeSet<>,
                                 policies::ArrayViewIndirection<SetPos, double>,
                                 policies::CompileTimeStride<SetPos, 3>>;
using MapArrayViewRT = slam::Map<double,
                                 slam::RangeSet<>,
                                 policies::ArrayViewIndirection<SetPos, double>,
                                 policies::RuntimeStride<SetPos>>;
using BMapArrayViewCT =
  slam::BivariateMap<double,
                     ConcreteProductSet,
                     policies::ArrayViewIndirection<typename ConcreteProductSet::PositionType, double>,
                     policies::CompileTimeStride<typename ConcreteProductSet::PositionType, 2>>;

// Force full instantiation (also exercises axom::Array's element requirements).
static_assert(sizeof(AliasArraySet) > 0, "ArraySet instantiates");
static_assert(sizeof(AliasCustomArraySet) > 0, "custom-position ArraySet instantiates");
static_assert(sizeof(AliasArrayViewSet) > 0, "ArrayViewSet instantiates");
static_assert(sizeof(AliasVarRelation) > 0, "VariableRelation instantiates");
static_assert(sizeof(AliasVarRelationView) > 0, "VariableRelationView instantiates");
static_assert(sizeof(AliasConstRelation) > 0, "ConstantRelation instantiates");
static_assert(sizeof(AliasConstRelationView) > 0, "ConstantRelationView instantiates");
static_assert(sizeof(AliasRuntimeConstRelation) > 0, "RuntimeConstantRelation instantiates");
static_assert(sizeof(AliasRuntimeConstRelationView) > 0, "RuntimeConstantRelationView instantiates");
static_assert(sizeof(MapArrayCT) > 0, "Array-backed Map instantiates");
static_assert(sizeof(MapArrayViewCT) > 0, "ArrayView Map instantiates");
static_assert(sizeof(MapArrayViewRT) > 0, "runtime-stride ArrayView Map instantiates");
static_assert(sizeof(BMapArrayViewCT) > 0, "ArrayView BivariateMap instantiates");

// Map/BivariateMap default to an axom::Array indirection (managing their own buffer).
static_assert(std::is_same_v<typename slam::Map<double>::OrderedMap, axom::Array<double>>,
              "Map's default indirection uses an axom::Array");
static_assert(
  std::is_same_v<typename slam::BivariateMap<double>::MapType::OrderedMap, axom::Array<double>>,
  "BivariateMap's default indirection uses an axom::Array");

// Check that the aliases match the expected concepts
static_assert(slam::is_set_like_v<AliasArraySet>, "ArraySet is set-like");
static_assert(std::is_same_v<AliasCustomArraySet, slam::ArrayIndirectionSet<SetPos, double>>,
              "ArraySet preserves the normal <Position, Element> set template order");
static_assert(slam::is_ordered_set_like_v<AliasArraySet>, "ArraySet is an ordered set");
static_assert(slam::is_ordered_set_like_v<AliasArrayViewSet>, "ArrayViewSet is an ordered set");
static_assert(slam::is_relation_like_v<AliasVarRelation>, "VariableRelation is relation-like");
static_assert(slam::is_relation_like_v<AliasVarRelationView>,
              "VariableRelationView is relation-like");
static_assert(slam::is_relation_like_v<AliasConstRelation>, "ConstantRelation is relation-like");
static_assert(slam::is_relation_like_v<AliasConstRelationView>,
              "ConstantRelationView is relation-like");
static_assert(slam::is_relation_like_v<AliasRuntimeConstRelation>,
              "RuntimeConstantRelation is relation-like");
static_assert(slam::is_relation_like_v<AliasRuntimeConstRelationView>,
              "RuntimeConstantRelationView is relation-like");
static_assert(slam::is_bivariate_set_like_v<ConcreteProductSet>,
              "a ProductSet is bivariate-set-like");
static_assert(slam::is_map_like_v<MapArrayCT>, "Array-backed Map is map-like");
static_assert(slam::is_map_like_v<MapArrayViewCT>, "ArrayView Map is map-like");
static_assert(slam::is_map_like_v<MapArrayViewRT>, "runtime-stride ArrayView Map is map-like");
static_assert(slam::is_map_over_v<MapArrayCT, slam::RangeSet<>>,
              "Array-backed Map is a map over its set");
static_assert(slam::is_map_over_v<MapArrayViewCT, slam::RangeSet<>>,
              "ArrayView Map is a map over its set");
static_assert(slam::is_map_over_v<MapArrayViewRT, slam::RangeSet<>>,
              "runtime-stride ArrayView Map is a map over its set");
static_assert(slam::is_map_like_v<BMapArrayViewCT>, "ArrayView BivariateMap is map-like");

// Check that the aliases match the Array/ArrayView make_* helpers.
using RSPos = slam::RangeSet<>::PositionType;
static_assert(
  std::is_same_v<AliasVarRelation,
                 decltype(slam::make_variable_relation(std::declval<slam::RangeSet<>*>(),
                                                       std::declval<slam::RangeSet<>*>(),
                                                       std::declval<axom::Array<RSPos>&>(),
                                                       std::declval<axom::Array<RSPos>&>()))>,
  "VariableRelation matches make_variable_relation's axom::Array overload");
static_assert(
  std::is_same_v<AliasVarRelationView,
                 decltype(slam::make_variable_relation(std::declval<slam::RangeSet<>*>(),
                                                       std::declval<slam::RangeSet<>*>(),
                                                       std::declval<axom::ArrayView<RSPos>>(),
                                                       std::declval<axom::ArrayView<RSPos>>()))>,
  "VariableRelationView matches make_variable_relation's ArrayView overload");
static_assert(
  std::is_same_v<AliasConstRelation,
                 decltype(slam::make_constant_relation_ct<3>(std::declval<slam::RangeSet<>*>(),
                                                             std::declval<slam::RangeSet<>*>(),
                                                             std::declval<axom::Array<RSPos>&>()))>,
  "ConstantRelation matches make_constant_relation_ct<N>'s axom::Array overload");
static_assert(
  std::is_same_v<AliasConstRelationView,
                 decltype(slam::make_constant_relation_ct<3>(std::declval<slam::RangeSet<>*>(),
                                                             std::declval<slam::RangeSet<>*>(),
                                                             std::declval<axom::ArrayView<RSPos>>()))>,
  "ConstantRelationView matches make_constant_relation_ct<N>'s ArrayView overload");
static_assert(
  std::is_same_v<AliasRuntimeConstRelation,
                 decltype(slam::make_constant_relation(std::declval<slam::RangeSet<>*>(),
                                                       std::declval<slam::RangeSet<>*>(),
                                                       std::declval<RSPos>(),
                                                       std::declval<axom::Array<RSPos>&>()))>,
  "RuntimeConstantRelation matches make_constant_relation's axom::Array overload");
static_assert(
  std::is_same_v<AliasRuntimeConstRelationView,
                 decltype(slam::make_constant_relation(std::declval<slam::RangeSet<>*>(),
                                                       std::declval<slam::RangeSet<>*>(),
                                                       std::declval<RSPos>(),
                                                       std::declval<axom::ArrayView<RSPos>>()))>,
  "RuntimeConstantRelationView matches make_constant_relation's ArrayView overload");
static_assert(
  std::is_same_v<AliasArrayViewSet,
                 decltype(slam::make_indirection_set(std::declval<axom::ArrayView<RSPos>>()))>,
  "ArrayViewSet matches make_indirection_set's ArrayView overload");
static_assert(std::is_same_v<MapArrayViewCT,
                             decltype(slam::make_map_ct<3>(std::declval<const slam::RangeSet<>*>(),
                                                           std::declval<axom::ArrayView<double>>()))>,
              "make_map_ct<N>'s ArrayView overload yields the ArrayView Map");
static_assert(
  std::is_same_v<MapArrayViewRT,
                 decltype(slam::make_map(std::declval<const slam::RangeSet<>*>(),
                                         std::declval<RSPos>(),
                                         std::declval<axom::ArrayView<double>>()))>,
  "make_map's runtime-stride ArrayView overload yields the runtime-stride ArrayView Map");

}  // anonymous namespace

//------------------------------------------------------------------------------
TEST(slam_static_asserts, compile_time_value_holds)
{
  // Checks that a constexpr value survives to run time.

  constexpr int wrapped = int(Mod5(7));
  EXPECT_EQ(wrapped, 2);

  constexpr int idx = flatIndex(Off3().offset(), Stride4().stride(), 2);
  EXPECT_EQ(idx, 11);
}
