// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_set_BivariateSet.cpp
 *
 * This file tests BivariateSet, ProductSet and RelationSet within Slam.
 * It uses a templated test fixture to test many different types of bivariate sets
 */

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"

#include <cstdint>
#include <limits>
#include <type_traits>
#include <sstream>
#include <iostream>
#include <vector>

namespace slam = axom::slam;
namespace policies = axom::slam::policies;

namespace
{
constexpr int SET_SIZE_1 = 5;
constexpr int SET_SIZE_2 = 16;

constexpr int SET_OFFSET_1 = 3;
constexpr int SET_OFFSET_2 = 2;

// Template aliases to simplify specifying some sets and relations
template <typename P, typename E, typename FromSet, typename ToSet>
using RelType =
  slam::StaticRelation<P,
                       E,
                       policies::VariableCardinality<P, policies::ArrayIndirection<P, P>>,
                       policies::ArrayIndirection<P, E>,
                       FromSet,
                       ToSet>;

template <typename SetType>
using PositionSetType =
  slam::PositionSet<typename SetType::PositionType, typename SetType::ElementType>;

template <typename SetType>
using RangeSetType = slam::RangeSet<typename SetType::PositionType, typename SetType::ElementType>;

template <typename SetType>
using IndirSetType =
  slam::VectorIndirectionSet<typename SetType::PositionType, typename SetType::ElementType>;

// Predicate to check if b%a is zero
// Note: Also returns true when a is zero to avoid division by zero
bool modCheck(int a, int b) { return (a == 0) || (b % a == 0); }

struct ZoneTag;
struct NodeTag;

template <typename EntityTag, typename IdType>
struct TypedHandle
{
  IdType id {};
};

}  // end anonymous namespace

// Test fixture for testing slam's BivariateSet class
//
// Since we test BivariateSet parameterized by different possible underlying Set
// types, we initialize all of the possibilities in this test fixture,
// and initialize the concrete instances based on the actual input types.
template <typename BSet>
class BivariateSetTester : public ::testing::Test
{
public:
  using BivariateSetType = BSet;
  using FirstSetType = typename BSet::FirstSetType;
  using SecondSetType = typename BSet::SecondSetType;

  using PositionType = typename BSet::PositionType;
  using ElementType = typename BSet::ElementType;
  using FirstElementType = typename FirstSetType::ElementType;
  using SecondElementType = typename SecondSetType::ElementType;
  using SecondPositionType = typename SecondSetType::PositionType;

  using FirstElementVec = std::vector<FirstElementType>;
  using SecondElementVec = std::vector<SecondElementType>;
  using PositionArray = axom::Array<PositionType>;
  using RelationElementArray = axom::Array<SecondPositionType>;
  using RelationType = ::RelType<PositionType, SecondPositionType, FirstSetType, SecondSetType>;

  using PSet1 = ::PositionSetType<FirstSetType>;
  using PSet2 = ::PositionSetType<SecondSetType>;

  using RSet1 = ::RangeSetType<FirstSetType>;
  using RSet2 = ::RangeSetType<SecondSetType>;

  using ISet1 = ::IndirSetType<FirstSetType>;
  using ISet2 = ::IndirSetType<SecondSetType>;

  BivariateSetTester()
    : m_set1(nullptr)
    , m_set2(nullptr)
    , m_pset1(nullptr)
    , m_pset2(nullptr)
    , m_rset1(nullptr)
    , m_rset2(nullptr)
    , m_iset1(nullptr)
    , m_iset2(nullptr)
  { }

  virtual void SetUp()
  {
    constructSets();

    m_set1 = getFirstSet<FirstSetType>();
    EXPECT_TRUE(m_set1->isValid(true));

    m_set2 = getSecondSet<SecondSetType>();
    EXPECT_TRUE(m_set2->isValid(true));

    constructRelation();
  }

  virtual void TearDown()
  {
    m_set1 = nullptr;
    m_set2 = nullptr;

    deleteSet(m_pset1);
    deleteSet(m_pset2);

    deleteSet(m_rset1);
    deleteSet(m_rset2);

    deleteSet(m_iset1);
    deleteSet(m_iset2);
  }

private:
  // Construct all the position, range and indirection sets that we might need
  void constructSets()
  {
    // Construct the position sets
    {
      m_pset1 = new PSet1(SET_SIZE_1);
      EXPECT_TRUE(m_pset1->isValid(true));

      m_pset2 = new PSet2(SET_SIZE_2);
      EXPECT_TRUE(m_pset2->isValid(true));
    }

    // Construct the range sets
    {
      m_rset1 = new RSet1(SET_OFFSET_1, SET_OFFSET_1 + SET_SIZE_1);
      EXPECT_TRUE(m_rset1->isValid(true));

      m_rset2 = new RSet2(SET_OFFSET_2, SET_OFFSET_2 + SET_SIZE_2);
      EXPECT_TRUE(m_rset2->isValid(true));
    }

    // Construct the first indirection sets elements are multiples of 2
    {
      setIndices1.resize(SET_SIZE_1);
      m_iset1 = new ISet1(SET_SIZE_1);
      m_iset1->ptr() = &setIndices1;
      for(auto idx : m_iset1->positions())
      {
        (*m_iset1)[idx] = 2 * idx;
      }

      EXPECT_TRUE(m_iset1->isValid(true));
    }

    // Construct the second indirection sets elements are multiples of 3
    {
      setIndices2.resize(SET_SIZE_2);
      m_iset2 = new ISet2(SET_SIZE_2);
      m_iset2->ptr() = &setIndices2;
      for(auto idx : m_iset2->positions())
      {
        (*m_iset2)[idx] = 3 * idx;
      }
    }
    EXPECT_TRUE(m_iset2->isValid(true));
  }

  template <typename SetType>
  void deleteSet(SetType*& s)
  {
    delete s;
    s = nullptr;
  }

  // Construct the relation for the RelationSet
  void constructRelation()
  {
    // Generate the mod relation on the two sets:
    // Add entry (outer,inner) to relation when (inner % outer == 0)
    relationBegins.clear();
    relationIndices.clear();

    for(int i = 0; i < m_set1->size(); ++i)
    {
      auto outer = m_set1->at(i);
      relationBegins.push_back(static_cast<PositionType>(relationIndices.size()));

      for(int j = 0; j < m_set2->size(); ++j)
      {
        auto inner = m_set2->at(j);
        if(modCheck(outer, inner))
        {
          relationIndices.push_back(j);
        }
      }
    }
    relationBegins.push_back(static_cast<PositionType>(relationIndices.size()));

    // Construct the relation using this data
    using RelationBuilder = typename RelationType::RelationBuilder;
    modRelation =
      RelationBuilder()
        .fromSet(m_set1)
        .toSet(m_set2)
        .begins(
          typename RelationBuilder::BeginsSetBuilder().size(relationBegins.size()).data(&relationBegins))
        .indices(typename RelationBuilder::IndicesSetBuilder()
                   .size(relationIndices.size())
                   .data(&relationIndices));

    EXPECT_TRUE(modRelation.isValid(true));
  }

private:
  // Template magic to access the appropriate first set based on
  // the FirstSetType. This is necessary since FirstSetType can
  // be either the abstract slam::Set<> or one of its derived types.
  // In the former case, we return PositionSet
  template <typename S>
  typename std::enable_if<std::is_abstract<S>::value, S*>::type getFirstSet()
  {
    return m_pset1;
  }

  template <typename S>
  typename std::enable_if<std::is_same<S, RSet1>::value, S*>::type getFirstSet()
  {
    return m_rset1;
  }

  template <typename S>
  typename std::enable_if<std::is_same<S, ISet1>::value, S*>::type getFirstSet()
  {
    return m_iset1;
  }

  template <typename S>
  typename std::enable_if<std::is_same<S, PSet1>::value, S*>::type getFirstSet()
  {
    return m_pset1;
  }

  // Template magic to access the appropriate second set based on
  // the SecondSetType. When SecondSetType is slam::Set, we return PositionSet
  template <typename S>
  typename std::enable_if<std::is_abstract<S>::value, S*>::type getSecondSet()
  {
    return m_pset2;
  }

  template <typename S>
  typename std::enable_if<std::is_same<S, RSet1>::value, S*>::type getSecondSet()
  {
    return m_rset2;
  }

  template <typename S>
  typename std::enable_if<std::is_same<S, ISet1>::value, S*>::type getSecondSet()
  {
    return m_iset2;
  }

  template <typename S>
  typename std::enable_if<std::is_same<S, PSet1>::value, S*>::type getSecondSet()
  {
    return m_pset2;
  }

protected:
  FirstSetType* m_set1;
  SecondSetType* m_set2;

  PSet1* m_pset1;
  PSet2* m_pset2;

  RSet1* m_rset1;
  RSet2* m_rset2;

  ISet1* m_iset1;
  ISet2* m_iset2;

  FirstElementVec setIndices1;
  SecondElementVec setIndices2;

  PositionArray relationBegins;
  RelationElementArray relationIndices;
  RelationType modRelation;
};

// Tests several types of BivariateSet that differ in their underlying set types
using MyTypes =
  ::testing::Types<slam::BivariateSet<>,
                   slam::BivariateSet<slam::PositionSet<>, slam::PositionSet<>>,
                   slam::BivariateSet<slam::PositionSet<>, slam::RangeSet<>>,
                   slam::BivariateSet<slam::RangeSet<>, slam::PositionSet<>>,
                   slam::BivariateSet<slam::RangeSet<>, slam::RangeSet<>>,
                   slam::BivariateSet<slam::Set<>, slam::RangeSet<>>,
                   slam::BivariateSet<slam::PositionSet<>, slam::Set<>>,
                   slam::BivariateSet<slam::VectorIndirectionSet<>, slam::VectorIndirectionSet<>>,
                   slam::BivariateSet<slam::VectorIndirectionSet<>>,
                   slam::BivariateSet<slam::VectorIndirectionSet<>, slam::RangeSet<>>,
                   slam::BivariateSet<slam::RangeSet<>, slam::VectorIndirectionSet<>>>;
TYPED_TEST_SUITE(BivariateSetTester, MyTypes);

//-----------------------------------------------------------------------------

TYPED_TEST(BivariateSetTester, smoke)
{
  using BSet = typename TestFixture::BivariateSetType;
  using S1 = typename TestFixture::FirstSetType;
  using S2 = typename TestFixture::SecondSetType;

  // Test a default-constructed BivariateSet
  {
    BSet* bset = new slam::ProductSet<S1, S2>();

    // We expect this object to be valid when the types of both S1 and S2 are
    // abstract slam::Sets, but not if either is specialized.
    bool isSet1Abstract = std::is_abstract<S1>::value;
    bool isSet2Abstract = std::is_abstract<S2>::value;

    EXPECT_EQ(isSet1Abstract && isSet2Abstract, bset->isValid(true));
    delete bset;
  }

  // Test a non-default-constructed BivariateSet
  {
    BSet* bset = new slam::ProductSet<S1, S2>(this->m_set1, this->m_set2);
    EXPECT_TRUE(bset->isValid(true));
    delete bset;
  }
}

//-----------------------------------------------------------------------------

// Tests BivariateSet::firstSetSize(), secondSetSize()
template <typename S1, typename S2>
void bSetSizesTest(slam::BivariateSet<S1, S2>* bset)
{
  const auto firstSetSize = bset->firstSetSize();
  const auto secondSetSize = bset->secondSetSize();

  EXPECT_EQ(SET_SIZE_1, firstSetSize);
  EXPECT_EQ(SET_SIZE_2, secondSetSize);

  EXPECT_EQ(firstSetSize, bset->getFirstSet()->size());
  EXPECT_EQ(secondSetSize, bset->getSecondSet()->size());
}

TYPED_TEST(BivariateSetTester, sizes)
{
  using S1 = typename TestFixture::FirstSetType;
  using S2 = typename TestFixture::SecondSetType;

  // Test the size functions for a ProductSet
  {
    using PSet = slam::ProductSet<S1, S2>;

    PSet pset = PSet(this->m_set1, this->m_set2);
    EXPECT_TRUE(pset.isValid(true));

    bSetSizesTest<S1, S2>(&pset);
  }

  // Test the size functions for a RelationSet
  {
    using RType = typename TestFixture::RelationType;
    using RSet = slam::RelationSet<RType, S1, S2>;

    RSet rset = RSet(&this->modRelation);
    EXPECT_TRUE(rset.isValid(true));

    bSetSizesTest<S1, S2>(&rset);
  }
}

//-----------------------------------------------------------------------------

TEST(slam_set_bivariate_set, concrete_interface_abstract_endpoint_sizes)
{
  using AbstractSetType = slam::Set<>;
  using ProductSetType = typename slam::ProductSet<AbstractSetType, AbstractSetType>::ConcreteSet;

  slam::RangeSet<> firstSet(SET_SIZE_1);
  slam::RangeSet<> secondSet(SET_SIZE_2);
  ProductSetType productSet(&firstSet, &secondSet);

  EXPECT_EQ(productSet.firstSetSize(), SET_SIZE_1);
  EXPECT_EQ(productSet.secondSetSize(), SET_SIZE_2);
  EXPECT_EQ(productSet.size(), SET_SIZE_1 * SET_SIZE_2);
}

//-----------------------------------------------------------------------------

// Tests BivariateSet::getFirstSet(), getSecondSet(), getElements(idx)
// and the ability to operate and iterate on the resulting sets
template <typename S1, typename S2, typename DerivedSetType>
void bSetTraverseTest(slam::BivariateSet<S1, S2>* bset, bool shouldCheckMod)
{
  const auto firstSetSize = bset->firstSetSize();

  const auto* firstSet = bset->getFirstSet();
  const auto* secondSet = bset->getSecondSet();

  int flatIdx = 0;
  for(int idx = 0; idx < firstSetSize; ++idx)
  {
    auto outer = firstSet->at(idx);

    std::stringstream sstr;
    sstr << idx << ": " << outer << " -> { ";

    auto elems = bset->getElements(idx);
    int sparseIdx = 0;
    for(auto innerIdx : elems)
    {
      auto inner = secondSet->at(innerIdx);

      if(shouldCheckMod)
      {
        EXPECT_TRUE(modCheck(outer, inner));
      }

      EXPECT_EQ(sparseIdx, bset->findElementIndex(idx, innerIdx));
      EXPECT_EQ(flatIdx, bset->findElementFlatIndex(idx, innerIdx));

      {
        auto sparse_opt = bset->findElementIndexOptional(idx, innerIdx);
        ASSERT_TRUE(sparse_opt.has_value());
        EXPECT_EQ(sparseIdx, *sparse_opt);
      }

      {
        auto flat_opt = bset->findElementFlatIndexOptional(idx, innerIdx);
        ASSERT_TRUE(flat_opt.has_value());
        EXPECT_EQ(flatIdx, *flat_opt);
      }

      EXPECT_EQ(idx, bset->flatToFirstIndex(flatIdx));
      EXPECT_EQ(innerIdx, bset->flatToSecondIndex(flatIdx));

      sstr << "(" << outer << "," << inner << ") ";
      sparseIdx++;
      flatIdx++;
    }
    SLIC_INFO(sstr.str() << "}");
  }

  SLIC_INFO(
    "Flat traversal through virtual bivariate set:\n       "
    "----------------------");

  {
    std::stringstream sstr;

    // Iterate through the bivariate set as a list of indexes (i, j) in {S1, S2}
    int flatIndex = 0;
    for(auto bsetElem = bset->begin(); bsetElem != bset->end(); ++bsetElem)
    {
      EXPECT_EQ(flatIndex, bsetElem.flatIndex());
      EXPECT_EQ(flatIndex, bset->findElementFlatIndex(bsetElem.firstIndex(), bsetElem.secondIndex()));

      {
        auto flat_opt =
          bset->findElementFlatIndexOptional(bsetElem.firstIndex(), bsetElem.secondIndex());
        ASSERT_TRUE(flat_opt.has_value());
        EXPECT_EQ(flatIndex, *flat_opt);
      }

      EXPECT_EQ(bsetElem.firstIndex(), bset->flatToFirstIndex(flatIndex));
      EXPECT_EQ(bsetElem.secondIndex(), bset->flatToSecondIndex(flatIndex));

      sstr << flatIndex << ": (" << firstSet->at(bsetElem.firstIndex()) << ","
           << secondSet->at(bsetElem.secondIndex()) << "), ";

      flatIndex++;
    }

    SLIC_INFO("{ " << sstr.str() << " }");
  }

  DerivedSetType* derivedSet = static_cast<DerivedSetType*>(bset);

  SLIC_INFO(
    "Flat traversal through derived bivariate set:\n       "
    "----------------------");

  {
    std::stringstream sstr;

    // Iterate through the bivariate set as a list of indexes (i, j) in {S1, S2}
    int flatIndex = 0;
    for(auto bsetElem = derivedSet->begin(); bsetElem != derivedSet->end(); ++bsetElem)
    {
      EXPECT_EQ(flatIndex, bsetElem.flatIndex());
      EXPECT_EQ(flatIndex,
                derivedSet->findElementFlatIndex(bsetElem.firstIndex(), bsetElem.secondIndex()));

      {
        auto flat_opt =
          derivedSet->findElementFlatIndexOptional(bsetElem.firstIndex(), bsetElem.secondIndex());
        ASSERT_TRUE(flat_opt.has_value());
        EXPECT_EQ(flatIndex, *flat_opt);
      }

      EXPECT_EQ(bsetElem.firstIndex(), derivedSet->flatToFirstIndex(flatIndex));
      EXPECT_EQ(bsetElem.secondIndex(), derivedSet->flatToSecondIndex(flatIndex));

      sstr << flatIndex << ": (" << firstSet->at(bsetElem.firstIndex()) << ","
           << secondSet->at(bsetElem.secondIndex()) << "), ";

      flatIndex++;
    }

    SLIC_INFO("{ " << sstr.str() << " }");
  }
}

template <typename SetType>
std::pair<typename SetType::PositionType, typename SetType::PositionType> find_missing_mod_pair(
  const SetType& bset)
{
  using PositionType = typename SetType::PositionType;
  const auto* firstSet = bset.getFirstSet();
  const auto* secondSet = bset.getSecondSet();

  for(PositionType i = 0; i < bset.firstSetSize(); ++i)
  {
    const auto outer = firstSet->at(i);
    if(outer == 0)
    {
      continue;
    }

    for(PositionType j = 0; j < bset.secondSetSize(); ++j)
    {
      const auto inner = secondSet->at(j);
      if(!modCheck(outer, inner))
      {
        return {i, j};
      }
    }
  }

  return {SetType::INVALID_POS, SetType::INVALID_POS};
}

TYPED_TEST(BivariateSetTester, traverse)
{
  using S1 = typename TestFixture::FirstSetType;
  using S2 = typename TestFixture::SecondSetType;

  // Test traversal functions for a ProductSet
  SLIC_INFO("Traversing product set" << "\n       ----------------------");
  {
    using PSet = slam::ProductSet<S1, S2>;

    PSet pset = PSet(this->m_set1, this->m_set2);
    EXPECT_TRUE(pset.isValid(true));

    bool checkMod = false;
    bSetTraverseTest<S1, S2, PSet>(&pset, checkMod);
  }

  // Test traversal functions for a RelationSet
  std::cout << std::endl;
  SLIC_INFO("Traversing relation set" << "\n       -----------------------");
  {
    using RType = typename TestFixture::RelationType;
    using RSet = slam::RelationSet<RType, S1, S2>;

    RSet rset = RSet(&this->modRelation);
    EXPECT_TRUE(rset.isValid(true));

    bool checkMod = true;
    bSetTraverseTest<S1, S2, RSet>(&rset, checkMod);
  }
}

TEST(slam_bivariate_set, product_set_empty_second_set_has_no_first_flat_index)
{
  using SetType = slam::PositionSet<>;
  using ProductSetType = slam::ProductSet<SetType, SetType>;

  SetType firstSet(3);
  SetType secondSet(0);
  ProductSetType productSet(&firstSet, &secondSet);

  ASSERT_TRUE(productSet.isValid(true));
  EXPECT_EQ(productSet.firstSetSize(), 3);
  EXPECT_EQ(productSet.secondSetSize(), 0);
  EXPECT_EQ(productSet.size(), 0);
  EXPECT_EQ(productSet.size(0), 0);
  EXPECT_EQ(productSet.getElements(0).size(), 0);

  EXPECT_EQ(productSet.findElementFlatIndex(0), ProductSetType::INVALID_POS);
  EXPECT_FALSE(productSet.findElementFlatIndexOptional(0).has_value());

  const slam::BivariateSet<SetType, SetType>* baseSet = &productSet;
  EXPECT_EQ(baseSet->findElementFlatIndex(0), ProductSetType::INVALID_POS);
  EXPECT_FALSE(baseSet->findElementFlatIndexOptional(0).has_value());
}

TEST(slam_bivariate_set, product_set_elements_are_coordinate_pairs)
{
  using PositionType = std::int32_t;
  using FirstSetType = slam::RangeSet<PositionType, double>;
  using SecondSetType = slam::RangeSet<PositionType, float>;
  using ProductSetType = slam::ProductSet<FirstSetType, SecondSetType>;

  static_assert(
    std::is_same_v<typename ProductSetType::ElementType, std::pair<PositionType, PositionType>>);

  FirstSetType firstSet(2);
  SecondSetType secondSet(3);
  ProductSetType productSet(&firstSet, &secondSet);

  ASSERT_TRUE(productSet.isValid(true));
  EXPECT_EQ(productSet.size(), 6);
  EXPECT_EQ(productSet.at(5), std::make_pair(PositionType {1}, PositionType {2}));
  EXPECT_EQ(productSet.flatToSecondIndex(5), PositionType {2});
}

TEST(slam_bivariate_set, product_set_preserves_heterogeneous_coordinate_types)
{
  using FirstPosition = std::int32_t;
  using SecondPosition = std::int64_t;
  using FirstSet = slam::RangeSet<FirstPosition, double>;
  using SecondSet = slam::RangeSet<SecondPosition, float>;
  using Product = slam::ProductSet<FirstSet, SecondSet>;

  static_assert(std::is_same_v<typename Product::FirstPositionType, FirstPosition>);
  static_assert(std::is_same_v<typename Product::SecondPositionType, SecondPosition>);
  static_assert(std::is_same_v<typename Product::PositionType, SecondPosition>);

  FirstSet firstSet(2);
  SecondSet secondSet(3);
  Product product(&firstSet, &secondSet);

  const auto coordinate = product.at(5);
  EXPECT_EQ(coordinate.first, FirstPosition {1});
  EXPECT_EQ(coordinate.second, SecondPosition {2});
  EXPECT_EQ(product.findElementFlatIndex(FirstPosition {1}, SecondPosition {2}), 5);
}

TEST(slam_bivariate_set, heterogeneous_coordinates_preserve_indices_above_int32)
{
  using FirstPosition = std::int32_t;
  using SecondPosition = std::int64_t;
  using FirstSet = slam::RangeSet<FirstPosition, double>;
  using SecondSet = slam::RangeSet<SecondPosition, float>;
  using Product = typename slam::ProductSet<FirstSet, SecondSet>::ConcreteSet;

  constexpr SecondPosition largeSecondPosition =
    static_cast<SecondPosition>(std::numeric_limits<FirstPosition>::max()) + 7;
  FirstSet firstSet(1);
  SecondSet secondSet(largeSecondPosition + 1);

  // RangeSet and ProductSet represent this coordinate without allocating one
  // entry per position.
  Product product(&firstSet, &secondSet);
  EXPECT_EQ(product.flatToSecondIndex(largeSecondPosition), largeSecondPosition);
  EXPECT_EQ(product.at(largeSecondPosition), std::make_pair(FirstPosition {0}, largeSecondPosition));

  std::vector<SecondPosition> begins {0, 1};
  std::vector<SecondPosition> indices {largeSecondPosition};
  auto relation = slam::make_variable_relation(firstSet, secondSet, begins, indices);
  using RelationSet = typename slam::RelationSet<decltype(relation)>::ConcreteSet;
  RelationSet connectivity(&relation);

  EXPECT_EQ(connectivity.flatToSecondIndex(0), largeSecondPosition);
  EXPECT_EQ(connectivity.at(0), std::make_pair(FirstPosition {0}, largeSecondPosition));
}

TEST(slam_bivariate_set, relation_set_projects_distinct_typed_handles)
{
  using ZonePosition = std::int32_t;
  using NodePosition = std::int64_t;
  using FlatPosition = std::common_type_t<ZonePosition, NodePosition>;
  using ZoneHandle = TypedHandle<ZoneTag, ZonePosition>;
  using NodeHandle = TypedHandle<NodeTag, NodePosition>;
  using ZoneSet = slam::VectorIndirectionSet<ZonePosition, ZoneHandle>;
  using NodeSet = slam::VectorIndirectionSet<NodePosition, NodeHandle>;

  static_assert(!std::is_same_v<ZoneHandle, NodeHandle>);

  std::vector<ZoneHandle> zoneHandles {{100}, {200}};
  std::vector<NodeHandle> nodeHandles {{10}, {20}, {30}, {40}};
  ZoneSet zones = slam::make_indirection_set<ZonePosition>(zoneHandles);
  NodeSet nodes = slam::make_indirection_set<NodePosition>(nodeHandles);

  std::vector<FlatPosition> begins {0, 3, 6};
  std::vector<NodePosition> nodePositions {0, 1, 2, 1, 2, 3};
  auto relation = slam::make_variable_relation(zones, nodes, begins, nodePositions);
  using Relation = decltype(relation);

  static_assert(slam::RelationLike<Relation>);
  static_assert(std::is_same_v<typename Relation::FromPositionType, ZonePosition>);
  static_assert(std::is_same_v<typename Relation::ToPositionType, NodePosition>);
  static_assert(std::is_same_v<typename Relation::FlatPositionType, FlatPosition>);
  static_assert(std::is_same_v<typename Relation::FromSetType::ElementType, ZoneHandle>);
  static_assert(std::is_same_v<typename Relation::ToSetType::ElementType, NodeHandle>);

  using ConnectivitySet = typename slam::RelationSet<Relation>::ConcreteSet;
  ConnectivitySet connectivity(&relation);
  static_assert(slam::BivariateSetLike<ConnectivitySet>);
  static_assert(
    std::is_same_v<typename ConnectivitySet::ElementType, std::pair<ZonePosition, NodePosition>>);

  using WeightMap = slam::BivariateMap<double, ConnectivitySet>;
  static_assert(slam::BivariateMapLike<WeightMap>);
  static_assert(slam::MapOver<WeightMap, ConnectivitySet>);
  WeightMap weights(connectivity, 0.0);

  NodePosition processedNodeIds = 0;
  for(ZonePosition zonePos = 0; zonePos < zones.size(); ++zonePos)
  {
    const auto row = relation[zonePos];
    for(FlatPosition local = 0; local < row.size(); ++local)
    {
      const NodePosition nodePos = row[local];
      const FlatPosition flat = connectivity.findElementFlatIndex(zonePos, nodePos);
      weights.flatValue(flat) = static_cast<double>(zones[zonePos].id + nodes[nodePos].id);
      processedNodeIds += nodes[nodePos].id;
    }
  }

  ASSERT_TRUE(relation.isValid(true));
  ASSERT_TRUE(connectivity.isValid(true));
  EXPECT_EQ(processedNodeIds, NodePosition {150});
  const auto* weight = weights.findValue(ZonePosition {1}, NodePosition {2});
  ASSERT_NE(weight, nullptr);
  EXPECT_DOUBLE_EQ(*weight, 230.0);

  const auto coordinate = connectivity.at(FlatPosition {4});
  EXPECT_EQ(zones[coordinate.first].id, ZonePosition {200});
  EXPECT_EQ(nodes[coordinate.second].id, NodePosition {30});

  const auto rowMap = weights(ZonePosition {1});
  EXPECT_EQ(rowMap.index(1), std::make_pair(ZonePosition {1}, NodePosition {2}));
}

TYPED_TEST(BivariateSetTester, optional_find_returns_empty_for_missing_relation_entries)
{
  using S1 = typename TestFixture::FirstSetType;
  using S2 = typename TestFixture::SecondSetType;
  using RType = typename TestFixture::RelationType;
  using RSet = slam::RelationSet<RType, S1, S2>;

  RSet rset = RSet(&this->modRelation);
  ASSERT_TRUE(rset.isValid(true));

  const auto [i, j] = find_missing_mod_pair(rset);
  ASSERT_NE(i, RSet::INVALID_POS);
  ASSERT_NE(j, RSet::INVALID_POS);

  auto sparse_opt = rset.findElementIndexOptional(i, j);
  EXPECT_FALSE(sparse_opt.has_value());

  auto flat_opt = rset.findElementFlatIndexOptional(i, j);
  EXPECT_FALSE(flat_opt.has_value());
}

//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  // create & initialize test logger. finalized when exiting main scope
  axom::slic::SimpleLogger logger;

  int result = RUN_ALL_TESTS();

  return result;
}
