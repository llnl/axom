// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_set_BivariateSet.cpp
 *
 * This file tests BivariateSet, ProductSet and RelationSet within Slam.
 * It uses a templated test fixture to test many different types of
 * bivariate sets
 */

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"

#include <type_traits>
#include <sstream>
#include <iostream>

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
  slam::StaticRelation<P, E, policies::VariableCardinality<E>, policies::STLVectorIndirection<P, E>, FromSet, ToSet>;

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

  using Vec = std::vector<PositionType>;
  using RelationType = ::RelType<PositionType, ElementType, FirstSetType, SecondSetType>;

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
      relationBegins.push_back(relationIndices.size());

      for(int j = 0; j < m_set2->size(); ++j)
      {
        auto inner = m_set2->at(j);
        if(modCheck(outer, inner))
        {
          relationIndices.push_back(j);
        }
      }
    }
    relationBegins.push_back(relationIndices.size());

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

  Vec setIndices1;
  Vec setIndices2;

  Vec relationBegins;
  Vec relationIndices;
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
      EXPECT_EQ(bsetElem.firstIndex(), derivedSet->flatToFirstIndex(flatIndex));
      EXPECT_EQ(bsetElem.secondIndex(), derivedSet->flatToSecondIndex(flatIndex));

      sstr << flatIndex << ": (" << firstSet->at(bsetElem.firstIndex()) << ","
           << secondSet->at(bsetElem.secondIndex()) << "), ";

      flatIndex++;
    }

    SLIC_INFO("{ " << sstr.str() << " }");
  }
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

//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  // create & initialize test logger. finalized when exiting main scope
  axom::slic::SimpleLogger logger;

  int result = RUN_ALL_TESTS();

  return result;
}
