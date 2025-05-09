// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_relation_StaticVariable.cpp
 *
 * \brief Unit tests for Slam's StaticRelation class
 *  configured with variable per-element cardinality
 */

#include <iostream>
#include <iterator>
#include <sstream>

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/slam/ModularInt.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/slam/Relation.hpp"

#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/StaticRelation.hpp"

namespace
{
namespace slam = axom::slam;
namespace policies = axom::slam::policies;

using SetPosition = slam::DefaultPositionType;
using SetElement = slam::DefaultElementType;

using RangeSetType = slam::RangeSet<SetPosition, SetElement>;
using RelationType = slam::Relation<SetPosition, SetElement>;

using IndexVec = std::vector<SetPosition>;

constexpr SetPosition FROMSET_SIZE = 7;
constexpr SetPosition TOSET_SIZE = 8;

using STLIndirection = policies::STLVectorIndirection<SetPosition, SetElement>;

using VariableCardinality = policies::VariableCardinality<SetPosition, STLIndirection>;

using StaticVariableRelationType =
  slam::StaticRelation<SetPosition, SetElement, VariableCardinality, STLIndirection, RangeSetType, RangeSetType>;

// Use a slam::ModularInt type for more interesting test data
using CTSize = policies::CompileTimeSize<int, TOSET_SIZE>;
using FixedModularInt = slam::ModularInt<CTSize>;

SetPosition elementCardinality(SetPosition fromPos) { return fromPos; }

SetPosition relationData(SetPosition fromPos, SetPosition toPos)
{
  return FixedModularInt(fromPos + toPos);
}

template <typename StrType, typename VecType>
void printVector(StrType const& msg, VecType const& vec)
{
  std::stringstream sstr;

  sstr << "\n** " << msg << "\n\t";
  sstr << "Array of size " << vec.size() << ": ";
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<SetPosition>(sstr, " "));

  SLIC_INFO(sstr.str());
}

template <typename VecType>
void generateIncrementingRelations(VecType* begins, VecType* indices)
{
  VecType& beginsVec = *begins;
  VecType& indicesVec = *indices;

  SetPosition curIdx = SetPosition();

  for(SetPosition i = 0; i < FROMSET_SIZE; ++i)
  {
    beginsVec.push_back(curIdx);
    for(SetPosition j = 0; j < elementCardinality(i); ++j)
    {
      indicesVec.push_back(relationData(i, j));
      ++curIdx;
    }
  }
  beginsVec.push_back(curIdx);
}

/**
 * \brief Traverses the relation's entities using the double subscript access API
 *
 * \note Assumes that the relation data has been set
 * using the generateIncrementingRelations() function
 */
template <typename RelationType>
void traverseRelation_doubleSubscript(RelationType& rel)
{
  SLIC_INFO("Traversing relation data using double subscript: ");
  for(SetPosition fromPos = 0; fromPos < rel.fromSet()->size(); ++fromPos)
  {
    const SetPosition fromSize = rel.size(fromPos);
    EXPECT_EQ(elementCardinality(fromPos), fromSize);

    for(int toPos = 0; toPos < fromSize; ++toPos)
    {
      SetPosition actualVal = rel[fromPos][toPos];

      EXPECT_EQ(relationData(fromPos, toPos), actualVal);
    }
  }
}

/**
 * \brief Traverses relation using separated subscript operators
 *
 * The first subscript operator gets the set of entities in the ToSet
 * that are mapped to the given element of the relation's FromSet
 *
 * \note Assumes that the relation data has been set using
 * the generateIncrementingRelations() function
 */
template <typename RelationType>
void traverseRelation_delayedSubscript(RelationType& rel)
{
  SLIC_INFO("Traversing relation data using delayed second subscript: ");
  for(SetPosition fromPos = 0; fromPos < rel.fromSet()->size(); ++fromPos)
  {
    const SetPosition fromSize = rel.size(fromPos);
    EXPECT_EQ(elementCardinality(fromPos), fromSize);

    typename RelationType::RelationSubset set = rel[fromPos];
    for(int toPos = 0; toPos < set.size(); ++toPos)
    {
      SetPosition actualVal = rel[fromPos][toPos];
      EXPECT_EQ(relationData(fromPos, toPos), actualVal);
    }
  }
}

/**
 * \brief Traverses relation using the iterator API (begin()/end() )
 *
 * \note The iterator API depends on C++11
 * \note Assumes that the relation data has been set
 * using the generateIncrementingRelations() function
 */
template <typename RelationType>
void iterateRelation_begin_end(RelationType& rel)
{
  using FromSet = typename RelationType::FromSetType;
  using FromSetIter = typename FromSet::iterator;
  using RelIter = typename RelationType::RelationIterator;

  SLIC_INFO("Traversing relation data using iterator begin()/end() functions");
  for(FromSetIter sIt = rel.fromSet()->begin(), sItEnd = rel.fromSet()->end(); sIt != sItEnd; ++sIt)
  {
    SetPosition actualSize = rel.size(*sIt);

    SetPosition fromSetEltNum = std::distance(rel.fromSet()->begin(), sIt);
    EXPECT_EQ(elementCardinality(fromSetEltNum), actualSize);

    RelIter toSetBegin = rel.begin(*sIt);
    RelIter toSetEnd = rel.end(*sIt);
    for(RelIter relIt = toSetBegin; relIt != toSetEnd; ++relIt)
    {
      SetPosition toSetEltNum = std::distance(toSetBegin, relIt);
      ASSERT_EQ(relationData(fromSetEltNum, toSetEltNum), *relIt);
    }
  }
}

/**
 * \brief Traverses relation using the iterator range API
 *
 * \note The iterator API depends on C++11
 * \note Assumes that the relation data has been set
 * using the generateIncrementingRelations() function
 */
template <typename RelationType>
void iterateRelation_range(RelationType& rel)
{
  using FromSet = typename RelationType::FromSetType;
  using FromSetIter = typename FromSet::iterator;
  using FromSetIterPair = typename FromSet::iterator_pair;

  using RelIter = typename RelationType::RelationIterator;
  using RelIterPair = typename RelationType::RelationIteratorPair;

  SLIC_INFO("Traversing relation data using iterator range() functions");
  FromSetIterPair itPair = rel.fromSet()->range();
  for(FromSetIter sIt = itPair.first; sIt != itPair.second; ++sIt)
  {
    SetPosition fromSetEltNum = std::distance(itPair.first, sIt);

    RelIterPair toSetItPair = rel.range(*sIt);
    for(RelIter relIt = toSetItPair.first; relIt != toSetItPair.second; ++relIt)
    {
      SetPosition toSetEltNum = std::distance(toSetItPair.first, relIt);
      ASSERT_EQ(relationData(fromSetEltNum, toSetEltNum), *relIt);
    }
  }
}

}  // end anonymous namespace

TEST(slam_static_variable_relation, construct_empty)
{
  SLIC_INFO("Testing empty relation.  isValid() should be false.");

  StaticVariableRelationType emptyRel;

  EXPECT_FALSE(emptyRel.isValid(true));
}

TEST(slam_static_variable_relation, construct_uninitialized)
{
  SLIC_INFO("Testing uninitialized relation.  isValid() should be false.");

  RangeSetType fromSet(FROMSET_SIZE);
  RangeSetType toSet(TOSET_SIZE);

  StaticVariableRelationType emptyRel(&fromSet, &toSet);

  EXPECT_FALSE(emptyRel.isValid(true));
}

TEST(slam_static_variable_relation, construct_relation)
{
  SLIC_INFO("Testing simple incrementing relation.  isValid() should be true.");

  RangeSetType fromSet(FROMSET_SIZE);
  RangeSetType toSet(TOSET_SIZE);

  IndexVec relOffsets;
  IndexVec relIndices;

  SLIC_INFO("Uninitialized relation data");
  printVector("begins vector", relOffsets);
  printVector("indices vector", relIndices);
  generateIncrementingRelations(&relOffsets, &relIndices);

  SLIC_INFO("Initialized relation data");
  printVector("begins vector", relOffsets);
  printVector("indices vector", relIndices);

  StaticVariableRelationType incrementingRel(&fromSet, &toSet);
  incrementingRel.bindBeginOffsets(fromSet.size(), &relOffsets);
  incrementingRel.bindIndices(relIndices.size(), &relIndices);

  EXPECT_TRUE(incrementingRel.isValid(true));

  // Test traversal of the relation data
  SCOPED_TRACE("Vector_simple_construct");
  traverseRelation_doubleSubscript(incrementingRel);
  traverseRelation_delayedSubscript(incrementingRel);
  iterateRelation_begin_end(incrementingRel);
  iterateRelation_range(incrementingRel);
}

TEST(slam_static_variable_relation, construct_builder)
{
  SLIC_INFO("Testing construction using builder interface.");

  RangeSetType fromSet(FROMSET_SIZE);
  RangeSetType toSet(TOSET_SIZE);

  IndexVec offsets;
  IndexVec relIndices;
  generateIncrementingRelations(&offsets, &relIndices);

  using RelationBuilder = StaticVariableRelationType::RelationBuilder;
  StaticVariableRelationType relation = RelationBuilder()
                                          .fromSet(&fromSet)
                                          .toSet(&toSet)
                                          .begins(RelationBuilder::BeginsSetBuilder()  //
                                                    .size(offsets.size())              //
                                                    .data(&offsets))
                                          .indices(RelationBuilder::IndicesSetBuilder()  //
                                                     .size(relIndices.size())            //
                                                     .data(&relIndices));
  EXPECT_TRUE(relation.isValid(true));

  // Test traversal of the relation data
  SCOPED_TRACE("Vector_builder");
  traverseRelation_doubleSubscript(relation);
  traverseRelation_delayedSubscript(relation);
  iterateRelation_begin_end(relation);
  iterateRelation_range(relation);
}

TEST(slam_static_variable_relation, empty_relation_out_of_bounds)
{
  StaticVariableRelationType emptyRel;

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode, so this test will only fail
  // in debug mode
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED(emptyRel[FROMSET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif
}

/// Tests for out-of-bounds access

TEST(slam_static_variable_relation, initialized_rel_out_of_bounds)
{
  SLIC_INFO("Testing simple incrementing relation.  isValid() should be true.");

  IndexVec relOffsets, relIndices;
  generateIncrementingRelations(&relOffsets, &relIndices);

  RangeSetType fromSet(FROMSET_SIZE), toSet(TOSET_SIZE);
  StaticVariableRelationType incrementingRel(&fromSet, &toSet);
  incrementingRel.bindBeginOffsets(fromSet.size(), &relOffsets);
  incrementingRel.bindIndices(relIndices.size(), &relIndices);

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode, so this test will only fail
  // in debug mode
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED(incrementingRel[FROMSET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
