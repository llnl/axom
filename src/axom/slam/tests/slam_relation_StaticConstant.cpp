// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_relation_StaticConstant.cpp
 *
 * \brief Unit tests for Slam's StaticRelation class
 * configured with constant per-element cardinality
 *
 * Exercises several different variants of the class.
 * Namely:
 * * with a runtime or compile time striding of offsets per element
 * * where the underlying indirection array uses STL vectors or arrays
 *
 */

#include <iostream>
#include <iterator>

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/slam/RangeSet.hpp"
#include "axom/slam/Relation.hpp"
#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/policies/CardinalityPolicies.hpp"
#include "axom/slam/ModularInt.hpp"

#include "axom/slam/policies/StridePolicies.hpp"
#include "axom/slam/StaticRelation.hpp"

namespace
{
namespace slam = axom::slam;
namespace policies = axom::slam::policies;

using SetPosition = slam::DefaultPositionType;
using SetElement = slam::DefaultElementType;

using RangeSetType = slam::RangeSet<SetPosition, SetElement>;

constexpr SetPosition FROMSET_SIZE = 5;
constexpr SetPosition TOSET_SIZE = 6;
constexpr SetPosition ELEM_STRIDE = 5;

using CTStride = policies::CompileTimeStride<SetPosition, ELEM_STRIDE>;
using RTStride = policies::RuntimeStride<SetPosition>;

using ConstantCardinalityCT = policies::ConstantCardinality<SetPosition, CTStride>;
using ConstantCardinalityRT = policies::ConstantCardinality<SetPosition, RTStride>;

using STLIndirection = policies::STLVectorIndirection<SetPosition, SetElement>;
using CArrayIndirection = policies::CArrayIndirection<SetPosition, SetElement>;

using IndexVec = std::vector<SetPosition>;

using StaticConstantRelationType =
  slam::StaticRelation<SetPosition, SetElement, ConstantCardinalityCT, STLIndirection, RangeSetType, RangeSetType>;

// Use a slam::ModularInt type for more interesting test data
using CTSize = policies::CompileTimeSize<int, ELEM_STRIDE>;
using FixedModularInt = slam::ModularInt<CTSize>;

template <typename StrType, typename VecType>
void printVector(StrType const& msg, VecType const& vec)
{
  std::stringstream sstr;

  sstr << "Array of size " << vec.size() << ": ";
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<SetElement>(sstr, " "));

  SLIC_INFO(msg << ": " << sstr.str());
}

SetPosition elementCardinality(SetPosition AXOM_UNUSED_PARAM(fromPos)) { return ELEM_STRIDE; }

/**
 * \brief Sets the value at relation element (i,j) to (i + j) % ELEM_SIZE using
 *  slam::ModularInt
 */
SetPosition relationData(SetPosition fromPos, SetPosition toPos)
{
  return FixedModularInt(fromPos + toPos);
}

template <typename VecType>
void generateIncrementingRelations(SetPosition stride, VecType* offsets)
{
  VecType& offsetsVec = *offsets;

  for(SetPosition i = 0; i < FROMSET_SIZE; ++i)
  {
    EXPECT_EQ(elementCardinality(i), stride);

    for(SetPosition j = 0; j < elementCardinality(i); ++j)
    {
      offsetsVec.push_back(relationData(i, j));
    }
  }
}

/**
 * \brief Traverses the relation's entities using the double subscript access
 *  API
 *
 * \note Expects cardinality and relation data for each element to match the
 *  results of the elementCardinality() and relationData() functions above,
 *  respectively.
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
 * The first subscript operator gets the set of entities in the ToSet that are
 * mapped to the given element of the relation's FromSet
 *
 * \note Expects cardinality and relation data for each element to match the
 *  results of the elementCardinality() and relationData() functions above,
 *  respectively.
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
 * \note Expects cardinality and relation data for each element to match the
 *  results of the elementCardinality() and relationData() functions above,
 *  respectively.
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
 * \note Expects cardinality and relation data for each element to match the
 *  results of the elementCardinality() and relationData() functions above,
 *  respectively.
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

TEST(slam_relation_static_constant, construct_empty)
{
  SLIC_INFO("Testing empty relation.  isValid() should be false.");

  StaticConstantRelationType emptyRel;
  EXPECT_FALSE(emptyRel.isValid());
}

TEST(slam_relation_static_constant, construct_uninitialized)
{
  SLIC_INFO("Testing uninitialized relation.  isValid() should be false.");

  RangeSetType fromSet(FROMSET_SIZE);
  RangeSetType toSet(TOSET_SIZE);

  StaticConstantRelationType emptyRel(&fromSet, &toSet);
  EXPECT_FALSE(emptyRel.isValid(true));
}

TEST(slam_relation_static_constant, construct_relation)
{
  SLIC_INFO("Testing simple incrementing relation.  isValid() should be true.");

  IndexVec relIndices;
  generateIncrementingRelations(ELEM_STRIDE, &relIndices);
  printVector("Relation indices vector ", relIndices);

  RangeSetType fromSet(FROMSET_SIZE);
  RangeSetType toSet(TOSET_SIZE);

  StaticConstantRelationType incrementingRel(&fromSet, &toSet);
  // intialize the begins set
  incrementingRel.bindBeginOffsets(fromSet.size(), ELEM_STRIDE);
  // initialize the relation indices set
  incrementingRel.bindIndices(relIndices.size(), &relIndices);

  EXPECT_TRUE(incrementingRel.isValid(true));

  // Test traversal of the relation data
  SCOPED_TRACE("CTStride_vector_simple_construct");
  traverseRelation_doubleSubscript(incrementingRel);
  traverseRelation_delayedSubscript(incrementingRel);
  iterateRelation_begin_end(incrementingRel);
  iterateRelation_range(incrementingRel);
}

TEST(slam_relation_static_constant, construct_builder)
{
  SLIC_INFO("Checking if we can instantiate a concrete StaticRelation (constant).");

  using RelationBuilder = StaticConstantRelationType::RelationBuilder;

  RangeSetType fromSet(FROMSET_SIZE);
  RangeSetType toSet(TOSET_SIZE);

  IndexVec offsets;
  generateIncrementingRelations(ELEM_STRIDE, &offsets);

  StaticConstantRelationType relation =
    RelationBuilder()
      .fromSet(&fromSet)
      .toSet(&toSet)
      .begins(RelationBuilder::BeginsSetBuilder().stride(ELEM_STRIDE))
      .indices(RelationBuilder::IndicesSetBuilder().size(offsets.size()).data(&offsets));

  EXPECT_TRUE(relation.isValid(true));

  // Test traversal of the relation data
  SCOPED_TRACE("CTStride_vector_builder_construct");
  traverseRelation_doubleSubscript(relation);
  traverseRelation_delayedSubscript(relation);
  iterateRelation_begin_end(relation);
  iterateRelation_range(relation);
}

/// Tests for out-of-bounds access

TEST(slam_relation_static_constant, out_of_bounds_empty)
{
  SLIC_INFO("Testing access on empty relation " << "-- code is expected to assert and die.");

  RangeSetType fromSet(FROMSET_SIZE);
  RangeSetType toSet(TOSET_SIZE);
  StaticConstantRelationType emptyRel(&fromSet, &toSet);

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode,
  // so this test will only fail in debug mode

  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED(emptyRel[FROMSET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure checks in release mode.");
#endif
}

TEST(slam_relation_static_constant, out_of_bounds_initialized)
{
  SLIC_INFO("Testing out of bounds access on initialized relation.  "
            << "Code is expected to assert and die.");

  RangeSetType fromSet(FROMSET_SIZE);
  RangeSetType toSet(TOSET_SIZE);

  IndexVec relIndices;
  generateIncrementingRelations(ELEM_STRIDE, &relIndices);

  StaticConstantRelationType incrementingRel(&fromSet, &toSet);
  // initialize the begins set
  incrementingRel.bindBeginOffsets(fromSet.size(), ELEM_STRIDE);
  // intialize the relation indices set
  incrementingRel.bindIndices(relIndices.size(), &relIndices);

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode,
  // so this test will only fail in debug mode

  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED(incrementingRel[FROMSET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure checks in release mode.");
#endif
}

/// Tests that vary the relation's policies (STL vs. array indirection and
// compile time vs. runtime offset striding)

TEST(slam_relation_static_constant, runtime_stride_STLIndirection)
{
  SLIC_INFO("Tests for Static Relation with runtime stride and STL Indirection");

  using StaticConstantRelation_RT_STL =
    slam::StaticRelation<SetPosition, SetElement, ConstantCardinalityRT, STLIndirection, RangeSetType, RangeSetType>;

  RangeSetType fromSet(FROMSET_SIZE);
  RangeSetType toSet(TOSET_SIZE);

  IndexVec relIndices;
  generateIncrementingRelations(ELEM_STRIDE, &relIndices);

  // --  Construct empty and uninitialized relation
  StaticConstantRelation_RT_STL emptyRel;
  EXPECT_FALSE(emptyRel.isValid(true));

  StaticConstantRelation_RT_STL uninitRel(&fromSet, &toSet);
  EXPECT_FALSE(uninitRel.isValid(true));

  // -- Simple relation construction
  StaticConstantRelation_RT_STL relation(&fromSet, &toSet);
  relation.bindBeginOffsets(fromSet.size(), ELEM_STRIDE);
  relation.bindIndices(relIndices.size(), &relIndices);

  EXPECT_TRUE(relation.isValid(true));

  // Test traversal of the relation data
  SCOPED_TRACE("RTStride_vector_simple_construct");
  traverseRelation_doubleSubscript(relation);
  traverseRelation_delayedSubscript(relation);
  iterateRelation_begin_end(relation);
  iterateRelation_range(relation);

  // -- Construction using set and relation builder objects
  using RelationBuilder = StaticConstantRelation_RT_STL::RelationBuilder;
  StaticConstantRelation_RT_STL builderRel = RelationBuilder()
                                               .fromSet(&fromSet)
                                               .toSet(&toSet)
                                               .begins(RelationBuilder::BeginsSetBuilder()  //
                                                         .stride(ELEM_STRIDE))
                                               .indices(RelationBuilder::IndicesSetBuilder()  //
                                                          .size(relIndices.size())            //
                                                          .data(&relIndices));
  EXPECT_TRUE(builderRel.isValid(true));

  // Test traversal of the relation data
  SCOPED_TRACE("RTStride_vector_builder construct");
  traverseRelation_doubleSubscript(builderRel);
  traverseRelation_delayedSubscript(builderRel);
  iterateRelation_begin_end(builderRel);
  iterateRelation_range(builderRel);
}

TEST(slam_relation_static_constant, runtime_stride_CArrayIndirection)
{
  SLIC_INFO("Tests for Static Relation " << " with runtime stride and C-array Indirection");

  using StaticConstantRelation_RT_Array =
    slam::StaticRelation<SetPosition, SetElement, ConstantCardinalityRT, CArrayIndirection, RangeSetType, RangeSetType>;

  RangeSetType fromSet(FROMSET_SIZE);
  RangeSetType toSet(TOSET_SIZE);

  IndexVec relIndices;
  generateIncrementingRelations(ELEM_STRIDE, &relIndices);
  SetPosition* data = &relIndices[0];  // Get a pointer to the data

  // --  Construct empty and uninitialized relation
  StaticConstantRelation_RT_Array emptyRel;
  EXPECT_FALSE(emptyRel.isValid());

  StaticConstantRelation_RT_Array uninitRel(&fromSet, &toSet);
  EXPECT_FALSE(uninitRel.isValid());

  // -- Simple relation construction
  StaticConstantRelation_RT_Array relation(&fromSet, &toSet);
  relation.bindBeginOffsets(fromSet.size(), ELEM_STRIDE);
  relation.bindIndices(relIndices.size(), data);

  EXPECT_TRUE(relation.isValid(true));

  // Test traversal of the relation data
  SCOPED_TRACE("RTStride_array_simple_construct");
  traverseRelation_doubleSubscript(relation);
  traverseRelation_delayedSubscript(relation);
  iterateRelation_begin_end(relation);
  iterateRelation_range(relation);

  // -- Construction using set and relation builder objects
  using RelationBuilder = StaticConstantRelation_RT_Array::RelationBuilder;
  StaticConstantRelation_RT_Array builderRel =
    RelationBuilder()
      .fromSet(&fromSet)
      .toSet(&toSet)
      .begins(RelationBuilder::BeginsSetBuilder().stride(ELEM_STRIDE))
      .indices(RelationBuilder::IndicesSetBuilder().size(relIndices.size()).data(data));
  EXPECT_TRUE(builderRel.isValid(true));

  // Test traversal of the relation data
  SCOPED_TRACE("RTStride_array_builder construct");
  traverseRelation_doubleSubscript(builderRel);
  traverseRelation_delayedSubscript(builderRel);
  iterateRelation_begin_end(builderRel);
  iterateRelation_range(builderRel);
}

TEST(slam_relation_static_constant, compileTime_stride_CArrayIndirection)
{
  SLIC_INFO("Tests for Static Relation with " << " runtime stride and C-array Indirection");

  using StaticConstantRelation_CT_Array =
    slam::StaticRelation<SetPosition, SetElement, ConstantCardinalityCT, CArrayIndirection, RangeSetType, RangeSetType>;

  RangeSetType fromSet(FROMSET_SIZE);
  RangeSetType toSet(TOSET_SIZE);

  IndexVec relIndices;
  generateIncrementingRelations(ELEM_STRIDE, &relIndices);
  SetPosition* data = &relIndices[0];  // Get a pointer to the data

  // --  Construct empty and uninitialized relation
  StaticConstantRelation_CT_Array emptyRel;
  EXPECT_FALSE(emptyRel.isValid());

  StaticConstantRelation_CT_Array uninitRel(&fromSet, &toSet);
  EXPECT_FALSE(uninitRel.isValid());

  // -- Simple relation construction
  StaticConstantRelation_CT_Array relation(&fromSet, &toSet);
  relation.bindIndices(relIndices.size(), data);

  // Note: Since the striding is a compile time constant,
  //       we don't need to set the striding
  EXPECT_TRUE(relation.isValid(true));

  //        .. but we can, if we'd like to
  relation.bindBeginOffsets(fromSet.size(), ELEM_STRIDE);
  EXPECT_TRUE(relation.isValid(true));

  // Test traversal of the relation data
  SCOPED_TRACE("RTStride_array_simple_construct");
  traverseRelation_doubleSubscript(relation);
  traverseRelation_delayedSubscript(relation);
  iterateRelation_begin_end(relation);
  iterateRelation_range(relation);

  // -- Construction using set and relation builder objects
  using RelationBuilder = StaticConstantRelation_CT_Array::RelationBuilder;
  StaticConstantRelation_CT_Array builderRel = RelationBuilder()
                                                 .fromSet(&fromSet)
                                                 .toSet(&toSet)
                                                 .begins(RelationBuilder::BeginsSetBuilder()  //
                                                           .stride(ELEM_STRIDE))
                                                 .indices(RelationBuilder::IndicesSetBuilder()  //
                                                            .size(relIndices.size())            //
                                                            .data(data));
  EXPECT_TRUE(builderRel.isValid(true));

  // Test traversal of the relation data
  SCOPED_TRACE("RTStride_array_builder construct");
  traverseRelation_doubleSubscript(builderRel);
  traverseRelation_delayedSubscript(builderRel);
  iterateRelation_begin_end(builderRel);
  iterateRelation_range(builderRel);

  // Similarly here, it is ok to omit the striding in the setup
  // due to the compile time striding policy
  StaticConstantRelation_CT_Array builderRel_implicitStride =
    RelationBuilder()
      .fromSet(&fromSet)
      .toSet(&toSet)
      //  .begins( RelationBuilder::BeginsSetBuilder()
      //          .stride(ELEM_STRIDE))
      .indices(RelationBuilder::IndicesSetBuilder()  //
                 .size(relIndices.size())            //
                 .data(data));
  EXPECT_TRUE(builderRel_implicitStride.isValid(true));

  // Test traversal of the relation data
  SCOPED_TRACE("RTStride_array_builder construct");
  traverseRelation_doubleSubscript(builderRel_implicitStride);
  traverseRelation_delayedSubscript(builderRel_implicitStride);
  iterateRelation_begin_end(builderRel_implicitStride);
  iterateRelation_range(builderRel_implicitStride);
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
