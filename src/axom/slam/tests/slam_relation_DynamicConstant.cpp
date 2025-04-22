// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_relation_DynamicConstant.cpp
 *
 * \brief Unit tests for Slam's DynamicConstantRelation class
 */

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
#include "axom/slam/DynamicConstantRelation.hpp"

namespace
{
namespace slam = axom::slam;
namespace policies = axom::slam::policies;

using PositionType = slam::DefaultPositionType;
using ElementType = slam::DefaultElementType;

using RangeSetType = slam::RangeSet<PositionType, ElementType>;

constexpr PositionType FROMSET_SIZE = 5;
constexpr PositionType TOSET_SIZE = 6;
constexpr PositionType ELEM_STRIDE = 6;

using CTStride = policies::CompileTimeStride<PositionType, ELEM_STRIDE>;
using RTStride = policies::RuntimeStride<PositionType>;

using ConstantCardinalityCT = policies::ConstantCardinality<PositionType, CTStride>;
using ConstantCardinalityRT = policies::ConstantCardinality<PositionType, RTStride>;

using STLIndirection = policies::STLVectorIndirection<PositionType, ElementType>;

using IndexVec = std::vector<PositionType>;
using RelationType = slam::DynamicConstantRelation<PositionType, ElementType, ConstantCardinalityCT>;

}  //end anonymous namespace

/** \brief Utility function to set the value at rel[i][j] to i */
template <typename DynamicSetType>
RelationType generateRelation(DynamicSetType& fromSet, DynamicSetType& toSet)
{
  RelationType rel(&fromSet, &toSet);

  const int fromSize = fromSet.size();
  const int toSize = toSet.size();

  for(int i = 0; i < fromSize; ++i)
  {
    for(int j = 0; j < toSize; ++j)
    {
      rel.modify(i, j, i);
    }
  }

  return rel;
}

TEST(slam_relation_dynamic_constant, construct_empty)
{
  SLIC_INFO("Testing empty relation.  isValid() should be false.");

  RelationType emptyRel;
  EXPECT_FALSE(emptyRel.isValid());
}

TEST(slam_relation_dynamic_constant, assignment)
{
  SLIC_INFO("Testing assignment of values to relation");

  slam::DynamicSet<PositionType, ElementType> fromSet(FROMSET_SIZE);
  slam::DynamicSet<PositionType, ElementType> toSet(TOSET_SIZE);

  RelationType rel = generateRelation(fromSet, toSet);

  EXPECT_TRUE(rel.isValid(true));

  const int fromSize = rel.size();
  EXPECT_EQ(FROMSET_SIZE, fromSize);

  // Check that accessed relation values are as expected
  for(int i = 0; i < fromSize; ++i)
  {
    EXPECT_EQ(rel[i].size(), rel.size(i));

    // Test double indirection
    {
      const int relSize = rel.size(i);
      EXPECT_EQ(ELEM_STRIDE, relSize);
      for(int j = 0; j < relSize; j++)
      {
        EXPECT_EQ(i, rel[i][j]);
      }
    }

    // Test RelationSubset
    {
      RelationType::RelationSubset set = rel[i];
      const int relSize = set.size();
      for(int j = 0; j < relSize; ++j)
      {
        EXPECT_EQ(i, set[j]);
      }

      EXPECT_EQ(set, rel.at(i));
    }
  }
}

TEST(slam_relation_dynamic_constant, iterators)
{
  SLIC_INFO("Testing iterator interface");

  // Add tests for relation iterators
  slam::DynamicSet<PositionType, ElementType> fromSet(FROMSET_SIZE);
  slam::DynamicSet<PositionType, ElementType> toSet(TOSET_SIZE);

  {
    RelationType rel = generateRelation(fromSet, toSet);
    const int fromSize = rel.size();

    for(int i = 0; i < fromSize; ++i)
    {
      ElementType val = i;
      PositionType idx = 0;
      for(RelationType::RelationIterator it = rel.begin(i), itEnd = rel.end(i); it != itEnd;
          ++it, ++idx)
      {
        EXPECT_EQ(val, *it);
        EXPECT_EQ(idx, it.index());
      }

      idx = 0;
      for(RelationType::RelationIteratorPair itPair = rel.range(i); itPair.first != itPair.second;
          ++itPair.first, ++idx)
      {
        EXPECT_EQ(val, *itPair.first);
        EXPECT_EQ(idx, itPair.first.index());
      }
    }
  }
}

TEST(slam_relation_dynamic_constant, const_iterators)
{
  SLIC_INFO("Testing iterator interface");

  // Add tests for relation iterators
  slam::DynamicSet<PositionType, ElementType> fromSet(FROMSET_SIZE);
  slam::DynamicSet<PositionType, ElementType> toSet(TOSET_SIZE);

  {
    const RelationType rel = generateRelation(fromSet, toSet);
    const int fromSize = rel.size();

    for(int i = 0; i < fromSize; ++i)
    {
      const ElementType val = i;
      PositionType idx = 0;
      for(RelationType::RelationConstIterator it = rel.begin(i), itEnd = rel.end(i); it != itEnd;
          ++it, ++idx)
      {
        EXPECT_EQ(val, *it);
        EXPECT_EQ(idx, it.index());
      }

      idx = 0;
      for(RelationType::RelationConstIteratorPair itPair = rel.range(i);
          itPair.first != itPair.second;
          ++itPair.first, ++idx)
      {
        EXPECT_EQ(val, *itPair.first);
        EXPECT_EQ(idx, itPair.first.index());
      }
    }
  }
}

TEST(slam_relation_dynamic_constant, remove)
{
  SLIC_INFO("Testing ability to remove relations");

  slam::DynamicSet<PositionType, ElementType> fromSet(FROMSET_SIZE);
  slam::DynamicSet<PositionType, ElementType> toSet(TOSET_SIZE);

  RelationType rel = generateRelation(fromSet, toSet);
  EXPECT_EQ(rel.size(), rel.numberOfValidEntries());

  const int fromSize = rel.size();
  int removeCount = 0;
  for(int i = 1; i < fromSize; i += 2)
  {
    rel.remove(i);
    ++removeCount;
  }

  // Check counts
  EXPECT_TRUE(rel.isValid());
  EXPECT_EQ(rel.size() - removeCount, rel.numberOfValidEntries());

  // Test individual relations
  for(int i = 0; i < fromSize; i += 2)
  {
    EXPECT_TRUE(rel.isValidEntry(i));
  }

  for(int i = 1; i < fromSize; i += 2)
  {
    EXPECT_FALSE(rel.isValidEntry(i));
  }
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  int result = RUN_ALL_TESTS();

  return result;
}
