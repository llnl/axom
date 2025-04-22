// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_AccessingRelationDataInMap.cpp
 *
 * \brief This file tests Sets, Relations and Maps working together.
 *
 * \note There are no maps.
 */

#include <iostream>
#include <iterator>
#include <sstream>

#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/slam/ModularInt.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/slam/Relation.hpp"

#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/policies/CardinalityPolicies.hpp"

#include "axom/slam/StaticRelation.hpp"

#include "axom/slam/Map.hpp"

namespace
{
namespace slam = axom::slam;
namespace policies = axom::slam::policies;

using SetPosition = slam::DefaultPositionType;
using SetElement = slam::DefaultElementType;

using RangeSetType = slam::RangeSet<SetPosition, SetElement>;
using RelationType = slam::Relation<SetPosition, SetElement>;

using IndexVec = std::vector<SetPosition>;

constexpr SetPosition FROMSET_SIZE = 10;
constexpr SetPosition TOSET_SIZE = 8;

using STLIndirection = policies::STLVectorIndirection<SetPosition, SetElement>;

using VariableCardinality = policies::VariableCardinality<SetPosition, STLIndirection>;

using StaticVariableRelationType =
  slam::StaticRelation<SetPosition, SetElement, VariableCardinality, STLIndirection, RangeSetType, RangeSetType>;

// Use a slam::ModularInt type for more interesting test data
using CTSize = policies::CompileTimeSize<SetPosition, TOSET_SIZE>;
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
void generateIncrementingRelations(VecType* begins, VecType* offsets)
{
  VecType& beginsVec = *begins;
  VecType& offsetsVec = *offsets;

  auto curIdx = SetPosition();

  for(auto i = 0; i < FROMSET_SIZE; ++i)
  {
    beginsVec.push_back(curIdx);
    for(auto j = 0; j < elementCardinality(i); ++j)
    {
      offsetsVec.push_back(relationData(i, j));
      ++curIdx;
    }
  }
  beginsVec.push_back(curIdx);
}

}  // namespace

TEST(slam_set_relation_map, access_pattern)
{
  SLIC_INFO("Testing accessing relation data.");

  IndexVec relOffsets, relIndices;
  generateIncrementingRelations(&relOffsets, &relIndices);

  RangeSetType fromSet(FROMSET_SIZE), toSet(TOSET_SIZE);
  StaticVariableRelationType incrementingRel(&fromSet, &toSet);
  incrementingRel.bindBeginOffsets(fromSet.size(), &relOffsets);
  incrementingRel.bindIndices(relIndices.size(), &relIndices);

  // Note: Nothing requires the relations elements to be unique
  //  -- the relation can still be valid with duplicates
  EXPECT_TRUE(incrementingRel.isValid(true));

  SLIC_INFO("-- Looking at relation's stored values...");
  for(auto fromPos = SetPosition(); fromPos < fromSet.size(); ++fromPos)
  {
    SLIC_INFO("--Inspecting element " << fromSet[fromPos] << " in position " << fromPos
                                      << " of first set.");

    for(auto idx = 0; idx < incrementingRel.size(fromPos); ++idx)
    {
      auto posInToSet_actual = incrementingRel[fromPos][idx];
      auto posInToSet_expected = relationData(fromPos, idx);
      EXPECT_EQ(posInToSet_expected, posInToSet_actual);

      SLIC_INFO("-- \t pos: " << idx << " ToSet position: " << incrementingRel[fromPos][idx]
                              << " ToSet element " << toSet[incrementingRel[fromPos][idx]]);
      ;
    }
  }
  SLIC_INFO("done.");
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
