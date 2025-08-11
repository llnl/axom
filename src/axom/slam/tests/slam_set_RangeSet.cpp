// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_set_rangeset.cpp
 *
 * \brief Unit tests for Slam's RangeSet
 *
 * A RangeSet is a contiguous range of integers from low to high.
 * It has a configurable offset and size but no indirection and a stride of +1
 *
 * Also tests GenericRangeSets, which provide a runtime size and offset,
 * and the user supplies templates for strides, indirection and subsetting.
 */

#include <iterator>
#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/slam/Set.hpp"
#include "axom/slam/RangeSet.hpp"

#include "axom/fmt.hpp"

namespace
{
using SetType = axom::slam::RangeSet<>;
using SetPosition = SetType::PositionType;
using SetElement = SetType::ElementType;

constexpr SetPosition MAX_SIZE = 20;
constexpr SetElement lowerIndex = static_cast<SetElement>(.3 * MAX_SIZE);
constexpr SetElement upperIndex = static_cast<SetElement>(.7 * MAX_SIZE);
constexpr SetElement range = upperIndex - lowerIndex;

}  // end anonymous namespace

TEST(slam_range_set, construct)
{
  SLIC_INFO("Testing construction RangeSets using constructors");

  // Empty set
  SetType s0;
  EXPECT_TRUE(s0.isValid());
  EXPECT_TRUE(s0.empty());
  EXPECT_EQ(0, s0.size());
  EXPECT_EQ(0, s0.offset());

  // Construct using only a size
  SetType s1(range);
  EXPECT_TRUE(s1.isValid());
  EXPECT_EQ(range, s1.size());
  EXPECT_EQ(0, s1.offset());

  // Construct a sized range starting at 0
  SetType s1_b(0, range);
  EXPECT_TRUE(s1_b.isValid());
  EXPECT_EQ(range, s1_b.size());
  EXPECT_EQ(0, s1_b.offset());

  EXPECT_EQ(s1, s1_b);

  // Construct using a lower and upper bound
  SetType s2(lowerIndex, upperIndex);
  EXPECT_TRUE(s2.isValid(true));
  EXPECT_EQ(range, s2.size());
  EXPECT_EQ(lowerIndex, s2.offset());

  EXPECT_EQ(s1.size(), s2.size());

  // Check that negative ranges are fine
  SetType s3(-upperIndex, -lowerIndex);
  EXPECT_TRUE(s3.isValid(true));
  EXPECT_EQ(range, s3.size());
  EXPECT_EQ(-upperIndex, s3.offset());
}

TEST(slam_range_set, set_builder)
{
  SLIC_INFO("Testing construction RangeSets using SetBuilders");

  using SetBuilder = SetType::SetBuilder;

  // Empty set, with and without offset
  SetType s0(SetBuilder().size(0));
  EXPECT_TRUE(s0.isValid());
  EXPECT_TRUE(s0.empty());
  EXPECT_EQ(0, s0.size());
  EXPECT_EQ(0, s0.offset());

  SetType s0_b(SetBuilder().offset(10));
  EXPECT_TRUE(s0_b.isValid());
  EXPECT_TRUE(s0_b.empty());
  EXPECT_EQ(0, s0_b.size());
  EXPECT_EQ(10, s0_b.offset());

  // Construct using only a size
  SetType s1(SetBuilder().size(range));
  EXPECT_TRUE(s1.isValid());
  EXPECT_EQ(range, s1.size());
  EXPECT_EQ(0, s1.offset());

  // Construct a sized range starting at 0
  SetType s1_b(SetBuilder().range(0, range));
  EXPECT_TRUE(s1_b.isValid());
  EXPECT_EQ(range, s1_b.size());
  EXPECT_EQ(0, s1_b.offset());

  EXPECT_EQ(s1, s1_b);

  // Construct using a lower and upper bound
  SetType s2(SetBuilder().range(lowerIndex, upperIndex));
  EXPECT_TRUE(s2.isValid(true));
  EXPECT_EQ(range, s2.size());
  EXPECT_EQ(lowerIndex, s2.offset());

  SetType s2_b(SetBuilder()    //
                 .size(range)  //
                 .offset(lowerIndex));
  EXPECT_TRUE(s2_b.isValid(true));
  EXPECT_EQ(range, s2_b.size());
  EXPECT_EQ(lowerIndex, s2_b.offset());

  EXPECT_EQ(s2, s2_b);

  // Check that negative ranges are fine
  SetType s3(SetBuilder().range(-upperIndex, -lowerIndex));
  EXPECT_TRUE(s3.isValid(true));
  EXPECT_EQ(range, s3.size());
  EXPECT_EQ(-upperIndex, s3.offset());

  SetType s3_b(SetBuilder()    //
                 .size(range)  //
                 .offset(-upperIndex));
  EXPECT_TRUE(s3_b.isValid(true));
  EXPECT_EQ(range, s3_b.size());
  EXPECT_EQ(-upperIndex, s3_b.offset());

  EXPECT_EQ(s3, s3_b);
}

TEST(slam_range_set, iterate)
{
  SetType s(lowerIndex, upperIndex);

  SLIC_INFO("Iterating through range set of size "
            << s.size() << "\n\twith lower element " << lowerIndex << " (included in set)"
            << "\n\twith upper element " << upperIndex << " (not included in set)");
  const SetPosition expectedSize = upperIndex - lowerIndex;
  EXPECT_EQ(expectedSize, s.size());

  SLIC_INFO("Testing random access -- operator[] and at() function");
  {
    std::stringstream sstr;
    for(SetPosition pos = SetPosition(); pos < s.size(); ++pos)
    {
      SetElement expected = pos + lowerIndex;
      EXPECT_EQ(expected, s[pos]);
      EXPECT_EQ(expected, s.at(pos));
      EXPECT_EQ(s[pos], s.at(pos));

      sstr << s[pos] << "\t";
    }
    SLIC_INFO(sstr.str());

    // same test, using range-for over OrderedSet::positions()
    for(auto pos : s.positions())
    {
      SetElement expected = pos + lowerIndex;
      EXPECT_EQ(expected, s[pos]);
      EXPECT_EQ(expected, s.at(pos));
      EXPECT_EQ(s[pos], s.at(pos));
    }
  }

  SLIC_INFO("Testing iterator access");
  {
    std::stringstream sstr;

    using SetIterator = SetType::iterator;

    EXPECT_EQ(SetIterator(s.offset(), s), s.begin());

    //iter++ access
    {
      int idx = 0;
      for(SetIterator it = s.begin(), itEnd = s.end(); it != itEnd; ++it, ++idx)
      {
        SetElement expected = idx + lowerIndex;
        EXPECT_EQ(expected, *it) << "Iterator dereference should be equal "
                                 << "to its translated position in the windowed range set";
        sstr << *it << "\t";
      }
    }

    //iter+n access
    {
      SetIterator beginIter = s.begin();
      for(int idx = 0; idx < s.size(); ++idx)
      {
        SetIterator iter = beginIter + idx;

        //SetPosition position = std::distance(s.begin(), iter);
        SetElement expected = idx + lowerIndex;
        EXPECT_EQ(expected, *iter) << "Iterator dereference should be equal "
                                   << "to its translated position in the windowed range set";
        sstr << *iter << "\t";
      }
    }

    //iter-n access
    {
      SetIterator endIter = s.end();
      for(int idx = 1; idx <= s.size(); ++idx)
      {
        SetIterator iter = endIter - idx;

        //SetPosition position = std::distance(s.begin(), iter);
        SetElement expected = upperIndex - idx;
        EXPECT_EQ(expected, *iter) << "Iterator dereference should be equal "
                                   << "to its translated position in the windowed range set";
        sstr << *iter << "\t";
      }
    }

    //iter+=n access
    {
      for(int idx = 0; idx < s.size(); ++idx)
      {
        SetIterator iter = s.begin();
        iter += idx;
        SetElement expected = idx + lowerIndex;
        EXPECT_EQ(expected, *iter) << "Iterator dereference should be equal "
                                   << "to its translated position in the windowed range set";
        sstr << *iter << "\t";
      }
    }

    //iter-=n access
    {
      for(int idx = 1; idx <= s.size(); ++idx)
      {
        SetIterator iter = s.end();
        iter -= idx;
        SetElement expected = upperIndex - idx;
        EXPECT_EQ(expected, *iter) << "Iterator dereference should be equal "
                                   << "to its translated position in the windowed range set";
        sstr << *iter << "\t";
      }
    }

    //iter[n] access
    {
      SetIterator beginIter = s.begin();
      for(int idx = 0; idx < s.size(); ++idx)
      {
        SetElement expected = idx + lowerIndex;
        EXPECT_EQ(expected, beginIter[idx])
          << "Iterator dereference should be equal "
          << "to its translated position in the windowed range set";
        sstr << beginIter[idx] << "\t";
      }
    }

    //std::distance calculation
    //iter+n access
    {
      int idx = 0;
      for(SetIterator it = s.begin(), itEnd = s.end(); it != itEnd; ++it, ++idx)
      {
        SetElement position = std::distance(s.begin(), it);
        SetElement expected = idx;
        EXPECT_EQ(expected, position) << "Iterator distance should be equal "
                                      << "to the number of times it's been incremented";
        sstr << idx << "\t";
      }
    }

    EXPECT_TRUE(s.isValid(true));

    SLIC_INFO(sstr.str());
  }
}

TEST(slam_range_set, out_of_range)
{
  SetType s(lowerIndex, upperIndex);

  SLIC_INFO("Using random access on invalid address "
            << "-- Note: We are testing for the expected failures.");
#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode,
  // so this test will only fail in debug mode
  EXPECT_DEATH_IF_SUPPORTED(s.at(upperIndex), "");
  EXPECT_DEATH_IF_SUPPORTED(s.at(MAX_SIZE), "");
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif
}

TEST(slam_generic_range_set, virtual_parent_set)
{
  namespace policies = axom::slam::policies;

  using GenericRangeSet =
    axom::slam::GenericRangeSet<SetPosition,
                                SetElement,
                                policies::RuntimeOffset<SetPosition>,
                                policies::StrideOne<SetPosition>,
                                policies::NoIndirection<SetPosition, SetElement>,
                                policies::VirtualParentSubset>;

  using SetBuilder = GenericRangeSet::SetBuilder;

  SLIC_INFO("Generating a parent set, and a subset and checking validity");
  GenericRangeSet parentSet(SetBuilder().size(MAX_SIZE));
  GenericRangeSet childSet(SetBuilder().range(lowerIndex, upperIndex).parent(&parentSet));

  SetType nonChildSet(SetType::SetBuilder().range(lowerIndex, upperIndex));

  EXPECT_TRUE(parentSet.isValid(true));
  EXPECT_TRUE(childSet.isValid(true));
  EXPECT_TRUE(nonChildSet.isValid(true));

  SLIC_INFO("Checking that the child is a subset, "
            << " but not the parent or the non-child windowed set.");
  EXPECT_FALSE(parentSet.isSubset());
  EXPECT_TRUE(childSet.isSubset());
  EXPECT_FALSE(nonChildSet.isSubset());

  SLIC_INFO("Checking that the child set's parent is equal to "
            << " the parent set (according to the equality operator==).");
  EXPECT_EQ(parentSet, *childSet.parentSet());

  // Since the parent is using the virtual subsetting policy,
  // it does not have an operator[]
  for(SetPosition pos = 0; pos > childSet.parentSet()->size(); ++pos)
  {
    EXPECT_EQ(parentSet[pos], childSet.parentSet()->at(pos));
  }

  // Note: Equality is based on Base class Set
  //-- it does not differentiate based on whether a set is a subset of another
  // set
  EXPECT_EQ(childSet, nonChildSet);
}

TEST(slam_generic_range_set, concrete_parent_set)
{
  namespace policies = axom::slam::policies;

  using ParentType = SetType;

  using GenericRangeSet =
    axom::slam::GenericRangeSet<SetPosition,
                                SetElement,
                                policies::RuntimeOffset<SetPosition>,
                                policies::StrideOne<SetPosition>,
                                policies::NoIndirection<SetPosition, SetElement>,
                                policies::ConcreteParentSubset<ParentType>>;

  using SetBuilder = GenericRangeSet::SetBuilder;

  SLIC_INFO("Generating a parent set, and a subset and checking validity");
  ParentType parentSet(ParentType::SetBuilder().size(MAX_SIZE));
  GenericRangeSet childSet(SetBuilder()                      //
                             .range(lowerIndex, upperIndex)  //
                             .parent(&parentSet));

  SetType nonChildSet(SetType::SetBuilder()  //
                        .range(lowerIndex, upperIndex));

  EXPECT_TRUE(parentSet.isValid(true));
  EXPECT_TRUE(childSet.isValid(true));
  EXPECT_TRUE(nonChildSet.isValid(true));

  SLIC_INFO("Checking that the child is a subset, "
            << "but not the parent or the non-child windowed set.");
  EXPECT_FALSE(parentSet.isSubset());
  EXPECT_TRUE(childSet.isSubset());
  EXPECT_FALSE(nonChildSet.isSubset());

  SLIC_INFO("Checking that the child set's parent is equal to "
            << "the parent set (according to the equality operator==).");
  EXPECT_EQ(parentSet, *childSet.parentSet());

  // Since the parent is a concrete OrderedSet, it should have an operator[]
  const ParentType& childParSet = *(childSet.parentSet());
  for(SetPosition pos = 0; pos > childParSet.size(); ++pos)
  {
    EXPECT_EQ(childParSet.at(pos), childParSet[pos]);
    EXPECT_EQ(parentSet[pos], childParSet[pos]);
  }

  // Note: Equality is based on Base class Set -- it does not differentiate based on whether a set is a subset of another set
  EXPECT_EQ(childSet, nonChildSet);
}

TEST(slam_generic_range_set, strided_range_set)
{
  namespace policies = axom::slam::policies;

  using GenericRangeSet = axom::slam::GenericRangeSet<SetPosition,
                                                      SetElement,
                                                      policies::RuntimeOffset<SetPosition>,
                                                      policies::RuntimeStride<SetPosition>>;
  using SetBuilder = GenericRangeSet::SetBuilder;

  for(SetPosition strideVal = -14; strideVal <= 14; ++strideVal)
  {
    const SetPosition localLower = strideVal > 0 ? lowerIndex : upperIndex;
    const SetPosition localUpper = strideVal > 0 ? upperIndex : lowerIndex;

    SLIC_INFO(axom::fmt::format("testing strided_range_set with: lower: {}, upper: {}, stride: {}",
                                localLower,
                                localUpper,
                                strideVal));

    GenericRangeSet stridedSet(SetBuilder()  //
                                 .range(localLower, localUpper)
                                 .stride(strideVal));

    // stride of zero is not valid; check and then continue to next case
    if(strideVal == 0)
    {
      EXPECT_FALSE(stridedSet.isValid());
      continue;
    }

    SLIC_INFO(
      axom::fmt::format("strided_range_set: "
                        "size: {}, stride: {}, offset: {}, first index: {}, last index: {}",
                        stridedSet.size(),
                        stridedSet.stride(),
                        stridedSet.offset(),
                        stridedSet[0],
                        stridedSet[stridedSet.size() - 1]));

    EXPECT_EQ(localLower, stridedSet.offset());
    EXPECT_EQ(strideVal, stridedSet.stride());
    if(axom::utilities::abs(localLower - localUpper) > axom::utilities::abs(strideVal))
    {
      EXPECT_GT(stridedSet.size(), 1);
    }
    EXPECT_GT(stridedSet.size(), 0);

    SetPosition idx = 0;
    for(auto it = stridedSet.begin(); it != stridedSet.end(); ++it, ++idx)
    {
      SLIC_INFO(axom::fmt::format("\tindex: {}, val: {}", it.index(), *it));

      EXPECT_EQ(idx, it.index());
      EXPECT_EQ(stridedSet[idx], stridedSet[it.index()]);

      const auto expVal = localLower + idx * strideVal;
      EXPECT_EQ(expVal, *it);

      if(strideVal > 0)
      {
        EXPECT_GE(*it, localLower);
        EXPECT_LE(*it, localUpper);
      }
      else
      {
        EXPECT_LE(*it, localLower);
        EXPECT_GE(*it, localUpper);
      }
    }
  }
}

//----------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
#ifdef AXOM_DEBUG
  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
#endif

  axom::slic::SimpleLogger logger;
  //axom::slic::debug::checksAreErrors = true;

  result = RUN_ALL_TESTS();

  return result;
}
