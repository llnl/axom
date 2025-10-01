// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_set_PositionSet.cpp
 *
 * This file tests PositionSets within Slam.  These are OrderedSets with a
 * runtime size, but no offsets, stride or indirection.
 */

#include <iterator>
#include <sstream>  // for std::stringstream

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/slam/RangeSet.hpp"  // for PositionSet
#include "axom/slam/Utilities.hpp"

namespace
{
using SetType = axom::slam::PositionSet<>;
using SetPosition = SetType::PositionType;
using SetElement = SetType::ElementType;

constexpr SetPosition MAX_SET_SIZE = 10;

}  // end anonymous namespace

TEST(slam_set_positionset, construct_valid)
{
  SLIC_INFO("Testing that we can construct a PositionSet from an integer");

  SetType s(MAX_SET_SIZE);

  EXPECT_TRUE(s.isValid());

  EXPECT_EQ(MAX_SET_SIZE, s.size());

  if(MAX_SET_SIZE > SetPosition())
  {
    EXPECT_FALSE(s.empty());
  }
}

TEST(slam_set_positionset, construct_empty)
{
  SLIC_INFO("Testing that we can construct PositionSets of size 0");

  SetType s1;
  EXPECT_TRUE(s1.isValid());
  EXPECT_EQ(0, s1.size());
  EXPECT_TRUE(s1.empty());

  SetType s2(0);
  EXPECT_TRUE(s2.isValid());
  EXPECT_EQ(0, s2.size());
  EXPECT_TRUE(s2.empty());

  EXPECT_EQ(s1, s2);
}

TEST(slam_set_positionset, construct_strictly_positive)
{
  SLIC_INFO("Testing that PositionSets with negative size are invalid");

  const int NEGATIVE_NUMBER = -2;
  SetType s(NEGATIVE_NUMBER);
  EXPECT_FALSE(s.isValid());
}

TEST(slam_set_positionset, construct_set_builder)
{
  SLIC_INFO("Testing construction of PositionSets using SetBuilders");

  using SetBuilder = SetType::SetBuilder;
  SetBuilder builder = SetBuilder().size(MAX_SET_SIZE);

  SetType s(builder);
  EXPECT_TRUE(s.isValid());
  EXPECT_EQ(MAX_SET_SIZE, s.size());

  const int ZERO_OFFSET = 0;
  const int DEFAULT_STRIDE = 1;
  SetBuilder ok_builder = SetBuilder()            //
                            .size(MAX_SET_SIZE)   //
                            .offset(ZERO_OFFSET)  //
                            .stride(DEFAULT_STRIDE);
  SetType s2(ok_builder);
  EXPECT_TRUE(s2.isValid());
  EXPECT_EQ(MAX_SET_SIZE, s2.size());

  // The two sets should be equal
  EXPECT_EQ(s, s2);

#ifdef AXOM_DEBUG
  // Using inappropriate SetBuilder features generates an assert failure
  SLIC_INFO("Cannot construct a PositionSets with an invalid SetBuilder");

  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  const int NON_ZERO_OFFSET = 3;
  EXPECT_DEATH_IF_SUPPORTED(SetType(SetBuilder()           //
                                      .size(MAX_SET_SIZE)  //
                                      .offset(NON_ZERO_OFFSET)),
                            "");
#endif
}

TEST(slam_set_positionset, construct_runtime)
{
  SLIC_INFO("Testing that PositionSets can be constructed with runtime sizes");

  for(int i = 10; i < 15; ++i)
  {
    // mark as volatile to ensure compiler doesn't see this at compile time
    volatile int v = i;
    SetType s(v);

    EXPECT_TRUE(s.isValid());
    EXPECT_EQ(i, s.size());
  }
}

TEST(slam_set_positionset, iterate)
{
  SetType s(MAX_SET_SIZE);

  SLIC_INFO("Iterating through set of size " << s.size());

  SLIC_INFO("Using random access -- operator[]");
  {
    std::stringstream sstr;
    for(SetPosition pos = SetPosition(); pos < s.size(); ++pos)
    {
      SetElement elt = static_cast<SetElement>(pos);
      EXPECT_EQ(elt, s[pos]);
      sstr << s[pos] << "\t";
    }
    SLIC_INFO("Element of slam set using operator[]:\n" << sstr.str());

    // same test, using range-for over OrderedSet::positions()
    int count = 0;
    for(auto pos : s.positions())
    {
      SetElement elt = static_cast<SetElement>(pos);
      EXPECT_EQ(elt, s[pos]);
      EXPECT_EQ(elt, s.at(pos));
      EXPECT_EQ(s[pos], s.at(pos));
      ++count;
    }
    EXPECT_EQ(s.size(), count);
  }

  SLIC_INFO("Using checked random access -- at()");
  {
    std::stringstream sstr;
    for(SetPosition pos = SetPosition(); pos < s.size(); ++pos)
    {
      SetElement elt = static_cast<SetElement>(pos);
      EXPECT_EQ(elt, s.at(pos));

      sstr << s[pos] << "\t";
    }
    SLIC_INFO("Element of slam set using at():\n" << sstr.str());
  }

  SLIC_INFO("Using iterators begin/end");
  {
    // also tests default constructor and operator=
    std::stringstream sstr;

    SetType::iterator it, itEnd;
    EXPECT_EQ(it, itEnd);

    int count = 0;
    for(it = s.begin(), itEnd = s.end(); it != itEnd; ++it)
    {
      EXPECT_EQ(std::distance(s.begin(), it), *it);
      ++count;

      sstr << *it << "\t";
    }
    EXPECT_EQ(s.size(), count);

    SLIC_INFO("Element of slam set using iterators:\n" << sstr.str());
  }
}

TEST(slam_set_positionset, out_of_bounds_at)
{
  SLIC_INFO("Testing out of bounds access using at() " << "-- code is expected to assert and die.");
  SetType s(MAX_SET_SIZE);

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode, so this test will only fail
  // in debug mode

  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED(s.at(MAX_SET_SIZE), "");
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif

  SLIC_INFO("done.");
}

TEST(slam_set_positionset, out_of_bounds_bracket)
{
  SLIC_INFO("Testing out of bounds access using operator[] "
            << "-- code is expected to assert and die.");

  SetType s(MAX_SET_SIZE);

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode,
  // so this test will only fail in debug mode
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED(s[MAX_SET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif

  SLIC_INFO("done.");
}

TEST(slam_set_positionset, awkward_resize)
{
  SLIC_INFO("This test shows that we can modify the underlying size of a "
            << " PositionSet through its RuntimeSizePolicy parent class");

  SetType s1;
  EXPECT_TRUE(s1.isValid());
  EXPECT_EQ(0, s1.size());
  EXPECT_TRUE(s1.empty());
  SLIC_INFO("Original size of set: " << s1.size());

  // Resize the set after construction by casting to its SizePolicyType base
  // class
  // Note that this is awkward, but the only way to do this using the existing
  // API
  static_cast<SetType::SizePolicyType&>(s1).size() = MAX_SET_SIZE;

  EXPECT_TRUE(s1.isValid());
  EXPECT_EQ(MAX_SET_SIZE, s1.size());
  EXPECT_FALSE(s1.empty());
  SLIC_INFO("Modified size of set: " << s1.size());
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
