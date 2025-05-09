// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_map_DynamicMap.cpp
 *
 * \brief Unit tests for slam's DynamicMap
 */

#include <iterator>
#include "gtest/gtest.h"

#include "axom/slam/Utilities.hpp"
#include "axom/slam/DynamicSet.hpp"
#include "axom/slam/DynamicMap.hpp"

#include "axom/slic.hpp"

namespace
{
namespace slam = axom::slam;

using SetPosition = slam::DefaultPositionType;
using SetElement = slam::DefaultElementType;
using SetType = slam::DynamicSet<SetPosition, SetElement>;

constexpr SetPosition MAX_SET_SIZE = 10;

using IntMap = slam::DynamicMap<SetType, int>;
using RealMap = slam::DynamicMap<SetType, double>;

}  // end anonymous namespace

TEST(slam_map, construct_empty_map)
{
  IntMap m;

  EXPECT_TRUE(m.isValid(true));
  EXPECT_EQ(nullptr, m.set());
}

TEST(slam_map, construct_from_set_int)
{
  SetType s(MAX_SET_SIZE);

  IntMap m(&s);
  EXPECT_TRUE(m.isValid(true));

  EXPECT_EQ(MAX_SET_SIZE, m.size());
  EXPECT_EQ(s.size(), m.size());
}

TEST(slam_map, construct_from_set_real)
{
  SetType s(MAX_SET_SIZE);

  const double DEFAULT_VAL = -1.0;
  RealMap m(&s, DEFAULT_VAL);
  EXPECT_TRUE(m.isValid(true));
  EXPECT_EQ(&s, m.set());

  EXPECT_EQ(MAX_SET_SIZE, m.size());
  EXPECT_EQ(s.size(), m.size());

  for(int i = 0; i < MAX_SET_SIZE; ++i)
  {
    EXPECT_EQ(DEFAULT_VAL, m[i]);
  }

  // Check that we can update the data
  const double UPDATE_VAL = 5.5;
  m[2] = UPDATE_VAL;
  EXPECT_EQ(UPDATE_VAL, m[2]);
}

// Check that resizing the map to a smaller size
// TODO: This shouldn't shrink the data
TEST(slam_map, modify_size)
{
  SetType s(MAX_SET_SIZE);

  IntMap m(&s);

  // Resizing to a larger size increases the map's size
  {
    const int newSize = 2 * MAX_SET_SIZE;
    m.resize(newSize);

    EXPECT_EQ(newSize, m.size());
  }

  // Resizing to a larger size increases the map's size
  {
    const int oldSize = m.size();
    const int newSize = MAX_SET_SIZE / 2;
    m.resize(newSize);

    EXPECT_GT(oldSize, newSize);
  }
}

// Check that resizing the map
// TODO: This should resize the underlying set
TEST(slam_map, modify_size_increases_set)
{
  SetType s(MAX_SET_SIZE);
  IntMap m(&s);

  const int newSize = 2 * MAX_SET_SIZE;
  m.resize(newSize);

  EXPECT_EQ(newSize, m.size());
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  int result = RUN_ALL_TESTS();

  return result;
}
