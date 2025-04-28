// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_map_SubMap.cpp
 *
 * \brief Unit tests for Slam's SubMap
 */

#include <iterator>
#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/slam/Utilities.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/slam/SubMap.hpp"

#include <type_traits>

namespace
{
namespace slam = axom::slam;

using SetBase = slam::Set<>;
using PositionType = SetBase::PositionType;
using ElementType = SetBase::ElementType;

using RangeSetType = slam::RangeSet<PositionType, ElementType>;

template <typename T>
using SuperMap = slam::Map<T, SetBase>;

using OrderedSetType =
  axom::slam::OrderedSet<PositionType,
                         ElementType,
                         slam::policies::RuntimeSize<PositionType>,
                         slam::policies::ZeroOffset<PositionType>,
                         slam::policies::StrideOne<PositionType>,
                         slam::policies::STLVectorIndirection<PositionType, ElementType>>;

constexpr double multFac = 1.0001;

PositionType const MAX_SET_SIZE = 10;

template <typename T>
T getValue(int idx)
{
  return static_cast<T>(idx * multFac);
}

//A struct to construct the Map for testing
template <typename T>
struct MapForTest
{
  MapForTest(int size)
    : set_data(size)
    , s(OrderedSetType::SetBuilder().size(size).data(&set_data))
    , m(&s)
  {
    SLIC_INFO("Initializing set of size " << s.size() << " and '"
                                          << slam::util::TypeToString<T>::to_string()
                                          << "' map on the set ");

    for(auto i : s.positions())
    {
      s[i] = 100 + i;
      m[i] = getValue<T>(i);
    }

    EXPECT_TRUE(s.isValid());
    EXPECT_TRUE(m.isValid());
  }

  std::vector<PositionType> set_data;
  OrderedSetType s;
  SuperMap<T> m;
};

}  // namespace

TEST(slam_map, construct_empty_subsetmap)
{
  slam::SubMap<SuperMap<int>, RangeSetType> m;
  EXPECT_TRUE(m.isValid(true));
}

template <typename T>
bool constructAndTestSubMap()
{
  using MapType = SuperMap<T>;

  MapForTest<T> mft(MAX_SET_SIZE);
  MapType& m = mft.m;
  OrderedSetType& s = mft.s;

  int submapOffset = 3;
  int submapSize = 5;
  {
    using SubMapType = slam::SubMap<MapType, RangeSetType>;

    SLIC_INFO("Creating the Subset.");
    RangeSetType ss(submapOffset, submapOffset + submapSize);
    SubMapType ssm(&m, ss);
    EXPECT_TRUE(ssm.isValid(true));

    SLIC_INFO("Checking the elements.");
    for(auto idx = 0; idx < ssm.size(); ++idx)
    {
      auto expVal = getValue<T>(submapOffset + idx);
      EXPECT_EQ(expVal, ssm[idx]);
    }

    SLIC_INFO("Checking the elements using SubMap iterator.");
    int cnt = 0;
    for(auto it = ssm.begin(); it != ssm.end(); ++it, ++cnt)
    {
      // Check iterator's .index() function
      {
        auto subMapElt = it.index();

        auto setIdx = ss[cnt];
        auto expSetElt = s[setIdx];
        EXPECT_EQ(expSetElt, subMapElt);
      }

      // Check iterator's value access functions
      {
        auto expectedValue = getValue<T>(submapOffset + cnt);
        EXPECT_EQ(expectedValue, *it);
        EXPECT_EQ(expectedValue, it[0]);

        auto expVal = getValue<T>(submapOffset + cnt);
        EXPECT_EQ(expVal, expectedValue);
      }
    }
  }

  {
    using SubMapType = slam::SubMap<MapType, OrderedSetType>;

    SLIC_INFO("Creating Subset 2");
    std::vector<PositionType> subset_indices_data(submapSize);
    for(int i = 0; i < submapSize; i++)
    {
      subset_indices_data[i] = i * 2;
    }
    OrderedSetType subset_indices = OrderedSetType::SetBuilder().size(5).data(&subset_indices_data);

    SubMapType ssm(&m, subset_indices);

    SLIC_INFO("Checking the elements.");
    for(PositionType idx = 0; idx < submapSize; ++idx)
    {
      EXPECT_EQ(ssm[idx], getValue<T>(subset_indices[idx]));
      EXPECT_EQ(ssm.value(idx), getValue<T>(subset_indices[idx]));
      EXPECT_EQ(ssm.index(idx), s[subset_indices[idx]]);
    }

    SLIC_INFO("Checking the elements using SubMap range iterator.");
    {
      int cnt = 0;
      for(auto it = ssm.set_begin(); it != ssm.set_end(); ++it, ++cnt)
      {
        auto expectedValue = getValue<T>(subset_indices[cnt]);
        EXPECT_EQ(expectedValue, (*it)[0]);
        EXPECT_EQ(expectedValue, it(0));
        EXPECT_EQ(expectedValue, it.value(0));
        EXPECT_EQ(it.index(), s[subset_indices[cnt]]);
      }
      EXPECT_EQ(cnt, subset_indices.size());
    }
  }
  return true;
}

TEST(slam_map, construct_int_submap) { EXPECT_TRUE(constructAndTestSubMap<int>()); }

TEST(slam_map, construct_double_submap) { EXPECT_TRUE(constructAndTestSubMap<double>()); }

template <typename T>
bool constructBySubMap()
{
  //This tests modifying the values in the original map via a Submap
  //Create a Map, then create a SubMap on the Map, negate all values covered by
  //the Submap, then check the values are negative in the SubMap range.
  using MapType = SuperMap<T>;
  using SubMapType = slam::SubMap<MapType, RangeSetType>;

  MapForTest<T> mft(MAX_SET_SIZE);
  MapType& m = mft.m;

  int submapOffset = 3;
  int submapSize = 5;

  SLIC_INFO("Creating the Subset.");
  RangeSetType ss(submapOffset, submapOffset + submapSize);
  SubMapType ssm(&m, ss);
  EXPECT_TRUE(m.isValid(true));

  SLIC_INFO("Negating elements");
  for(PositionType idx = 0; idx < submapSize; ++idx)
  {
    ssm[idx] = -getValue<T>(submapOffset + idx);
  }

  SLIC_INFO("Checking the elements.");
  for(PositionType idx = 0; idx < m.size(); ++idx)
  {
    T val = getValue<T>(idx);
    if(idx >= submapOffset && idx < submapOffset + submapSize)
    {
      val = -getValue<T>(idx);
    }
    EXPECT_EQ(m[idx], val);
  }

  SLIC_INFO("Checking elements");

  return true;
}

TEST(slam_map, construct_with_int_submap) { EXPECT_TRUE(constructBySubMap<int>()); }

TEST(slam_map, construct_with_double_submap) { EXPECT_TRUE(constructBySubMap<double>()); }

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
#ifdef AXOM_DEBUG
  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
#endif

  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  int result = RUN_ALL_TESTS();

  return result;
}
