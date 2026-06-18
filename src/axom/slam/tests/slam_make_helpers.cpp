// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file slam_make_helpers.cpp
 *
 * \brief Unit tests for the make* set/relation construction helpers,
 *        which deduce a SLAM type's full policy stack from a buffer or range.
 */

#include "gtest/gtest.h"

#include "axom/core/Array.hpp"
#include "axom/slam/SetBuilders.hpp"
#include "axom/slam/RelationBuilders.hpp"

#include <type_traits>
#include <vector>

namespace
{
namespace slam = axom::slam;
using Pos = slam::DefaultPositionType;
}  // anonymous namespace

TEST(slam_make_helpers, make_range_set_size)
{
  auto s = slam::make_range_set(5);
  EXPECT_EQ(s.size(), 5);
  EXPECT_EQ(s[0], 0);
  EXPECT_EQ(s[4], 4);

  static_assert(std::is_same_v<typename decltype(s)::PositionType, Pos>,
                "make_range_set uses DefaultPositionType by default");
  static_assert(std::is_same_v<typename decltype(s)::ElementType, slam::DefaultElementType>,
                "make_range_set uses DefaultElementType by default");

  // The deduced type is the blessed RangeSet alias.
  static_assert(std::is_same_v<decltype(s), slam::RangeSet<Pos, slam::DefaultElementType>>,
                "make_range_set(size) yields RangeSet");
}

TEST(slam_make_helpers, make_range_set_bounds)
{
  auto s = slam::make_range_set(3, 8);
  EXPECT_EQ(s.size(), 5);
  EXPECT_EQ(s[0], 3);
  EXPECT_EQ(s[4], 7);

  static_assert(std::is_same_v<typename decltype(s)::PositionType, Pos>,
                "make_range_set(lower, upper) uses DefaultPositionType by default");
}

TEST(slam_make_helpers, make_range_set_explicit_position_type)
{
  auto s = slam::make_range_set<std::int64_t>(5);
  static_assert(std::is_same_v<typename decltype(s)::PositionType, std::int64_t>,
                "explicit PosType flows through make_range_set");
  EXPECT_EQ(s.size(), 5);
}

TEST(slam_make_helpers, make_array_view_set_deduces_element_type)
{
  axom::Array<double> data {10., 20., 30., 40.};
  axom::ArrayView<double> view = data.view();

  auto s = slam::make_array_view_set(view);

  // Element type double was deduced from the view.
  static_assert(std::is_same_v<decltype(s), slam::ArrayViewIndirectionSet<Pos, double>>,
                "make_array_view_set deduces ArrayViewIndirectionSet<.., double>");
  static_assert(std::is_same_v<typename decltype(s)::PositionType, Pos>,
                "make_array_view_set uses DefaultPositionType by default");

  ASSERT_EQ(s.size(), 4);
  EXPECT_TRUE(s.isValid());
  EXPECT_DOUBLE_EQ(s[0], 10.0);
  EXPECT_DOUBLE_EQ(s[3], 40.0);
}

TEST(slam_make_helpers, make_vector_set_deduces_element_type)
{
  std::vector<int> vec {7, 8, 9};
  auto s = slam::make_vector_set(vec);

  static_assert(std::is_same_v<decltype(s), slam::VectorIndirectionSet<Pos, int>>,
                "make_vector_set deduces VectorIndirectionSet<.., int>");
  static_assert(std::is_same_v<typename decltype(s)::PositionType, Pos>,
                "make_vector_set uses DefaultPositionType by default");

  ASSERT_EQ(s.size(), 3);
  EXPECT_TRUE(s.isValid());
  EXPECT_EQ(s[0], 7);
  EXPECT_EQ(s[2], 9);
}

TEST(slam_make_helpers, make_carray_set)
{
  int data[3] = {2, 4, 6};
  auto s = slam::make_c_array_set(data, 3);

  static_assert(std::is_same_v<decltype(s), slam::CArrayIndirectionSet<Pos, int>>,
                "make_c_array_set deduces CArrayIndirectionSet<.., int>");
  static_assert(std::is_same_v<typename decltype(s)::PositionType, Pos>,
                "make_c_array_set uses DefaultPositionType by default");

  ASSERT_EQ(s.size(), 3);
  EXPECT_TRUE(s.isValid());
  EXPECT_EQ(s[1], 4);
}

TEST(slam_make_helpers, make_array_view_set_explicit_position_type)
{
  axom::Array<float> data {1.5f, 2.5f};
  auto view = data.view();

  // Position type can be supplied explicitly as the leading template argument.
  auto s = slam::make_array_view_set<std::int64_t>(view);
  static_assert(std::is_same_v<decltype(s), slam::ArrayViewIndirectionSet<std::int64_t, float>>,
                "explicit PosType flows through");
  EXPECT_EQ(s.size(), 2);
}

TEST(slam_make_helpers, make_variable_relation)
{
  // from-set of 3 elements, to-set of 5 elements.
  auto fromSet = slam::make_range_set(3);
  auto toSet = slam::make_range_set(5);

  // CSR layout: element 0 -> {1,2}, element 1 -> {3}, element 2 -> {0,4}
  std::vector<Pos> begins {0, 2, 3, 5};  // size == fromSet.size() + 1
  std::vector<Pos> indices {1, 2, 3, 0, 4};

  auto rel = slam::make_variable_relation(&fromSet, &toSet, begins, indices);

  EXPECT_TRUE(rel.isValid());
  EXPECT_EQ(rel.size(0), 2);
  EXPECT_EQ(rel.size(1), 1);
  EXPECT_EQ(rel.size(2), 2);

  // Spot-check the related indices for from-element 2.
  auto rel2 = rel[2];
  ASSERT_EQ(rel2.size(), 2);
  EXPECT_EQ(rel2[0], 0);
  EXPECT_EQ(rel2[1], 4);
}

TEST(slam_make_helpers, make_constant_relation_runtime_stride)
{
  auto fromSet = slam::make_range_set(3);
  auto toSet = slam::make_range_set(5);

  //from { 0 -> {1,2}; 1 -> {3,4}; from 2 -> {0,2} }
  std::vector<Pos> indices {1, 2, 3, 4, 0, 2};

  auto rel = slam::make_constant_relation(&fromSet, &toSet, Pos {2}, indices);

  EXPECT_TRUE(rel.isValid());
  EXPECT_EQ(rel.size(0), 2);
  EXPECT_EQ(rel.size(1), 2);
  EXPECT_EQ(rel.size(2), 2);

  auto r1 = rel[1];
  ASSERT_EQ(r1.size(), 2);
  EXPECT_EQ(r1[0], 3);
  EXPECT_EQ(r1[1], 4);
}

TEST(slam_make_helpers, make_constant_relation_compile_time_stride)
{
  auto fromSet = slam::make_range_set(2);
  auto toSet = slam::make_range_set(5);

  Pos indices[4] = {0, 4, 1, 3};

  auto rel = slam::make_constant_relation_ct<2>(&fromSet, &toSet, indices, Pos {4});

  EXPECT_TRUE(rel.isValid());
  EXPECT_EQ(rel.size(0), 2);
  EXPECT_EQ(rel.size(1), 2);

  auto r0 = rel[0];
  ASSERT_EQ(r0.size(), 2);
  EXPECT_EQ(r0[0], 0);
  EXPECT_EQ(r0[1], 4);
}
