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

#include "axom/slic.hpp"
#include "axom/core/Array.hpp"
#include "axom/slam/SetBuilders.hpp"
#include "axom/slam/RelationBuilders.hpp"
#include "axom/slam/MapBuilders.hpp"
#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/policies/StridePolicies.hpp"

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

  auto s = slam::make_indirection_set(view);

  // Element type double was deduced from the view.
  static_assert(std::is_same_v<decltype(s), slam::ArrayViewIndirectionSet<Pos, double>>,
                "make_indirection_set(ArrayView) deduces ArrayViewIndirectionSet<.., double>");
  static_assert(std::is_same_v<typename decltype(s)::PositionType, Pos>,
                "make_indirection_set(ArrayView) uses DefaultPositionType by default");

  ASSERT_EQ(s.size(), 4);
  EXPECT_TRUE(s.isValid());
  EXPECT_DOUBLE_EQ(s[0], 10.0);
  EXPECT_DOUBLE_EQ(s[3], 40.0);
}

TEST(slam_make_helpers, make_vector_set_deduces_element_type)
{
  std::vector<int> vec {7, 8, 9};
  auto s = slam::make_indirection_set(vec);

  static_assert(std::is_same_v<decltype(s), slam::VectorIndirectionSet<Pos, int>>,
                "make_indirection_set(vector) deduces VectorIndirectionSet<.., int>");
  static_assert(std::is_same_v<typename decltype(s)::PositionType, Pos>,
                "make_indirection_set(vector) uses DefaultPositionType by default");

  ASSERT_EQ(s.size(), 3);
  EXPECT_TRUE(s.isValid());
  EXPECT_EQ(s[0], 7);
  EXPECT_EQ(s[2], 9);
}

TEST(slam_make_helpers, make_array_set_deduces_element_type)
{
  axom::Array<double> data {1.0, 2.0, 3.0};
  auto s = slam::make_indirection_set(data);

  static_assert(std::is_same_v<decltype(s), slam::ArrayViewIndirectionSet<Pos, double>>,
                "make_indirection_set(Array) deduces ArrayViewIndirectionSet<.., double>");
  static_assert(std::is_same_v<typename decltype(s)::PositionType, Pos>,
                "make_indirection_set(Array) uses DefaultPositionType by default");

  ASSERT_EQ(s.size(), 3);
  EXPECT_TRUE(s.isValid());
  EXPECT_DOUBLE_EQ(s[0], 1.0);
  EXPECT_DOUBLE_EQ(s[2], 3.0);

  data[1] = 20.0;
  EXPECT_DOUBLE_EQ(s[1], 20.0);

  s[2] = 30.0;
  EXPECT_DOUBLE_EQ(data[2], 30.0);
}

TEST(slam_make_helpers, make_multidimensional_array_set_uses_flat_storage)
{
  const axom::StackArray<axom::IndexType, 2> shape {2, 3};
  axom::Array<int, 2> data(shape);

  for(axom::IndexType i = 0; i < data.size(); ++i)
  {
    data.flatIndex(i) = static_cast<int>(10 + i);
  }

  auto s = slam::make_indirection_set(data);

  static_assert(std::is_same_v<decltype(s), slam::ArrayViewIndirectionSet<Pos, int>>,
                "multidimensional Array helper returns flat ArrayViewIndirectionSet");

  ASSERT_EQ(s.size(), 6);
  EXPECT_TRUE(s.isValid());
  EXPECT_EQ(s[0], 10);
  EXPECT_EQ(s[4], data.flatIndex(4));

  s[5] = 42;
  EXPECT_EQ(data.flatIndex(5), 42);
}

TEST(slam_make_helpers, make_carray_set)
{
  int data[3] = {2, 4, 6};
  auto s = slam::make_indirection_set(data, 3);

  static_assert(std::is_same_v<decltype(s), slam::CArrayIndirectionSet<Pos, int>>,
                "make_indirection_set(C array) deduces CArrayIndirectionSet<.., int>");
  static_assert(std::is_same_v<typename decltype(s)::PositionType, Pos>,
                "make_indirection_set(C array) uses DefaultPositionType by default");

  ASSERT_EQ(s.size(), 3);
  EXPECT_TRUE(s.isValid());
  EXPECT_EQ(s[1], 4);
}

TEST(slam_make_helpers, make_array_view_set_explicit_position_type)
{
  axom::Array<float> data {1.5f, 2.5f};
  auto view = data.view();

  // Position type can be supplied explicitly as the leading template argument.
  auto s = slam::make_indirection_set<std::int64_t>(view);
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

TEST(slam_make_helpers, make_variable_relation_carray_buffers)
{
  auto fromSet = slam::make_range_set(3);
  auto toSet = slam::make_range_set(5);

  Pos begins[4] = {0, 2, 3, 5};
  Pos indices[5] = {1, 2, 3, 0, 4};

  auto rel = slam::make_variable_relation(&fromSet, &toSet, begins, Pos {4}, indices, Pos {5});

  static_assert(std::is_same_v<typename decltype(rel)::CardinalityPolicy::BeginsIndirectionPolicy,
                               slam::policies::CArrayIndirection<Pos, Pos>>,
                "raw-pointer begins remain C-array backed");
  static_assert(std::is_same_v<typename decltype(rel)::IndicesIndirectionPolicy,
                               slam::policies::CArrayIndirection<Pos, Pos>>,
                "raw-pointer indices remain C-array backed");

  EXPECT_TRUE(rel.isValid());
  EXPECT_EQ(rel.size(0), 2);
  EXPECT_EQ(rel.size(1), 1);
  EXPECT_EQ(rel.size(2), 2);

  auto rel0 = rel[0];
  ASSERT_EQ(rel0.size(), 2);
  EXPECT_EQ(rel0[0], 1);
  EXPECT_EQ(rel0[1], 2);
}

TEST(slam_make_helpers, make_variable_relation_carray_rejects_short_begins_size)
{
  auto fromSet = slam::make_range_set(3);
  auto toSet = slam::make_range_set(5);

  // begins claims 3 offsets for a size-3 from-set, but needs 4
  // make_variable_relation asserts this invariant at construction in debug builds
  // in release the check compiles out and the malformed relation is instead caught by isValid().
  Pos begins[4] = {0, 2, 3, 3};
  Pos indices[3] = {1, 2, 3};

#ifdef AXOM_DEBUG
  EXPECT_DEATH_IF_SUPPORTED(
    slam::make_variable_relation(&fromSet, &toSet, begins, Pos {3}, indices, Pos {3}),
    "");
#else
  auto rel = slam::make_variable_relation(&fromSet, &toSet, begins, Pos {3}, indices, Pos {3});
  EXPECT_FALSE(rel.isValid());
#endif
}

TEST(slam_make_helpers, make_constant_relation_rejects_undersized_indices)
{
  auto fromSet = slam::make_range_set(3);
  auto toSet = slam::make_range_set(5);

  // A stride-2 constant relation over a size-3 from-set needs 6 indices but we supply 4 here.
  // make_constant_relation asserts the exact size at construction in debug builds
  // the check compiles out in release builds.
  Pos indices[4] = {0, 1, 2, 3};

#ifdef AXOM_DEBUG
  EXPECT_DEATH_IF_SUPPORTED(slam::make_constant_relation(&fromSet, &toSet, Pos {2}, indices, Pos {4}),
                            "");
#else
  SLIC_INFO("Skipped constant-relation size assertion check in release mode.");
#endif
}

TEST(slam_make_helpers, make_variable_relation_axom_array_buffers)
{
  auto fromSet = slam::make_range_set(3);
  auto toSet = slam::make_range_set(5);

  axom::Array<Pos> begins {0, 2, 3, 5};
  axom::Array<Pos> indices {1, 2, 3, 0, 4};

  auto rel = slam::make_variable_relation(&fromSet, &toSet, begins, indices);

  EXPECT_TRUE(rel.isValid());
  EXPECT_EQ(rel.size(0), 2);
  EXPECT_EQ(rel.size(1), 1);
  EXPECT_EQ(rel.size(2), 2);
}

TEST(slam_make_helpers, make_variable_relation_indices_are_to_set_positions)
{
  // toSet has element type != position type; relation indices should still be positions.
  auto fromSet = slam::make_range_set(2);

  axom::Array<double> values {10., 20., 30.};
  auto toSet = slam::make_indirection_set(values.view());  // ArrayViewIndirectionSet<Pos, double>

  std::vector<Pos> begins {0, 2, 3};
  std::vector<Pos> indices {0, 2, 1};  // positions into toSet

  auto rel = slam::make_variable_relation(&fromSet, &toSet, begins, indices);

  static_assert(std::is_same_v<typename decltype(rel)::SetElement, Pos>,
                "relation element type defaults to ToSet::PositionType");

  auto r0 = rel[0];
  ASSERT_EQ(r0.size(), 2);
  EXPECT_DOUBLE_EQ(toSet[r0[0]], 10.0);
  EXPECT_DOUBLE_EQ(toSet[r0[1]], 30.0);
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

TEST(slam_make_helpers, make_constant_relation_runtime_stride_array_view)
{
  auto fromSet = slam::make_range_set(2);
  auto toSet = slam::make_range_set(5);

  axom::Array<Pos> indices {0, 4, 1, 3};
  auto rel = slam::make_constant_relation(&fromSet, &toSet, Pos {2}, indices.view());

  EXPECT_TRUE(rel.isValid());
  EXPECT_EQ(rel.size(0), 2);
  EXPECT_EQ(rel.size(1), 2);
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

TEST(slam_make_helpers, make_constant_relation_compile_time_stride_vector_backed)
{
  auto fromSet = slam::make_range_set(2);
  auto toSet = slam::make_range_set(5);

  std::vector<Pos> indices {0, 4, 1, 3};
  auto rel = slam::make_constant_relation_ct<2>(&fromSet, &toSet, indices);

  static_assert(std::is_same_v<typename decltype(rel)::IndicesIndirectionPolicy,
                               slam::policies::STLVectorIndirection<Pos, Pos>>,
                "make_constant_relation_ct(vector) keeps STLVectorIndirection backing");

  EXPECT_TRUE(rel.isValid());
  EXPECT_EQ(rel.size(0), 2);
  EXPECT_EQ(rel.size(1), 2);

  auto r1 = rel[1];
  ASSERT_EQ(r1.size(), 2);
  EXPECT_EQ(r1[0], 1);
  EXPECT_EQ(r1[1], 3);
}

TEST(slam_make_helpers, make_constant_relation_compile_time_stride_axom_array_backed)
{
  auto fromSet = slam::make_range_set(2);
  auto toSet = slam::make_range_set(5);

  axom::Array<Pos> indices {0, 4, 1, 3};
  auto rel = slam::make_constant_relation_ct<2>(&fromSet, &toSet, indices);

  static_assert(std::is_same_v<typename decltype(rel)::IndicesIndirectionPolicy,
                               slam::policies::ArrayIndirection<Pos, Pos>>,
                "make_constant_relation_ct(axom::Array) keeps ArrayIndirection backing");

  EXPECT_TRUE(rel.isValid());
  EXPECT_EQ(rel.size(0), 2);
  EXPECT_EQ(rel.size(1), 2);

  auto r0 = rel[0];
  ASSERT_EQ(r0.size(), 2);
  EXPECT_EQ(r0[0], 0);
  EXPECT_EQ(r0[1], 4);
}

TEST(slam_make_helpers, make_constant_relation_runtime_stride_carray)
{
  // tests make_constant_relation helper function for C-style arrays

  auto fromSet = slam::make_range_set(3);
  auto toSet = slam::make_range_set(5);

  // from { 0 -> {1,2}; 1 -> {3,4}; 2 -> {0,2} }
  Pos indices[6] = {1, 2, 3, 4, 0, 2};

  // C-array backing: raw pointer + element count.
  auto rel = slam::make_constant_relation(&fromSet, &toSet, Pos {2}, indices, Pos {6});

  EXPECT_TRUE(rel.isValid());
  EXPECT_EQ(rel.size(0), 2);
  EXPECT_EQ(rel.size(1), 2);
  EXPECT_EQ(rel.size(2), 2);

  auto r2 = rel[2];
  ASSERT_EQ(r2.size(), 2);
  EXPECT_EQ(r2[0], 0);
  EXPECT_EQ(r2[1], 2);
}

TEST(slam_make_helpers, make_constant_relation_runtime_stride_axom_array)
{
  // tests make_constant_relation helper function for axom::Array

  auto fromSet = slam::make_range_set(3);
  auto toSet = slam::make_range_set(5);

  // from { 0 -> {1,2}; 1 -> {3,4}; 2 -> {0,2} }
  axom::Array<Pos> indices {1, 2, 3, 4, 0, 2};

  // axom::Array backing (by reference, distinct from the ArrayView overload).
  auto rel = slam::make_constant_relation(&fromSet, &toSet, Pos {2}, indices);

  EXPECT_TRUE(rel.isValid());
  EXPECT_EQ(rel.size(0), 2);
  EXPECT_EQ(rel.size(1), 2);
  EXPECT_EQ(rel.size(2), 2);

  auto r0 = rel[0];
  ASSERT_EQ(r0.size(), 2);
  EXPECT_EQ(r0[0], 1);
  EXPECT_EQ(r0[1], 2);
}

//------------------------------------------------------------------------------
// Test for make_map helpers
//------------------------------------------------------------------------------

TEST(slam_make_helpers, make_map_stride_one_array_view)
{
  auto set = slam::make_range_set(4);
  using SetT = decltype(set);

  axom::Array<double> data {10., 20., 30., 40.};
  auto m = slam::make_map(&set, data.view());

  using Indirection = slam::policies::ArrayViewIndirection<SetT::PositionType, double>;
  using Stride = slam::policies::StrideOne<SetT::PositionType>;
  static_assert(std::is_same_v<decltype(m), slam::Map<double, SetT, Indirection, Stride>>,
                "make_map(set, ArrayView) yields a stride-one ArrayView-backed Map");

  EXPECT_TRUE(m.isValid());
  EXPECT_EQ(m.size(), set.size());
  EXPECT_EQ(m.stride(), 1);
  EXPECT_DOUBLE_EQ(m[0], 10.);
  EXPECT_DOUBLE_EQ(m[3], 40.);

  // The map is a non-owning view: writes round-trip through the backing storage.
  m[2] = 33.;
  EXPECT_DOUBLE_EQ(data[2], 33.);
  data[1] = 22.;
  EXPECT_DOUBLE_EQ(m[1], 22.);
}

TEST(slam_make_helpers, make_map_runtime_stride_array_view)
{
  auto set = slam::make_range_set(2);

  // 2 elements, each with stride/numComp 3 => 6 backing values, laid out element-major:
  //   element 0 -> {1, 2, 3}, element 1 -> {4, 5, 6}
  axom::Array<double> data {1., 2., 3., 4., 5., 6.};
  auto m = slam::make_map(&set, Pos {3}, data.view());

  EXPECT_TRUE(m.isValid());
  EXPECT_EQ(m.size(), 2);
  EXPECT_EQ(m.stride(), 3);

  // operator()(elem, comp) reads flat index elem*stride + comp.
  EXPECT_DOUBLE_EQ(m(0, 0), 1.);
  EXPECT_DOUBLE_EQ(m(0, 2), 3.);
  EXPECT_DOUBLE_EQ(m(1, 0), 4.);
  EXPECT_DOUBLE_EQ(m(1, 2), 6.);

  // operator[] is the flat accessor over size()*numComp() values.
  EXPECT_DOUBLE_EQ(m[4], 5.);
}

TEST(slam_make_helpers, make_map_ct_compile_time_stride_array_view)
{
  auto set = slam::make_range_set(2);
  using SetT = decltype(set);

  axom::Array<double> data {10., 11., 20., 21.};
  auto m = slam::make_map_ct<2>(&set, data.view());

  using Indirection = slam::policies::ArrayViewIndirection<SetT::PositionType, double>;
  using Stride = slam::policies::CompileTimeStride<SetT::PositionType, 2>;
  static_assert(std::is_same_v<decltype(m), slam::Map<double, SetT, Indirection, Stride>>,
                "make_map_ct yields an ArrayView-backed compile-time-stride Map");

  EXPECT_TRUE(m.isValid());
  EXPECT_EQ(m.size(), 2);
  EXPECT_EQ(m.stride(), 2);
  EXPECT_DOUBLE_EQ(m(0, 0), 10.);
  EXPECT_DOUBLE_EQ(m(1, 1), 21.);
}

TEST(slam_make_helpers, make_map_raw_pointer_overloads_size_the_view)
{
  auto set = slam::make_range_set(3);

  // stride-one raw pointer overload: view sized to set->size().
  double one[3] = {7., 8., 9.};
  auto m1 = slam::make_map(&set, one);
  EXPECT_TRUE(m1.isValid());
  EXPECT_EQ(m1.stride(), 1);
  EXPECT_DOUBLE_EQ(m1[2], 9.);

  // strided raw pointer overload: view sized to set->size() * stride.
  double strided[6] = {1., 2., 3., 4., 5., 6.};
  auto m2 = slam::make_map(&set, Pos {2}, strided);
  EXPECT_TRUE(m2.isValid());
  EXPECT_EQ(m2.stride(), 2);
  EXPECT_DOUBLE_EQ(m2(2, 1), 6.);

  // compile-time-stride raw pointer overload.
  auto m3 = slam::make_map_ct<2>(&set, strided);
  EXPECT_TRUE(m3.isValid());
  EXPECT_EQ(m3.stride(), 2);
  EXPECT_DOUBLE_EQ(m3(0, 1), 2.);
}

TEST(slam_make_helpers, make_map_undersized_view_is_a_precondition_violation)
{
  // The ArrayView-taking make_map overloads require data.size() == set->size() * stride.
  // An undersized view is checked in debug builds (no-op in release).
  auto set = slam::make_range_set(4);

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode, so these checks are skipped there.
  axom::Array<double> tooSmall {1., 2., 3.};  // need 4 for stride one
  EXPECT_DEATH_IF_SUPPORTED(slam::make_map(&set, tooSmall.view()), "");

  axom::Array<double> tooSmallStrided {1., 2., 3., 4., 5.};  // need 8 for stride two
  EXPECT_DEATH_IF_SUPPORTED(slam::make_map(&set, Pos {2}, tooSmallStrided.view()), "");
  EXPECT_DEATH_IF_SUPPORTED(slam::make_map_ct<2>(&set, tooSmallStrided.view()), "");
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

  // Construction-precondition tests below use death tests in debug builds.
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  result = RUN_ALL_TESTS();

  return result;
}
