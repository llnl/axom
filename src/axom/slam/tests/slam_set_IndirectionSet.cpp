// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_set_IndirectionSet.cpp
 *
 * Most unit tests in this file are parametrized on an Indirection set type
 * using the IndirectionSetTester class
 */

#include "axom/slic.hpp"
#include "axom/slam.hpp"

#include "gtest/gtest.h"

#include <iterator>
#include <sstream>
#include <algorithm>
#include <iterator>

namespace slam = axom::slam;

namespace
{
constexpr int MAX_SET_SIZE = 10;
}

/**
 * Simple test class for IndirectionSets
 * Initializes three vectors of data
 *
 * \note Uses specialized variants of getDataBuffer() to access
 * the data in an appropriated way for the indirectionSet \a TheSet
 */
template <typename TheSet>
class IndirectionSetTester : public ::testing::Test
{
public:
  using SetType = TheSet;
  using IndirectionType = typename SetType::IndirectionPolicyType;
  using BufferType = typename IndirectionType::IndirectionBufferType;
  using BufferPtrType = typename IndirectionType::IndirectionPtrType;

  using PosType = typename SetType::PositionType;
  using ElemType = typename SetType::ElementType;

  virtual void SetUp()
  {
    // Allocate and initialize three incrementing vectors
    // vecs 0 and 1 have the same size; vec 2 is shorter

    mDataVec[0].reserve(MAX_SET_SIZE);
    mDataVec[1].reserve(MAX_SET_SIZE);
    mDataVec[2].reserve(MAX_SET_SIZE / 2);

    for(int i = 0; i < MAX_SET_SIZE; ++i)
    {
      mDataVec[0].push_back(static_cast<ElemType>(i));
      mDataVec[1].push_back(static_cast<ElemType>(i));
    }
    for(int i = 0; i < MAX_SET_SIZE / 2; ++i)
    {
      mDataVec[2].push_back(static_cast<ElemType>(i));
    }

    for(int i = 0; i < 3; ++i)
    {
      mDataArr[i].resize(mDataVec[i].size());
      std::copy(mDataVec[i].begin(), mDataVec[i].end(), mDataArr[i].begin());
    }
  }

  /**
   * Permute the data buffer with index \a index
   * \pre index must be 0, 1 or 2
   */
  void permuteData(int index = 0)
  {
    checkIndex(index);

    for(int i = 0; i < 1000; ++i)
    {
      std::next_permutation(mDataVec[index].begin(), mDataVec[index].end());
      std::next_permutation(mDataArr[index].begin(), mDataArr[index].end());
    }
  }

  /**
   * Get a pointer to the data buffer
   * Specialization for CArrayIndirection
   * \pre index must be 0, 1 or 2
   */
  void getDataBuffer(ElemType*& ptr, int index = 0)
  {
    checkIndex(index);
    ptr = mDataVec[index].data();
  }

  /**
   * Get a pointer to the data buffer
   * Specialization for VectorIndirection
   * \pre index must be 0, 1 or 2
   */
  void getDataBuffer(std::vector<ElemType>*& ptr, int index = 0)
  {
    checkIndex(index);
    ptr = &mDataVec[index];
  }

  /**
   * Get a pointer to the data buffer
   * Specialization for ArrayIndirection
   * \pre index must be 0, 1 or 2
   */
  void getDataBuffer(axom::Array<ElemType>*& ptr, int index = 0)
  {
    checkIndex(index);
    ptr = &mDataArr[index];
  }

  /**
   * Get a pointer to the data buffer
   * Specialization for ArrayViewIndirection
   * \pre index must be 0, 1 or 2
   */
  void getDataBuffer(axom::ArrayView<ElemType>& ptr, int index = 0)
  {
    checkIndex(index);
    ptr = axom::ArrayView<ElemType>(mDataVec[index].data(), mDataVec[index].size());
  }

private:
  /// Check that the index is valid; fail test if not
  void checkIndex(int index)
  {
    ASSERT_GE(index, 0);
    ASSERT_LT(index, 3);
  }

private:
  std::vector<ElemType> mDataVec[3];
  axom::Array<ElemType> mDataArr[3];
};

template <typename ElemType>
bool compareData(ElemType* a, ElemType* b)
{
  return a == b;
}

template <typename ElemType>
bool compareData(std::vector<ElemType>* a, std::vector<ElemType>* b)
{
  return a == b;
}

template <typename ElemType>
bool compareData(axom::ArrayView<ElemType> a, axom::ArrayView<ElemType> b)
{
  return (a.data() == b.data()) && (a.size() == b.size());
}

// Tests several types of indirection sets
using MyTypes = ::testing::Types<slam::CArrayIndirectionSet<std::int32_t, std::int64_t>,
                                 slam::VectorIndirectionSet<std::int32_t, std::int64_t>,
                                 slam::ArrayIndirectionSet<std::int32_t, std::int64_t>,
                                 slam::ArrayViewIndirectionSet<std::int32_t, std::int64_t>>;

TYPED_TEST_SUITE(IndirectionSetTester, MyTypes);

TYPED_TEST(IndirectionSetTester, constuct)
{
  SLIC_INFO("Testing simple indirection Set constructors");

  using SetType = typename TestFixture::SetType;

  // Empty set
  {
    SetType s_empty;
    EXPECT_TRUE(s_empty.empty());
    EXPECT_EQ(0, s_empty.size());
    EXPECT_TRUE(s_empty.isValid());
    EXPECT_FALSE(s_empty.hasIndirection());
  }

  // Incrementally build the set
  {
    SetType s(MAX_SET_SIZE);

    // set is originally empty
    EXPECT_FALSE(s.empty());
    EXPECT_EQ(MAX_SET_SIZE, s.size());
    EXPECT_FALSE(s.hasIndirection());
    EXPECT_FALSE(s.isValid()) << "IndirectionSet not valid until we set the indirection buffer";

    // Add data to the set
    this->getDataBuffer(s.ptr());

    EXPECT_FALSE(s.empty());
    EXPECT_EQ(MAX_SET_SIZE, s.size());
    EXPECT_TRUE(s.hasIndirection());
    EXPECT_TRUE(s.isValid());
  }
}

TYPED_TEST(IndirectionSetTester, set_builder)
{
  SLIC_INFO("Testing construction of IndirectionSet using SetBuilders");

  using SetType = typename TestFixture::SetType;
  using BufferPtrType = typename TestFixture::BufferPtrType;
  using SetBuilder = typename SetType::SetBuilder;

  BufferPtrType bufPtr {};
  this->getDataBuffer(bufPtr);

  bool bVerbose = false;

  // Only set the size, but not the data
  {
    SetType s_empty(SetBuilder().size(MAX_SET_SIZE));

    EXPECT_FALSE(s_empty.empty());
    EXPECT_EQ(MAX_SET_SIZE, s_empty.size());
    EXPECT_FALSE(s_empty.hasIndirection());
    EXPECT_FALSE(s_empty.isValid(bVerbose));
  }

  // Set the data but not the size
  {
    SetType s_empty2(SetBuilder().data(bufPtr));

    EXPECT_TRUE(s_empty2.empty());
    EXPECT_EQ(0, s_empty2.size());
    EXPECT_TRUE(s_empty2.hasIndirection());
    EXPECT_TRUE(s_empty2.isValid(bVerbose));
  }

  // Set the data and the size
  {
    SetType s(SetBuilder().size(MAX_SET_SIZE).data(bufPtr));

    EXPECT_FALSE(s.empty());
    EXPECT_EQ(MAX_SET_SIZE, s.size());
    EXPECT_EQ(0, s.offset());
    EXPECT_EQ(1, s.stride());
    EXPECT_TRUE(s.hasIndirection());
    EXPECT_TRUE(s.isValid(bVerbose));
  }

  // Set the data and the size
  // also set the offset to 0 and the stride to 1
  {
    SetType s(SetBuilder().size(MAX_SET_SIZE).offset(0).stride(1).data(bufPtr));

    EXPECT_FALSE(s.empty());
    EXPECT_EQ(MAX_SET_SIZE, s.size());
    EXPECT_EQ(0, s.offset());
    EXPECT_EQ(1, s.stride());
    EXPECT_TRUE(s.hasIndirection());
    EXPECT_TRUE(s.isValid(bVerbose));
  }
}

TYPED_TEST(IndirectionSetTester, iterate)
{
  using SetType = typename TestFixture::SetType;
  using BufferPtrType = typename TestFixture::BufferPtrType;
  using SetPosition = typename TestFixture::PosType;
  using SetElement = typename TestFixture::ElemType;

  BufferPtrType bufPtr {};
  this->getDataBuffer(bufPtr);

  SetType s(typename SetType::SetBuilder().size(MAX_SET_SIZE).data(bufPtr));

  SLIC_INFO("Iterating through set of size " << s.size());
  EXPECT_EQ(MAX_SET_SIZE, s.size());
  EXPECT_TRUE(s.isValid());

  SLIC_INFO("Using random access -- operator[] and at()");
  {
    std::stringstream sstr;
    for(SetPosition pos = SetPosition(); pos < s.size(); ++pos)
    {
      SetElement exp = static_cast<SetElement>(pos);
      EXPECT_EQ(exp, s[pos]);
      EXPECT_EQ(exp, s.at(pos));

      sstr << s[pos] << "\t";
    }
    SLIC_INFO("Data using operator[]:\t" << sstr.str());

    // using range for over OrderedSet::positions()
    int count = 0;
    for(auto pos : s.positions())
    {
      SetElement exp = static_cast<SetElement>(pos);
      EXPECT_EQ(exp, s[pos]);
      EXPECT_EQ(exp, s.at(pos));
      ++count;
    }
    EXPECT_EQ(s.size(), count);
  }

  SLIC_INFO("Using iterator interface");
  {
    // also tests default .ctor and operator=()

    std::stringstream sstr;
    int count = 0;
    typename SetType::iterator it, itEnd;
    EXPECT_EQ(it, itEnd);

    for(it = s.begin(), itEnd = s.end(); it != itEnd; ++it)
    {
      EXPECT_EQ(std::distance(s.begin(), it), *it);
      ++count;

      sstr << *it << "\t";
    }
    EXPECT_EQ(s.size(), count);
    SLIC_INFO("Data from iterator:\t" << sstr.str());
  }

  const int NUM_PERMUTATIONS = 3;
  SLIC_INFO("Printing a few permutations of the data.");
  for(int i = 0; i < NUM_PERMUTATIONS; ++i)
  {
    this->permuteData();
    std::stringstream sstr;
    std::copy(s.begin(), s.end(), std::ostream_iterator<SetPosition>(sstr, "\t"));
    SLIC_INFO("Permutation " << i << ":\t" << sstr.str());
  }
}

TYPED_TEST(IndirectionSetTester, equality)
{
  SLIC_INFO("Testing equality/inequality for several sets");

  using SetType = typename TestFixture::SetType;
  using BufferPtrType = typename TestFixture::BufferPtrType;
  using SetPosition = typename TestFixture::PosType;
  using SetElement = typename TestFixture::ElemType;

  // Get pointers to the data
  BufferPtrType bufPtr[3] = {{}, {}, {}};
  this->getDataBuffer(bufPtr[0], 0);
  this->getDataBuffer(bufPtr[1], 1);
  this->getDataBuffer(bufPtr[2], 2);
  EXPECT_FALSE(compareData(bufPtr[0], bufPtr[1]));
  EXPECT_FALSE(compareData(bufPtr[0], bufPtr[2]));
  EXPECT_FALSE(compareData(bufPtr[1], bufPtr[2]));

  // Initialize the first set
  SetType s0(typename SetType::SetBuilder()  //
               .size(MAX_SET_SIZE)           //
               .data(bufPtr[0]));

  EXPECT_TRUE(s0.isValid());

  // Test self-equality
  EXPECT_EQ(s0, s0);

  // Test equality against a set with the same array
  SetType s0_b(typename SetType::SetBuilder()  //
                 .size(MAX_SET_SIZE)           //
                 .data(s0.ptr()));
  EXPECT_TRUE(s0_b.isValid());
  EXPECT_EQ(s0, s0_b);
  EXPECT_EQ(s0.ptr(), s0_b.ptr());

  // Test equality against a set with a different array (same values)

  SetType s1(typename SetType::SetBuilder()  //
               .size(MAX_SET_SIZE)           //
               .data(bufPtr[1]));
  EXPECT_TRUE(s1.isValid());
  EXPECT_EQ(s0, s1);
  EXPECT_EQ(s1, s0);

  // Check that sets are not equal if we permute the data
  this->permuteData(1);
  EXPECT_NE(s0, s1);

  // Test against a smaller set
  SetType s2(typename SetType::SetBuilder()  //
               .size(MAX_SET_SIZE / 2)       //
               .data(bufPtr[2]));
  EXPECT_TRUE(s2.isValid());
  EXPECT_NE(s0, s2);
  EXPECT_NE(s1, s2);

  // Check against a PositionSet with the same elements
  slam::PositionSet<SetPosition, SetElement> posSet(MAX_SET_SIZE);
  EXPECT_TRUE(posSet.isValid());
  EXPECT_EQ(posSet, s0);
}

TYPED_TEST(IndirectionSetTester, out_of_bounds)
{
  SLIC_INFO("Testing out of bounds access on initialized set"
            << "-- code is expected to assert and die.");

  using SetType = typename TestFixture::SetType;
  using BufferPtrType = typename TestFixture::BufferPtrType;

  BufferPtrType bufPtr {};
  this->getDataBuffer(bufPtr);

  SetType s(typename SetType::SetBuilder()  //
              .size(MAX_SET_SIZE)           //
              .data(bufPtr));

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode,
  // so this test will only fail in debug mode
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED(s[MAX_SET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif
}

TEST(slam_set_indirectionset, compare_array_and_vector)
{
  SLIC_INFO("Comparing Array and Vector indirection Sets");

  using P = int;
  using E = int;

  using VecSet = slam::VectorIndirectionSet<P, E>;
  using CArrSet = slam::CArrayIndirectionSet<P, E>;
  using PosSet = slam::PositionSet<P, E>;

  const int SZ = MAX_SET_SIZE;
  std::vector<E> vVals(SZ);
  std::vector<E> aVals(SZ);

  PosSet ps(SZ);

  // Create vector and array sets without initializing their elements
  VecSet vs(VecSet::SetBuilder()  //
              .size(SZ)           //
              .data(&vVals));

  CArrSet as(CArrSet::SetBuilder()  //
               .size(SZ)            //
               .data(aVals.data()));

  EXPECT_NE(vs, ps);
  EXPECT_NE(as, ps);

  // Use position set to initialize vector and array set elements
  for(auto i : ps)
  {
    vs[i] = i;
    as[i] = i;
  }

  EXPECT_EQ(vs, as);
  EXPECT_EQ(as, ps);
  EXPECT_EQ(ps, vs);
}

TEST(slam_set_indirectionset, negative_stride)
{
  SLIC_INFO("Testing negative strides on an indirection set");

  namespace policies = axom::slam::policies;

  using SetPosition = slam::DefaultPositionType;
  using SetElement = slam::DefaultElementType;

  using SizePol = policies::RuntimeSize<SetPosition>;
  ;
  using OffPol = policies::RuntimeOffset<SetPosition>;
  using StridePol = policies::RuntimeStride<SetPosition>;
  using VecIndPol = policies::STLVectorIndirection<SetPosition, SetElement>;
  using CArrIndPol = policies::CArrayIndirection<SetPosition, SetElement>;

  using VecSet = slam::OrderedSet<SetPosition, SetElement, SizePol, OffPol, StridePol, VecIndPol>;

  using CArrSet = slam::OrderedSet<SetPosition, SetElement, SizePol, OffPol, StridePol, CArrIndPol>;

  // Set up data -- an array of incrementing integers
  std::vector<SetElement> intVec(MAX_SET_SIZE);
  for(auto i : slam::PositionSet<>(MAX_SET_SIZE))
  {
    intVec[i] = i;
  }

  const int setSize = MAX_SET_SIZE / 2;
  const int setOffset = MAX_SET_SIZE - 1;
  const int setStride = -2;
  const bool bVerbose = false;

  // Setup the VectorIndirectionSet and test basic functionality
  VecSet vSet(VecSet::SetBuilder()  //
                .size(setSize)      //
                .offset(setOffset)  //
                .stride(setStride)  //
                .data(&intVec));
  {
    EXPECT_TRUE(vSet.isValid());
    EXPECT_FALSE(vSet.empty());
    EXPECT_EQ(setSize, vSet.size());
    EXPECT_TRUE(vSet.hasIndirection());

    SLIC_INFO("Ordered vector set has:" << "{ size: " << vSet.size() << ", stride: " << vSet.stride()
                                        << ", offset: " << vSet.offset() << ", first elt: " << vSet[0]
                                        << ", last elt: " << vSet[vSet.size() - 1] << "}");

    // Test the elements
    EXPECT_EQ(intVec[setOffset], vSet[0]);
    for(int i = 0; i < vSet.size(); ++i)
    {
      EXPECT_EQ(setOffset + setStride * i, vSet[i]);
    }

    /// Several checks that sets with bad offsets and strides are invalid
    SLIC_DEBUG_IF(
      bVerbose,
      "--- Checking isValid() on several sets with " << "bad sizes, offsets and strides.");
    VecSet noDataVSet(VecSet::SetBuilder()  // Note: Missing a data pointer
                        .size(setSize)
                        .offset(setOffset)
                        .stride(setStride));
    EXPECT_FALSE(noDataVSet.isValid(bVerbose));

    VecSet outOfBoundsVSet(VecSet::SetBuilder()
                             .size(setSize + 1)  // Note: last index is out of bounds
                             .offset(setOffset)
                             .stride(setStride)
                             .data(&intVec));
    EXPECT_FALSE(outOfBoundsVSet.isValid(bVerbose));

    VecSet outOfBoundsVSet2(VecSet::SetBuilder()
                              .size(setSize)
                              .offset(-1)  // Note: first index is out of bounds
                              .stride(1)
                              .data(&intVec));
    EXPECT_FALSE(outOfBoundsVSet2.isValid(bVerbose));

    VecSet zeroStrideVSet(VecSet::SetBuilder()
                            .size(setSize)
                            .offset(setOffset)
                            .stride(0)  // Note: A stride of zero is not valid
                            .data(&intVec));
    EXPECT_FALSE(zeroStrideVSet.isValid(bVerbose));

    SLIC_DEBUG_IF(bVerbose, "--- Done.");
  }

  // Setup the CArrayIndirectionSet and test basic functionality
  CArrSet cSet(CArrSet::SetBuilder()  //
                 .size(setSize)       //
                 .offset(setOffset)   //
                 .stride(setStride)   //
                 .data(intVec.data()));
  {
    EXPECT_TRUE(cSet.isValid());
    EXPECT_FALSE(cSet.empty());
    EXPECT_EQ(setSize, cSet.size());
    EXPECT_TRUE(cSet.hasIndirection());

    EXPECT_EQ(intVec[setOffset], cSet[0]);

    SLIC_INFO("Ordered array set has:" << "{ size: " << cSet.size() << ", stride: " << cSet.stride()
                                       << ", offset: " << cSet.offset() << ", first elt: " << cSet[0]
                                       << ", last elt: " << cSet[cSet.size() - 1] << "}");

    for(int i = 0; i < cSet.size(); ++i)
    {
      EXPECT_EQ(setOffset + setStride * i, cSet[i]);
    }

    /// Several checks that sets with bad offsets and strides are invalid
    SLIC_DEBUG_IF(
      bVerbose,
      "--- Checking isValid() on several sets with " << "bad sizes, offsets and strides.");

    CArrSet noDataCSet(CArrSet::SetBuilder()  // Note: Missing a data pointer
                         .size(setSize)
                         .offset(setOffset)
                         .stride(setStride));
    EXPECT_FALSE(noDataCSet.isValid(bVerbose));

    CArrSet outOfBoundsCSet1(CArrSet::SetBuilder()
                               .size(setSize + 1)  // Note: last index is out of bounds
                               .offset(setOffset)
                               .stride(setStride)
                               .data(intVec.data()));
    EXPECT_FALSE(outOfBoundsCSet1.isValid(bVerbose));

    CArrSet outOfBoundsCSet2(CArrSet::SetBuilder()
                               .size(setSize)
                               .offset(-1)  // Note: first index is out of bounds
                               .stride(1)
                               .data(intVec.data()));
    EXPECT_FALSE(outOfBoundsCSet2.isValid(bVerbose));

    CArrSet zeroStrideCSet(CArrSet::SetBuilder()
                             .size(setSize)
                             .offset(setOffset)
                             .stride(0)  // Note: stride of zero is not valid
                             .data(intVec.data()));
    EXPECT_FALSE(zeroStrideCSet.isValid(bVerbose));

    SLIC_DEBUG_IF(bVerbose, "--- Done.");
  }

  // check that vset and cset are equivalent
  EXPECT_EQ(vSet, cSet);
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  int result = RUN_ALL_TESTS();

  return result;
}
