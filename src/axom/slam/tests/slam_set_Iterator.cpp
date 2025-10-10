// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_set_Iterator.cpp
 *
 * \brief Unit tests for Slam's OrderedSet iterator class
 *
 * Most of the tests in this file are templated by Set type
 * on the OrderedSetIteratorTester<SetType> test class.
 *
 * Note: To test a new OrderedSet class, you need to
 * a) Add a specialization of generateSet()
 * b) Add your type to the MyTypes alias below
 */

#include <iterator>
#include <vector>

#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/slam/OrderedSet.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/slam/IndirectionSet.hpp"

namespace
{
namespace slam = axom::slam;

using SetPosition = slam::DefaultPositionType;
using SetElement = slam::DefaultElementType;

using SetBase = slam::Set<SetPosition, SetElement>;

using PositionSet = slam::PositionSet<SetPosition, SetElement>;
using RangeSet = slam::RangeSet<SetPosition, SetElement>;
using VectorSet = slam::VectorIndirectionSet<SetPosition, SetElement>;
using CArraySet = slam::CArrayIndirectionSet<SetPosition, SetElement>;

constexpr int SET_SIZE = 10;

/// Utility function to initialize a set for the test data
/// Specialized for each set type
template <typename S, typename E>
S generateSet(std::vector<E>& vec, int size);

/// Specialization of \a generateSet for \a PositionSet
template <>
PositionSet generateSet(std::vector<SetElement>&, int size)
{
  return PositionSet(size);
}

/// Specialization of \a generateSet for \a RangeSet
template <>
RangeSet generateSet(std::vector<SetElement>&, int size)
{
  auto lower = 10;
  auto upper = lower + size;

  return RangeSet(lower, upper);
}

/// Specialization of \a generateSet for \a VectorSet
template <>
VectorSet generateSet(std::vector<SetElement>& vec, int size)
{
  return VectorSet(VectorSet::SetBuilder().size(size).data(&vec));
}

/// Specialization of \a generateSet for \a CArraySet
template <>
CArraySet generateSet(std::vector<SetElement>& vec, int size)
{
  return CArraySet(CArraySet::SetBuilder().size(size).data(vec.data()));
}

}  // end anonymous namespace

/**
 * Simple test class for iterators on ordered sets
 * Templated on a Set type (see MyTypes alias below)
 */
template <typename TheSet>
class OrderedSetIteratorTester : public ::testing::Test
{
public:
  using SetType = TheSet;

  using PosType = typename SetType::PositionType;
  using ElemType = typename SetType::ElementType;

  virtual void SetUp()
  {
    // initialize data for indirection-based sets
    // The values don't really matter for this test.
    mVec = std::vector<ElemType>(SET_SIZE);

    for(int i = 0; i < SET_SIZE; ++i)
    {
      mVec[i] = (SET_SIZE) / 2 - i;
    }

    // generate the set using specialized functions (defined above)
    mSet = generateSet<SetType>(mVec, SET_SIZE);
  }

  SetType& getSet() { return mSet; }

private:
  SetType mSet;
  std::vector<ElemType> mVec;
};

// Tests several types of sets
using MyTypes = ::testing::Types<PositionSet, RangeSet, VectorSet, CArraySet>;
TYPED_TEST_SUITE(OrderedSetIteratorTester, MyTypes);

TYPED_TEST(OrderedSetIteratorTester, basic_operations)
{
  using SetType = typename TestFixture::SetType;
  using iterator = typename SetType::iterator;

  SetType set = this->getSet();

  // equal to self
  {
    iterator it = set.begin();
    EXPECT_EQ(it, set.begin());
    EXPECT_NE(it, set.end());
    EXPECT_EQ(it, it);
  }

  // simple equality
  {
    iterator it1 = set.begin();
    iterator it2 = set.begin();
    EXPECT_EQ(it1, it2);
    EXPECT_TRUE(it1 == it2);

    EXPECT_EQ(*it1, *it2);
  }

  // test assignment
  {
    iterator it1 = set.begin();
    iterator it2 = set.end();
    EXPECT_NE(it1, it2);

    it2 = it1;
    EXPECT_EQ(it1, it2);
    EXPECT_EQ(*it1, *it2);
  }

  // test increment/decrement
  {
    iterator it1 = set.begin() + set.size() / 2;
    iterator it2 = it1;

    EXPECT_EQ(it1, it2);

    // post-increment
    EXPECT_EQ(it1, it2++);
    EXPECT_NE(it1, it2);
    EXPECT_EQ(it1 + 1, it2);

    // post-decrement
    it2 = it1;
    EXPECT_EQ(it1, it2--);
    EXPECT_NE(it1, it2);
    EXPECT_EQ(it1 - 1, it2);

    // pre-increment
    it2 = it1;
    EXPECT_NE(it1, ++it2);
    EXPECT_EQ(it1 + 1, it2);

    // pre-decrement
    it2 = it1;
    EXPECT_NE(it1, --it2);
    EXPECT_EQ(it1 - 1, it2);

    // pre- and post- differ
    it2 = it1;
    EXPECT_NE(it1++, ++it2);
    EXPECT_NE(it1--, --it2);
  }

  // test relational ops
  {
    iterator it1 = set.begin();
    iterator it2 = it1 + 1;

    EXPECT_NE(it2, it1);

    EXPECT_GT(it2, it1);
    EXPECT_GE(it2, it1);

    EXPECT_LT(it1, it2);
    EXPECT_LE(it1, it2);
  }

  // test arithmetic ops and random access
  {
    iterator it1 = set.begin();
    iterator it2 = it1 + 2;

    EXPECT_EQ(*it2, it2[0]);

    EXPECT_NE(it1, it2);
    EXPECT_EQ(it1 + 2, it2);
    EXPECT_EQ(it1, it2 - 2);

    EXPECT_EQ(it1[2], *it2);
    EXPECT_EQ(*it1, it2[-2]);
    EXPECT_EQ(it1[1], it2[-1]);
  }
}

TYPED_TEST(OrderedSetIteratorTester, const_set)
{
  using SetType = typename TestFixture::SetType;

  const SetType set = this->getSet();

  // cannot create (non-const) iterator from const Set
  {
    // Note: Uncommenting the following line generates a compiler error
    //   error: no viable conversion from 'OrderedSetIterator<[...], true>'
    //                                 to 'OrderedSetIterator<[...], false>
    //  (where the bool param indicates whether the iterator is const or not)

    // typename SetType::iterator invalid = set.begin();
  }

  // check equality of const iterators
  {
    typename SetType::const_iterator b1 = set.begin();
    typename SetType::const_iterator b2 = set.begin();
    typename SetType::const_iterator e1 = set.end();

    EXPECT_EQ(b1, b2);
    EXPECT_EQ(b2, b1);

    EXPECT_NE(b1, e1);
    EXPECT_NE(e1, b1);
  }

  // apply some simple operators
  {
    auto b1 = set.begin();
    auto b2 = set.begin();
    auto e1 = set.end();

    EXPECT_EQ(b1, b2);
    EXPECT_EQ(b2, b1);

    EXPECT_NE(b1, e1);
    EXPECT_NE(e1, b1);

    EXPECT_EQ(*b1, *b2);
    EXPECT_EQ(b1[1], b2[1]);

    EXPECT_NE(b1 + 1, b2);

    b2++;
    EXPECT_EQ(b1 + 1, b2);
  }
}

TYPED_TEST(OrderedSetIteratorTester, equality_const_and_non_const)
{
  using SetType = typename TestFixture::SetType;
  using iterator = typename SetType::iterator;
  using const_iterator = typename SetType::const_iterator;

  SetType set = this->getSet();

  // check equality of non-const iterators
  {
    iterator b1 = set.begin();
    iterator b2 = set.begin();
    iterator e1 = set.end();

    EXPECT_EQ(b1, b2);
    EXPECT_EQ(b2, b1);

    EXPECT_NE(b1, e1);
    EXPECT_NE(e1, b1);
  }

  // check equality of const iterators
  {
    const_iterator b1 = set.begin();
    const_iterator b2 = set.begin();
    const_iterator e1 = set.end();

    EXPECT_EQ(b1, b2);
    EXPECT_EQ(b2, b1);

    EXPECT_NE(b1, e1);
    EXPECT_NE(e1, b1);
  }

  // mix const and non-const iterators to non-const set
  // TODO: Re-enable this
  {
    iterator b1 = set.begin();
    iterator e1 = set.end();

    const_iterator b2 = set.begin();
    const_iterator b3 = set.begin();

    //EXPECT_EQ(b1, b2);
    //EXPECT_EQ(b2, b1);
    EXPECT_EQ(b2, b3);

    EXPECT_NE(b1, e1);
    EXPECT_NE(e1, b1);

    //EXPECT_NE(b2, e1);
    //EXPECT_NE(e1, b2);
  }
}

TEST(slam_set_iterator, modify_with_iterators)
{
  using SetType = VectorSet;
  using PositionType = typename SetType::PositionType;
  using ElementType = typename SetType::ElementType;
  using const_iterator = SetType::const_iterator;

  int vecSize = 10;
  VectorSet::IndirectionBufferType vec1(vecSize);
  VectorSet::IndirectionBufferType vec2(vecSize);
  VectorSet::IndirectionBufferType vec3(vecSize);

  int setSize = 5;
  SetType set1 = generateSet<SetType>(vec1, setSize);
  SetType set2 = generateSet<SetType>(vec2, setSize);
  SetType set3 = generateSet<SetType>(vec3, setSize);

  /// Modify the set values through their iterators

  // set1's values will be incrementing integers
  {
    int counter = 0;
    for(auto& el : set1)
    {
      el = counter++;  // set values using element reference
    }
  }

  // set2's values will be incrementing ints
  // but will use the previous element's value, accessed by iterator
  {
    for(auto it = set2.begin(); it != set2.end(); ++it)
    {
      *it = (it.index() == 0) ? 0 : it[-1] + 1;
    }
  }

  // set3's values will use a simple formula:
  //   set2[i] = 2 N - i,  where N is the set size
  {
    const auto sz = set3.size();
    for(auto it = set3.begin(); it != set3.end(); ++it)
    {
      // Get the element index using the iterator distance function
      PositionType pos = it - set3.begin();

      *it = (2 * sz) - pos;  // set value using iterator dereference
    }
  }

  /// Print the set values
  {
    SLIC_INFO("set 1 := {0,1,2,..., n-1}," << " where n is the set's size (" << set1.size() << ")");
    for(auto& el : set1)
    {
      SLIC_INFO(el);
    }

    SLIC_INFO("set 2 := {0,1,2,..., n-1}," << " where n is the set's size (" << set2.size() << ")");
    for(auto& el : set2)
    {
      SLIC_INFO(el);
    }

    SLIC_INFO("set 3 := { 2n -i }, " << " where n is the set's size (" << set3.size() << ")"
                                     << " and i is the element index");

    for(auto& el : set3)
    {
      SLIC_INFO(el);
    }
  }

  /// Check the values with const iterators
  {
    const_iterator it1 = set1.begin();
    const_iterator it2 = set2.begin();
    for(auto p : set1.positions())
    {
      EXPECT_EQ(p, it1[p]);
      EXPECT_EQ(*it1, *it2);
    }
  }

  /// Check the values with const iterators
  /// TODO: use explicit const_iterator
  {
    const auto sz = set3.size();

    for(auto it = set3.begin(); it != set3.end(); ++it)
    {
      auto pos = it.index();
      ElementType exp = 2 * sz - pos;

      EXPECT_EQ(exp, *it);
    }
  }
}

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
