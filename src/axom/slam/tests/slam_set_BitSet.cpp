// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 *  This file tests the BitSet class of slam
 *  Most tests are parameterized by the values in the testSizes() function
 */

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/slam/BitSet.hpp"

#include "gtest/gtest.h"

namespace slam = axom::slam;

namespace
{
// Utility function to initialize a bitset of size \a size
// where every stride^th bit is set, starting with offset \a offset
// Note: Use only after establishing that the constructor and set(idx)
//       functions work properly
slam::BitSet generateBitset(int size, int stride = 1, int offset = 0)
{
  using Index = slam::BitSet::Index;

  slam::BitSet bitset(size);
  for(Index i = offset; i < size; i += stride)
  {
    bitset.set(i);
  }

  return bitset;
}

// The test sizes for bitsets
std::vector<int> testSizes()
{
  std::vector<int> vals;

  vals.push_back(0);    // empty bitset
  vals.push_back(23);   // less than one word
  vals.push_back(31);   // one bit less than a 4-byte word
  vals.push_back(32);   // exactly one 4-byte word
  vals.push_back(33);   // one bit more than a 4-byte word
  vals.push_back(63);   // one bit less than an 8-byte word
  vals.push_back(64);   // exactly one 8-byte word
  vals.push_back(65);   // one bit more than an 8-byte word
  vals.push_back(127);  // one bit less than two words
  vals.push_back(128);  // two words
  vals.push_back(129);  // more than two words
  vals.push_back(153);
  vals.push_back(1547);  // large bitset
  //vals.push_back(1234567);  // large bitset

  return vals;
}
}  // namespace

class SlamBitSet : public ::testing::TestWithParam<int>
{ };

TEST_P(SlamBitSet, checkInitEmpty)
{
  const int NBITS = GetParam();
  SLIC_INFO("Testing bitset construction (" << NBITS << " bits)");

  slam::BitSet bitset(NBITS);

  EXPECT_TRUE(bitset.isValid());
  EXPECT_EQ(NBITS, bitset.size());
  EXPECT_EQ(0, bitset.count());

  // Test that each bit is off
  for(int i = 0; i < NBITS; ++i)
  {
    EXPECT_FALSE(bitset.test(i));
  }

  // It is always ok to test bits outside the size
  {
    EXPECT_FALSE(bitset.test(-1));
    EXPECT_FALSE(bitset.test(bitset.size() + 1));
  }
}

TEST_P(SlamBitSet, setClearFlipAllBits)
{
  const int NBITS = GetParam();
  SLIC_INFO("Testing bitset set, clear and flip (" << NBITS << " bits)");

  slam::BitSet bitset(NBITS);

  // count should be 0
  EXPECT_EQ(0, bitset.count());
  EXPECT_TRUE(bitset.isValid());

  // clear all bits, when already empty
  bitset.clear();
  EXPECT_EQ(0, bitset.count());
  EXPECT_TRUE(bitset.isValid());

  // set all bits
  bitset.set();
  EXPECT_EQ(NBITS, bitset.count());
  EXPECT_TRUE(bitset.isValid());

  // clear all bits
  bitset.clear();
  EXPECT_EQ(0, bitset.count());
  EXPECT_TRUE(bitset.isValid());

  // toggle all bits
  bitset.flip();
  EXPECT_EQ(NBITS, bitset.count());
  EXPECT_TRUE(bitset.isValid());
}

TEST_P(SlamBitSet, setAndTestIndividualBits)
{
  const int NBITS = GetParam();
  SLIC_INFO("Testing bitset set and test individual bits (" << NBITS << " bits)");

  const int STRIDE = 5;

  slam::BitSet bitset(NBITS);

  int nSetBits = 0;
  for(int i = 0; i < NBITS; i += STRIDE)
  {
    bitset.set(i);
    ++nSetBits;
  }
  EXPECT_TRUE(bitset.isValid());

  for(int i = 0; i < NBITS; ++i)
  {
    const bool exp = (i % STRIDE == 0) ? true : false;
    EXPECT_EQ(exp, bitset.test(i));
  }

  // count should be equal to number of manually set bits
  EXPECT_EQ(nSetBits, bitset.count());
}

TEST_P(SlamBitSet, flipAndTestIndividualBits)
{
  int const NBITS = GetParam();
  SLIC_INFO("Testing bitset flip for individual bits (" << NBITS << " bits)");

  slam::BitSet bitset(NBITS);

  // Set fizzbuzz bits
  const int STRIDE1 = 3;
  for(int i = 0; i < NBITS; i += STRIDE1)
  {
    bitset.flip(i);
  }

  const int STRIDE2 = 5;
  for(int i = 0; i < NBITS; i += STRIDE2)
  {
    bitset.flip(i);
  }

  // Check single-bit flip()
  for(int i = 0; i < NBITS; ++i)
  {
    const bool exp = (i % STRIDE1 == 0) ^ (i % STRIDE2 == 0) ? true : false;
    EXPECT_EQ(exp, bitset.test(i));
  }
  EXPECT_TRUE(bitset.isValid());

  // Flip all bits and recheck all-bits flip()
  bitset.flip();
  for(int i = 0; i < NBITS; ++i)
  {
    const bool exp = (i % STRIDE1 == 0) ^ (i % STRIDE2 == 0) ? false  // note: oppposite of last test
                                                             : true;
    EXPECT_EQ(exp, bitset.test(i));
  }
  EXPECT_TRUE(bitset.isValid());
}

TEST_P(SlamBitSet, copyAssign)
{
  using Index = slam::BitSet::Index;
  const int NBITS = GetParam();
  const Index STRIDE = 5;
  const Index OFFSET = 2;

  SLIC_INFO("Testing bitset copy and assignment (" << NBITS << " bits)");

  // Construct the first bitset and check validity
  slam::BitSet bitset1 = generateBitset(NBITS, STRIDE, OFFSET);
  EXPECT_TRUE(bitset1.isValid());
  EXPECT_EQ(NBITS, bitset1.size());

  // Test copy constructor and check validity
  slam::BitSet bitset2(bitset1);
  EXPECT_TRUE(bitset2.isValid());

  EXPECT_EQ(bitset1.size(), bitset2.size());
  EXPECT_EQ(bitset1.count(), bitset2.count());
  EXPECT_EQ(bitset1, bitset2);

  // Test the assignment operator and check validity
  slam::BitSet bitset3(NBITS);
  EXPECT_EQ(bitset1.size(), bitset3.size());

  bitset3 = bitset1;
  EXPECT_EQ(bitset1.size(), bitset3.size());
  EXPECT_EQ(bitset1.count(), bitset3.count());
  EXPECT_EQ(bitset1, bitset3);
}

TEST_P(SlamBitSet, iterator)
{
  const int NBITS = GetParam();
  SLIC_INFO("Testing bitset iterator interface with " << NBITS << " bits");

  using Index = slam::BitSet::Index;
  const Index STRIDE = 5;
  const Index OFFSET = 2;

  const Index npos = slam::BitSet::npos;
  slam::BitSet bitset = generateBitset(NBITS, STRIDE, OFFSET);

  if(NBITS < OFFSET)
  {
    EXPECT_EQ(npos, bitset.find_first());
  }
  else
  {
    Index startIdx = bitset.find_first();
    EXPECT_EQ(OFFSET, startIdx);

    int numFound = 0;
    for(Index idx = startIdx; idx != slam::BitSet::npos; idx = bitset.find_next(idx))
    {
      EXPECT_EQ(OFFSET, idx % STRIDE);
      ++numFound;
    }

    EXPECT_EQ(bitset.count(), numFound);
  }
}

TEST_P(SlamBitSet, atomicOps)
{
  using Index = slam::BitSet::Index;
  const int NBITS = GetParam();
  const Index STRIDE = 5;
  const Index OFFSET = 2;

  SLIC_INFO("Testing atomic ops on bitset (" << NBITS << " bits)");

  // Construct the first bitset and check validity
  slam::BitSet bitset1 = generateBitset(NBITS, STRIDE, OFFSET);
  EXPECT_TRUE(bitset1.isValid());
  EXPECT_EQ(NBITS, bitset1.size());

  // Construct the second bitset from the first
  slam::BitSet bitset2(NBITS);
  for(Index i = 0; i < NBITS; ++i)
  {
    if(bitset1.test(i))
    {
      bitset2.atomicSet(i);
    }
    else
    {
      bitset2.atomicClear(i);
    }
  }
  EXPECT_TRUE(bitset2.isValid());
  EXPECT_EQ(NBITS, bitset2.size());

  // Check for equality
  EXPECT_EQ(bitset1.count(), bitset2.count());
  EXPECT_EQ(bitset1, bitset2);

  // Apply several atomic operations w/ checks
  for(int i = 0; i < NBITS; ++i)
  {
    EXPECT_EQ(bitset1.test(i), bitset2.test(i));

    bitset2.atomicFlip(i);
    EXPECT_NE(bitset1.test(i), bitset2.test(i));

    bitset1.flip(i);
    EXPECT_EQ(bitset1.test(i), bitset2.test(i));

    bitset2.atomicClear(i);
    EXPECT_FALSE(bitset2.test(i));

    bitset2.atomicSet(i);
    EXPECT_TRUE(bitset2.test(i));

    bitset2.atomicClear(i);
    EXPECT_FALSE(bitset2.test(i));
  }

  bitset1.clear();
  EXPECT_EQ(bitset1.count(), bitset2.count());
  EXPECT_EQ(bitset1, bitset2);
}

TEST_P(SlamBitSet, unionOperator)
{
  const int NBITS = GetParam();
  SLIC_INFO("Testing bitset union operator (" << NBITS << " bits");

  using Index = slam::BitSet::Index;
  const Index STRIDE1 = 3;
  const Index STRIDE2 = 5;

  slam::BitSet bitset1 = generateBitset(NBITS, STRIDE1);
  slam::BitSet bitset2 = generateBitset(NBITS, STRIDE2);

  // Apply and test union operator
  slam::BitSet bitset3 = bitset1 | bitset2;

  for(int i = 0; i < NBITS; ++i)
  {
    const bool exp = (i % STRIDE1 == 0) || (i % STRIDE2 == 0) ? true : false;
    EXPECT_EQ(exp, bitset3.test(i));
  }
  EXPECT_TRUE(bitset3.isValid());

  // Test union-assignment operator
  bitset1 |= bitset2;
  EXPECT_EQ(bitset1, bitset3);
  EXPECT_TRUE(bitset1.isValid());
}

TEST_P(SlamBitSet, intersectOperator)
{
  const int NBITS = GetParam();
  SLIC_INFO("Testing bitset intersect operator (" << NBITS << " bits");

  using Index = slam::BitSet::Index;
  const Index STRIDE1 = 3;
  const Index STRIDE2 = 5;

  slam::BitSet bitset1 = generateBitset(NBITS, STRIDE1);
  slam::BitSet bitset2 = generateBitset(NBITS, STRIDE2);

  // Apply and test intersection operator
  slam::BitSet bitset3 = bitset1 & bitset2;

  for(int i = 0; i < NBITS; ++i)
  {
    const bool exp = (i % STRIDE1 == 0) && (i % STRIDE2 == 0) ? true : false;
    EXPECT_EQ(exp, bitset3.test(i));
  }
  EXPECT_TRUE(bitset3.isValid());

  // Test intersection-assignment operator
  bitset1 &= bitset2;
  EXPECT_EQ(bitset1, bitset3);
  EXPECT_TRUE(bitset1.isValid());
}

TEST_P(SlamBitSet, xorOperator)
{
  const int NBITS = GetParam();
  SLIC_INFO("Testing bitset xor operator (" << NBITS << " bits");

  using Index = slam::BitSet::Index;
  const Index STRIDE1 = 3;
  const Index STRIDE2 = 5;

  slam::BitSet bitset1 = generateBitset(NBITS, STRIDE1);
  slam::BitSet bitset2 = generateBitset(NBITS, STRIDE2);

  // Apply and test xor operator
  slam::BitSet bitset3 = bitset1 ^ bitset2;

  for(int i = 0; i < NBITS; ++i)
  {
    const bool exp = (i % STRIDE1 == 0) ^ (i % STRIDE2 == 0) ? true : false;
    EXPECT_EQ(exp, bitset3.test(i));
  }
  EXPECT_TRUE(bitset3.isValid());

  // Test xor-assignment operator
  bitset1 ^= bitset2;
  EXPECT_EQ(bitset1, bitset3);
  EXPECT_TRUE(bitset1.isValid());
}

TEST_P(SlamBitSet, differenceOperator)
{
  const int NBITS = GetParam();
  SLIC_INFO("Testing bitset difference operator (" << NBITS << " bits");

  using Index = slam::BitSet::Index;
  const Index STRIDE1 = 3;
  const Index STRIDE2 = 5;

  slam::BitSet bitset1 = generateBitset(NBITS, STRIDE1);
  slam::BitSet bitset2 = generateBitset(NBITS, STRIDE2);

  // Apply and test difference operator
  slam::BitSet bitset3 = bitset1 - bitset2;

  for(int i = 0; i < NBITS; ++i)
  {
    const bool exp = (i % STRIDE1 == 0) && !(i % STRIDE2 == 0) ? true : false;
    EXPECT_EQ(exp, bitset3.test(i));
  }
  EXPECT_TRUE(bitset3.isValid());

  // Test difference-assignment operator
  bitset1 -= bitset2;
  EXPECT_EQ(bitset1, bitset3);
  EXPECT_TRUE(bitset1.isValid());
}

INSTANTIATE_TEST_SUITE_P(SlamBitSetParam, SlamBitSet, ::testing::ValuesIn(testSizes()));

TEST(slam_set_bitset, settingOutOfRange)
{
  SLIC_INFO("Testing assert/invalid when setting out-of-range bit");

  // set a bit on a zero-sized bitset
  {
    slam::BitSet bitset(0);

#ifdef AXOM_DEBUG
    SLIC_INFO("** Expecting error ** ");
    EXPECT_DEATH_IF_SUPPORTED(bitset.set(0), "");
#else
    bitset.set(0);
    EXPECT_FALSE(bitset.isValid());
#endif
  }

  // set a bit on a bitset that is less than one word long
  {
    int const NBITS = 20;
    slam::BitSet bitset(NBITS);

#ifdef AXOM_DEBUG
    SLIC_INFO("** Expecting error ** ");
    EXPECT_DEATH_IF_SUPPORTED(bitset.set(NBITS), "");
#else
    bitset.set(NBITS);
    EXPECT_FALSE(bitset.isValid());
#endif
  }

  // set a bit on a bitset that is more than one word long
  {
    int const NBITS = 120;
    slam::BitSet bitset(NBITS);

#ifdef AXOM_DEBUG
    SLIC_INFO("** Expecting error ** ");
    EXPECT_DEATH_IF_SUPPORTED(bitset.set(200), "");
#else
    bitset.set(NBITS);
    EXPECT_FALSE(bitset.isValid());

    // Note: The following will cause the program to access memory
    // that is out-of-bounds, and will cause an assert in debug runs
    // but will not fail the isValid() test
    //bitset.set(200);
#endif
  }
}

TEST(slam_set_bitset, moreIterators)
{
  SLIC_INFO("More testing of BitSet iteration interface (first_bit(), next_bit())");

  using Index = slam::BitSet::Index;

  const Index npos = slam::BitSet::npos;

  // Test empty bitset
  {
    slam::BitSet bitset(0);

    EXPECT_EQ(npos, bitset.find_first());
    EXPECT_EQ(npos, bitset.find_next(npos));
  }

  // Test non-empty bitset with one word
  {
    const int NBITS = 21;
    slam::BitSet bitset(NBITS);
    bitset.set(10);

    EXPECT_EQ(10, bitset.find_first());
    EXPECT_EQ(10, bitset.find_next(0));
    EXPECT_EQ(10, bitset.find_next(1));
    EXPECT_EQ(10, bitset.find_next(9));
    EXPECT_EQ(npos, bitset.find_next(10));

    bitset.set(15);
    EXPECT_EQ(15, bitset.find_next(10));

    // Set last bit
    bitset.set(20);
    EXPECT_EQ(20, bitset.find_next(15));
    EXPECT_EQ(npos, bitset.find_next(20));

    // Find next bit after npos
    EXPECT_EQ(npos, bitset.find_next(npos));
  }
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

#if AXOM_DEBUG
  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
#endif

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
