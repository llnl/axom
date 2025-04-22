// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/NumericLimits.hpp"

#include "axom/spin/MortonIndex.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/slic.hpp"

#include <cstdlib>

// Uncomment the line below for true randomized points
#ifndef MORTON_TESTER_SHOULD_SEED
//  #define MORTON_TESTER_SHOULD_SEED
#endif

#ifdef MORTON_TESTER_SHOULD_SEED
  #include <ctime>  // for time() used by srand()
#endif

using axom::primal::Point;

namespace
{
constexpr int MAX_ITER = 10000;

// Generate a random integer in the range [beg, end)
template <typename CoordType>
CoordType randomInt(CoordType beg, CoordType end)
{
  CoordType range = end - beg;

  if(range == 0)
  {
    range = axom::numeric_limits<CoordType>::max();
  }

  return (std::rand() % range) + beg;
}

template <typename CoordType, int DIM>
Point<CoordType, DIM> randomPoint(CoordType beg, CoordType end)
{
  Point<CoordType, DIM> pt;
  for(int i = 0; i < DIM; ++i)
  {
    pt[i] = randomInt(beg, end);
  }

  return pt;
}
}  // end anonymous namespace

TEST(spin_morton, test_max_set_bit)
{
  SLIC_INFO("Checks that MortonBase's maxSetBit function works properly");

  using CoordType = int;
  using MortonIndexType = std::size_t;

  axom::spin::Mortonizer<CoordType, MortonIndexType, 2> morton2;
  EXPECT_EQ(morton2.maxSetBit(0), 0);

  int maxBit = axom::numeric_limits<CoordType>::digits;
  for(int i = 0; i <= maxBit; ++i)
  {
    int val = 1 << i;
    EXPECT_EQ(morton2.maxSetBit(val), i);

    for(int j = 0; j < MAX_ITER; ++j)
    {
      int randVal = randomInt<int>(0, val);
      EXPECT_EQ(morton2.maxSetBit(val + randVal), i);
    }
  }
}

TEST(spin_morton, test_mortonizer)
{
  using namespace axom::spin;

  SLIC_INFO("Testing Morton conversion on some simple points");

  axom::slic::setLoggingMsgLevel(axom::slic::message::Debug);
  using MortonIndexType = std::size_t;

  Point<int, 2> pt2(2);  // (0b10, 0b10)
  Mortonizer<int, MortonIndexType, 2> morton2;
  MortonIndexType mIdx2 = morton2.mortonize(pt2);       // interleaved bits
  EXPECT_EQ(mIdx2, static_cast<MortonIndexType>(0xC));  // 0b1100

  MortonIndexType mIdx2Alt = morton2.mortonize(pt2[0], pt2[1]);
  EXPECT_EQ(mIdx2, mIdx2Alt);
  EXPECT_EQ(morton2.demortonize(mIdx2Alt), pt2);

  // Example from doxygen documentation Point(6,3) has MortonIndex 29
  Point<int, 2> ptDoxygenExample;
  ptDoxygenExample[0] = 6;
  EXPECT_EQ(ptDoxygenExample[0], 0x6);  // 0b0110

  ptDoxygenExample[1] = 3;
  EXPECT_EQ(ptDoxygenExample[1], 0x3);  // 0b0011

  MortonIndexType mEx = morton2.mortonize(ptDoxygenExample);
  EXPECT_EQ(mEx, static_cast<MortonIndexType>(30));  // interleaved bits
  EXPECT_EQ(30, 0x1e);                               // 0b 0001 1110

  Point<int, 3> pt3(2);  // (0b10, 0b10, 0b10)
  Mortonizer<int, MortonIndexType, 3> morton3;
  MortonIndexType mIdx3 = morton3.mortonize(pt3);        // interleaved
                                                         // bits
  EXPECT_EQ(mIdx3, static_cast<MortonIndexType>(0x38));  // 0b 0011 1000

  MortonIndexType mIdx3Alt = morton3.mortonize(pt3[0], pt3[1], pt3[2]);
  EXPECT_EQ(mIdx3, mIdx3Alt);
  EXPECT_EQ(morton3.demortonize(mIdx3Alt), pt3);

  axom::slic::setLoggingMsgLevel(axom::slic::message::Info);

  // The following will not compile -- static_assert CoordType must be integral
  //Mortonizer<double,2> mortonD;
}

template <typename CoordType, typename MortonIndexType, int DIM>
void testMortonizer()
{
  using namespace axom::spin;

  using GridPoint = Point<CoordType, DIM>;

  int maxBits = axom::spin::Mortonizer<CoordType, MortonIndexType, DIM>::maxBitsPerCoord();
  SLIC_INFO("\tMax bits per dimension: " << axom::numeric_limits<CoordType>::digits);
  SLIC_INFO("\tMax unique bits per dimension: " << maxBits);

  int maxIter = std::min(1 << (maxBits - 1), MAX_ITER);

  SLIC_DEBUG("Testing " << maxIter << " random points");
  for(int i = 0; i < maxIter; ++i)
  {
    GridPoint origPt = randomPoint<CoordType, DIM>(0, 1 << maxBits);
    SLIC_DEBUG("\tOriginal point: " << origPt);

    MortonIndexType mortonIdx = convertPointToMorton<MortonIndexType>(origPt);
    SLIC_DEBUG("\tMorton index: " << mortonIdx);

    GridPoint convertedPt = convertMortonToPoint<CoordType, DIM>(mortonIdx);
    SLIC_DEBUG("\tConverted point: " << convertedPt << "\n..");

    EXPECT_EQ(origPt, convertedPt);

    MortonIndexType convertedMortonIdx = convertPointToMorton<MortonIndexType>(convertedPt);
    EXPECT_EQ(mortonIdx, convertedMortonIdx);
  }
}

template <int DIM>
void testIntegralTypes()
{
  SLIC_INFO("Testing char in " << DIM << "d -- ");
  testMortonizer<std::int8_t, std::uint8_t, DIM>();
  testMortonizer<std::int8_t, std::uint16_t, DIM>();
  testMortonizer<std::int8_t, std::uint32_t, DIM>();
  testMortonizer<std::int8_t, std::uint64_t, DIM>();

  SLIC_INFO("Testing uchar in " << DIM << "d -- ");
  testMortonizer<std::uint8_t, std::uint8_t, DIM>();
  testMortonizer<std::uint8_t, std::uint16_t, DIM>();
  testMortonizer<std::uint8_t, std::uint32_t, DIM>();
  testMortonizer<std::uint8_t, std::uint64_t, DIM>();

  // --
  SLIC_INFO("Testing short in " << DIM << "d -- ");
  testMortonizer<std::int16_t, std::uint8_t, DIM>();
  testMortonizer<std::int16_t, std::uint16_t, DIM>();
  testMortonizer<std::int16_t, std::uint32_t, DIM>();
  testMortonizer<std::int16_t, std::uint64_t, DIM>();

  SLIC_INFO("Testing ushort in " << DIM << "d -- ");
  testMortonizer<std::uint16_t, std::uint8_t, DIM>();
  testMortonizer<std::uint16_t, std::uint16_t, DIM>();
  testMortonizer<std::uint16_t, std::uint32_t, DIM>();
  testMortonizer<std::uint16_t, std::uint64_t, DIM>();

  // --
  SLIC_INFO("Testing int in " << DIM << "d -- ");
  testMortonizer<std::int32_t, std::uint8_t, DIM>();
  testMortonizer<std::int32_t, std::uint16_t, DIM>();
  testMortonizer<std::int32_t, std::uint32_t, DIM>();
  testMortonizer<std::int32_t, std::uint64_t, DIM>();

  SLIC_INFO("Testing uint in " << DIM << "d -- ");
  testMortonizer<std::uint32_t, std::uint8_t, DIM>();
  testMortonizer<std::uint32_t, std::uint16_t, DIM>();
  testMortonizer<std::uint32_t, std::uint32_t, DIM>();
  testMortonizer<std::uint32_t, std::uint64_t, DIM>();

  // --
  SLIC_INFO("Testing long long in " << DIM << "d -- ");
  testMortonizer<std::int64_t, std::uint8_t, DIM>();
  testMortonizer<std::int64_t, std::uint16_t, DIM>();
  testMortonizer<std::int64_t, std::uint32_t, DIM>();
  testMortonizer<std::int64_t, std::uint64_t, DIM>();

  SLIC_INFO("Testing ull in " << DIM << "d -- ");
  testMortonizer<std::uint64_t, std::uint8_t, DIM>();
  testMortonizer<std::uint64_t, std::uint16_t, DIM>();
  testMortonizer<std::uint64_t, std::uint32_t, DIM>();
  testMortonizer<std::uint64_t, std::uint64_t, DIM>();
}

TEST(spin_morton, test_integral_types_2D)
{
  SLIC_INFO(
    "*** Testing morton indexing in 2D with different coord and Morton index "
    "types");

  const int DIM = 2;
  testIntegralTypes<DIM>();
}

TEST(spin_morton, test_integral_types_3D)
{
  SLIC_INFO(
    "*** Testing morton indexing in 3D with different coord and Morton index "
    "types");

  const int DIM = 3;
  testIntegralTypes<DIM>();
}

TEST(spin_morton, test_point_hasher)
{
  using namespace axom::spin;

  SLIC_INFO(
    "** Here we test the point hasher which can be used e.g. in an "
    "unordered_map");

  axom::slic::setLoggingMsgLevel(axom::slic::message::Debug);

  using CoordType = int;
  PointHash<CoordType> ptHash;
  Point<CoordType, 1> p1(2);
  std::size_t exp = 0x2;  // 0b10
  EXPECT_EQ(ptHash(p1), exp);

  Point<CoordType, 2> p2(2);
  exp = 0xC;  // 0b1100
  EXPECT_EQ(ptHash(p2), exp);
  p2[0] = 0x38;  // 0b 0011 1000
  p2[1] = 0x0C;  // 0b 0000 1100
  exp = 0x5E0;   // 0b 0101 1110 0000
  EXPECT_EQ(ptHash(p2), exp);

  Point<CoordType, 3> p3(2);
  exp = 0x38;  // 0b 0011 1000
  EXPECT_EQ(ptHash(p3), exp);
  p3[0] = 0x1F;  // 0b 0001 1111
  p3[1] = 0X04;  // 0b 0000 0100
  p3[2] = 0x10;  // 0b 0001 0000
  exp = 051311;  // in octal (read bits bottom up, left to right)
  EXPECT_EQ(ptHash(p3), exp);

  Point<CoordType, 4> p4(2);
  exp = 0XF0;  // 0b 1111 0000
  EXPECT_EQ(ptHash(p4), exp);
  p4[0] = 0x3F;    // 0b 0011 1111
  p4[1] = 0x08;    // 0b 0000 1000
  p4[2] = 0x10;    // 0b 0001 0000
  p4[3] = 0x20;    // 0b 0010 0000
  exp = 0x953111;  // in hex (read bits bottom up, left to right)
  EXPECT_EQ(ptHash(p4), exp);

  axom::slic::setLoggingMsgLevel(axom::slic::message::Info);
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger(axom::slic::message::Info);

#ifdef MORTON_TESTER_SHOULD_SEED
  std::srand(std::time(0));
#else
  std::srand(105);
#endif

  result = RUN_ALL_TESTS();

  return result;
}
