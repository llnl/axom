// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core/Types.hpp"
#include "axom/core/NumericLimits.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/core/utilities/BitUtilities.hpp"

#include <limits>

namespace
{
constexpr int NRAND_SAMPLES = 1000;

// Returns a random integer value between 0 and the max value
template <typename T>
T random_int()
{
  static_assert(std::is_integral<T>::value, "T must be an integral type");

  constexpr T max_int = axom::numeric_limits<T>::max();
  constexpr double max_d = static_cast<double>(max_int);

  const auto val = axom::utilities::random_real(0., max_d);
  return static_cast<T>(val);
}

std::uint64_t shifted(int bit) { return static_cast<std::uint64_t>(1) << bit; }

}  // namespace

TEST(core_bit_utilities, trailingZeroes)
{
  constexpr std::uint64_t ZERO = std::uint64_t(0);
  constexpr int BITS = axom::utilities::BitTraits<std::uint64_t>::BITS_PER_WORD;
  ASSERT_EQ(64, BITS);

  // Axom's countr_zero will return 64 when given 0
  {
    EXPECT_EQ(BITS, axom::utilities::countr_zero(ZERO));
  }

  // Test with a known trailing bit
  for(int i = 0; i < BITS; ++i)
  {
    std::uint64_t val = ::shifted(i);
    EXPECT_EQ(i, axom::utilities::countr_zero(val));

    // Value doesn't change when you set bits to left of trailing bit
    for(int j = i + 1; j < BITS; ++j)
    {
      std::uint64_t val2 = ::shifted(i) + ::shifted(j);
      EXPECT_EQ(axom::utilities::countr_zero(val), axom::utilities::countr_zero(val2));
    }
  }

  // Test trailing zeros on some random numbers
  for(int n = 0; n < NRAND_SAMPLES; ++n)
  {
    std::uint64_t rand_val = ::random_int<std::uint64_t>();

    // Explicitly find first bit and check that it matches
    int bit = 0;
    for(; bit < BITS; ++bit)
    {
      if(rand_val & shifted(bit))
      {
        break;
      }
    }
    EXPECT_EQ(bit, axom::utilities::countr_zero(rand_val));
  }
}

TEST(core_bit_utilities, popcount)
{
  constexpr std::uint64_t ZERO = std::uint64_t(0);
  constexpr int BITS = axom::utilities::BitTraits<std::uint64_t>::BITS_PER_WORD;
  ASSERT_EQ(64, BITS);

  // Test pop count when zero bits are set
  {
    EXPECT_EQ(0, axom::utilities::popcount(ZERO));
    EXPECT_EQ(BITS, axom::utilities::popcount(~ZERO));
  }

  // Test pop count when one bit is set
  for(int i = 0; i < BITS; ++i)
  {
    std::uint64_t val = ::shifted(i);
    EXPECT_EQ(1, axom::utilities::popcount(val));
    EXPECT_EQ(BITS - 1, axom::utilities::popcount(~val));
  }

  // Test pop count when two bits are set
  for(int i = 0; i < BITS; ++i)
  {
    for(int j = 0; j < i; ++j)
    {
      std::uint64_t val = shifted(i) + shifted(j);
      EXPECT_EQ(2, axom::utilities::popcount(val));
      EXPECT_EQ(BITS - 2, axom::utilities::popcount(~val));
    }
  }

  // Test pop count when three bits are set
  for(int i = 0; i < BITS; ++i)
  {
    for(int j = 0; j < i; ++j)
    {
      for(int k = 0; k < j; ++k)
      {
        std::uint64_t val = shifted(i) + shifted(j) + shifted(k);
        EXPECT_EQ(3, axom::utilities::popcount(val));
        EXPECT_EQ(BITS - 3, axom::utilities::popcount(~val));
      }
    }
  }

  // Test pop count on some random numbers
  for(int n = 0; n < NRAND_SAMPLES; ++n)
  {
    auto val = ::random_int<std::uint64_t>();

    int bits = 0;
    for(int i = 0; i < BITS; ++i)
    {
      if(val & shifted(i))
      {
        ++bits;
      }
    }

    EXPECT_EQ(bits, axom::utilities::popcount(val));
  }
}

TEST(core_bit_utilities, countl_zero)
{
  constexpr std::int32_t ZERO = std::int32_t(0);
  constexpr int BITS = axom::utilities::BitTraits<std::uint32_t>::BITS_PER_WORD;
  ASSERT_EQ(32, BITS);

  // Axom's countl_zero will return 32 when given 0
  EXPECT_EQ(BITS, axom::utilities::countl_zero(ZERO));

  for(int i = 0; i < BITS; ++i)
  {
    std::int32_t val = ::shifted(i);
    EXPECT_EQ(BITS - i - 1, axom::utilities::countl_zero(val));

    // Value doesn't change if you set bits to right of leading zero
    for(int j = 0; j < i; ++j)
    {
      std::int32_t val2 = ::shifted(i) + ::shifted(j);
      EXPECT_EQ(axom::utilities::countl_zero(val), axom::utilities::countl_zero(val2));
    }
  }

  // Test leading zeros on some random numbers
  for(int n = 0; n < NRAND_SAMPLES; ++n)
  {
    std::int32_t rand_val = ::random_int<std::int32_t>();

    // Explicitly find first bit and check that it matches
    int bit = 0;
    for(; bit < BITS; ++bit)
    {
      if(rand_val & shifted(BITS - bit - 1))
      {
        break;
      }
    }
    EXPECT_EQ(bit, axom::utilities::countl_zero(rand_val));
  }
}

TEST(core_bit_utilities, setbit_bitisset)
{
  const std::uint32_t pattern = 0xaaaaaaaa;
  for(size_t bit = 0; bit < axom::utilities::BitTraits<std::uint32_t>::BITS_PER_WORD; bit++)
  {
    EXPECT_EQ(axom::utilities::bitIsSet(pattern, bit), ((bit & 1) == 1));
  }

  const bool bitvals[] = {false, true, true, false, true, true, false, true};
  std::uint8_t value = 0;
  for(size_t i = 0; i < axom::utilities::BitTraits<std::uint8_t>::BITS_PER_WORD; i++)
  {
    axom::utilities::setBit(value, i, bitvals[i]);
    for(size_t b = 0; b <= i; b++)
    {
      EXPECT_EQ(axom::utilities::bitIsSet(value, b), bitvals[b]);
    }
  }

  for(size_t i = 0; i < axom::utilities::BitTraits<std::uint8_t>::BITS_PER_WORD; i++)
  {
    axom::utilities::setBit(value, i, false);
    EXPECT_EQ(axom::utilities::bitIsSet(value, i), false);
  }
}
