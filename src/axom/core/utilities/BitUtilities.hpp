// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *
 * \file BitUtilities.hpp
 *
 * \brief Header file containing bitwise utility functions.
 *
 */

#ifndef AXOM_BIT_UTILITIES_HPP
#define AXOM_BIT_UTILITIES_HPP

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"

// CUDA intrinsics: https://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__INTRINSIC__INT.html
// HIP intrinsics: https://rocm.docs.amd.com/projects/HIP/en/latest/reference/kernel_language.html

// Check for and setup defines for platform-specific intrinsics
// Note: `__GNUC__` is defined for the gnu, clang and intel compilers
#if defined(AXOM_USE_CUDA)
  // Intrinsics included implicitly

#elif defined(AXOM_USE_HIP)
  #include <hip/hip_runtime.h>

#elif defined(_WIN64) && (_MSC_VER >= 1600)
  #define _AXOM_CORE_USE_INTRINSICS_MSVC
  #include <intrin.h>
#elif defined(__x86_64__) && defined(__GNUC__)
  #define _AXOM_CORE_USE_INTRINSICS_GCC
  #include <x86intrin.h>
#elif defined(__powerpc64__) && (defined(__GNUC__) || defined(__ibmxl__))
  #define _AXOM_CORE_USE_INTRINSICS_PPC
#endif

namespace axom
{
namespace utilities
{
/**
 * Helper traits template for determining the number of bits and bytes
 * required to represent a given type \a T.
 */
template <typename T>
struct BitTraits;

template <>
struct BitTraits<std::int64_t>
{
  constexpr static int NUM_BYTES = 8;
  constexpr static int BITS_PER_WORD = NUM_BYTES << 3;
  constexpr static int LG_BITS_PER_WORD = 6;
};

template <>
struct BitTraits<std::int32_t>
{
  constexpr static int NUM_BYTES = 4;
  constexpr static int BITS_PER_WORD = NUM_BYTES << 3;
  constexpr static int LG_BITS_PER_WORD = 5;
};

template <>
struct BitTraits<std::int16_t>
{
  constexpr static int NUM_BYTES = 2;
  constexpr static int BITS_PER_WORD = NUM_BYTES << 3;
  constexpr static int LG_BITS_PER_WORD = 4;
};

template <>
struct BitTraits<std::int8_t>
{
  constexpr static int NUM_BYTES = 1;
  constexpr static int BITS_PER_WORD = NUM_BYTES << 3;
  constexpr static int LG_BITS_PER_WORD = 3;
};

template <>
struct BitTraits<std::uint64_t>
{
  constexpr static int NUM_BYTES = 8;
  constexpr static int BITS_PER_WORD = NUM_BYTES << 3;
  constexpr static int LG_BITS_PER_WORD = 6;
};

template <>
struct BitTraits<std::uint32_t>
{
  constexpr static int NUM_BYTES = 4;
  constexpr static int BITS_PER_WORD = NUM_BYTES << 3;
  constexpr static int LG_BITS_PER_WORD = 5;
};

template <>
struct BitTraits<std::uint16_t>
{
  constexpr static int NUM_BYTES = 2;
  constexpr static int BITS_PER_WORD = NUM_BYTES << 3;
  constexpr static int LG_BITS_PER_WORD = 4;
};

template <>
struct BitTraits<std::uint8_t>
{
  constexpr static int NUM_BYTES = 1;
  constexpr static int BITS_PER_WORD = NUM_BYTES << 3;
  constexpr static int LG_BITS_PER_WORD = 3;
};

/**
 * \brief Counts the number of trailing zeros in \a word
 * \accelerated
 * \return The number of zeros to the right of the first set bit in \word,
 * starting with the least significant bit, or 64 if \a word == 0.
 */
AXOM_HOST_DEVICE inline int countr_zero(std::uint64_t word) noexcept
{
  /* clang-format off */
#if defined(__CUDA_ARCH__) && defined(AXOM_USE_CUDA)
  return word != std::uint64_t(0) ? __ffsll(word) - 1 : 64;
#elif defined(__HIP_DEVICE_COMPILE__) && defined(AXOM_USE_HIP)
  return word != std::uint64_t(0) ? __ffsll(static_cast<unsigned long long int>(word)) - 1 : 64;
#elif defined(_AXOM_CORE_USE_INTRINSICS_MSVC)
  unsigned long cnt;
  return _BitScanForward64(&cnt, word) ? cnt : 64;
#elif defined(_AXOM_CORE_USE_INTRINSICS_GCC) || defined(_AXOM_CORE_USE_INTRINSICS_PPC)
  return word != std::uint64_t(0) ? __builtin_ctzll(word) : 64;
#else
  // Explicit implementation adapted from bit twiddling hacks
  // https://graphics.stanford.edu/~seander/bithacks.html#ZerosOnRightParallel
  // and modified for 64 bits:

  // cnt tracks the number of zero bits on the right of the first set bit
  int cnt = BitTraits<std::uint64_t>::BITS_PER_WORD;

  word &= -static_cast<std::int64_t>(word);
  if (word)                      { cnt--;     }
  if (word & 0x00000000FFFFFFFF) { cnt -= 32; }
  if (word & 0x0000FFFF0000FFFF) { cnt -= 16; }
  if (word & 0x00FF00FF00FF00FF) { cnt -=  8; }
  if (word & 0x0F0F0F0F0F0F0F0F) { cnt -=  4; }
  if (word & 0x3333333333333333) { cnt -=  2; }
  if (word & 0x5555555555555555) { cnt -=  1; }

  return cnt;
#endif
  /* clang-format on */
}

/*!
 * \brief Counts the number of set bits in \a word
 * \accelerated
 * \return number of bits in \a word that are set to 1
 */
AXOM_HOST_DEVICE inline int popcount(std::uint64_t word) noexcept
{
  /* clang-format off */
#if defined(__CUDA_ARCH__) && defined(AXOM_USE_CUDA)
  // Use CUDA intrinsic for popcount
  return __popcll(word);
#elif defined(__HIP_DEVICE_COMPILE__) && defined(AXOM_USE_HIP)
  // Use HIP intrinsic for popcount
  return __popcll(word);
#elif defined(_AXOM_CORE_USE_INTRINSICS_MSVC)
  return static_cast<int>(__popcnt64(word));
#elif defined(_AXOM_CORE_USE_INTRINSICS_GCC) || defined(_AXOM_CORE_USE_INTRINSICS_PPC)
  return __builtin_popcountll(word);
#else
  // 64 bit popcount implementation from:
  // http://chessprogramming.wikispaces.com/Population+Count#SWARPopcount

  using Word = std::uint64_t;

  const Word masks[] = {
    0x5555555555555555,  //  0 --  -1/3
    0x3333333333333333,  //  1 --  -1/5
    0x0F0F0F0F0F0F0F0F,  //  2 --  -1/17
    0x0101010101010101,  //  3 --  -1/255
  };

  // aggregate counts of 2-,4-,8- and 16-bits
  word = word - ((word >> 1) & masks[0]);
  word = (word & masks[1]) + ((word >> 2) & masks[1]);
  word = (word + (word >> 4)) & masks[2];
  return static_cast<int>((word * masks[3]) >> 56);
#endif
  /* clang-format off */
}

// gpu_macros_example_start
/*!
 * \brief Counts the number of leading zeros in \a word
 * \accelerated
 * \return The number of zeros to the left of the first set bit in \word,
 * starting with the least significant bit.
 */
AXOM_HOST_DEVICE inline std::int32_t countl_zero(std::int32_t word) noexcept
{
  /* clang-format off */
#if defined(__CUDA_ARCH__) && defined(AXOM_USE_CUDA)
  // Use CUDA intrinsic for count leading zeros
  return __clz(word);
#elif defined(__HIP_DEVICE_COMPILE__) && defined(AXOM_USE_HIP)
  // Use HIP intrinsic for count leading zeros
  return __clz(word);
#elif defined(_AXOM_CORE_USE_INTRINSICS_MSVC)
  unsigned long cnt;
  return _BitScanReverse(&cnt, word) ? 31 - cnt : 32;
#elif defined(_AXOM_CORE_USE_INTRINSICS_GCC) || defined(_AXOM_CORE_USE_INTRINSICS_PPC)
  return word != std::int32_t(0) ? __builtin_clz(word) : 32;
#else
  std::int32_t y {};
  std::int32_t n = 32;
  y = word >> 16; if(y != 0) { n -= 16; word = y;}
  y = word >>  8; if(y != 0) { n -=  8; word = y;}
  y = word >>  4; if(y != 0) { n -=  4; word = y;}
  y = word >>  2; if(y != 0) { n -=  2; word = y;}
  y = word >>  1; if(y != 0) { return std::int32_t(n - 2); }
  return std::int32_t(n - word);
#endif
  /* clang-format off */
}
// gpu_macros_example_end

/*!
 * \brief Determine whether the \a bit is set in \a flags.
 * \accelerated
 *
 * \tparam FlagType The integer type that contains the flags.
 * \tparam BitType  The index type that stores the bit index.
 *
 * \param[in] flags The flags whose bit we're testing.
 * \param[in] bit   The bit we're testing.
 *
 * \return True if the bit is set; false otherwise.
 */
template <typename FlagType, typename BitType>
AXOM_HOST_DEVICE
constexpr bool bitIsSet(FlagType flags, BitType bit)
{
  assert(static_cast<int>(bit) < BitTraits<FlagType>::BITS_PER_WORD << 3);
  return (flags & (1 << bit)) > 0;
}

/*!
 * \brief Set the bit in flags.
 * \accelerated
 *
 * \tparam FlagType The integer type that contains the flags.
 * \tparam BitType  The index type that stores the bit index.
 *
 * \param[inout] flags The flags whose bit we're setting.
 * \param[in] bit   The bit we're setting.
 * \param[in] value The value we're setting into the bit.
 */
template <typename FlagType, typename BitType>
AXOM_HOST_DEVICE
constexpr void setBit(FlagType &flags, BitType bit, bool value = true)
{
  assert(static_cast<int>(bit) < BitTraits<FlagType>::BITS_PER_WORD << 3);
  const auto mask = 1 << bit;
  flags = value ? (flags | mask) : (flags & ~mask);
}

/*!
 * \brief Set the bit in flags.
 * \accelerated
 *
 * \tparam FlagType The integer type that contains the flags.
 * \tparam BitType  The index type that stores the bit index.
 *
 * \param[inout] flags The flags whose bit we're setting.
 * \param[in] bit   The bit we're setting.
 */
template <typename FlagType, typename BitType>
AXOM_HOST_DEVICE
constexpr void setBitOn(FlagType &flags, BitType bit)
{
  assert(static_cast<int>(bit) < BitTraits<FlagType>::BITS_PER_WORD << 3);
  flags |= (1 << bit);
}

}  // namespace utilities
}  // namespace axom


// Undefine intrinsic defines from the top of the file
#undef _AXOM_CORE_USE_INTRINSICS_MSVC
#undef _AXOM_CORE_USE_INTRINSICS_GCC
#undef _AXOM_CORE_USE_INTRINSICS_PPC

#endif  // AXOM_BIT_UTILITIES_HPP
