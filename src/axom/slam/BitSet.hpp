// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file BitSet.hpp
 *
 * \brief Contains a BitSet class for manipulating ordered sequences of bits.
 */

#ifndef SLAM_BITSET_H_
#define SLAM_BITSET_H_

#include "axom/config.hpp"
#include "axom/core/Array.hpp"
#include "axom/core/execution/atomics.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/core/utilities/BitUtilities.hpp"
#include "axom/slic.hpp"

#include <vector>

namespace axom
{
namespace slam
{
// Forward declare the BitSet class and some operator functions

class BitSet;

/**
 * \brief Union operator for two bit sets.
 *
 *
 * \param lhs The first bitset
 * \param rhs The second bitset
 * \return The union of the two input bitsets, where the i^th bit is set
 * it is set in either lhs or rhs (e.g. b[i] = lhs[i] || rhs[i]).
 *
 * \pre lhs.size() == rhs.size()
 */
BitSet operator|(const BitSet& lhs, const BitSet& rhs);

/**
 * \brief Intersection operator for two bit sets.
 *
 *
 * \param lhs The first bitset
 * \param rhs The second bitset
 * \return The intersection of the two input bitsets, where the i^th bit is set
 * it is set in both lhs and rhs (e.g. b[i] = lhs[i] && rhs[i]).
 *
 * \pre lhs.size() == rhs.size()
 */
BitSet operator&(const BitSet& lhs, const BitSet& rhs);

/**
 * \brief Exclusive or (xor) operator for two bit sets.
 *
 *
 * \param lhs The first bitset
 * \param rhs The second bitset
 * \return The xor of the two input bitsets, where the i^th bit is set
 * it is set in either lhs or rhs but not both (e.g. b[i] = lhs[i] ^ rhs[i]).
 *
 * \pre lhs.size() == rhs.size()
 */
BitSet operator^(const BitSet& lhs, const BitSet& rhs);

/**
 * \brief Set difference operator for two bit sets.
 *
 *
 * \param lhs The first bitset
 * \param rhs The second bitset
 * \return The set difference of the two input bitsets, where the i^th bit is
 * set it is set in lhs and not rhs (e.g. b[i] = lhs[i] && ! rhs[i]).
 *
 * \pre lhs.size() == rhs.size()
 */
BitSet operator-(const BitSet& lhs, const BitSet& rhs);

/**
 * \class BitSet
 * \brief A bitset class
 *
 * This class supports bitwise manipulation operations (e.g. set intersection,
 * union and difference) on an ordered set of bits. The class has a similar
 * interface to std::bitset and boost::dynamic_bitset, but with the following
 * differences:
 *   - The size of the bitset is supplied at runtime in the class constructor.
 *   However, BitSet does not currently support changing the size after
 *   construction.
 *   - We do not support the random access operation ( operator[](index) ).
 *   The value of individual bits can be checked using the test function
 *   (e.g. bset.test(i) )
 *   - There is no support for directly initializing the bits in the bitset
 *   (e.g. via strings).
 *
 * The individual bits in the bitset are packed into a contiguous array of
 * Words (an unsigned integer type). Many bitset operations such as count()
 * and flip() are performed at the granularity of a word.
 *
 * BitSet uses the same interface as boost::dynamic_bitset to enumerate the
 * set bits. Specifically, find_first() returns the index of the first set bit.
 * and find_next(idx) returns the next bit that is set after bit index idx.
 * BitSet::npos is used as a sentinel to indicate no more set bits.
 */
class BitSet
{
public:
  using Index = int;

#if (defined(__x86_64__) && defined(__GNUC__)) || (defined(_WIN64) && (_MSC_VER >= 1600))
  using Word = std::uint64_t;
#else
  using Word = std::uint32_t;
#endif

  // TODO: update using a policy
  using ArrayType = axom::Array<Word, 1>;

  static constexpr Index npos = -2;
  static constexpr int BitsPerWord = axom::utilities::BitTraits<Word>::BITS_PER_WORD;

private:
  static constexpr int LG_BITS_PER_WORD = axom::utilities::BitTraits<Word>::LG_BITS_PER_WORD;

public:
  /**
   * \brief BitSet class constructor
   *
   * \param numBits The number of bits in the bitset.
   * \pre numBits must be non-negative
   * \post bset.size() == numBits
   * \post All bits will be off
   */
  explicit BitSet(int numBits = 0, int allocatorID = axom::getDefaultAllocatorID())
  {
    SLIC_ASSERT_MSG(numBits >= 0, "slam::BitSet must be initialized with a non-zero number of bits");

    m_numBits = axom::utilities::max(numBits, 0);
    axom::IndexType numWords = (m_numBits == 0) ? 1 : 1 + (m_numBits - 1) / BitsPerWord;

    m_data = ArrayType(axom::ArrayOptions::Uninitialized {}, numWords, numWords, allocatorID);
    m_data.fill(0);
  }

  /** \brief Equality operator for two bitsets */
  bool operator==(const BitSet& other) const;

  /** \brief Inequality operator for two bitsets */
  bool operator!=(const BitSet& other) const { return !(operator==(other)); }

  /*!
   * \brief Gets the underlying data of a BitSet.
   */
  AXOM_HOST_DEVICE const Word* data() const { return m_data.data(); }

public:
  /// \name Bitset bitwise assignment operators
  /// @{

  /**
   * \brief BitSet union-assignment operator
   *
   * \param other The other bitset
   * \pre other.size() == size()
   * \return A reference to the modified bitset
   *
   * Set union of the current instance and other. The i^th
   * bit will be set if it is set in this bitset or in \a other
   */
  BitSet& operator|=(const BitSet& other);

  /**
   * \brief BitSet intersection-assignment operator
   *
   * \param other The other bitset
   * \pre other.size() == size()
   * \return A reference to the modified bitset
   *
   * Set intersection of the current instance and other. The i^th
   * bit will be set if it is set in this bitset and in \a other
   */
  BitSet& operator&=(const BitSet& other);

  /**
   * \brief BitSet exclusive-or-assignment operator
   *
   * \param other The other bitset
   * \pre other.size() == size()
   * \return A reference to the modified bitset
   *
   * Set exclusive-or of the current instance and other. The i^th
   * bit will be set if it is set in this bitset or in \a other,
   * but not both.
   */
  BitSet& operator^=(const BitSet& other);

  /**
   * \brief BitSet difference-assignment operator
   *
   * \param other The other bitset
   * \pre other.size() == size()
   * \return A reference to the modified bitset
   *
   * Set difference of the current instance and other. The i^th
   * bit will be set if it is set in this bitset and not in \a other
   */
  BitSet& operator-=(const BitSet& other);

  /// @}

public:
  /// \name Bitset iteration interface
  /// @{

  /**
   * \brief Finds the index of the first bit that is set in the bitset
   *
   * \return The index of the first set bit,
   * or BitSet::npos if no bits are set
   */
  Index find_first() const;

  /**
   * \brief Finds the index of the next set bit in the bitset after \a idx
   *
   * \param idx The starting index
   * \return The index of the first set bit after index \a idx
   * or BitSet::npos if none can be found
   * \note Will also return BitSet::npos if \a idx is BitSet::npos
   */
  Index find_next(Index idx) const;

  /// @}

public:
  /// \name Operations that affect all bits in the bitset
  /// @{

  /** \brief Returns the cardinality of the bitset */
  AXOM_HOST_DEVICE int size() const { return m_numBits; }

  /** \brief Returns the number of bits that are set */
  int count() const;

  /** \brief Clears all bits in the bitset */
  void clear();

  /** \brief Sets all bits in the bitset */
  void set();

  /** \brief Toggles all bits in the bitset */
  void flip();

  /**
   * \brief Checks if the bitset instance is valid
   *
   * \return True if the bitset is valid, false otherwise
   *
   * A bitset is valid if it has sufficient storage for size() bits.
   * If we have storage for more than size() bits, none of these additional
   * bits are set.
   */
  bool isValid() const;

  /// @}

public:
  /// \name Operations that affect a single bits in the bitset
  /// @{

  /**
   * \brief Clears bit at index \a idx
   *
   * \pre \a idx must be between 0 and bitset.size()
   */
  void clear(Index idx) { getWord(idx) &= ~mask(idx); }

  /**
   * \brief Sets bit at index \a idx
   *
   * \pre \a idx must be between 0 and bitset.size()
   */
  void set(Index idx) { getWord(idx) |= mask(idx); }

  /**
   * \brief Toggles bit at index \a idx
   *
   * \pre \a idx must be between 0 and bitset.size()
   */
  void flip(Index idx) { getWord(idx) ^= mask(idx); }

  /**
   * \brief Tests the bit at index \a idx
   * \return True if \a idx is valid and its bit is set, false otherwise
   */
  bool test(Index idx) const
  {
    return (idx >= 0 && idx < m_numBits) &&   // idx is in range
      (getWord(idx) & mask(idx)) != Word(0);  // and its bit is set
  }

  /// @}

  /// \name Atomic versions of single-bit operations
  /// @{

  /**
   * \brief Clears bit at index \a idx
   *
   * \pre \a idx must be between 0 and bitset.size()
   */
  void atomicClear(Index idx) { axom::atomicAnd<axom::auto_atomic>(&getWord(idx), ~mask(idx)); }

  /**
   * \brief Sets bit at index \a idx
   *
   * \pre \a idx must be between 0 and bitset.size()
   */
  AXOM_HOST_DEVICE void atomicSet(Index idx)
  {
    axom::atomicOr<axom::auto_atomic>(&getWord(idx), mask(idx));
  }

  /**
   * \brief Toggles bit at index \a idx
   *
   * \pre \a idx must be between 0 and bitset.size()
   */
  void atomicFlip(Index idx) { axom::atomicXor<axom::auto_atomic>(&getWord(idx), mask(idx)); }

  /// @}
private:
  /**
   * \brief Gets the index of the word containing bit at index \a idx
   *
   * \param idx The index of the word containing the desired bit
   * \param checkIndexValid Option to enable bounds checking to ensure
   * that \a idx is within range [0, size() )
   */
  AXOM_HOST_DEVICE
  Word& getWord(Index idx, bool checkIndexValid = true)
  {
    if(checkIndexValid)
    {
      checkValidIndex(idx);
    }

    const Index wIdx = idx / BitsPerWord;
    return m_data[wIdx];
  }

  /**
   * \brief Const implementation of getWord()
   * \sa getWord()
   */
  AXOM_HOST_DEVICE
  const Word& getWord(Index idx, bool checkIndexValid = true) const
  {
    if(checkIndexValid)
    {
      checkValidIndex(idx);
    }

    const Index wIdx = idx / BitsPerWord;
    return m_data[wIdx];
  }

  /**
   * \brief Returns a bitmask for the desired bit within the word
   * containing \a idx
   *
   * \param idx The index of the desired bit
   */
  AXOM_HOST_DEVICE
  Word mask(Index idx) const
  {
    const Index wOffset = idx % BitsPerWord;
    return Word(1) << wOffset;
  }

  /**
   * \brief Returns a bitmask for the bits in the last word of the bitset.
   *
   * This function is only valid when isLasteWordFull() is false
   * \sa isLastWordFull()
   */
  Word lastWordMask() const { return mask(m_numBits) - 1; }

  /**
   * \brief Checks if index \idx corresponds to a valid index
   *
   * \note This function is a no-op in Release builds
   */
  AXOM_HOST_DEVICE
  void checkValidIndex(Index idx) const
  {
    AXOM_UNUSED_VAR(idx);
    SLIC_ASSERT_MSG(idx >= 0 && idx < m_numBits,
                    "slam::Bitset attempted to out of range bit " << idx << ". Valid range is [0, "
                                                                  << m_numBits << ").");
  }

  /**
   * \brief Predicate to determine if we need special processing
   * for the final word of the bitset
   *
   * The last word is full when the bitset has exactly
   * m_words * BitsPerWord bits
   */
  bool isLastWordFull() const
  {
    const int lg = (1 << (LG_BITS_PER_WORD)) - 1;
    return (m_numBits & lg) == 0;
  }

private:
  ArrayType m_data;

  int m_numBits;
};

}  // end namespace slam
}  // end namespace axom

#endif  //  SLAM_BITSET_H_
