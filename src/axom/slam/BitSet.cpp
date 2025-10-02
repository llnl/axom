// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/BitUtilities.hpp"

#include "axom/slam/BitSet.hpp"

namespace axom
{
namespace slam
{
void BitSet::clear()
{
  for(int i = 0; i < m_data.size(); ++i)
  {
    m_data[i] = 0;
  }
}

void BitSet::set()
{
  if(m_numBits == 0)
  {
    return;
  }

  const Word ones = ~Word(0);
  for(int i = 0; i < m_data.size() - 1; ++i)
  {
    m_data[i] = ones;
  }

  // Handle last word
  m_data[m_data.size() - 1] = isLastWordFull() ? ones : lastWordMask();
}

void BitSet::flip()
{
  if(m_numBits == 0)
  {
    return;
  }

  const Word ones = ~Word(0);
  for(int i = 0; i < m_data.size() - 1; ++i)
  {
    m_data[i] ^= ones;
  }

  // Handle last word
  m_data[m_data.size() - 1] ^= isLastWordFull() ? ones : lastWordMask();
}

int BitSet::count() const
{
  int ctr = 0;

  for(int i = 0; i < m_data.size(); ++i)
  {
    ctr += axom::utilities::popcount(m_data[i]);
  }
  return ctr;
}

bool BitSet::isValid() const
{
  bool valid = true;

  if(m_numBits < 0 || m_data.size() < 0)
  {
    valid = false;
  }

  if(m_numBits == 0)
  {
    if(m_data.size() != 1 || m_data[0] != Word(0))
    {
      valid = false;
    }
  }
  else
  {
    // check num words vs. num bits
    int expWords = (m_numBits - 1) / BitsPerWord + 1;
    if(expWords != m_data.size())
    {
      valid = false;
    }

    // check that highest bits are not set
    if(!isLastWordFull())
    {
      const Word& lastWord = getWord(m_numBits, false);

      Word upperMask = ~lastWordMask();
      if((lastWord & upperMask) != 0)
      {
        valid = false;
      }
    }
  }

  return valid;
}

BitSet::Index BitSet::find_first() const { return find_next(-1); }

BitSet::Index BitSet::find_next(Index idx) const
{
  // Handle boundary cases
  if(idx == npos || (idx + 1 >= m_numBits))
  {
    return npos;
  }

  Index startWordIdx = 0;
  if(idx >= 0)
  {
    checkValidIndex(idx);

    const Index startIdx = idx + 1;
    startWordIdx = startIdx / BitsPerWord;

    // Check for next set bit in current word
    const Index startOffset = startIdx % BitsPerWord;

    const Word startWord = m_data[startWordIdx] >> startOffset;
    if(startWord != Word(0))
    {
      return (startWordIdx * BitsPerWord) + axom::utilities::countr_zero(startWord << startOffset);
    }

    ++startWordIdx;
  }

  // If not in current word, check remaining words
  for(int i = startWordIdx; i < m_data.size(); ++i)
  {
    const Word& w = m_data[i];
    if(w != Word(0))
    {
      return (i * BitsPerWord) + axom::utilities::countr_zero(w);
    }
  }
  return BitSet::npos;
}

BitSet& BitSet::operator|=(const BitSet& other)
{
  SLIC_ASSERT_MSG(size() == other.size(),
                  "slam::BitSet Sizes must be the same for bit set operators."
                    << " In operator|=(), BitSet has size " << size()
                    << " and other BitSet has size " << other.size() << ".");

  for(int i = 0; i < m_data.size(); ++i)
  {
    m_data[i] |= other.m_data[i];
  }

  return *this;
}

BitSet& BitSet::operator&=(const BitSet& other)
{
  SLIC_ASSERT_MSG(size() == other.size(),
                  "slam::BitSet Sizes must be the same for bit set operators."
                    << " In operator&=(), BitSet has size " << size()
                    << " and other BitSet has size " << other.size() << ".");

  for(int i = 0; i < m_data.size(); ++i)
  {
    m_data[i] &= other.m_data[i];
  }

  return *this;
}

BitSet& BitSet::operator^=(const BitSet& other)
{
  SLIC_ASSERT_MSG(size() == other.size(),
                  "slam::BitSet Sizes must be the same for bit set operators."
                    << " In operator^=(), BitSet has size " << size()
                    << " and other BitSet has size " << other.size() << ".");

  for(int i = 0; i < m_data.size(); ++i)
  {
    m_data[i] ^= other.m_data[i];
  }

  return *this;
}

BitSet& BitSet::operator-=(const BitSet& other)
{
  SLIC_ASSERT_MSG(size() == other.size(),
                  "slam::BitSet Sizes must be the same for bit set operators."
                    << " In operator-=(), BitSet has size " << size()
                    << " and other BitSet has size " << other.size() << ".");

  for(int i = 0; i < m_data.size(); ++i)
  {
    m_data[i] &= ~other.m_data[i];
  }

  return *this;
}

BitSet operator|(const BitSet& lhs, const BitSet& rhs)
{
  BitSet s(lhs);
  s |= rhs;
  return s;
}

BitSet operator&(const BitSet& lhs, const BitSet& rhs)
{
  BitSet s(lhs);
  s &= rhs;
  return s;
}

BitSet operator^(const BitSet& lhs, const BitSet& rhs)
{
  BitSet s(lhs);
  s ^= rhs;
  return s;
}

BitSet operator-(const BitSet& lhs, const BitSet& rhs)
{
  BitSet s(lhs);
  s -= rhs;
  return s;
}

bool BitSet::operator==(const BitSet& other) const
{
  if(size() != other.size() || m_numBits != other.m_numBits)
  {
    return false;
  }

  return m_data == other.m_data;
}

}  // end namespace slam
}  // end namespace axom
