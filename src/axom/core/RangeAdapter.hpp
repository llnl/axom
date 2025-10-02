// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef Axom_Core_RangeAdapter_HPP
#define Axom_Core_RangeAdapter_HPP

namespace axom
{
/*!
 * \class RangeAdapter
 *
 * \brief Simple adapter class that converts a pair of iterators into a range
 *  that can be iterated over in a range-based for-loop.
 */
template <typename IteratorType>
class RangeAdapter
{
public:
  RangeAdapter(IteratorType begin, IteratorType end) : m_begin(begin), m_end(end) { }
  /// \brief Returns an iterator to the beginning of the range.
  IteratorType begin() const { return m_begin; }
  /// \brief Returns an iterator to the end of the range.
  IteratorType end() const { return m_end; }

private:
  IteratorType m_begin;
  IteratorType m_end;
};

}  // namespace axom

#endif
