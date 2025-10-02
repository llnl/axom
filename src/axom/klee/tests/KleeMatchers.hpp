// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_KLEE_MATCHERS_HPP_
#define AXOM_KLEE_MATCHERS_HPP_

#include "axom/core.hpp"
#include "axom/primal.hpp"

#include "axom/klee/GeometryOperators.hpp"

#include "gmock/gmock.h"

// Note: After updating to the GMock in blt@0.5.3, the clang+nvcc compiler on LLNL's blueos
// complained about using GMock's MATCHER_P macros, as originally implemented.
// This required explicitly generating equivalent classes for the Matcher functionality.

namespace axom
{
namespace klee
{
namespace test
{
constexpr double tolerance = 1e-10;

/// Create a GMock Matcher for an axom::numerics:Matrix
template <typename T>
class AlmostEqMatrixMatcher
{
public:
  using is_gtest_matcher = void;

  explicit AlmostEqMatrixMatcher(const axom::numerics::Matrix<T>& mat) : m_mat(mat) { }

  bool MatchAndExplain(const axom::numerics::Matrix<T>& other, std::ostream*) const
  {
    if(other.getNumRows() != m_mat.getNumRows() || other.getNumColumns() != m_mat.getNumColumns())
    {
      return false;
    }

    for(axom::IndexType row = 0; row < other.getNumRows(); ++row)
    {
      for(axom::IndexType column = 0; column < other.getNumColumns(); ++column)
      {
        if(!::testing::Matches(::testing::DoubleNear(m_mat(row, column), tolerance))(
             other(row, column)))
        {
          return false;
        }
      }
    }

    return true;
  }

  void DescribeTo(std::ostream*) const { }
  void DescribeNegationTo(std::ostream*) const { }

  friend bool operator==(const AlmostEqMatrixMatcher& lhs, const axom::numerics::Matrix<T>& rhs)
  {
    return lhs.MatchAndExplain(rhs, nullptr);
  }

  friend bool operator==(const axom::numerics::Matrix<T>& lhs, const AlmostEqMatrixMatcher& rhs)
  {
    return rhs == lhs;
  }

private:
  const axom::numerics::Matrix<T> m_mat;
};

template <typename T>
inline auto AlmostEqMatrix(const axom::numerics::Matrix<T>& mat)
{
  return AlmostEqMatrixMatcher<T>(mat);
}

// ----------------------------------------------------------------------------

/// Create a GMock Matcher for (numeric array) types that have
/// a dimension() and operator[] such as primal::Point and primal::Vector
template <typename T>
class AlmostEqArrMatcher
{
public:
  using is_gtest_matcher = void;

  explicit AlmostEqArrMatcher(const T& arr) : m_arr(arr) { }

  bool MatchAndExplain(const T& other, std::ostream*) const
  {
    for(int i = 0; i < m_arr.dimension(); ++i)
    {
      if(!::testing::Matches(::testing::DoubleNear(m_arr[i], tolerance))(other[i]))
      {
        return false;
      }
    }
    return true;
  }

  void DescribeTo(std::ostream*) const { }
  void DescribeNegationTo(std::ostream*) const { }

  friend bool operator==(const AlmostEqArrMatcher& lhs, const T& rhs)
  {
    return lhs.MatchAndExplain(rhs, nullptr);
  }

  friend bool operator==(const T& lhs, const AlmostEqArrMatcher& rhs) { return rhs == lhs; }

  explicit operator const T() const { return m_arr; }

private:
  const T m_arr;
};

template <typename T, int DIM>
inline auto AlmostEqVector(const primal::Vector<T, DIM>& vec)
{
  return AlmostEqArrMatcher<primal::Vector<T, DIM>>(vec);
}

template <typename T, int DIM>
inline auto AlmostEqPoint(const primal::Point<T, DIM>& pt)
{
  return AlmostEqArrMatcher<primal::Point<T, DIM>>(pt);
}

// ----------------------------------------------------------------------------

/// Create a GMock Matcher for a klee::SliceOperator
class AlmostEqSliceMatcher
{
public:
  using is_gtest_matcher = void;

  explicit AlmostEqSliceMatcher(const klee::SliceOperator& slice) : m_slice(slice) { }

  bool MatchAndExplain(const klee::SliceOperator& other, std::ostream*) const
  {
    return ::testing::Matches(AlmostEqPoint(m_slice.getOrigin()))(other.getOrigin()) &&
      ::testing::Matches(AlmostEqVector(m_slice.getNormal()))(other.getNormal()) &&
      ::testing::Matches(AlmostEqVector(m_slice.getUp()))(other.getUp()) &&
      ::testing::Matches(::testing::Eq(m_slice.getStartProperties()))(other.getStartProperties());
  }

  void DescribeTo(std::ostream*) const { }
  void DescribeNegationTo(std::ostream*) const { }

  friend bool operator==(const AlmostEqSliceMatcher& lhs, const klee::SliceOperator& rhs)
  {
    return lhs.MatchAndExplain(rhs, nullptr);
  }

  friend bool operator==(const klee::SliceOperator& lhs, const AlmostEqSliceMatcher& rhs)
  {
    return rhs == lhs;
  }

private:
  const klee::SliceOperator m_slice;
};

inline auto AlmostEqSlice(const klee::SliceOperator& slice) { return AlmostEqSliceMatcher(slice); }

}  // namespace test
}  // namespace klee
}  // namespace axom

#endif  // AXOM_KLEE_MATCHERS_HPP_
