// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "gtest/gtest.h"

#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/core/ItemCollection.hpp"

#include <algorithm>
#include <complex>
#include <functional>
#include <iterator>
#include <numeric>
#include <type_traits>
#include <utility>

namespace
{
// Detects whether subscripting a const iterator is a valid expression.
template <typename Iter, typename = void>
struct has_subscript : std::false_type
{ };

template <typename Iter>
struct has_subscript<Iter,
                     std::void_t<decltype(std::declval<const Iter&>()[std::declval<
                       typename std::iterator_traits<Iter>::difference_type>()])>> : std::true_type
{ };
}  // namespace

//------------------------------------------------------------------------------
TEST(core_array_iterator, random_access_contract_is_static)
{
  using ArrayIter = axom::Array<double>::ArrayIterator;
  using ConstArrayIter = axom::Array<double>::ConstArrayIterator;
  using ViewIter = decltype(std::declval<axom::ArrayView<double>&>().begin());

  static_assert(std::is_same<typename std::iterator_traits<ArrayIter>::iterator_category,
                             std::random_access_iterator_tag>::value,
                "Array's iterator advertises random access");
  static_assert(std::is_same<typename std::iterator_traits<ConstArrayIter>::iterator_category,
                             std::random_access_iterator_tag>::value,
                "Array's const iterator advertises random access");
  static_assert(std::is_same<typename std::iterator_traits<ViewIter>::iterator_category,
                             std::random_access_iterator_tag>::value,
                "ArrayView's iterator advertises random access");

  static_assert(has_subscript<ArrayIter>::value, "Array iterator needs i[n]");
  static_assert(has_subscript<ConstArrayIter>::value, "Array const_iterator needs i[n]");
  static_assert(has_subscript<ViewIter>::value, "ArrayView iterator needs i[n]");

  static_assert(std::is_same<decltype(std::declval<const ArrayIter&>()[0]), double&>::value,
                "Array iterator subscripting must preserve mutability");
  static_assert(std::is_same<decltype(std::declval<const ConstArrayIter&>()[0]), const double&>::value,
                "Array const_iterator subscripting must return a const reference");

  // Note: The subscript is defined on ArrayIteratorBase, not on the common IteratorBase
  // used by derived forward-only iterators (like ItemCollection)
  using ForwardIter = axom::ItemCollection<double>::iterator;
  static_assert(!has_subscript<ForwardIter>::value,
                "ItemCollection's forward iterator must not acquire random access");

  SUCCEED();
}

//------------------------------------------------------------------------------
TEST(core_array_iterator, subscript_matches_offset_dereference)
{
  axom::Array<int> arr(5);
  std::iota(arr.begin(), arr.end(), 10);

  auto it = arr.begin();
  for(axom::IndexType n = 0; n < arr.size(); ++n)
  {
    EXPECT_EQ(it[n], *(it + n));
    EXPECT_EQ(it[n], arr[n]);
  }

  // Subscript is relative to the iterator's position, not to begin().
  auto mid = arr.begin() + 2;
  EXPECT_EQ(mid[0], arr[2]);
  EXPECT_EQ(mid[2], arr[4]);
  EXPECT_EQ(mid[-2], arr[0]);
}

//------------------------------------------------------------------------------
TEST(core_array_iterator, subscript_is_a_mutable_reference)
{
  axom::Array<int> arr(3);
  arr.fill(0);

  const auto it = arr.begin();
  it[1] = 42;
  EXPECT_EQ(arr[1], 42);
}

//------------------------------------------------------------------------------
// Regression: libc++ 22 rewrote std::__sift_down to index the iterator rather
// than dereference it, so std::sort over an axom::Array failed to compile.
// See numerics::solve_polynomial_durand_kerner, which sorts an
// axom::Array<std::complex<double>>.
//------------------------------------------------------------------------------
TEST(core_array_iterator, sort_over_array_of_complex)
{
  using Complex = std::complex<double>;

  axom::Array<Complex> roots;
  roots.push_back(Complex {3.0, 0.0});
  roots.push_back(Complex {1.0, 2.0});
  roots.push_back(Complex {1.0, -2.0});
  roots.push_back(Complex {2.0, 0.0});

  std::sort(roots.begin(), roots.end(), [](const Complex& lhs, const Complex& rhs) {
    if(lhs.real() != rhs.real())
    {
      return lhs.real() < rhs.real();
    }
    return lhs.imag() < rhs.imag();
  });

  EXPECT_EQ(roots[0], Complex(1.0, -2.0));
  EXPECT_EQ(roots[1], Complex(1.0, 2.0));
  EXPECT_EQ(roots[2], Complex(2.0, 0.0));
  EXPECT_EQ(roots[3], Complex(3.0, 0.0));
}

//------------------------------------------------------------------------------
TEST(core_array_iterator, heap_algorithms_over_array_view)
{
  axom::Array<int> arr(6);
  const int values[6] = {5, 1, 4, 2, 6, 3};
  for(axom::IndexType i = 0; i < arr.size(); ++i)
  {
    arr[i] = values[i];
  }

  axom::ArrayView<int> view(arr);

  // make_heap/sort_heap route through __sift_down, which requires the subscript operator.
  std::make_heap(view.begin(), view.end());
  EXPECT_TRUE(std::is_heap(view.begin(), view.end()));

  std::sort_heap(view.begin(), view.end());
  EXPECT_TRUE(std::is_sorted(view.begin(), view.end()));
  EXPECT_EQ(view[0], 1);
  EXPECT_EQ(view[5], 6);

  // partial_sort also reaches __sift_down.
  std::partial_sort(view.begin(), view.begin() + 3, view.end(), std::greater<int> {});
  EXPECT_EQ(view[0], 6);
  EXPECT_EQ(view[1], 5);
  EXPECT_EQ(view[2], 4);
}
