// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/core/utilities/Sorting.hpp"

#include <algorithm>
#include <random>
#include <numeric>

TEST(core_utilities, insertion_sort_int)
{
  constexpr int NUM_ITERS = 100;
  constexpr int NUM_INTS = 32;
  // Create some RNG state used by std::shuffle
  std::random_device rd;
  std::mt19937 mersenne_twister(rd());

  // Run this for a few iterations to test that sorting works on different
  // random shuffles.
  for(int iter = 0; iter < NUM_ITERS; iter++)
  {
    // Create an array of 32 integers
    int data[NUM_INTS];

    // Fill it with uniform numbers from 0...31
    std::iota(data, data + NUM_INTS, 1);

    // Randomly shuffle the data.
    std::shuffle(data, data + NUM_INTS, mersenne_twister);

    // Sort the data.
    axom::utilities::insertionSort(data, NUM_INTS);

    // The result should be our ordered range 0...31
    for(int i = 0; i < NUM_INTS; i++)
    {
      ASSERT_EQ(data[i], i + 1);
    }
  }
}

TEST(core_utilities, insertion_sort_doubles)
{
  constexpr int NUM_ITERS = 100;
  constexpr int NUM_DBLS = 64;
  // Create some RNG state used by std::shuffle
  std::random_device rd;
  std::mt19937 mersenne_twister(rd());

  // Run this for a few iterations to test that sorting works on different
  // random shuffles.
  for(int iter = 0; iter < NUM_ITERS; iter++)
  {
    // Create an array of 64 doubles
    int data[NUM_DBLS];

    // Fill it with random doubles in the range of [0, 1024)
    std::generate_n(data, NUM_DBLS, []() -> double {
      return axom::utilities::random_real(0.0, 1024.0);
    });

    // Randomly shuffle the data.
    std::shuffle(data, data + NUM_DBLS, mersenne_twister);

    // Sort the data.
    axom::utilities::insertionSort(data, NUM_DBLS);

    // Check that the range is now sorted.
    for(int i = 0; i < NUM_DBLS - 1; i++)
    {
      ASSERT_LE(data[i], data[i + 1]);
    }
  }
}

TEST(core_utilities, qsort_sort_double)
{
  constexpr int NUM_ITERS = 1000;
  constexpr int MAX_SIZE = 100;

  // Run this for a few iterations to test that sorting works on different random shuffles.
  for(int n = 1; n < MAX_SIZE; ++n)
  {
    std::vector<double> vec;
    vec.resize(n);

    for(int iter = 0; iter < NUM_ITERS; iter++)
    {
      // Fill it with random doubles in the range of [0, 1024)
      std::generate_n(vec.data(), n, []() -> double {
        return axom::utilities::random_real(0., 1024.);
      });

      // sort the first two iterations in descending/ascending order
      if(iter == 0)
      {
        std::sort(vec.begin(), vec.end(), [](int a, int b) { return a > b; });
      }
      else if(iter == 1)
      {
        std::sort(vec.begin(), vec.end(), [](int a, int b) { return a < b; });
      }

      axom::utilities::Sorting<double, MAX_SIZE>::qsort(vec.data(), n);

      // Check that the results are sorted
      for(int i = 1; i < n; i++)
      {
        EXPECT_TRUE(vec[i - 1] <= vec[i]) << "data not sorted!";
      }
    }
  }
}

/*!
 * \brief Examine a container and check that all of its elements are increasing order.
 * \param arr The array to examine.
 * \return True if all elements are in increasing order; False otherwise.
 */
template <typename ArrayType>
bool is_increasing(const ArrayType &arr)
{
  bool retval = true;
  for(size_t i = 1; i < arr.size(); i++)
  {
    retval &= (arr[i] >= arr[i - 1]);
  }
  return retval;
}

/*!
 * \brief Examine a container and check that all of its elements are decreasing order.
 * \param arr The array to examine.
 * \return True if all elements are in decreasing order; False otherwise.
 */
template <typename ArrayType>
bool is_decreasing(const ArrayType &arr)
{
  bool retval = true;
  for(size_t i = 1; i < arr.size(); i++)
  {
    retval &= (arr[i] <= arr[i - 1]);
  }
  return retval;
}

TEST(core_utilities, sort_multiple_values)
{
  constexpr axom::IndexType NUM_ITERS = 500;
  constexpr axom::IndexType MAX_SIZE = 100;

  // Run this for a few iterations to test that sorting works on different random shuffles.
  for(axom::IndexType n = 1; n < MAX_SIZE; ++n)
  {
    std::vector<double> vec;
    vec.resize(n);

    for(axom::IndexType iter = 0; iter < NUM_ITERS; iter++)
    {
      // Fill it with random doubles in the range of [0, 1024)
      std::generate_n(vec.data(), n, []() -> double {
        return axom::utilities::random_real(0., 1024.);
      });

      // sort the first two iterations in descending/ascending order
      if(iter == 0)
      {
        std::sort(vec.begin(), vec.end(), [](int a, int b) { return a > b; });
      }
      else if(iter == 1)
      {
        std::sort(vec.begin(), vec.end(), [](int a, int b) { return a < b; });
      }

      // Make a second array.
      std::vector<int> vec2(n);
      std::iota(vec2.begin(), vec2.end(), 0);

      auto vecCopy(vec);
      auto vec2Copy(vec2);

      // Sort values in vec and vec2 based on vec.
      axom::utilities::sort_multiple(vec.data(), vec2.data(), n);
      EXPECT_TRUE(is_increasing(vec));

      // Sort back based on vec2.
      axom::utilities::sort_multiple(vec2.data(), vec.data(), n);
      EXPECT_TRUE(is_increasing(vec2));

      // Make sure both vectors are the same as the initial vectors.
      EXPECT_EQ(vec, vecCopy);
      EXPECT_EQ(vec2, vec2Copy);

      // Reverse values based on vec.
      axom::utilities::reverse_sort_multiple(vec.data(), vec2.data(), n);
      EXPECT_TRUE(is_decreasing(vec));

      // Reverse back based on vec2.
      axom::utilities::reverse_sort_multiple(vec2.data(), vec.data(), n);
      EXPECT_TRUE(is_decreasing(vec2));
    }
  }
}

TEST(core_utilities, sign_of)
{
  using axom::utilities::sign_of;
  EXPECT_EQ(sign_of(1234), 1);
  EXPECT_EQ(sign_of(-1234), -1);
  EXPECT_EQ(sign_of(0), 0);

  EXPECT_EQ(sign_of(0.1234), 1);
  EXPECT_EQ(sign_of(-0.1234), -1);
  EXPECT_EQ(sign_of(0.0), 0);

  EXPECT_EQ(sign_of(1e-9), 1);
  EXPECT_EQ(sign_of(-1e-9), -1);

  EXPECT_EQ(sign_of(1e-9, 1e-12), 1);
  EXPECT_EQ(sign_of(-1e-9, 1e-12), -1);

  EXPECT_EQ(sign_of(1e-9, 1e-8), 0);
  EXPECT_EQ(sign_of(-1e-9, 1e-8), 0);
}
