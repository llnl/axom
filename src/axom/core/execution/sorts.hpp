// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_CORE_EXECUTION_SORTS_HPP_
#define AXOM_CORE_EXECUTION_SORTS_HPP_

#include "axom/config.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"

#if defined(AXOM_USE_RAJA)
  #include "RAJA/RAJA.hpp"
#else
  #include "axom/core/utilities/Sorting.hpp"
  #include <algorithm>
  #include <numeric>
  #include <vector>
#endif

namespace axom
{

/*!
 * \brief Sort an array.
 *
 * \tparam ExecSpace The execution space where the sort occurs.
 * \tparam T The type of data to sort.
 *
 * \param input The data array to sort.
 * \param size The number of elements to sort.
 */
template <typename ExecSpace, typename T>
inline void sort(T *input, axom::IndexType size)
{
#if defined(AXOM_USE_RAJA)
  // Sort using RAJA
  using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
  RAJA::sort<loop_policy>(RAJA::make_span(input, size));
#else
  constexpr bool is_serial = std::is_same<ExecSpace, SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);
  axom::utilities::Sorting<T>::sort(input, size);
#endif
}

/*!
 * \brief Sort a container.
 *
 * \tparam ExecSpace The execution space where the sort occurs.
 * \tparam ContiguousMemoryContainer Container type for the data to sort.
 *
 * \param input The container to sort.
 */
template <typename ExecSpace, typename ContiguousMemoryContainer>
inline void sort(ContiguousMemoryContainer &input)
{
  sort<ExecSpace>(input.data(), input.size());
}

/*!
 * \brief Sort a pair of containers using the first container's elements as the
 *        values to sort. The second container is sorted the same way.
 *
 * \tparam ExecSpace The execution space where the sort occurs.
 * \tparam Container1 The container type to sort.
 * \tparam Container2 The second container type to sort.
 *
 * \param input1 The container to sort (used as sorting key values).
 * \param input2 A second container to sort (according to input1's sort order).
 */
template <typename ExecSpace, typename Container1, typename Container2>
inline void sort_pairs(Container1 &input1, Container2 &input2)
{
  assert(input1.size() == input2.size());

#if defined(AXOM_USE_RAJA)
  // Sort using RAJA
  using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
  RAJA::sort_pairs<loop_policy>(RAJA::make_span(input1.data(), input1.size()),
                                RAJA::make_span(input2.data(), input2.size()));

#else
  constexpr bool is_serial = std::is_same<ExecSpace, SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);
  axom::utilities::sort_multiple(input1, input2);
#endif
}

/*!
 * \brief Sort a pair of arrays using the first array's elements as the
 *        values to sort. The second array is sorted the same way. This sort
 *        is stable.
 *
 * \tparam ExecSpace The execution space where the sort occurs.
 * \tparam T Type for the first data array to sort.
 * \tparam U Type for the second data array to sort.
 *
 * \param input1 The data array to sort (used as sorting key values).
 * \param input2 A second array to sort (according to input1's sort order).
 * \param size The number of elements in input1 and input2.
 */
template <typename ExecSpace, typename T, typename U>
inline void stable_sort_pairs(T *input1, U *input2, axom::IndexType size)
{
#if defined(AXOM_USE_RAJA)
  // Sort using RAJA
  using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
  RAJA::stable_sort_pairs<loop_policy>(RAJA::make_span(input1, size), RAJA::make_span(input2, size));

#else
  constexpr bool is_serial = std::is_same<ExecSpace, SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);

  // Do stable sort of indices using input1 as the sort key.
  std::vector<axom::IndexType> indices(size);
  std::iota(indices.begin(), indices.end(), 0);
  std::stable_sort(indices.begin(), indices.end(), [&](axom::IndexType index1, axom::IndexType index2) {
    return input1[index1] < input1[index2];
  });

  // Store the values back into the input containers in sort order.
  std::vector<T> input1_copy(input1, input1 + size);
  std::vector<U> input2_copy(input2, input2 + size);
  for(axom::IndexType i = 0; i < size; i++)
  {
    input1[i] = input1_copy[indices[i]];
    input2[i] = input2_copy[indices[i]];
  }
#endif
}

/*!
 * \brief Sort a pair of containers using the first container's elements as the
 *        values to sort. The second container is sorted the same way. This sort
 *        is stable.
 *
 * \tparam ExecSpace The execution space where the sort occurs.
 * \tparam Container1 The container type to sort.
 * \tparam Container2 The second container type to sort.
 *
 * \param input1 The data container to sort (used as sorting key values).
 * \param input2 A second container to sort (according to input1's sort order).
 * \param size The number of elements in input1 and input2.
 */
template <typename ExecSpace, typename Container1, typename Container2>
inline void stable_sort_pairs(Container1 &input1, Container2 &input2)
{
  assert(input1.size() == input2.size());
  stable_sort_pairs<ExecSpace>(input1.data(), input2.data(), input1.size());
}

}  // namespace axom
#endif
