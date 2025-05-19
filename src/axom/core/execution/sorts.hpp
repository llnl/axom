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
#endif

namespace axom
{

/*!
 * \brief Sort a pair of containers using the first container's elements as the
 *        values to sort. The second container is sorted the same way.
 *
 * \tparam ExecSpace The execution space where the sort occurs.
 * \tparam ContiguousMemoryContainer Container type for the data to sort.
 *
 * \brief input1 The container to sort (used as sorting key values).
 * \brief input2 A second container to sort (according to input1's sort order).
 */
template <typename ExecSpace, typename ContiguousMemoryContainer>
inline void sort_pairs(ContiguousMemoryContainer &input1, ContiguousMemoryContainer &input2)
{
  assert(input1.size() == input2.size());

#if defined(AXOM_USE_RAJA)
  // Sort using RAJA
  using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
  RAJA::sort_pairs<loop_policy>(RAJA::make_span(input1.data(), input1.size()),
                                RAJA::make_span(input2.data(), input2.size()));

#else
  AXOM_STATIC_ASSERT(std::is_same<ExecSpace, SEQ_EXEC>::value);
  axom::utilities::sort_multiple(input1, input2);
#endif
}

/*!
 * \brief Sort a pair of containers using the first container's elements as the
 *        values to sort. The second container is sorted the same way. This sort
 *        is stable.
 *
 * \tparam ExecSpace The execution space where the sort occurs.
 * \tparam ContiguousMemoryContainer Container type for the data to sort.
 *
 * \brief input1 The container to sort (used as sorting key values).
 * \brief input2 A second container to sort (according to input1's sort order).
 */
template <typename ExecSpace, typename ContiguousMemoryContainer>
inline void stable_sort_pairs(ContiguousMemoryContainer &input1, ContiguousMemoryContainer &input2)
{
  assert(input1.size() == input2.size());

#if defined(AXOM_USE_RAJA)
  // Sort using RAJA
  using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
  RAJA::stable_sort_pairs<loop_policy>(RAJA::make_span(input1.data(), input1.size()),
                                       RAJA::make_span(input2.data(), input2.size()));

#else
  AXOM_STATIC_ASSERT(std::is_same<ExecSpace, SEQ_EXEC>::value);

  // Do stable sort of indices using input1 as the sort key.
  std::vector<axom::IndexType> indices(input1.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::stable_sort(indices.begin(), indices.end(), [&](axom::IndexType index1, axom::IndexType index2) {
    return input1[index1] < input1[index2];
  });

  // Store the values back into the input containers in sort order.
  ContiguousMemoryContainer input1_copy(input1), input2_copy(input2);
  const auto n = input1.size();
  for(axom::IndexType i = 0; i < n; i++)
  {
    input1[i] = input1_copy[indices[i]];
    input2[i] = input2_copy[indices[i]];
  }
#endif
}

}  // namespace axom
#endif
