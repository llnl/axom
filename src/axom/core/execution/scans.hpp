// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_CORE_EXECUTION_SCANS_HPP_
#define AXOM_CORE_EXECUTION_SCANS_HPP_

#include "axom/config.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"

// C/C++ includes
#include <type_traits>
#include <utility>

namespace axom
{
/// \name Scans
/// @{

/*!
 * \brief Performs exclusive scan over \a input view and stores result in \a output.
 *
 * \param [in] input The input container to be scanned.
 * \param [out] output The container that will contain the output scan data. This 
 *                     must have the same number of elements as \a input.
 *
 * \tparam ExecSpace the execution space where to run the supplied kernel
 * \tparam ContiguousMemoryContainer The container type that holds the data
 *
 * \see axom::execution_space
 *
 * Usage Example:
 * \code
 *
 *    axom::ArrayView<int> sizesView = sizes.view();
 *    axom::ArrayView<int> offsetsView = offsets.view();
 *
 *    // Compute the scan for all elements in sizesView, store scan in offsetsView.
 *    axom::exclusive_scan<ExecSpace>(sizesView, offsetsView);
 *
 * \endcode
 *
 */
template <typename ExecSpace, typename ContiguousMemoryContainer>
inline void exclusive_scan(const ContiguousMemoryContainer &input, ContiguousMemoryContainer &output)
{
  assert(input.size() == output.size());

#ifdef AXOM_USE_RAJA

  using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
  RAJA::exclusive_scan<loop_policy>(RAJA::make_span(input.data(), input.size()),
                                    RAJA::make_span(output.data(), output.size()));

#else
  constexpr bool is_serial = std::is_same<ExecSpace, SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);

  typename ContiguousMemoryContainer::value_type total {0};
  for(IndexType i = 0; i < input.size(); ++i)
  {
    output[i] = total;
    total += input[i];
  }
#endif
}

/*!
 * \brief Performs inclusive scan over \a input view and stores result in \a output.
 *
 * \param [in] input The input container to be scanned.
 * \param [out] output The container that will contain the output scan data. This 
 *                     must have the same number of elements as \a input.
 *
 * \tparam ExecSpace the execution space where to run the supplied kernel
 * \tparam ContiguousMemoryContainer The container type that holds the data
 *
 * \see axom::execution_space
 *
 * Usage Example:
 * \code
 *
 *    axom::ArrayView<int> sizesView = sizes.view();
 *    axom::ArrayView<int> totalView = totals.view();
 *
 *    // Compute the scan for all elements in sizesView, store scan in totalView.
 *    axom::inclusive_scan<ExecSpace>(sizesView, totalView);
 *
 * \endcode
 *
 */
template <typename ExecSpace, typename ContiguousMemoryContainer>
inline void inclusive_scan(const ContiguousMemoryContainer &input, ContiguousMemoryContainer &output)
{
  assert(input.size() == output.size());

#ifdef AXOM_USE_RAJA

  using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
  RAJA::inclusive_scan<loop_policy>(RAJA::make_span(input.data(), input.size()),
                                    RAJA::make_span(output.data(), output.size()));

#else
  constexpr bool is_serial = std::is_same<ExecSpace, SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);

  typename ContiguousMemoryContainer::value_type total {0};
  for(IndexType i = 0; i < input.size(); ++i)
  {
    total += input[i];
    output[i] = total;
  }
#endif
}
/// @}

}  // namespace axom

#endif  // AXOM_CORE_EXECUTION_FOR_ALL_HPP_
