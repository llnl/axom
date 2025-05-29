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

#if defined(AXOM_USE_RAJA)
  #include "RAJA/RAJA.hpp"
#endif

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
 * \tparam Container1 The container type that holds the input data.
 * \tparam Container2 The container type that holds the output data.
 *
 * \note  The container types Container1 and Container2 work with Axom's Array,
 *        ArrayView types as well as std::vector or any other type that provides
 *        a similar interface with these characteristics: 1) size() method,
 *        2) data() method to return container memory pointer, 3) operator[] to
 *        return element access, 4) a "value_type" sub-type that indicates the
 *        type of data stored in the container.
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
template <typename ExecSpace, typename Container1, typename Container2>
inline void exclusive_scan(const Container1 &input, Container2 &&output)
{
  assert(input.size() == output.size());

#if defined(AXOM_USE_RAJA)
  #if defined(AXOM_USE_OPENMP) && defined(__INTEL_LLVM_COMPILER)
  // NOTE: This workaround was brought to this central location instead of
  //       replicating throughout Axom.
  // Intel oneAPI compiler workaround for OpenMP RAJA scan
  using exec_space =
    typename std::conditional<std::is_same<ExecSpace, axom::OMP_EXEC>::value, axom::SEQ_EXEC, ExecSpace>::type;
  #else
  using exec_space = ExecSpace;
  #endif
  using loop_policy = typename axom::execution_space<exec_space>::loop_policy;
  RAJA::exclusive_scan<loop_policy>(RAJA::make_span(input.data(), input.size()),
                                    RAJA::make_span(output.data(), output.size()));

#else
  constexpr bool is_serial = std::is_same<ExecSpace, SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);

  typename std::remove_const<typename Container1::value_type>::type total {0};
  for(IndexType i = 0; i < input.size(); ++i)
  {
    output[i] = total;
    total += input[i];
  }
#endif
}

/*!
 * \brief Performs exclusive scan over \a input view and stores result also in \a input.
 *
 * \param [inout] input The container to be scanned.
 *
 * \tparam ExecSpace the execution space where to run the supplied kernel
 * \tparam Container The container type that holds the data
 *
 * \note  The Container type works with Axom's Array,
 *        ArrayView types as well as std::vector or any other type that provides
 *        a similar interface with these characteristics: 1) size() method,
 *        2) data() method to return container memory pointer, 3) operator[] to
 *        return element access, 4) a "value_type" sub-type that indicates the
 *        type of data stored in the container.
 */
template <typename ExecSpace, typename Container>
inline void exclusive_scan_inplace(Container &&input)
{
#if defined(AXOM_USE_RAJA)
  using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
  RAJA::exclusive_scan_inplace<loop_policy>(RAJA::make_span(input.data(), input.size()));
#else
  constexpr bool is_serial = std::is_same<ExecSpace, SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);

  typename std::remove_const<typename Container::value_type>::type total {0};
  for(IndexType i = 0; i < input.size(); ++i)
  {
    const auto tmp = input[i];
    input[i] = total;
    total += tmp;
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
 * \tparam Container1 The container type that holds the input data
 * \tparam Container2 The container type that holds the output data
 *
 * \note  The container types Container1 and Container2 work with Axom's Array,
 *        ArrayView types as well as std::vector or any other type that provides
 *        a similar interface with these characteristics: 1) size() method,
 *        2) data() method to return container memory pointer, 3) operator[] to
 *        return element access, 4) a "value_type" sub-type that indicates the
 *        type of data stored in the container.
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
template <typename ExecSpace, typename Container1, typename Container2>
inline void inclusive_scan(const Container1 &input, Container2 &&output)
{
  assert(input.size() == output.size());

#if defined(AXOM_USE_RAJA)
  using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
  RAJA::inclusive_scan<loop_policy>(RAJA::make_span(input.data(), input.size()),
                                    RAJA::make_span(output.data(), output.size()));

#else
  constexpr bool is_serial = std::is_same<ExecSpace, SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);

  typename std::remove_const<typename Container1::value_type>::type total {0};
  for(IndexType i = 0; i < input.size(); ++i)
  {
    total += input[i];
    output[i] = total;
  }
#endif
}

/*!
 * \brief Performs inclusive scan over \a input view and stores result in \a input.
 *
 * \param [inout] input The container to be scanned.
 *
 * \tparam ExecSpace the execution space where to run the supplied kernel
 * \tparam Container The container type that holds the data
 *
 * \note  The Container type works with Axom's Array,
 *        ArrayView types as well as std::vector or any other type that provides
 *        a similar interface with these characteristics: 1) size() method,
 *        2) data() method to return container memory pointer, 3) operator[] to
 *        return element access, 4) a "value_type" sub-type that indicates the
 *        type of data stored in the container.
 */
template <typename ExecSpace, typename Container>
inline void inclusive_scan_inplace(Container &&input)
{
#if defined(AXOM_USE_RAJA)
  using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
  RAJA::inclusive_scan_inplace<loop_policy>(RAJA::make_span(input.data(), input.size()));
#else
  constexpr bool is_serial = std::is_same<ExecSpace, SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);

  typename std::remove_const<typename Container::value_type>::type total {0};
  for(IndexType i = 0; i < input.size(); ++i)
  {
    total += input[i];
    input[i] = total;
  }
#endif
}

/// @}

}  // namespace axom

#endif  // AXOM_CORE_EXECUTION_FOR_ALL_HPP_
