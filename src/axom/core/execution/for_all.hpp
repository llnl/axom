// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_CORE_EXECUTION_FOR_ALL_HPP_
#define AXOM_CORE_EXECUTION_FOR_ALL_HPP_

#include "axom/config.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"

// C/C++ includes
#include <type_traits>
#include <utility>

namespace axom
{
/// \name Generic Loop Traversal Functions
/// @{

/*!
 * \brief Loops over a specified contiguous range, I:[begin,end-1].
 *
 * \param [in] begin start index of the iteration.
 * \param [in] end length of the iteration space.
 * \param [in] kernel user-supplied kernel, i.e., a lambda or functor.
 *
 * \tparam ExecSpace the execution space where to run the supplied kernel
 * \tparam KernelType the type of the supplied kernel (detected by the compiler)
 *
 * \see axom::execution_space
 *
 * Usage Example:
 * \code
 *
 *    double* A = ...
 *    double* B = ...
 *    double* C = ...
 *
 *    // compute C[ idx ] for all entries in [100-499]
 *    axom::for_all< axom::OMP_EXEC >( 100, 500, AXOM_LAMBDA( IndexType idx ) {
 *      C[ idx ] = A[ idx ] + B[ idx ];
 *    } );
 *
 * \endcode
 *
 */
template <typename ExecSpace, typename KernelType>
inline void for_all(const IndexType& begin, const IndexType& end, KernelType&& kernel) noexcept
{
  AXOM_STATIC_ASSERT(execution_space<ExecSpace>::valid());

#ifdef AXOM_USE_RAJA

  using loop_exec = typename execution_space<ExecSpace>::loop_policy;
  RAJA::forall<loop_exec>(RAJA::RangeSegment(begin, end), std::forward<KernelType>(kernel));

#else

  constexpr bool is_serial = std::is_same<ExecSpace, SEQ_EXEC>::value;
  AXOM_STATIC_ASSERT(is_serial);
  for(IndexType i = begin; i < end; ++i)
  {
    kernel(i);
  }

#endif
}

/*!
 * \brief Loops over the contiguous range, I:[0,N-1], given by its length, N.
 *
 * \param [in] N the length of the contiguous range.
 * \param [in] kernel user-supplied kernel, i.e., a lambda or functor.
 *
 * \tparam ExecSpace the execution space where to run the supplied kernel
 * \tparam KernelType the type of the supplied kernel (detected by the compiler)
 *
 * \see axom::execution_space
 *
 * Usage Example:
 * \code
 *
 *    double* A = ...
 *    double* B = ...
 *    double* C = ...
 *
 *    axom::for_all< axom::OMP_EXEC >( 500, AXOM_LAMBDA( IndexType idx ) {
 *      C[ idx ] = A[ idx ] + B[ idx ];
 *    } );
 *
 * \endcode
 *
 */
template <typename ExecSpace, typename KernelType>
inline void for_all(const IndexType& N, KernelType&& kernel) noexcept
{
  AXOM_STATIC_ASSERT(execution_space<ExecSpace>::valid());
  for_all<ExecSpace>(0, N, std::forward<KernelType>(kernel));
}

/// @}

}  // namespace axom

#endif  // AXOM_CORE_EXECUTION_FOR_ALL_HPP_
