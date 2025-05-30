// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_CORE_EXECUTION_FOR_ALL_HPP_
#define AXOM_CORE_EXECUTION_FOR_ALL_HPP_

#include "axom/config.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/nested_for_exec.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"
#include "axom/core/StackArray.hpp"

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
inline void for_all(const IndexType &begin, const IndexType &end, KernelType &&kernel) noexcept
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
inline void for_all(const IndexType &N, KernelType &&kernel) noexcept
{
  AXOM_STATIC_ASSERT(execution_space<ExecSpace>::valid());
  for_all<ExecSpace>(0, N, std::forward<KernelType>(kernel));
}

/*!
 * \brief Loops over a 2D specified contiguous range.
 *
 * \param [in] iRange The start/end values for the inner loop.
 * \param [in] jRange The start/end values for the outer loop.
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
 *    // compute C[ idx ] for all entries in i:[0-99],j:[0-299]
 *    IndexType NX = 100;
 *    IndexType NY = 300;
 *    StackArray<IndexType, 2> iRange{{0, NX}};
 *    StackArray<IndexType, 2> jRange{{0, NY}};
 *    axom::for_all< axom::OMP_EXEC >( iRange, jRange, AXOM_LAMBDA( IndexType i, IndexType j ) {
 *      const auto idx = j * NX + i;
 *      C[ idx ] = A[ idx ] + B[ idx ];
 *    } );
 *
 * \endcode
 *
 */
template <typename ExecSpace, typename KernelType>
inline void for_all(const axom::StackArray<IndexType, 2> &iRange,
                    const axom::StackArray<IndexType, 2> &jRange,
                    KernelType &&kernel) noexcept
{
  AXOM_STATIC_ASSERT(execution_space<ExecSpace>::valid());
  assert(iRange[1] >= iRange[0] && jRange[1] >= jRange[0]);

#if defined(AXOM_USE_RAJA)
  RAJA::RangeSegment jRangeSeg(jRange[0], jRange[1]);
  RAJA::RangeSegment iRangeSeg(iRange[0], iRange[1]);
  using EXEC_POL = typename axom::internal::nested_for_exec<ExecSpace>::loop2d_policy;
  RAJA::kernel<EXEC_POL>(RAJA::make_tuple(iRangeSeg, jRangeSeg), std::forward<KernelType>(kernel));
#else
  for(IndexType j = jRange[0]; j < jRange[1]; j++)
  {
    for(IndexType i = iRange[0]; i < iRange[1]; i++)
    {
      kernel(i, j);
    }
  }
#endif
}

/*!
 * \brief Loops over a 2D specified contiguous range.
 *
 * \param [in] shape 2 values that indicate the range of the inner, outer loops, respectively.
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
 *    // compute C[ idx ] for all entries in i:[0-99],j:[0-299]
 *    IndexType NX = 100;
 *    IndexType NY = 300;
 *    StackArray<IndexType, 2> shape{{NX, NY}};
 *    axom::for_all< axom::OMP_EXEC >( shape, AXOM_LAMBDA( IndexType i, IndexType j ) {
 *      const auto idx = j * NX + i;
 *      C[ idx ] = A[ idx ] + B[ idx ];
 *    } );
 *
 * \endcode
 *
 */

template <typename ExecSpace, typename KernelType>
inline void for_all(const axom::StackArray<IndexType, 2> &shape, KernelType &&kernel) noexcept
{
  for_all<ExecSpace>(axom::StackArray<IndexType, 2> {{0, shape[0]}},
                     axom::StackArray<IndexType, 2> {{0, shape[1]}},
                     std::forward<KernelType>(kernel));
}

/*!
 * \brief Loops over a 3D specified contiguous range.
 *
 * \param [in] iRange The start/end values for the inner loop.
 * \param [in] jRange The start/end values for the middle loop.
 * \param [in] kRange The start/end values for the outer loop.
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
 *    // compute C[ idx ] for all entries in i:[0-99],j:[0-299],k:[0-49]
 *    IndexType NX = 100;
 *    IndexType NY = 300;
 *    IndexType NZ = 50;
 *    StackArray<IndexType, 2> iRange{{0, NX}};
 *    StackArray<IndexType, 2> kRange{{0, NY}};
 *    StackArray<IndexType, 2> kRange{{0, NZ}};
 *    axom::for_all< axom::OMP_EXEC >(iRange, jRange, kRange, AXOM_LAMBDA( IndexType i, IndexType j, IndexType k ) {
 *      const auto idx = (k * NX * NY) + (j * NX) + i;
 *      C[ idx ] = A[ idx ] + B[ idx ];
 *    } );
 *
 * \endcode
 *
 */
template <typename ExecSpace, typename KernelType>
inline void for_all(const axom::StackArray<IndexType, 2> &iRange,
                    const axom::StackArray<IndexType, 2> &jRange,
                    const axom::StackArray<IndexType, 2> &kRange,
                    KernelType &&kernel) noexcept
{
  AXOM_STATIC_ASSERT(execution_space<ExecSpace>::valid());
  assert(iRange[1] >= iRange[0] && jRange[1] >= jRange[0] && kRange[1] >= kRange[0]);

#if defined(AXOM_USE_RAJA)
  RAJA::RangeSegment iR(iRange[0], iRange[1]);
  RAJA::RangeSegment jR(jRange[0], jRange[1]);
  RAJA::RangeSegment kR(kRange[0], kRange[1]);
  using EXEC_POL = typename axom::internal::nested_for_exec<ExecSpace>::loop3d_policy;
  RAJA::kernel<EXEC_POL>(RAJA::make_tuple(iR, jR, kR), std::forward<KernelType>(kernel));
#else
  for(IndexType k = kRange[0]; k < kRange[1]; k++)
  {
    for(IndexType j = jRange[0]; j < jRange[1]; j++)
    {
      for(IndexType i = iRange[0]; i < iRange[1]; i++)
      {
        kernel(i, j, k);
      }
    }
  }
#endif
}

/*!
 * \brief Loops over a 3D specified contiguous range.
 *
 * \param [in] shape An array containing the x,y,z sizes of the loops.
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
 *    // compute C[ idx ] for all entries in i:[0-99],j:[0-299],k:[0-49]
 *    IndexType NX = 100;
 *    IndexType NY = 300;
 *    IndexType NZ = 50;
 *    StackArray<IndexType, 3> shape{{NX, NY, NZ}};
 *    axom::for_all< axom::OMP_EXEC >(shape, AXOM_LAMBDA( IndexType i, IndexType j, IndexType k ) {
 *      const auto idx = (k * NX * NY) + (j * NX) + i;
 *      C[ idx ] = A[ idx ] + B[ idx ];
 *    } );
 *
 * \endcode
 *
 */
template <typename ExecSpace, typename KernelType>
inline void for_all(const StackArray<IndexType, 3> &shape, KernelType &&kernel) noexcept
{
  for_all<ExecSpace>(axom::StackArray<IndexType, 2> {{0, shape[0]}},
                     axom::StackArray<IndexType, 2> {{0, shape[1]}},
                     axom::StackArray<IndexType, 2> {{0, shape[2]}},
                     std::forward<KernelType>(kernel));
}

/// @}

}  // namespace axom

#endif  // AXOM_CORE_EXECUTION_FOR_ALL_HPP_
