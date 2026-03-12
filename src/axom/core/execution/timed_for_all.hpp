// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_CORE_EXECUTION_TIMED_FOR_ALL_HPP_
#define AXOM_CORE_EXECUTION_TIMED_FOR_ALL_HPP_

#include "axom/config.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/AnnotationMacros.hpp"

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP)
  #include "axom/core/Array.hpp"
  #include <iostream>
  #include <chrono>
  #include <omp.h>
#endif

namespace axom
{
/// \name Generic Timed Loop Traversal Functions
/// @{

namespace detail
{
/*!
 * \brief Default implementation for timing axom::for_all.
 */
template <typename ExecSpace, typename KernelType>
struct TimedForAll
{
  /*!
   * \brief Execute the for_all.
   *
   * \param name The name of the loop being timed.
   * \param n The number it items in the loop.
   * \param kernel The kernel to execute.
   */
  static void execute([[maybe_unused]] const std::string &name, axom::IndexType n, KernelType &&kernel)
  {
    AXOM_ANNOTATE_SCOPE(name);
    axom::for_all<ExecSpace>(n, std::forward<KernelType>(kernel));
  }
};

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP)
/*!
 * \brief A specialization for OpenMP that times each thread and prints the start/end times.
 */
template <typename KernelType>
struct TimedForAll<axom::OMP_EXEC, KernelType>
{
  using ExecSpace = axom::OMP_EXEC;

  /*!
   * \brief Execute the for_all using OpenMP.
   *
   * \param name The name of the loop being timed.
   * \param n The number it items in the loop.
   * \param kernel The kernel to execute.
   */
  static void execute(const std::string &name, axom::IndexType n, KernelType &&kernel)
  {
    AXOM_ANNOTATE_BEGIN(name);
    const auto now1 =
      std::chrono::duration<double>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    const int nthreads = omp_get_max_threads() + 1;  // make an extra slot.
    axom::Array<double> ompStart(nthreads, nthreads, allocatorID);
    axom::Array<double> ompEnd(nthreads, nthreads, allocatorID);
    auto ompStartView = ompStart.view();
    auto ompEndView = ompEnd.view();

    ompStart.fill(-1.);
    // Save the start time as the last element in the array.
    ompStart[nthreads - 1] = now1;

    auto outer = [&](axom::IndexType i) {
      // Save the start time.
      double &start = ompStartView[omp_get_thread_num()];
      if(start < 0.)
      {
        start =
          std::chrono::duration<double>(std::chrono::high_resolution_clock::now().time_since_epoch())
            .count();
      }

      kernel(i);

      // Save the end time.
      ompEndView[omp_get_thread_num()] =
        std::chrono::duration<double>(std::chrono::high_resolution_clock::now().time_since_epoch())
          .count();
    };

    // Run the outer kernel to gather timings and run the kernel.
    axom::for_all<ExecSpace>(n, outer);

    // Save the end time as the last element in the array.
    ompEnd[nthreads - 1] =
      std::chrono::duration<double>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    AXOM_ANNOTATE_END(name);

    std::cout << name << ":\n";
    std::cout << "\tn: " << n << "\n";
    std::cout << "\tnthreads: " << (nthreads - 1) << "\n";
    std::cout << std::setprecision(20) << "\tstart=" << ompStart << std::endl;
    std::cout << "\tend=" << ompEnd << std::endl;
  }
};
#endif
}  // end namespace detail

/*!
 * \brief Execute axom::for_all and add a caliper timer (if enabled). Certain ExecSpace types may output additional timing information.
 *
 * \param name The name of the loop being timed.
 * \param n The number it items in the loop.
 * \param kernel The kernel to execute.
 */
template <typename ExecSpace, typename KernelType>
void timed_for_all(const std::string &name, axom::IndexType n, KernelType &&kernel)
{
  detail::TimedForAll<ExecSpace, KernelType>::execute(name, n, std::forward<KernelType>(kernel));
}

/// @}

}  // namespace axom

#endif  // AXOM_CORE_EXECUTION_TIMED_FOR_ALL_HPP_
