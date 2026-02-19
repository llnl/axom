// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_CORE_EXECUTION_SYNCHRONIZE_HPP_
#define AXOM_CORE_EXECUTION_SYNCHRONIZE_HPP_

#include "axom/config.hpp"                         /* for compile time defs. */
#include "axom/core/Macros.hpp"                    /* for AXOM_STATIC_ASSERT */
#include "axom/core/execution/execution_space.hpp" /* execution_space traits */

namespace axom
{
// _gpu_synchronize_start
/*!
 * \brief Synchronizes all execution threads when using an ASYNC policy with
 *  the specified execution space.
 *
 * \tparam ExecSpace the execution space
 */
template <typename ExecSpace>
inline void synchronize() noexcept
{
  AXOM_STATIC_ASSERT(execution_space<ExecSpace>::valid());

#ifdef AXOM_USE_RAJA
  using sync_policy = typename execution_space<ExecSpace>::sync_policy;
  RAJA::synchronize<sync_policy>();
#endif
}
// _gpu_synchronize_end

template <>
inline void synchronize<SEQ_EXEC>() noexcept
{ }

}  // namespace axom

#endif /* AXOM_CORE_EXECUTION_SYNCHRONIZE_HPP_ */
