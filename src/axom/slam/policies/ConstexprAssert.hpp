// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file ConstexprAssert.hpp
 *
 * \brief A constexpr-friendly assertion for slam's compile-time arithmetic core.
 *
 * slam's debug assertion (SLIC_ASSERT) expands, on the host debug path, 
 * to a block containing a `std::ostringstream` (a non-literal type)
 * whih is ill-formed before C++23 -- so SLIC_ASSERT cannot appear in a function
 * that must also be usable in a constant expression. 
 *
 * SLAM_CONSTEXPR_ASSERT bridges the gap:
 *   - When the condition holds, it is a true no-op and is fully usable inside a
 *     constant expression.
 *   - At run time with AXOM_DEBUG enabled, a violated condition logs via SLIC
 *     and aborts, matching SLIC_ASSERT's runtime behavior.
 *   - During constant evaluation, a violated condition calls a non-constexpr
 *     function, which is a hard compile error pointing at the offending constant
 *   - When AXOM_DEBUG is off, it collapses to a no-op.
 */

#ifndef SLAM_POLICIES_CONSTEXPR_ASSERT_H_
#define SLAM_POLICIES_CONSTEXPR_ASSERT_H_

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"

#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)
  #include "axom/slic/interface/slic_macros.hpp"
  #include "axom/core/utilities/Utilities.hpp"
#endif

namespace axom::slam::detail
{
#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)

/*!
 * \brief Runtime handler for a failed constexpr assert (host debug only).
 *
 * Not marked constexpr: reaching this during constant evaluation is the
 * mechanism by which a compile-time invariant violation becomes a hard error.
 */
[[noreturn]] inline void constexprAssertFail(const char* expr, const char* /*file*/, int /*line*/)
{
  SLIC_ERROR("Failed constexpr assert: " << expr);
  // SLIC_ERROR may not abort if abort-on-error is disabled
  axom::utilities::processAbort();
}

/*!
 * \brief constexpr-safe assertion.
 * \param cond the invariant that must hold
 * \param expr stringized form of \a cond (for diagnostics)
 */
constexpr void constexprAssert(bool cond, const char* expr, const char* file, int line)
{
  cond ? void(0) : constexprAssertFail(expr, file, line);
}

#else  // release, or device code: no-op

constexpr void constexprAssert(bool, const char*, const char*, int) { }

#endif

}  // namespace axom::slam::detail

/*!
 * \def SLAM_CONSTEXPR_ASSERT(EXP)
 * \brief Assert \a EXP in a way that is valid inside constexpr functions.
 */
#define SLAM_CONSTEXPR_ASSERT(EXP) \
  ::axom::slam::detail::constexprAssert((EXP), #EXP, __FILE__, __LINE__)

#endif  // SLAM_POLICIES_CONSTEXPR_ASSERT_H_
