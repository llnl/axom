// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file ConstexprAssert.hpp
 *
 * \brief A constexpr-friendly assertion for host/device code.
 *
 * This header provides the low-level implementation routine
 * `axom::detail::constexprAssert(bool, const char*, const char*, int)`.
 *
 * Most Axom code should use the convenience macro `AXOM_CONSTEXPR_ASSERT(EXP)`
 * (defined in `axom/core/Macros.hpp`), which forwards expression text and source location
 * to `axom::detail::constexprAssert(...)`.
 *
 * \note Rationale: Some debug assertion mechanisms use non-literal types (e.g. iostreams)
 * and therefore cannot appear inside functions that must be usable in constant expressions under C++17.
 * This facility provides an assertion-like hook that is valid in constexpr-capable code and remains
 * device-compilable across Axom's supported backends.
 *
 * Semantics:
 * - If the condition is false during constant evaluation, compilation fails.
 * - If the condition is false at run time and AXOM_DEBUG is enabled (host code), the process aborts.
 * - In non-debug runtime builds, it is a no-op.
 * - In device compilation, it is a no-op (kernels cannot throw/abort portably).
 */

#ifndef AXOM_CORE_UTILITIES_CONSTEXPR_ASSERT_HPP_
#define AXOM_CORE_UTILITIES_CONSTEXPR_ASSERT_HPP_

#include "axom/config.hpp"
#include "axom/core/utilities/Abort.hpp"

// This header is included by `axom/core/Macros.hpp`, so redefine the necessary macro(s) or spell them out
//
// The following is equivalent to AXOM_HOST_DEVICE
#if defined(__CUDACC__) || defined(__HIPCC__)
  #define AXOM_DETAIL_CONSTEXPR_ASSERT_HOST_DEVICE __host__ __device__
#else
  #define AXOM_DETAIL_CONSTEXPR_ASSERT_HOST_DEVICE
#endif

#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
  #include <cassert>
#endif

namespace axom::detail
{
// This conditional is equivalent to !defined(AXOM_DEVICE_CODE)
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)

[[noreturn]] inline void constexprAssertFail(const char* /*expr*/, const char* /*file*/, int /*line*/)
{
  #if defined(AXOM_DEBUG)
  // Provide a debugger-friendly trap site in debug builds.
  assert(false && "Failed AXOM_CONSTEXPR_ASSERT");
  #endif
  axom::utilities::processAbort();
}

AXOM_DETAIL_CONSTEXPR_ASSERT_HOST_DEVICE constexpr void constexprAssert(bool cond,
                                                                        const char* expr,
                                                                        const char* file,
                                                                        int line)
{
  if(!cond)
  {
  #if defined(__clang__) || defined(__GNUC__)
    if(__builtin_is_constant_evaluated())
    {
      // Not constexpr: reaching this in constant evaluation is a hard error.
      constexprAssertFail(expr, file, line);
    }
  #endif

  #if defined(AXOM_DEBUG)
    constexprAssertFail(expr, file, line);
  #else
    static_cast<void>(expr);
    static_cast<void>(file);
    static_cast<void>(line);
  #endif
  }
}

#else  // device compilation: no-op

AXOM_DETAIL_CONSTEXPR_ASSERT_HOST_DEVICE constexpr void constexprAssert(bool,
                                                                        const char*,
                                                                        const char*,
                                                                        int)
{ }

#endif
}  // namespace axom::detail

#undef AXOM_DETAIL_CONSTEXPR_ASSERT_HOST_DEVICE

#endif  // AXOM_CORE_UTILITIES_CONSTEXPR_ASSERT_HPP_
