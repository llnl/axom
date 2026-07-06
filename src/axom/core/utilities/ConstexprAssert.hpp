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
  #include <type_traits>
#endif

// Detect if we are in a constant-evaluation context:
//   - std::is_constant_evaluated() is a C++20 library function.
//   - __builtin_is_constant_evaluated() is a compiler builtin available in C++17
//   - otherwise, AXOM_DETAIL_HAS_IS_CONSTANT_EVALUATED is 0
#if defined(__cpp_lib_is_constant_evaluated)
  #define AXOM_DETAIL_IS_CONSTANT_EVALUATED() (std::is_constant_evaluated())
  #define AXOM_DETAIL_HAS_IS_CONSTANT_EVALUATED 1
#elif defined(__has_builtin)
  #if __has_builtin(__builtin_is_constant_evaluated)
    #define AXOM_DETAIL_IS_CONSTANT_EVALUATED() (__builtin_is_constant_evaluated())
    #define AXOM_DETAIL_HAS_IS_CONSTANT_EVALUATED 1
  #endif
#elif defined(_MSC_VER) && _MSC_VER >= 1926
  // MSVC provides the builtin but historically did not implement __has_builtin.
  #define AXOM_DETAIL_IS_CONSTANT_EVALUATED() (__builtin_is_constant_evaluated())
  #define AXOM_DETAIL_HAS_IS_CONSTANT_EVALUATED 1
#endif

#if !defined(AXOM_DETAIL_HAS_IS_CONSTANT_EVALUATED)
  #define AXOM_DETAIL_HAS_IS_CONSTANT_EVALUATED 0
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
  #if AXOM_DETAIL_HAS_IS_CONSTANT_EVALUATED
    if(AXOM_DETAIL_IS_CONSTANT_EVALUATED())
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
#undef AXOM_DETAIL_IS_CONSTANT_EVALUATED
#undef AXOM_DETAIL_HAS_IS_CONSTANT_EVALUATED

#endif  // AXOM_CORE_UTILITIES_CONSTEXPR_ASSERT_HPP_
