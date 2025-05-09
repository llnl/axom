// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_CORE_EXECUTION_RUNTIME_POLICY_HPP_
#define AXOM_CORE_EXECUTION_RUNTIME_POLICY_HPP_

#include "axom/config.hpp" /* for compile time defs. */
#include <string>
#include "axom/fmt/format.h" /* for axom::fmt */

#include <map>

/*!
  @file runtime_policy.hpp

  @brief Define runtime policies symbols for selecting.

  The policies are enums corresponding to
  \a axom::execution_space template parameters.
  The difference is that runtime policies are selected at
  run time while \a axom::execution_space is specialized
  at build time.
  @see execution_space.hpp.

  The possible runtime parameters are
  - @c seq: sequential execution on the host
  - @c omp: OpenMP execution
  - @c cuda: GPU execution via CUDA
  - @c hip: GPU execution via HIP

  The available policies depend on how Axom is configured.
  RAJA is required for using OpenMP, CUDA and HIP.
  UMPIRE is required for using CUDA and HIP.
  Sequential execution on host is always available.

  These macros are defined to indicate available non-sequential
  policies.
  - @c AXOM_RUNTIME_POLICY_USE_OPENMP
  - @c AXOM_RUNTIME_POLICY_USE_CUDA
  - @c AXOM_RUNTIME_POLICY_USE_HIP
*/

// Helper preprocessor defines for using OPENMP, CUDA, and HIP policies.
#if defined(AXOM_USE_RAJA)
  #ifdef AXOM_USE_OPENMP
    #define AXOM_RUNTIME_POLICY_USE_OPENMP
  #endif
  #if defined(__CUDACC__) && defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
    #define AXOM_RUNTIME_POLICY_USE_CUDA
  #endif
  #if defined(__HIPCC__) && defined(AXOM_USE_HIP) && defined(AXOM_USE_UMPIRE)
    #define AXOM_RUNTIME_POLICY_USE_HIP
  #endif
#endif

namespace axom
{
namespace runtime_policy
{
/// Execution policies.  The supported set depends on Axom's configuration.
enum class Policy
{
  seq = 0
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  ,
  omp = 1
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  ,
  cuda = 2
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  ,
  hip = 3
#endif
};

//! @brief Mapping from policy name to policy enum.
// clang-format off
  static const std::map<std::string, Policy> s_nameToPolicy
  {
    {"seq", Policy::seq}
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
    , {"omp", Policy::omp}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
    , {"cuda", Policy::cuda}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
    , {"hip", Policy::hip}
#endif
  };

  //! @brief Mapping from policy enum to policy name.
  static const std::map<Policy, std::string> s_policyToName
  {
    {Policy::seq, "seq"}
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
    , {Policy::omp, "omp"}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
    , {Policy::cuda, "cuda"}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
    , {Policy::hip, "hip"}
#endif
  };
// clang-format on

inline Policy nameToPolicy(const std::string &name) { return s_nameToPolicy.find(name)->second; }

inline std::string policyToName(Policy policy) { return s_policyToName.find(policy)->second; }

/// Utility function to allow formating a Policy
static inline auto format_as(Policy pol) { return axom::fmt::underlying(pol); }

}  // end namespace runtime_policy
}  // end namespace axom

#endif /* AXOM_CORE_EXECUTION_RUNTIME_POLICY_HPP_ */
