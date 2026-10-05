// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Instantiates the 2D MarchingCubes implementation for the HIP policy.
// Each policy and dimension has its own translation unit to limit compile times.

#include "axom/quest/detail/MarchingCubesSingleDomainPolicy.hpp"

#if defined(AXOM_RUNTIME_POLICY_USE_HIP)

namespace axom::quest::detail::marching_cubes
{
std::unique_ptr<MarchingCubesSingleDomain::ImplBase> MarchingCubesSingleDomain::newMarchingCubesHipImpl(
  std::integral_constant<int, 2>)
{
  return newMarchingCubesPolicyImpl<2, axom::HIP_EXEC<256>, axom::HIP_EXEC<1>>();
}

}  // namespace axom::quest::detail::marching_cubes

#endif  // AXOM_RUNTIME_POLICY_USE_HIP
