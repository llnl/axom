// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Instantiates the 3D MarchingCubes implementation for the CUDA policy.
// Each policy and dimension has its own translation unit to limit compile times.

#include "axom/quest/detail/MarchingCubesSingleDomainPolicy.hpp"

#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)

namespace axom::quest::detail::marching_cubes
{
std::unique_ptr<MarchingCubesSingleDomain::ImplBase>
MarchingCubesSingleDomain::newMarchingCubesCudaImpl(std::integral_constant<int, 3>)
{
  return newMarchingCubesPolicyImpl<3, axom::CUDA_EXEC<256>, axom::CUDA_EXEC<1>>();
}

}  // namespace axom::quest::detail::marching_cubes

#endif  // AXOM_RUNTIME_POLICY_USE_CUDA
