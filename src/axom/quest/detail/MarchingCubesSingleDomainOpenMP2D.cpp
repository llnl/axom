// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Instantiates the 2D MarchingCubes implementation for the OpenMP policy.
// Each policy and dimension has its own translation unit to limit compile times.

#include "axom/quest/detail/MarchingCubesSingleDomainPolicy.hpp"

#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)

namespace axom::quest::detail::marching_cubes
{
std::unique_ptr<MarchingCubesSingleDomain::ImplBase>
MarchingCubesSingleDomain::newMarchingCubesOpenMPImpl(std::integral_constant<int, 2>)
{
  return newMarchingCubesPolicyImpl<2, axom::OMP_EXEC, axom::SEQ_EXEC>();
}

}  // namespace axom::quest::detail::marching_cubes

#endif  // AXOM_RUNTIME_POLICY_USE_OPENMP
