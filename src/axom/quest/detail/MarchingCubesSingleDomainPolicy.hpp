// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/*!
 * \file MarchingCubesSingleDomainPolicy.hpp
 *
 * \brief Defines the factory that creates a MarchingCubes implementation
 *        for one execution policy and dimension.
 *
 * Only the MarchingCubesSingleDomain<Policy><DIM>D.cpp sources include this header,
 * so each implementation is instantiated in exactly one translation unit.
 * The header is internal to the quest library and is not installed.
 */

#include "axom/config.hpp"

#ifndef AXOM_USE_CONDUIT
  #error "MarchingCubesSingleDomainPolicy.hpp requires conduit"
#endif

#include "axom/core/execution/execution_space.hpp"
#include "axom/quest/detail/MarchingCubesSingleDomain.hpp"
#include "axom/quest/detail/MarchingCubesImpl.hpp"
#if defined(AXOM_USE_BUMP)
  #include "axom/quest/detail/MarchingCubesBumpImpl.hpp"
#endif

#include <memory>
#include <type_traits>

namespace axom::quest::detail::marching_cubes
{
template <int DIM, typename ExecSpace, typename SequentialExecSpace>
std::unique_ptr<MarchingCubesSingleDomain::ImplBase> MarchingCubesSingleDomain::newMarchingCubesPolicyImpl()
{
  static_assert(DIM == 2 || DIM == 3, "MarchingCubes supports only 2D and 3D meshes");

#if defined(AXOM_USE_BUMP)
  if(m_mc.m_useBumpBackend)
  {
  #if defined(_WIN32)
    if constexpr(std::is_same<ExecSpace, axom::SEQ_EXEC>::value)
    {
      return std::make_unique<MarchingCubesBumpImpl<DIM, ExecSpace>>(m_mc.m_allocatorID);
    }
    else
    {
      SLIC_ERROR(
        "MarchingCubes bump backend is not enabled for this runtime policy "
        "on Windows shared-library builds.");
      // With a non-aborting error handler for SLIC_ERROR, we could
      // fall through to the common return below and silently hand back the
      // structured-only legacy kernel for what may be an unstructured mesh.
      return nullptr;
    }
  #else
    return std::make_unique<MarchingCubesBumpImpl<DIM, ExecSpace>>(m_mc.m_allocatorID);
  #endif
  }
#endif

  return std::make_unique<MarchingCubesImpl<DIM, ExecSpace, SequentialExecSpace>>(
    m_mc.m_allocatorID,
    m_mc.m_caseIdsFlat,
    m_mc.m_crossingFlags,
    m_mc.m_scannedFlags,
    m_mc.m_facetIncrs);
}

}  // namespace axom::quest::detail::marching_cubes
