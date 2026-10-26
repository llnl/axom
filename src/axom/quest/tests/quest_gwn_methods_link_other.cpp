// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/mint.hpp"
#include "axom/quest/GWNMethods.hpp"

auto setup_gwn_mesh_address() -> decltype(&axom::quest::setup_gwn_mesh)
{
  return &axom::quest::setup_gwn_mesh;
}

auto compute_field_stats_address() -> decltype(&axom::quest::compute_field_stats)
{
  return &axom::quest::compute_field_stats;
}

auto compute_integrals_address() -> decltype(&axom::quest::compute_integrals)
{
  return &axom::quest::compute_integrals;
}
