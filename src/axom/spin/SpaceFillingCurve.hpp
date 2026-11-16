// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

namespace axom::spin
{
/// \brief Different space-filling curves for spatial ordering.
enum class SpaceFillingCurve
{
  Morton,
  Hilbert
};
}  // namespace axom::spin
