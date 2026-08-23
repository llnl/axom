// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/**
 * \file Ranges.hpp
 *
 * \brief Host-only integration between SLAM types and the C++20 ranges library.
 *
 * This header is intentionally excluded from axom/slam.hpp so device-facing translation units
 * do not acquire a dependency on the host standard library's ranges implementation.
 */

#include "axom/slam/RangeSet.hpp"

#include <ranges>

namespace std::ranges
{
/**
 * A RangeSet iterator stores the complete concrete range state by value,
 * so it remains valid when the RangeSet object that produced it is destroyed.
 */
template <typename PositionType, typename ElementType, typename InterfacePolicy>
inline constexpr bool enable_borrowed_range<
  axom::slam::GenericRangeSet<PositionType,
                              ElementType,
                              axom::slam::policies::RuntimeOffset<PositionType>,
                              axom::slam::policies::StrideOne<PositionType>,
                              axom::slam::policies::NoIndirection<PositionType, ElementType>,
                              axom::slam::policies::NoSubset,
                              InterfacePolicy>> = true;
}  // namespace std::ranges
