// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/ArrayView.hpp"
#include "axom/core/utilities/Checksum.hpp"

TEST(core_checksum, scalar_matches_singleton_view_and_scale_factor)
{
  const double value = 3.25;
  axom::ArrayView<const double> singletonView(&value, 1);

  const auto scalarChecksum = axom::utilities::checksum(value);

  EXPECT_EQ(scalarChecksum, axom::utilities::checksum(singletonView));
  EXPECT_EQ(scalarChecksum * 2.5L, axom::utilities::checksum(value, 2.5));
}

TEST(core_checksum, empty_view_is_zero)
{
  axom::ArrayView<const int> emptyView(nullptr, 0);

  EXPECT_EQ(0.0L, axom::utilities::checksum(emptyView));
}

TEST(core_checksum, array_order_and_values_change_checksum)
{
  const double ordered[] = {1.0, 2.0, 3.0, 4.0};
  const double reordered[] = {4.0, 3.0, 2.0, 1.0};
  const double modified[] = {1.0, 2.0, 3.0, 5.0};

  const auto orderedChecksum =
    axom::utilities::checksum(axom::ArrayView<const double>(ordered, 4));

  EXPECT_NE(
    orderedChecksum,
    axom::utilities::checksum(axom::ArrayView<const double>(reordered, 4)));
  EXPECT_NE(
    orderedChecksum,
    axom::utilities::checksum(axom::ArrayView<const double>(modified, 4)));
}
