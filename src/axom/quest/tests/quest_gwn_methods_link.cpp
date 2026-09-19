// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/mint.hpp"
#include "axom/quest/GWNMethods.hpp"

#include "gtest/gtest.h"

#include <cmath>

// These functions are defined in a separate translation unit that also includes
// GWNMethods.hpp. The helpers must link and have one address across both files.
auto setup_gwn_mesh_address() -> decltype(&axom::quest::setup_gwn_mesh);
auto compute_field_stats_address() -> decltype(&axom::quest::compute_field_stats);
auto compute_integrals_address() -> decltype(&axom::quest::compute_integrals);

TEST(quest_gwn_methods_link, shared_function_addresses)
{
  EXPECT_EQ(&axom::quest::setup_gwn_mesh, setup_gwn_mesh_address());
  EXPECT_EQ(&axom::quest::compute_field_stats, compute_field_stats_address());
  EXPECT_EQ(&axom::quest::compute_integrals, compute_integrals_address());
}

TEST(quest_gwn_methods_link, helper_results_on_cartesian_mesh)
{
  mfem::DataCollection dc("gwn_link_test");
  auto* mesh =
    new mfem::Mesh(mfem::Mesh::MakeCartesian2D(1, 1, mfem::Element::QUADRILATERAL, true, 2.0, 3.0));
  setup_gwn_mesh_address()(dc, mesh, 1);

  auto* winding = dc.GetField("winding");
  ASSERT_NE(winding, nullptr);
  ASSERT_NE(dc.GetField("inout"), nullptr);
  *winding = 2.0;

  const auto field_stats = compute_field_stats_address()(*winding);
  EXPECT_DOUBLE_EQ(field_stats.min, 2.0);
  EXPECT_DOUBLE_EQ(field_stats.max, 2.0);
  EXPECT_NEAR(field_stats.l2, std::sqrt(24.0), 1e-12);

  const auto integrals = compute_integrals_address()(*winding);
  EXPECT_NEAR(integrals.integral, 12.0, 1e-12);
  EXPECT_NEAR(integrals.domain_volume, 6.0, 1e-12);
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;
  return RUN_ALL_TESTS();
}
