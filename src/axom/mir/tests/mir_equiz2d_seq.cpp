// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/bump/tests/blueprint_testing_helpers.hpp"
#include "axom/mir/tests/mir_equiz2d_impl.hpp"

axom::blueprint::testing::TestApplication TestApp;

TEST(mir_equiz2d_seq, miralgorithm)
{
  axom::mir::MIRAlgorithm* m = nullptr;
  EXPECT_EQ(m, nullptr);
}

TEST(mir_equiz2d_seq, materialinformation)
{
  conduit::Node matset;
  matset["material_map/a"] = 1;
  matset["material_map/b"] = 2;
  matset["material_map/c"] = 0;

  auto mi = axom::bump::views::materials(matset);
  EXPECT_EQ(mi.size(), 3);
  EXPECT_EQ(mi[0].m_number, 1);
  EXPECT_EQ(mi[0].m_name, "a");

  EXPECT_EQ(mi[1].m_number, 2);
  EXPECT_EQ(mi[1].m_name, "b");

  EXPECT_EQ(mi[2].m_number, 0);
  EXPECT_EQ(mi[2].m_name, "c");
}

TEST(mir_equiz2d_seq, equiz_uniform_unibuffer_seq)
{
  AXOM_ANNOTATE_SCOPE("equiz_uniform_unibuffer_seq");
  test_equiz_uniform_unibuffer<seq_exec>();
}

TEST(mir_equiz2d_seq, equiz_polygonal_unibuffer_seq)
{
  AXOM_ANNOTATE_SCOPE("equiz_polygonal_unibuffer_seq");
  test_Polygonal_MIR<seq_exec>::test("equiz_polygonal_unibuffer");
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return TestApp.execute(argc, argv);
}
