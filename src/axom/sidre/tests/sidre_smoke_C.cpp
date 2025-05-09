// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/sidre/interface/sidre.h"

//------------------------------------------------------------------------------

TEST(C_sidre_smoke, create_datastore)
{
  SIDRE_DataStore ds_buf;

  SIDRE_DataStore* ds = SIDRE_DataStore_new(&ds_buf);
  SIDRE_DataStore_delete(ds);
  EXPECT_TRUE(true);
}

//------------------------------------------------------------------------------

TEST(sidre_smoke, valid_invalid)
{
  SIDRE_DataStore ds_buf;
  SIDRE_Group root_buf;

  SIDRE_DataStore* ds = SIDRE_DataStore_new(&ds_buf);

  SIDRE_IndexType idx = 3;
  EXPECT_TRUE(idx != SIDRE_InvalidIndex);

  const char* name = "foo";
  EXPECT_TRUE(SIDRE_name_is_valid(name));

  SIDRE_Group* root = SIDRE_DataStore_get_root(ds, &root_buf);

  const char* gp_name = SIDRE_Group_get_group_name(root, idx);
  EXPECT_TRUE(gp_name == NULL);
  EXPECT_TRUE(gp_name == SIDRE_InvalidName);
  EXPECT_FALSE(SIDRE_name_is_valid(gp_name));
  EXPECT_TRUE(SIDRE_Group_get_group_index(root, name) == SIDRE_InvalidIndex);

  SIDRE_DataStore_delete(ds);
}
