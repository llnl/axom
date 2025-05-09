// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"
#include "axom/sidre.hpp"

using axom::sidre::DataStore;
using axom::sidre::Group;
using axom::sidre::indexIsValid;
using axom::sidre::IndexType;
using axom::sidre::InvalidIndex;
using axom::sidre::InvalidName;
using axom::sidre::nameIsValid;

//------------------------------------------------------------------------------

TEST(sidre_smoke, create_datastore)
{
  DataStore* ds = new DataStore();
  delete ds;
  EXPECT_TRUE(true);
}

//------------------------------------------------------------------------------

TEST(sidre_smoke, valid_invalid)
{
  DataStore* ds = new DataStore();

  IndexType idx = 3;
  EXPECT_TRUE(idx != InvalidIndex);

  std::string name("foo");
  EXPECT_TRUE(nameIsValid(name));

  Group* root = ds->getRoot();

  EXPECT_TRUE(root->getGroupName(idx) == InvalidName);
  EXPECT_TRUE(root->getGroupIndex(name) == InvalidIndex);

  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_smoke, conduit_in_sidre_smoke)
{
  // make sure we are linking with conduit ok.
  conduit::Node n;
  n["field"] = 100;
  EXPECT_EQ(n["field"].to_index_t(), 100u);
}
