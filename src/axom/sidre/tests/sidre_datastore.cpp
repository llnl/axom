// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"
#include "axom/sidre.hpp"

using axom::sidre::Buffer;
using axom::sidre::DataStore;
using axom::sidre::Group;
using axom::sidre::View;

using axom::sidre::DOUBLE_ID;
using axom::sidre::IndexType;
using axom::sidre::INT_ID;

#include <iostream>

//------------------------------------------------------------------------------

TEST(sidre_datastore, destroy_buffer)
{
  DataStore* ds = new DataStore();
  Buffer* dbuff = ds->createBuffer();

  IndexType bufferIndex = ds->getFirstValidBufferIndex();
  EXPECT_EQ(dbuff->getIndex(), bufferIndex);

  ds->destroyBuffer(bufferIndex);
  // should be no buffers
  EXPECT_TRUE(ds->getFirstValidBufferIndex() == -1);

  // After destroy, that buffer index should be available again, and have been
  // re-used..
  Buffer* dbuff2 = ds->createBuffer();
  (void)dbuff2;
  std::cerr << ds->getFirstValidBufferIndex() << std::endl;
  std::cerr << bufferIndex << std::endl;
  EXPECT_TRUE(ds->getFirstValidBufferIndex() == bufferIndex);

  // check error condition
  IndexType badBufferIndex = 9999;
  ds->destroyBuffer(badBufferIndex);

  delete ds;
}

TEST(sidre_datastore, buffer_info)
{
  IndexType num_buffers_chk = 0;
  IndexType num_buffers_referenced_chk = 0;
  IndexType num_buffers_detached_chk = 0;
  IndexType num_bytes_allocated_chk = 0;

  DataStore* ds = new DataStore();

  //
  // Ha-ha check
  //
  conduit::Node n0;
  ds->getBufferInfo(n0);

  IndexType num_buffers = n0["num_buffers"].value();
  IndexType num_buffers_referenced = n0["num_buffers_referenced"].value();
  IndexType num_buffers_detached = n0["num_buffers_detached"].value();
  IndexType num_bytes_allocated = n0["num_bytes_allocated"].value();

  num_buffers_detached_chk = ds->getNumBuffers() - num_buffers_referenced_chk;

  EXPECT_EQ(num_buffers_chk, ds->getNumBuffers());
  EXPECT_EQ(num_buffers_referenced_chk, num_buffers_referenced);
  EXPECT_EQ(num_buffers_detached_chk, num_buffers_detached);
  EXPECT_EQ(num_bytes_allocated_chk, num_bytes_allocated);

  //
  // More complex checks...
  //
  ds->createBuffer(INT_ID, 5);
  num_bytes_allocated_chk += 0;

  Buffer* buff1 = ds->createBuffer(INT_ID, 10)->allocate();
  num_bytes_allocated_chk += buff1->getTotalBytes();

  Buffer* buff2 = ds->createBuffer(DOUBLE_ID, 10)->allocate();
  num_bytes_allocated_chk += buff2->getTotalBytes();

  Group* root = ds->getRoot();
  Group* group_a = root->createGroup("a");
  Group* group_b = root->createGroup("b");

  View* view_a = group_a->createViewAndAllocate("fielda", INT_ID, 5);
  num_buffers_referenced_chk++;
  num_bytes_allocated_chk += view_a->getTotalBytes();

  View* view_b = group_b->createView("fieldb");
  view_b->attachBuffer(DOUBLE_ID, 5, buff2);
  num_buffers_referenced_chk++;

  conduit::Node n;
  ds->getBufferInfo(n);

  num_buffers = n["num_buffers"].value();
  num_buffers_referenced = n["num_buffers_referenced"].value();
  num_buffers_detached = n["num_buffers_detached"].value();
  num_bytes_allocated = n["num_bytes_allocated"].value();

  num_buffers_detached_chk = ds->getNumBuffers() - num_buffers_referenced_chk;

  EXPECT_EQ(num_buffers, ds->getNumBuffers());
  EXPECT_EQ(num_buffers_referenced_chk, num_buffers_referenced);
  EXPECT_EQ(num_buffers_detached_chk, num_buffers_detached);
  EXPECT_EQ(num_bytes_allocated_chk, num_bytes_allocated);

  delete ds;
}
