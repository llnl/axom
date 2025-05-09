// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"
#include "axom/sidre.hpp"

using axom::sidre::Buffer;
using axom::sidre::DataStore;
using axom::sidre::Group;
using axom::sidre::indexIsValid;
using axom::sidre::IndexType;
using axom::sidre::InvalidIndex;
using axom::sidre::View;

#include <map>

//------------------------------------------------------------------------------

void verifyEmptyGroupNamed(Group* dg, std::string name)
{
  EXPECT_EQ(name, dg->getName());

  EXPECT_EQ(0UL, dg->getNumGroups());
  EXPECT_FALSE(dg->hasGroup(-1));
  EXPECT_FALSE(dg->hasGroup(0));
  EXPECT_FALSE(dg->hasGroup(1));
  EXPECT_FALSE(dg->hasGroup("some_name"));
  EXPECT_EQ(InvalidIndex, dg->getGroupIndex("some_other_name"));
  EXPECT_EQ(InvalidIndex, dg->getFirstValidGroupIndex());
  EXPECT_EQ(InvalidIndex, dg->getNextValidGroupIndex(0));
  EXPECT_EQ(InvalidIndex, dg->getNextValidGroupIndex(4));

  EXPECT_EQ(0UL, dg->getNumViews());
  EXPECT_FALSE(dg->hasView(-1));
  EXPECT_FALSE(dg->hasView(0));
  EXPECT_FALSE(dg->hasView(1));
  EXPECT_FALSE(dg->hasView("some_name"));
  EXPECT_EQ(InvalidIndex, dg->getViewIndex("some_other_name"));
  EXPECT_EQ(InvalidIndex, dg->getFirstValidViewIndex());
  EXPECT_EQ(InvalidIndex, dg->getNextValidViewIndex(0));
  EXPECT_EQ(InvalidIndex, dg->getNextValidViewIndex(4));
}

void verifyBufferIdentity(DataStore* ds, std::map<IndexType, Buffer*>& bs)
{
  int bufcount = bs.size();

  // Does ds contain the number of buffers we expect?
  EXPECT_EQ(bufcount, static_cast<int>(ds->getNumBuffers()));

  // Does ds contain the buffer IDs and pointers we expect?
  int iteratedCount = 0;
  IndexType idx = ds->getFirstValidBufferIndex();
  while(indexIsValid(idx) && iteratedCount < bufcount)
  {
    EXPECT_TRUE(bs.count(idx) == 1);
    if(bs.count(idx) == 1)
    {
      EXPECT_EQ(bs[idx], ds->getBuffer(idx));
    }
    idx = ds->getNextValidBufferIndex(idx);
    iteratedCount += 1;
  }

  // Have we iterated over exactly the number of buffers we expect, finishing on
  // InvalidIndex?
  EXPECT_EQ(bufcount, iteratedCount);
  EXPECT_EQ(InvalidIndex, idx);
}

TEST(sidre_datastore, default_ctor)
{
  DataStore* ds = new DataStore();

  // After construction, the DataStore should contain no buffers.

  EXPECT_EQ(0, static_cast<int>(ds->getNumBuffers()));
  EXPECT_FALSE(ds->hasBuffer(-15));
  EXPECT_FALSE(ds->hasBuffer(-1));
  EXPECT_FALSE(ds->hasBuffer(0));
  EXPECT_FALSE(ds->hasBuffer(1));
  EXPECT_FALSE(ds->hasBuffer(8));

  EXPECT_EQ(InvalidIndex, ds->getFirstValidBufferIndex());
  EXPECT_EQ(InvalidIndex, ds->getNextValidBufferIndex(0));
  EXPECT_EQ(InvalidIndex, ds->getNextValidBufferIndex(4));

  // The new DataStore should contain exactly one group, the root group.
  // The root group should be named "/" and should contain no views and no
  // groups.

  Group* dg = ds->getRoot();

  EXPECT_FALSE(nullptr == dg);
  EXPECT_EQ(dg, dg->getParent());
  EXPECT_EQ(ds, dg->getDataStore());

  verifyEmptyGroupNamed(dg, "");

  delete ds;
}

// The dtor destroys all buffers and deletes the root group.
// An outside tool (valgrind) is used to check for proper memory cleanup.

TEST(sidre_datastore, create_destroy_buffers_basic)
{
  DataStore* ds = new DataStore();
  EXPECT_EQ(0, static_cast<int>(ds->getNumBuffers()));

  // Basic tests

  Buffer* dbuff = ds->createBuffer();
  EXPECT_EQ(1, static_cast<int>(ds->getNumBuffers()));

  IndexType bufferIndex = ds->getFirstValidBufferIndex();
  EXPECT_EQ(0, dbuff->getIndex());
  EXPECT_EQ(dbuff->getIndex(), bufferIndex);
  EXPECT_EQ(InvalidIndex, ds->getNextValidBufferIndex(bufferIndex));

  // Do we get the buffer we expect?
  EXPECT_EQ(dbuff, ds->getBuffer(bufferIndex));
  IndexType badBufferIndex = 9999;
  EXPECT_EQ(static_cast<void*>(nullptr), ds->getBuffer(badBufferIndex));

  ds->destroyBuffer(bufferIndex);
  // should be no buffers
  EXPECT_EQ(0, static_cast<int>(ds->getNumBuffers()));
  EXPECT_EQ(InvalidIndex, ds->getFirstValidBufferIndex());
  EXPECT_FALSE(ds->hasBuffer(bufferIndex));
  EXPECT_EQ(static_cast<void*>(nullptr), ds->getBuffer(bufferIndex));
  EXPECT_EQ(static_cast<void*>(nullptr), ds->getBuffer(badBufferIndex));

  delete ds;
}

TEST(sidre_datastore, create_destroy_buffers_order)
{
  DataStore* ds = new DataStore();
  EXPECT_EQ(0, static_cast<int>(ds->getNumBuffers()));

  Buffer* dbuff = ds->createBuffer();
  EXPECT_EQ(1, static_cast<int>(ds->getNumBuffers()));

  IndexType bufferIndex = ds->getFirstValidBufferIndex();
  ds->destroyBuffer(dbuff);

  // After destroy, test that buffer index should be available again for reuse.
  Buffer* dbuff2 = ds->createBuffer(axom::sidre::FLOAT32_ID, 16);
  IndexType d2Index = dbuff2->getIndex();
  EXPECT_EQ(bufferIndex, ds->getFirstValidBufferIndex());
  EXPECT_EQ(bufferIndex, d2Index);
  EXPECT_TRUE(ds->hasBuffer(d2Index));
  EXPECT_EQ(dbuff2, ds->getBuffer(bufferIndex));

  Buffer* dbuff3 = ds->createBuffer();
  IndexType d3Index = dbuff3->getIndex();
  EXPECT_EQ(2, static_cast<int>(ds->getNumBuffers()));
  EXPECT_TRUE(ds->hasBuffer(d3Index));

  // Try destroying the first valid buffer; see if we have the correct count and
  // indices
  ds->destroyBuffer(bufferIndex);
  EXPECT_EQ(1, static_cast<int>(ds->getNumBuffers()));
  EXPECT_FALSE(ds->hasBuffer(d2Index));
  EXPECT_TRUE(ds->hasBuffer(d3Index));

  // Add some more buffers, then try destroying the second one; see if we have
  // the correct count and indices
  Buffer* dbuff4 = ds->createBuffer();
  Buffer* dbuff5 = ds->createBuffer();
  IndexType d4Index = dbuff4->getIndex();
  IndexType d5Index = dbuff5->getIndex();
  EXPECT_EQ(3, static_cast<int>(ds->getNumBuffers()));
  EXPECT_EQ(d4Index, bufferIndex);  // dbuff4 should have recycled index 0
  EXPECT_TRUE(ds->hasBuffer(d3Index));
  EXPECT_TRUE(ds->hasBuffer(d4Index));
  EXPECT_TRUE(ds->hasBuffer(d5Index));

  ds->destroyBuffer(d3Index);  // Destroy dbuff3 (not dbuff4) because we already
                               // tested destroying index 0
  EXPECT_EQ(2, static_cast<int>(ds->getNumBuffers()));
  EXPECT_FALSE(ds->hasBuffer(d3Index));
  EXPECT_TRUE(ds->hasBuffer(d4Index));
  EXPECT_TRUE(ds->hasBuffer(d5Index));

  // Can we destroy all buffers?
  ds->destroyAllBuffers();
  EXPECT_EQ(0, static_cast<int>(ds->getNumBuffers()));
  EXPECT_FALSE(ds->hasBuffer(d2Index));
  EXPECT_FALSE(ds->hasBuffer(d4Index));
  EXPECT_FALSE(ds->hasBuffer(d5Index));

  delete ds;
}

TEST(sidre_datastore, create_destroy_buffers_views)
{
  DataStore* ds = new DataStore();
  EXPECT_EQ(0, static_cast<int>(ds->getNumBuffers()));

  Buffer* dbuff3 = ds->createBuffer();
  IndexType d3Index = dbuff3->getIndex();
  Buffer* dbuff4 = ds->createBuffer();
  IndexType d4Index = dbuff4->getIndex();
  Buffer* dbuff5 = ds->createBuffer();
  IndexType d5Index = dbuff5->getIndex();
  Buffer* dbuff6 = ds->createBuffer();
  IndexType d6Index = dbuff6->getIndex();

  EXPECT_EQ(dbuff3, ds->getBuffer(d3Index));
  EXPECT_EQ(dbuff4, ds->getBuffer(d4Index));
  EXPECT_EQ(dbuff5, ds->getBuffer(d5Index));
  EXPECT_EQ(dbuff6, ds->getBuffer(d6Index));

  // Create and verify views referencing buffers
  View* vA = ds->getRoot()->createView("vA", dbuff3);
  View* vB = ds->getRoot()->createView("vB", dbuff3);
  View* vC = ds->getRoot()->createView("vC", dbuff4);
  View* vD = ds->getRoot()->createView("vD", dbuff6);
  View* vE = ds->getRoot()->createView("vE", dbuff6);
  EXPECT_EQ(dbuff3, vA->getBuffer());
  EXPECT_EQ(dbuff3, vB->getBuffer());
  EXPECT_EQ(2, dbuff3->getNumViews());
  EXPECT_EQ(dbuff4, vC->getBuffer());
  EXPECT_EQ(dbuff6, vD->getBuffer());
  EXPECT_EQ(dbuff6, vE->getBuffer());
  EXPECT_EQ(2, dbuff6->getNumViews());

  // Destroying a buffer should detach it from the view
  ds->destroyBuffer(d4Index);
  EXPECT_EQ(3, static_cast<int>(ds->getNumBuffers()));
  EXPECT_TRUE(ds->hasBuffer(d3Index));
  EXPECT_FALSE(ds->hasBuffer(d4Index));
  EXPECT_TRUE(ds->hasBuffer(d5Index));
  EXPECT_TRUE(ds->hasBuffer(d6Index));
  EXPECT_EQ(dbuff3, vA->getBuffer());
  EXPECT_EQ(dbuff3, vB->getBuffer());
  EXPECT_FALSE(vC->hasBuffer());
  EXPECT_EQ(dbuff6, vD->getBuffer());
  EXPECT_EQ(dbuff6, vE->getBuffer());

  // Destroying a buffer should detach it from all of its views
  ds->destroyBuffer(d3Index);
  EXPECT_EQ(2, static_cast<int>(ds->getNumBuffers()));
  EXPECT_FALSE(ds->hasBuffer(d3Index));
  EXPECT_FALSE(ds->hasBuffer(d4Index));
  EXPECT_TRUE(ds->hasBuffer(d5Index));
  EXPECT_TRUE(ds->hasBuffer(d6Index));
  EXPECT_FALSE(vA->hasBuffer());
  EXPECT_FALSE(vB->hasBuffer());
  EXPECT_FALSE(vC->hasBuffer());
  EXPECT_EQ(dbuff6, vD->getBuffer());
  EXPECT_EQ(dbuff6, vE->getBuffer());

  // Destroying all buffers should detach them from all of their views
  dbuff3 = ds->createBuffer();
  dbuff4 = ds->createBuffer();
  vA->attachBuffer(dbuff3);
  vB->attachBuffer(dbuff3);
  vC->attachBuffer(dbuff4);
  ds->destroyAllBuffers();
  EXPECT_EQ(0, static_cast<int>(ds->getNumBuffers()));
  EXPECT_FALSE(vA->hasBuffer());
  EXPECT_FALSE(vB->hasBuffer());
  EXPECT_FALSE(vC->hasBuffer());
  EXPECT_FALSE(vD->hasBuffer());
  EXPECT_FALSE(vE->hasBuffer());

  // More tests will be found in sidre_group_unit.cpp.

  delete ds;
}

int psrand(int min, int max)
{
  // Returns a pseudorandom int in [min, max].  Note the closed interval.
  std::int64_t range = max - min + 1;
  return min + int(range * rand() / (RAND_MAX + 1.0));
}

int irhall(int n)
{
  // Approximates the Irwin-Hall distribution, a sum of uniformly-distributed
  // samples that approaches Gaussian distribution
  int retval = 0;

  for(int i = 0; i < n; ++i)
  {
    retval += psrand(-5, 5);
  }

  return retval;
}

// Test iteration through buffers, as well as proper index and buffer behavior
// while buffers are created and deleted
TEST(sidre_datastore, iterate_buffers_basic)
{
  DataStore* ds = new DataStore();
  EXPECT_EQ(0, static_cast<int>(ds->getNumBuffers()));

  IndexType badBufferIndex = 9999;
  // Do we get sidre::InvalidIndex for several queries with no buffers?
  EXPECT_EQ(InvalidIndex, ds->getFirstValidBufferIndex());
  EXPECT_EQ(InvalidIndex, ds->getNextValidBufferIndex(0));
  EXPECT_EQ(InvalidIndex, ds->getNextValidBufferIndex(badBufferIndex));
  EXPECT_EQ(InvalidIndex, ds->getNextValidBufferIndex(InvalidIndex));

  // Create one data buffer, verify its index is zero, and that iterators behave
  // as expected
  Buffer* initial = ds->createBuffer();
  EXPECT_EQ(0, initial->getIndex());
  EXPECT_EQ(1, static_cast<int>(ds->getNumBuffers()));
  EXPECT_EQ(0, ds->getFirstValidBufferIndex());
  EXPECT_EQ(InvalidIndex, ds->getNextValidBufferIndex(0));
  // Destroy the data buffer, verify that iterators behave as expected
  ds->destroyBuffer(initial);
  EXPECT_EQ(0, static_cast<int>(ds->getNumBuffers()));
  EXPECT_EQ(InvalidIndex, ds->getFirstValidBufferIndex());
  EXPECT_EQ(InvalidIndex, ds->getNextValidBufferIndex(0));

  delete ds;
}

// Test a few buffers: can we create them and iterate?
TEST(sidre_datastore, iterate_buffers_simple)
{
  DataStore* ds = new DataStore();
  EXPECT_EQ(0, static_cast<int>(ds->getNumBuffers()));

  std::map<IndexType, Buffer*> bs;
  int bufcount = 20;

  for(int i = 0; i < bufcount; ++i)
  {
    Buffer* b = ds->createBuffer(axom::sidre::FLOAT64_ID, 400 * i);
    IndexType idx = b->getIndex();
    bs[idx] = b;
  }

  verifyBufferIdentity(ds, bs);
  delete ds;
}

TEST(sidre_datastore, iterate_buffers_iterators)
{
  DataStore* ds = new DataStore();
  EXPECT_EQ(0, static_cast<int>(ds->getNumBuffers()));

  std::map<IndexType, Buffer*> bs;
  int bufcount = 20;

  for(int i = 0; i < bufcount; ++i)
  {
    Buffer* b = ds->createBuffer(axom::sidre::FLOAT64_ID, 400 * i);
    IndexType idx = b->getIndex();
    bs[idx] = b;
  }

  // non-const iterator on non--const datastore
  {
    int foundBuffers = 0;
    for(auto& buff : ds->buffers())
    {
      const IndexType idx = buff.getIndex();
      ++foundBuffers;
      EXPECT_TRUE(indexIsValid(idx));
      EXPECT_EQ(ds->getBuffer(idx), &buff);
      EXPECT_EQ(bs[idx], &buff);
    }
    EXPECT_EQ(bufcount, foundBuffers);
  }

  // const iterator on non-const datastore
  {
    int foundBuffers = 0;
    for(const auto& buff : ds->buffers())
    {
      const IndexType idx = buff.getIndex();
      ++foundBuffers;
      EXPECT_TRUE(indexIsValid(idx));
      EXPECT_EQ(ds->getBuffer(idx), &buff);
      EXPECT_EQ(bs[idx], &buff);
    }
    EXPECT_EQ(bufcount, foundBuffers);
  }

  // const iterator on const datastore
  {
    const DataStore* cds = ds;
    int foundBuffers = 0;
    for(auto& buff : cds->buffers())
    {
      const IndexType idx = buff.getIndex();
      ++foundBuffers;
      EXPECT_TRUE(indexIsValid(idx));
      EXPECT_EQ(cds->getBuffer(idx), &buff);
      EXPECT_EQ(bs[idx], &buff);
    }
    EXPECT_EQ(bufcount, foundBuffers);
  }

  delete ds;
}

// Test creating and allocating buffers, then destroying several of them
TEST(sidre_datastore, create_delete_buffers_iterate)
{
  DataStore* ds = new DataStore();
  EXPECT_EQ(0, static_cast<int>(ds->getNumBuffers()));

  std::map<IndexType, Buffer*> bs;
  int bufcount = 50;  // Arbitrary number of buffers

  // Initially, create some buffers of varying size
  for(int i = 0; i < bufcount; ++i)
  {
    Buffer* b = ds->createBuffer(axom::sidre::FLOAT64_ID, 400 * i % 10000)->allocate();
    IndexType idx = b->getIndex();
    bs[idx] = b;
  }

  int i = 0;
  std::map<IndexType, Buffer*> nbs;
  std::map<IndexType, Buffer*>::iterator bsit = bs.begin(), bsend = bs.end();
  for(; bsit != bsend; ++bsit)
  {
    // Eliminate some buffers (arbitrarily chosen)
    if(i % 5 && i % 7)
    {
      nbs[bsit->first] = bsit->second;
    }
    else
    {
      ds->destroyBuffer(bsit->first);
    }
    i += 1;
  }

  verifyBufferIdentity(ds, nbs);
  delete ds;
}

TEST(sidre_datastore, iterate_buffers_with_delete_iterators)
{
  DataStore* ds = new DataStore();
  EXPECT_EQ(0, static_cast<int>(ds->getNumBuffers()));

  const int initBuffCount = 22;
  for(int i = 0; i < initBuffCount; ++i)
  {
    ds->createBuffer(axom::sidre::FLOAT64_ID, 400 * i);
  }
  EXPECT_EQ(initBuffCount, ds->getNumBuffers());

  // remove a few buffers
  for(auto idx : {5, 10, 15, 20, 25})
  {
    if(ds->hasBuffer(idx))
    {
      ds->destroyBuffer(idx);
    }
  }
  int expBuffCount = 18;
  EXPECT_EQ(expBuffCount, ds->getNumBuffers());

  // Add a buffer
  {
    EXPECT_FALSE(ds->hasBuffer(5));
    auto* buff = ds->createBuffer(axom::sidre::FLOAT64_ID, 10);
    IndexType idx = buff->getIndex();
    EXPECT_TRUE(idx < initBuffCount);
    EXPECT_TRUE(ds->hasBuffer(idx));

    expBuffCount = 19;
    EXPECT_EQ(expBuffCount, ds->getNumBuffers());
  }

  // remove a few more
  for(auto idx : {3, 6, 9, 12, 15, 18, 21})
  {
    if(ds->hasBuffer(idx))
    {
      auto* buff = ds->getBuffer(idx);
      ds->destroyBuffer(buff);
    }
  }
  expBuffCount = 13;
  EXPECT_EQ(expBuffCount, ds->getNumBuffers());

  // iterate using 'range-for' syntax
  {
    int foundBuffers = 0;
    for(auto& buff : ds->buffers())
    {
      const IndexType idx = buff.getIndex();
      ++foundBuffers;
      EXPECT_TRUE(indexIsValid(idx));
      EXPECT_EQ(ds->getBuffer(idx), &buff);
    }
    EXPECT_EQ(expBuffCount, foundBuffers);
  }

  // iterate using begin/end iterator syntax
  {
    int foundBuffers = 0;
    for(auto it = ds->buffers().begin(), itEnd = ds->buffers().end(); it != itEnd; ++it)
    {
      auto& buff = *it;
      const IndexType idx = buff.getIndex();
      ++foundBuffers;
      EXPECT_TRUE(indexIsValid(idx));
      EXPECT_EQ(ds->getBuffer(idx), &buff);
    }
    EXPECT_EQ(expBuffCount, foundBuffers);
  }

  // iterate
  {
    int foundBuffers = 0;
    for(const auto& buff : ds->buffers())
    {
      const IndexType idx = buff.getIndex();
      ++foundBuffers;
      EXPECT_TRUE(indexIsValid(idx));
      EXPECT_EQ(ds->getBuffer(idx), &buff);
    }
    EXPECT_EQ(expBuffCount, foundBuffers);
  }

  delete ds;
}

// Test creating+allocating buffers, then destroying several of them, repeatedly
TEST(sidre_datastore, loop_create_delete_buffers_iterate)
{
  DataStore* ds = new DataStore();
  EXPECT_EQ(0, static_cast<int>(ds->getNumBuffers()));

  std::map<IndexType, Buffer*> bs;
  std::vector<int> idxlist;
  int initbufcount = 50;  // Arbitrary number of buffers

  // Initially, create some buffers of varying size
  for(int i = 0; i < initbufcount; ++i)
  {
    Buffer* b = ds->createBuffer(axom::sidre::FLOAT64_ID, 400 * i % 10000)->allocate();
    IndexType idx = b->getIndex();
    bs[idx] = b;
    idxlist.push_back(idx);
  }

  int totalrounds = 100;  // Arbitrary number of rounds
  for(int round = 0; round < totalrounds; ++round)
  {
    SCOPED_TRACE(round);

    // In each round, choose a random number of buffers to delete or create
    {
      int delta = irhall(5);  // Arbitrary argument to Irwin-Hall distribution
      SCOPED_TRACE(delta);

      if(delta < 0)
      {
        int rmvcount = abs(delta);
        if(rmvcount > static_cast<int>(bs.size()))
        {
          rmvcount = bs.size();
        }

        for(int i = 0; i < rmvcount; ++i)
        {
          int rmvidx = psrand(0, idxlist.size() - 1);
          int rmvid = idxlist[rmvidx];
          EXPECT_TRUE(ds->hasBuffer(rmvid));
          EXPECT_TRUE(bs.count(rmvid) == 1);
          ds->destroyBuffer(rmvid);
          bs.erase(rmvid);
          idxlist.erase(idxlist.begin() + rmvidx);
          EXPECT_FALSE(ds->hasBuffer(rmvid));
          EXPECT_FALSE(bs.count(rmvid) == 1);
        }
      }
      else if(delta > 0)
      {
        int addcount = delta;
        for(int i = 0; i < addcount; ++i)
        {
          Buffer* buf = ds->createBuffer(axom::sidre::FLOAT64_ID, 400)->allocate();
          int addid = buf->getIndex();
          EXPECT_TRUE(ds->hasBuffer(addid));
          EXPECT_TRUE(bs.count(addid) < 1);
          bs[addid] = buf;
          idxlist.push_back(addid);
        }
      }
    }

    verifyBufferIdentity(ds, bs);
  }

  delete ds;
}
