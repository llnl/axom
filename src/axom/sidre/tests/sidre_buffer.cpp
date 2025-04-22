// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/sidre.hpp"

using axom::sidre::Buffer;
using axom::sidre::DataStore;
using axom::sidre::DataType;
using axom::sidre::Group;
using axom::sidre::IndexType;
using axom::sidre::INT_ID;
using axom::sidre::View;

constexpr IndexType ZERO = 0;

//------------------------------------------------------------------------------

// id's are created in sequence.
// After detaching or destroying a buffer, the next createBuffer will recycle
// the id.

TEST(sidre_buffer, create_buffers)
{
  DataStore* ds = new DataStore();
  EXPECT_EQ(0u, ds->getNumBuffers());

  // Create two buffers
  Buffer* dbuff_0 = ds->createBuffer();
  EXPECT_EQ(1u, ds->getNumBuffers());
  EXPECT_EQ(0, dbuff_0->getIndex());

  Buffer* dbuff_1 = ds->createBuffer();
  EXPECT_EQ(2u, ds->getNumBuffers());
  EXPECT_EQ(1, dbuff_1->getIndex());

  //----------
  // destroy by index
  ds->destroyBuffer(ZERO);
  EXPECT_EQ(1u, ds->getNumBuffers());

  Buffer* dbuff_0b = ds->createBuffer();
  EXPECT_EQ(2u, ds->getNumBuffers());
  EXPECT_EQ(0, dbuff_0b->getIndex());

  //----------
  // destroy by pointer
  ds->destroyBuffer(dbuff_1);
  EXPECT_EQ(1u, ds->getNumBuffers());

  Buffer* dbuff_1b = ds->createBuffer();
  EXPECT_EQ(2u, ds->getNumBuffers());
  EXPECT_EQ(1, dbuff_1b->getIndex());

  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_buffer, create_buffer_with_description)
{
  DataStore* ds = new DataStore();
  Buffer* dbuff = ds->createBuffer(INT_ID, 10);
  dbuff->allocate();

  EXPECT_EQ(dbuff->getTypeID(), INT_ID);
  EXPECT_EQ(dbuff->getNumElements(), 10u);
  EXPECT_EQ(dbuff->getTotalBytes(), static_cast<IndexType>(sizeof(int) * 10));

  int* data_ptr = dbuff->getData();

  for(int i = 0; i < 10; i++)
  {
    data_ptr[i] = i * i;
  }

  dbuff->print();

  ds->print();

  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_buffer, alloc_buffer_for_int_array)
{
  DataStore* ds = new DataStore();
  Buffer* dbuff = ds->createBuffer();

  dbuff->allocate(INT_ID, 10);
  // Should be a warning and no-op, buffer is already allocated, we don't want
  // to re-allocate and leak memory.
  dbuff->allocate();

  EXPECT_EQ(dbuff->getTypeID(), INT_ID);
  EXPECT_EQ(dbuff->getNumElements(), 10u);
  EXPECT_EQ(dbuff->getTotalBytes(), static_cast<IndexType>(sizeof(int) * 10));

  //  int * data_ptr = static_cast<int *>(dbuff->getData());
  int* data_ptr = dbuff->getData();

  for(int i = 0; i < 10; i++)
  {
    data_ptr[i] = i * i;
  }

  dbuff->print();

  ds->print();

  dbuff->deallocate();

  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_buffer, init_buffer_for_int_array)
{
  DataStore* ds = new DataStore();
  Buffer* dbuff = ds->createBuffer();

  dbuff->allocate(INT_ID, 10);

  EXPECT_EQ(dbuff->getTypeID(), INT_ID);
  EXPECT_EQ(dbuff->getNumElements(), 10u);
  EXPECT_EQ(dbuff->getTotalBytes(), static_cast<IndexType>(sizeof(int) * 10));

  int* data_ptr = static_cast<int*>(dbuff->getData());

  for(int i = 0; i < 10; i++)
  {
    data_ptr[i] = i * i;
  }

  dbuff->print();

  ds->print();
  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_buffer, realloc_buffer)
{
  DataStore* ds = new DataStore();
  Buffer* dbuff = ds->createBuffer();

  dbuff->allocate(INT_ID, 5);

  EXPECT_EQ(dbuff->getTypeID(), INT_ID);
  EXPECT_EQ(dbuff->getNumElements(), 5u);
  EXPECT_EQ(dbuff->getTotalBytes(), static_cast<IndexType>(sizeof(int) * 5));

  int* data_ptr = static_cast<int*>(dbuff->getVoidPtr());

  for(int i = 0; i < 5; i++)
  {
    data_ptr[i] = 5;
  }

  for(int i = 0; i < 5; i++)
  {
    EXPECT_EQ(data_ptr[i], 5);
  }

  //  dbuff->print();

  dbuff->reallocate(10);

  std::cerr << dbuff->getTypeID() << std::endl;
  EXPECT_EQ(dbuff->getTypeID(), INT_ID);
  EXPECT_EQ(dbuff->getNumElements(), 10u);
  EXPECT_EQ(dbuff->getTotalBytes(), static_cast<IndexType>(sizeof(int) * 10));

  // data buffer changes
  data_ptr = static_cast<int*>(dbuff->getVoidPtr());

  for(int i = 5; i < 10; i++)
  {
    data_ptr[i] = 10;
  }

  for(int i = 0; i < 10; i++)
  {
    int value = 5;
    if(i > 4)
    {
      value = 10;
    }
    EXPECT_EQ(data_ptr[i], value);
  }

  //  dbuff->print();

  //  ds->print();
  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_buffer, create_buffer_view)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  const IndexType len = 11;
  const int ndims = 1;
  IndexType shape[] = {len};

  Buffer* buff = ds->createBuffer(INT_ID, len)->allocate();

  int* idata = buff->getData();

  for(int ii = 0; ii < len; ++ii)
  {
    idata[ii] = ii;
  }

  for(unsigned int i = 0; i < 8; i++)
  {
    View* view;

    switch(i)
    {
    case 0:
      view = root->createView("data0", INT_ID, len, buff);
      break;
    case 1:
      view = root->createView("data1", INT_ID, len)->attachBuffer(buff);
      break;
    case 2:
      view = root->createView("data2")->attachBuffer(INT_ID, len, buff);
      break;
    case 3:
      view = root->createView("data3", buff)->apply(INT_ID, len);
      break;

    case 4:
      view = root->createViewWithShape("data4", INT_ID, ndims, shape, buff);
      break;
    case 5:
      view = root->createViewWithShape("data5", INT_ID, ndims, shape)->attachBuffer(buff);
      break;
    case 6:
      view = root->createView("data6")->attachBuffer(INT_ID, ndims, shape, buff);
      break;
    case 7:
      view = root->createView("data7", buff)->apply(INT_ID, ndims, shape);
      break;
    }

    EXPECT_EQ(root->getNumViews(), i + 1);

    EXPECT_TRUE(view->isDescribed());
    EXPECT_TRUE(view->isAllocated());
    EXPECT_TRUE(view->isApplied());

    EXPECT_TRUE(view->hasBuffer());
    EXPECT_EQ(buff, view->getBuffer());

    EXPECT_EQ(view->getTypeID(), INT_ID);
    EXPECT_EQ(view->getNumElements(), len);

    view->print();

    int* idata_chk = view->getData();
    for(int ii = 0; ii < len; ++ii)
    {
      EXPECT_EQ(idata_chk[ii], idata[ii]);
    }
  }
  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_buffer, with_multiple_views)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Buffer* dbuff;
  View *dv1, *dv2;

  dbuff = ds->createBuffer();
  IndexType idx = dbuff->getIndex();

  dv1 = root->createView("e1", dbuff);
  dv2 = root->createView("e2", dbuff);
  EXPECT_TRUE(dv1->hasBuffer());
  EXPECT_TRUE(dv2->hasBuffer());
  EXPECT_EQ(ds->getNumBuffers(), 1u);
  EXPECT_EQ(dbuff->getNumViews(), 2);

  // Detach buffer from first view will not detach from datastore.
  dv1->attachBuffer(nullptr);
  EXPECT_FALSE(dv1->hasBuffer());
  EXPECT_EQ(ds->getNumBuffers(), 1u);
  EXPECT_EQ(dbuff->getNumViews(), 1);

  // Detach buffer from second view will detach from datastore.
  dv2->attachBuffer(nullptr);
  EXPECT_FALSE(dv2->hasBuffer());
  EXPECT_EQ(ds->getNumBuffers(), 0u);

  // Buffer has been destroyed since there are no more attached views
  EXPECT_TRUE(ds->getBuffer(idx) == nullptr);

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_buffer, move_buffer)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Buffer *dbuff, *dbuff2;
  View *dv1, *dv2;

  dbuff = ds->createBuffer();

  dv1 = root->createView("e1", dbuff);

  dbuff2 = dv1->detachBuffer();
  EXPECT_FALSE(dv1->hasBuffer());
  EXPECT_EQ(dbuff2, dbuff);

  dv2 = root->createView("e2");
  dv2->attachBuffer(dbuff2);

  EXPECT_TRUE(dv2->hasBuffer());
  EXPECT_EQ(1u, ds->getNumBuffers());
  EXPECT_EQ(1, dbuff->getNumViews());

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_buffer, destroy_all_buffers)
{
  DataStore* ds = new DataStore();

  ds->createBuffer();
  ds->createBuffer();
  ds->createBuffer();
  ds->createBuffer();
  EXPECT_EQ(4u, ds->getNumBuffers());

  ds->destroyAllBuffers();
  EXPECT_EQ(0u, ds->getNumBuffers());

  delete ds;
}

#ifdef AXOM_USE_UMPIRE

class UmpireTest : public ::testing::TestWithParam<int>
{
public:
  void SetUp() override { allocID = GetParam(); }

  void TearDown() override
  {
    int allocID = axom::getUmpireResourceAllocatorID(umpire::resource::Host);
    axom::setDefaultAllocator(allocID);
  }

  DataStore ds;
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  int allocID;
};

//------------------------------------------------------------------------------
TEST_P(UmpireTest, allocate)
{
  constexpr int SIZE = 100;

  {
    Buffer* buff = ds.createBuffer();
    buff->describe(INT_ID, SIZE);
    buff->allocate(allocID);

    ASSERT_EQ(allocID, rm.getAllocator(buff->getVoidPtr()).getId());
    buff->deallocate();
  }

  {
    Buffer* buff = ds.createBuffer();
    buff->allocate(INT_ID, SIZE, allocID);

    ASSERT_EQ(allocID, rm.getAllocator(buff->getVoidPtr()).getId());
    buff->deallocate();
  }

  {
    Buffer* buff = ds.createBuffer(INT_ID, SIZE);
    buff->allocate(allocID);

    ASSERT_EQ(allocID, rm.getAllocator(buff->getVoidPtr()).getId());
    buff->deallocate();
  }
}

TEST_P(UmpireTest, allocate_default)
{
  constexpr int SIZE = 100;

  axom::setDefaultAllocator(allocID);

  {
    Buffer* buff = ds.createBuffer();
    buff->describe(INT_ID, SIZE);
    buff->allocate();

    ASSERT_EQ(allocID, rm.getAllocator(buff->getVoidPtr()).getId());
    buff->deallocate();
  }

  {
    Buffer* buff = ds.createBuffer();
    buff->allocate(INT_ID, SIZE);

    ASSERT_EQ(allocID, rm.getAllocator(buff->getVoidPtr()).getId());
    buff->deallocate();
  }

  {
    Buffer* buff = ds.createBuffer(INT_ID, SIZE);
    buff->allocate();

    ASSERT_EQ(allocID, rm.getAllocator(buff->getVoidPtr()).getId());
    buff->deallocate();
  }
}

TEST_P(UmpireTest, reallocate)
{
  constexpr int SIZE = 100;

  #if defined(AXOM_USE_GPU) && defined(UMPIRE_ENABLE_CONST)
  if(allocID == axom::getUmpireResourceAllocatorID(umpire::resource::Constant))
  {
    return;
  }
  #endif

  {
    Buffer* buff = ds.createBuffer();
    buff->describe(INT_ID, SIZE);
    buff->allocate(allocID);
    buff->reallocate(2 * SIZE);

    ASSERT_EQ(allocID, rm.getAllocator(buff->getVoidPtr()).getId());

    buff->deallocate();
  }

  {
    Buffer* buff = ds.createBuffer();
    buff->allocate(INT_ID, SIZE, allocID);
    buff->reallocate(2 * SIZE);

    ASSERT_EQ(allocID, rm.getAllocator(buff->getVoidPtr()).getId());
    buff->deallocate();
  }

  {
    Buffer* buff = ds.createBuffer(INT_ID, SIZE);
    buff->allocate(allocID);
    buff->reallocate(2 * SIZE);

    ASSERT_EQ(allocID, rm.getAllocator(buff->getVoidPtr()).getId());
    buff->deallocate();
  }
}

TEST_P(UmpireTest, reallocate_zero)
{
  constexpr int SIZE = 100;

  #if defined(AXOM_USE_GPU) && defined(UMPIRE_ENABLE_CONST)
  if(allocID == axom::getUmpireResourceAllocatorID(umpire::resource::Constant))
  {
    return;
  }
  #endif

  // set the default allocator
  axom::setDefaultAllocator(allocID);

  {
    Buffer* buff = ds.createBuffer();
    buff->describe(INT_ID, SIZE);
    buff->allocate(allocID);
    buff->reallocate(0);
    buff->reallocate(SIZE);

    ASSERT_EQ(axom::getDefaultAllocatorID(), rm.getAllocator(buff->getVoidPtr()).getId());

    buff->deallocate();
  }

  {
    Buffer* buff = ds.createBuffer();
    buff->allocate(INT_ID, SIZE, allocID);
    buff->reallocate(0);
    buff->reallocate(SIZE);

    ASSERT_EQ(axom::getDefaultAllocatorID(), rm.getAllocator(buff->getVoidPtr()).getId());
    buff->deallocate();
  }

  {
    Buffer* buff = ds.createBuffer(INT_ID, SIZE);
    buff->allocate(allocID);
    buff->reallocate(0);
    buff->reallocate(SIZE);

    ASSERT_EQ(axom::getDefaultAllocatorID(), rm.getAllocator(buff->getVoidPtr()).getId());
    buff->deallocate();
  }
}

const int allocators[] = {axom::getUmpireResourceAllocatorID(umpire::resource::Host)
  #ifdef AXOM_USE_GPU

    #ifdef UMPIRE_ENABLE_PINNED
                            ,
                          axom::getUmpireResourceAllocatorID(umpire::resource::Pinned)
    #endif

    #ifdef UMPIRE_ENABLE_DEVICE
                            ,
                          axom::getUmpireResourceAllocatorID(umpire::resource::Device)
    #endif

    #ifdef UMPIRE_ENABLE_CONST
                            ,
                          axom::getUmpireResourceAllocatorID(umpire::resource::Constant)
    #endif

    #ifdef UMPIRE_ENABLE_UM
                            ,
                          axom::getUmpireResourceAllocatorID(umpire::resource::Unified)
    #endif

  #endif /* defined(AXOM_USE_GPU) */
};

INSTANTIATE_TEST_SUITE_P(sidre_buffer, UmpireTest, ::testing::ValuesIn(allocators));

#endif  // AXOM_USE_UMPIRE
