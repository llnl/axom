// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"
#include "axom/core/Types.hpp"
#include "axom/sidre/core/ConduitMemory.hpp"
#include "axom/slic.hpp"
#include "axom/sidre.hpp"

using axom::sidre::Buffer;
using axom::sidre::CHAR8_STR_ID;
using axom::sidre::DataStore;
using axom::sidre::DOUBLE_ID;
using axom::sidre::Group;
using axom::sidre::IndexType;
using axom::sidre::INT_ID;
using axom::sidre::NO_TYPE_ID;
using axom::sidre::TypeID;
using axom::sidre::View;

using namespace conduit;

#define BLEN 10

// Match ViewBuffer states
enum State
{
  EMPTY,
  BUFFER,
  EXTERNAL,
  TUPLE,
  STRING,
  NOTYPE
};

// Check state of View based on booleans since m_state is private.

static State getState(View* view)
{
  if(view->isEmpty())
  {
    return EMPTY;
  }
  else if(view->hasBuffer())
  {
    return BUFFER;
  }
  else if(view->isExternal())
  {
    return EXTERNAL;
  }
  else if(view->isScalar())
  {
    return TUPLE;
  }
  else if(view->isTuple())
  {
    return TUPLE;
  }
  else if(view->isString())
  {
    return STRING;
  }
  else
  {
    return NOTYPE;
  }
}

//
// Test various paths to creating a buffer.
// Assume all are int[len]
//
// Since this function is called in several places, any test failure
// will be ambiguous. To use, wrap in the EXPECT_TRUE macro.
//    EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));
// Any failures will have two reports, one from this routine and one
// from the call site.
//
static bool checkViewValues(View* view,
                            State state,
                            bool isDescribed,
                            bool isAllocated,
                            bool isApplied,
                            IndexType len)
{
  bool rv = true;
  IndexType dims[2];

  if(getState(view) != state)
  {
    EXPECT_EQ(getState(view), state);
    rv = false;
  }

  if(view->isDescribed() != isDescribed)
  {
    EXPECT_EQ(isDescribed, view->isDescribed());
    rv = false;
  }
  if(view->isAllocated() != isAllocated)
  {
    EXPECT_EQ(isAllocated, view->isAllocated());
    rv = false;
  }
  if(view->isApplied() != isApplied)
  {
    EXPECT_EQ(isApplied, view->isApplied());
    rv = false;
  }

  if(view->getNumElements() != len)
  {
    EXPECT_EQ(len, view->getNumElements());
    rv = false;
  }

  if(view->isDescribed())
  {
    if(view->getTypeID() != INT_ID)
    {
      EXPECT_EQ(INT_ID, view->getTypeID());
      rv = false;
    }
    if(view->getNumDimensions() != 1)
    {
      EXPECT_EQ(1, view->getNumDimensions());
      rv = false;
    }
    if(view->getShape(1, dims) != 1 || dims[0] != len)
    {
      EXPECT_TRUE(view->getShape(1, dims) == 1 && dims[0] == len);
      rv = false;
    }
    if(view->getTotalBytes() != static_cast<IndexType>(sizeof(int) * len))
    {
      EXPECT_EQ(view->getTotalBytes(), static_cast<IndexType>(sizeof(int) * len));
      rv = false;
    }
  }
  else
  {
    TypeID id = view->getTypeID();
    if(id != NO_TYPE_ID)
    {
      EXPECT_EQ(NO_TYPE_ID, id);
      rv = false;
    }
  }

  if(view->isAllocated() && view->isApplied())
  {
    // Fill with data to help the print function
    int* data_ptr = view->getData();
    int nelem = view->getNumElements();

    for(int i = 0; i < nelem; i++)
    {
      data_ptr[i] = i * i;
    }
  }

  //  view->print();

  return rv;
}

#if 0
//------------------------------------------------------------------------------

TEST(sidre_view,create_views)
{
  DataStore* ds   = new DataStore();
  Group* root = ds->getRoot();

  View* dv_0 = root->createViewAndAllocate("field0", INT_ID, 1);
  View* dv_1 = root->createViewAndAllocate("field1", INT_ID, 1);

  EXPECT_EQ(0, dv_0->getIndex());
  EXPECT_EQ(1, dv_1->getIndex());

  Buffer* db_0 = dv_0->getBuffer();
  Buffer* db_1 = dv_1->getBuffer();

  EXPECT_EQ(0, db_0->getIndex());
  EXPECT_EQ(1, db_1->getIndex());
  delete ds;
}
#endif

//------------------------------------------------------------------------------

TEST(sidre_view, get_path_name)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  View* v1 = root->createView("test/a/b/v1");
  View* v2 = root->createView("test/v2");
  View* v3 = root->createView("v3");
  View* v4 = root->createView("v4");

  EXPECT_EQ(std::string("v1"), v1->getName());
  EXPECT_EQ(std::string("test/a/b"), v1->getPath());
  EXPECT_EQ(std::string("test/a/b/v1"), v1->getPathName());

  EXPECT_EQ(std::string("v2"), v2->getName());
  EXPECT_EQ(std::string("test"), v2->getPath());
  EXPECT_EQ(std::string("test/v2"), v2->getPathName());

  EXPECT_EQ(0, v3->getIndex());
  EXPECT_EQ(std::string("v3"), v3->getName());
  EXPECT_EQ(std::string(""), v3->getPath());
  EXPECT_EQ(std::string("v3"), v3->getPathName());

  EXPECT_EQ(1, v4->getIndex());

  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_view, create_view_from_path)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  // Verify create works when groups must be created on demand.
  View* baz = root->createView("foo/bar/baz");
  // Groups should have been created.
  EXPECT_TRUE(root->hasGroup("foo"));
  EXPECT_TRUE(root->getGroup("foo")->hasGroup("bar"));

  Group* bar = root->getGroup("foo")->getGroup("bar");
  EXPECT_TRUE(bar->hasView("baz"));
  EXPECT_EQ(bar->getView("baz"), baz);

  (void)baz;

  delete ds;

#if 0
  ds = new DataStore();
  root = ds->getRoot();

  // Verify create works when groups already exist.
  baz = root->createView("foo/bar/baz");
  EXPECT_TRUE( baz != nullptr );

  EXPECT_TRUE( root->hasGroup("foo") );

  std::cerr << "HERE2" << std::endl;

  foo = root->getGroup("foo");
  EXPECT_TRUE( foo->hasGroup("bar") );

  std::cerr << "HERE3" << std::endl;

  bar = foo->getGroup("bar");
  EXPECT_TRUE( bar->hasView("baz"));
  EXPECT_EQ( baz, bar->getView("baz"));

  std::cerr << "HERE4" << std::endl;

  delete ds;
#endif
}

//------------------------------------------------------------------------------

static void checkScalarValues(View* view,
                              State state,
                              bool isDescribed,
                              bool isAllocated,
                              bool isApplied,
                              TypeID type,
                              IndexType len)
{
  IndexType dims[2];

  SCOPED_TRACE(view->getName());

  EXPECT_EQ(getState(view), state);

  EXPECT_EQ(view->isDescribed(), isDescribed);
  EXPECT_EQ(view->isAllocated(), isAllocated);
  EXPECT_EQ(view->isApplied(), isApplied);

  EXPECT_EQ(view->getTypeID(), type);
  EXPECT_EQ(view->getNumElements(), len);
  EXPECT_EQ(view->getNumDimensions(), 1);
  EXPECT_TRUE(view->getShape(1, dims) == 1 && dims[0] == len);
}

TEST(sidre_view, tuple_view)
{
  DataStore ds;
  Group* root = ds.getRoot();

  constexpr size_t N = 5;

  axom::Array<int> iArray(N, N);
  for(int i = 0; i < iArray.size(); ++i)
  {
    iArray[i] = i + 100;
  }

  View* iView = root->createView("i")->setTuple(iArray.view());
  EXPECT_EQ(iView->getNumElements(), iArray.size());
  EXPECT_TRUE(iView->isTuple());
  EXPECT_FALSE(iView->isScalar() || N == 1);
  EXPECT_FALSE(iView->isEmpty());
  EXPECT_FALSE(iView->isExternal());

  int* iPtr = static_cast<int*>(iView->getVoidPtr());
  EXPECT_NE(iPtr, iArray.data());

  for(int i = 0; i < iArray.size(); ++i)
  {
    EXPECT_EQ(iArray[i], iPtr[i]);
  }
}

TEST(sidre_view, scalar_view)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  int i;
  const char* s;

  View* i0view = root->createView("i0")->setScalar(1);
  checkScalarValues(i0view, TUPLE, true, true, true, INT_ID, 1);
  i = i0view->getScalar();
  EXPECT_EQ(1, i);

  View* i1view = root->createViewScalar("i1", 2);
  checkScalarValues(i1view, TUPLE, true, true, true, INT_ID, 1);
  i = i1view->getScalar();
  EXPECT_EQ(2, i);

  View* s0view = root->createView("s0")->setString("I am a string");
  checkScalarValues(s0view, STRING, true, true, true, CHAR8_STR_ID, 14);
  s = s0view->getString();
  EXPECT_TRUE(strcmp(s, "I am a string") == 0);

  View* s1view = root->createViewString("s1", "I too am a string");
  checkScalarValues(s1view, STRING, true, true, true, CHAR8_STR_ID, 18);
  s = s1view->getString();
  EXPECT_TRUE(strcmp(s, "I too am a string") == 0);

  // Check illegal operations
  i0view->apply(INT_ID, 1);
  i0view->allocate();
  i0view->deallocate();

  s0view->apply(INT_ID, 1);
  s0view->allocate();
  s0view->deallocate();

  View* empty = root->createView("empty");
  const char* svalue = empty->getString();
  EXPECT_EQ(nullptr, svalue);

  delete ds;
}

//------------------------------------------------------------------------------

// Most tests deallocate via the DataStore destructor
// This is an explicit deallocate test

TEST(sidre_view, dealloc)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Buffer* dbuff;
  View* dv;

  //----------  EMPTY(F,F,F)
  dv = root->createView("e1");
  EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));

  dv->deallocate();
  EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));

  //----------  EMPTY(T,F,F)
  dv = root->createView("e2", INT_ID, BLEN);
  EXPECT_TRUE(checkViewValues(dv, EMPTY, true, false, false, BLEN));

  dv->deallocate();
  EXPECT_TRUE(checkViewValues(dv, EMPTY, true, false, false, BLEN));

  //----------  BUFFER(F,F,F)
  dv = root->createView("b1");
  dbuff = ds->createBuffer();
  dv->attachBuffer(dbuff);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, false, false, false, 0));

  dv->deallocate();
  EXPECT_TRUE(checkViewValues(dv, BUFFER, false, false, false, 0));
  EXPECT_FALSE(dv->getBuffer()->isAllocated());

  //---------- BUFFER(T,F,F)
  dv = root->createView("b2", INT_ID, BLEN);
  dbuff = ds->createBuffer()->describe(INT_ID, BLEN);
  dv->attachBuffer(dbuff);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, false, false, BLEN));

  dv->deallocate();
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, false, false, BLEN));
  EXPECT_FALSE(dv->getBuffer()->isAllocated());

  //---------- BUFFER(T,T,T)
  dv = root->createView("b3", INT_ID, BLEN);
  dbuff = ds->createBuffer()->allocate(INT_ID, BLEN);
  dv->attachBuffer(dbuff);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));

  dv->deallocate();
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, false, false, BLEN));
  EXPECT_FALSE(dv->getBuffer()->isAllocated());

  delete ds;
}

//------------------------------------------------------------------------------

// allocate/reallocate with zero items results in allocated (yet zero).

TEST(sidre_view, alloc_zero_items)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  // Allocate View with zero items
  {
    bool expDesc = false;
    bool expAlloc = false;
    bool expAppl = false;

    View* dv = root->createView("z0");
    EXPECT_TRUE(checkViewValues(dv, EMPTY, expDesc, expAlloc, expAppl, 0));

    // Allocate with negative value leaves the view unallocated
    SLIC_INFO("Attempting View::allocate(-1).  Warning expected.");
    dv->allocate(INT_ID, -1);
    EXPECT_TRUE(checkViewValues(dv, EMPTY, expDesc, expAlloc, expAppl, 0));

    // Allocate with size zero
    dv->allocate(INT_ID, 0);
    expDesc = expAppl = true;
    expAlloc = true;
    EXPECT_TRUE(checkViewValues(dv, BUFFER, expDesc, expAlloc, expAppl, 0));
    EXPECT_TRUE(dv->getBuffer()->isAllocated());
    EXPECT_EQ(0, dv->getBuffer()->getNumElements());
  }

  // Reallocate View with zero items
  {
    bool expDesc = false;
    bool expAlloc = false;
    bool expAppl = false;

    View* dv = root->createView("z1");
    EXPECT_TRUE(checkViewValues(dv, EMPTY, expDesc, expAlloc, expAppl, 0));

    dv->allocate(INT_ID, BLEN);
    expDesc = expAlloc = expAppl = true;
    EXPECT_TRUE(checkViewValues(dv, BUFFER, expDesc, expAlloc, expAppl, BLEN));
    EXPECT_TRUE(dv->getBuffer()->isAllocated());
    EXPECT_EQ(BLEN, dv->getBuffer()->getNumElements());

    // Try to reallocate with negative size; no-op
    SLIC_INFO("Attempting View::reallocate(-1).  Warning expected.");
    dv->reallocate(-1);
    expDesc = expAlloc = expAppl = true;
    EXPECT_TRUE(checkViewValues(dv, BUFFER, expDesc, expAlloc, expAppl, BLEN));
    EXPECT_TRUE(dv->getBuffer()->isAllocated());
    EXPECT_EQ(BLEN, dv->getBuffer()->getNumElements());

    // Call realloc(0); view and associated buffer are resized
    dv->reallocate(0);
    expDesc = expAppl = true;
    expAlloc = true;
    EXPECT_TRUE(checkViewValues(dv, BUFFER, expDesc, expAlloc, expAppl, 0));
    EXPECT_TRUE(dv->getBuffer()->isAllocated());
    EXPECT_EQ(0, dv->getBuffer()->getNumElements());

    // Deallocate and then allocate() when described with zero items
    dv->deallocate();
    expDesc = true;
    expAlloc = expAppl = false;
    EXPECT_TRUE(checkViewValues(dv, BUFFER, expDesc, expAlloc, expAppl, 0));
    expAlloc = expAppl = true;
    dv->allocate();
    EXPECT_TRUE(checkViewValues(dv, BUFFER, expDesc, expAlloc, expAppl, 0));

    // Deallocate and then allocate(0) when described with zero items
    dv->deallocate();
    dv->allocate(INT_ID, 0);
    EXPECT_TRUE(checkViewValues(dv, BUFFER, expDesc, expAlloc, expAppl, 0));

    // Deallocate and then reallocate(), size zero
    dv->deallocate();
    dv->reallocate(0);
    EXPECT_TRUE(checkViewValues(dv, BUFFER, expDesc, expAlloc, expAppl, 0));

    // Deallocate and then reallocate(), non-zero sized
    dv->deallocate();
    dv->reallocate(BLEN);
    expDesc = expAlloc = expAppl = true;
    EXPECT_TRUE(checkViewValues(dv, BUFFER, expDesc, expAlloc, expAppl, BLEN));
  }

  // Allocate View of size zero into non-empty buffer
  {
    bool expDesc = true;
    bool expAlloc = true;
    bool expAppl = true;

    Buffer* nonEmptyBuf = ds->createBuffer(INT_ID, BLEN)->allocate();

    View* dv = root->createView("z_nonEmptyBuf_attach_apply");
    dv->attachBuffer(nonEmptyBuf)->apply(0);
    EXPECT_TRUE(dv->getBuffer()->isAllocated());
    EXPECT_TRUE(checkViewValues(dv, BUFFER, expDesc, expAlloc, expAppl, 0));
    EXPECT_EQ(BLEN, dv->getBuffer()->getNumElements());

    // reallocate view to have non-zero size and check that buffer is resized
    int sz = BLEN / 2;
    dv->reallocate(sz);
    EXPECT_TRUE(checkViewValues(dv, BUFFER, expDesc, expAlloc, expAppl, sz));
    EXPECT_EQ(sz, dv->getBuffer()->getNumElements());

    // reallocate view to have size zero and check that buffer is resized
    dv->reallocate(0);
    expDesc = expAppl = true;
    expAlloc = true;
    EXPECT_TRUE(checkViewValues(dv, BUFFER, expDesc, expAlloc, expAppl, 0));
    EXPECT_EQ(0, dv->getBuffer()->getNumElements());

    // reallocate view to have size BLEN and check that buffer is resized
    dv->reallocate(BLEN);
    expDesc = expAlloc = expAppl = true;
    EXPECT_TRUE(checkViewValues(dv, BUFFER, expDesc, expAlloc, expAppl, BLEN));
    EXPECT_EQ(BLEN, dv->getBuffer()->getNumElements());

    dv = root->createView("z_nonEmptyBuf_described_attach", INT_ID, 0);
    dv->attachBuffer(nonEmptyBuf);
    EXPECT_TRUE(checkViewValues(dv, BUFFER, expDesc, expAlloc, expAppl, 0));
    EXPECT_EQ(BLEN, dv->getBuffer()->getNumElements());
  }

  // Allocate View of size 0 into empty buffer
  {
    bool expDesc = true;
    bool expAlloc = true;
    bool expAppl = true;

    Buffer* emptyBuf = ds->createBuffer(INT_ID, 0)->allocate();
    View* dv = root->createView("z_emptyBuf_attach_apply");
    dv->attachBuffer(emptyBuf)->allocate()->apply(0);
    EXPECT_TRUE(dv->getBuffer()->isAllocated());
    EXPECT_TRUE(checkViewValues(dv, BUFFER, expDesc, expAlloc, expAppl, 0));

    // reallocate buffer with non-empty size; view should still be zero
    emptyBuf->reallocate(BLEN);
    expDesc = expAlloc = expAppl = true;
    EXPECT_TRUE(checkViewValues(dv, BUFFER, expDesc, expAlloc, expAppl, 0));
    EXPECT_TRUE(dv->getBuffer()->isAllocated());
    EXPECT_EQ(BLEN, dv->getBuffer()->getNumElements());

    // reallocate view to have size zero; view and buffer should be resized
    dv->reallocate(0);
    expDesc = expAlloc = expAppl = true;
    EXPECT_TRUE(checkViewValues(dv, BUFFER, expDesc, expAlloc, expAppl, 0));
    EXPECT_TRUE(dv->getBuffer()->isAllocated());
    EXPECT_EQ(0, dv->getBuffer()->getNumElements());

    // attach a second view with size zero
    auto* v2 = root->createView("z_emptyBuf_described_attach", INT_ID, 0);
    v2->attachBuffer(emptyBuf);
    expDesc = true;
    expAlloc = expAppl = true;
    EXPECT_TRUE(checkViewValues(dv, BUFFER, expDesc, expAlloc, expAppl, 0));
    EXPECT_EQ(0, dv->getBuffer()->getNumElements());
  }

  delete ds;
}

TEST(sidre_view, save_empty_view_non_empty_buffer)
{
  DataStore ds;
  Group* root = ds.getRoot();

  // allocate non-empty buffer
  const int SZ = 10;
  Buffer* buf = ds.createBuffer(INT_ID, SZ)->allocate();

  // attach zero-sized array to it
  root->createView("a")->attachBuffer(buf)->apply(0);

  // fill the buffer
  int* data = buf->getData();
  for(int i = 0; i < SZ; ++i)
  {
    data[i] = i;
  }

  root->save("empty_view_non_empty_buffer.sidre.json", "sidre_json");  // ok
#ifdef AXOM_USE_HDF5
  root->save("empty_view_non_empty_buffer.sidre.hdf5", "sidre_hdf5");  // ok
#endif
}

TEST(sidre_view, save_empty_view)
{
  DataStore ds;
  Group* root = ds.getRoot();

  root->createView("a", INT_ID, 0)->allocate()->apply();  // create View and
                                                          // allocate view of
                                                          // size 0

  root->save("empty_view.sidre.json", "sidre_json");  // this is ok
  // root->save("empty_view.sidre.hdf5", "sidre_hdf5");   // <-- problem here
  // (in conduit 0.2.1)
}

//------------------------------------------------------------------------------

// Test allocate, reallocate, and deallocate when there are multiple views
// attached to a buffer.  All operations should be no-ops.

TEST(sidre_view, alloc_and_dealloc_multiview)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Buffer* dbuff;
  View *dv1, *dv2;
  void* baddr;

  //---------- allocate
  dbuff = ds->createBuffer()->describe(INT_ID, BLEN);
  EXPECT_FALSE(dbuff->isAllocated());
  EXPECT_EQ(dbuff->getNumElements(), BLEN);
  baddr = dbuff->getVoidPtr();

  dv1 = root->createView("dv1alloc", dbuff);
  EXPECT_TRUE(checkViewValues(dv1, BUFFER, false, false, false, 0));
  dv2 = root->createView("dv2alloc", dbuff);
  EXPECT_TRUE(checkViewValues(dv2, BUFFER, false, false, false, 0));
  EXPECT_EQ(dbuff->getNumViews(), 2);

  dv1->allocate(INT_ID, BLEN + 10);
  EXPECT_TRUE(checkViewValues(dv2, BUFFER, false, false, false, 0));

  // buffer is unchanged
  EXPECT_FALSE(dbuff->isAllocated());
  EXPECT_EQ(dbuff->getNumElements(), BLEN);
  EXPECT_EQ(dbuff->getVoidPtr(), baddr);

  //---------- reallocate
  dbuff = ds->createBuffer()->allocate(INT_ID, BLEN);
  EXPECT_TRUE(dbuff->isAllocated());
  EXPECT_EQ(dbuff->getNumElements(), BLEN);
  baddr = dbuff->getVoidPtr();

  dv1 = root->createView("dv1realloc", dbuff);
  EXPECT_TRUE(checkViewValues(dv1, BUFFER, false, false, false, 0));
  dv2 = root->createView("dv2realloc", dbuff);
  EXPECT_TRUE(checkViewValues(dv1, BUFFER, false, false, false, 0));
  EXPECT_EQ(dbuff->getNumViews(), 2);

  dv1->reallocate(BLEN + 10);
  EXPECT_TRUE(checkViewValues(dv2, BUFFER, false, false, false, 0));

  // buffer is unchanged
  EXPECT_TRUE(dbuff->isAllocated());
  EXPECT_EQ(dbuff->getNumElements(), BLEN);
  EXPECT_EQ(dbuff->getVoidPtr(), baddr);

  //---------- deallocate
  dbuff = ds->createBuffer()->allocate(INT_ID, BLEN);
  EXPECT_TRUE(dbuff->isAllocated());
  EXPECT_EQ(dbuff->getNumElements(), BLEN);
  baddr = dbuff->getVoidPtr();

  dv1 = root->createView("dv1dealloc", dbuff);
  EXPECT_TRUE(checkViewValues(dv1, BUFFER, false, false, false, 0));
  dv2 = root->createView("dv2dealloc", dbuff);
  EXPECT_TRUE(checkViewValues(dv1, BUFFER, false, false, false, 0));
  EXPECT_EQ(dbuff->getNumViews(), 2);

  dv1->deallocate();
  EXPECT_TRUE(checkViewValues(dv2, BUFFER, false, false, false, 0));

  // buffer is unchanged
  EXPECT_TRUE(dbuff->isAllocated());
  EXPECT_EQ(dbuff->getNumElements(), BLEN);
  EXPECT_EQ(dbuff->getVoidPtr(), baddr);

  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_view, int_alloc_view)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  View* dv;
  IndexType shape[] = {BLEN};

  dv = root->createView("u0");
  EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));
  dv->allocate(INT_ID, BLEN);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));
#if 0
  dv = root->createView("u1");
  EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));
  dv->allocate(INT_ID, 1, shape);  // XXX - missing overload
  EXPECT_TRUE(checkViewValues(dv, BUFFEER, true, true, true, BLEN));
#endif
  dv = root->createView("u2");
  EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));
  dv->allocate(DataType::c_int(BLEN));
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));

  dv = root->createView("v0", INT_ID, 10);
  EXPECT_TRUE(checkViewValues(dv, EMPTY, true, false, false, BLEN));
  dv->allocate();
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));

  dv = root->createViewWithShape("v1", INT_ID, 1, shape);
  EXPECT_TRUE(checkViewValues(dv, EMPTY, true, false, false, BLEN));
  dv->allocate();
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));

  dv = root->createView("v2", DataType::c_int(BLEN));
  EXPECT_TRUE(checkViewValues(dv, EMPTY, true, false, false, BLEN));
  dv->allocate();
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));

  dv = root->createViewAndAllocate("a0", INT_ID, BLEN);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));

  dv = root->createViewWithShapeAndAllocate("a1", INT_ID, 1, shape);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));

  dv = root->createViewAndAllocate("a2", DataType::c_int(BLEN));
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));

  delete ds;
}

//------------------------------------------------------------------------------

// Test some association/state transitions.
// All of these tests have only one view per buffer.

TEST(sidre_view, int_buffer_view)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Buffer *dbuff, *otherbuffer;
  View* dv;
  IndexType bindex;

  //     view                        buffer
  //                  undescribed  described  allocated
  //     undescribed      1           2           3
  //     described        4           5           6
  //
  //  buffer description incomptable with view description
  //                               described  allocated
  //     described                    7           8

  //---------- 1
  // Attach undescribed buffer to undescribed view
  dv = root->createView("u1");
  EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));

  // no-op, attach nullptr buffer to EMPTY view
  dv->attachBuffer(nullptr);
  EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));

  dbuff = ds->createBuffer();
  dv->attachBuffer(dbuff);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, false, false, false, 0));
  EXPECT_EQ(dv->getBuffer(), dbuff);  // sanity check

  dv->allocate();  // no-op, no description
  EXPECT_TRUE(checkViewValues(dv, BUFFER, false, false, false, 0));

  dv->allocate(INT_ID, 10);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));

  dv->reallocate(BLEN + 5);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN + 5));

  // After deallocate, description is intact.
  dv->deallocate();
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, false, false, BLEN + 5));

  //---------- 2
  // Attach described buffer to undescribed view
  dv = root->createView("u2");
  EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));

  dbuff = ds->createBuffer()->describe(INT_ID, BLEN);
  dv->attachBuffer(dbuff);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, false, false, false, 0));

  //---------- 3
  // Attach allocated buffer to undescribed view
  dv = root->createView("u3");
  EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));

  dbuff = ds->createBuffer()->allocate(INT_ID, BLEN);
  dv->attachBuffer(dbuff);
  bindex = dbuff->getIndex();
  EXPECT_TRUE(checkViewValues(dv, BUFFER, false, false, false, 0));
  EXPECT_EQ(dbuff->getNumViews(), 1);

  // no-op, attaching a buffer to a view which already has a buffer
  otherbuffer = ds->createBuffer();
  dv->attachBuffer(otherbuffer);
  EXPECT_EQ(dv->getBuffer(), dbuff);  // same buffer as before
  EXPECT_EQ(otherbuffer->getNumViews(), 0);

  // The current attached buffer will be destroyed.
  dv->attachBuffer(nullptr);
  EXPECT_TRUE(checkViewValues(dv, EMPTY, false, false, false, 0));
  EXPECT_TRUE(ds->getBuffer(bindex) == nullptr);

  //---------- 4
  // Attach undescribed buffer to described view
  dv = root->createView("u4", INT_ID, BLEN);
  EXPECT_TRUE(checkViewValues(dv, EMPTY, true, false, false, BLEN));

  dbuff = ds->createBuffer();
  dv->attachBuffer(dbuff);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, false, false, BLEN));

  //---------- 5
  // Attach described buffer to described view
  dv = root->createView("u5", INT_ID, BLEN);
  EXPECT_TRUE(checkViewValues(dv, EMPTY, true, false, false, BLEN));

  dbuff = ds->createBuffer()->describe(INT_ID, BLEN);
  dv->attachBuffer(dbuff);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, false, false, BLEN));

  //---------- 6
  // Attach allocated buffer to described view
  dv = root->createView("u6", INT_ID, BLEN);
  EXPECT_TRUE(checkViewValues(dv, EMPTY, true, false, false, BLEN));

  dbuff = ds->createBuffer()->allocate(INT_ID, BLEN);
  dv->attachBuffer(dbuff);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));

  // Deallocate the buffer which will update the view.
  dbuff->deallocate();
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, false, false, BLEN));

  // Allocate the buffer via the view.
  dv->allocate();
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, true, BLEN));
  EXPECT_TRUE(dbuff->isAllocated());

  //---------- 7
  // Attach incompatable described buffer to described view
  dv = root->createView("u7", INT_ID, BLEN + 5);
  EXPECT_TRUE(checkViewValues(dv, EMPTY, true, false, false, BLEN + 5));

  dbuff = ds->createBuffer()->describe(INT_ID, BLEN);
  dv->attachBuffer(dbuff);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, false, false, BLEN + 5));

  //---------- 8
  // Attach incompatable allocated buffer to described view
  dv = root->createView("u8", INT_ID, BLEN + 5);
  EXPECT_TRUE(checkViewValues(dv, EMPTY, true, false, false, BLEN + 5));

  dbuff = ds->createBuffer()->allocate(INT_ID, BLEN);
  dv->attachBuffer(dbuff);
  EXPECT_TRUE(checkViewValues(dv, BUFFER, true, true, false, BLEN + 5));  // XXX -
                                                                          // how
                                                                          // is
  // isAllocated
  // useful

  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_view, int_array_strided_views)
{
  const IndexType num_elts = 10;
  const IndexType num_view_elts = 5;
  const IndexType elt_offset_e = 0;
  const IndexType elt_offset_o = 1;
  const IndexType elt_stride = 2;
  const IndexType elt_bytes = sizeof(int);

  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Buffer* dbuff = ds->createBuffer(INT_ID, num_elts);

  dbuff->allocate();
  int* data_ptr = dbuff->getData();

  for(int i = 0; i < num_elts; i++)
  {
    data_ptr[i] = i;
  }

  dbuff->print();

  EXPECT_EQ(dbuff->getTotalBytes(), num_elts * elt_bytes);

  View* dv_e = root->createView("even", dbuff);
  View* dv_o = root->createView("odd", dbuff);
  EXPECT_EQ(dbuff->getNumViews(), 2);

  // Set up views through conduit DataType interface
  // c_int(num_elems, offset [in bytes], stride [in bytes])
  dv_e->apply(DataType::c_int(num_view_elts, elt_offset_e * elt_bytes, elt_stride * elt_bytes));
  dv_o->apply(DataType::c_int(num_view_elts, elt_offset_o * elt_bytes, elt_stride * elt_bytes));

  // check that the view base pointers match that of the buffer
  EXPECT_EQ(dbuff->getVoidPtr(), dv_e->getVoidPtr());
  EXPECT_EQ(dbuff->getVoidPtr(), dv_o->getVoidPtr());

  dv_e->print();
  dv_o->print();

  // Check base pointer case:
  int* v_e_ptr = dv_e->getData();
  int* v_o_ptr = dv_o->getData();

  EXPECT_EQ(v_e_ptr, static_cast<int*>(dv_e->getVoidPtr()) + dv_e->getOffset());
  EXPECT_EQ(v_o_ptr, static_cast<int*>(dv_o->getVoidPtr()) + dv_o->getOffset());
  for(int i = 0; i < num_elts; i += 2)
  {
    std::cout << "idx:" << i << " e:" << v_e_ptr[i] << " o:" << v_o_ptr[i]
              << " em:" << v_e_ptr[i] % 2 << " om:" << v_o_ptr[i] % 2 << std::endl;

    EXPECT_EQ(v_e_ptr[i] % 2, 0);
    EXPECT_EQ(v_o_ptr[i] % 2, 1);
  }

  // Check Conduit mem-map struct case:
  int_array dv_e_ptr = dv_e->getData();
  int_array dv_o_ptr = dv_o->getData();
  for(int i = 0; i < 5; ++i)
  {
    std::cout << "idx:" << i << " e:" << dv_e_ptr[i] << " o:" << dv_o_ptr[i]
              << " em:" << dv_e_ptr[i] % 2 << " om:" << dv_o_ptr[i] % 2 << std::endl;

    EXPECT_EQ(dv_e_ptr[i] % 2, 0);
    EXPECT_EQ(dv_o_ptr[i] % 2, 1);
  }

  // Run similar test to above with different view apply method
  // Set up views through Sidre's native interface
  View* dv_e1 = root->createView("even1", dbuff);
  View* dv_o1 = root->createView("odd1", dbuff);
  EXPECT_EQ(dbuff->getNumViews(), 4);

  // (num_elems, offset [in # elems], stride [in # elems])
  dv_e1->apply(INT_ID, num_view_elts, elt_offset_e, elt_stride);
  dv_o1->apply(INT_ID, num_view_elts, elt_offset_o, elt_stride);

  // check that the view base pointers match that of the buffer
  EXPECT_EQ(dbuff->getVoidPtr(), dv_e1->getVoidPtr());
  EXPECT_EQ(dbuff->getVoidPtr(), dv_o1->getVoidPtr());

  dv_e1->print();
  dv_o1->print();

  // Check base pointer case:
  int* v_e1_ptr = dv_e1->getData();
  int* v_o1_ptr = dv_o1->getData();
  EXPECT_EQ(v_e1_ptr, static_cast<int*>(dv_e1->getVoidPtr()) + dv_e1->getOffset());
  EXPECT_EQ(v_o1_ptr, static_cast<int*>(dv_o1->getVoidPtr()) + dv_o1->getOffset());

  for(int i = 0; i < num_elts; i += 2)
  {
    std::cout << "idx:" << i << " e1:" << v_e1_ptr[i] << " oj:" << v_o1_ptr[i]
              << " em1:" << v_e1_ptr[i] % 2 << " om1:" << v_o1_ptr[i] % 2 << std::endl;

    EXPECT_EQ(v_e1_ptr[i], v_e_ptr[i]);
    EXPECT_EQ(v_o1_ptr[i], v_o_ptr[i]);
  }

  // Check Conduit mem-map struct case:
  int_array dv_e1_ptr = dv_e1->getData();
  int_array dv_o1_ptr = dv_o1->getData();
  for(int i = 0; i < 5; i++)
  {
    std::cout << "idx:" << i << " e1:" << dv_e1_ptr[i] << " o1:" << dv_o1_ptr[i]
              << " em1:" << dv_e1_ptr[i] % 2 << " om1:" << dv_o1_ptr[i] % 2 << std::endl;

    EXPECT_EQ(dv_e1_ptr[i], dv_e_ptr[i]);
    EXPECT_EQ(dv_o1_ptr[i], dv_o_ptr[i]);
  }

  ds->print();

  // Delete buffer and make sure views are no longer allocated.
  dbuff->deallocate();
  EXPECT_FALSE(dv_e->isAllocated());
  EXPECT_FALSE(dv_o->isAllocated());
  EXPECT_FALSE(dv_e1->isAllocated());
  EXPECT_FALSE(dv_o1->isAllocated());

  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_view, int_array_depth_view)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  const IndexType depth_nelems = 10;
  Buffer* dbuff = ds->createBuffer(INT_ID, 4 * depth_nelems);

  // Allocate buffer to hold data for 4 "depth" views
  dbuff->allocate();
  int* data_ptr = dbuff->getData();

  for(size_t i = 0; i < 4 * depth_nelems; ++i)
  {
    data_ptr[i] = i / depth_nelems;
  }

  dbuff->print();

  EXPECT_EQ(dbuff->getNumElements(), 4 * depth_nelems);

  // create 4 "depth" views and apply offsets into buffer
  View* views[4];
  std::string view_names[4] = {"depth_0", "depth_1", "depth_2", "depth_3"};

  for(int id = 0; id < 2; ++id)
  {
    views[id] = root->createView(view_names[id], dbuff)->apply(depth_nelems, id * depth_nelems);
  }
  //
  // call path including type
  for(int id = 2; id < 4; ++id)
  {
    views[id] =
      root->createView(view_names[id], dbuff)->apply(INT_ID, depth_nelems, id * depth_nelems);
  }
  EXPECT_EQ(dbuff->getNumViews(), 4);

  // Verify that the pointers for each view match that of the buffer
  for(int id = 0; id < 4; ++id)
  {
    void* vptr = views[id]->getVoidPtr();
    IndexType off = views[id]->getOffset();

    EXPECT_EQ(dbuff->getVoidPtr(), vptr);
    EXPECT_EQ(views[id]->getData<int*>(), static_cast<int*>(vptr) + off);
  }

  // print depth views...
  for(int id = 0; id < 4; ++id)
  {
    views[id]->print();
  }

  // check values in depth views...
  for(int id = 0; id < 4; ++id)
  {
    int* dv_ptr = views[id]->getData();
    for(IndexType i = 0; i < depth_nelems; ++i)
    {
      EXPECT_EQ(dv_ptr[i], id);
    }
  }

  ds->print();
  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_view, view_offset_and_stride)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  const IndexType nelems = 15;
  Buffer* dbuff = ds->createBuffer(DOUBLE_ID, nelems)->allocate();

  double* data_ptr = dbuff->getData();
  for(int i = 0; i < nelems; ++i)
  {
    data_ptr[i] = 1.01 * i;
  }

  // Array layout:
  //   [-,a,b,c,a,b,c,a,b,c,d,d,d,d,d]
  const int NUM_GROUPS = 4;
  int size[NUM_GROUPS] = {3, 3, 3, 5};
  int stride[NUM_GROUPS] = {3, 3, 3, 1};
  int offset[NUM_GROUPS] = {1, 2, 3, 10};
  std::string names[NUM_GROUPS] = {"a", "b", "c", "d"};
  const int ELT_SIZE = sizeof(double);

  // -- Test offsets and strides on buffer-based views

  Group* bufferGroup = root->createGroup("buffer");

  // Create the views off of the buffer
  for(int i = 0; i < NUM_GROUPS; ++i)
  {
    bufferGroup->createView(names[i], dbuff)->apply(size[i], offset[i], stride[i]);
  }

  // Test the sizes, offsets and strides of the views
  for(int i = 0; i < NUM_GROUPS; ++i)
  {
    View* view = bufferGroup->getView(names[i]);
    EXPECT_EQ(size[i], view->getNumElements());
    EXPECT_EQ(offset[i], view->getOffset());
    EXPECT_EQ(stride[i], view->getStride());
    EXPECT_EQ(ELT_SIZE, view->getBytesPerElement());
  }

  // test the offsets and strides applied to the data
  // Note: View::getData() already incorporates the offset
  //       of the view, View::getVoidPtr() does not
  for(int i = 0; i < NUM_GROUPS; ++i)
  {
    View* view = bufferGroup->getView(names[i]);
    void* raw_vp = view->getVoidPtr();
    double* offset_ptr = static_cast<double*>(raw_vp) + view->getOffset();
    double* getData_ptr = view->getData();

    EXPECT_EQ(data_ptr, raw_vp);
    EXPECT_EQ(offset_ptr, getData_ptr);
    EXPECT_EQ(stride[i], view->getStride());
  }

  // -- Test offsets and strides on external pointer based views
  Group* extGroup = root->createGroup("ext");

  // Create the views off of external data pointer
  for(int i = 0; i < NUM_GROUPS; ++i)
  {
    extGroup->createView(names[i], data_ptr)->apply(DOUBLE_ID, size[i], offset[i], stride[i]);
  }

  // Test the sizes, offsets and strides of the views
  for(int i = 0; i < NUM_GROUPS; ++i)
  {
    View* view = extGroup->getView(names[i]);
    EXPECT_EQ(size[i], view->getNumElements());
    EXPECT_EQ(offset[i], view->getOffset());
    EXPECT_EQ(stride[i], view->getStride());
    EXPECT_EQ(ELT_SIZE, view->getBytesPerElement());
  }

  // test the offsets and strides applied to the data
  // Note: View::getData() already incorporates the offset
  //       of the view, View::getVoidPtr() does not
  for(int i = 0; i < NUM_GROUPS; ++i)
  {
    View* view = extGroup->getView(names[i]);
    void* raw_vp = view->getVoidPtr();
    double* offset_ptr = static_cast<double*>(raw_vp) + view->getOffset();
    double* getData_ptr = view->getData();

    EXPECT_EQ(data_ptr, raw_vp);
    EXPECT_EQ(offset_ptr, getData_ptr);
    EXPECT_EQ(stride[i], view->getStride());
  }

  // -- Test offset and stride on the other view types:
  //          string, scalar, empty, opaque
  Group* othersGroup = root->createGroup("others");

  using ViewVec = std::vector<View*>;
  ViewVec views;
  std::uint8_t ui8 = 3;
  std::uint16_t ui16 = 4;
  std::uint32_t ui32 = 5;
#ifndef AXOM_NO_INT46_T
  std::uint64_t ui64 = 6;
#endif

  std::int8_t i8 = -3;
  std::int16_t i16 = -4;
  std::int32_t i32 = -5;
#ifndef AXOM_NO_INT46_T
  std::int64_t i64 = -6;
#endif

  axom::float32 f32 = 7.7f;
  axom::float64 f64 = 8.8;

  views.push_back(othersGroup->createView("key_empty"));
  // key_opaque is not described.
  views.push_back(othersGroup->createView("key_opaque", data_ptr));
  views.push_back(othersGroup->createViewString("key_string", "string_value"));

  views.push_back(othersGroup->createViewScalar("key_uint8", ui8));
  views.push_back(othersGroup->createViewScalar("key_uint16", ui16));
  views.push_back(othersGroup->createViewScalar("key_uint32", ui32));
#ifndef AXOM_NO_INT46_T
  views.push_back(othersGroup->createViewScalar("key_uint64", ui64));
#endif

  views.push_back(othersGroup->createViewScalar("key_int8", i8));
  views.push_back(othersGroup->createViewScalar("key_int16", i16));
  views.push_back(othersGroup->createViewScalar("key_int32", i32));
#ifndef AXOM_NO_INT46_T
  views.push_back(othersGroup->createViewScalar("key_int64", i64));
#endif

  views.push_back(othersGroup->createViewScalar("key_float32", f32));
  views.push_back(othersGroup->createViewScalar("key_float64", f64));

  std::map<std::string, int> sizeMap;
  sizeMap["key_empty"] = 0;
  sizeMap["key_opaque"] = 0;
  sizeMap["key_string"] = 1;
  sizeMap["key_uint8"] = sizeof(ui8);
  sizeMap["key_uint16"] = sizeof(ui16);
  sizeMap["key_uint32"] = sizeof(ui32);
#ifndef AXOM_NO_INT46_T
  sizeMap["key_uint64"] = sizeof(ui64);
#endif
  sizeMap["key_int8"] = sizeof(i8);
  sizeMap["key_int16"] = sizeof(i16);
  sizeMap["key_int32"] = sizeof(i32);
#ifndef AXOM_NO_INT46_T
  sizeMap["key_int64"] = sizeof(i64);
#endif
  sizeMap["key_float32"] = sizeof(f32);
  sizeMap["key_float64"] = sizeof(f64);

  for(ViewVec::iterator it = views.begin(); it != views.end(); ++it)
  {
    View* view = *it;
    EXPECT_EQ(0, view->getOffset());
    EXPECT_EQ(1, view->getStride());
    EXPECT_EQ(sizeMap[view->getName()], view->getBytesPerElement());
  }

  // -- cleanup

  ds->print();
  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_view, int_array_view_attach_buffer)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  const IndexType field_nelems = 10;

  // create 2 "field" views with type and # elems
  IndexType elem_count = 0;
  View* field0 = root->createView("field0", INT_ID, field_nelems);
  elem_count += field0->getNumElements();
  View* field1 = root->createView("field1", INT_ID, field_nelems);
  elem_count += field1->getNumElements();
  EXPECT_EQ(elem_count, 2 * field_nelems);

  // create buffer to hold data for fields and allocate
  Buffer* dbuff = ds->createBuffer()->allocate(INT_ID, elem_count);
  EXPECT_EQ(dbuff->getNumElements(), elem_count);

  // Initialize buffer data for testing below.
  int* b_ptr = dbuff->getData();
  for(IndexType i = 0; i < elem_count; ++i)
  {
    b_ptr[i] = i / field_nelems;
  }

  dbuff->print();

  // attach field views to buffer and apply offsets into buffer
  field0->attachBuffer(dbuff)->apply(field_nelems, 0 * field_nelems);
  field1->attachBuffer(dbuff)->apply(field_nelems, 1 * field_nelems);
  EXPECT_EQ(dbuff->getNumViews(), 2);

  // print field views...
  field0->print();
  field1->print();

  // check values in field views...
  int* f0_ptr = field0->getData();
  for(IndexType i = 0; i < field_nelems; ++i)
  {
    EXPECT_EQ(f0_ptr[i], 0);
  }
  int* f1_ptr = field1->getData();
  for(IndexType i = 0; i < field_nelems; ++i)
  {
    EXPECT_EQ(f1_ptr[i], 1);
  }

  ds->print();
  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_view, int_array_multi_view_resize)
{
  ///
  /// This example creates a 4 * 10 buffer of ints,
  /// and 4 views that point the 4 sections of 10 ints
  ///
  /// We then create a new buffer to support 4*12 ints
  /// and 4 views that point into them
  ///
  /// after this we use the old buffers to copy the values
  /// into the new views
  ///

  // create our main data store
  DataStore* ds = new DataStore();
  // get access to our root data Group
  Group* root = ds->getRoot();

  // create a group to hold the "old" or data we want to copy
  Group* r_old = root->createGroup("r_old");
  // create a view to hold the base buffer and allocate
  View* base_old = r_old->createViewAndAllocate("base_data", DataType::c_int(40));

  // we will create 4 sub views of this array
  int* data_ptr = base_old->getData();

  // init the buff with values that align with the
  // 4 subsections.
  for(int i = 0; i < 10; i++)
  {
    data_ptr[i] = 1;
  }
  for(int i = 10; i < 20; i++)
  {
    data_ptr[i] = 2;
  }
  for(int i = 20; i < 30; i++)
  {
    data_ptr[i] = 3;
  }
  for(int i = 30; i < 40; i++)
  {
    data_ptr[i] = 4;
  }

  /// setup our 4 views
  Buffer* buff_old = base_old->getBuffer();
  buff_old->print();
  View* r0_old = r_old->createView("r0", buff_old);
  View* r1_old = r_old->createView("r1", buff_old);
  View* r2_old = r_old->createView("r2", buff_old);
  View* r3_old = r_old->createView("r3", buff_old);

  // each view is offset by 10 * the # of bytes in a int
  // c_int(num_elems, offset)
  index_t offset = 0;
  r0_old->apply(DataType::c_int(10, offset));

  offset += sizeof(int) * 10;
  r1_old->apply(DataType::c_int(10, offset));

  offset += sizeof(int) * 10;
  r2_old->apply(DataType::c_int(10, offset));

  offset += sizeof(int) * 10;
  r3_old->apply(DataType::c_int(10, offset));

  /// check that our views actually point to the expected data
  //
  int* r0_ptr = r0_old->getData();
  for(int i = 0; i < 10; i++)
  {
    EXPECT_EQ(r0_ptr[i], 1);
    // check pointer relation
    EXPECT_EQ(&r0_ptr[i], &data_ptr[i]);
  }

  int* r3_ptr = r3_old->getData();
  for(int i = 0; i < 10; i++)
  {
    EXPECT_EQ(r3_ptr[i], 4);
    // check pointer relation
    EXPECT_EQ(&r3_ptr[i], &data_ptr[i + 30]);
  }

  // create a group to hold the "old" or data we want to copy into
  Group* r_new = root->createGroup("r_new");
  // create a view to hold the base buffer and allocate
  View* base_new = r_new->createViewAndAllocate("base_data", DataType::c_int(4 * 12));

  int* base_new_data = base_new->getData();
  for(int i = 0; i < 4 * 12; ++i)
  {
    base_new_data[i] = 0;
  }

  Buffer* buff_new = base_new->getBuffer();
  buff_new->print();

  // create the 4 sub views of this array
  View* r0_new = r_new->createView("r0", buff_new);
  View* r1_new = r_new->createView("r1", buff_new);
  View* r2_new = r_new->createView("r2", buff_new);
  View* r3_new = r_new->createView("r3", buff_new);

  // apply views to r0,r1,r2,r3
  // each view is offset by 12 * the # of bytes in a int

  // c_int(num_elems, offset)
  offset = 0;
  r0_new->apply(DataType::c_int(12, offset));

  offset += sizeof(int) * 12;
  r1_new->apply(DataType::c_int(12, offset));

  offset += sizeof(int) * 12;
  r2_new->apply(DataType::c_int(12, offset));

  offset += sizeof(int) * 12;
  r3_new->apply(DataType::c_int(12, offset));

  /// update r2 as an example first
  buff_new->print();
  r2_new->print();

  /// copy the subset of value
  r2_new->getNode().update(r2_old->getNode());
  r2_new->getNode().print();
  buff_new->print();

  /// check pointer values
  int* r2_new_ptr = r2_new->getData();

  for(int i = 0; i < 10; i++)
  {
    EXPECT_EQ(r2_new_ptr[i], 3);
  }

  for(int i = 10; i < 12; i++)
  {
    EXPECT_EQ(r2_new_ptr[i], 0);  // assumes zero-ed alloc
  }

  /// update the other views
  r0_new->getNode().update(r0_old->getNode());
  r1_new->getNode().update(r1_old->getNode());
  r3_new->getNode().update(r3_old->getNode());

  buff_new->print();

  ds->print();
  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_view, int_array_realloc)
{
  ///
  /// info
  ///

  // create our main data store
  DataStore* ds = new DataStore();
  // get access to our root data Group
  Group* root = ds->getRoot();

  // create a view to hold the base buffer
  View* a1 = root->createViewAndAllocate("a1", DataType::c_float(5));
  View* a2 = root->createViewAndAllocate("a2", DataType::c_int(5));

  float* a1_ptr = a1->getData();
  int* a2_ptr = a2->getData();

  for(int i = 0; i < 5; i++)
  {
    a1_ptr[i] = 5.0;
    a2_ptr[i] = -5;
  }

  EXPECT_EQ(a1->getTotalBytes(), static_cast<IndexType>(sizeof(float) * 5));
  EXPECT_EQ(a2->getTotalBytes(), static_cast<IndexType>(sizeof(int) * 5));

  a1->reallocate(DataType::c_float(10));
  a2->reallocate(DataType::c_int(15));

  a1_ptr = a1->getData();
  a2_ptr = a2->getData();

  for(int i = 0; i < 5; i++)
  {
    EXPECT_EQ(a1_ptr[i], 5.0);
    EXPECT_EQ(a2_ptr[i], -5);
  }

  for(int i = 5; i < 10; i++)
  {
    a1_ptr[i] = 10.0;
    a2_ptr[i] = -10;
  }

  for(int i = 10; i < 15; i++)
  {
    a2_ptr[i] = -15;
  }

  EXPECT_EQ(a1->getTotalBytes(), static_cast<IndexType>(sizeof(float) * 10));
  EXPECT_EQ(a2->getTotalBytes(), static_cast<IndexType>(sizeof(int) * 15));

  // Try some errors
  // XXX  a1->reallocate(DataType::c_int(20));

  ds->print();
  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_view, simple_opaque)
{
  // create our main data store
  DataStore* ds = new DataStore();
  // get access to our root data Group
  Group* root = ds->getRoot();
  int* src_data = new int[1];

  src_data[0] = 42;

  void* src_ptr = (void*)src_data;

  View* opq_view = root->createView("my_opaque", src_ptr);

  // External pointers are held in the view, should not have a buffer.
  EXPECT_EQ(ds->getNumBuffers(), 0u);

  EXPECT_TRUE(opq_view->isExternal());
  EXPECT_TRUE(!opq_view->isApplied());
  EXPECT_TRUE(opq_view->isOpaque());

  void* opq_ptr = opq_view->getVoidPtr();
  EXPECT_EQ(src_ptr, opq_ptr);

  int* out_data = (int*)opq_ptr;
  EXPECT_EQ(opq_ptr, src_ptr);
  EXPECT_EQ(out_data[0], 42);

  ds->print();
  delete ds;
  delete[] src_data;
}

//------------------------------------------------------------------------------
TEST(sidre_view, rename_view)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  View* v1 = root->createView("v_a");
  View* v2 = root->createView("v_b");
  View* v3 = root->createView("v_c");

  bool success = v1->rename("v_r");
  EXPECT_TRUE(success);
  EXPECT_TRUE(v1->getName() == "v_r");
  EXPECT_TRUE(root->hasView("v_r"));
  EXPECT_FALSE(root->hasView("v_a"));

  success = v2->rename("fields/v_s");
  EXPECT_FALSE(success);
  EXPECT_TRUE(v2->getName() == "v_b");

  success = v3->rename("v_b");
  EXPECT_FALSE(success);
  EXPECT_TRUE(v3->getName() == "v_c");

  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_datastore, destroy_buffer)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Buffer* dbuff1 = ds->createBuffer()->allocate(INT_ID, BLEN);
  View* view1a = root->createView("view1a", INT_ID, BLEN)->attachBuffer(dbuff1);
  View* view1b = root->createView("view1b", INT_ID, BLEN)->attachBuffer(dbuff1);

  EXPECT_TRUE(checkViewValues(view1a, BUFFER, true, true, true, BLEN));
  EXPECT_TRUE(checkViewValues(view1b, BUFFER, true, true, true, BLEN));

  // destroyBuffer will detach from views.
  IndexType bindex = dbuff1->getIndex();
  ds->destroyBuffer(bindex);
  EXPECT_TRUE(ds->getBuffer(bindex) == nullptr);

  // views no longer have buffers (but retain descriptions)
  EXPECT_TRUE(checkViewValues(view1a, EMPTY, true, false, false, BLEN));
  EXPECT_TRUE(checkViewValues(view1b, EMPTY, true, false, false, BLEN));

  EXPECT_FALSE(view1a->hasBuffer());
  EXPECT_FALSE(view1b->hasBuffer());

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_view, value_from_uninited_view)
{
  // Note: This test relies on re-wiring conduit error handlers
  DataStore::setConduitSLICMessageHandlers();

  DataStore ds;
  View* view = ds.getRoot()->createView("empty");

  // check getScalar
  int val = view->getScalar();
  EXPECT_EQ(val, 0);

  // check getArray
  int* aval_ptr = view->getArray();
  EXPECT_TRUE(aval_ptr == nullptr);

  int aval = view->getArray();
  EXPECT_EQ(aval, 0);

  // check getData
  int* dval_ptr = view->getData();
  EXPECT_TRUE(dval_ptr == nullptr);

  int dval = view->getData();
  EXPECT_EQ(dval, 0);

  // restore conduit default errors
  DataStore::setConduitDefaultMessageHandlers();
}

//------------------------------------------------------------------------------
TEST(sidre_view, import_array_node)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  //Import Node holding int array
  std::vector<int> int_vec;
  for(int i = 0; i < 10; ++i)
  {
    int_vec.push_back(i * i + 3);
  }

  conduit::Node n_ints;
  n_ints.set(int_vec);

  View* v1 = root->createView("v1");
  v1->importArrayNode(n_ints);

  EXPECT_TRUE(v1->hasBuffer());
  EXPECT_TRUE(v1->isAllocated());
  EXPECT_TRUE(v1->isApplied());
  EXPECT_EQ(v1->getNumElements(), 10);

  int* v_ints = v1->getData();

  for(int i = 0; i < 10; ++i)
  {
    EXPECT_EQ(v_ints[i], i * i + 3);
  }

  //Import Node holding double array
  std::vector<double> dbl_vec;
  for(int i = 0; i < 8; ++i)
  {
    dbl_vec.push_back(static_cast<double>(i) / 2.0 + 1.1);
  }

  conduit::Node n_dbls;
  n_dbls.set(dbl_vec);

  View* v2 = root->createView("v2");
  v2->importArrayNode(n_dbls);

  EXPECT_TRUE(v2->hasBuffer());
  EXPECT_TRUE(v2->isAllocated());
  EXPECT_TRUE(v2->isApplied());
  EXPECT_EQ(v2->getNumElements(), 8);

  double* v_dbls = v2->getData();

  for(int i = 0; i < 8; ++i)
  {
    EXPECT_NEAR(v_dbls[i], static_cast<double>(i) / 2.0 + 1.1, 1.0e-12);
  }

  //Attempt to import a non-array node, will not work.
  conduit::Node n_obj;
  n_obj.set(conduit::DataType::object());

  View* v3 = root->createView("v3");
  v3->importArrayNode(n_obj);

  EXPECT_TRUE(v3->isEmpty());

  //Attempt to import into View that is already a string, will not work.
  View* v4 = root->createView("v4");
  v4->setString("string_view");

  v4->importArrayNode(n_ints);
  EXPECT_TRUE(v4->isString());

  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_view, clear_view)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  // Create an empty view.
  {
    View* view = root->createView("v_empty");
    EXPECT_TRUE(checkViewValues(view, EMPTY, false, false, false, 0));
    view->clear();
    EXPECT_TRUE(checkViewValues(view, EMPTY, false, false, false, 0));
  }

  // Describe an empty view.
  {
    View* view = root->createView("v_described", INT_ID, BLEN);
    EXPECT_TRUE(checkViewValues(view, EMPTY, true, false, false, BLEN));
    view->clear();
    EXPECT_TRUE(checkViewValues(view, EMPTY, false, false, false, 0));
  }

  // scalar view.
  {
    View* view = root->createViewScalar("v_scalar", 1);
    EXPECT_TRUE(checkViewValues(view, TUPLE, true, true, true, 1));
    view->clear();
    EXPECT_TRUE(checkViewValues(view, EMPTY, false, false, false, 0));
  }

  // string view.
  {
    View* view = root->createViewString("v_string", "string-test");
    view->clear();
    EXPECT_TRUE(checkViewValues(view, EMPTY, false, false, false, 0));
  }

  // Allocated view, Buffer will be released.
  {
    int nbuf = ds->getNumBuffers();
    View* view = root->createViewAndAllocate("v_allocated", INT_ID, BLEN);
    EXPECT_TRUE(checkViewValues(view, BUFFER, true, true, true, BLEN));
    view->clear();
    EXPECT_TRUE(checkViewValues(view, EMPTY, false, false, false, 0));
    EXPECT_EQ(ds->getNumBuffers(), nbuf);
  }

  // Undescribed Buffer.
  {
    int nbuf = ds->getNumBuffers();
    Buffer* dbuff = ds->createBuffer();
    View* view = root->createView("v_undescribed_buffer", dbuff);
    EXPECT_TRUE(checkViewValues(view, BUFFER, false, false, false, 0));
    view->clear();
    EXPECT_TRUE(checkViewValues(view, EMPTY, false, false, false, 0));
    EXPECT_EQ(ds->getNumBuffers(), nbuf);
  }

  // Explicit Buffer attached to two views.
  // Buffer will not be released.
  {
    Buffer* dbuff = ds->createBuffer()->allocate(INT_ID, BLEN);
    int nbuf = ds->getNumBuffers();
    EXPECT_EQ(dbuff->getNumViews(), 0);

    View* vother = root->createView("v_other", INT_ID, BLEN);
    View* view = root->createView("v_buffer", INT_ID, BLEN);
    vother->attachBuffer(dbuff);
    EXPECT_EQ(dbuff->getNumViews(), 1);
    view->attachBuffer(dbuff);
    EXPECT_EQ(dbuff->getNumViews(), 2);

    EXPECT_TRUE(checkViewValues(view, BUFFER, true, true, true, BLEN));
    view->clear();
    EXPECT_TRUE(checkViewValues(view, EMPTY, false, false, false, 0));

    EXPECT_EQ(ds->getNumBuffers(), nbuf);
    EXPECT_EQ(dbuff->getNumViews(), 1);
  }

  // External view.
  {
    int extdata[BLEN];
    View* view = root->createView("v_external", INT_ID, BLEN, &extdata);
    EXPECT_TRUE(checkViewValues(view, EXTERNAL, true, true, true, BLEN));
    view->clear();
    EXPECT_TRUE(checkViewValues(view, EMPTY, false, false, false, 0));
  }

  // Opaque view.
  {
    int extdata[BLEN];
    View* view = root->createView("v_opaque", &extdata);
    EXPECT_TRUE(checkViewValues(view, EXTERNAL, false, true, false, 0));
    view->clear();
    EXPECT_TRUE(checkViewValues(view, EMPTY, false, false, false, 0));
  }

  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_view, deep_copy_shape)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Group* groupA = root->createGroup("groupA");
  Group* groupB = root->createGroup("groupB");

  {
    const IndexType shapeA[3] = {3, 2, 5};
    View* viewA = groupA->createViewWithShape("array3d", INT_ID, 3, shapeA);
    View* viewB = groupB->deepCopyView(viewA);
    IndexType shapeB[3] = {-1, -1, -1};
    viewB->getShape(3, shapeB);
    EXPECT_EQ(viewB->getNumDimensions(), 3);
    EXPECT_EQ(shapeB[0], shapeA[0]);
    EXPECT_EQ(shapeB[1], shapeA[1]);
    EXPECT_EQ(shapeB[2], shapeA[2]);
  }

  delete ds;
}

// Return vector of known allocator ids.
std::vector<int> getKnownAllocIds()
{
  std::vector<int> allocIds(1, axom::MALLOC_ALLOCATOR_ID);
#ifdef AXOM_USE_UMPIRE
  allocIds.push_back(axom::detail::getAllocatorID<axom::MemorySpace::Host>());
  #ifdef AXOM_USE_GPU
  allocIds.push_back(axom::detail::getAllocatorID<axom::MemorySpace::Device>());
  allocIds.push_back(axom::detail::getAllocatorID<axom::MemorySpace::Unified>());
  // Does it make sense to check Pinned and Constant memory spaces?
  #endif
#endif
  return allocIds;
}

//------------------------------------------------------------------------------

TEST(sidre_view, deep_copy_to_conduit)
{
  std::vector<int> allocIds = getKnownAllocIds();

  DataStore ds;

  constexpr int N = 5;
  auto dtypeDouble = conduit::DataType::float64(1);
  auto dtypeIntArray = conduit::DataType::int32(N);

  const double doubleValue = 10.13456;
  axom::Array<std::int32_t> intArray(N, N);
  for(int i = 0; i < N; ++i)
  {
    intArray[i] = 1001.5 + i;
  }
  std::string stringValue = "a string";

  double tmpDoubleValue = doubleValue;
  axom::Array<std::int32_t> tmpIntArray(N, N);
  axom::Array<char> tmpCharArray;

  // For each combination of srcArrayAllocId and dstAllocId,
  // allocate src and copy to dst.

  for(size_t si = 0; si < allocIds.size(); ++si)
  {
    int srcArrayAllocId = allocIds[si];
    int srcTupleAllocId = allocIds[allocIds.size() - 1 - si];

    Group* src = ds.getRoot()->createGroup("src");
    src->setDefaultArrayAllocator(srcArrayAllocId);
    src->setDefaultTupleAllocator(srcArrayAllocId);

    //
    // Create, initialize and sanity-check src objects to test.
    //

    View* srcString = src->createViewString("aString", stringValue);
    const char* srcStringPtr = static_cast<char*>(srcString->getNode().data_ptr());
    if(srcArrayAllocId == axom::MALLOC_ALLOCATOR_ID)
    {
      // Skip this check if srcArrayAllocId is non-MALLOC_ALLOCATOR_ID,
      // because a current sidre issue results in always
      // placing strings on Conduit's default alloc id.
      EXPECT_EQ(axom::getAllocatorIDFromPointer(srcStringPtr), srcArrayAllocId);
    }

    View* srcScalar = src->createViewScalar("aDouble", dtypeDouble);
    double* srcScalarPtr = static_cast<double*>(srcScalar->getNode().data_ptr());
    EXPECT_EQ(axom::getAllocatorIDFromPointer(srcScalarPtr), srcArrayAllocId);
    axom::copy(srcScalarPtr, &doubleValue, sizeof(double));

    View* srcArray = src->createViewAndAllocate("aIntArray", dtypeIntArray);
    std::int32_t* srcArrayPtr = static_cast<std::int32_t*>(srcArray->getNode().data_ptr());
    EXPECT_EQ(axom::getAllocatorIDFromPointer(srcArrayPtr), srcArrayAllocId);
    axom::copy(srcArrayPtr, intArray.data(), N * sizeof(std::int32_t));

    if(axom::execution_space<axom::SEQ_EXEC>::usesAllocId(srcArrayAllocId))
    {
      std::cout << "src group:" << std::endl;
      src->print();
    }

    //
    // Copy the source into destinations of different alloc ids
    // Conduit doesn't support separate allocators for arrays and tuples
    // so there's only one destination allocator.
    //

    for(size_t di = 0; di < allocIds.size(); ++di)
    {
      int dstAllocId = allocIds[di];

      std::cout << "Testing copying allocator ids " << srcArrayAllocId << " and " << srcTupleAllocId
                << " to " << dstAllocId << std::endl;

      const auto& idConverter = axom::sidre::ConduitMemory::instanceForAxomId(dstAllocId);
      auto dstAllocIdConduit = idConverter.conduitId();

      conduit::Node dst;
      dst.set_allocator(dstAllocIdConduit);

      src->deepCopyToConduit(dst);

      //
      // Check pointers.  Copy data to temporary host buffers and check data.
      //
      int checkAllocId;

      double* dstScalarPtr =
        static_cast<double*>(dst.fetch_existing(srcScalar->getName()).data_ptr());
      EXPECT_NE(dstScalarPtr, nullptr);
      EXPECT_NE(dstScalarPtr, srcScalarPtr);
      checkAllocId = axom::getAllocatorIDFromPointer(dstScalarPtr);
      EXPECT_EQ(checkAllocId, dstAllocId);
      axom::copy(&tmpDoubleValue, dstScalarPtr, sizeof(double));
      EXPECT_EQ(tmpDoubleValue, doubleValue);

      std::int32_t* dstArrayPtr =
        static_cast<std::int32_t*>(dst.fetch_existing(srcArray->getName()).data_ptr());
      EXPECT_NE(dstArrayPtr, nullptr);
      EXPECT_NE(dstArrayPtr, srcArrayPtr);
      checkAllocId = axom::getAllocatorIDFromPointer(dstArrayPtr);
      EXPECT_EQ(checkAllocId, dstAllocId);
      axom::copy(tmpIntArray.data(), srcArrayPtr, N * sizeof(std::int32_t));
      for(int i = 0; i < N; ++i)
      {
        EXPECT_EQ(tmpIntArray[i], intArray[i]);
      }

      char* dstStringPtr = static_cast<char*>(dst.fetch_existing(srcString->getName()).data_ptr());
      EXPECT_NE(dstStringPtr, nullptr);
      EXPECT_NE(dstStringPtr, srcStringPtr);
      checkAllocId = axom::getAllocatorIDFromPointer(dstStringPtr);
      EXPECT_EQ(checkAllocId, dstAllocId);
      tmpCharArray.resize(srcString->getNumElements());
      axom::copy(tmpCharArray.data(), srcStringPtr, srcString->getNumElements());
      for(int i = 0; i < N; ++i)
      {
        EXPECT_EQ(tmpCharArray[i], stringValue[i]);
      }
    }

    ds.getRoot()->destroyGroup("src");
  }
}

//------------------------------------------------------------------------------

TEST(sidre_view, deep_copy_to_conduit_with_allocid)
{
  using axom::sidre::ConduitMemory;

  std::vector<int> allocIds = getKnownAllocIds();
  for(std::size_t si = 0; si < allocIds.size(); ++si)
  {
    auto aId = allocIds[si];
    auto cId = ConduitMemory::axomAllocIdToConduit(aId);
    std::cout << "si=" << si << "  axomAllocId=" << aId << "  conduitAllocId=" << cId << std::endl;
  }
  if(allocIds.size() < 2)
  {
    std::cout << "Skipping test deep_copy_to_conduit_with_allocid."
                 "At least 2 allocators are needed to test."
              << std::endl;
    return;
  }

  DataStore ds;

  Group* srcGrp = ds.getRoot()->createGroup("src");
  View* srcView = srcGrp->createView("a");
  srcView->setScalar(int(15));

  // Copy srcView to conduit destinations, varying the allocator.
  for(std::size_t si = 0; si < allocIds.size(); ++si)
  {
    // dst1 uses default Conduit allocator
    conduit::Node dst1;
    auto defaultConduitAllocId = dst1.allocator();
    srcView->deepCopyToConduit(dst1);
    EXPECT_EQ(dst1.allocator(), defaultConduitAllocId);

    // dst2 sets and uses allocator allocIds[si]
    conduit::Node dst2;
    dst2.set_allocator(ConduitMemory::axomAllocIdToConduit(allocIds[si]));
    srcView->deepCopyToConduit(dst2);
    EXPECT_EQ(ConduitMemory::conduitAllocIdToAxom(dst2.allocator()), allocIds[si]);

    // dst3 sets allocator allocIds[si] but overrides it in the copy
    conduit::Node dst3;
    dst3.set_allocator(ConduitMemory::axomAllocIdToConduit(allocIds[si]));
    srcView->deepCopyToConduit(dst3, allocIds[0]);
    EXPECT_EQ(ConduitMemory::conduitAllocIdToAxom(dst3.allocator()), allocIds[0]);
  }

  ds.getRoot()->destroyGroup("src");
}

//------------------------------------------------------------------------------

TEST(sidre_view, reshape_array)
{
  using namespace axom;
  DataStore ds;
  Group* root = ds.getRoot();

  constexpr int DMAX = 4;
  IndexType shapeOutput[DMAX];
  int nDim = -1;

  // A view with array in a buffer.
  auto* viewA = root->createView("viewA", sidre::detail::SidreTT<int32_t>::id, 12);

  {
    nDim = viewA->getShape(DMAX, shapeOutput);
    EXPECT_EQ(viewA->getNumElements(), 12);
    EXPECT_EQ(nDim, 1);
    EXPECT_EQ(shapeOutput[0], 12);
  }

  {
    // Attempt to reshape empty View should be no-op.
    EXPECT_TRUE(viewA->isEmpty());
    IndexType badShape[] = {1, 2, 3};
    SLIC_INFO(
      "Next warning about reshaping non-array view is expected and can be "
      "ignored.");
    viewA->reshapeArray(3, badShape);
    nDim = viewA->getShape(DMAX, shapeOutput);
    EXPECT_EQ(viewA->getNumElements(), 12);
    EXPECT_EQ(nDim, 1);
    EXPECT_EQ(shapeOutput[0], 12);
  }

  viewA->allocate();
  EXPECT_FALSE(viewA->isEmpty());
  EXPECT_FALSE(viewA->isScalar());
  EXPECT_TRUE(viewA->hasBuffer());

  {
    // Attempt to change size should be no-op.
    IndexType badShape[] = {1, 2, 3};
    SLIC_INFO(
      "Next warning about changing number of elemnents is expected and can be "
      "ignored.");
    viewA->reshapeArray(3, badShape);
    nDim = viewA->getShape(DMAX, shapeOutput);
    EXPECT_EQ(viewA->getNumElements(), 12);
    EXPECT_EQ(nDim, 1);
    EXPECT_EQ(shapeOutput[0], 12);
  }

  {
    IndexType shape2d[] = {3, 4};
    viewA->reshapeArray(2, shape2d);
    nDim = viewA->getShape(DMAX, shapeOutput);
    EXPECT_EQ(viewA->getNumElements(), 12);
    EXPECT_EQ(nDim, 2);
    EXPECT_EQ(shapeOutput[0], 3);
    EXPECT_EQ(shapeOutput[1], 4);
  }

  viewA->deallocate();
  EXPECT_FALSE(viewA->isAllocated());

  {
    // Reshaping an unallocated array is allowed, as long as it is described.
    IndexType shape2d[] = {2, 6};
    EXPECT_TRUE(viewA->isDescribed());
    viewA->reshapeArray(2, shape2d);
    nDim = viewA->getShape(DMAX, shapeOutput);
    EXPECT_EQ(viewA->getNumElements(), 12);
    EXPECT_EQ(nDim, 2);
    EXPECT_EQ(shapeOutput[0], 2);
    EXPECT_EQ(shapeOutput[1], 6);
  }

  viewA->allocate();
  EXPECT_TRUE(viewA->isAllocated());

  {
    IndexType shape3d[] = {3, 2, 2};
    viewA->reshapeArray(3, shape3d);
    nDim = viewA->getShape(DMAX, shapeOutput);
    EXPECT_EQ(viewA->getNumElements(), 12);
    EXPECT_EQ(nDim, 3);
    EXPECT_EQ(shapeOutput[0], 3);
    EXPECT_EQ(shapeOutput[1], 2);
    EXPECT_EQ(shapeOutput[2], 2);
  }

  {
    IndexType shape1d = viewA->getNumElements();
    viewA->reshapeArray(1, &shape1d);
    nDim = viewA->getShape(DMAX, shapeOutput);
    EXPECT_EQ(viewA->getNumElements(), shape1d);
    EXPECT_EQ(nDim, 1);
    EXPECT_EQ(shapeOutput[0], shape1d);
    // Test a valid reshape that doesn't change the shape.
    viewA->reshapeArray(1, &shape1d);
    nDim = viewA->getShape(DMAX, shapeOutput);
    EXPECT_EQ(viewA->getNumElements(), shape1d);
    EXPECT_EQ(nDim, 1);
    EXPECT_EQ(shapeOutput[0], shape1d);
  }

  // A view with external array data.
  auto* viewB = root->createView("viewB");

  std::int32_t extData[24];

  {
    viewB->setExternalDataPtr(sidre::detail::SidreTT<std::int32_t>::id, 24, extData);
    EXPECT_TRUE(viewB->isExternal());
    nDim = viewB->getShape(DMAX, shapeOutput);
    EXPECT_EQ(viewB->getNumElements(), 24);
    EXPECT_EQ(nDim, 1);
    EXPECT_EQ(shapeOutput[0], 24);
  }

  {
    IndexType shape4d[] = {1, 2, 3, 4};
    viewB->reshapeArray(4, shape4d);
    nDim = viewB->getShape(DMAX, shapeOutput);
    EXPECT_EQ(viewB->getNumElements(), 24);
    EXPECT_EQ(nDim, 4);
    EXPECT_EQ(shapeOutput[0], 1);
    EXPECT_EQ(shapeOutput[1], 2);
    EXPECT_EQ(shapeOutput[2], 3);
    EXPECT_EQ(shapeOutput[3], 4);
  }
}

//------------------------------------------------------------------------------

TEST(sidre_view, apply_and_shape_behavior)
{
  using namespace axom;
  DataStore ds;
  Group* root = ds.getRoot();

  const int nfoo = 10;
  int int2d[nfoo * 2];
  IndexType shape[] = {nfoo, 2};

  View* view1 = root->createViewWithShape("int2d", INT_ID, 2, shape, int2d);

  IndexType shapeOutput[] = {0, 0};

  view1->getShape(2, shapeOutput);
  EXPECT_EQ(shapeOutput[0], nfoo);
  EXPECT_EQ(shapeOutput[1], 2);
  EXPECT_EQ(view1->getNumDimensions(), 2);
  EXPECT_EQ(view1->getStride(), 1);

  // 20 elements, offset of 0, stride of 2
  // Note - dimensions are flattened to 1D
  view1->apply(nfoo * 2, 0, 2);
  view1->getShape(2, shapeOutput);
  EXPECT_EQ(shapeOutput[0], nfoo * 2);
  EXPECT_EQ(shapeOutput[1], 0);
  EXPECT_EQ(view1->getNumDimensions(), 1);
  EXPECT_EQ(view1->getStride(), 2);

  // Reshape back to [10,2]
  view1->reshapeArray(2, shape);
  view1->getShape(2, shapeOutput);
  EXPECT_EQ(shapeOutput[0], nfoo);
  EXPECT_EQ(shapeOutput[1], 2);
  EXPECT_EQ(view1->getNumDimensions(), 2);
  EXPECT_EQ(view1->getStride(), 1);

  // apply() (no args) - no change
  view1->apply();
  view1->getShape(2, shapeOutput);
  EXPECT_EQ(shapeOutput[0], nfoo);
  EXPECT_EQ(shapeOutput[1], 2);
  EXPECT_EQ(view1->getNumDimensions(), 2);
  EXPECT_EQ(view1->getStride(), 1);
}

//------------------------------------------------------------------------------

#ifdef AXOM_USE_UMPIRE

class UmpireTest : public ::testing::TestWithParam<int>
{
public:
  void SetUp() override
  {
    allocID = GetParam();
    root = ds.getRoot();
  }

  void TearDown() override
  {
    const int allocID = axom::getUmpireResourceAllocatorID(umpire::resource::Host);
    axom::setDefaultAllocator(allocID);
  }

  static constexpr int SIZE = 100;
  DataStore ds;
  Group* root;
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  int allocID;
};

//------------------------------------------------------------------------------
TEST_P(UmpireTest, allocate)
{
  {
    View* view = root->createView("v");
    view->allocate(INT_ID, SIZE, allocID);

    ASSERT_EQ(allocID, rm.getAllocator(view->getVoidPtr()).getId());
    root->destroyViewAndData("v");
  }

  {
    View* view = root->createView("v");
    DataType dtype = conduit::DataType::default_dtype(INT_ID);
    dtype.set_number_of_elements(SIZE);
    view->allocate(dtype, allocID);

    ASSERT_EQ(allocID, rm.getAllocator(view->getVoidPtr()).getId());
    root->destroyViewAndData("v");
  }
}

//------------------------------------------------------------------------------
TEST_P(UmpireTest, allocate_default)
{
  root->setDefaultAllocator(allocID);

  {
    View* view = root->createView("v");
    view->allocate(INT_ID, SIZE);

    ASSERT_EQ(allocID, rm.getAllocator(view->getVoidPtr()).getId());
    root->destroyViewAndData("v");
  }

  {
    View* view = root->createView("v");
    DataType dtype = conduit::DataType::default_dtype(INT_ID);
    dtype.set_number_of_elements(SIZE);
    view->allocate(dtype);

    ASSERT_EQ(allocID, rm.getAllocator(view->getVoidPtr()).getId());
    root->destroyViewAndData("v");
  }
}

//------------------------------------------------------------------------------
TEST_P(UmpireTest, reallocate)
{
  #if defined(AXOM_USE_GPU) && defined(UMPIRE_ENABLE_CONST)
  if(allocID == axom::getUmpireResourceAllocatorID(umpire::resource::Constant))
  {
    return;
  }
  #endif

  {
    View* view = root->createView("v");
    view->allocate(INT_ID, SIZE, allocID);
    view->reallocate(2 * SIZE);

    ASSERT_EQ(allocID, rm.getAllocator(view->getVoidPtr()).getId());
    root->destroyViewAndData("v");
  }

  {
    View* view = root->createView("v");
    DataType dtype = conduit::DataType::default_dtype(INT_ID);
    dtype.set_number_of_elements(SIZE);
    view->allocate(dtype, allocID);
    view->reallocate(2 * SIZE);

    ASSERT_EQ(allocID, rm.getAllocator(view->getVoidPtr()).getId());
    root->destroyViewAndData("v");
  }
}

//------------------------------------------------------------------------------
TEST_P(UmpireTest, reallocate_zero)
{
  #if defined(AXOM_USE_GPU) && defined(UMPIRE_ENABLE_CONST)
  if(allocID == axom::getUmpireResourceAllocatorID(umpire::resource::Constant))
  {
    return;
  }
  #endif

  axom::setDefaultAllocator(allocID);

  {
    View* view = root->createView("v");
    view->allocate(INT_ID, SIZE, allocID);
    view->reallocate(0);
    view->reallocate(SIZE);

    ASSERT_EQ(axom::getDefaultAllocatorID(), rm.getAllocator(view->getVoidPtr()).getId());

    root->destroyViewAndData("v");
  }

  {
    View* view = root->createView("v");
    DataType dtype = conduit::DataType::default_dtype(INT_ID);
    dtype.set_number_of_elements(SIZE);
    view->allocate(dtype, allocID);
    view->reallocate(0);
    view->reallocate(SIZE);

    ASSERT_EQ(axom::getDefaultAllocatorID(), rm.getAllocator(view->getVoidPtr()).getId());

    root->destroyViewAndData("v");
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

INSTANTIATE_TEST_SUITE_P(sidre_view, UmpireTest, ::testing::ValuesIn(allocators));

#endif  // AXOM_USE_UMPIRE

TEST(sidre_view, isUpdateableFrom)
{
  constexpr int SIZE = 100;
  DataStore ds;
  Group* root = ds.getRoot();
  View* host = root->createViewAndAllocate("v", INT_ID, SIZE);

  // Cannot update from an empty View
  {
    View* v = root->createView("v0");
    ASSERT_FALSE(host->isUpdateableFrom(v));
    ASSERT_FALSE(v->isUpdateableFrom(host));
  }

  // Cannot update from a scalar View
  {
    View* v = root->createViewScalar("v1", 5);
    ASSERT_FALSE(host->isUpdateableFrom(v));
    ASSERT_FALSE(v->isUpdateableFrom(host));
  }

  // Cannot update from a string View
  {
    View* v = root->createViewString("v2", "dummy");
    ASSERT_FALSE(host->isUpdateableFrom(v));
    ASSERT_FALSE(v->isUpdateableFrom(host));
  }

  // Cannot update from a View with a different number of bytes
  {
    View* v = root->createViewAndAllocate("v3", INT_ID, SIZE + 1);
    ASSERT_FALSE(host->isUpdateableFrom(v));
    ASSERT_FALSE(v->isUpdateableFrom(host));
  }

  // Cannot update from a View with a non-unit stride.
  {
    View* v = root->createViewAndAllocate("v4", INT_ID, SIZE);
    v->apply(SIZE / 2, 0, 2);
    ASSERT_FALSE(host->isUpdateableFrom(v));
    ASSERT_FALSE(v->isUpdateableFrom(host));
  }

  // Can update from a simlar view.
  {
    View* v = root->createViewAndAllocate("v5", INT_ID, SIZE);
    ASSERT_TRUE(host->isUpdateableFrom(v));
    ASSERT_TRUE(v->isUpdateableFrom(host));
  }

  // Can update from a simlar view with offset.
  {
    View* v = root->createViewAndAllocate("v6", INT_ID, SIZE + 10);
    v->apply(SIZE, 10);
    ASSERT_TRUE(host->isUpdateableFrom(v));
    ASSERT_TRUE(v->isUpdateableFrom(host));
  }

  // Can update from a view with a different type but same number of bytes.
  {
    View* v = root->createViewAndAllocate("v7", TypeID::INT8_ID, SIZE * sizeof(int));
    ASSERT_TRUE(host->isUpdateableFrom(v));
    ASSERT_TRUE(v->isUpdateableFrom(host));
  }
}

class UpdateTest : public ::testing::TestWithParam<::testing::tuple<std::string, std::string, int>>
{
public:
  void SetUp() override
  {
    src_string = ::testing::get<0>(GetParam());
    dst_string = ::testing::get<1>(GetParam());
    offset = ::testing::get<2>(GetParam());
    int size = SIZE - offset;

    Group* root = ds.getRoot();
    host = root->createViewAndAllocate("host", INT_ID, size);

    if(src_string == "NEW")
    {
      m_src_array = new int[SIZE];
      src = root->createView("src")->setExternalDataPtr(m_src_array)->apply(INT_ID, size, offset);
    }
    else if(src_string == "MALLOC")
    {
      m_src_array = static_cast<int*>(std::malloc(SIZE * sizeof(int)));
      src = root->createView("src")->setExternalDataPtr(m_src_array)->apply(INT_ID, size, offset);
    }
    else if(src_string == "STATIC")
    {
      src =
        root->createView("src")->setExternalDataPtr(m_static_src_array)->apply(INT_ID, size, offset);
    }
#ifdef AXOM_USE_UMPIRE
    else
    {
      int src_alloc_id = rm.getAllocator(src_string).getId();
      src = root->createViewAndAllocate("src", INT_ID, SIZE, src_alloc_id)->apply(size, offset);
    }
#endif

    if(dst_string == "NEW")
    {
      m_dst_array = new int[SIZE];
      dst = root->createView("dst")->setExternalDataPtr(m_dst_array)->apply(INT_ID, size, offset);
    }
    else if(dst_string == "MALLOC")
    {
      m_dst_array = static_cast<int*>(std::malloc(SIZE * sizeof(int)));
      dst = root->createView("dst")->setExternalDataPtr(m_dst_array)->apply(INT_ID, size, offset);
    }
    else if(dst_string == "STATIC")
    {
      dst =
        root->createView("dst")->setExternalDataPtr(m_static_dst_array)->apply(INT_ID, size, offset);
    }
#ifdef AXOM_USE_UMPIRE
    else
    {
      int dst_alloc_id = rm.getAllocator(dst_string).getId();
      dst = root->createViewAndAllocate("dst", INT_ID, SIZE, dst_alloc_id)->apply(size, offset);
    }
#endif
  }

  void TearDown() override
  {
    if(src_string == "NEW")
    {
      delete[] m_src_array;
    }
    else if(src_string == "MALLOC")
    {
      std::free(m_src_array);
    }

    if(dst_string == "NEW")
    {
      delete[] m_dst_array;
    }
    else if(dst_string == "MALLOC")
    {
      std::free(m_dst_array);
    }
  }

  static constexpr int SIZE = 100;

#ifdef AXOM_USE_UMPIRE
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
#endif

  std::string src_string;
  std::string dst_string;
  int offset;

  DataStore ds;
  View* host;
  View* src;
  View* dst;

private:
  int m_static_src_array[SIZE];
  int m_static_dst_array[SIZE];
  int* m_src_array = nullptr;
  int* m_dst_array = nullptr;
};

TEST_P(UpdateTest, updateFrom)
{
  std::cout << "SRC = " << src_string << ", DST = " << dst_string << ", OFFSET = " << offset
            << std::endl;

  ASSERT_TRUE(src->isUpdateableFrom(host));
  ASSERT_TRUE(dst->isUpdateableFrom(src));
  ASSERT_TRUE(host->isUpdateableFrom(dst));

  int* host_ptr = host->getData();
  const int size = host->getNumElements();
  ASSERT_EQ(size, SIZE - offset);

  for(int i = 0; i < size; ++i)
  {
    host_ptr[i] = i;
  }

  src->updateFrom(host);

  for(int i = 0; i < size; ++i)
  {
    host_ptr[i] = -i;
  }

  dst->updateFrom(src);
  host->updateFrom(dst);

  for(int i = 0; i < size; ++i)
  {
    ASSERT_EQ(host_ptr[i], i);
  }
}

const std::string copy_locations[] = {"NEW",
                                      "MALLOC",
                                      "STATIC"
#if defined(AXOM_USE_UMPIRE)
                                      ,
                                      "HOST"
  #if defined(UMPIRE_ENABLE_DEVICE)
                                      ,
                                      "DEVICE"
  #endif
  #if defined(UMPIRE_ENABLE_UM)
                                      ,
                                      "UM"
  #endif
  #if defined(UMPIRE_ENABLE_PINNED)
                                      ,
                                      "PINNED"
  #endif
#endif
};

const int offsets[] = {0, 29};

INSTANTIATE_TEST_SUITE_P(sidre_view,
                         UpdateTest,
                         ::testing::Combine(::testing::ValuesIn(copy_locations),
                                            ::testing::ValuesIn(copy_locations),
                                            ::testing::ValuesIn(offsets)));

//----------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();

  return result;
}
