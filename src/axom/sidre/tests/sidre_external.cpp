// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core/Types.hpp"
#include "axom/sidre.hpp"

#include "gtest/gtest.h"

#include <vector>

using axom::sidre::DataStore;
using axom::sidre::DOUBLE_ID;
using axom::sidre::Group;
using axom::sidre::IndexType;
using axom::sidre::INT_ID;
using axom::sidre::View;

/* This test code contains snippets used in the Sidre Sphinx documentation.
 * They begin and end with comments.
 *
 * external_save_load_start
 * external_save_load_end
 */

//------------------------------------------------------------------------------
// Test Group::createView() -- external
//------------------------------------------------------------------------------
TEST(sidre_external, create_external_view)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  const IndexType len = 11;
  const int ndims = 1;
  IndexType shape[] = {len};

  int* idata = new int[len];

  for(int ii = 0; ii < len; ++ii)
  {
    idata[ii] = ii;
  }

  for(unsigned int i = 0; i < 8; i++)
  {
    View* view = nullptr;

    switch(i)
    {
    case 0:
      view = root->createView("data0", INT_ID, len, idata);
      break;
    case 1:
      view = root->createView("data1", INT_ID, len)->setExternalDataPtr(idata);
      break;
    case 2:
      view = root->createView("data2")->setExternalDataPtr(INT_ID, len, idata);
      break;
    case 3:
      view = root->createView("data3", idata)->apply(INT_ID, len);
      break;

    case 4:
      view = root->createViewWithShape("data4", INT_ID, ndims, shape, idata);
      break;
    case 5:
      view = root->createViewWithShape("data5", INT_ID, ndims, shape)->setExternalDataPtr(idata);
      break;
    case 6:
      view = root->createView("data6")->setExternalDataPtr(INT_ID, ndims, shape, idata);
      break;
    case 7:
      view = root->createView("data7", idata)->apply(INT_ID, ndims, shape);
      break;
    }

    EXPECT_TRUE(view != nullptr);
    EXPECT_EQ(root->getNumViews(), i + 1);

    EXPECT_TRUE(view->isDescribed());
    EXPECT_TRUE(view->isAllocated());
    EXPECT_TRUE(view->isApplied());

    EXPECT_TRUE(view->isExternal());
    EXPECT_FALSE(view->isOpaque());

    EXPECT_EQ(view->getTypeID(), INT_ID);
    EXPECT_EQ(view->getNumElements(), len);

    EXPECT_EQ(idata, view->getVoidPtr());

    view->print();

    int* idata_chk = view->getData();
    for(int ii = 0; ii < len; ++ii)
    {
      EXPECT_EQ(idata_chk[ii], idata[ii]);
    }
  }
  delete ds;
  delete[] idata;
}

//------------------------------------------------------------------------------
// Test Group::createView() -- external, null pointer
//------------------------------------------------------------------------------
TEST(sidre_external, create_external_view_null)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  int* idata = nullptr;

  for(int i = 0; i < 2; i++)
  {
    View* view = nullptr;

    switch(i)
    {
    case 0:
      view = root->createView("null0", idata)->apply(INT_ID, 0);
      break;
    case 1:
      view = root->createView("null1", INT_ID, 0)->setExternalDataPtr(idata);
      break;
    case 2:
      //    view = root->createView("null2", INT_ID, 0, idata);
      break;
    }

    EXPECT_TRUE(view != nullptr);
    EXPECT_TRUE(view->isDescribed());
    EXPECT_FALSE(view->isAllocated());
    EXPECT_FALSE(view->isApplied());

    EXPECT_TRUE(view->isEmpty());
    EXPECT_FALSE(view->isExternal());
    EXPECT_FALSE(view->isOpaque());

    EXPECT_EQ(view->getTypeID(), INT_ID);
    EXPECT_EQ(view->getNumElements(), 0u);
    EXPECT_EQ(view->getTotalBytes(), 0u);

    void* ptr = view->getVoidPtr();
    EXPECT_EQ(static_cast<void*>(nullptr), ptr);

    // getData will not work since the address is NULL
    //  int * idata_chk = view->getData();
    //EXPECT_EQ(idata_chk, nullptr);

    view->print();
  }
  delete ds;
}

//------------------------------------------------------------------------------
// Test EXTERNAL-> EMPTY
//------------------------------------------------------------------------------
TEST(sidre_external, transition_external_view_to_empty)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  const IndexType len = 11;
  int idata[len];

  for(int ii = 0; ii < len; ++ii)
  {
    idata[ii] = ii;
  }

  View* view = root->createView("data0", INT_ID, len)->setExternalDataPtr(idata);
  EXPECT_TRUE(view != nullptr);
  EXPECT_TRUE(view->isExternal());
  EXPECT_EQ(idata, view->getVoidPtr());

  // Transition from EXTERNAL to EMPTY
  view->setExternalDataPtr(nullptr);

  EXPECT_TRUE(view->isDescribed());
  EXPECT_FALSE(view->isAllocated());
  EXPECT_FALSE(view->isApplied());

  EXPECT_TRUE(view->isEmpty());
  EXPECT_FALSE(view->isExternal());
  EXPECT_FALSE(view->isOpaque());

  EXPECT_EQ(view->getTypeID(), INT_ID);
  EXPECT_EQ(view->getNumElements(), len);

  void* ptr = view->getVoidPtr();
  EXPECT_EQ(static_cast<void*>(nullptr), ptr);

  view->print();

  delete ds;
}

TEST(sidre_external, verify_external_layout)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  const int SZ = 11;
  int extData[SZ];

  SLIC_INFO("Tests that the hierarchy generated by Group::"
            << "createExternalLayout() consists only of groups and views"
            << " associated with external pointers (described or undescribed).");

  // Create some internal views
  root->createViewAndAllocate("int/desc/bufferview", axom::sidre::INT_ID, SZ);
  root->createViewScalar("int/scalar/scalarview", SZ);
  root->createViewString("int/string/stringview", "A string");

  // Sanity check on internal views
  EXPECT_FALSE(root->getView("int/desc/bufferview")->isExternal());
  EXPECT_TRUE(root->getView("int/desc/bufferview")->isDescribed());
  EXPECT_FALSE(root->getView("int/scalar/scalarview")->isExternal());
  EXPECT_TRUE(root->getView("int/scalar/scalarview")->isDescribed());
  EXPECT_FALSE(root->getView("int/string/stringview")->isExternal());
  EXPECT_TRUE(root->getView("int/string/stringview")->isDescribed());

  // Check that the generated external layout is empty
  {
    axom::sidre::Node emptyNode;
    root->createExternalLayout(emptyNode);

    SLIC_INFO("External node's json before adding external data." << " (should be empty): \n\t"
                                                                  << emptyNode.to_json());

    EXPECT_EQ(0, emptyNode.number_of_children());
  }

  // Create some external views
  root->createView("ext/desc/external_desc", axom::sidre::INT_ID, SZ, extData);
  root->createView("ext/undesc/external_opaque")->setExternalDataPtr(extData);

  // Sanity check on the external views
  EXPECT_TRUE(root->getView("ext/desc/external_desc")->isExternal());
  EXPECT_TRUE(root->getView("ext/desc/external_desc")->isDescribed());

  EXPECT_TRUE(root->getView("ext/undesc/external_opaque")->isExternal());
  EXPECT_FALSE(root->getView("ext/undesc/external_opaque")->isDescribed());

  // Initialize the data so we can test it later
  int* bufData = root->getView("int/desc/bufferview")->getData();
  for(int i = 0; i < SZ; ++i)
  {
    extData[i] = i;
    bufData[i] = 100 + i;
  }

  // Check that the generated external layout matches expectation
  {
    /// Copy the external layout into a conduit node and test the layout
    axom::sidre::Node node;
    root->createExternalLayout(node);

    SLIC_INFO("External node's json: \n\t" << node.to_json());

    EXPECT_TRUE(node.number_of_children() == 1);
    EXPECT_TRUE(node.has_path("ext"));
    EXPECT_TRUE(node["ext"].number_of_children() == 2);

    // Described external views are present and we can access the data
    EXPECT_TRUE(node.has_path("ext/desc"));
    EXPECT_TRUE(node.has_path("ext/desc/external_desc"));
    int* extLayoutData = node["ext/desc/external_desc"].value();
    for(int i = 0; i < SZ; ++i)
    {
      EXPECT_EQ(extData[i], extLayoutData[i]);
    }

    // Opaque views are not present, but their containing groups are present
    EXPECT_TRUE(node.has_path("ext/undesc"));
    EXPECT_FALSE(node.has_path("ext/undesc/external_opaque"));

    // Buffer and scalar views are not present in the external layout
    EXPECT_FALSE(node.has_path("int"));
  }

  delete ds;
}

//------------------------------------------------------------------------------
// Test Group::save(), Group::load() with described external views
//------------------------------------------------------------------------------
TEST(sidre_external, save_load_external_view)
{
#ifndef AXOM_USE_HDF5
  SUCCEED() << "sidre::Group::loadExternalData() is only implemented "
               "for the 'sidre_hdf5' protocol";
  return;
#endif

  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  constexpr IndexType len = 11;

  std::array<int, len> idata;
  std::array<double, len> ddata;

  for(int ii = 0; ii < len; ++ii)
  {
    idata[ii] = ii;
    ddata[ii] = idata[ii] * 2.0;
  }

  root->createView("idata", idata.data())->apply(INT_ID, len);
  root->createView("ddata", ddata.data())->apply(DOUBLE_ID, len);
  EXPECT_EQ(root->getNumViews(), 2u);

  ds->getRoot()->save("sidre_external_save_load_external_view", "sidre_hdf5");

  DataStore* ds2 = new DataStore();
  Group* load_group = ds2->getRoot();

  // _external_save_load_start
  // Load from file, the Views with external data will be described but
  // have no pointer to data
  load_group->load("sidre_external_save_load_external_view");

  // Verify load_group has external Views named "idata" and "ddata".
  // idata describes an int array of length len.
  // ddata describes a double array of length len.
  // These Views do not yet have valid data.
  EXPECT_TRUE(load_group->hasView("idata"));
  EXPECT_TRUE(load_group->hasView("ddata"));
  View* load_idata = load_group->getView("idata");
  View* load_ddata = load_group->getView("ddata");
  EXPECT_TRUE(load_idata->isExternal());
  EXPECT_TRUE(load_ddata->isExternal());
  EXPECT_TRUE(load_idata->getNumElements() == len);
  EXPECT_TRUE(load_ddata->getNumElements() == len);
  EXPECT_TRUE(load_idata->getTypeID() == INT_ID);
  EXPECT_TRUE(load_ddata->getTypeID() == DOUBLE_ID);

  // Create arrays that will serve as locations for external data
  std::array<int, len> new_idata;
  std::array<double, len> new_ddata;

  // Set the new arrays' pointers into the Views
  load_idata->setExternalDataPtr(new_idata.data());
  load_ddata->setExternalDataPtr(new_ddata.data());

  // Load external data; values located in the file will be loaded into
  // the storage identified by the external pointers.
  load_group->loadExternalData("sidre_external_save_load_external_view");
  // _external_save_load_end

  // The pointer retrieved from each View is the same address as the new
  // std::arrays.
  int* idata_chk = load_group->getView("idata")->getData();
  EXPECT_EQ(idata_chk, new_idata.data());

  // idata_chk has been loaded with the values from the file, which must
  // be the same as the original idata array that was saved
  for(int ii = 0; ii < len; ++ii)
  {
    EXPECT_EQ(idata_chk[ii], idata[ii]);
  }

  double* ddata_chk = load_group->getView("ddata")->getData();
  EXPECT_EQ(ddata_chk, new_ddata.data());
  for(int ii = 0; ii < len; ++ii)
  {
    EXPECT_EQ(ddata_chk[ii], ddata[ii]);
  }

  delete ds;
  delete ds2;
}
