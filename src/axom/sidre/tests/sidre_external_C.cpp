// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/sidre/interface/sidre.h"

#include <stdlib.h>

//------------------------------------------------------------------------------
// Test Group::create_external_view()
//------------------------------------------------------------------------------
TEST(C_sidre_external, create_external_view)
{
  SIDRE_DataStore ds_buf;
  SIDRE_Group root_buf;
  SIDRE_View iview_buf, dview_buf;

  SIDRE_DataStore* ds = SIDRE_DataStore_new(&ds_buf);
  SIDRE_Group* root = SIDRE_DataStore_get_root(ds, &root_buf);

  const int len = 11;

  int* idata = (int*)malloc(sizeof(int) * len);
  double* ddata = (double*)malloc(sizeof(double) * len);

  for(int ii = 0; ii < len; ++ii)
  {
    idata[ii] = ii;
    ddata[ii] = idata[ii] * 2.0;
  }

  SIDRE_View* iview = SIDRE_Group_create_view_external(root, "idata", idata, &iview_buf);
  SIDRE_View_apply_type_nelems(iview, SIDRE_INT_ID, len);
  SIDRE_View* dview = SIDRE_Group_create_view_external(root, "ddata", ddata, &dview_buf);
  SIDRE_View_apply_type_nelems(dview, SIDRE_DOUBLE_ID, len);
  EXPECT_EQ(SIDRE_Group_get_num_views(root), 2u);

  SIDRE_View_print(iview);
  SIDRE_View_print(dview);

  int* idata_chk = (int*)SIDRE_View_get_void_ptr(iview);
  for(int ii = 0; ii < len; ++ii)
  {
    EXPECT_EQ(idata_chk[ii], idata[ii]);
  }

  double* ddata_chk = (double*)SIDRE_View_get_void_ptr(dview);
  for(int ii = 0; ii < len; ++ii)
  {
    EXPECT_EQ(ddata_chk[ii], ddata[ii]);
  }

  SIDRE_DataStore_delete(ds);
  free(idata);
  free(ddata);
}

#if 0
//------------------------------------------------------------------------------
// Test Group::save(), Group::load() with external buffers
//------------------------------------------------------------------------------
TEST(C_sidre_external, save_load_external_view)
{
  SIDRE_DataStore ds_buf;
  SIDRE_Group root_buf;

  SIDRE_DataStore* ds = SIDRE_DataStore_new(&ds_buf);
  SIDRE_Group* root = SIDRE_DataStore_get_root(ds, &root_buf);

  const int len = 11;

  int* idata = (int*) malloc(sizeof(int) * len);
  double* ddata = (double*) malloc(sizeof(double) * len);

  for (int ii = 0 ; ii < len ; ++ii)
  {
    idata[ii] = ii;
    ddata[ii] = idata[ii] * 2.0;
  }

  SIDRE_View* iview =
    SIDRE_Group_create_view_external(root, "idata", idata);
  SIDRE_View_apply_type_nelems(iview, SIDRE_INT_ID, len);
  SIDRE_View* dview =
    SIDRE_Group_create_view_external(root, "ddata", ddata);
  SIDRE_View_apply_type_nelems(dview, SIDRE_DOUBLE_ID, len);

  EXPECT_EQ(SIDRE_Group_get_num_views(root), 2u);

  SIDRE_View_print(iview);
  SIDRE_View_print(dview);

  SIDRE_Group_save(root, "out_sidre_external_save_restore_external_view",
                   "conduit");

  SIDRE_DataStore_print(ds);


  SIDRE_DataStore* ds2 = SIDRE_DataStore_new();
  SIDRE_Group* root2 = SIDRE_DataStore_get_root(ds);

  SIDRE_Group_load(root, "out_sidre_external_save_restore_external_view",
                   "conduit");

  SIDRE_DataStore_print(ds2);

  SIDRE_View* iview2 = SIDRE_Group_get_view_from_name(root2, "idata");
  SIDRE_View* dview2 = SIDRE_Group_get_view_from_name(root2, "ddata");

  EXPECT_EQ(SIDRE_Group_get_num_views(root2), 2u);

  int* idata_chk = (int*) SIDRE_View_get_void_ptr(iview2);
  for (int ii = 0 ; ii < len ; ++ii)
  {
    EXPECT_EQ(idata_chk[ii], idata[ii]);
  }

  double* ddata_chk = (double*) SIDRE_View_get_void_ptr(dview2);
  for (int ii = 0 ; ii < len ; ++ii)
  {
    EXPECT_EQ(ddata_chk[ii], ddata[ii]);
  }

  SIDRE_DataStore_delete(ds);
  SIDRE_DataStore_delete(ds2);
  free(idata);
  free(ddata);
}

#endif
