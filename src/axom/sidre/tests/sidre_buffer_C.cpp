// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <stdint.h>

#include "gtest/gtest.h"

#include "axom/sidre/interface/sidre.h"

//------------------------------------------------------------------------------

TEST(C_sidre_buffer, create_buffers)
{
  SIDRE_DataStore ds_buf;
  SIDRE_Buffer dbuff_0_buf, dbuff_1_buf, dbuff_3_buf;

  SIDRE_DataStore* ds = SIDRE_DataStore_new(&ds_buf);
  SIDRE_Buffer* dbuff_0 = SIDRE_DataStore_create_buffer_empty(ds, &dbuff_0_buf);
  SIDRE_Buffer* dbuff_1 = SIDRE_DataStore_create_buffer_empty(ds, &dbuff_1_buf);

  EXPECT_EQ(SIDRE_Buffer_get_index(dbuff_0), 0);
  EXPECT_EQ(SIDRE_Buffer_get_index(dbuff_1), 1);
  SIDRE_DataStore_destroy_buffer(ds, 0);

  SIDRE_Buffer* dbuff_3 = SIDRE_DataStore_create_buffer_empty(ds, &dbuff_3_buf);
  EXPECT_EQ(SIDRE_Buffer_get_index(dbuff_3), 0);

  SIDRE_DataStore_print(ds);
  SIDRE_DataStore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_buffer, alloc_buffer_for_int_array)
{
  SIDRE_DataStore ds_buf;
  SIDRE_Buffer dbuff_buf;

  SIDRE_DataStore* ds = SIDRE_DataStore_new(&ds_buf);
  SIDRE_Buffer* dbuff = SIDRE_DataStore_create_buffer_from_type(ds, SIDRE_INT_ID, 10, &dbuff_buf);

  //  SIDRE_Buffer_declare(dbuff, SIDRE_INT_ID, 10);
  SIDRE_Buffer_allocate_existing(dbuff);

  EXPECT_EQ(SIDRE_Buffer_get_type_id(dbuff), SIDRE_INT_ID);
  EXPECT_EQ(SIDRE_Buffer_get_num_elements(dbuff), 10u);
  EXPECT_EQ(SIDRE_Buffer_get_bytes_per_element(dbuff), sizeof(int));
  EXPECT_EQ(SIDRE_Buffer_get_total_bytes(dbuff), sizeof(int) * 10);

  int* data_ptr = (int*)SIDRE_Buffer_get_void_ptr(dbuff);

  for(int i = 0; i < 10; i++)
  {
    data_ptr[i] = i * i;
  }

  SIDRE_Buffer_print(dbuff);

  SIDRE_DataStore_print(ds);
  SIDRE_DataStore_delete(ds);
}

//------------------------------------------------------------------------------

TEST(C_sidre_buffer, init_buffer_for_int_array)
{
  SIDRE_DataStore ds_buf;
  SIDRE_Buffer dbuff_buf;

  SIDRE_DataStore* ds = SIDRE_DataStore_new(&ds_buf);
  SIDRE_Buffer* dbuff = SIDRE_DataStore_create_buffer_empty(ds, &dbuff_buf);

  SIDRE_Buffer_allocate_from_type(dbuff, SIDRE_INT_ID, 10);

  EXPECT_EQ(SIDRE_Buffer_get_type_id(dbuff), SIDRE_INT_ID);
  EXPECT_EQ(SIDRE_Buffer_get_num_elements(dbuff), 10u);
  EXPECT_EQ(SIDRE_Buffer_get_bytes_per_element(dbuff), sizeof(int));
  EXPECT_EQ(SIDRE_Buffer_get_total_bytes(dbuff), sizeof(int) * 10);

  int* data_ptr = (int*)SIDRE_Buffer_get_void_ptr(dbuff);

  for(int i = 0; i < 10; i++)
  {
    data_ptr[i] = i * i;
  }

  SIDRE_Buffer_print(dbuff);

  SIDRE_DataStore_print(ds);
  SIDRE_DataStore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_buffer, realloc_buffer)
{
  SIDRE_DataStore ds_buf;
  SIDRE_Buffer dbuff_buf;

  SIDRE_DataStore* ds = SIDRE_DataStore_new(&ds_buf);
  SIDRE_Buffer* dbuff = SIDRE_DataStore_create_buffer_from_type(ds, SIDRE_LONG_ID, 5, &dbuff_buf);

  SIDRE_Buffer_allocate_existing(dbuff);

  EXPECT_EQ(SIDRE_Buffer_get_type_id(dbuff), SIDRE_LONG_ID);
  EXPECT_EQ(SIDRE_Buffer_get_num_elements(dbuff), 5u);
  EXPECT_EQ(SIDRE_Buffer_get_bytes_per_element(dbuff), sizeof(long));
  EXPECT_EQ(SIDRE_Buffer_get_total_bytes(dbuff), sizeof(long) * 5);

  long* data_ptr = (long*)SIDRE_Buffer_get_void_ptr(dbuff);

  for(int i = 0; i < 5; i++)
  {
    data_ptr[i] = 5;
  }

  SIDRE_Buffer_print(dbuff);

  SIDRE_Buffer_reallocate(dbuff, 10);

  EXPECT_EQ(SIDRE_Buffer_get_type_id(dbuff), SIDRE_LONG_ID);
  EXPECT_EQ(SIDRE_Buffer_get_num_elements(dbuff), 10u);
  EXPECT_EQ(SIDRE_Buffer_get_bytes_per_element(dbuff), sizeof(long));
  EXPECT_EQ(SIDRE_Buffer_get_total_bytes(dbuff), sizeof(long) * 10);

  // data buffer changes
  data_ptr = (long*)SIDRE_Buffer_get_void_ptr(dbuff);

  for(int i = 0; i < 5; i++)
  {
    EXPECT_EQ(data_ptr[i], 5);
  }

  for(int i = 5; i < 10; i++)
  {
    data_ptr[i] = 10;
  }

  SIDRE_Buffer_print(dbuff);

  SIDRE_DataStore_print(ds);
  SIDRE_DataStore_delete(ds);
}
