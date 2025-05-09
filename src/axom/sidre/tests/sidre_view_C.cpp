// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/sidre/interface/sidre.h"

//------------------------------------------------------------------------------

TEST(C_sidre_view, create_views)
{
  SIDRE_DataStore ds_buf;
  SIDRE_Group root_buf;
  SIDRE_View dv_0_buf, dv_1_buf;
  SIDRE_Buffer db_0_buf, db_1_buf;

  SIDRE_DataStore* ds = SIDRE_DataStore_new(&ds_buf);
  SIDRE_Group* root = SIDRE_DataStore_get_root(ds, &root_buf);

  SIDRE_View* dv_0 =
    SIDRE_Group_create_view_and_allocate_nelems(root, "field0", SIDRE_INT_ID, 1, &dv_0_buf);
  SIDRE_View* dv_1 =
    SIDRE_Group_create_view_and_allocate_nelems(root, "field1", SIDRE_INT_ID, 1, &dv_1_buf);

  SIDRE_Buffer* db_0 = SIDRE_View_get_buffer(dv_0, &db_0_buf);
  SIDRE_Buffer* db_1 = SIDRE_View_get_buffer(dv_1, &db_1_buf);

  EXPECT_EQ(SIDRE_Buffer_get_index(db_0), 0);
  EXPECT_EQ(SIDRE_Buffer_get_index(db_1), 1);
  SIDRE_DataStore_delete(ds);
}

//------------------------------------------------------------------------------

TEST(C_sidre_view, int_buffer_from_view)
{
  SIDRE_DataStore ds_buf;
  SIDRE_Group root_buf;
  SIDRE_View dv_view;

  SIDRE_DataStore* ds = SIDRE_DataStore_new(&ds_buf);
  SIDRE_Group* root = SIDRE_DataStore_get_root(ds, &root_buf);

  SIDRE_View* dv =
    SIDRE_Group_create_view_and_allocate_nelems(root, "u0", SIDRE_INT_ID, 10, &dv_view);

  EXPECT_EQ(SIDRE_View_get_type_id(dv), SIDRE_INT_ID);
  int* data_ptr = (int*)SIDRE_View_get_void_ptr(dv);

  for(int i = 0; i < 10; i++)
  {
    data_ptr[i] = i * i;
  }

  SIDRE_View_print(dv);

  EXPECT_EQ(SIDRE_View_get_bytes_per_element(dv), sizeof(int));
  EXPECT_EQ(SIDRE_View_get_total_bytes(dv), sizeof(int) * 10);
  SIDRE_DataStore_delete(ds);
}

//------------------------------------------------------------------------------

TEST(C_sidre_view, int_buffer_from_view_conduit_value)
{
  SIDRE_DataStore ds_buf;
  SIDRE_Group root_buf;
  SIDRE_View dv_buf;

  SIDRE_DataStore* ds = SIDRE_DataStore_new(&ds_buf);
  SIDRE_Group* root = SIDRE_DataStore_get_root(ds, &root_buf);

  SIDRE_View* dv = SIDRE_Group_create_view_and_allocate_nelems(root, "u0", SIDRE_INT_ID, 10, &dv_buf);
  int* data_ptr = (int*)SIDRE_View_get_void_ptr(dv);

  for(int i = 0; i < 10; i++)
  {
    data_ptr[i] = i * i;
  }

  SIDRE_View_print(dv);

  EXPECT_EQ(SIDRE_View_get_bytes_per_element(dv), sizeof(int));
  EXPECT_EQ(SIDRE_View_get_total_bytes(dv), sizeof(int) * 10);
  SIDRE_DataStore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_view, int_array_strided_views)
{
  const unsigned int num_elts = 10;
  const unsigned int num_view_elts = 5;
  const unsigned int offset_even = 0;
  const unsigned int offset_odd = 1;
  const unsigned int stride = 2;
  const unsigned int elt_bytes = sizeof(int);

  SIDRE_DataStore ds_buf;
  SIDRE_Group root_buf;
  SIDRE_Buffer dbuff_buf;
  SIDRE_View dv_e_buf, dv_o_buf;
  SIDRE_View dv_e1_buf, dv_o1_buf;

  SIDRE_DataStore* ds = SIDRE_DataStore_new(&ds_buf);
  SIDRE_Group* root = SIDRE_DataStore_get_root(ds, &root_buf);
  SIDRE_Buffer* dbuff =
    SIDRE_DataStore_create_buffer_from_type(ds, SIDRE_INT_ID, num_elts, &dbuff_buf);

  SIDRE_Buffer_allocate_existing(dbuff);
  int* data_ptr = (int*)SIDRE_Buffer_get_void_ptr(dbuff);

  for(int i = 0; i < (int)num_elts; i++)
  {
    data_ptr[i] = i;
  }

  SIDRE_Buffer_print(dbuff);

  EXPECT_EQ(num_elts * elt_bytes, SIDRE_Buffer_get_total_bytes(dbuff));

  SIDRE_View* dv_e = SIDRE_Group_create_view_into_buffer(root, "even", dbuff, &dv_e_buf);
  SIDRE_View* dv_o = SIDRE_Group_create_view_into_buffer(root, "odd", dbuff, &dv_o_buf);
  EXPECT_TRUE(dv_e != NULL);
  EXPECT_TRUE(dv_o != NULL);
  EXPECT_EQ(SIDRE_Buffer_get_num_views(dbuff), 2u);

  SIDRE_View_apply_nelems_offset_stride(dv_e, num_view_elts, offset_even, stride);
  SIDRE_View_apply_nelems_offset_stride(dv_o, num_view_elts, offset_odd, stride);

  // Test that the void_ptr matches the buffer's pointer
  EXPECT_EQ(data_ptr, SIDRE_View_get_void_ptr(dv_e));
  EXPECT_EQ(data_ptr, SIDRE_View_get_void_ptr(dv_o));

  // Test offsets and strides for dv_e and dv_o
  EXPECT_EQ(num_view_elts, SIDRE_View_get_num_elements(dv_e));
  EXPECT_EQ(offset_even, SIDRE_View_get_offset(dv_e));
  EXPECT_EQ(stride, SIDRE_View_get_stride(dv_e));

  EXPECT_EQ(num_view_elts, SIDRE_View_get_num_elements(dv_o));
  EXPECT_EQ(offset_odd, SIDRE_View_get_offset(dv_o));
  EXPECT_EQ(stride, SIDRE_View_get_stride(dv_o));

  SIDRE_View_print(dv_e);
  SIDRE_View_print(dv_o);

  // Note: Must manually take care of offsets (and strides) since templated
  //       function View::getData() is not exposed to C API for arrays
  int* dv_e_ptr = (int*)SIDRE_View_get_void_ptr(dv_e) + offset_even;
  int* dv_o_ptr = (int*)SIDRE_View_get_void_ptr(dv_o) + offset_odd;

  for(int i = 0; i < (int)num_elts; i += 2)
  {
    EXPECT_EQ(dv_e_ptr[i], i);
    EXPECT_EQ(dv_o_ptr[i], i + 1);
  }

  // Run similar test to above with different view apply method
  SIDRE_View* dv_e1 = SIDRE_Group_create_view_into_buffer(root, "even1", dbuff, &dv_e1_buf);
  SIDRE_View* dv_o1 = SIDRE_Group_create_view_into_buffer(root, "odd1", dbuff, &dv_o1_buf);
  EXPECT_TRUE(dv_e1 != NULL);
  EXPECT_TRUE(dv_o1 != NULL);
  EXPECT_EQ(SIDRE_Buffer_get_num_views(dbuff), 4u);

  SIDRE_View_apply_type_nelems_offset_stride(dv_e1, SIDRE_INT_ID, num_view_elts, offset_even, stride);
  SIDRE_View_apply_type_nelems_offset_stride(dv_o1, SIDRE_INT_ID, num_view_elts, offset_odd, stride);

  // Test that the void_ptr matches the buffer's pointer
  EXPECT_EQ(data_ptr, SIDRE_View_get_void_ptr(dv_e1));
  EXPECT_EQ(data_ptr, SIDRE_View_get_void_ptr(dv_o1));

  // Test offsets and strides for dv_e and dv_o
  EXPECT_EQ(num_view_elts, SIDRE_View_get_num_elements(dv_e1));
  EXPECT_EQ(offset_even, SIDRE_View_get_offset(dv_e1));
  EXPECT_EQ(stride, SIDRE_View_get_stride(dv_e1));

  EXPECT_EQ(num_view_elts, SIDRE_View_get_num_elements(dv_o1));
  EXPECT_EQ(offset_odd, SIDRE_View_get_offset(dv_o1));
  EXPECT_EQ(stride, SIDRE_View_get_stride(dv_o1));

  SIDRE_View_print(dv_e1);
  SIDRE_View_print(dv_o1);

  // Note: Must manually take care of offsets (and strides).
  int* dv_e1_ptr = (int*)SIDRE_View_get_void_ptr(dv_e1) + offset_even;
  int* dv_o1_ptr = (int*)SIDRE_View_get_void_ptr(dv_o1) + offset_odd;
  for(int i = 0; i < (int)num_elts; i += 2)
  {
    EXPECT_EQ(dv_e1_ptr[i], i);
    EXPECT_EQ(dv_o1_ptr[i], i + 1);
  }

  SIDRE_DataStore_print(ds);
  SIDRE_DataStore_delete(ds);
}
//------------------------------------------------------------------------------
TEST(C_sidre_view, int_array_depth_view)
{
  SIDRE_DataStore ds_buf;
  SIDRE_Group root_buf;
  SIDRE_Buffer dbuff_buf;

  SIDRE_DataStore* ds = SIDRE_DataStore_new(&ds_buf);
  SIDRE_Group* root = SIDRE_DataStore_get_root(ds, &root_buf);

  const size_t depth_nelems = 10;
  SIDRE_Buffer* dbuff =
    SIDRE_DataStore_create_buffer_from_type(ds, SIDRE_INT_ID, 4 * depth_nelems, &dbuff_buf);

  SIDRE_Buffer_allocate_existing(dbuff);
  int* data_ptr = (int*)SIDRE_Buffer_get_void_ptr(dbuff);

  for(size_t i = 0; i < 4 * depth_nelems; i++)
  {
    data_ptr[i] = i / depth_nelems;
  }

  SIDRE_Buffer_print(dbuff);

  EXPECT_EQ(SIDRE_Buffer_get_num_elements(dbuff), 4 * depth_nelems);

  // create 4 "depth" views and apply offsets into buffer
  SIDRE_View views_buf[4];
  SIDRE_View* views[4];
  const char* view_names[4] = {"depth_0", "depth_1", "depth_2", "depth_3"};
  unsigned int view_offsets[4];

  for(int id = 0; id < 4; ++id)
  {
    views[id] = SIDRE_Group_create_view_into_buffer(root, view_names[id], dbuff, views_buf + id);
    view_offsets[id] = id * depth_nelems;

    SIDRE_View_apply_nelems_offset(views[id], depth_nelems, view_offsets[id]);
  }
  EXPECT_EQ(SIDRE_Buffer_get_num_views(dbuff), 4u);

  // print depth views...
  for(int id = 0; id < 4; ++id)
  {
    SIDRE_View_print(views[id]);
  }

  // check values in depth views...
  for(int id = 0; id < 4; ++id)
  {
    EXPECT_EQ(data_ptr, SIDRE_View_get_void_ptr(views[id]));

    // Note: offsets and striding must be handled manually when
    //       accessing the data with the get_void_ptr() function.
    //       This is the only array access function exposed by the Sidre C API
    unsigned int offset = SIDRE_View_get_offset(views[id]);
    EXPECT_EQ(view_offsets[id], offset);
    int* dv_ptr = (int*)SIDRE_View_get_void_ptr(views[id]) + offset;
    for(size_t i = 0; i < depth_nelems; ++i)
    {
      EXPECT_EQ(dv_ptr[i], id);
    }
  }

  SIDRE_DataStore_print(ds);
  SIDRE_DataStore_delete(ds);
}
//------------------------------------------------------------------------------
TEST(sidre_view, int_array_view_attach_buffer)
{
  SIDRE_DataStore ds_buf;
  SIDRE_Group root_buf;
  SIDRE_View field0_buf, field1_buf;
  SIDRE_Buffer dbuff_buf;

  SIDRE_DataStore* ds = SIDRE_DataStore_new(&ds_buf);
  SIDRE_Group* root = SIDRE_DataStore_get_root(ds, &root_buf);

  const size_t field_nelems = 10;

  // create 2 "field" views with type and # elems
  size_t elem_count = 0;
  SIDRE_View* field0 =
    SIDRE_Group_create_view_from_type(root, "field0", SIDRE_INT_ID, field_nelems, &field0_buf);
  elem_count += SIDRE_View_get_num_elements(field0);
  SIDRE_View* field1 =
    SIDRE_Group_create_view_from_type(root, "field1", SIDRE_INT_ID, field_nelems, &field1_buf);
  elem_count += SIDRE_View_get_num_elements(field1);
  EXPECT_EQ(elem_count, 2 * field_nelems);

  // create buffer to hold data for all fields and allocate
  SIDRE_Buffer* dbuff =
    SIDRE_DataStore_create_buffer_from_type(ds, SIDRE_INT_ID, 2 * field_nelems, &dbuff_buf);
  SIDRE_Buffer_allocate_existing(dbuff);
  EXPECT_EQ(SIDRE_Buffer_get_num_elements(dbuff), elem_count);

  // Initilize buffer data for testing below.
  int* b_ptr = (int*)SIDRE_Buffer_get_void_ptr(dbuff);
  for(size_t i = 0; i < elem_count; ++i)
  {
    b_ptr[i] = i / field_nelems;
  }

  SIDRE_Buffer_print(dbuff);

  const int offset0 = 0 * field_nelems;
  const int offset1 = 1 * field_nelems;

  // attach field views to buffer and apply offsets into buffer
  SIDRE_View_attach_buffer_only(field0, dbuff);
  SIDRE_View_apply_nelems_offset(field0, field_nelems, offset0);
  SIDRE_View_attach_buffer_only(field1, dbuff);
  SIDRE_View_apply_nelems_offset(field1, field_nelems, offset1);

  EXPECT_EQ(SIDRE_Buffer_get_num_views(dbuff), 2u);
  EXPECT_EQ(b_ptr, SIDRE_View_get_void_ptr(field0));
  EXPECT_EQ(b_ptr, SIDRE_View_get_void_ptr(field1));

  // print field views...
  SIDRE_View_print(field0);
  SIDRE_View_print(field1);

  // check values in field views...
  int* f0_ptr = (int*)SIDRE_View_get_void_ptr(field0) + offset0;
  for(size_t i = 0; i < field_nelems; ++i)
  {
    EXPECT_EQ(f0_ptr[i], 0);
  }
  int* f1_ptr = (int*)SIDRE_View_get_void_ptr(field1) + offset1;
  for(size_t i = 0; i < field_nelems; ++i)
  {
    EXPECT_EQ(f1_ptr[i], 1);
  }

  SIDRE_DataStore_print(ds);
  SIDRE_DataStore_delete(ds);
}
//------------------------------------------------------------------------------

TEST(C_sidre_view, int_array_multi_view_resize)
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
  SIDRE_DataStore ds_buf;
  SIDRE_Group root_buf;
  SIDRE_Group r_old_buf;
  SIDRE_View base_old_buf;

  // create our main data store
  SIDRE_DataStore* ds = SIDRE_DataStore_new(&ds_buf);
  // get access to our root data Group
  SIDRE_Group* root = SIDRE_DataStore_get_root(ds, &root_buf);

  // create a group to hold the "old" or data we want to copy
  SIDRE_Group* r_old = SIDRE_Group_create_group(root, "r_old", &r_old_buf);
  // create a view to hold the base buffer
  SIDRE_View* base_old =
    SIDRE_Group_create_view_and_allocate_nelems(r_old, "base_data", SIDRE_INT_ID, 40, &base_old_buf);

  int* data_ptr = (int*)SIDRE_View_get_void_ptr(base_old);

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

#ifdef XXX
  /// setup our 4 views
  SIDRE_Buffer* buff_old = SIDRE_View_get_buffer(base_old);
  buff_old->getNode().print();
  SIDRE_View* r0_old = SIDRE_View_create_view(r_old, "r0", buff_old);
  SIDRE_View* r1_old = SIDRE_View_create_view(r_old, "r1", buff_old);
  SIDRE_View* r2_old = SIDRE_View_create_view(r_old, "r2", buff_old);
  SIDRE_View* r3_old = SIDRE_View_create_view(r_old, "r3", buff_old);

  // each view is offset by 10 * the # of bytes in a uint32
  // uint32(num_elems, offset)
  index_t offset = 0;
  r0_old->apply(r0_old, DataType::uint32(10, offset));

  offset += sizeof(int) * 10;
  r1_old->apply(r1_old, DataType::uint32(10, offset));

  offset += sizeof(int) * 10;
  r2_old->apply(r2_old, DataType::uint32(10, offset));

  offset += sizeof(int) * 10;
  r3_old->apply(r3_old, DataType::uint32(10, offset));

  /// check that our views actually point to the expected data
  //
  uint32* r0_ptr = r0_old->getNode().as_uint32_ptr();
  for(int i = 0; i < 10; i++)
  {
    EXPECT_EQ(r0_ptr[i], 1u);
    // check pointer relation
    EXPECT_EQ(&r0_ptr[i], &data_ptr[i]);
  }

  uint32* r3_ptr = r3_old->getNode().as_uint32_ptr();
  for(int i = 0; i < 10; i++)
  {
    EXPECT_EQ(r3_ptr[i], 4u);
    // check pointer relation
    EXPECT_EQ(&r3_ptr[i], &data_ptr[i + 30]);
  }

  // create a group to hold the "old" or data we want to copy into
  SIDRE_Group* r_new = SIDRE_Group_create_group(root, "r_new");
  // create a view to hold the base buffer
  SIDRE_View* base_new = SIDRE_Group_create_view_and_buffer_simple(r_new, "base_data");

  // alloc our buffer
  // create a buffer to hold larger subarrays
  base_new->allocate_from_type(base_new, DataType::uint32(4 * 12));
  int* base_new_data = (int*)SIDRE_Buffer_det_data(base_new);
  for(int i = 0; i < 4 * 12; ++i)
  {
    base_new_data[i] = 0;
  }

  SIDRE_Buffer* buff_new = SIDRE_View_get_buffer(base_new);
  buff_new->getNode().print();

  // create the 4 sub views of this array
  SIDRE_View* r0_new = SIDRE_Group_create_view(r_new, "r0", buff_new);
  SIDRE_View* r1_new = SIDRE_Group_create_view(r_new, "r1", buff_new);
  SIDRE_View* r2_new = SIDRE_Group_create_view(r_new, "r2", buff_new);
  SIDRE_View* r3_new = SIDRE_Group_create_view(r_new, "r3", buff_new);

  // apply views to r0,r1,r2,r3
  // each view is offset by 12 * the # of bytes in a uint32

  // uint32(num_elems, offset)
  offset = 0;
  r0_new->apply(DataType::uint32(12, offset));

  offset += sizeof(int) * 12;
  r1_new->apply(DataType::uint32(12, offset));

  offset += sizeof(int) * 12;
  r2_new->apply(DataType::uint32(12, offset));

  offset += sizeof(int) * 12;
  r3_new->apply(DataType::uint32(12, offset));

  /// update r2 as an example first
  buff_new->getNode().print();
  r2_new->getNode().print();

  /// copy the subset of value
  r2_new->getNode().update(r2_old->getNode());
  r2_new->getNode().print();
  buff_new->getNode().print();

  /// check pointer values
  int* r2_new_ptr = (int*)SIDRE_View_get_void_ptr(r2_new);

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

  buff_new->getNode().print();
#endif

  SIDRE_DataStore_print(ds);
  SIDRE_DataStore_delete(ds);
}

//------------------------------------------------------------------------------

TEST(C_sidre_view, int_array_realloc)
{
  ///
  /// info
  ///
  SIDRE_DataStore ds_buf;
  SIDRE_Group root_buf;
  SIDRE_View a1_buf, a2_buf;

  // create our main data store
  SIDRE_DataStore* ds = SIDRE_DataStore_new(&ds_buf);
  // get access to our root data Group
  SIDRE_Group* root = SIDRE_DataStore_get_root(ds, &root_buf);

  // create a view to hold the base buffer
  SIDRE_View* a1 =
    SIDRE_Group_create_view_and_allocate_nelems(root, "a1", SIDRE_FLOAT_ID, 5, &a1_buf);
  SIDRE_View* a2 = SIDRE_Group_create_view_and_allocate_nelems(root, "a2", SIDRE_INT_ID, 5, &a2_buf);

  float* a1_ptr = (float*)SIDRE_View_get_void_ptr(a1);
  int* a2_ptr = (int*)SIDRE_View_get_void_ptr(a2);
  for(int i = 0; i < 5; i++)
  {
    a1_ptr[i] = 5.0;
    a2_ptr[i] = -5;
  }

  EXPECT_EQ(SIDRE_View_get_bytes_per_element(a1), sizeof(float));
  EXPECT_EQ(SIDRE_View_get_total_bytes(a1), sizeof(float) * 5);

  EXPECT_EQ(SIDRE_View_get_bytes_per_element(a2), sizeof(int));
  EXPECT_EQ(SIDRE_View_get_total_bytes(a2), sizeof(int) * 5);

  SIDRE_View_reallocate(a1, 10);
  SIDRE_View_reallocate(a2, 15);

  a1_ptr = (float*)SIDRE_View_get_void_ptr(a1);
  a2_ptr = (int*)SIDRE_View_get_void_ptr(a2);

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

  EXPECT_EQ(SIDRE_View_get_bytes_per_element(a1), sizeof(float));
  EXPECT_EQ(SIDRE_View_get_total_bytes(a1), sizeof(float) * 10);

  EXPECT_EQ(SIDRE_View_get_bytes_per_element(a2), sizeof(int));
  EXPECT_EQ(SIDRE_View_get_total_bytes(a2), sizeof(int) * 15);

  SIDRE_DataStore_print(ds);
  SIDRE_DataStore_delete(ds);
}

//------------------------------------------------------------------------------

TEST(C_sidre_view, simple_opaque)
{
  SIDRE_DataStore ds_buf;
  SIDRE_Group root_buf;
  SIDRE_View opq_view_buf;

  // create our main data store
  SIDRE_DataStore* ds = SIDRE_DataStore_new(&ds_buf);
  // get access to our root data Group
  SIDRE_Group* root = SIDRE_DataStore_get_root(ds, &root_buf);
  int* src_data = (int*)malloc(sizeof(int));

  src_data[0] = 42;

  void* src_ptr = (void*)src_data;

  SIDRE_View* opq_view = SIDRE_Group_create_view_external(root, "my_opaque", src_ptr, &opq_view_buf);

  // External pointers are held in the view, no buffer should be present.
  EXPECT_EQ(SIDRE_DataStore_get_num_buffers(ds), 0u);

  EXPECT_TRUE(SIDRE_View_is_external(opq_view));
  EXPECT_TRUE(!SIDRE_View_is_applied(opq_view));
  EXPECT_TRUE(SIDRE_View_is_opaque(opq_view));

  // Opaque views do not know their sizes
  EXPECT_EQ(0u, SIDRE_View_get_bytes_per_element(opq_view));
  EXPECT_EQ(0u, SIDRE_View_get_total_bytes(opq_view));
  EXPECT_EQ(0u, SIDRE_View_get_offset(opq_view));
  EXPECT_EQ(1u, SIDRE_View_get_stride(opq_view));

  void* opq_ptr = SIDRE_View_get_void_ptr(opq_view);

  int* out_data = (int*)opq_ptr;
  EXPECT_EQ(src_ptr, opq_ptr);
  EXPECT_EQ(42, out_data[0]);

  SIDRE_DataStore_print(ds);
  SIDRE_DataStore_delete(ds);
  free(src_data);
}
