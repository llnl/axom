! Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
! other Axom Project Developers. See the top-level LICENSE file for details.
!
! SPDX-License-Identifier: (BSD-3-Clause)

module sidre_view_test
  use iso_c_binding
  use fruit
  use axom_sidre
  implicit none

  integer, parameter :: &
       EMPTYVIEW = 1, &
       BUFFERVIEW = 2, &
       EXTERNALVIEW = 3, &
       SCALARVIEW = 4, &
       STRINGVIEW = 5, &
       NOTYPE = 6

contains
!------------------------------------------------------------------------------

  function get_state(view) result(state)
    type(SidreView), intent(IN) :: view
    integer state

    if (view%is_empty()) then
       state = EMPTYVIEW
    else if (view%has_buffer()) then
       state = BUFFERVIEW
    else if (view%is_external()) then
       state = EXTERNALVIEW
    else if (view%is_scalar()) then
       state = SCALARVIEW
    else if (view%is_string()) then
       state = STRINGVIEW
    else
       state = NOTYPE
    endif
  end function get_state

  subroutine check_view_values(view, state, is_described, is_allocated, &
       is_applied, length)

    type(SidreView), intent(IN) :: view
    integer, intent(IN) :: state
    logical, intent(IN) :: is_described, is_allocated, is_applied
    integer, intent(IN) :: length
    character(30) name

    integer(SIDRE_IndexType) dims(2)

    name = view%get_name()

    call assert_equals(get_state(view), state)
    call assert_equals(view%is_described(), is_described, trim(name) // " is_described")
    call assert_equals(view%is_allocated(), is_allocated, trim(name) // "is_allocated")
    call assert_equals(view%is_applied(), is_applied, trim(name) // " is_applied")

    call assert_equals(int(view%get_num_elements(), kind(length)), length, trim(name) // " get_num_elements")
    if (view%is_described()) then
       call assert_equals(view%get_num_dimensions(), 1, trim(name) // " get_num_dimensions")
       call assert_true(view%get_shape(1, dims) == 1, trim(name) // " get_shape")
       call assert_true(dims(1) == length, trim(name) // " dims(1)")
    endif
  end subroutine check_view_values
  
!------------------------------------------------------------------------------

  subroutine create_views()
    type(SidreDataStore) ds
    type(SidreGroup) root
    type(SidreView) dv_0, dv_1
    type(SidreBuffer) db_0, db_1

    call set_case_name("create_views")

    ds = SidreDataStore()
    root = ds%get_root()

    dv_0 = root%create_view_and_allocate("field0", SIDRE_INT_ID, 1)
    dv_1 = root%create_view_and_allocate("field1", SIDRE_INT_ID, 1)

    db_0 = dv_0%get_buffer()
    db_1 = dv_1%get_buffer()

    call assert_true(db_0%get_index() == 0, "db_0%get_index(), 0")
    call assert_true(db_1%get_index() == 1, "db_1%get_index(), 1")
    call ds%delete()
  end subroutine create_views

!------------------------------------------------------------------------------

  subroutine get_path_name
    type(SidreDataStore) ds
    type(SidreGroup) root
    type(SidreView) v1, v2, v3

    call set_case_name("get_path_name")

    ds = SidreDataStore()
    root = ds%get_root()
    v1 = root%create_view("test/a/b/v1")
    v2 = root%create_view("test/v2")
    v3 = root%create_view("v3")

    call assert_true(v1%get_name() == "v1")
    call assert_true(v1%get_path() == "test/a/b")
    call assert_true(v1%get_path_name() == "test/a/b/v1")

    call assert_true(v2%get_name() == "v2")
    call assert_true(v2%get_path() == "test")
    call assert_true(v2%get_path_name() == "test/v2")

    call assert_true(v3%get_name() == "v3")
    call assert_true(v3%get_path() == " ")
    call assert_true(v3%get_path_name() == "v3")

    call ds%delete()
  end subroutine get_path_name

!------------------------------------------------------------------------------

  subroutine scalar_view
    type(SidreDataStore) ds
    type(SidreGroup) root
    type(SidreView) i0view, i1view, s0view, s1view
    integer(C_INT) i1, i2
    character(80) s1, s2

    call set_case_name("scalar_view")

    ds = SidreDataStore()
    root = ds%get_root()

    i1 = 1
    i0view = root%create_view("i0")
    call i0view%set_scalar(i1)
    call check_scalar_values(i0view, SCALARVIEW, .true., .true., .true., SIDRE_INT_ID, 1)
    i2 = i0view%get_data_int()
    call assert_equals(i1, i2)

    i1 = 2
    i1view = root%create_view_scalar("i1", i1)
    call check_scalar_values(i1view, SCALARVIEW, .true., .true., .true., SIDRE_INT_ID, 1)
    i2 = i1view%get_data_int()
    call assert_equals(i1, i2)

    ! TODO: passing len_trim to account for non-existent NULL
    s1 = "i am a string"
    s0view = root%create_view("s0")
    call s0view%set_string(trim(s1))
    call check_scalar_values(s0view, STRINGVIEW, .true., .true., .true., &
         SIDRE_CHAR8_STR_ID, len_trim(s1) + 1)
    call s0view%get_string(s2)
    call assert_equals(s1, s2)

    s1 = "i too am a string"
    s1view = root%create_view_string("s1", trim(s1))
    call check_scalar_values(s1view, STRINGVIEW, .true., .true., .true., &
         SIDRE_CHAR8_STR_ID, len_trim(s1) + 1)
    call s1view%get_string(s2)
    call assert_equals(s1, s2)

  ! check illegal operations
!  call i0view%apply(int_id, 1)
!  call i0view%allocate()
!  call i0view%deallocate()

!  call s0view%apply(int_id, 1)
!  call s0view%allocate()
!  call s0view%deallocate()

!type(SidreView) empty
!empty = root%create_view("empty")
!#if 0
!  try
!  {
!type(int) j
!j = empty%get_scalar()
!    !int j = empty%get_scalar()
!    call assert_equals(0, *j)
!  }
!  catch ( conduit::error e)
!  {}
!#endif
!  const char * svalue = empty%get_string()
!  call assert_equals(null, svalue)

    ! Test group access to scalars
    i1 = 100
    i2 = 0
    call root%set_scalar("i0", i1)
    call root%get_scalar("i0", i2)
    call assert_equals(i1, i2)

    s1 = "Replacement string value"
    s2 = " "
    call root%set_string("s0", s1)
    call root%get_string("s0", s2)
    call assert_equals(s1, s2)

    call ds%delete()

    contains

      subroutine check_scalar_values(view, state, is_described, is_allocated, &
           is_applied, type, length)

        type(SidreView), intent(IN) :: view
        integer, intent(IN) :: state
        logical, intent(IN) :: is_described, is_allocated, is_applied
        integer(TypeID), intent(IN) :: type
        integer, intent(IN) :: length
        character(30) name

        integer(SIDRE_IndexType) dims(2)

        name = view%get_name()

        call assert_equals(get_state(view), state)
        call assert_equals(view%is_described(), is_described, trim(name) // " is_described")
        call assert_equals(view%is_allocated(), is_allocated, trim(name) // "is_allocated")
        call assert_equals(view%is_applied(), is_applied, trim(name) // " is_applied")

        call assert_true(view%get_type_id() == type, trim(name) // " get_type_id")
        call assert_equals(int(view%get_num_elements(), kind(length)), length, trim(name) // " get_num_elements")
        call assert_equals(view%get_num_dimensions(), 1, trim(name) // " get_num_dimensions")
        call assert_true(view%get_shape(1, dims) == 1, trim(name) // " get_shape")
        call assert_true(dims(1) == length, trim(name) // " dims(1)")
      end subroutine check_scalar_values

  end subroutine scalar_view

!------------------------------------------------------------------------------

  subroutine int_buffer_from_view()
    type(SidreDataStore) ds
    type(SidreGroup) root
    type(SidreView) dv
    integer(C_INT), pointer :: data(:)
    integer(C_INT) i
    integer int_size, elem_count

    int_size = c_sizeof(i)
    elem_count = 10

    call set_case_name("int_buffer_from_view")

    ds = SidreDataStore()
    root = ds%get_root()

    dv = root%create_view_and_allocate("u0", SIDRE_INT_ID, elem_count)
    call assert_true(dv%get_type_id() == SIDRE_INT_ID, "dv%get_type_id(), SIDRE_INT_ID")
    call dv%get_data(data)

    do i = 1, elem_count
       data(i) = i * i
    enddo

    call dv%print()

    call assert_true(dv%get_num_elements() == elem_count)
    call assert_true(dv%get_bytes_per_element() == int_size)
    call assert_true(dv%get_total_bytes() == int_size * elem_count)

    call ds%delete()
  end subroutine int_buffer_from_view

!------------------------------------------------------------------------------

  subroutine int_buffer_from_view_conduit_value()
    type(SidreDataStore) ds
    type(SidreGroup) root
    type(SidreView) dv
    integer(C_INT), pointer :: data(:)
    integer i

    call set_case_name("int_buffer_from_view_conduit")

    ds = SidreDataStore()
    root = ds%get_root()

    dv = root%create_view_and_allocate("u0", SIDRE_INT_ID, 10_8)
    call dv%get_data(data)

    do i = 1, 10
       data(i) = i * i
    enddo

    call dv%print()

!--    EXPECT_EQ(SIDRE_view_get_total_bytes(dv), sizeof(int) * 10)
    call ds%delete()
  end subroutine int_buffer_from_view_conduit_value

!------------------------------------------------------------------------------

  subroutine int_array_multi_view()
    type(SidreDataStore) ds
    type(SidreGroup) root
    type(SidreBuffer) dbuff
    type(SidreView) dv_e, dv_o 
    type(C_PTR) data_ptr
    integer(C_INT), pointer :: data(:)
    integer i

    call set_case_name("int_array_multi_view")

    ds = SidreDataStore()
    root = ds%get_root()
    dbuff = ds%create_buffer(SIDRE_INT_ID, 10_8)

    call dbuff%allocate()
    data_ptr = dbuff%get_void_ptr()
    call c_f_pointer(data_ptr, data, [ 10 ])

    do i = 1, 10
       data(i) = i
    enddo

    call dbuff%print()

!--#ifdef XXX
!--    EXPECT_EQ(dbuff->getNode().schema().total_bytes(),
!--              dbuff->getSchema().total_bytes())
!--#endif

    dv_e = root%create_view("even", dbuff)
    dv_o = root%create_view("odd", dbuff)
    call assert_true(dbuff%get_num_views() == 2, "dbuff%get_num_views() == 2")

!--#ifdef XXX
!--  dv_e->apply(DataType::uint32(5,0,8))
!--
!--  dv_o->apply(DataType::uint32(5,4,8))
!--
!--  call dv_e%print()
!--  call dv_o%print()
!--
!--  uint32_array dv_e_ptr = dv_e->getNode().as_uint32_array()
!--  uint32_array dv_o_ptr = dv_o->getNode().as_uint32_array()
!--  for(int i=0  i<5  i++)
!--  {
!--    std::cout << "idx:" <<  i
!--              << " e:" << dv_e_ptr[i]
!--              << " o:" << dv_o_ptr[i]
!--              << " em:" << dv_e_ptr[i]  % 2
!--              << " om:" << dv_o_ptr[i]  % 2
!--              << std::endl
!--
!--    EXPECT_EQ(dv_e_ptr[i] % 2, 0u)
!--    EXPECT_EQ(dv_o_ptr[i] % 2, 1u)
!--  }
!--#endif
    call ds%print()
    call ds%delete()
  end subroutine int_array_multi_view

!------------------------------------------------------------------------------

  subroutine init_int_array_multi_view()
    type(SidreDataStore) ds
    type(SidreGroup) root
    type(SidreBuffer) dbuff
    type(SidreView) dv_e, dv_o 
    type(C_PTR) data_ptr
    integer, pointer :: data(:)
    integer i
    
    call set_case_name("init_int_array_multi_view")

    ds = SidreDataStore()
    root = ds%get_root()
    dbuff = ds%create_buffer()
    
    call dbuff%allocate(SIDRE_INT_ID, 10_8)
    data_ptr = dbuff%get_void_ptr()
    call c_f_pointer(data_ptr, data, [ 10 ])

    do i = 1, 10
       data(i) = i
    enddo

    call dbuff%print()

!--#ifdef XXX
!--  EXPECT_EQ(dbuff->getNode().schema().total_bytes(),
!--            dbuff->getSchema().total_bytes())
!--#endif

    dv_e = root%create_view("even", dbuff)
    dv_o = root%create_view("odd", dbuff)

!--#ifdef XXX
!--  ! uint32(num_elems, offset, stride)
!--  dv_e->apply(DataType::uint32(5,0,8))
!--
!--
!--  ! uint32(num_elems, offset, stride)
!--  dv_o->apply(DataType::uint32(5,4,8))
!--
!--
!--  call dv_e%print()
!--  call dv_o%print()
!--
!--  uint32_array dv_e_ptr = dv_e->getNode().as_uint32_array()
!--  uint32_array dv_o_ptr = dv_o->getNode().as_uint32_array()
!--  for(int i=0  i<5  i++)
!--  {
!--    std::cout << "idx:" <<  i
!--              << " e:" << dv_e_ptr[i]
!--              << " o:" << dv_o_ptr[i]
!--              << " em:" << dv_e_ptr[i]  % 2
!--              << " om:" << dv_o_ptr[i]  % 2
!--              << std::endl
!--
!--    EXPECT_EQ(dv_e_ptr[i] % 2, 0u)
!--    EXPECT_EQ(dv_o_ptr[i] % 2, 1u)
!--  }
!--#endif

    call ds%print()
    call ds%delete()
  end subroutine init_int_array_multi_view

!------------------------------------------------------------------------------

  subroutine int_array_depth_view()
    type(SidreDataStore) ds
    type(SidreGroup) root
    type(SidreBuffer) dbuff
    type(SidreView) view0
    type(SidreView) view1
    type(SidreView) view2
    type(SidreView) view3
    integer(C_INT), pointer :: data(:)
    type(C_PTR) data_ptr
    integer i
    integer(SIDRE_IndexType) depth_nelems
    integer(SIDRE_IndexType) total_nelems

    call set_case_name("int_array_depth_view")

    ! create our main data store
    ds = SidreDataStore()

    depth_nelems = 10 
    total_nelems = 4 * depth_nelems

    dbuff = ds%create_buffer(SIDRE_INT_ID, total_nelems)

    ! get access to our root data Group
    root = ds%get_root()

    ! Allocate buffer to hold data for 4 "depth" views
    call dbuff%allocate()

    data_ptr = dbuff%get_void_ptr()
    call c_f_pointer(data_ptr, data, [ total_nelems ])

    do i = 1, total_nelems
       data(i) = (i - 1) / depth_nelems
    enddo

    call dbuff%print()

    call assert_true( dbuff%get_num_elements() == 4 * depth_nelems, "dbuff%get_num_elements() == 4 * depth_nelems" )

    ! create 4 "depth" views and apply offsets into buffer
    view0 = root%create_view_into_buffer("depth_view_0", dbuff)
    call view0%apply_nelems_offset(depth_nelems, 0 * depth_nelems)

    view1 = root%create_view_into_buffer("depth_view_1", dbuff)
    call view1%apply_nelems_offset(depth_nelems, 1 * depth_nelems)

    view2 = root%create_view_into_buffer("depth_view_2", dbuff)
    call view2%apply_nelems_offset(depth_nelems, 2 * depth_nelems)

    view3 = root%create_view_into_buffer("depth_view_3", dbuff)
    call view3%apply_nelems_offset(depth_nelems, 3 * depth_nelems)

    call assert_true( dbuff%get_num_views() == 4, "dbuff%get_num_views() == 4" )

    call view0%print()
    call view1%print()
    call view2%print()
    call view3%print()

    ! check values in depth views...
    call view0%get_data(data)
    call assert_true(all(data == 0), "depth 0 does not compare")

    call view1%get_data(data)
    call assert_true(all(data == 1), "depth 1 does not compare")
 
    call view2%get_data(data)
    call assert_true(all(data == 2), "depth 2 does not compare")

    call view3%get_data(data)
    call assert_true(all(data == 3), "depth 3 does not compare")

    call ds%print()
    call ds%delete()

  end subroutine int_array_depth_view

!------------------------------------------------------------------------------

  subroutine int_array_view_attach_buffer()
    type(SidreDataStore) ds
    type(SidreGroup) root
    type(SidreBuffer) dbuff
    type(SidreView) field0
    type(SidreView) field1
    integer(C_INT), pointer :: data(:)
    type(C_PTR) data_ptr
    integer i
    integer(SIDRE_IndexType) field_nelems
    integer(SIDRE_IndexType) elem_count
    integer(SIDRE_IndexType) offset0, offset1

    call set_case_name("int_array_view_attach_buffer")

    ! create our main data store
    ds = SidreDataStore()

    ! get access to our root data Group
    root = ds%get_root()

    field_nelems = 10 

    ! create 2 "field" views with type and # elems
    elem_count = 0
    field0 = root%create_view("field0", SIDRE_INT_ID, field_nelems)
    elem_count = elem_count + field0%get_num_elements()
    print*,"elem_count field0",elem_count
    field1 = root%create_view("field1", SIDRE_INT_ID, field_nelems)
    elem_count = elem_count + field1%get_num_elements()
    print*,"elem_count field1",elem_count

    call assert_true( elem_count == 2 * field_nelems, "elem_count == 2 * field_nelems" )

    ! Create buffer to hold data for all fields and allocate
    dbuff = ds%create_buffer(SIDRE_INT_ID, elem_count)
    call dbuff%allocate()

    call assert_true( dbuff%get_num_elements() == elem_count, "dbuff%get_num_elements() == elem_count" ) 

    ! Initilize buffer data for testing below
    data_ptr = dbuff%get_void_ptr()
    call c_f_pointer(data_ptr, data, [ elem_count ])

    do i = 1, elem_count
       data(i) = (i - 1) / field_nelems
    enddo

    call dbuff%print()

    ! attach field views to buffer and apply offsets into buffer
    offset0 = 0 * field_nelems
    call field0%attach_buffer(dbuff)
    call field0%apply_nelems_offset(field_nelems, offset0)

    offset1 = 1 * field_nelems
    call field1%attach_buffer(dbuff)
    call field1%apply_nelems_offset(field_nelems, offset1)

    call assert_true( dbuff%get_num_views() == 2, "dbuff%get_num_views() == 2" )

    ! print field views...
    call field0%print()
    call field1%print()

    ! check values in field views...
    call field0%get_data(data)
    call assert_true(size(data) == field_nelems, &
         "depth 0 is incorrect size")
    call assert_true(all(data == 0), "depth 0 does not compare")
    call assert_true(offset0 == field0%get_offset())

    call field1%get_data(data)
    call assert_true(size(data) == field_nelems, &
         "depth 1 is incorrect size")
    call assert_true(all(data == 1), "depth 0 does not compare")
    call assert_true(offset1 == field1%get_offset())

    call ds%print()
    call ds%delete()

  end subroutine int_array_view_attach_buffer

!------------------------------------------------------------------------------

  subroutine int_array_offset_stride()
    type(SidreDataStore) ds
    type(SidreGroup) root, other
    type(SidreBuffer) dbuff
    type(SidreView) field0, view1, view2, view3

    real(C_DOUBLE), pointer :: data(:)
    type(C_PTR) data_ptr
    integer(C_INT) i
    integer(SIDRE_IndexType) field_nelems
    integer(SIDRE_IndexType) v1_nelems, v1_stride, v1_offset
    integer(SIDRE_IndexType) v2_nelems, v2_stride, v2_offset
    integer(SIDRE_INdexType) v3_nelems, v3_stride, v3_offset
    real(C_DOUBLE), pointer :: data1(:), data2(:), data3(:)
    real(C_DOUBLE) :: elem
    integer int_size, double_size
    integer(C_INT), pointer :: int_data

    int_size = c_sizeof(i)
    double_size = c_sizeof(elem)

    call set_case_name("int_array_offset_stride")

    ! create our main data store
    ds = SidreDataStore()

    ! get access to our root data Group
    root = ds%get_root()

    field_nelems = 20
    field0 = root%create_view_and_allocate("field0", SIDRE_DOUBLE_ID, field_nelems)
    call assert_true(field0%get_num_elements() == field_nelems)
    call assert_true(field0%get_bytes_per_element() == double_size)
    call assert_true(field0%get_total_bytes() == double_size * field_nelems)
    call assert_true(field0%get_offset() == 0)
    call assert_true(field0%get_stride() == 1)

    dbuff = field0%get_buffer()

    ! Initilize buffer data for testing below
    data_ptr = dbuff%get_void_ptr()
    call c_f_pointer(data_ptr, data, [ field_nelems ])

    do i = 1, field_nelems
       data(i) = 1.001 * i
    enddo

    call dbuff%print()

    ! create two more views into field0's buffer and test stride and offset
    v1_nelems = 3
    v1_stride = 3
    v1_offset = 2
    view1 = root%create_view_into_buffer("offset_stride_1", dbuff)
    call view1%apply_nelems_offset_stride(v1_nelems, v1_offset, v1_stride)
    call view1%get_data(data1)

    v2_nelems = 3
    v2_stride = 3
    v2_offset = 3
    view2 = root%create_view_into_buffer("offset_stride_2", dbuff)
    call view2%apply_nelems_offset_stride(v2_nelems, v2_offset, v2_stride)
    call view2%get_data(data2)

    v3_nelems = 5
    v3_stride = 1
    v3_offset = 12
    view3 = root%create_view_into_buffer("offset_stride_3", dbuff)
    call view3%apply_nelems_offset_stride(v3_nelems, v3_offset, v3_stride)
    call view3%get_data(data3)

    call assert_true(view1%get_num_elements() == v1_nelems)
    call assert_true(view1%get_bytes_per_element() == double_size)
    call assert_true(view1%get_offset() == v1_offset)
    call assert_true(view1%get_stride() == v1_stride)
    call assert_true(view1%get_total_bytes() == double_size * (1 + (v1_stride * (v1_nelems-1))))
    call assert_true(data1(0) == data(v1_offset))

    call assert_true(view2%get_num_elements() == v2_nelems)
    call assert_true(view2%get_bytes_per_element() == double_size, "view 2 byes per elt")
    call assert_true(view2%get_offset() == v2_offset)
    call assert_true(view2%get_stride() == v2_stride)
    call assert_true(view2%get_total_bytes() == double_size * (1 + (v2_stride * (v2_nelems-1))))
    call assert_true(data2(0) == data(v2_offset))

    call assert_true(view3%get_num_elements() == v3_nelems)
    call assert_true(view3%get_bytes_per_element() == double_size)
    call assert_true(view3%get_offset() == v3_offset)
    call assert_true(view3%get_stride() == v3_stride)
    call assert_true(view3%get_total_bytes() == double_size * (1 + (v3_stride * (v3_nelems-1))))
    call assert_true(data3(0) == data(v3_offset))


    ! test stride and offset against other types of  views
    other = root%create_group("other")
    view1 = other%create_view("key_empty")
    call assert_true(view1%get_offset() == 0)
    call assert_true(view1%get_stride() == 1)
    call assert_true(view1%get_num_elements() == 0)
    call assert_true(view1%get_bytes_per_element() == 0)

    view1 = other%create_view("key_opaque", data_ptr) ! opaque -- not described
    call assert_true(view1%get_offset() == 0)
    call assert_true(view1%get_stride() == 1)
    call assert_true(view1%get_num_elements() == 0)
    call assert_true(view1%get_bytes_per_element() == 0)
    call assert_true(view1%get_total_bytes() == 0)

    view1 = other%create_view_string("key_str", "val_str")
    call assert_true(view1%get_offset() == 0)
    call assert_true(view1%get_stride() == 1)
    call assert_true(view1%get_bytes_per_element() == 1)

    view1 = other%create_view_scalar_int("key_int", 5)
    call view1%get_data(int_data)
    call assert_true(int_data == 5)
    call assert_true(view1%get_offset() == 0)
    call assert_true(view1%get_stride() == 1)
    call assert_true(view1%get_num_elements() == 1)
    call assert_true(view1%get_bytes_per_element() == int_size)
    call assert_true(view1%get_total_bytes() == int_size)

    ! cleanup
    call ds%print()
    call ds%delete()

  end subroutine int_array_offset_stride

!------------------------------------------------------------------------------

  subroutine int_array_multi_view_resize()
     !
     ! This example creates a 4 * 10 buffer of ints,
     ! and 4 views that point the 4 sections of 10 ints
     !
     ! We then create a new buffer to support 4*12 ints
     ! and 4 views that point into them
     !
     ! after this we use the old buffers to copy the values
     ! into the new views
     !
    type(SidreDataStore) ds
    type(SidreGroup) root, r_old
    type(SidreView) base_old
    integer(C_INT), pointer :: data(:)
    integer i

    call set_case_name("int_array_multi_view_resize")

    ! create our main data store
    ds = SidreDataStore()

    ! get access to our root data Group
    root = ds%get_root()

    ! create a group to hold the "old" or data we want to copy
    r_old = root%create_group("r_old")

    ! create a view to hold the base buffer and allocate
    ! we will create 4 sub views of this array
    base_old = r_old%create_view_and_allocate("base_data", SIDRE_INT_ID, 40)

    call base_old%get_data(data)

    ! init the buff with values that align with the
    ! 4 subsections.
    data( 1:10) = 1
    data(11:20) = 2
    data(21:30) = 3
    data(31:40) = 4

!--#ifdef XXX
!--  ! setup our 4 views
!--  SIDRE_buffer * buff_old = SIDRE_view_get_buffer(base_old)
!--  call buff_old%print()
!--  SIDRE_view * r0_old = SIDRE_view_create_view(r_old, "r0",buff_old)
!--  SIDRE_view * r1_old = SIDRE_view_create_view(r_old, "r1",buff_old)
!--  SIDRE_view * r2_old = SIDRE_view_create_view(r_old, "r2",buff_old)
!--  SIDRE_view * r3_old = SIDRE_view_create_view(r_old, "r3",buff_old)
!--
!--  ! each view is offset by 10 * the # of bytes in a uint32
!--  ! uint32(num_elems, offset)
!--  index_t offset =0
!--  r0_old->apply(r0_old, DataType::uint32(10,offset))
!--
!--  offset += sizeof(int) * 10
!--  r1_old->apply(r1_old, DataType::uint32(10,offset))
!--
!--  offset += sizeof(int) * 10
!--  r2_old->apply(r2_old, DataType::uint32(10,offset))
!--
!--  offset += sizeof(int) * 10
!--  r3_old->apply(r3_old, DataType::uint32(10,offset))
!--
!--  ! check that our views actually point to the expected data
!--  !
!--  uint32 * r0_ptr = r0_old->getNode().as_uint32_ptr()
!--  for(int i=0  i<10  i++)
!--  {
!--    EXPECT_EQ(r0_ptr[i], 1u)
!--    ! check pointer relation
!--    EXPECT_EQ(&r0_ptr[i], &data_ptr[i])
!--  }
!--
!--  uint32 * r3_ptr = r3_old->getNode().as_uint32_ptr()
!--  for(int i=0  i<10  i++)
!--  {
!--    EXPECT_EQ(r3_ptr[i], 4u)
!--    ! check pointer relation
!--    EXPECT_EQ(&r3_ptr[i], &data_ptr[i+30])
!--  }
!--
!--  ! create a group to hold the "old" or data we want to copy into
!--  SIDRE_group * r_new = SIDRE_group_create_group(root, "r_new")
!--  ! create a view to hold the base buffer
!--  SIDRE_view * base_new = SIDRE_group_create_view_and_buffer(r_new, "base_data")
!--
!--  ! alloc our buffer
!--  ! create a buffer to hold larger subarrays
!--  base_new->allocate(base_new, DataType::uint32(4 * 12))
!--  int* base_new_data = (int *) SIDRE_buffer_det_data(base_new)
!--  for (int i = 0 i < 4 * 12 ++i) 
!--  {
!--     base_new_data[i] = 0
!--  } 
!--
!--  SIDRE_buffer * buff_new = SIDRE_view_get_buffer(base_new)
!--  call buff_new%print()
!--
!--  ! create the 4 sub views of this array
!--  SIDRE_view * r0_new = SIDRE_group_create_view(r_new, "r0",buff_new)
!--  SIDRE_view * r1_new = SIDRE_group_create_view(r_new, "r1",buff_new)
!--  SIDRE_view * r2_new = SIDRE_group_create_view(r_new, "r2",buff_new)
!--  SIDRE_view * r3_new = SIDRE_group_create_view(r_new, "r3",buff_new)
!--
!--  ! apply views to r0,r1,r2,r3
!--  ! each view is offset by 12 * the # of bytes in a uint32
!--
!--  ! uint32(num_elems, offset)
!--  offset =0
!--  r0_new->apply(DataType::uint32(12,offset))
!--
!--  offset += sizeof(int) * 12
!--  r1_new->apply(DataType::uint32(12,offset))
!--
!--  offset += sizeof(int) * 12
!--  r2_new->apply(DataType::uint32(12,offset))
!--
!--  offset += sizeof(int) * 12
!--  r3_new->apply(DataType::uint32(12,offset))
!--
!--  ! update r2 as an example first
!--  call buff_new%print()
!--  call r2_new%print()
!--
!--  ! copy the subset of value
!--  r2_new->getNode().update(r2_old->getNode())
!--  call r2_new%print()
!--  call buff_new%print()
!--
!--
!--  ! check pointer values
!--  int * r2_new_ptr = (int *) SIDRE_view_get_data_pointer(r2_new)
!--
!--  for(int i=0  i<10  i++)
!--  {
!--    EXPECT_EQ(r2_new_ptr[i], 3)
!--  }
!--
!--  for(int i=10  i<12  i++)
!--  {
!--    EXPECT_EQ(r2_new_ptr[i], 0)     ! assumes zero-ed alloc
!--  }
!--
!--
!--  ! update the other views
!--  r0_new->getNode().update(r0_old->getNode())
!--  r1_new->getNode().update(r1_old->getNode())
!--  r3_new->getNode().update(r3_old->getNode())
!--
!--  call buff_new%print()
!--#endif

    call ds%print()
    call ds%delete()

  end subroutine int_array_multi_view_resize

!------------------------------------------------------------------------------

  subroutine int_array_realloc()
    !
    ! info
    !
    type(SidreDataStore) ds
    type(SidreGroup) root
    type(SidreView) a1, a2
!    type(C_PTR) a1_ptr, a2_ptr
    real(C_FLOAT), pointer :: a1_data(:)
    integer(C_INT), pointer :: a2_data(:)
    integer i

    call set_case_name("int_array_realloc")

    ! create our main data store
    ds = SidreDataStore()

    ! get access to our root data Group
    root = ds%get_root()

    ! create a view to hold the base buffer
    a1 = root%create_view_and_allocate("a1", SIDRE_FLOAT_ID, 5)
    a2 = root%create_view_and_allocate("a2", SIDRE_FLOAT_ID, 5)

    call a1%get_data(a1_data)
    call a2%get_data(a2_data)

    call assert_true(size(a1_data) == 5, &
         "a1_data is incorrect size")
    call assert_true(size(a2_data) == 5, &
         "a2_data is incorrect size")

!--  EXPECT_EQ(SIDRE_view_get_total_bytes(a1), sizeof(float)*5)
!--  EXPECT_EQ(SIDRE_view_get_total_bytes(a2), sizeof(int)*5)

    a1_data(1:5) =  5.0
    a2_data(1:5) = -5

    call a1%reallocate(10)
    call a2%reallocate(15)

    call a1%get_data(a1_data)
    call a2%get_data(a2_data)

    call assert_true(size(a1_data) == 10, &
         "a1_data is incorrect size after realloc")
    call assert_true(size(a2_data) == 15, &
         "a2_data is incorrect size after realloc")

    call assert_true(all(a1_data(1:5) == 5.0), &
         "a1_data does not compare after realloc")
    call assert_true(all(a2_data(1:5) == -5), &
         "a2_data does not compare after realloc")

    a1_data(6:10) = 10.0
    a2_data(6:10) = -10

    a2_data(11:15) = -15

!--  EXPECT_EQ(SIDRE_view_get_total_bytes(a1), sizeof(float)*10)
!--  EXPECT_EQ(SIDRE_view_get_total_bytes(a2), sizeof(int)*15)

    call ds%print()
    call ds%delete()

  end subroutine int_array_realloc

!------------------------------------------------------------------------------

  subroutine simple_opaque()
    type(SidreDataStore) ds
    type(SidreGroup) root
    type(SidreView) opq_view
    integer(C_INT), target :: src_data
    integer(C_INT), pointer :: out_data
    type(C_PTR) src_ptr, opq_ptr

    call set_case_name("simple_opaque")

    ! create our main data store
    ds = SidreDataStore()

    ! get access to our root data Group
    root = ds%get_root()

    src_data = 42
   
    src_ptr = c_loc(src_data)

    opq_view = root%create_view_external("my_opaque", src_ptr)

    ! External pointers are held in the view, not buffer
    call assert_true(ds%get_num_buffers() == 0)

    call assert_true(opq_view%is_external(), "opq_view%is_external()")
    call assert_false(opq_view%is_applied(), "opq_view%is_applied()")
    call assert_true(opq_view%is_opaque(), "opq_view%is_opaque()")

    opq_ptr = opq_view%get_void_ptr()
    call c_f_pointer(opq_ptr, out_data)

    call assert_true(c_associated(opq_ptr, src_ptr), "c_associated(opq_ptr,src_ptr)")
    call assert_equals(out_data, 42, "out_data, 42")

    call ds%print()
    call ds%delete()
!--  free(src_data)
  end subroutine simple_opaque

!------------------------------------------------------------------------------

  subroutine clear_view
    type(SidreDataStore) ds
    type(SidreGroup) root
    integer, parameter :: BLEN = 10

    ds = SidreDataStore()
    root = ds%get_root()

    call v_empty
    call v_described
    call v_scalar
    call v_string
    call v_allocated
    call v_undescribed_buffer
    call v_buffer
    call v_external
    call v_opaque

    call ds%delete()
    
  contains 
    ! Create an empty view.
    subroutine v_empty
      type(SidreView) view

      view = root%create_view("v_empty")
      call check_view_values(view, EMPTYVIEW, .false., .false., .false., 0)
      call view%clear
      call check_view_values(view, EMPTYVIEW, .false., .false., .false., 0)
    end subroutine v_empty

    ! Describe an empty view.
    subroutine v_described
      type(SidreView) view

      view = root%create_view("v_described", SIDRE_INT_ID, BLEN)
      call check_view_values(view, EMPTYVIEW, .true., .false., .false., BLEN)
      call view%clear
      call check_view_values(view, EMPTYVIEW, .false., .false., .false., 0)
    end subroutine v_described

    ! scalar view.
    subroutine v_scalar
      type(SidreView) view

      view = root%create_view_scalar("v_scalar", 1)
      call check_view_values(view, SCALARVIEW, .true., .true., .true., 1)
      call view%clear
      call check_view_values(view, EMPTYVIEW, .false., .false., .false., 0)
    end subroutine v_scalar

    ! string view.
    subroutine v_string
      type(SidreView) view

      view = root%create_view_string("v_string", "string-test")
      call view%clear
      call check_view_values(view, EMPTYVIEW, .false., .false., .false., 0)
    end subroutine v_string

    ! Allocated view, Buffer will be released.
    subroutine v_allocated
      type(SidreView) view
      integer nbuf

      nbuf = ds%get_num_buffers()
      view = root%create_view_and_allocate("v_allocated", SIDRE_INT_ID, BLEN)
      call check_view_values(view, BUFFERVIEW, .true., .true., .true., BLEN)
      call view%clear
      call check_view_values(view, EMPTYVIEW, .false., .false., .false., 0)
      call assert_true(ds%get_num_buffers() == nbuf)
    end subroutine v_allocated

    ! Undescribed Buffer.
    subroutine v_undescribed_buffer
      type(SidreView) view
      type(SidreBuffer) dbuff
      integer nbuf

      nbuf = ds%get_num_buffers()
      dbuff = ds%create_buffer()
      view = root%create_view("v_undescribed_buffer", dbuff)
      call check_view_values(view, BUFFERVIEW, .false., .false., .false., 0)
      call view%clear
      call check_view_values(view, EMPTYVIEW, .false., .false., .false., 0)
      call assert_true(ds%get_num_buffers() == nbuf)
    end subroutine v_undescribed_buffer

    ! Explicit Buffer attached to two views.
    ! Buffer will not be released.
    subroutine v_buffer
      type(SidreView) view, vother
      type(SidreBuffer) dbuff
      integer nbuf

      dbuff = ds%create_buffer()
      call dbuff%allocate(SIDRE_INT_ID, BLEN)
      nbuf = ds%get_num_buffers()
      call assert_true(dbuff%get_num_views() == 0)

      vother = root%create_view("v_other", SIDRE_INT_ID, BLEN)
      view = root%create_view("v_buffer", SIDRE_INT_ID, BLEN)
      call vother%attach_buffer(dbuff)
      call assert_true(dbuff%get_num_views() == 1)
      call view%attach_buffer(dbuff)
      call assert_true(dbuff%get_num_views() == 2)
      
      call check_view_values(view, BUFFERVIEW, .true., .true., .true., BLEN)
      call view%clear
      call check_view_values(view, EMPTYVIEW, .false., .false., .false., 0)
      
      call assert_true(ds%get_num_buffers() == nbuf)
      call assert_true(dbuff%get_num_views() == 1)
    end subroutine v_buffer

    ! External view.
    subroutine v_external
      type(SidreView) view
      integer(C_INT), target :: extdata(BLEN)
      type(C_PTR) ptr

      ptr = c_loc(extdata)
      view = root%create_view("v_external", SIDRE_INT_ID, BLEN, ptr)
      call check_view_values(view, EXTERNALVIEW, .true., .true., .true., BLEN)
      call view%clear
      call check_view_values(view, EMPTYVIEW, .false., .false., .false., 0)
    end subroutine v_external

    ! Opaque view.
    subroutine v_opaque
      type(SidreView) view
      integer(C_INT), target :: extdata(BLEN)
      type(C_PTR) ptr

      ptr = c_loc(extdata)
      view = root%create_view("v_opaque", ptr)
      call check_view_values(view, EXTERNALVIEW, .false., .true., .false., 0)
      call view%clear
      call check_view_values(view, EMPTYVIEW, .false., .false., .false., 0)
    end subroutine v_opaque
  end subroutine clear_view
  
!----------------------------------------------------------------------
end module sidre_view_test
!----------------------------------------------------------------------

program fortran_test
  use fruit
  use sidre_view_test
  implicit none
  logical ok

  call init_fruit

  call create_views
  call get_path_name
! create_view_from_path
  call scalar_view
! dealloc
! alloc_zero_items
! alloc_and_dealloc_multiview
  call int_buffer_from_view
  call int_buffer_from_view_conduit_value
  call int_array_multi_view
  call init_int_array_multi_view
  call int_array_depth_view
  call int_array_view_attach_buffer
  call int_array_offset_stride
  call int_array_multi_view_resize
  call int_array_realloc
  call simple_opaque
  call clear_view

  call fruit_summary
  call fruit_finalize

  call is_all_successful(ok)
  if (.not. ok) then
     call exit(1)
  endif
end program fortran_test

