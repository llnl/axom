! Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
! other Axom Project Developers. See the top-level LICENSE file for details.
!
! SPDX-License-Identifier: (BSD-3-Clause)

module sidre_buffer_test
  use iso_c_binding
  use fruit
  use axom_sidre
  implicit none

contains
!------------------------------------------------------------------------------

  subroutine create_buffers
    type(SidreDataStore) ds
    type(SidreBuffer) dbuff_0, dbuff_1, dbuff_3

    call set_case_name("create_buffers")

    ds = SidreDataStore()

    dbuff_0 = ds%create_buffer()
    dbuff_1 = ds%create_buffer()

    call assert_true(dbuff_0%get_index() == 0, "dbuff_0%get_index()")
    call assert_true(dbuff_1%get_index() == 1, "dbuff_1%get_index()")
    call ds%destroy_buffer(0)

    dbuff_3 = ds%create_buffer()
    call assert_true(dbuff_3%get_index() == 0, "dbuff_3%get_index()")

    call ds%print()
    call ds%delete()
  end subroutine create_buffers

!------------------------------------------------------------------------------

  subroutine alloc_buffer_for_int_array
    type(SidreDataStore) ds
    type(SidreBuffer) dbuff
    type(C_PTR) data_ptr
    integer(C_INT), pointer :: data(:)
    integer(C_INT) i
    integer int_size, elem_count

    int_size = c_sizeof(i)
    elem_count = 10
    
    call set_case_name("alloc_buffer_for_int_array")

    ds = SidreDataStore()
    dbuff = ds%create_buffer()

    call dbuff%allocate(SIDRE_INT_ID, elem_count)

    call assert_true(dbuff%get_type_id() == SIDRE_INT_ID, "dbuff%get_typeid()")
    call assert_true(dbuff%get_num_elements() == elem_count, "dbuff%get_num_elements()")
    call assert_true(dbuff%get_bytes_per_element() == int_size)
    call assert_true(dbuff%get_total_bytes() == int_size * elem_count)

    data_ptr = dbuff%get_void_ptr()
    call c_f_pointer(data_ptr, data, [ elem_count ])

    do i = 1, elem_count
       data(i) = i * i
    enddo

    call dbuff%print()

    call ds%print()
    call ds%delete()
  end subroutine alloc_buffer_for_int_array

!------------------------------------------------------------------------------

  subroutine init_buffer_for_int_array
    type(SidreDataStore) ds
    type(SidreBuffer) dbuff
    type(C_PTR) data_ptr
    integer(C_INT), pointer :: data(:)
    integer(C_INT) i
    integer int_size, elem_count

    int_size = c_sizeof(i)
    elem_count = 10

    call set_case_name("init_buffer_for_int_array")

    ds = SidreDataStore()
    dbuff = ds%create_buffer()

    call dbuff%allocate(SIDRE_INT_ID, elem_count)

    call assert_true(dbuff%get_type_id() == SIDRE_INT_ID, "dbuff%get_type_id()")
    call assert_true(dbuff%get_num_elements() == elem_count, "dbuff%get_num_elements")
    call assert_true(dbuff%get_bytes_per_element() == int_size)
    call assert_true(dbuff%get_total_bytes() == int_size * elem_count)

    data_ptr = dbuff%get_void_ptr()
    call c_f_pointer(data_ptr, data, [ elem_count ])

    do i = 1, elem_count
       data(i) = i * i
    enddo

    call dbuff%print()

    call ds%print()
    call ds%delete()
  end subroutine init_buffer_for_int_array

!------------------------------------------------------------------------------

  subroutine realloc_buffer
    type(SidreDataStore) ds
    type(SidreBuffer) dbuff
    type(C_PTR) data_ptr
    integer(C_LONG), pointer :: data(:)
    integer i
    integer(C_LONG) elem
    integer long_size, orig_elem_count, mod_elem_count

    long_size = c_sizeof(elem)
    orig_elem_count = 5
    mod_elem_count = 10

    call set_case_name("realloc_buffer")

    ds = SidreDataStore()

    dbuff = ds%create_buffer()

    call dbuff%allocate(SIDRE_LONG_ID, orig_elem_count)

    call assert_true(dbuff%get_type_id() == SIDRE_LONG_ID, "dbuff%get_type_id()")
    call assert_true(dbuff%get_num_elements() == orig_elem_count, "dbuff%get_num_elements()")
    call assert_true(dbuff%get_bytes_per_element() == long_size)
    call assert_true(dbuff%get_total_bytes() == long_size * orig_elem_count)

    data_ptr = dbuff%get_void_ptr()
    call c_f_pointer(data_ptr, data, [ orig_elem_count ])

    data(:) = orig_elem_count

    call dbuff%print()
  
    call dbuff%reallocate( mod_elem_count )

    call assert_true(dbuff%get_type_id() == SIDRE_LONG_ID, "dbuff%get_type_id() after realloc")
    call assert_true(dbuff%get_num_elements() == mod_elem_count, "dbuff%get_num_elements() after realloc")
    call assert_true(dbuff%get_bytes_per_element() == long_size)
    call assert_true(dbuff%get_total_bytes() == long_size * mod_elem_count)

    ! data buffer changes
    data_ptr = dbuff%get_void_ptr()
    call c_f_pointer(data_ptr, data, [ mod_elem_count ])

    call assert_true(size(data) == mod_elem_count, "data_ptr is wrong size")
    call assert_true(all(data(1:5) == orig_elem_count), "data has wrong values")

    data(6:10) = mod_elem_count

    call dbuff%print()

    call ds%print()
    call ds%delete()
  end subroutine realloc_buffer

!----------------------------------------------------------------------
end module sidre_buffer_test
!----------------------------------------------------------------------

program fortran_test
  use fruit
  use sidre_buffer_test
  implicit none
  logical ok

  call init_fruit

  call create_buffers
  call alloc_buffer_for_int_array
  call init_buffer_for_int_array
  call realloc_buffer

  call fruit_summary
  call fruit_finalize

  call is_all_successful(ok)
  if (.not. ok) then
     call exit(1)
  endif
end program fortran_test
