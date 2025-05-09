! Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
! other Axom Project Developers. See the top-level LICENSE file for details.
!
! SPDX-License-Identifier: (BSD-3-Clause)

!
! Test allocatable arrays as meta-buffers
!

module sidre_allocatable_test
  use iso_c_binding
  use fruit
  use axom_sidre
  implicit none

contains

! Allocate array via Fortran
! Register with datastore then 
! Query metadata using datastore API.
!----------------------------------------------------------------------

  subroutine external_allocatable_int
    integer, allocatable :: iarray(:)
    integer, pointer :: ipointer(:)

    type(SidreDataStore) ds
    type(SidreGroup) root
    type(SidreView)  view
    integer type
    integer num_elements
    integer i
    integer rank
    integer(SIDRE_IndexType) extents(7)

    call set_case_name("external_allocatable_int")

    ds = SidreDataStore()
    root = ds%get_root()

    allocate(iarray(10))

    do i=1,10
       iarray(i) = i
    enddo

    view = root%create_array_view("iarray", iarray)

    call assert_true(view%is_external())

    type = view%get_type_id()
    call assert_true(type == SIDRE_INT_ID)

    num_elements = view%get_num_elements()
    call assert_equals(num_elements, size(iarray))

    rank = view%get_num_dimensions()
    call assert_equals(rank, 1)

    rank = view%get_shape(7, extents)
    call assert_equals(rank, 1)
    call assert_true(extents(1) == size(iarray, 1))

    ! get array via a pointer
    call view%get_data(ipointer)
    call assert_true(all(iarray.eq.ipointer))

    call ds%delete()

    deallocate(iarray)

  end subroutine external_allocatable_int

!----------------------------------------------------------------------

  subroutine external_allocatable_int_3d
    integer, allocatable :: iarray(:,:,:)
    integer, pointer :: ipointer(:,:,:)

    type(SidreDataStore) ds
    type(SidreGroup) root
    type(SidreView)  view
    integer type
    integer num_elements
    integer i, j, k
    integer rank
    integer(SIDRE_IndexType) extents(7)

    call set_case_name("external_allocatable_int_3d")

    ds = SidreDataStore()
    root = ds%get_root()

    allocate(iarray(2,3,4))

    do i=1,2
       do j=1,3
          do k=1,4
             iarray(i,j,k) = i*100 + j*1-0 + k
          enddo
       enddo
    enddo

    view = root%create_array_view("iarray", iarray)

    call assert_true(view%is_external())

    type = view%get_type_id()
    call assert_true(type == SIDRE_INT_ID)

    num_elements = view%get_num_elements()
    call assert_equals(num_elements, size(iarray))

    rank = view%get_num_dimensions()
    call assert_equals(rank, 3)

    rank = view%get_shape(7, extents)
    call assert_equals(rank, 3)
    call assert_true(extents(1) == size(iarray, 1))
    call assert_true(extents(2) == size(iarray, 2))
    call assert_true(extents(3) == size(iarray, 3))

    ! get array via a pointer
    call view%get_data(ipointer)
    call assert_true(all(iarray.eq.ipointer))

    call ds%delete()

    deallocate(iarray)

  end subroutine external_allocatable_int_3d

!----------------------------------------------------------------------
!
! register a static (non-allocatable) array with the datastore as external view

  subroutine external_static_int
    integer :: iarray(10)
    integer, pointer :: ipointer(:)

    type(SidreDataStore) ds
    type(SidreGroup) root
    type(SidreView)  view
    integer(TypeID) type
    integer num_elements
    integer i

    call set_case_name("external_static_int")

    ds = SidreDataStore()
    root = ds%get_root()

    do i=1,10
       iarray(i) = i
    enddo

    view = root%create_array_view("iarray", iarray)

    type = view%get_type_id()
    call assert_true(type == SIDRE_INT_ID)

    num_elements = view%get_num_elements()
    call assert_equals(num_elements, 10)

    ! get array via a pointer
    call view%get_data(ipointer)
    call assert_true(all(iarray.eq.ipointer))

    call ds%delete()

  end subroutine external_static_int

!----------------------------------------------------------------------
!--- check other types

  subroutine external_allocatable_double
    real(C_DOUBLE), allocatable :: darray(:)
    real(C_DOUBLE), pointer :: dpointer(:)

    type(SidreDataStore) ds
    type(SidreGroup) root
    type(SidreView)  view
    integer num_elements
    integer(TypeID) type
    integer i

    call set_case_name("external_allocatable_double")

    ds = SidreDataStore()
    root = ds%get_root()

    allocate(darray(10))

    do i=1,10
       darray(i) = i + 0.5d0
    enddo

    view = root%create_array_view("darray", darray)

    type = view%get_type_id()
    call assert_true(type == SIDRE_DOUBLE_ID)

    num_elements = view%get_num_elements()
    call assert_equals(num_elements, 10)

    ! get array via a pointer
    call view%get_data(dpointer)
    call assert_true(all(abs(darray-dpointer).lt..0005))

    call ds%delete()

    deallocate(darray)

  end subroutine external_allocatable_double

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Datastore owns a multi-dimension array.

  subroutine datastore_int_3d
    integer, pointer :: ipointer(:,:,:)

    type(SidreDataStore) ds
    type(SidreGroup) root
    type(SidreView)  view
    integer(TypeID) type
    integer num_elements
    integer i, j, k
    integer rank
    integer(SIDRE_IndexType) extents_in(3), extents(7)

    call set_case_name("datastore_int_3d")

    extents_in(1) = 2
    extents_in(2) = 3
    extents_in(3) = 4

    ds = SidreDataStore()
    root = ds%get_root()

    view = root%create_view_with_shape_and_allocate("iarray", SIDRE_INT_ID, 3, extents_in)

    call view%get_data(ipointer)

    type = view%get_type_id()
    call assert_true(type == SIDRE_INT_ID)

    num_elements = view%get_num_elements()
    call assert_equals(num_elements, size(ipointer))

    rank = view%get_num_dimensions()
    call assert_equals(rank, 3)

    rank = view%get_shape(7, extents)
    call assert_equals(rank, 3)
    call assert_true(extents(1) == size(ipointer, 1))
    call assert_true(extents(2) == size(ipointer, 2))
    call assert_true(extents(3) == size(ipointer, 3))

    ! reshape as 1-d using shape
    extents_in(1) = size(ipointer)
    call view%apply(SIDRE_INT_ID, 1, extents_in(1:1))
    num_elements = view%get_num_elements()
    call assert_equals(num_elements, size(ipointer))

    ! reshape as 1-d using length
    call view%apply(SIDRE_INT_ID, extents_in(1))
    num_elements = view%get_num_elements()
    call assert_equals(num_elements, size(ipointer))

    call ds%delete()

  end subroutine datastore_int_3d

!----------------------------------------------------------------------
!----------------------------------------------------------------------

end module sidre_allocatable_test


program fortran_test
  use fruit
  use sidre_allocatable_test
  implicit none
  logical ok

  call init_fruit

  call external_allocatable_int
  call external_allocatable_int_3d
  call external_static_int
  call external_allocatable_double
  call datastore_int_3d

  call fruit_summary
  call fruit_finalize

  call is_all_successful(ok)
  if (.not. ok) then
     call exit(1)
  endif
end program fortran_test
