! Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
! other Axom Project Developers. See the top-level LICENSE file for details.
!
! SPDX-License-Identifier: (BSD-3-Clause)

module sidre_external_test
  use iso_c_binding
  use fruit
  use axom_sidre
  implicit none

contains

!------------------------------------------------------------------------------
! Test Group::createExternalView()
!------------------------------------------------------------------------------
  subroutine create_external_view
    type(SidreDataStore) ds
    type(SidreGroup) root
    type(SidreView) iview, dview
    integer(C_INT), allocatable, target :: idata(:)
    real(C_DOUBLE), allocatable, target :: ddata(:)
    integer(C_INT), pointer :: idata_chk(:)
    real(C_DOUBLE), pointer :: ddata_chk(:)
    integer(SIDRE_IndexType), parameter :: len = 11
    integer ii

    call set_case_name("create_external_view")

    ds = SidreDataStore()
    root = ds%get_root()

    allocate(idata(len))
    allocate(ddata(len))

    do ii = 1, len
       idata(ii) = ii
       ddata(ii) = idata(ii) * 2.0
    enddo

    iview = root%create_view_external("idata", c_loc(idata))
    call iview%apply(SIDRE_INT_ID, len)
    dview = root%create_view_external("ddata", c_loc(ddata))
    call dview%apply(SIDRE_DOUBLE_ID, len)
    call assert_true(root%get_num_views() .eq. 2)

    call iview%print()
    call dview%print()

    call iview%get_data(idata_chk)
    call assert_true(size(idata_chk) == len, "idata_chk is wrong size")
    call assert_true(all(idata_chk == idata), "idata_chk != idata")

    call dview%get_data(ddata_chk)
    call assert_true(size(ddata_chk) == len, "ddata_chk is wrong size")
    call assert_true(all(ddata_chk == ddata), "ddata_chk != ddata")

    call ds%delete()
    deallocate(idata)
    deallocate(ddata)
  end subroutine create_external_view

!------------------------------------------------------------------------------
! Test Group::save(), Group::load() with external buffers
!------------------------------------------------------------------------------
  subroutine save_load_external_view
    type(SidreDataStore) ds, ds2
    type(SidreGroup) root, root2
    type(SidreView)  iview, dview, iview2, dview2
    type(SidreBuffer) tmpbuff
    integer(C_INT), allocatable, target :: idata(:)
    real(C_DOUBLE), allocatable, target :: ddata(:)
    integer(C_INT), pointer :: idata_chk(:)
    real(C_DOUBLE), pointer :: ddata_chk(:)
    integer(SIDRE_IndexType), parameter :: len = 11
    integer ii

    call set_case_name("save_load_external_view")

    ds = SidreDataStore()
    root = ds%get_root()

    allocate(idata(len))
    allocate(ddata(len))

    do ii = 1, len
       idata(ii) = ii
       ddata(ii) = idata(ii) * 2.0
    enddo

    iview = root%create_view_external("idata", c_loc(idata))
    call iview%apply_type_nelems(SIDRE_INT_ID, len)
    dview = root%create_view_external("ddata", c_loc(ddata))
    call dview%apply_type_nelems(SIDRE_DOUBLE_ID, len)

    call assert_true(root%get_num_views() .eq. 2)

    call iview%print()
    call dview%print()

!  TODO - fix wrapping, change to datastore save call
    !call root%save("out_sidre_external_save_restore_external_view", "conduit")

    call ds%print()


    ds2 = SidreDataStore()
    root2 = ds2%get_root()

! TODO - fix wrapping change to datastore load call
    !call root2%load("out_sidre_external_save_restore_external_view","conduit")

    call ds2%print()

    call assert_true(root2%get_num_views() .eq. 2)

    iview = root2%get_view("idata")
    dview = root2%get_view("ddata")

    call iview%get_data(idata_chk)
    call assert_true(size(idata_chk) == len, "idata_chk is wrong size")
    call assert_true(all(idata_chk == idata), "idata_chk != idata")

    call dview%get_data(ddata_chk)
    call assert_true(size(ddata_chk) == len, "ddata_chk is wrong size")
    call assert_true(all(ddata_chk == ddata), "ddata_chk != ddata")

    call ds%delete()
    call ds2%delete()
    deallocate(idata)
    deallocate(ddata)
  end subroutine save_load_external_view

!----------------------------------------------------------------------
end module sidre_external_test
!----------------------------------------------------------------------

program fortran_test
  use fruit
  use sidre_external_test
  implicit none
  logical ok

  call init_fruit

  call create_external_view
!  Needs to be fixed
!  call save_load_external_view


  call fruit_summary
  call fruit_finalize

  call is_all_successful(ok)
  if (.not. ok) then
     call exit(1)
  endif
end program fortran_test
