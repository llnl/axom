! Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
! other Axom Project Developers. See the top-level LICENSE file for details.
!
! SPDX-License-Identifier: (BSD-3-Clause)

! API coverage tests
! Each test should be documented with the interface functions being tested

#include "axom/config.hpp"

module sidre_group_test
  use iso_c_binding
  use fruit
  use axom_sidre
  implicit none

#ifdef AXOM_USE_HDF5
  ! Test protocols
  integer, parameter :: nprotocols = 3
  character(12) :: protocols(nprotocols) = [ "sidre_json  ", "sidre_hdf5  ", "json        " ]
#else
  ! Test protocols
  integer, parameter :: nprotocols = 2
  character(12) :: protocols(nprotocols) = [ "sidre_json  ", "json        " ]
#endif

contains

  !------------------------------------------------------------------------------
  ! get_name()
  !------------------------------------------------------------------------------
  subroutine get_name
    type(SidreDataStore) ds
    type(SidreGroup) root, grp, grp2

    call set_case_name("get_name")

    ds = datastore_new()
    root = ds%get_root()
    grp = root%create_group("test")

    call assert_true(grp%get_name() == "test" )

    grp2 = root%get_group("foo")
    call assert_false( grp2%associated() )
    
    call ds%delete()
  end subroutine get_name

  !------------------------------------------------------------------------------
  ! get_path_name()
  !------------------------------------------------------------------------------
  subroutine get_path_name
    type(SidreDataStore) ds
    type(SidreGroup) root, group, grp2, grp3

    call set_case_name("get_path_name")

    ds = datastore_new()
    root = ds%get_root()
    group = root%create_group("test/a/b/c")
    grp2 = root%get_group("test/a")
    grp3 = root%get_group("test")

    call assert_true(root%get_name() == " " )
    call assert_true(root%get_path() == " " )
    call assert_true(root%get_path_name() == " " )

    call assert_true(grp2%get_name() == "a" )
    call assert_true(grp2%get_path() == "test" )
    call assert_true(grp2%get_path_name() == "test/a" )

    call assert_true(grp3%get_name() == "test" )
    call assert_true(grp3%get_path() == " " )
    call assert_true(grp3%get_path_name() == "test" )

    call assert_true(group%get_name() == "c" )
    call assert_true(group%get_path() == "test/a/b" )
    call assert_true(group%get_path_name() == "test/a/b/c" )

    call ds%delete()
  end subroutine get_path_name

  !------------------------------------------------------------------------------
  ! get_parent()
  !------------------------------------------------------------------------------
  subroutine get_parent
    type(SidreDataStore) ds
    type(SidreGroup) root, parent, child

    call set_case_name("get_parent")

    ds = datastore_new()
    root = ds%get_root()
    parent = root%create_group("parent")
    child = parent%create_group("child")

    call assert_true( child%get_parent() == parent )

    call ds%delete()
  end subroutine get_parent

!------------------------------------------------------------------------------
! Verify get_data_store()
!------------------------------------------------------------------------------
  subroutine get_datastore
    type(SidreDataStore) ds, const_ds
    type(SidreGroup) root, grp

    call set_case_name("get_datastore")

    ds = datastore_new()
    root = ds%get_root()
    grp = root%create_group("parent")

    call assert_true( grp%get_data_store() == ds )

    const_ds = grp%get_data_store()
    call assert_true( const_ds == ds )

    call ds%delete()
  end subroutine get_datastore

  !------------------------------------------------------------------------------
  ! Verify get_group()
  !------------------------------------------------------------------------------
  subroutine get_group
    type(SidreDataStore) ds
    type(SidreGroup) root, parent, child, child1, child2, errgrp

    call set_case_name("has_group")

    ds = datastore_new()
    root = ds%get_root()

    parent = root%create_group("parent")
    child = parent%create_group("child")
    call assert_true( child%get_parent() == parent )

    child1 = parent%get_group("child")
    call assert_true( child == child1 )

    child2 = parent%get_group(0)
    call assert_true( child == child2 )

    ! check error condition
    errgrp = parent%get_group("non-existent group")
    call assert_false( c_associated(errgrp%get_instance()) )

    call ds%delete()
  end subroutine get_group

!------------------------------------------------------------------------------
! Verify get_view()
!------------------------------------------------------------------------------
  subroutine get_view
    type(SidreDataStore) ds
    type(SidreGroup) root,parent
    type(SidreView) view1, view2, view3

    call set_case_name("get_view")

    ds = datastore_new()
    root = ds%get_root()

    parent = root%create_group("parent")
    view1 = parent%create_view("view")

    view2 = parent%get_view("view")
    call assert_true( view1 == view2 )

    view3 = parent%get_view(0)
    call assert_true( view1 == view3 )

    ! check error condition
    view2 = parent%get_view("non-existent view")
    call assert_false( view2%associated() )

    call ds%delete()
  end subroutine get_view

!------------------------------------------------------------------------------
! Verify get_view_name(), get_view_index()
!------------------------------------------------------------------------------
  subroutine get_view_name_index
    type(SidreDataStore) ds
    type(SidreGroup) root, parent
    type(SidreView) view1, view2
    integer(SIDRE_IndexType) idx1, idx2, idx3
    character(len=MAXNAMESIZE) name1, name2, name3

    call set_case_name("get_view_name_index")

    ds = datastore_new()
    root = ds%get_root()

    parent = root%create_group("parent")
    view1 = parent%create_view("view1")
    view2 = parent%create_view("view2")

    call assert_true(parent%get_num_views() == 2)

    idx1 = parent%get_view_index("view1")
    idx2 = parent%get_view_index("view2")

    name1 = parent%get_view_name(idx1)
    name2 = parent%get_view_name(idx2)

    call assert_equals(name1, "view1")
    call assert_equals(view1%get_name(), name1)
    
    call assert_equals(name2, "view2")
    call assert_equals(view2%get_name(), name2)

    idx3 = parent%get_view_index("view3")
    call assert_true(idx3 == invalid_index)

    name3 = parent%get_view_name(idx3)
    call assert_true(name3 == " ")
    call assert_false(name_is_valid(name3))

    call ds%delete()
  end subroutine get_view_name_index

!------------------------------------------------------------------------------
! Verify getFirstValidGroupIndex, getNextValidGroupIndex
!------------------------------------------------------------------------------
  subroutine get_first_and_next_group_index
    type(SidreDataStore) ds
    type(SidreGroup) root, parent, emptygrp
    type(SidreGroup) group1, group2, group1out, group2out
    integer(SIDRE_IndexType) idx1, idx2, idx3, badidx1, badidx2

    call set_case_name("get_first_and_next_group_index")

    ds = datastore_new()
    root = ds%get_root()

    parent = root%create_group("parent")
    group1 = parent%create_group("group1")
    group2 = parent%create_group("group2")
    call assert_true(parent%get_num_groups() == 2)

    idx1 = parent%get_first_valid_group_index()
    idx2 = parent%get_next_valid_group_index(idx1)
    idx3 = parent%get_next_valid_group_index(idx2)
    call assert_true(0 == idx1)
    call assert_true(1 == idx2)
    call assert_true(invalid_index == idx3)

    group1out = parent%get_group(idx1)
    group2out = parent%get_group(idx2)
    call assert_true(group1 == group1out)
    call assert_true(group2 == group2out)

    ! check error conditions
    emptygrp = root%create_group("emptyGroup")
    badidx1 = emptygrp%get_first_valid_group_index()
    badidx2 = emptygrp%get_next_valid_group_index(badidx1)

    call assert_true(badidx1 == invalid_index)
    call assert_true(badidx2 == invalid_index)

    call ds%delete()
  end subroutine get_first_and_next_group_index

!------------------------------------------------------------------------------
! Verify getFirstValidViewIndex, getNextValidViewIndex
!------------------------------------------------------------------------------
  subroutine get_first_and_next_view_index
    type(SidreDataStore) ds
    type(SidreGroup) root, parent, emptygrp
    type(SidreView) view1, view2, view1out, view2out
    integer(SIDRE_IndexType) idx1, idx2, idx3, badidx1, badidx2

    call set_case_name("get_first_and_next_view_index")

    ds = datastore_new()
    root = ds%get_root()

    parent = root%create_group("parent")
    view1 = parent%create_view("view1")
    view2 = parent%create_view("view2")
    call assert_true(parent%get_num_views() == 2)

    idx1 = parent%get_first_valid_view_index()
    idx2 = parent%get_next_valid_view_index(idx1)
    idx3 = parent%get_next_valid_view_index(idx2)
    call assert_true(0 == idx1)
    call assert_true(1 == idx2)
    call assert_true(invalid_index == idx3)

    view1out = parent%get_view(idx1)
    view2out = parent%get_view(idx2)
    call assert_true(view1 == view1out)
    call assert_true(view2 == view2out)

    ! check error conditions
    emptygrp = root%create_group("emptyGroup")
    badidx1 = emptygrp%get_first_valid_view_index()
    badidx2 = emptygrp%get_next_valid_view_index(badidx1)

    call assert_true(badidx1 == invalid_index)
    call assert_true(badidx2 == invalid_index)

    call ds%delete()
  end subroutine get_first_and_next_view_index

!------------------------------------------------------------------------------
! Verify get_group_name(), get_group_index()
!------------------------------------------------------------------------------
  subroutine get_group_name_index
    type(SidreDataStore) ds
    type(SidreGroup) root, parent, grp1, grp2
    integer(SIDRE_IndexType) idx1, idx2, idx3
    character(len=MAXNAMESIZE) name1, name2, name3

    call set_case_name("get_group_name_index")

    ds = datastore_new()
    root = ds%get_root()

    parent = root%create_group("parent")
    grp1 = parent%create_group("grp1")
    grp2 = parent%create_group("grp2")

!--    call assert_equals(parent%get_num_groups(), 2)  ! size_t
    call assert_true(parent%get_num_groups() == 2)

    idx1 = parent%get_group_index("grp1")
    idx2 = parent%get_group_index("grp2")

    name1 = parent%get_group_name(idx1)
    name2 = parent%get_group_name(idx2)

    call assert_equals(name1, "grp1")
    call assert_equals(grp1%get_name(), name1)

    call assert_equals(name2, "grp2")
    call assert_equals(grp2%get_name(), name2)

    idx3 = parent%get_group_index("grp3")
    call assert_true(idx3 == invalid_index)

    name3 = parent%get_group_name(idx3)
    call assert_true(name3 == " ")
    call assert_false(name_is_valid(name3))

    call ds%delete()
  end subroutine get_group_name_index

  !------------------------------------------------------------------------------
  ! create_view_empty()
  ! destroy_view()
  ! has_view()
  !------------------------------------------------------------------------------
  subroutine create_destroy_has_view
    type(SidreDataStore) ds
    type(SidreGroup) root,grp
    type(SidreView) view0, view1, view2

    call set_case_name("create_destroy_has_view")

    ds = datastore_new()
    root = ds%get_root()
    grp = root%create_group("parent")

    view0 = grp%create_view("view0")
    call assert_true( grp%get_parent() == root )
    call assert_false( view0%has_buffer() )

    call assert_true( grp%has_view("view0") )
    ! try creating view again, should be a no-op.
    view1 = grp%create_view("view0")
    call assert_false( view1%associated() )

    call grp%destroy_view("view0")

    call assert_false( grp%has_view("view0"), 'grp%has_view("view0")')

    ! try api call that specifies specific type and length
    view1 = grp%create_view_and_allocate( "viewWithLength1", SIDRE_FLOAT_ID, 50 )

    ! error condition check - try again with duplicate name, should be a no-op
!XXX    view1 = grp%create_view_and_allocate( "viewWithLength1", SIDRE_FLOAT64_ID, 50 )
!XXX    call assert_true( view1%associated() )
    call grp%destroy_view_and_data("viewWithLength1")
    call assert_false( grp%has_view("viewWithLength1"), &
         'grp%has_view("viewWithLength1"' )

    view1 = grp%create_view_and_allocate( "viewWithLengthBadLen", SIDRE_FLOAT64_ID, -1 )
    call assert_false( view1%associated(), 'view1%associated()' )

    ! try api call that specifies data type in another way
    view1 = grp%create_view_and_allocate( "viewWithLength2", SIDRE_FLOAT64_ID, 50 )
    view2 = grp%create_view_and_allocate( "viewWithLength2", SIDRE_FLOAT64_ID, 50 )
    call assert_false( view2%associated(), 'view2%associated()')
    ! destroy this view using index
    call grp%destroy_view_and_data( grp%get_first_valid_view_index() )

    call ds%delete()

  end subroutine create_destroy_has_view

  !------------------------------------------------------------------------------
  ! create_group()
  ! destroy_group()
  ! has_group()
  !------------------------------------------------------------------------------
  subroutine create_destroy_has_group
    type(SidreDataStore) ds
    type(SidreGroup) root, grp, grp2

    call set_case_name("create_destroy_has_group")

    ds = datastore_new()
    root = ds%get_root()
    grp = root%create_group("grp")
    call assert_true( grp%get_parent() == root )

    call assert_true( root%has_group("grp") )

    call root%destroy_group("grp")
    call assert_false( root%has_group("grp") )

    grp2 = root%create_group("grp2")
    call root%destroy_group( root%get_first_valid_group_index() )

    call ds%delete()
  end subroutine create_destroy_has_group

  !------------------------------------------------------------------------------
  subroutine group_name_collisions
    type(SidreDataStore) ds
    type(SidreGroup) root, flds, badgrp
    type(SidreView) view

    call set_case_name("group_name_collisions")

    ds = datastore_new()
    root = ds%get_root()
    flds = root%create_group("fields")
    view = flds%create_view("a")

    call assert_true(flds%has_view("a"))

    ! attempt to create duplicate group name

    badgrp = root%create_group("fields")
    call assert_false( badgrp%associated() )

    ! check error condition
    ! attempt to create duplicate view name.
    view = flds%create_view("a")
    call assert_false( view%associated())

    call ds%delete()
  end subroutine group_name_collisions

  !------------------------------------------------------------------------------
! Restore this after copy_move is working. ATK-667
  subroutine view_copy_move
    type(SidreDataStore) ds
    type(SidreGroup) root, flds, subgrp
    type(SidreView) i0_view, f0_view, d0_view, tmpview

    call set_case_name("view_copy_move")

    ds = datastore_new()
    root = ds%get_root()
    flds = root%create_group("fields")

    i0_view = flds%create_view_scalar_int("i0", 1)
    f0_view = flds%create_view_scalar_float("f0", 100.0)
    d0_view = flds%create_view_scalar_double("d0", 3000.0d0)

    call assert_true(flds%has_view("i0"))
    call assert_true(flds%has_view("f0"))
    call assert_true(flds%has_view("d0"))

    ! test moving a view from flds to sub
!    subgrp = flds%create_group("sub")
!    tmpview = subgrp%move_view(flds%get_view("d0"))
!    call flds%print()
!    call assert_false(flds%has_view("d0"))
!    call assert_true(flds%has_group("sub"))
!    call assert_true(subgrp%has_view("d0"))

    ! check the data value
!    call assert_equals(tmpview%get_data_double(), 3000.0_C_DOUBLE)

    ! test copying a view from flds to sub
!    tmpview = subgrp%copy_view(flds%get_view("i0"))

!    call flds%print()

!    call assert_true(flds%has_view("i0"))
!    call assert_true(subgrp%has_view("i0"))

    ! we expect the actual data  pointers to be the same
!--    call assert_equals(flds%get_view("i0")%getBuffer(),
!--    call flds%get_group("sub")%get_view("i0")%getBuffer())
    call ds%delete()
  end subroutine view_copy_move

  !------------------------------------------------------------------------------
  subroutine groups_move_copy
    type(SidreDataStore) ds
    type(SidreGroup) root, flds, ga, gb, gc
    type(SidreGroup) subgrp, tmpgrp
    type(SidreView) i0_view, f0_view, d0_view, tmpview

    call set_case_name("groups_move_copy")

    ds = datastore_new()
    root = ds%get_root()
    flds = root%create_group("fields")

    ga = flds%create_group("a")
    gb = flds%create_group("b")
    gc = flds%create_group("c")

    i0_view = ga%create_view_scalar_int("i0", 1)
    f0_view = gb%create_view_scalar_float("f0", 100.0)
    d0_view = gc%create_view_scalar_double("d0", 3000.0d0)

    ! check that all sub groups exist
    call assert_true(flds%has_group("a"))
    call assert_true(flds%has_group("b"))
    call assert_true(flds%has_group("c"))

    !move "b" to a child of "sub"
    subgrp = flds%create_group("sub")
    tmpgrp = subgrp%move_group(gb)

    call flds%print()

    call assert_true(flds%has_group("a"))
    call assert_true(flds%has_group("sub"))
    call assert_true(flds%has_group("c"))
    
    call assert_true(tmpgrp == gb)
    
    call ds%delete()
  end subroutine groups_move_copy

  !------------------------------------------------------------------------------
  subroutine create_destroy_view_and_data
    type(SidreDataStore) ds
    type(SidreGroup) root, grp
    type(SidreView) view1, view2
    type(SidreBuffer) tmpbuf
!XX    type(SidreBuffer) buffer1, tmpbuf
    integer(SIDRE_IndexType) bufferid1
    character(len=30) view_name1, view_name2
!XX    logical buffvalid

    call set_case_name("create_destroy_view_and_data")

    ds = datastore_new()
    root = ds%get_root()
    grp = root%create_group("grp")

    view_name1 = "viewBuffer1"
    view_name2 = "viewBuffer2"

    view1 = grp%create_view_and_allocate(view_name1, SIDRE_INT_ID, 1_8)
    view2 = grp%create_view_and_allocate(view_name2, SIDRE_INT_ID, 1_8)

    call assert_true(grp%has_view(view_name1))
    call assert_true(grp%get_view(view_name1) == view1)

    call assert_true(grp%has_view(view_name2))
    call assert_true(grp%get_view(view_name2) == view2)

    tmpbuf = view1%get_buffer()
    bufferid1 = tmpbuf%get_index()

    call grp%destroy_view_and_data(view_name1)


    call assert_false(grp%has_view(view_name1))
    call assert_true(ds%get_num_buffers() == 1)

!XX    buffer1 = ds%get_buffer(bufferId1)
!XX   buffvalid = .true.
!--    if( buffer1 == nullptr ) buffValid = .false.
!XX    call assert_false(buffValid)

    call ds%delete()
  end subroutine create_destroy_view_and_data


  !------------------------------------------------------------------------------
  subroutine create_destroy_alloc_view_and_data
    type(SidreDataStore) ds
    type(SidreGroup) root, grp
    type(SidreView) view1
    integer i
    character(len=30) view_name1, view_name2
    integer(C_INT), pointer :: v1_vals(:)
!--    real(C_DOUBLE), pointer :: v2_vals(:)

    call set_case_name("create_destroy_alloc_view_and_data")

    ds = datastore_new()
    root = ds%get_root()
    grp = root%create_group("grp")

    view_name1 = "viewBuffer1"
    view_name2 = "viewBuffer2"

    ! use create + alloc convenience methods
    ! this one is the DataType & method
    view1 = grp%create_view_and_allocate(view_name1, SIDRE_INT_ID, 10)

    call assert_true(grp%has_view(view_name1))
    call assert_true(grp%get_view(view_name1) == view1)

    call view1%get_data(v1_vals)
  
    do i = 1, 10
       v1_vals(i) = i
    enddo

    call assert_true(view1%get_num_elements() == 10)
!--    call assert_equals(view1%get_total_bytes(), 10 * sizeof(int))

    call grp%destroy_view_and_data(view_name1)

    call ds%delete()
  end subroutine create_destroy_alloc_view_and_data

  !------------------------------------------------------------------------------
  subroutine create_view_of_buffer_with_datatype
    type(SidreDataStore) ds
    type(SidreGroup) root
    type(SidreView) base
    type(SidreBuffer) base_buff
    integer(C_INT), pointer :: base_vals(:)
!--    integer i
!    integer(C_INT), pointer :: sub_a_vals(:)

    call set_case_name("create_view_of_buffer_with_datatype")

    ds = datastore_new()
    root = ds%get_root()

    ! use create + alloc convenience methods
    ! this one is the DataType & method
    base =  root%create_view_and_allocate("base", SIDRE_INT_ID, 10)
    call base%get_data(base_vals)

    base_vals(1:5) = 10
    base_vals(6:10) = 20

    base_buff = base%get_buffer()
    ! create two views into this buffer
    ! view for the first 5 values
!--    root%createView("sub_a", base_buff, DataType::c_int(5))
    ! view for the second 5 values

!--    int * sub_a_vals = root%get_view("sub_a")%get_data(sub_a_vals)

!--    do i = 1, 5
!--       call assert_equals(sub_a_vals(i), 10)
!--    enddo

    call ds%delete()
  end subroutine create_view_of_buffer_with_datatype

  !------------------------------------------------------------------------------
  subroutine save_restore_empty_datastore
    character(24) :: file_path_base = "F_sidre_empty_datastore_"
    character(80) file_path
    integer i
    type(SidreDataStore) ds1, ds2
    type(SidreGroup) root1, root2

    call set_case_name("save_restore_empty_datastore")

    ds1 = datastore_new()
    root1 = ds1%get_root()

    do i = 1, nprotocols
       file_path = file_path_base //  protocols(i)
       call root1%save(file_path, protocols(i))
    enddo

    call ds1%delete()

    do i = 1, nprotocols
       ! Only restore sidre_hdf5 protocol
       if(protocols(i) .ne. "sidre_hdf5") then 
          cycle
       endif

       file_path = file_path_base //  protocols(i)

       ds2 = datastore_new()
       root2 = ds2%get_root()

       call root2%load(file_path, protocols(i))

       call assert_true(ds2%get_num_buffers() == 0 )
       call assert_true(root2%get_num_groups() == 0 )
       call assert_true(root2%get_num_views() == 0 )

       call ds2%delete()
    enddo
  end subroutine save_restore_empty_datastore

  !------------------------------------------------------------------------------

  subroutine save_restore_scalars_and_strings
    character(33) :: file_path_base = "F_sidre_save_scalars_and_strings_"
    character(80) file_path
    integer i
    type(SidreDataStore) ds1, ds2
    type(SidreGroup) root1, root2
    type(SidreView) view

    integer(C_INT) :: i0
    real(C_FLOAT) :: f0
    real(C_DOUBLE) :: d0
    character(80) :: s0

    call set_case_name("save_restore_scalars_and_strings")

    ds1 = datastore_new()
    root1 = ds1%get_root()

    view = root1%create_view_scalar_int("i0", 1)
    view = root1%create_view_scalar_float("f0", 1.e0)
    view = root1%create_view_scalar_double("d0", 10.d0)
    view = root1%create_view_string("s0", "I am a string")

    do i = 1, nprotocols
       file_path = file_path_base // protocols(i)
       call root1%save(file_path, protocols(i))
    enddo

    ! Restore datastore
    do i = 1, nprotocols
       ! Only restore sidre_hdf5 protocol
       if(protocols(i) .ne. "sidre_hdf5") then 
          cycle
       endif
       
       file_path = file_path_base //  protocols(i)

       ds2 = datastore_new()
       root2 = ds2%get_root()

       call root2%load(file_path, protocols(i))

       call assert_true( root1%is_equivalent_to( root2 ))

       view = root2%get_view("i0")
       i0 = view%get_data_int()
       call assert_equals(i0, 1)

       view = root2%get_view("f0")
       f0 = view%get_data_float()
       call assert_equals(f0, 1.0)

       view = root2%get_view("d0")
       d0 = view%get_data_double()
       call assert_equals(d0, 10.0d0)

       view = root2%get_view("s0")
       call view%get_string(s0)
       call assert_equals(s0, "I am a string")

       call ds2%delete()
    enddo

    call ds1%delete()
  end subroutine save_restore_scalars_and_strings

  !------------------------------------------------------------------------------
  subroutine save_restore_external_data
    character(22) :: file_path_base = "F_sidre_save_external_"
    character(80) file_path
    integer i, j
    integer, parameter :: nfoo = 10
    integer(C_INT), target :: foo1(nfoo), foo2(nfoo)
!    integer, pointer :: foo3(:) => null()
    integer(C_INT), target :: foo4(nfoo)
    integer int2d1(nfoo,2), int2d2(nfoo,2)
    integer, pointer :: int2d3(:,:)
    integer rank
    integer(SIDRE_IndexType) extents(7)
    type(SidreDataStore) ds1, ds2
    type(SidreGroup) root1, root2
    type(SidreView) view1, view2, view3, view4

    call set_case_name("save_restore_external_data")

    do i = 1, nfoo
       foo1(i) = i - 1   ! -1 to match C++
       foo2(i) = 0
       foo4(i) = i
       int2d1(i, 1) = i - 1
       int2d1(i, 2) = (i - 1) + nfoo
    enddo
    int2d2(:,:) = 0

    ds1 = datastore_new()
    root1 = ds1%get_root()

    view1 = root1%create_array_view("external_array", foo1)
    ! unable to use create_array_view since foo3 is not associated.
    view2 = root1%create_view("empty_array")
    call view2%apply_type_nelems(SIDRE_INT_ID, 0_SIDRE_IndexType)
    view3 = root1%create_view("external_undescribed")
    call view3%set_external_data_ptr(C_LOC(foo4))
    view4 = root1%create_array_view("int2d", int2d1)

    do i = 1, nprotocols
       file_path = file_path_base //  protocols(i)
       call root1%save(file_path, protocols(i))
    enddo

    call ds1%delete()
    
    ! Now load back in; 
    do i = 1, nprotocols
       ! Only restore sidre_hdf5 protocol
       if(protocols(i) .ne. "sidre_hdf5") then 
          cycle
       endif

       file_path = file_path_base //  protocols(i)

       ds2 = datastore_new()
       root2 = ds2%get_root()

       call root2%load(file_path, protocols(i))

       ! load has set the type and size of the view.
       ! now set the external address before calling load_external.
       view1 = root2%get_view("external_array")
       call assert_true(view1%is_external(), "external_array is_external")
       call assert_true(view1%is_described(), "external_array is_described")
       call assert_true(view1%get_type_id() == SIDRE_INT_ID, "external_array get_type_id")
       call assert_true(view1%get_num_elements() == nfoo, "external_array get_num_elements")
       call view1%set_array_data_ptr(foo2)

       view2 = root2%get_view("empty_array")
       call assert_true(view2%is_empty(), "empty_array is_empty")
       call assert_true(view2%is_described(), "empty_array is_described")
       call assert_true(view2%get_type_id() == SIDRE_INT_ID, "empty_array get_type_id")
       call assert_true(view2%get_num_elements() == 0, "empty_array get_num_elements")
!       call view2%set_array_data_ptr(foo3)
!       call root2%set_array_data_ptr("empty_array", foo3)

       view3 = root2%get_view("external_undescribed")
       call assert_true(view3%is_empty(), "external_undescribed is_empty")
       call assert_false(view3%is_described(), "external_undescribed is_described")

       view4 = root2%get_view("int2d")
       call assert_true(view4%is_external(), "int2d is_external")
       call assert_true(view4%is_described(), "int2d is_described")
       call assert_true(view4%get_type_id() == SIDRE_INT_ID, "int2d get_type_id")
       call assert_true(view4%get_num_elements() == nfoo*2, "int2d get_num_elements")
       call assert_equals(view4%get_num_dimensions(), 2, "int2d get_num_dimensions")
       rank = view4%get_shape(7, extents)
       call assert_equals(rank, 2, "int2d rank")
       call assert_true(extents(1) == nfoo .and. extents(2) == 2, "int2d extents")
       call view4%set_array_data_ptr(int2d2)

       ! read external data into views
       call root2%load_external_data(file_path)

       call assert_true( all(foo1 == foo2), "compare foo1 foo2" )

       ! make sure shape translates into fortran
       call view4%get_data(int2d3)
       call assert_equals(size(int2d3,1), nfoo, "int2d fortran size(1)")
       call assert_equals(size(int2d3,2), 2, "int2d fortran size(2)")
       call assert_true( all(int2d1 == int2d3), "compare ind2d1==ind2d3" )

       call ds2%delete()
    enddo
   
   end subroutine save_restore_external_data

  !------------------------------------------------------------------------------

  subroutine save_restore_other
    character(24) :: file_path_base = "F_sidre_empty_other_"
    character(80) file_path
    integer i
    integer, parameter :: ndata = 10
    integer(SIDRE_IndexType) :: shape1(2), shape2(7)
    integer rank
    type(SidreDataStore) ds1, ds2
    type(SidreGroup) root1, root2
    type(SidreView) view1, view2, view3, view4

    call set_case_name("save_restore_other")

    ds1 = datastore_new()
    root1 = ds1%get_root()

    shape1 = [ndata, 2]
    view1 = root1%create_view("empty_view")
    view2 = root1%create_view("empty_described", SIDRE_INT_ID, ndata)
    view3 = root1%create_view_with_shape("empty_shape", SIDRE_INT_ID, 2, shape1)

    view4 = root1%create_view_with_shape_and_allocate("buffer_shape", &
         SIDRE_INT_ID, 2, shape1)

    do i = 1, nprotocols
       file_path = file_path_base //  protocols(i)
       call root1%save(file_path, protocols(i))
    enddo

    call ds1%delete()

    ! Now load back in; 
    do i = 1, nprotocols
       ! Only restore sidre_hdf5 protocol
       if(protocols(i) .ne. "sidre_hdf5") then 
          cycle
       endif

       file_path = file_path_base //  protocols(i)

       ds2 = datastore_new()
       root2 = ds2%get_root()

       call root2%load(file_path, protocols(i))

       view1 = root2%get_view("empty_view")
       call assert_true(view1%is_empty(), "empty_view is_empty")
       call assert_false(view1%is_described(), "empty_view is_described")

       view2 = root2%get_view("empty_described")
       call assert_true(view2%is_empty(), "empty_described is_empty")
       call assert_true(view2%is_described(), "empty_described is_described")
       call assert_true(view2%get_type_id() == SIDRE_INT_ID, "empty_described get_type_id")
       call assert_true(view2%get_num_elements() == ndata, "empty_described get_num_elements")

       view3 = root2%get_view("empty_shape")
       call assert_true(view3%is_empty(), "empty_shape is_empty")
       call assert_true(view3%is_described(), "empty_shape is_described")
       call assert_true(view3%get_type_id() == SIDRE_INT_ID, "empty_shape get_type_id")
       call assert_true(view3%get_num_elements() == ndata*2, "empty_shape get_num_elements")
       shape2(:) = 0
       rank = view3%get_shape(7, shape2)
       call assert_equals(rank, 2, "empty_shape rank")
       call assert_true(shape2(1) == ndata .and. shape2(2) == 2, "empty_shape get_shape")

       view4 = root2%get_view("buffer_shape")
       call assert_true(view4%has_buffer(), "buffer_shape has_buffer")
       call assert_true(view4%is_described(), "buffer_shape is_described")
       call assert_true(view4%get_type_id() == SIDRE_INT_ID, "buffer_shape get_type_id")
       call assert_true(view4%get_num_elements() == ndata*2, "buffer_shape get_num_elements")
       shape2(:) = 0
       rank = view4%get_shape(7, shape2)
       call assert_equals(rank, 2, "buffer_shape rank")
       call assert_true(shape2(1) == ndata .and. shape2(2) == 2, "buffer_shape get_shape")

       call ds2%delete()
    enddo
  end subroutine save_restore_other

  !------------------------------------------------------------------------------
  subroutine rename_group
    type(SidreDataStore) ds
    type(SidreGroup) root, child1, child2, child3

    call set_case_name("rename_group")

    ds = datastore_new()
    root = ds%get_root()
    child1 = root%create_group("g_a")
    child2 = root%create_group("g_b")
    child3 = root%create_group("g_c")

    call assert_true(child1%rename("g_r"))
    call assert_equals(child1%get_name(), "g_r")
    call assert_true(root%has_group("g_r"))
    call assert_false(root%has_group("g_a"))

    call assert_false(child2%rename("fields/g_s"))
    call assert_equals(child2%get_name(), "g_b")

    call assert_false(child3%rename("g_b"))
    call assert_equals(child3%get_name(), "g_c")
  end subroutine rename_group

  !------------------------------------------------------------------------------
  subroutine save_restore_simple
    type(SidreDataStore) ds, ds2
    type(SidreGroup) root, root2, flds, ga
    type(SidreView) i0_view

    call set_case_name("save_restore_simple")

    ds = datastore_new()
    root = ds%get_root()
    flds = root%create_group("fields")

    ga = flds%create_group("a")

    i0_view = ga%create_view_scalar_int("i0", 1)

    call assert_true(root%has_group("fields"))
    call assert_true(flds%has_group("a"))
    call assert_true(ga%has_view("i0"))

!   TODO - fix wrapping, change to datastore save call, pass in group pointer.
!    call root%save("F_out_sidre_group_save_restore_simple","conduit")

    call ds%print()

!    ds2 = datastore_new()
!    root2 = ds2%get_root()

!    call root2%load("F_out_sidre_group_save_restore_simple","conduit")

!    call ds2%print()

!    flds = root2%get_group("fields")
    ! check that all sub groups exist
!    call assert_true(flds%has_group("a"))
!    ga = flds%get_group("a")
!    i0_view = ga%get_view("i0")
!    call assert_equals(i0_view%get_data_int(), 1)

!    call ds2%print()
    
!    call ds%delete()
!    call ds2%delete()
  end subroutine save_restore_simple

  !------------------------------------------------------------------------------
  subroutine save_restore_complex
    type(SidreDataStore) ds, ds2
    type(SidreGroup) root, flds, root2
    type(SidreGroup) ga, gb, gc
    type(SidreView) i0_view, f0_view, d0_view
    
    call set_case_name("save_restore_complex")

    ds = datastore_new()
    root = ds%get_root()
    flds = root%create_group("fields")

    ga = flds%create_group("a")
    gb = flds%create_group("b")
    gc = flds%create_group("c")

    i0_view = ga%create_view_scalar_int("i0", 1)
    f0_view = gb%create_view_scalar_float("f0", 100.0)
    d0_view = gc%create_view_scalar_double("d0", 3000.0d0)

    ! check that all sub groups exist
    call assert_true(flds%has_group("a"))
    call assert_true(flds%has_group("b"))
    call assert_true(flds%has_group("c"))

    call ds%print()
! TODO - fix wrapping change to datastore save call ( pass in group ptr )
    !call root%save("F_out_sidre_group_save_restore_complex","conduit")

!    ds2 = datastore_new()
!    root2 = ds2%get_root()

!    call root2%load("F_out_sidre_group_save_restore_complex","conduit")

!    flds = root2%get_group("fields")
    ! check that all sub groups exist
!    call assert_true(flds%has_group("a"))
!    call assert_true(flds%has_group("b"))
!    call assert_true(flds%has_group("c"))
    
!    ga = flds%get_group("a")
!    gb = flds%get_group("b")
!    gc = flds%get_group("c")

!    i0_view = ga%get_view("i0")
!    f0_view = gb%get_view("f0")
!    d0_view = gc%get_view("d0")

!    call assert_equals(i0_view%get_data_int(), 1)
!    call assert_equals(f0_view%get_data_float(), 100.0)
!    call assert_equals(d0_view%get_data_double(), 3000.0d0)

!    call ds2%print()

!    call ds%delete()
!    call ds2%delete()
  end subroutine save_restore_complex

  !------------------------------------------------------------------------------
  subroutine save_load_preserve_contents
    character(38) :: file_path_base0 = "F_sidre_save_preserve_contents_tree0_"
    character(38) :: file_path_base1 = "F_sidre_save_preserve_contents_tree1_"
    character(80) :: file_path0, file_path1

    type(SidreDataStore) ds, ds2, dsload
    type(SidreGroup) root, tree0, tree1, ldroot, ldtree0
    type(SidreGroup) ga, gb, gc, gx, gy, gz
    type(SidreView) i0_view, f0_view, s0_view, i10_view
    type(SidreView) i20_view, i1_view, f1_view
    character(80) :: s0
    integer(C_INT), pointer :: v1_vals(:), v2_vals(:)
    integer(C_INT) i, j, i0, i1
    real(C_FLOAT) :: f0, f1

    call set_case_name("save_load_preserve_contents")

    ds = datastore_new()
    root = ds%get_root()
    tree0 = root%create_group("tree0")

    ga = tree0%create_group("a")
    gb = tree0%create_group("b")
    gc = tree0%create_group("c")

    i0_view = ga%create_view_scalar_int("i0", 100)
    f0_view = ga%create_view_scalar_float("f0", 3000.0)
    s0_view = gb%create_view_string("s0", "foo")
    i10_view = gc%create_view_and_allocate("int10", SIDRE_INT_ID, 10)

    call i10_view%get_data(v1_vals)

    do i = 1, 10
       v1_vals(i) = i
    enddo

    do j = 1, nprotocols
       file_path0 = file_path_base0 //  protocols(j)
       call tree0%save(file_path0, protocols(j) )

       tree1 = root%create_group("tree1")

       gx = tree1%create_group("x")
       gy = tree1%create_group("y")
       gz = tree1%create_group("z")

       i20_view = gx%create_view_and_allocate("int20", SIDRE_INT_ID, 20)

       call i20_view%get_data(v2_vals)

       do i = 1, 10
          v2_vals(i) = 2*i
       enddo

       i1_view = gy%create_view_scalar_int("i1", 400)
       f1_view = gz%create_view_scalar_float("f1", 17.0)

       file_path1 =  file_path_base1 // protocols(j) 
       call tree1%save(file_path1, protocols(j) )

       call ds%print()

       dsload = datastore_new()
       ldroot = dsload%get_root()

       ldtree0 = ldroot%create_group("tree0")
       call ldtree0%load(file_path0, protocols(j) )
       call ldtree0%load(file_path1, protocols(j), .true.)

       call dsload%print()

       i0_view = ldroot%get_view("tree1/a/i0")
       i0 = i0_view%get_data_int()
       call assert_equals(i0, 100)

       f0_view = ldroot%get_view("tree1/a/f0")
       f0 = f0_view%get_data_float()
       call assert_equals(f0, 3000.0)

       s0_view = ldroot%get_view("tree1/b/s0")
       call s0_view%get_string(s0)
       call assert_equals(s0, "foo")

       i1_view = ldroot%get_view("tree1/y/i1")
       i1 = i1_view%get_data_int()
       call assert_equals(i1, 400)

       f1_view = ldroot%get_view("tree1/z/f1")
       f1 = f1_view%get_data_float()
       call assert_equals(f1, 17.0)

       i10_view = ldroot%get_view("tree1/c/int10")
       i20_view = ldroot%get_view("tree1/c/int20")

       call i10_view%get_data(v1_vals)
       call i20_view%get_data(v2_vals)

       do i = 1, 10
         call assert_equals(v1_vals(i), i)
       enddo

       do i = 1, 20
         call assert_equals(v2_vals(i), 2*i)
       enddo
       
       ! Delete the group so it is ready to use by the next protocol
       call root%destroy_group("tree1")
    enddo

  end subroutine save_load_preserve_contents



!----------------------------------------------------------------------
end module sidre_group_test
!----------------------------------------------------------------------

program fortran_test
  use fruit
  use sidre_group_test
  implicit none
  logical ok

  call init_fruit

  call get_name
  call get_path_name
  call get_parent
  call get_datastore
  call get_group
  call get_view
  call get_view_name_index
  call get_first_and_next_group_index
  call get_first_and_next_view_index
  call get_group_name_index
  call create_destroy_has_view
  call create_destroy_has_group
  call group_name_collisions
  call view_copy_move
  call groups_move_copy
  call create_destroy_view_and_data
  call create_destroy_alloc_view_and_data
  call create_view_of_buffer_with_datatype
  call save_restore_empty_datastore
  call save_restore_scalars_and_strings
  call save_restore_external_data
  call save_restore_other
! TODO - redo these, the C++ tests were heavily rewritten
!  call save_restore_simple
!  call save_restore_complex

  call fruit_summary
  call fruit_finalize

  call is_all_successful(ok)
  if (.not. ok) then
     call exit(1)
  endif
end program fortran_test
