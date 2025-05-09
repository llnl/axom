! Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
! other Axom Project Developers. See the top-level LICENSE file for details.
!
! SPDX-License-Identifier: (BSD-3-Clause)

!------------------------------------------------------------------------------
! Some simple types and functions used in tests
! (included in namespace to prevent clashes)
!------------------------------------------------------------------------------

module sidre_opaque_util_test
  use iso_c_binding

  ! Centering
  integer, parameter :: &
       Zone_Centering = 1, &
       Node_Centering = 2, &
       Unknown_Centering = 3

  ! DType
  integer, parameter :: &
       Double_Type = 1, &
       Int_Type = 2, &
       Unknown_Type = 3

  type Extent
     integer(C_INT) m_ilo
     integer(C_INT) m_ihi
   contains
     procedure :: getNumPts => extent_getNumPts
  end type Extent

  type MeshVar
     integer(C_INT) m_cent  ! Centering
     integer(C_INT) m_type  ! DType
     integer(C_INT) m_depth
   contains
     procedure :: getNumVals => MeshVar_getNumVals
  end type MeshVar

contains

!  function new_Extent(ilo, ihi)
!    integer, intent(in) :: ilo, ihi
!    type(Extent) new_extent
!    
!    allocate(new_extent)
!    new_extent%m_ilo = ilo
!    new_extent%m_ihi = ihi
!  end function new_base

  function extent_getNumPts(self, cent) result(retval)
    class(Extent) self
    integer(C_INT) cent
    integer(C_INT) retval

    select case ( cent )
    case(Zone_Centering)
       retval = (self%m_ihi - self%m_ilo + 1)
    case(Node_Centering)
       retval = (self%m_ihi - self%m_ilo + 2)
    case default
       retval = -1         ! I know magic numbers are bad. So sue me.
    end select
    return
  end function extent_getNumPts

  function MeshVar_getNumVals(self, ext) result(rv)
    class(MeshVar) self
    type(Extent) ext
    integer(C_INT) rv
    
    rv = ext%getNumPts(self%m_cent) * self%m_depth
  end function MeshVar_getNumVals
end module sidre_opaque_util_test

module sidre_opaque_test
  use iso_c_binding
  use fruit
  use axom_sidre
  implicit none

contains

!------------------------------------------------------------------------------
!
! Simple test that adds an opaque data object, retrieves it and checks if
! the retrieved object is in the expected state.
!
  subroutine basic_inout
    use sidre_opaque_util_test
    type(SidreDataStore) ds
    type(SidreGroup) root, problem_gp
    type(SidreView) ext_view
    type(SidreView) ext2_view
    integer, parameter :: ihi_val = 9
    integer(C_INT) test_ihi
    integer(C_INT) test_ihi2
    type(Extent), pointer :: ext, test_extent, ext2, test_extent2
    type(C_PTR) :: test_extent_ptr, test_extent2_ptr

    call set_case_name("basic_inout")

    ds = SidreDataStore()
    root = ds%get_root()

    problem_gp = root%create_group("problem")

    allocate(ext)
    ext = Extent(0, ihi_val)

    ext_view = problem_gp%create_view_external("ext", c_loc(ext))

    call assert_equals(ext_view%is_external(), .true.)
    call assert_equals(ext_view%is_applied(), .false.)
    call assert_equals(ext_view%is_opaque(), .true.)

    test_extent_ptr = ext_view%get_void_ptr()
    call c_f_pointer(test_extent_ptr, test_extent)

    test_ihi = test_extent%m_ihi
    call assert_equals(test_extent%m_ihi, ihi_val)

    ! Similar test with different view methods
    
    allocate(ext2)
    ext2 = Extent(0, 2 * ihi_val)

    ext2_view = problem_gp%create_view_empty("ext2")
    call ext2_view%set_external_data_ptr( c_loc(ext2) )

    call assert_equals(ext2_view%is_opaque(), .true.)

    test_extent2_ptr = ext2_view%get_void_ptr()
    call c_f_pointer(test_extent2_ptr, test_extent2)

    test_ihi2 = test_extent2%m_ihi
    call assert_equals(test_extent2%m_ihi, 2 * ihi_val)

    ! clean up...
    deallocate(ext)
    deallocate(ext2)
    call ds%delete()
  end subroutine basic_inout

  !------------------------------------------------------------------------------
  !
  ! Test that adds "MeshVars" as opaque data objects, creates views for their
  ! data on each of two domains, allocates their data (based on centering,
  ! domain size, and depth), and then checks to if the allocated data
  ! lengths match the expected values.
  !
  subroutine meshvar_test
    use sidre_opaque_util_test

    type(SidreDataStore) ds
    type(SidreGroup) root, problem_gp, meshvar_gp
    type(SidreView) zone_mv_view, node_mv_view

    type(SidreGroup) dom_gp
    type(SidreView) dom_zone_view, dom_node_view, ext_view

    type(SidreGroup) tmpgroup
    type(SidreView) tmpview

    integer num_zone_vals, num_node_vals
    integer test_num_zone_vals, test_num_node_vals
    integer idom

    integer(C_INT), parameter :: ilo_val(2) = [0, 10]
    integer(C_INT), parameter :: ihi_val(2) = [9, 21]
    character(7), parameter :: dom_name(2) = [ "domain0", "domain1" ]
    integer(C_INT), parameter :: zone_var_depth = 1
    integer(C_INT), parameter :: node_var_depth = 2

    type(MeshVar), allocatable, target :: zone_mv, node_mv
    type(Extent), pointer :: dom_ext
    type(Extent), allocatable, target :: domains(:)

    type(MeshVar), pointer :: zonemv, nodemv
    type(C_PTR) zonemv_ptr, nodemv_ptr, dom_ext_ptr

    call set_case_name("meshvar_test")

    ds   = SidreDataStore()
    root = ds%get_root()

    problem_gp = root%create_group("problem")

    ! Add two different mesh vars to mesh var group
    meshvar_gp = problem_gp%create_group("mesh_var")
    allocate(zone_mv)
    zone_mv = MeshVar(Zone_Centering, Int_Type, zone_var_depth)
    zone_mv_view = meshvar_gp%create_view_external("zone_mv", c_loc(zone_mv))
    allocate(node_mv)
    node_mv = MeshVar(Node_Centering, Double_Type, node_var_depth)
    node_mv_view = meshvar_gp%create_view_external("node_mv", c_loc(node_mv))

    !
    ! Create domain groups, add extents
    ! Create data views for mesh var data on domains and allocate
    !
    allocate(domains(2))
    do idom = 1, 2
       dom_gp = problem_gp%create_group(dom_name(idom))
       dom_ext => domains(idom)
       dom_ext = Extent(ilo_val(idom), ihi_val(idom))
       ext_view = dom_gp%create_view_external("ext", c_loc(dom_ext))

       zonemv_ptr = zone_mv_view%get_void_ptr()
       call c_f_pointer(zonemv_ptr, zonemv)

       dom_zone_view = dom_gp%create_view_empty("zone_data")
       call dom_zone_view%allocate(SIDRE_INT_ID, zone_mv%getNumVals(dom_ext))

       nodemv_ptr = node_mv_view%get_void_ptr()
       call c_f_pointer(nodemv_ptr, nodemv)

       dom_node_view = dom_gp%create_view_empty("node_data")
       call dom_node_view%allocate(SIDRE_DOUBLE_ID, nodemv%getNumVals(dom_ext))
    enddo

!
!  Print datastore contents to see what's going on.
!
!  ds%print()


    !
    ! Check data lengths
    !
    do idom = 1, 2
       dom_gp = problem_gp%get_group(dom_name(idom))

       tmpview = dom_gp%get_view("ext")
       dom_ext_ptr = tmpview%get_void_ptr()
       call c_f_pointer(dom_ext_ptr, dom_ext)

       zonemv_ptr = zone_mv_view%get_void_ptr()
       call c_f_pointer(zonemv_ptr, zonemv)

       nodemv_ptr = node_mv_view%get_void_ptr()
       call c_f_pointer(nodemv_ptr, nodemv)

       num_zone_vals = zonemv%getNumVals(dom_ext)
       tmpview = dom_gp%get_view("zone_data")
       test_num_zone_vals = tmpview%get_num_elements()
       call assert_equals(num_zone_vals, test_num_zone_vals)

       num_node_vals = nodemv%getNumVals(dom_ext)
       tmpview = dom_gp%get_view("node_data")
       test_num_node_vals = tmpview%get_num_elements()
       call assert_equals(num_node_vals, test_num_node_vals)
    enddo

    ! clean up...
    deallocate(zone_mv)
    deallocate(node_mv)
    deallocate(domains)
    call ds%delete()
  end subroutine meshvar_test

!----------------------------------------------------------------------
end module sidre_opaque_test
!----------------------------------------------------------------------

program fortran_test
  use fruit
  use sidre_opaque_test
  implicit none
  logical ok

  call init_fruit

  call basic_inout
  call meshvar_test

  call fruit_summary
  call fruit_finalize

  call is_all_successful(ok)
  if (.not. ok) then
     call exit(1)
  endif
end program fortran_test

