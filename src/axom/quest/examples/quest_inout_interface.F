! Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
! other Axom Project Developers. See the top-level LICENSE file for details.
!
! SPDX-License-Identifier: (BSD-3-Clause)

!-------------------------------------------------------------------------------
! Simple Fortran example using the quest_inout Fortran interface
!
! Usage: Pass in a single command line argument with the path to a surface mesh
!-------------------------------------------------------------------------------

module quest_inout_interface_ex

#include "axom/config.hpp"

 use iso_c_binding, only : C_DOUBLE

#if defined AXOM_USE_MPI && ! defined AXOM_USE_MPIF_HEADER
  use mpi
#endif
  use axom_quest
  implicit none

#if defined AXOM_USE_MPI && defined AXOM_USE_MPIF_HEADER
  include 'mpif.h'
#endif

  integer :: num_pts = 10, rank, nprocs, ierr
  character(LEN=:), allocatable :: FileName
  integer :: dim = 3

contains

  !-------------------------------------------------------------------------------
  ! random_double( lb, ub )
  !
  ! Returns a random double, d, between lowerbound lb and upperbound ub
  !-------------------------------------------------------------------------------
  function random_double(lb,ub ) result(d)
    implicit none

    real(C_DOUBLE), value, intent(in) ::  lb, ub
    real(C_DOUBLE) :: d

    call random_number(d)
    d =  (ub-lb) * d + lb

  end function random_double

  subroutine init_mpi()
    !...initialize MPI
#ifdef AXOM_USE_MPI
    call mpi_init( ierr )
    call mpi_comm_rank( MPI_COMM_WORLD, rank, ierr )
    call mpi_comm_size( MPI_COMM_WORLD, nprocs, ierr )
#endif
  end subroutine init_mpi

  subroutine finalize_mpi()
    !... finalize MPI
#ifdef AXOM_USE_MPI
    call mpi_finalize( ierr )
#endif
  end subroutine finalize_mpi

  !-------------------------------------------------------------------------------
  ! run_inout_queries( rank )
  !
  ! Exercises quest's inout containment query functionality on a given mpi rank
  !-------------------------------------------------------------------------------
  subroutine run_inout_queries()
    implicit none

    integer :: i, rc
    logical :: ins, is_verbose = .true.
    real(C_DOUBLE)  :: x, y, z
    real(C_DOUBLE)  :: EPS = 1E-7
    real(C_DOUBLE) :: mesh_bb_min(3), mesh_bb_max(3), mesh_center_of_mass(3)
    character(LEN=60) :: Format = "(I4, I4, F8.2, F8.2, F8.2, L3)"

    !...set parameters for the inout query
    rc = quest_inout_set_verbose(is_verbose)
    rc = quest_inout_set_vertex_weld_threshold(EPS)

    if (rank == 0) then
      write(*,*) "** Initializing quest inout query"
    end if

    !...initialize inout query
#ifdef AXOM_USE_MPI
    rc = quest_inout_init( FileName, MPI_COMM_WORLD )
#else
    rc = quest_inout_init( FileName )
#endif

    !...apply some basic error checking
    if ( rc /= 0 ) then 
      error stop "quest_inout_init() failed!"
    endif

    if (.not. quest_inout_initialized() ) then 
      error stop "quest_inout not initialized after call quest_inout_init()"
    endif

    !...access properties of the query
    rc = quest_inout_mesh_min_bounds(mesh_bb_min)
    rc = quest_inout_mesh_max_bounds(mesh_bb_max)
    rc = quest_inout_mesh_center_of_mass(mesh_center_of_mass)

    if (rank == 0) then
      write(*,*) "** Running containment queries using quest fortran interface"
      write(*,*) "[rank, id, (  x,  y,  z),   inside?]"
    end if

    !...run the inout query at some fixed and random points
    do i=1,num_pts + 1
      if (i == 1) then
        x = mesh_center_of_mass(1)
        y = mesh_center_of_mass(2)
        z = mesh_center_of_mass(3)
      else
        x = random_double( mesh_bb_min(1), mesh_bb_max(1) )
        y = random_double( mesh_bb_min(2), mesh_bb_max(2) )
        z = random_double( mesh_bb_min(3), mesh_bb_max(3) )
      end if

      ins = quest_inout_evaluate(x,y,z)

      write(*,Format) rank, i-1, x, y, z, ins
    end do

    !...finalize quest_inout query
    rc = quest_inout_finalize()
  end subroutine run_inout_queries

end module quest_inout_interface_ex


!-------------------------------------------------------------------------------
! Main Program
!-------------------------------------------------------------------------------
program quest_inout_interface_fortran_ex

  use quest_inout_interface_ex

  implicit none
  integer len_name, arg_count

  call init_mpi()

  !...parse command line arguments
  arg_count = command_argument_count()
  if (arg_count /= 1) then 
    error stop "Usage: [mpiriun -np N] ./quest_inout_interface_F_ex <stl_file>"
  end if

  call get_command_argument(1, length=len_name)
  allocate(character(len_name):: FileName)
  call get_command_argument(1, FileName)

  !...exercise the inout query functions in quest module
  call run_inout_queries()

  deallocate(FileName)
  call finalize_mpi()

end program
