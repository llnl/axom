! Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
! other Axom Project Developers. See the top-level LICENSE file for details.
!
! SPDX-License-Identifier: (BSD-3-Clause)

!...include axom/config.hpp for AXOM_USE_MPI
#include "axom/config.hpp"

!------------------------------------------------------------------------------
! Helper Methods
!------------------------------------------------------------------------------
module helper_methods

  use iso_c_binding, only : C_DOUBLE

contains

  function random_double( lo, hi ) result(d)
    implicit none      
    real(C_DOUBLE), value, intent(in) ::  lo, hi
    real(C_DOUBLE) :: d

    call random_number(d)
    d =  (hi-lo) * d + lo 
  end function random_double

end module helper_methods

!------------------------------------------------------------------------------
! Program main
!------------------------------------------------------------------------------
program quest_signed_distance_interface

  use iso_c_binding, only : C_DOUBLE      

!...include helper methods and quest
  use helper_methods
  use axom_quest 

!...include mpi stuff         
#if defined (AXOM_USE_MPI) && defined (AXOM_USE_MPIF_HEADER)
  include 'mpif.h'
#elif defined (AXOM_USE_MPI)  
  use mpi
#endif


!...variables  
  integer :: ierr, rc, i, file_name_length
  integer :: ndims = 3, N=10
  character(LEN=:), allocatable :: file_name
  logical :: is_verbose = .true.
  logical :: is_watertight = .true.
  logical :: compute_signs = .true.
  real(C_DOUBLE) :: mesh_bb_min(3), mesh_bb_max(3), pt(3), cp(3), normal(3), phi, phi_cp, eps, radius

!...initialize mpi  
#ifdef AXOM_USE_MPI
  call mpi_init( ierr )
#endif

!...get filename from the command line
  call get_command_argument(1, length=file_name_length)       
  allocate(character(file_name_length):: file_name)
  call get_command_argument(1, file_name) 

!...set various parameters for the signed distance query
  call quest_signed_distance_set_dimension( ndims )
  call quest_signed_distance_set_closed_surface( is_watertight )
  call quest_signed_distance_set_verbose( is_verbose )
  call quest_signed_distance_set_compute_signs( compute_signs )

!...initialize the signed distance query
#ifdef AXOM_USE_MPI
  rc = quest_signed_distance_init( file_name, MPI_COMM_WORLD )
#else
  rc = quest_signed_distance_init( file_name )
#endif

  if ( rc /= 0 ) then 
    stop 'signed_distance_init() failed!'
  endif

!...compute the mesh bounds
  call quest_signed_distance_get_mesh_bounds( mesh_bb_min, mesh_bb_max ) 

!...evaluate the signed distance at random points within the mesh bounding box  
  do i=1,N

     do j=1,ndims
       pt(j) = random_double( mesh_bb_min(j), mesh_bb_max(j) )
     end do ! end for all dimensions

     ! evaluate signed distance (simple)
     phi = quest_signed_distance_evaluate( pt(1), pt(2), pt(3) )

     ! evaluate signed distance with closest point and normal at that point
     phi_cp = quest_signed_distance_evaluate( pt(1), pt(2), pt(3), &
                                              cp(1), cp(2), cp(3), &
                                              normal(1), normal(2), normal(3) )
     
     radius = sqrt(cp(1)**2 + cp(2)**2 + cp(3)**2)

     ! check that the two calls gave approximately the same answer
     if (abs(phi-phi_cp) >= epsilon(phi)) then
        print*, 'Got different values when calling overloads of quest_signed_distance_evaluate! ', phi, phi_cp
     end if

     print '("phi: ",F8.3,"; closest_point: ",F8.3, F8.3, F8.3, "; normal: ",F8.3,F8.3,F8.3,"; radius @ cp: ", F8.3)',  &
            phi, cp(1), cp(2), cp(3), normal(1), normal(2), normal(3), radius

  end do ! end for all points

!...finalize the signed distance query
  call quest_signed_distance_finalize( )

!...return all dynamically allocated memory to memory manager  
  deallocate(file_name)

!...finalize mpi  
#ifdef AXOM_USE_MPI
  call mpi_finalize( ierr )
#endif

end program
