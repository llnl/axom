! Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
! other Axom Project Developers. See the top-level LICENSE file for details.
!
! SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

program spio_external_write_read

  use iso_c_binding
  use axom_sidre
  use axom_spio
  implicit none

  include 'mpif.h'

  integer, parameter :: nvals = 10
  integer(C_INT), target :: orig_vals1(nvals), restored_vals1(nvals)
  integer(C_INT), target :: orig_vals2(nvals), restored_vals2(nvals)
  integer my_rank, num_ranks, mpierr
  integer num_output, num_files
  integer testvalue1, testvalue2
  integer return_val
  integer i

  type(SidreDataStore) ds1, ds2
  type(SidreGroup) root1, root2
  type(SidreGroup) flds, flds2
  type(SidreGroup) ga, gb
  type(SidreView)  view1, view2

  type(IOManager) writer, reader

! Loading external data is only currently valid for 'sidre_hdf5' protocol
#ifdef AXOM_USE_HDF5
  call mpi_init(mpierr)

  call mpi_comm_rank(MPI_COMM_WORLD, my_rank, mpierr)
  call mpi_comm_size(MPI_COMM_WORLD, num_ranks, mpierr)

  num_output = num_ranks / 2 
  if (num_output == 0) then
     num_output = 1
  endif

  do i = 0, 9
     orig_vals1(i+1) = (i+10) * (404-my_rank-i)
     orig_vals2(i+1) = (i+10) * (404-my_rank-i) + 20
  enddo

  ! create a datastore and give it a small hierarchy of groups and views.
  !
  ! the views are filled with repeatable nonsense data that will vary based
  ! on rank.

  ds1 = SidreDataStore()
  root1 = ds1%get_root()

  flds = root1%create_group("fields")
  flds2 = root1%create_group("fields2")

  ga = flds%create_group("a")
  gb = flds2%create_group("b")
  view1 = ga%create_array_view("external_array", orig_vals1)
  view1 = gb%create_view("external_undescribed")
  call view1%set_external_data_ptr(c_loc(orig_vals2))

 ! contents of the datastore written to files with io_manager.
  num_files = num_output
  writer = IOManager(MPI_COMM_WORLD)

  call writer%write(root1, num_files, "F_out_spio_external_write_read", "sidre_hdf5")

  ! create another datastore than holds nothing but the root group.
  ds2 = SidreDataStore()
  root2 = ds2%get_root()

  ! read from the files that were written above.
  reader = IOManager(MPI_COMM_WORLD)

  call reader%read(root2, "F_out_spio_external_write_read.root")

  restored_vals1(:) = -1
  restored_vals2(:) = -1
  view1 = root2%get_view("fields/a/external_array")
  call view1%set_array_data_ptr(restored_vals1)

  view2 = root2%get_view("fields2/b/external_undescribed")
  call view2%set_external_data_ptr(c_loc(restored_vals2))

  call reader%load_external_data(root2, "F_out_spio_external_write_read.root") 


  ! verify that the contents of ds2 match those written from ds1.
  return_val = 0
  if (.not. root2%is_equivalent_to(root1)) then
    return_val = 1
  else if (view1%get_num_elements() .ne. nvals) then
     return_val = 1
  else if (any(orig_vals1 .ne. restored_vals1)) then
     return_val = 1
  else if (any(-1 .ne. restored_vals2)) then
     ! external_undescribed was not written to disk (since it is undescribed)
     ! make sure it was not read in.
     return_val = 1
  endif

  call ds1%delete()
  call ds2%delete()

  call mpi_finalize(mpierr)
#else
  return_val = 0
  print *, "Loading external data only valid for 'sidre_hdf5' format &
           & and when Axom is configured with hdf5. Skipping test."
#endif

  if (return_val .ne. 0) then
     call exit(1)
  endif
end program spio_external_write_read

