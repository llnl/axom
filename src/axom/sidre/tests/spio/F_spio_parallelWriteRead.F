! Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
! other Axom Project Developers. See the top-level LICENSE file for details.
!
! SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

program spio_parallel_write_read
  use axom_sidre
  use axom_spio
  implicit none

  include 'mpif.h'

  integer i
  integer mpierr
  integer num_files
  integer testvalue1, testvalue2
  integer num_elems1, num_elems2
  integer return_val
  integer my_rank, num_ranks, num_output
  integer, pointer :: i1_vals(:), i2_vals(:)

  type(SidreDataStore) ds1, ds2
  type(SidreGroup) root1, root2
  type(SidreGroup) flds, flds2
  type(SidreGroup) ga, gb
  type(SidreView)  view1, view2

  type(SidreDataStore) dsextra
  type(SidreGroup) extra_root, extra, child
  type(SidreView) view

  type(IOManager) writer, reader

#ifdef AXOM_USE_HDF5
  character(18) :: protocol = "sidre_hdf5        "
  character(18) :: root_ext = ".root             "
#else
  character(18) :: protocol = "sidre_conduit_json"
  character(18) :: root_ext = ".root"
#endif  

  call mpi_init(mpierr)

  call mpi_comm_rank(MPI_COMM_WORLD, my_rank, mpierr)
  call mpi_comm_size(MPI_COMM_WORLD, num_ranks, mpierr)

#ifdef AXOM_USE_HDF5
  num_output = num_ranks / 2 
  if (num_output == 0) then
     num_output = 1
  endif
#else
  num_output = num_ranks 
#endif

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
  
  view1 = ga%create_view_scalar_int("i0", 101*my_rank)
  view2 = gb%create_view_and_allocate("i1", SIDRE_INT_ID, 10)

  call view2%get_data(i1_vals)

  do i = 0, 9
     i1_vals(i+1) = (i+10) * (404-my_rank-i)
  enddo

  ! contents of the datastore written to files with io_manager.
  num_files = num_output
  writer = IOManager(MPI_COMM_WORLD)

  call writer%write(root1, num_files, "F_out_spio_parallel_write_read", protocol )

  ! Extra stuff to exercise writeGroupToRootFile for 'sidre_hdf5' protocol
  call mpi_barrier(MPI_COMM_WORLD, mpierr)
  if (my_rank == 0 .and. protocol .eq. "sidre_hdf5") then
     dsextra = SidreDataStore()
     extra_root = dsextra%get_root()
     extra = extra_root%create_group("extra")
     view = extra%create_View_Scalar("dval", 1.1d0)
     child = extra%create_group("child")
     view = child%create_view_scalar("ival", 7)
     view = child%create_view_string("word0", "hello")
     view = child%create_view_string("word1", "world")

     call writer%write_group_to_root_file(extra, "F_out_spio_parallel_write_read" // root_ext )

     call dsextra%delete()
  endif
  call mpi_barrier(MPI_COMM_WORLD, mpierr)

  ! create another datastore that holds nothing but the root group.
  ds2 = SidreDataStore()

  ! read from the files that were written above.
  reader = IOManager(MPI_COMM_WORLD)

  root2 = ds2%get_root()
  call reader%read(root2, "F_out_spio_parallel_write_read" // root_ext )

  ! verify that the contents of ds2 match those written from ds.
  return_val = 0
  if (.not. root2%is_equivalent_to(root1)) then
     return_val = 1 
  endif

  view1 = root1%get_view("fields/a/i0")
  view2 = root2%get_view("fields/a/i0")
  testvalue1 = view1%get_data_int()
  testvalue2 = view1%get_data_int()

  if (testvalue1 .ne. testvalue2) then
     return_val = 1
  endif

  view1 = root1%get_view("fields2/b/i1")
  view2 = root2%get_view("fields2/b/i1")

  num_elems1 = view1%get_num_elements()
  num_elems2 = view2%get_num_elements()
  if (num_elems1 .ne. num_elems2) then
     return_val = 1
  else
     call view1%get_data(i1_vals)
     call view2%get_data(i2_vals)
     if (any(i1_vals.ne.i2_vals)) then
        return_val = 1
     endif
  endif

  call ds1%delete()
  call ds2%delete()

  call mpi_finalize(mpierr)

  if (return_val .ne. 0) then
     call exit(1)
  endif
end program spio_parallel_write_read

