! Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
! other Axom Project Developers. See the top-level LICENSE file for details.
!
! SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

program spio_irregularWriteRead
  use axom_sidre
  use axom_spio
  implicit none

  include 'mpif.h'

  integer i, f, g
  integer mpierr
  integer num_files
  integer return_val
  integer my_rank, num_ranks, num_output
  integer num_fields, num_subgroups
  integer num_elems1, num_elems2
  integer testvalue1, testvalue2
  integer, pointer :: vals1(:), vals2(:)
  character(80) name

  type(SidreDataStore) ds1, ds2
  type(SidreGroup) root1, root2
  type(SidreGroup) flds1, flds2
  type(SidreGroup) sg1, sg2
  type(SidreView)  view1, view2

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

  num_fields = my_rank + 2

  do f = 0, num_fields-1
     write(name, "(a,i0)") "fields", f
     flds1 = root1%create_group(name)

     num_subgroups = mod(f+my_rank,3) + 1
     do g = 0, num_subgroups-1
        write(name, "(a,i0)") "subgroup", g
        sg1 = flds1%create_group(name)

        write(name, "(a,i0)") "view", g
        if (mod(g, 2) .ne. 0) then
           view1 = sg1%create_view_and_allocate(name, SIDRE_INT_ID, 10+my_rank)
           call view1%get_data(vals1)
           do i = 0, 9+my_rank
              vals1(i+1) = (i+10) * (404-my_rank-i-g-f)
           enddo
        else
           view1 = sg1%create_view_scalar_int(name, 101*my_rank*(f+g+1))
        endif
     enddo
  enddo

  ! contents of the datastore written to files with io_manager.
  num_files = num_output
  writer = IOManager(MPI_COMM_WORLD)

  call writer%write(root1, num_files, "F_out_spio_irregular_write_read", protocol )

  ! create another datastore that holds nothing but the root group.
  ds2 = SidreDataStore()

  ! read from the files that were written above.
  reader = IOManager(MPI_COMM_WORLD)

  root2 = ds2%get_root()
  call reader%read(root2, "F_out_spio_irregular_write_read" // root_ext )

  ! verify that the contents of ds2 match those written from ds.
  return_val = 0
  if (.not. root2%is_equivalent_to(root1)) then
     return_val = 1
  endif

  do f = 0, num_fields-1
     write(name, "(a,i0)") "fields", f
     flds1 = root1%get_group(name)
     flds2 = root2%get_group(name)

     num_subgroups = mod(f+my_rank,3) + 1
     do g = 0, num_subgroups-1
        write(name, "(a,i0)") "subgroup", g
        sg1 = flds1%get_group(name)
        sg2 = flds2%get_group(name)

        write(name, "(a,i0)") "view", g
        if (mod(g, 2) .ne. 0) then
           view1 = sg1%get_view(name)
           view2 = sg2%get_view(name)

           num_elems1 = view1%get_num_elements()
           num_elems2 = view2%get_num_elements()
           if (num_elems1 .ne. num_elems2) then
              return_val = 1 
           else
              call view1%get_data(vals1)
              call view2%get_data(vals2)
              if (any(vals1.ne.vals2)) then
                 return_val = 1
              endif
           endif
        else
           view1 = sg1%get_view(name)
           view2 = sg2%get_view(name)
           testvalue1 = view1%get_data_int()
           testvalue2 = view1%get_data_int()
           if (testvalue1 .ne. testvalue2) then
              return_val = 1
           endif
        endif
     enddo
  enddo

  call ds1%delete()
  call ds2%delete()

  call mpi_finalize(mpierr)

  if (return_val .ne. 0) then
     call exit(1)
  endif
end program spio_irregularWriteRead

