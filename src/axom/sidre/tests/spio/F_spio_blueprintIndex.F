! Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
! other Axom Project Developers. See the top-level LICENSE file for details.
!
! SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

program spio_blueprint_index
  use axom_sidre
  use axom_spio
  use iso_c_binding
  implicit none

  include 'mpif.h'

  integer i
  integer mpierr
  integer return_val
  integer my_rank, num_ranks
  integer node_count, elt_count
  integer, pointer :: node_vals(:), conn_vals(:)
  real(C_DOUBLE), pointer :: x_vals(:), y_vals(:), z_vals(:), elt_vals(:)


  type(SidreDataStore) ds1, cds
  type(SidreGroup) cdroot, mroot
  type(SidreGroup) coords, topos, fields
  type(SidreGroup) connmesh, coordval, elts
  type(SidreGroup) nodefield, eltfield
  type(SidreGroup) nodegrp, eltgrp
  type(SidreView) connectivity, coordtype
  type(SidreView) eltshape, meshtype, topocoords
  type(SidreView) xview, yview, zview
  type(SidreView) nvrtx, nsclr, nmesh, nvals
  type(SidreView) evrtx, esclr, emesh, evals

  type(IOManager) writer

  call mpi_init(mpierr)

  call mpi_comm_rank(MPI_COMM_WORLD, my_rank, mpierr)
  call mpi_comm_size(MPI_COMM_WORLD, num_ranks, mpierr)

  return_val = 0
  node_count = 12
  elt_count = 2

  cds = SidreDataStore()
  cdroot = cds%get_root()
  mroot = cdroot%create_group("domain_data/domain")
  coords = mroot%create_group("mesh/coordsets/coords")
  topos = mroot%create_group("mesh/topologies")
  fields = mroot%create_group("mesh/fields")

  coordtype = coords%create_view_string("type", "explicit")
  coordval = coords%create_group("values")
  xview = coordval%create_view_and_allocate("x", SIDRE_DOUBLE_ID, node_count)
  yview = coordval%create_view_and_allocate("y", SIDRE_DOUBLE_ID, node_count)
  zview = coordval%create_view_and_allocate("z", SIDRE_DOUBLE_ID, node_count)

  call xview%get_data(x_vals)
  call yview%get_data(y_vals)
  call zview%get_data(z_vals)

  do i = 1, node_count
    x_vals(i) = real(modulo(((i+1) / 2),2))
    y_vals(i) = real(modulo((i / 2),2))
    z_vals(i) = real(i / 4)
  enddo

  connmesh = topos%create_group("mesh")
  meshtype = connmesh%create_view_string("type", "unstructured")
  topocoords = connmesh%create_view_string("coordset", "coords");
  elts = connmesh%create_group("elements");
  eltshape = elts%create_view_string("shape", "hex");

  connectivity = elts%create_view_and_allocate("connectivity", SIDRE_INT_ID, 16)

  call connectivity%get_data(conn_vals);
  do i = 1, 7
    conn_vals(i) = i
  enddo
  do i = 8, 16
    conn_vals(i) = i - 4
  enddo

  nodegrp = fields%create_group("nodefield")
  nvrtx =  nodegrp%create_view_string("association", "vertex")
  nsclr = nodegrp%create_view_string("type", "scalar")
  nmesh = nodegrp%create_view_string("topology", "mesh")
  nvals = nodegrp%create_view_and_allocate("values", SIDRE_INT_ID, node_count)

  call nvals%get_data(node_vals)

  do i = 1, node_count
    node_vals(i) = int(x_vals(i) + y_vals(i) + z_vals(i))
  enddo

  eltgrp = fields%create_group("eltfield")
  evrtx = eltgrp%create_view_string("association", "vertex")
  esclr = eltgrp%create_view_string("type", "scalar")
  emesh = eltgrp%create_view_string("topology", "mesh")
  evals = eltgrp%create_view_and_allocate("values", SIDRE_DOUBLE_ID, elt_count)

  call evals%get_data(elt_vals)
  elt_vals(1) = 2.65;
  elt_vals(2) = 1.96;

  writer = IOManager(MPI_COMM_WORLD)

  call writer%write(mroot, 1, "fbp", "sidre_hdf5")

  call writer%write_blueprint_index_to_root_file(cds, "domain_data/domain/mesh", "fbp.root", "mesh")

  call mpi_finalize(mpierr)

  if (return_val .ne. 0) then
     call exit(1)
  endif

end program spio_blueprint_index


