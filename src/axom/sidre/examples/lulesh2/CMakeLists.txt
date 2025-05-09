# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Sidre's implementation of Lulesh
#------------------------------------------------------------------------------

set(lulesh_headers
    lulesh.h
    lulesh_tuple.h )

set(lulesh_sources
    lulesh.cc
    lulesh-comm.cc
    lulesh-init.cc
    lulesh-util.cc
    lulesh-viz.cc )

axom_add_executable(
    NAME       sidre_lulesh2_ex
    SOURCES    ${lulesh_sources}
    HEADERS    ${lulesh_headers}
    OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON sidre
    FOLDER     axom/sidre/examples )

if(AXOM_ENABLE_TESTS)
    if(AXOM_ENABLE_MPI)
        axom_add_test(NAME sidre_lulesh2
                      COMMAND sidre_lulesh2_ex -s 5
                      NUM_MPI_TASKS 1
                      NUM_OMP_THREADS 2 )
    else()
        axom_add_test(NAME sidre_lulesh2 -s 5
                      COMMAND sidre_lulesh2_ex
                      NUM_OMP_THREADS 2 )
    endif()
endif()
