# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Slam version of lulesh 2.0
#------------------------------------------------------------------------------


set(lulesh_headers
    lulesh.hpp
    # lulesh_tuple.hpp
    )

set(lulesh_sources
    lulesh.cpp
    lulesh-comm.cpp
    lulesh-init.cpp
    lulesh-util.cpp
    lulesh-viz.cpp
    )


# Ignore omp pragmas when openMP is not defined
if ( NOT AXOM_ENABLE_OPENMP )
    set_source_files_properties(${lulesh_sources}
        PROPERTIES COMPILE_FLAGS "${AXOM_DISABLE_OMP_PRAGMA_WARNINGS}" )
endif()

axom_add_executable(
    NAME        slam_lulesh_ex
    SOURCES     ${lulesh_sources}
    HEADERS     ${lulesh_headers}
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  slam
    FOLDER      axom/slam/examples )

if(AXOM_ENABLE_TESTS)
    if(AXOM_ENABLE_MPI)

        # CMake inserts Fortran flag that forces one OpenMP thread when
        # compiling with MPI and xlf
        if(ENABLE_FORTRAN)
          list(REMOVE_ITEM CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "xlomp_ser")
        endif()

        axom_add_test(NAME        slam_lulesh
                      COMMAND     slam_lulesh_ex
                      NUM_MPI_TASKS   8
                      NUM_OMP_THREADS 1 )
    else()
        axom_add_test(NAME        slam_lulesh
                      COMMAND     slam_lulesh_ex
                      NUM_OMP_THREADS 4 )
    endif()
endif()
