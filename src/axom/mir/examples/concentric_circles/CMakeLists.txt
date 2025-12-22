# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

set( mir_example_dependencies
     core
     slic
     mir
     bump
   )

# Speed up compilation by separating out the execspaces into different runMIR_xxx.cpp files.
axom_add_executable(
    NAME       mir_concentric_circles
    SOURCES    mir_concentric_circles.cpp MIRApplication.cpp runMIR_seq.cpp runMIR_omp.cpp runMIR_cuda.cpp runMIR_hip.cpp
    OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON ${mir_example_dependencies}
    FOLDER     axom/mir/examples
    )

if(AXOM_ENABLE_MPI AND CONDUIT_RELAY_MPI_ENABLED)
  # Add a parallel version of the executable.
  axom_add_executable(
      NAME       mir_concentric_circles_mpi
      SOURCES    mir_concentric_circles_mpi.cpp MIRApplication.cpp runMIR_seq.cpp runMIR_omp.cpp runMIR_cuda.cpp runMIR_hip.cpp
      OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
      DEPENDS_ON ${mir_example_dependencies} conduit_relay_mpi_io
      FOLDER     axom/mir/examples
    )
endif()

if(AXOM_ENABLE_TESTS)
    set (_policies "seq")
    if(RAJA_FOUND AND UMPIRE_FOUND)
        blt_list_append(TO _policies ELEMENTS "omp" IF AXOM_ENABLE_OPENMP)
        blt_list_append(TO _policies ELEMENTS "cuda" IF AXOM_ENABLE_CUDA)
        blt_list_append(TO _policies ELEMENTS "hip" IF AXOM_ENABLE_HIP)
    endif()

    foreach(_policy ${_policies})
        # sets the number of threads for the omp policy; leave empty for non-omp policies
        set(_num_threads)
        if(_policy STREQUAL "omp")
            set(_num_threads ${AXOM_TEST_NUM_OMP_THREADS})
        endif()
            
        set(_testname "mir_concentric_circles_${_policy}")

        axom_add_test(
            NAME ${_testname}
            COMMAND mir_concentric_circles
                --gridsize 100
                --numcircles 5
                --policy ${_policy}
                --disable-write
            NUM_OMP_THREADS ${_num_threads})

        set_tests_properties(${_testname} PROPERTIES
            PASS_REGULAR_EXPRESSION  "Material interface reconstruction time:")
    endforeach()
endif()
