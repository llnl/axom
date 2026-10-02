# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

set( mir_example_dependencies
     core
     slic
     mir
     bump
   )

set(mir_concentric_circles_common_sources MIRApplication.cpp runMIR_seq.cpp)
blt_list_append(TO mir_concentric_circles_common_sources ELEMENTS runMIR_omp.cpp  IF AXOM_ENABLE_OPENMP)
blt_list_append(TO mir_concentric_circles_common_sources ELEMENTS runMIR_cuda.cpp IF AXOM_ENABLE_CUDA)
blt_list_append(TO mir_concentric_circles_common_sources ELEMENTS runMIR_hip.cpp  IF AXOM_ENABLE_HIP)

# Speed up compilation by separating out the execspaces into different runMIR_xxx.cpp files.
axom_add_executable(
    NAME       mir_concentric_circles
    SOURCES    mir_concentric_circles.cpp ${mir_concentric_circles_common_sources}
    OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON ${mir_example_dependencies}
    FOLDER     axom/mir/examples
    )

if(AXOM_ENABLE_MPI AND CONDUIT_RELAY_MPI_ENABLED)
  # Add a parallel version of the executable.
  axom_add_executable(
      NAME       mir_concentric_circles_mpi
      SOURCES    mir_concentric_circles_mpi.cpp ${mir_concentric_circles_common_sources}
      OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
      DEPENDS_ON ${mir_example_dependencies} conduit_relay_mpi_io
      FOLDER     axom/mir/examples
    )
endif()

if(AXOM_ENABLE_TESTS)
    foreach(_policy ${AXOM_EXECUTION_POLICIES})
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
