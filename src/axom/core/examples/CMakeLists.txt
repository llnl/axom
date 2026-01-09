# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Core containers, numerics, and utilities examples
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# List of single source file examples
#------------------------------------------------------------------------------
set(example_sources
    core_numerics.cpp
    core_containers.cpp 
    core_acceleration.cpp)

#------------------------------------------------------------------------------
# Add targets and tests for serial examples
#------------------------------------------------------------------------------
foreach(example_source ${example_sources})
    get_filename_component(exe_name ${example_source} NAME_WE)
    axom_add_executable(
        NAME        ${exe_name}_ex
        SOURCES     ${example_source}
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  core
        FOLDER      axom/core/examples )

    if(AXOM_ENABLE_TESTS)
        axom_add_test(
            NAME    ${exe_name}
            COMMAND ${exe_name}_ex )
    endif()
endforeach()

if(AXOM_ENABLE_CUDA)
    axom_add_executable(
          NAME        core_acceleration_cuda_on_ex
          SOURCES     core_acceleration.cpp
          OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
          DEPENDS_ON  core
          FOLDER      axom/core/examples )

    axom_add_executable(
          NAME        core_containers_cuda_on_ex
          SOURCES     core_containers.cpp
          OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
          DEPENDS_ON  core
          FOLDER      axom/core/examples )
endif()

axom_add_executable(
    NAME        core_utilities_ex
    SOURCES     core_utilities.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  core
    FOLDER      axom/core/examples )

if(AXOM_ENABLE_TESTS)
    axom_add_test(
        NAME    core_utilities
        COMMAND core_utilities_ex ${CMAKE_CURRENT_SOURCE_DIR} )
endif()

axom_add_executable(
    NAME        core_array_perf_ex
    SOURCES     core_array_perf.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  core
    FOLDER      axom/core/examples )

if(AXOM_ENABLE_TESTS)
    # Run the core array performance example on N ranks for each enabled policy
    set(_policies "seq")
    if(RAJA_FOUND)
        blt_list_append(TO _policies ELEMENTS "omp"  IF AXOM_ENABLE_OPENMP)
        blt_list_append(TO _policies ELEMENTS "cuda" IF AXOM_ENABLE_CUDA)
        blt_list_append(TO _policies ELEMENTS "hip"  IF AXOM_ENABLE_HIP)
    endif()
    foreach(_pol ${_policies})
        # sets the number of threads for the omp policy; leave empty for non-omp policies
        set(_num_threads)
        if(_pol STREQUAL "omp")
            set(_num_threads ${AXOM_TEST_NUM_OMP_THREADS})
        endif()

        axom_add_test(
            NAME            core_array_perf_1d_${_pol}
            COMMAND         core_array_perf_ex --policy ${_pol} --shape 1200000 -r 10
            NUM_OMP_THREADS ${_num_threads})
        axom_add_test(
            NAME            core_array_perf_2d_${_pol}
            COMMAND         core_array_perf_ex --policy ${_pol} --shape 1000 1200 -r 10
            NUM_OMP_THREADS ${_num_threads})
        axom_add_test(
            NAME            core_array_perf_3d_${_pol}
            COMMAND         core_array_perf_ex --policy ${_pol} --shape 100 100 120 -r 10
            NUM_OMP_THREADS ${_num_threads})
        if(NOT RAJA_FOUND)
            # RAJA provides support for up to 3D tiled nested loops,
            # so run 4D test only for non-RAJA builds.
            axom_add_test(
                NAME            core_array_perf_4d_${_pol}
                COMMAND         core_array_perf_ex --policy ${_pol} --shape 50 20 30 6
                NUM_OMP_THREADS ${_num_threads})
        endif()
    endforeach()
endif()

axom_add_executable(
    NAME        core_flatmap_perf_ex
    SOURCES     core_flatmap_perf.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  core
    FOLDER      axom/core/examples )
