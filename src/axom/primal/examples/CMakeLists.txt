# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Primal examples
#------------------------------------------------------------------------------

set(primal_examples 
    primal_introduction.cpp
    )

set(primal_example_depends
    primal
    fmt
    )

foreach ( example ${primal_examples} )

    get_filename_component( example_name ${example} NAME_WE )

    axom_add_executable(
        NAME        ${example_name}_ex
        SOURCES     ${example}
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${primal_example_depends}
        FOLDER      axom/primal/examples
        )

endforeach()

if (RAJA_FOUND AND UMPIRE_FOUND)

    axom_add_executable(
        NAME        hex_tet_volume_ex
        SOURCES     hex_tet_volume.cpp
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${primal_example_depends}
        FOLDER      axom/primal/examples
        )

    # Add unit tests
    if(AXOM_ENABLE_TESTS)

        # Run the hex_tet_volume_ex example with different raja policies
        # to check for completion

        set (_policies "seq")
        blt_list_append(TO _policies ELEMENTS "omp" IF AXOM_ENABLE_OPENMP)
        blt_list_append(TO _policies ELEMENTS "cuda" IF AXOM_ENABLE_CUDA)
        blt_list_append(TO _policies ELEMENTS "hip" IF AXOM_ENABLE_HIP)

        foreach(_policy ${_policies})
            # sets the number of threads for the omp policy; leave empty for non-omp policies
            set(_num_threads)
            if(_policy STREQUAL "omp")
                set(_num_threads ${AXOM_TEST_NUM_OMP_THREADS})
            endif()

            set(_testname "primal_hex_tet_volume_ex_${_policy}")
            axom_add_test(
              NAME           ${_testname}
              COMMAND        hex_tet_volume_ex --policy ${_policy}
              NUM_OMP_THREADS ${_num_threads})

            set_tests_properties(${_testname} PROPERTIES
                PASS_REGULAR_EXPRESSION  "Difference between sums")
        endforeach()
    endif()
endif()
