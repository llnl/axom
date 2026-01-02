# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

set( mir_example_dependencies
     core
     slic
     mir
   )

# Break the different MIR cases into multiple files to speed up compilation.
set(mir_tutorial_simple_sources
    mir_tutorial_simple.cpp
    runMIR_seq.cpp
    runMIR_seq_tri.cpp
    runMIR_seq_quad.cpp
    runMIR_seq_hex.cpp
    runMIR_omp.cpp
    runMIR_omp_tri.cpp
    runMIR_omp_quad.cpp
    runMIR_omp_hex.cpp
    runMIR_cuda.cpp
    runMIR_cuda_tri.cpp
    runMIR_cuda_quad.cpp
    runMIR_cuda_hex.cpp
    runMIR_hip.cpp
    runMIR_hip_tri.cpp
    runMIR_hip_quad.cpp
    runMIR_hip_hex.cpp
    )

axom_add_executable(
    NAME       mir_tutorial_simple
    SOURCES    ${mir_tutorial_simple_sources}
    OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON ${mir_example_dependencies}
    FOLDER     axom/mir/examples
    )

if(AXOM_ENABLE_TESTS)
    set(_test_numbers 1 2 3 4 5)
    set (_policies "seq")
    if(RAJA_FOUND AND UMPIRE_FOUND)
        blt_list_append(TO _policies ELEMENTS "omp" IF AXOM_ENABLE_OPENMP)
        blt_list_append(TO _policies ELEMENTS "cuda" IF AXOM_ENABLE_CUDA)
        blt_list_append(TO _policies ELEMENTS "hip" IF AXOM_ENABLE_HIP)
    endif()

    foreach(_policy ${_policies})
        set(_num_threads)
        if(_policy STREQUAL "omp")
            set(_num_threads ${AXOM_TEST_NUM_OMP_THREADS})
        endif()

        foreach(_testnum ${_test_numbers})
            set(_testname "mir_tutorial_simple_${_policy}_${_testnum}")
            axom_add_test(
                NAME ${_testname}
                COMMAND mir_tutorial_simple
                    --test-case ${_testnum}
                    --policy ${_policy}
                    --disable-write
                NUM_OMP_THREADS ${_num_threads})

            set_tests_properties(${_testname} PROPERTIES
                PASS_REGULAR_EXPRESSION  "Reconstruction time:")
        endforeach()
    endforeach()
endif()
