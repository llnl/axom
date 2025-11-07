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
    NAME       mir_heavily_mixed
    SOURCES    mir_heavily_mixed.cpp HMApplication.cpp runMIR_seq.cpp runMIR_omp.cpp runMIR_cuda.cpp runMIR_hip.cpp
    OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON ${mir_example_dependencies}
    FOLDER     axom/mir/examples
    )

if(AXOM_ENABLE_TESTS)
    set (_policies "seq")
    if(RAJA_FOUND AND UMPIRE_FOUND)
        blt_list_append(TO _policies ELEMENTS "omp" IF AXOM_ENABLE_OPENMP)
        blt_list_append(TO _policies ELEMENTS "cuda" IF AXOM_ENABLE_CUDA)
        blt_list_append(TO _policies ELEMENTS "hip" IF AXOM_ENABLE_HIP)
    endif()

    foreach(_policy ${_policies})
        set(_testname "mir_heavily_mixed_2D_${_policy}")
        axom_add_test(
            NAME ${_testname}
            COMMAND mir_heavily_mixed
                --dimension 2
                --dims 50,50
                --refinement 10
                --nmats 10
                --policy ${_policy}
                --disable-write
            )

        set_tests_properties(${_testname} PROPERTIES
            PASS_REGULAR_EXPRESSION  "Material interface reconstruction time:")

        set(_testname "mir_heavily_mixed_3D_${_policy}")
        axom_add_test(
            NAME ${_testname}
            COMMAND mir_heavily_mixed
                --dimension 3
                --dims 50,50,10
                --refinement 10
                --nmats 10
                --policy ${_policy}
                --disable-write
            )

        set_tests_properties(${_testname} PROPERTIES
            PASS_REGULAR_EXPRESSION  "Material interface reconstruction time:")
    endforeach()
endif()
