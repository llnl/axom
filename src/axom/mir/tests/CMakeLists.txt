# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Mir unit tests
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Specify list of tests
#------------------------------------------------------------------------------

set(gtest_mir_tests
    mir_elvira2d.cpp
    mir_elvira3d.cpp)

set(mir_tests_depends_on 
    slic
    mir
    bump
    gtest)

#-------------------------------------------------------------------------------------
# Policy-based tests are separated to reduce compile/link time and memory requirements
#-------------------------------------------------------------------------------------
set(policy_mir_tests
    mir_coupled
    mir_coupled3d
    mir_equiz2d
    mir_equiz3d)

foreach(test ${policy_mir_tests})
   blt_list_append(TO gtest_mir_tests ELEMENTS ${test}_seq.cpp)
   blt_list_append(TO gtest_mir_tests ELEMENTS ${test}_omp.cpp  IF AXOM_ENABLE_OPENMP)
   blt_list_append(TO gtest_mir_tests ELEMENTS ${test}_cuda.cpp IF AXOM_ENABLE_CUDA)
   blt_list_append(TO gtest_mir_tests ELEMENTS ${test}_hip.cpp  IF AXOM_ENABLE_HIP)
endforeach()

#------------------------------------------------------------------------------
# Add gtest based tests
#------------------------------------------------------------------------------
foreach(test ${gtest_mir_tests})
    get_filename_component( test_name ${test} NAME_WE )
    axom_add_executable( NAME ${test_name}_test
                         SOURCES ${test}
                         OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                         DEPENDS_ON ${mir_tests_depends_on}
                         FOLDER axom/mir/tests )

    axom_add_test( NAME ${test_name}
                   COMMAND ${test_name}_test )

    if(ENABLE_HIP AND test MATCHES "^mir_(coupled_hip|equiz2d_hip|elvira2d)\\.cpp$")
        axom_force_release_for_target(${test_name}_test)
    endif()
endforeach()

