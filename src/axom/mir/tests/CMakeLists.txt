# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Mir unit tests
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Specify list of tests
#------------------------------------------------------------------------------

set(gtest_mir_tests
    mir_coupled.cpp
    mir_coupled3d.cpp
    mir_elvira2d.cpp
    mir_elvira3d.cpp
    mir_equiz2d.cpp
    mir_equiz3d.cpp
    )

set(mir_tests_depends_on 
    slic
    mir
    bump
    gtest)

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
endforeach()

