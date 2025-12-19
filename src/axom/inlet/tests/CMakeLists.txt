# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Tests for Inlet component
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Specify list of tests
#------------------------------------------------------------------------------

axom_add_library(NAME        inlet_test_utils
                 SOURCES     inlet_test_utils.cpp
                 HEADERS     inlet_test_utils.hpp
                 DEPENDS_ON  inlet gtest
                 FOLDER      axom/inlet/tests )

# Add Serial GTests based tests
set(gtest_inlet_tests
    inlet_Reader.cpp
    inlet_restart.cpp
    inlet_errors.cpp )

blt_list_append(TO gtest_inlet_tests ELEMENTS inlet_function.cpp IF SOL_FOUND)
blt_list_append(TO gtest_inlet_tests ELEMENTS inlet_jsonschema_writer.cpp IF JSONSCHEMA_FOUND)

foreach(test ${gtest_inlet_tests})
    get_filename_component( test_name ${test} NAME_WE )
    axom_add_executable(NAME       ${test_name}_test
                        SOURCES    ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON inlet_test_utils
                        FOLDER     axom/inlet/tests )
    axom_add_test( NAME    ${test_name} 
                   COMMAND ${test_name}_test )
endforeach()

if(ENABLE_GMOCK)
    # Add Serial GMock based tests
    set(gmock_inlet_tests
        inlet_Inlet.cpp 
        inlet_object.cpp )

    foreach(test ${gmock_inlet_tests})
        get_filename_component( test_name ${test} NAME_WE )
        axom_add_executable(NAME       ${test_name}_test
                            SOURCES    ${test}
                            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                            DEPENDS_ON inlet_test_utils gmock
                            FOLDER     axom/inlet/tests )
        axom_add_test( NAME    ${test_name} 
                       COMMAND ${test_name}_test )
    endforeach()
endif()
