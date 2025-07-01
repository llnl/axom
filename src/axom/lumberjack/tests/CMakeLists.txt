# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Tests for Lumberjack component
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Add serial GTests based tests
#------------------------------------------------------------------------------
set(lumberjack_serial_tests
    lumberjack_Lumberjack.hpp
    lumberjack_Message.hpp
    lumberjack_TextEqualityCombiner.hpp
    lumberjack_TextTagCombiner.hpp
    lumberjack_LineFileTagCombiner.hpp )

axom_add_executable(NAME       lumberjack_serial_tests
                    SOURCES    lumberjack_serial_main.cpp
                    HEADERS    ${lumberjack_serial_tests}
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON lumberjack gtest 
                    FOLDER     axom/lumberjack/tests )

foreach(test_suite ${lumberjack_serial_tests})
    get_filename_component(test_suite ${test_suite} NAME_WE)
    axom_add_test(NAME    ${test_suite}
                  COMMAND lumberjack_serial_tests --gtest_filter=${test_suite}* )
endforeach()


#------------------------------------------------------------------------------
# Add MPI GTest based tests
#------------------------------------------------------------------------------
set(lumberjack_mpi_tests
    lumberjack_BinaryCommunicator.hpp
    lumberjack_RootCommunicator.hpp
    lumberjack_NonCollectiveRootCommunicator.hpp )

axom_add_executable(NAME       lumberjack_mpi_tests
                    SOURCES    lumberjack_mpi_main.cpp
                    HEADERS    ${lumberjack_mpi_tests}
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON lumberjack gtest 
                    FOLDER     axom/lumberjack/tests )

foreach(test_suite ${lumberjack_mpi_tests})
    get_filename_component( test_suite ${test_suite} NAME_WE )
    foreach(task_count RANGE 1 8)
        axom_add_test(NAME          ${test_suite}_${task_count}
                      COMMAND       lumberjack_mpi_tests --gtest_filter=${test_suite}*
                      NUM_MPI_TASKS ${task_count} )
    endforeach()
endforeach()

#------------------------------------------------------------------------------
# Single source file tests
#------------------------------------------------------------------------------
axom_add_executable(NAME       lumberjack_speed_test
                    SOURCES    lumberjack_speedTest.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON lumberjack
                    FOLDER     axom/lumberjack/tests )

set(lumberjack_sample_input_dir ${PROJECT_SOURCE_DIR}/axom/lumberjack/tests/sampleInput)
axom_add_test(NAME          lumberjack_speedTest_binary
              COMMAND       lumberjack_speed_test b 10 ${lumberjack_sample_input_dir}/loremIpsum02
              NUM_MPI_TASKS 4)
              
axom_add_test(NAME          lumberjack_speedTest_root
              COMMAND       lumberjack_speed_test r 10 ${lumberjack_sample_input_dir}/loremIpsum02
              NUM_MPI_TASKS 4)

#------------------------------------------------------------------------------
# Performance test without combiners
#------------------------------------------------------------------------------
axom_add_executable(NAME       lumberjack_speed_test_no_combiners
                    SOURCES    lumberjack_speedTestNoCombiners.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON lumberjack
                    FOLDER     axom/lumberjack/tests )

axom_add_test(NAME          lumberjack_speedTestNoCombiners
              COMMAND       lumberjack_speed_test_no_combiners
              NUM_MPI_TASKS 1)
