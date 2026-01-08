# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Unit tests for Slic component
#------------------------------------------------------------------------------

#
# Serial Unit Tests
#
set(serial_slic_tests
    slic_asserts.cpp
    slic_fmt.cpp
    slic_interface.cpp
    slic_macros.cpp
    slic_scoped_abort.cpp
    slic_uninit.cpp )

set(slic_tests_depends slic gtest fmt)

foreach(test ${serial_slic_tests})
    get_filename_component( test_name ${test} NAME_WE )
    axom_add_executable(NAME       ${test_name}_test
                        SOURCES    ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${slic_tests_depends}
                        FOLDER     axom/slic/tests )
    axom_add_test(NAME    ${test_name}
                  COMMAND ${test_name}_test )
endforeach()


#
# Parallel Unit Tests
#
if (AXOM_ENABLE_MPI AND AXOM_ENABLE_LUMBERJACK)
    set(parallel_slic_tests
        slic_log_stream_status_monitor.cpp
        slic_macros_parallel.cpp
        slic_flush_perf.cpp)

    foreach(test ${parallel_slic_tests})
        get_filename_component( test_name ${test} NAME_WE )
        axom_add_executable(NAME       ${test_name}_test
                            SOURCES    ${test}
                            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                            DEPENDS_ON ${slic_tests_depends}
                            FOLDER     axom/slic/tests )
        axom_add_test(NAME          ${test_name}
                      COMMAND       ${test_name}_test 
                      NUM_MPI_TASKS 4)
    endforeach()
endif()

#
# Benchmarks
#
if (ENABLE_BENCHMARKS)
    set(slic_benchmarks
        slic_benchmark_asserts.cpp )

    foreach(test ${slic_benchmarks})
        get_filename_component( test_name ${test} NAME_WE )

        axom_add_executable(NAME       ${test_name}
                            SOURCES    ${test}
                            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                            DEPENDS_ON slic gbenchmark
                            FOLDER     axom/slic/benchmarks )

        blt_add_benchmark(
            NAME    ${test_name}
            COMMAND ${test_name} )
    endforeach()
endif()


#
# Fortran tests
#
if (ENABLE_FORTRAN)
    set(fruit_slic_tests
        slic_interface_F.f )

    foreach(test ${fruit_slic_tests})
        get_filename_component( test_name ${test} NAME_WE )
        axom_add_executable(NAME ${test_name}_test
                            SOURCES ${test}
                            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                            DEPENDS_ON slic fruit
                            FOLDER axom/slic/tests )
        axom_add_test(NAME    ${test_name}
                      COMMAND ${test_name}_test )
    endforeach()
endif()
