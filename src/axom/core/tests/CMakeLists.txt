# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Tests for Axom's Core component
#
# Note: The tests in this file are consolidated into a few executables.
# In addition to adding files to the lists below, you'll also have  to #include 
# your headers in the appropriate 'main' cpp file, e.g. `core_serial_main.cpp`
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Specify list of serial tests
#------------------------------------------------------------------------------

set(core_serial_tests
    core_about.hpp
    core_array.hpp
    core_array_mapping.hpp
    core_array_for_all.hpp
    core_utilities.hpp
    core_bit_utilities.hpp
    core_device_hash.hpp
    core_execution_for_all.hpp
    core_execution_space.hpp
    core_map.hpp
    core_flatmap.hpp
    core_flatmap_for_all.hpp
    core_memory_management.hpp
    core_numeric_array.hpp
    core_numeric_limits.hpp
    core_Path.hpp
    core_stack_array.hpp
    core_static_array.hpp

    numerics_determinants.hpp
    numerics_eigen_solve.hpp
    numerics_eigen_sort.hpp
    numerics_floating_point_limits.hpp
    numerics_jacobi_eigensolve.hpp
    numerics_quadrature.hpp
    numerics_linear_solve.hpp
    numerics_lu.hpp
    numerics_matrix.hpp
    numerics_matvecops.hpp
    numerics_polynomial_solvers.hpp

    utils_endianness.hpp
    utils_fileUtilities.hpp
    utils_locale.hpp
    utils_stringUtilities.hpp
    utils_system.hpp
    utils_Timer.hpp
    utils_utilities.hpp
    )

set(core_test_depends
    core
    gtest
    fmt
    )

#------------------------------------------------------------------------------
# Serial GTests
#------------------------------------------------------------------------------

# Only enable this test in serial mode if AXOM_ENABLE_MPI is off
if (NOT AXOM_ENABLE_MPI)
  list(APPEND core_serial_tests core_types.hpp)
endif()

axom_add_executable(NAME       core_serial_tests
                    SOURCES    core_serial_main.cpp
                    HEADERS    ${core_serial_tests}
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON ${core_test_depends}
                    FOLDER     axom/core/tests )

foreach(test_suite ${core_serial_tests})
    set(_num_threads)
    if(test_suite MATCHES "(for_all|device)" AND AXOM_ENABLE_OPENMP)
        set(_num_threads ${AXOM_TEST_NUM_OMP_THREADS})
    endif()

    get_filename_component( test_suite ${test_suite} NAME_WE )
    axom_add_test( NAME            ${test_suite}_serial
                   COMMAND         core_serial_tests --gtest_filter=${test_suite}*
                   NUM_OMP_THREADS ${_num_threads})
endforeach()

#------------------------------------------------------------------------------
# Test for ItemCollection classes
#------------------------------------------------------------------------------

axom_add_executable(NAME core_collections_test
                    SOURCES core_collections.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON ${core_test_depends}
                    FOLDER axom/core/tests
                    )

axom_add_test( NAME    core_collections
               COMMAND core_collections_test
             )


#------------------------------------------------------------------------------
# MPI GTests
#------------------------------------------------------------------------------
if (AXOM_ENABLE_MPI)
  set( core_mpi_tests
       core_types.hpp
       core_shared_memory.hpp
       )

  axom_add_executable(NAME       core_mpi_tests
                      SOURCES    core_mpi_main.cpp
                      HEADERS    ${core_mpi_tests}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${core_test_depends} mpi
                      FOLDER     axom/core/tests )

  foreach( test_suite ${core_mpi_tests} )
    get_filename_component( test_suite ${test_suite} NAME_WE )
    axom_add_test( NAME          ${test_suite}_mpi
                   COMMAND       core_mpi_tests --gtest_filter=${test_suite}*
                   NUM_MPI_TASKS 2 )
  endforeach()
endif()

#------------------------------------------------------------------------------
# OpenMP GTests
#------------------------------------------------------------------------------
if (AXOM_ENABLE_OPENMP AND RAJA_FOUND)
  set( core_openmp_tests
       core_openmp_map.hpp
       )

  axom_add_executable(NAME       core_openmp_tests
                      SOURCES    core_openmp_main.cpp
                      HEADERS    ${core_mpi_tests}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${core_test_depends}
                      FOLDER     axom/core/tests )

  foreach( test_suite ${core_openmp_tests} )
    get_filename_component( test_suite ${test_suite} NAME_WE )
    axom_add_test( NAME            ${test_suite}_openmp
                   COMMAND         core_openmp_tests --gtest_filter=${test_suite}*
                   NUM_OMP_THREADS ${AXOM_TEST_NUM_OMP_THREADS} )
  endforeach()
endif()

#------------------------------------------------------------------------------
# Add a test to query properties of the configuration
#------------------------------------------------------------------------------
set(utils_config_test_depends core gtest )

blt_list_append( TO utils_config_test_depends ELEMENTS mfem IF MFEM_FOUND )

axom_add_executable(NAME       utils_config_test
                    SOURCES    utils_config.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON ${utils_config_test_depends}
                    FOLDER     axom/core/tests )

axom_add_test( NAME    utils_config
               COMMAND utils_config_test )

# Add flag to ignore unknown openmp pragmas in utils_config
set_property(TARGET utils_config_test
    APPEND_STRING PROPERTY COMPILE_FLAGS "${AXOM_DISABLE_OMP_PRAGMA_WARNINGS}")

# Have to ignore unused parameters for mfem
if(MFEM_FOUND)
    set_property(TARGET utils_config_test
        APPEND_STRING PROPERTY COMPILE_FLAGS "${AXOM_DISABLE_UNUSED_PARAMETER_WARNINGS}")
endif()

#------------------------------------------------------------------------------
# Test different annotation modes
#------------------------------------------------------------------------------
axom_add_executable(NAME       utils_annotations_test
                    SOURCES    utils_annotations.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON core gtest fmt cli11
                    FOLDER     axom/core/tests )

foreach(_mode none report counts gputx)
  set(_test_name utils_annotations_${_mode}_test)

  set(_num_ranks 1)
  if(AXOM_ENABLE_MPI)
    set(_num_ranks 2)
  endif()

  axom_add_test(NAME    ${_test_name}
                COMMAND utils_annotations_test -m ${_mode}
                NUM_MPI_TASKS ${_num_ranks})

  set_property(TEST ${_test_name}
               APPEND PROPERTY ENVIRONMENT  "CALI_LOG_VERBOSITY=2")
endforeach()

#------------------------------------------------------------------------------
# Benchmarks
#------------------------------------------------------------------------------
if (ENABLE_BENCHMARKS)

  set(core_benchmarks
      core_benchmark_array.cpp )

  foreach(test ${core_benchmarks})
    get_filename_component(test_name ${test} NAME_WE)

    axom_add_executable(NAME       ${test_name}
                        SOURCES    ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON fmt core slic gbenchmark cli11
                        FOLDER     axom/core/benchmarks)

    blt_add_benchmark(NAME    ${test_name}
                      COMMAND ${test_name} --benchmark_min_time=0.0001s)
  endforeach()
endif()

