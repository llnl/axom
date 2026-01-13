# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Unit tests for Klee component
#------------------------------------------------------------------------------

#
# Specify list of tests
#

set(gtest_klee_tests
    klee_config.cpp
    klee_dimensions.cpp
    klee_geometry.cpp
    klee_geometry_operators.cpp
    klee_geometry_operators_io.cpp
    klee_io.cpp
    klee_io_util.cpp
    klee_shape.cpp
    klee_shape_set.cpp
    klee_units.cpp
   )

axom_add_library(
        NAME        klee_test_utils
        SOURCES     KleeTestUtils.cpp
        HEADERS     KleeTestUtils.hpp KleeMatchers.hpp
        DEPENDS_ON  klee gtest gmock
        FOLDER      axom/klee/tests)

#
# Add gtest based tests
#
foreach(test ${gtest_klee_tests})
    get_filename_component( test_name ${test} NAME_WE )
    axom_add_executable(NAME        ${test_name}_test
                        SOURCES     ${test}
                        OUTPUT_DIR  ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON  klee_test_utils
                        FOLDER      axom/klee/tests )
    axom_add_test(NAME ${test_name}
                  COMMAND ${test_name}_test )
endforeach()

