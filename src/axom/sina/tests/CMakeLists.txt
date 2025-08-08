# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Sina unit tests
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Add gtest C++ tests
#------------------------------------------------------------------------------
set(gtest_sina_tests
    sina_File.cpp
   )

set(sina_gtests_depends_on gtest sina conduit::conduit)

foreach(test ${gtest_sina_tests})
    get_filename_component( test_name ${test} NAME_WE )
    axom_add_executable(NAME ${test_name}_test
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${sina_gtests_depends_on}
                        FOLDER axom/sina/tests
                        )

    axom_add_test( NAME    ${test_name}
                   COMMAND ${test_name}_test
                   )
endforeach()


if (ENABLE_GMOCK)

    #------------------------------------------------------------------------------
    # Create test utilities library for Sina
    #------------------------------------------------------------------------------

    # Define Sina test utility sources and headers
    set(sina_test_utils_sources
        SinaMatchers.cpp
        TestRecord.cpp
        )
    set(sina_test_utils_headers
        SinaMatchers.hpp
        TestRecord.hpp
        )

    set(sina_test_utils_depends_on gmock ${sina_gtests_depends_on})

    # Create a library for the test utilities so they can be used
    axom_add_library(
        NAME       sina_test_utils
        SOURCES    ${sina_test_utils_sources}
        HEADERS    ${sina_test_utils_headers}
        DEPENDS_ON ${sina_test_utils_depends_on}
        FOLDER     axom/sina/tests)

    #------------------------------------------------------------------------------
    # Add gmock C++ tests
    #------------------------------------------------------------------------------

    set(gmock_sina_tests
        sina_ConduitUtil.cpp
        sina_Curve.cpp
        sina_CurveSet.cpp
        sina_DataHolder.cpp
        sina_Datum.cpp
        sina_Document.cpp
        sina_ID.cpp
        sina_Record.cpp
        sina_Relationship.cpp
        sina_Run.cpp
        )

    set(sina_gmock_depends_on sina_test_utils)

    # Add tests using Adiak if necessary and Adiak dependency
    blt_list_append( TO gmock_sina_tests ELEMENTS sina_AdiakWriter.cpp IF AXOM_USE_ADIAK)
    blt_list_append( TO sina_gmock_depends_on ELEMENTS adiak::adiak IF AXOM_USE_ADIAK )

    foreach(test ${gmock_sina_tests})
        get_filename_component( test_name ${test} NAME_WE )
        axom_add_executable(NAME ${test_name}_test
                            SOURCES ${test}
                            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                            DEPENDS_ON ${sina_gmock_depends_on}
                            FOLDER axom/sina/tests
                            )

        axom_add_test( NAME    ${test_name}
                       COMMAND ${test_name}_test
                       )
    endforeach()
endif()

#------------------------------------------------------------------------------
# Add fortran integration test
#------------------------------------------------------------------------------
if (ENABLE_FORTRAN AND DEFINED Python_EXECUTABLE)
    configure_file(
        ${CMAKE_SOURCE_DIR}/axom/sina/tests/test_fortran_integration.py
        ${TEST_OUTPUT_DIRECTORY}/test_fortran_integration.py
        COPYONLY
    )
    configure_file(
        ${CMAKE_SOURCE_DIR}/axom/sina/interface/sina_schema.json
        ${TEST_OUTPUT_DIRECTORY}/sina_schema.json
        COPYONLY
    )

        # Define Python-friendly boolean for HDF5 guarding in Fortran Test.
        if(AXOM_USE_HDF5)
            set(AXOM_USE_HDF5_PY "True")
        else()
            set(AXOM_USE_HDF5_PY "False")
        endif()

        # Generate config.py in the test output directory
        configure_file(
            ${CMAKE_CURRENT_SOURCE_DIR}/config.py.in
            ${TEST_OUTPUT_DIRECTORY}/config.py
            @ONLY
        )

    axom_add_test( NAME     sina_fortran_integration_test
                COMMAND  ${Python_EXECUTABLE} ${TEST_OUTPUT_DIRECTORY}/test_fortran_integration.py -bd ${PROJECT_BINARY_DIR}
    )
endif()
