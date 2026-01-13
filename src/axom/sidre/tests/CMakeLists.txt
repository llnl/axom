# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Sidre unit tests
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Specify list of tests
#------------------------------------------------------------------------------

set(gtest_sidre_tests
   sidre_smoke.cpp
   sidre_buffer.cpp
   sidre_buffer_unit.cpp
   sidre_class.cpp
   sidre_datastore.cpp
   sidre_datastore_unit.cpp
   sidre_external.cpp
   sidre_group.cpp
   sidre_opaque.cpp
   sidre_view.cpp
   sidre_native_layout.cpp
   sidre_attribute.cpp
   sidre_mcarray.cpp
   sidre_read_write_userdefined_data.cpp
   )

set(gtest_sidre_C_tests
   sidre_smoke_C.cpp
   sidre_buffer_C.cpp
   sidre_external_C.cpp
   sidre_group_C.cpp
   sidre_opaque_C.cpp
   sidre_types_C.cpp
   sidre_view_C.cpp
   )

set(fruit_sidre_tests
   sidre_smoke_F.f
   sidre_group_F.F
   sidre_buffer_F.f
   sidre_view_F.f
   sidre_external_F.f
   sidre_opaque_F.f
   sidre_allocatable_F.f
   )

set(python_sidre_tests
   sidre_smoke_Py.py
   sidre_datastore_unit_Py.py
   sidre_group_Py.py
   sidre_buffer_Py.py
   sidre_view_Py.py
   sidre_external_Py.py
   sidre_attribute_Py.py
   )

set(sidre_gtests_depends_on sidre fmt gtest)

#------------------------------------------------------------------------------
# Add gtest C++ tests
#------------------------------------------------------------------------------
foreach(test ${gtest_sidre_tests})
    get_filename_component( test_name ${test} NAME_WE )
    axom_add_executable(NAME ${test_name}_test
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${sidre_gtests_depends_on}
                        FOLDER axom/sidre/tests
                        )

    axom_add_test( NAME    ${test_name}
                   COMMAND ${test_name}_test
                   )
endforeach()

#------------------------------------------------------------------------------
# Add gtest C tests
#------------------------------------------------------------------------------
if (ENABLE_FORTRAN)
    # Sidre's C API is only built to provide a Fortran interface.
    foreach(test ${gtest_sidre_C_tests})
        get_filename_component( test_name ${test} NAME_WE )
        axom_add_executable(NAME ${test_name}_test
                            SOURCES ${test}
                            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                            DEPENDS_ON sidre gtest
                            FOLDER axom/sidre/tests
                            )
        axom_add_test( NAME    ${test_name}
                       COMMAND ${test_name}_test
                       )
    endforeach()
endif()

#------------------------------------------------------------------------------
# Add Fortran tests
#------------------------------------------------------------------------------
if (ENABLE_FORTRAN)
    foreach(test ${fruit_sidre_tests})
        get_filename_component( test_name ${test} NAME_WE )
        axom_add_executable(NAME ${test_name}_test
                            SOURCES ${test}
                            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                            DEPENDS_ON sidre fruit
                            FOLDER axom/sidre/tests
                            )
        axom_add_test( NAME    ${test_name}
                       COMMAND ${test_name}_test
                       )
    endforeach()
endif()

#------------------------------------------------------------------------------
# Add Python tests
#------------------------------------------------------------------------------
if(NANOBIND_FOUND AND AXOM_ENABLE_PYTHON_TESTS)
    foreach(test ${python_sidre_tests})
        get_filename_component( test_name ${test} NAME_WE )

        # Copy python test file to build
        axom_configure_file ("${test}"
                             "${TEST_OUTPUT_DIRECTORY}/${test}" COPYONLY)

        # Run unit test with pytest
        # "-p no:cacheprovider" disables caching
        #
        # Need to invoke with "python3 -m pytest" to get generated nanobind library added to search path:
        # https://docs.pytest.org/en/latest/explanation/pythonpath.html#invoking-pytest-versus-python-m-pytest
        add_test (NAME    ${test_name}
          COMMAND ${Python_EXECUTABLE} -m pytest -s -p no:cacheprovider ${TEST_OUTPUT_DIRECTORY}/${test}
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/lib
        )
    endforeach()
endif()

#------------------------------------------------------------------------------
# Add MFEMSidreDataCollection tests
#------------------------------------------------------------------------------
if (MFEM_FOUND AND AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION)
    axom_add_executable(NAME sidre_mfem_datacollection_test
                        SOURCES sidre_mfem_datacollection.cpp
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON sidre gtest
                        FOLDER axom/sidre/tests
                        )
    if (AXOM_ENABLE_MPI)
        if(MFEM_USE_MPI)
            set(_nranks 2)
        else()
            set(_nranks 1)
        endif()

        axom_add_test( NAME sidre_mfem_datacollection
                       COMMAND sidre_mfem_datacollection_test
                       NUM_MPI_TASKS ${_nranks}
                       )
    else()
        axom_add_test( NAME    sidre_mfem_datacollection
                       COMMAND sidre_mfem_datacollection_test
                       )
    endif()
endif()


if (AXOM_ENABLE_MPI)
    add_subdirectory(spio)
endif()
