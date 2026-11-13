# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
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
   sidre_checksum.cpp
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
   sidre_lifetime_Py.py
   sidre_import_Py.py
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
    set(_sidre_python_test_names "")
    foreach(test ${python_sidre_tests})
        get_filename_component( test_name ${test} NAME_WE )

        axom_add_python_test( NAME ${test_name}
                              SOURCE ${test}
                              OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                              )
        list(APPEND _sidre_python_test_names ${test_name})
    endforeach()

    if(AXOM_ENABLE_MPI)
        axom_add_python_test( NAME sidre_spio_Py
                              SOURCE sidre_spio_Py.py
                              OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                              NUM_MPI_TASKS 4
                              )
        list(APPEND _sidre_python_test_names sidre_spio_Py)
    endif()

    # CCE 21.0.0 can crash in libcraymp initialization when Python loads _sidre
    # after the ROCm profiling runtime. Preload CMake's detected Cray OpenMP
    # runtime to avoid this initializer ordering. This workaround is limited
    # to Sidre Python tests; direct/installed Python use still needs LD_PRELOAD.
    # Preserve existing preload entries from the configure-time environment.
    if(AXOM_ENABLE_HIP AND AXOM_ENABLE_OPENMP
       AND CMAKE_CXX_COMPILER_ID STREQUAL "CrayClang"
       AND CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL "21.0.0")
        if(NOT EXISTS "${OpenMP_craymp_LIBRARY}")
            message(FATAL_ERROR
                    "CCE 21.0.0 HIP/OpenMP Sidre Python tests require OpenMP_craymp_LIBRARY")
        endif()

        set(_cray_preload "${OpenMP_craymp_LIBRARY}")
        if(NOT "$ENV{LD_PRELOAD}" STREQUAL "")
            string(APPEND _cray_preload ":$ENV{LD_PRELOAD}")
        endif()
        set_property(TEST ${_sidre_python_test_names} APPEND PROPERTY ENVIRONMENT
                     "LD_PRELOAD=${_cray_preload}")
        unset(_cray_preload)
    endif()
    unset(_sidre_python_test_names)
endif()

#------------------------------------------------------------------------------
# Add MFEMSidreDataCollection tests
#------------------------------------------------------------------------------
if (MFEM_FOUND)
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
