# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#------------------------------------------------------------------------------
# This directory contains source code for utility executables that utilize
# Axom or that support Axom
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# convert_protocol is a utility that converts
# a datastore from the sidre_hdf5 protocol
# to any of the other supported protocols.
# It requires slam, sidre, mpi and hdf5.
#------------------------------------------------------------------------------
if( AXOM_ENABLE_SLAM AND AXOM_ENABLE_SIDRE AND AXOM_ENABLE_MPI AND HDF5_FOUND )
    set(datastore_converter_sources
        convert_sidre_protocol.cpp
    )
    set(datastore_converter_depends
        slam
        sidre
        cli11
        fmt)

    blt_list_append(TO datastore_converter_depends ELEMENTS scr IF SCR_FOUND)

    axom_add_executable(
        NAME       convert_sidre_protocol
        SOURCES    ${datastore_converter_sources}
        OUTPUT_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
        DEPENDS_ON ${datastore_converter_depends}
        FOLDER     axom/examples
        )
    install(TARGETS              convert_sidre_protocol
            DESTINATION          bin
            )
endif()

#------------------------------------------------------------------------------
# mesh_tester is a utility that tests a mesh for a variety of problems.
#------------------------------------------------------------------------------
if(AXOM_ENABLE_QUEST)
    set(mesh_tester_sources
        mesh_tester.cpp
    )
    set(mesh_tester_depends
        quest
        )

    blt_list_append(TO mesh_tester_depends ELEMENTS umpire IF UMPIRE_FOUND)
    blt_list_append(TO mesh_tester_depends ELEMENTS RAJA IF RAJA_FOUND)

    axom_add_executable(
        NAME        mesh_tester
        SOURCES     ${mesh_tester_sources}
        OUTPUT_DIR  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${mesh_tester_depends}
        FOLDER      axom/examples
        )
    install(TARGETS              mesh_tester
            DESTINATION          bin
            )

    # Add separable compilation flag
    if (AXOM_ENABLE_HIP)
        blt_add_target_compile_flags(TO mesh_tester FLAGS "-fgpu-rdc")
        blt_add_target_link_flags(TO mesh_tester FLAGS "-fgpu-rdc")
    endif()

    # Add unit tests
    if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR AND RAJA_FOUND AND UMPIRE_FOUND)

        # Run the mesh_tester with the plane_simp_problems.stl example and
        # different spatial indexes, raja policies

        set(plane_data_dir "${AXOM_DATA_DIR}/quest/plane_simp_problems.stl")

        set(_methods "bvh" "implicit" "uniform")

        set (_policies "raja_seq")
        blt_list_append(TO _policies ELEMENTS "raja_omp" IF AXOM_ENABLE_OPENMP)
        blt_list_append(TO _policies ELEMENTS "raja_cuda" IF AXOM_ENABLE_CUDA)
        blt_list_append(TO _policies ELEMENTS "raja_hip" IF AXOM_ENABLE_HIP)

        foreach(_method ${_methods})
            foreach(_policy ${_policies})

                set(_testname "mesh_tester_${_method}_${_policy}")
                axom_add_test(
                  NAME ${_testname}
                  COMMAND mesh_tester
                          --infile ${plane_data_dir}
                          --method ${_method}
                          --policy ${_policy}
                )

                set_tests_properties(${_testname} PROPERTIES
                    PASS_REGULAR_EXPRESSION  "5 intersecting tri pairs")
            endforeach()
        endforeach()
    endif()

endif()

#------------------------------------------------------------------------------
# data_collection_util is a utility that generates blueprint-conforming sidre
# output files from an input mfem mesh
#------------------------------------------------------------------------------
if(AXOM_ENABLE_SIDRE AND AXOM_ENABLE_PRIMAL AND MFEM_FOUND)

    set(data_collection_util
        data_collection_util.cpp
    )
    set(data_collection_util_depends
        sidre
        primal
        quest
        cli11
        fmt
        )

    blt_list_append(TO data_collection_util_depends ELEMENTS mpi IF AXOM_ENABLE_MPI)

    axom_add_executable(
        NAME        data_collection_util
        SOURCES     ${data_collection_util}
        OUTPUT_DIR  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${data_collection_util_depends}
        FOLDER      axom/examples
        )

    set(_nranks 1)
    if(AXOM_ENABLE_MPI)
        set(_nranks 2)
    endif()
        
    axom_add_test(
        NAME    data_collection_util_box2D
        COMMAND data_collection_util --min -1 -1 --max 1 1 --resolution 16 16 -p 1
        NUM_MPI_TASKS ${_nranks}
    )
    axom_add_test(
        NAME    data_collection_util_box3D
        COMMAND data_collection_util --min -1 -1 -1 --max 1 1 1 --resolution 16 16 16 -p 1
        NUM_MPI_TASKS ${_nranks}
    )

    install(TARGETS              data_collection_util
            DESTINATION          bin
            )

endif()

#------------------------------------------------------------------------------
# mesh_converter is a utility that converts between mesh formats
#------------------------------------------------------------------------------
if(AXOM_ENABLE_QUEST)
    set(mesh_converter_sources
        mesh_converter.cpp
    )
    set(mesh_converter_depends
        axom
        )

    axom_add_executable(
        NAME        mesh_converter
        SOURCES     ${mesh_converter_sources}
        OUTPUT_DIR  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${mesh_converter_depends}
        FOLDER      axom/examples
        )
    install(TARGETS              mesh_converter
            DESTINATION          bin
            )

    # Add unit tests
    if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR)
        axom_add_test(
            NAME    mesh_converter_cup
            COMMAND mesh_converter -i "${AXOM_DATA_DIR}/quest/cup.proe" -o cup.stl)
        set_tests_properties(mesh_converter_cup 
            PROPERTIES PASS_REGULAR_EXPRESSION  "Retained 272 boundary faces")
    endif()
endif()


if(NANOBIND_FOUND)
    #--------------------------------------------------------------------------
    # run_python_with_axom is a utility that allows running python with
    # nanobind-generated extensions in the PYTHONPATH.
    # Based on Conduit's run_python_with_conduit.sh.in script.
    #--------------------------------------------------------------------------
    if(WIN32)
        set(_run_python_with_axom_bat_template "${CMAKE_CURRENT_SOURCE_DIR}/run_python_with_axom.bat.in")
        set(_run_python_with_axom_sh_template "${CMAKE_CURRENT_SOURCE_DIR}/run_python_with_axom_gitbash.sh.in")
        set(_run_python_with_axom_launcher_template "${CMAKE_CURRENT_SOURCE_DIR}/run_python_with_axom.py.in")
        set(_run_python_with_axom_build "${PROJECT_BINARY_DIR}/bin/run_python_with_axom.bat")
        set(_run_python_with_axom_sh_build "${PROJECT_BINARY_DIR}/bin/run_python_with_axom.sh")
        set(_run_python_with_axom_launcher_build "${PROJECT_BINARY_DIR}/bin/run_python_with_axom.py")
        set(_run_python_with_axom_install "${CMAKE_INSTALL_PREFIX}/bin/run_python_with_axom.bat")
        set(_run_python_with_axom_sh_install "${CMAKE_INSTALL_PREFIX}/bin/run_python_with_axom.sh")
        set(_run_python_with_axom_launcher_install "${CMAKE_INSTALL_PREFIX}/bin/run_python_with_axom.py")
    else()
        set(_run_python_with_axom_sh_template "${CMAKE_CURRENT_SOURCE_DIR}/run_python_with_axom.sh.in")
        set(_run_python_with_axom_build "${PROJECT_BINARY_DIR}/bin/run_python_with_axom.sh")
        set(_run_python_with_axom_install "${CMAKE_INSTALL_PREFIX}/bin/run_python_with_axom.sh")
    endif()

    # gen python helper to build directory
    set(_PYEXT_DIR ${PROJECT_BINARY_DIR}/lib)
    set(_PYDLL_DIR ${PROJECT_BINARY_DIR}/bin)
    set(_PYLAUNCHER ${_run_python_with_axom_launcher_build})

    if(WIN32)
        axom_configure_file ("${_run_python_with_axom_launcher_template}"
                             "${_run_python_with_axom_launcher_build}" @ONLY)
        axom_configure_file ("${_run_python_with_axom_bat_template}"
                             "${_run_python_with_axom_build}" @ONLY)
        axom_configure_file ("${_run_python_with_axom_sh_template}"
                             "${_run_python_with_axom_sh_build}" @ONLY)
    else()
        axom_configure_file ("${_run_python_with_axom_sh_template}"
                             "${_run_python_with_axom_build}" @ONLY)
    endif()

    # gen python helper to install directory
    set(_PYEXT_DIR ${CMAKE_INSTALL_PREFIX}/lib)
    set(_PYDLL_DIR ${CMAKE_INSTALL_PREFIX}/bin)
    set(_PYLAUNCHER ${_run_python_with_axom_launcher_install})

    if(WIN32)
        axom_configure_file ("${_run_python_with_axom_launcher_template}"
                             "${_run_python_with_axom_launcher_install}" @ONLY)
        axom_configure_file ("${_run_python_with_axom_bat_template}"
                             "${_run_python_with_axom_install}" @ONLY)
        axom_configure_file ("${_run_python_with_axom_sh_template}"
                             "${_run_python_with_axom_sh_install}" @ONLY)
    else()
        axom_configure_file ("${_run_python_with_axom_sh_template}"
                             "${_run_python_with_axom_install}" @ONLY)
    endif()

    unset(_PYEXT_DIR)
    unset(_PYDLL_DIR)
    unset(_PYLAUNCHER)
    unset(_run_python_with_axom_bat_template)
    unset(_run_python_with_axom_sh_template)
    unset(_run_python_with_axom_launcher_template)
    unset(_run_python_with_axom_install)
    unset(_run_python_with_axom_sh_install)
    unset(_run_python_with_axom_launcher_install)
    unset(_run_python_with_axom_sh_build)
    unset(_run_python_with_axom_launcher_build)

    # Smoke tests for the script.
    # The wrapper provides the runtime environment only:
    # Axom extensions, conduit (Node interop), numpy (ndarray returns), plus mpi4py in MPI configurations. 
    # nanobind is a build-time dependency (statically linked into the extensions).
    # pytest/pluggy/iniconfig are test-harness dependencies.
    if(AXOM_ENABLE_TESTS AND AXOM_ENABLE_SIDRE)
        axom_run_python_with_axom_command(_run_python_with_axom_test)
        axom_add_test(
            NAME    run_python_with_axom_build
            COMMAND ${_run_python_with_axom_test} -c "import pysidre, conduit, numpy")

        if(AXOM_ENABLE_PYTHON_TESTS)
            # The pytest harness is provided per-test via the ENVIRONMENT property;
            # verify that the wrapper + injected ENVIRONMENT combination resolves.
            axom_add_test(
                NAME    run_python_with_axom_pytest_harness
                COMMAND ${_run_python_with_axom_test} -c "import pytest")
            axom_python_test_environment(_py_test_env)
            if(_py_test_env)
                set_property(TEST run_python_with_axom_pytest_harness
                             APPEND
                             PROPERTY ENVIRONMENT "${_py_test_env}")
            endif()
            unset(_py_test_env)
        endif()
        unset(_run_python_with_axom_test)
    endif()

    #--------------------------------------------------------------------------
    # convert_sidre_protocol.py is a python version of the
    # convert_sidre_protocol.cpp utility
    #--------------------------------------------------------------------------
    axom_configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/convert_sidre_protocol.py"
                         "${PROJECT_BINARY_DIR}/bin/convert_sidre_protocol.py" COPYONLY)

    axom_configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/convert_sidre_protocol.py"
                         "${CMAKE_INSTALL_PREFIX}/bin/convert_sidre_protocol.py" COPYONLY)

    if(AXOM_ENABLE_MPI AND AXOM_ENABLE_TESTS AND AXOM_DATA_DIR)

        set(box_dir "${AXOM_DATA_DIR}/quest/box_2D_r3.root")

        set(_testname "convert_sidre_protocol_py")
        axom_run_python_with_axom_command(_run_python_with_axom_convert)
        axom_add_test(
          NAME ${_testname}
          COMMAND ${_run_python_with_axom_convert}
                  ${PROJECT_BINARY_DIR}/bin/convert_sidre_protocol.py
                  --input ${box_dir}
                  --output csp_output
                  --protocol json
                  --verbose
        )

        set_tests_properties(${_testname} PROPERTIES
            PASS_REGULAR_EXPRESSION  "Writing out datastore")
        unset(_run_python_with_axom_convert)
    endif()
    unset(_run_python_with_axom_build)
endif()
