# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
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
if( AXOM_ENABLE_SIDRE AND AXOM_ENABLE_PRIMAL AND
    MFEM_FOUND AND AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION)

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
        COMMAND data_collection_util --min -1 -1 --max 1 1 --res 16 16 -p 1
        NUM_MPI_TASKS ${_nranks}
    )
    axom_add_test(
        NAME    data_collection_util_box3D
        COMMAND data_collection_util --min -1 -1 -1 --max 1 1 1 --res 16 16 16 -p 1
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

#------------------------------------------------------------------------------
# run_python_with_axom.sh.in is a utility that allows running python with
# nanobind-generated extensions in the PYTHONPATH.
# Based on Conduit's run_python_with_conduit.sh.in script.
#------------------------------------------------------------------------------
if(NANOBIND_FOUND)

    # gen python helper to build directory
    set(_PYEXT_DIR ${PROJECT_BINARY_DIR}/lib)

    axom_configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/run_python_with_axom.sh.in"
                         "${PROJECT_BINARY_DIR}/bin/run_python_with_axom.sh" @ONLY)

    # gen python helper to install directory
    set(_PYEXT_DIR ${CMAKE_INSTALL_PREFIX}/lib)

    axom_configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/run_python_with_axom.sh.in"
                         "${CMAKE_INSTALL_PREFIX}/bin/run_python_with_axom.sh" @ONLY)

    unset(_PYEXT_DIR)

    # Add unit test to check script
    if(AXOM_ENABLE_TESTS AND AXOM_ENABLE_SIDRE)
        axom_add_test(
            NAME    run_python_with_axom_build
            COMMAND ${PROJECT_BINARY_DIR}/bin/run_python_with_axom.sh -c "import pysidre")
    endif()
endif()
