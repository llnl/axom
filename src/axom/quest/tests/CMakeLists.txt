# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Quest tests
#------------------------------------------------------------------------------

set(quest_tests
    quest_all_nearest_neighbors.cpp
    quest_inout_octree.cpp
    quest_inout_quadtree.cpp
    quest_linearize_curves.cpp
    quest_signed_distance.cpp
    quest_discretize.cpp
    quest_pro_e_reader.cpp
    quest_stl_reader.cpp
    quest_stl_writer.cpp
    quest_vertex_weld.cpp
   )

blt_list_append(TO       quest_tests
                IF       C2C_FOUND
                ELEMENTS quest_c2c_reader.cpp)

# Optionally, add tests that require AXOM_DATA_DIR
blt_list_append(TO       quest_tests
                IF       AXOM_DATA_DIR
                ELEMENTS quest_meshtester.cpp)

set(quest_tests_depends
    quest
    gtest
   )

foreach(test ${quest_tests})
    get_filename_component( test_name ${test} NAME_WE )

    axom_add_executable(
        NAME       ${test_name}_test
        SOURCES    ${test}
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON ${quest_tests_depends}
        FOLDER     axom/quest/tests
        )

    # sets the number of threads for the quest_discretize test
    set(_num_threads)
    if(AXOM_ENABLE_OPENMP AND "${test_name}" STREQUAL "quest_discretize")
        set(_num_threads ${AXOM_TEST_NUM_OMP_THREADS})
    endif()

    axom_add_test(
        NAME ${test_name}
        COMMAND ${test_name}_test
        NUM_OMP_THREADS ${_num_threads})
endforeach()

#------------------------------------------------------------------------------
# Tests that use Conduit when available
#------------------------------------------------------------------------------

if(CONDUIT_FOUND AND AXOM_DATA_DIR)
    axom_add_executable(
            NAME        quest_mesh_view_util_test
            SOURCES     quest_mesh_view_util.cpp
            OUTPUT_DIR  ${TEST_OUTPUT_DIRECTORY}
            DEPENDS_ON  ${quest_tests_depends} conduit
            FOLDER      axom/quest/tests
            )
    set(quest_data_dir  ${AXOM_DATA_DIR}/quest)

    # These examples currently segfault on windows without MPI when loading the test data into conduit
    if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR AND NOT WIN32)
        axom_add_test(
            NAME    quest_mesh_view_util
            COMMAND quest_mesh_view_util_test --verbose)
        unset(_nranks)
        unset(_test)
    endif()

endif()

#------------------------------------------------------------------------------
# Tests that use MFEM when available
#------------------------------------------------------------------------------

if(MFEM_FOUND AND AXOM_ENABLE_SIDRE AND AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION)
    function(add_quest_mfem_test test_name)
        axom_add_executable(
            NAME       ${test_name}
            SOURCES    ${test_name}.cpp
            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
            DEPENDS_ON ${quest_tests_depends} mfem
            FOLDER      axom/quest/tests
            )

        string(STRIP "${AXOM_DISABLE_UNUSED_PARAMETER_WARNINGS}" MFEM_COMPILE_FLAGS)

        if (AXOM_ENABLE_CUDA)
            set(MFEM_COMPILE_FLAGS "-Xcompiler=${MFEM_COMPILE_FLAGS}")
        endif()

        blt_add_target_compile_flags( TO ${test_name} FLAGS "${MFEM_COMPILE_FLAGS}" )

        axom_add_test(
            NAME ${test_name}
            COMMAND ${test_name}
            )
    endfunction()

    add_quest_mfem_test(quest_point_in_cell_mfem)
    if(AXOM_DATA_DIR)
        add_quest_mfem_test(quest_mfem_reader)
    endif()
endif()


#------------------------------------------------------------------------------
# Tests that use MPI when available
#------------------------------------------------------------------------------

set(quest_mpi_tests
    quest_initialize.cpp
    quest_signed_distance_interface.cpp
    )


# Tests depend on MPI
blt_list_append(TO quest_mpi_tests
                IF AXOM_ENABLE_MPI
                ELEMENTS quest_pro_e_reader_parallel.cpp)

# Optionally, add tests that require AXOM_DATA_DIR
blt_list_append(TO       quest_mpi_tests
                IF       AXOM_DATA_DIR
                ELEMENTS quest_inout_interface.cpp )

foreach(test ${quest_mpi_tests})
    get_filename_component( test_name ${test} NAME_WE )
    axom_add_executable(
        NAME       ${test_name}_test
        SOURCES    ${test_name}.cpp
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON ${quest_tests_depends}
        FOLDER     axom/quest/tests
    )

    # sets the number of MPI ranks; quest_inout_interface can only run w/ 1 rank
    set(_numMPITasks 1)
    if (AXOM_ENABLE_MPI)
        set(_numMPITasks 4)
    endif()
    if(${test_name} STREQUAL "quest_inout_interface")
        set(_numMPITasks 1)
    endif()

    axom_add_test(
        NAME ${test_name}
        COMMAND ${test_name}_test
        NUM_MPI_TASKS ${_numMPITasks})
endforeach()

#------------------------------------------------------------------------------
# These tests require MPI+MFEM+C2C+Axom data+...
#------------------------------------------------------------------------------
if(ENABLE_MPI AND MFEM_FOUND AND MFEM_USE_MPI
              AND AXOM_ENABLE_SIDRE AND AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION
              AND AXOM_ENABLE_KLEE AND AXOM_DATA_DIR
              AND C2C_FOUND)
    set(quest_mpi_mfem_c2c_data_tests
        quest_intersection_shaper.cpp
        quest_intersection_shaper_dynamic.cpp
        quest_sampling_shaper.cpp
       )
    foreach(test ${quest_mpi_mfem_c2c_data_tests})
        get_filename_component( test_name ${test} NAME_WE )
        axom_add_executable(
            NAME ${test_name}_test
            SOURCES ${test_name}.cpp
            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
            DEPENDS_ON ${quest_tests_depends} mfem
            FOLDER axom/quest/tests
        )
    endforeach()

    # Invoke the intersection shaper test program as many tests that run a subset of the tests.
    set(_policies "seq")
    if(RAJA_FOUND)
      blt_list_append(TO _policies ELEMENTS "omp" IF AXOM_ENABLE_OPENMP)
      if(AXOM_USE_UMPIRE)
        blt_list_append(TO _policies ELEMENTS "cuda" IF AXOM_ENABLE_CUDA)
        blt_list_append(TO _policies ELEMENTS "hip" IF AXOM_ENABLE_HIP)
      endif()
    endif()

    foreach(_pol ${_policies})
      # sets the number of threads
      set(_num_threads)
      if(_pol STREQUAL "omp")
          set(_num_threads ${AXOM_TEST_NUM_OMP_THREADS})
      endif()

      # For seq, limit the cases since the tests are somewhat slow.
      set(caseNumbers 1 2 3 4 5)
      if(_pol STREQUAL "seq")
        set(caseNumbers 5)
      endif()

      foreach(caseNumber ${caseNumbers})
        axom_add_test(NAME quest_intersection_shaper_${caseNumber}_${_pol}
          COMMAND quest_intersection_shaper_test --casenumber ${caseNumber} --policy ${_pol}
          NUM_MPI_TASKS 1
          NUM_OMP_THREADS ${_num_threads})
      endforeach()
    endforeach()

    # quest_intersection_shaper_dynamic does not take a policy argument, 
    # but uses OpenMP when enabled
    set(_dynamic_num_threads)
    if(AXOM_ENABLE_OPENMP)
        set(_dynamic_num_threads ${AXOM_TEST_NUM_OMP_THREADS})
    endif()

    axom_add_test(NAME quest_intersection_shaper_dynamic
        COMMAND quest_intersection_shaper_dynamic_test
        NUM_MPI_TASKS 1
        NUM_OMP_THREADS ${_dynamic_num_threads}
    )

    axom_add_test(NAME quest_sampling_shaper
        COMMAND quest_sampling_shaper_test
        NUM_MPI_TASKS 1
    )

endif()

#------------------------------------------------------------------------------
# Geometry clipping tests require Conduit and RAJA
#------------------------------------------------------------------------------
if(CONDUIT_FOUND AND RAJA_FOUND AND AXOM_ENABLE_SIDRE)
    axom_add_executable(
        NAME        quest_mesh_clipper_test
        SOURCES     quest_mesh_clipper.cpp
        OUTPUT_DIR  ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${quest_tests_depends}
        FOLDER      axom/quest/tests
        )

    # 3D shaping tests.  No 2D shaping in this feature, at least for now.
    set(_nranks 1)

    # Run the geometry clipping test on with each enabled policy
    set(_policies "seq")
    if(RAJA_FOUND)
      blt_list_append(TO _policies ELEMENTS "omp"  IF AXOM_ENABLE_OPENMP)
      blt_list_append(TO _policies ELEMENTS "cuda" IF AXOM_ENABLE_CUDA)
      blt_list_append(TO _policies ELEMENTS "hip"  IF AXOM_ENABLE_HIP)
    endif()

    set(_testgeoms "plane" "tet" "hex" "sphere" "plane,hex,tet")
    set(_meshTypes "bpSidre" "bpConduit")
    foreach(_meshType ${_meshTypes})
      foreach(_policy ${_policies})
        # sets the number of threads for the omp policy; leave empty for non-omp policies
        set(_num_threads)
        if(_policy STREQUAL "omp")
            set(_num_threads ${AXOM_TEST_NUM_OMP_THREADS})
        endif()

        foreach(_testgeom ${_testgeoms})
          set(_testname "quest_mesh_clipper_${_meshType}_${_policy}_${_testgeom}")
          axom_add_test(
            NAME ${_testname}
            COMMAND quest_mesh_clipper_test
                        --policy ${_policy}
                        --verbose
                        --testGeom ${_testgeom}
                        -o ${_testname}
                        --refinements 5
                        --scale .99 .99 .99
                        --dir 8 4 2
                        --meshType ${_meshType}
                        inline_mesh --min -2 -2 -2 --max 2 2 2 --res 16 16 16
            NUM_MPI_TASKS ${_nranks}
            NUM_OMP_THREADS ${_num_threads})
        endforeach()
      endforeach()
    endforeach()

endif()

#------------------------------------------------------------------------------
# Regression tests for quest signed distance and inout queries
#
# Note: By default we only run a small subset of the quest regression tests.
#       To run additional regression tests, configure the code
#       with AXOM_QUEST_ENABLE_EXTRA_REGRESSION_TESTS set to ON.
#------------------------------------------------------------------------------

cmake_dependent_option(
    AXOM_QUEST_ENABLE_EXTRA_REGRESSION_TESTS
    "Enable extra regression tests for quest" OFF
    "AXOM_ENABLE_MPI;AXOM_ENABLE_SIDRE;AXOM_DATA_DIR;HDF5_FOUND" OFF )
mark_as_advanced(AXOM_QUEST_ENABLE_EXTRA_REGRESSION_TESTS)

if (AXOM_ENABLE_MPI AND AXOM_ENABLE_SIDRE AND HDF5_FOUND)

    set(quest_regression_depends
        ${quest_tests_depends}
        cli11
        )

    blt_list_append(TO quest_regression_depends ELEMENTS hdf5 IF HDF5_FOUND)
    blt_list_append(TO quest_regression_depends ELEMENTS scr IF SCR_FOUND)

    set(test_name quest_regression_test)

    axom_add_executable(
        NAME       ${test_name}
        SOURCES    quest_regression.cpp
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON ${quest_regression_depends}
        FOLDER     axom/quest/tests
        )

    if ( NOT AXOM_ENABLE_OPENMP )
        set(extra_compile_flags "${AXOM_DISABLE_OMP_PRAGMA_WARNINGS}")

        if ( AXOM_ENABLE_CUDA )
            string(STRIP "${extra_compile_flags}" extra_compile_flags)
            set(extra_compile_flags "-Xcompiler=${extra_compile_flags}")
        endif()

        blt_add_target_compile_flags( TO ${test_name} FLAGS "${extra_compile_flags}" )

    endif()

    # Add resolution tests for each dataset and resolution
    if(AXOM_DATA_DIR)
        set(quest_data_dir  ${AXOM_DATA_DIR}/quest)

        # Datasets contain a mesh name and resolution delimited by a colon
        set(quest_regression_datasets sphere:16
                                      tetrahedron:31)

        # Add additional regression datasets when its config variable is defined
        blt_list_append(
            TO          quest_regression_datasets
            ELEMENTS    aatbase_3_binary:16
                        aatbase_3_binary:32
                        boxedSphere:40
                        cap:16
                        cap:32
                        fin0:32
                        fin1:32
                        fin2:32
                        naca0012:16
                        naca0012:32
                        plane_binary:16
                        plane_binary:32
                        plane_simp:16
                        plane_simp:32
                        plane_simp:120_120_60
                        sphere:32
                        tetrahedron:100
            IF          AXOM_QUEST_ENABLE_EXTRA_REGRESSION_TESTS
           )

        foreach(el ${quest_regression_datasets})
            # Extract mesh and resolution from the dataset string
            string(REPLACE ":" ";" _dataset_res ${el})
            list(GET _dataset_res 0 _dataset)
            list(GET _dataset_res 1 _res)

            axom_add_test(
                NAME quest_regression_${_dataset}_${_res}
                COMMAND ${test_name}
                        --mesh ${quest_data_dir}/${_dataset}.stl
                        --baseline ${quest_data_dir}/regression/${_dataset}_${_res}_baseline.root
                NUM_MPI_TASKS 1
            )
        endforeach()
    endif()
endif()
