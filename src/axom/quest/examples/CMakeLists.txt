# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Quest examples
#------------------------------------------------------------------------------

set(quest_example_depends quest fmt cli11)

blt_list_append(TO quest_example_depends ELEMENTS RAJA IF RAJA_FOUND)

set(quest_data_dir  ${AXOM_DATA_DIR}/quest)
set(shaping_data_dir  ${AXOM_DATA_DIR}/shaping)

# In/out octree containment example -------------------------------------------
axom_add_executable(
    NAME        quest_containment_driver_ex
    SOURCES     containment_driver.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends}
    FOLDER      axom/quest/examples
    )

if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR AND CALIPER_FOUND)
    set(_input_file "${AXOM_DATA_DIR}/quest/sphere_binary.stl")
    
    # this test requires MPI when available (but only 1 rank)
    axom_add_test(
        NAME    quest_containment_sphere_ex
        COMMAND quest_containment_driver_ex
                --input ${_input_file}
                --caliper report
        NUM_MPI_TASKS 1
    )
endif()

# BVH two pass example --------------------------------------------------------
axom_add_executable(
    NAME        quest_bvh_two_pass_ex
    SOURCES     quest_bvh_two_pass.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends}
    FOLDER      axom/quest/examples
    )

# Add unit tests
if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR)

    # Run the quest_bvh_two_pass example with different raja policies

    set(input_file "${AXOM_DATA_DIR}/quest/unit_cube.stl")

    set (_policies "seq")
    if(RAJA_FOUND)
        blt_list_append(TO _policies ELEMENTS "omp" IF AXOM_ENABLE_OPENMP)
        blt_list_append(TO _policies ELEMENTS "cuda" IF AXOM_ENABLE_CUDA)
        blt_list_append(TO _policies ELEMENTS "hip" IF AXOM_ENABLE_HIP)
    endif()

    foreach(_policy ${_policies})
        # sets the number of threads for the omp policy; leave empty for non-omp policies
        set(_num_threads)
        if(_policy STREQUAL "omp")
            set(_num_threads ${AXOM_TEST_NUM_OMP_THREADS})
        endif()

        set(_testname "quest_bvh_two_pass_example_${_policy}")
        axom_add_test(
          NAME ${_testname}
          COMMAND quest_bvh_two_pass_ex
                  --file ${input_file}
                  --exec_space ${_policy}
          NUM_OMP_THREADS ${_num_threads})

        set_tests_properties(${_testname} PROPERTIES
            PASS_REGULAR_EXPRESSION  "Found 4 actual collisions")
    endforeach()
endif()

# Read ProE example -----------------------------------------------------------
axom_add_executable(
    NAME        quest_proe_bbox_ex
    SOURCES     quest_proe_bbox.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends}
    FOLDER      axom/quest/examples
    )

# BVH silo example ------------------------------------------------------------
if (CONDUIT_FOUND AND UMPIRE_FOUND)
    axom_add_executable(
        NAME        quest_candidates_example_ex
        SOURCES     quest_candidates_example.cpp
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${quest_example_depends} conduit::conduit
        FOLDER      axom/quest/examples
        )

    # Add unit tests; the (current) input meshes are in the hdf5 format
    if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR AND HDF5_FOUND)

        # Run the candidates example with the different spatial indices and raja policies

        # Use same file for input and query
        set(input_file "${AXOM_DATA_DIR}/quest/ucart10.cycle_000000.root")
        set(query_file "${AXOM_DATA_DIR}/quest/ucart10_shifted.cycle_000000.root")

        set(_methods "bvh" "implicit")

        set (_policies "seq")
        if(RAJA_FOUND)
          blt_list_append(TO _policies ELEMENTS "omp" IF AXOM_ENABLE_OPENMP)
          blt_list_append(TO _policies ELEMENTS "cuda" IF AXOM_ENABLE_CUDA)
          blt_list_append(TO _policies ELEMENTS "hip" IF AXOM_ENABLE_HIP)
        endif()

        foreach(_method ${_methods})
            foreach(_policy ${_policies})
                # sets the number of threads for the omp policy; leave empty for non-omp policies
                set(_num_threads)
                if(_policy STREQUAL "omp")
                    set(_num_threads ${AXOM_TEST_NUM_OMP_THREADS})
                endif()

                set(_testname "quest_candidates_example_${_method}_${_policy}")
                axom_add_test(
                  NAME ${_testname}
                  COMMAND quest_candidates_example_ex
                          --infile ${input_file}
                          --queryfile ${query_file}
                          --method ${_method}
                          --policy ${_policy}
                  NUM_OMP_THREADS ${_num_threads})

                # Match either one comma or none for portability
                set_tests_properties(${_testname} PROPERTIES
                    PASS_REGULAR_EXPRESSION  "Mesh had 6[,]?859 candidates pairs")
            endforeach()
        endforeach()
    endif()

endif()

# Shaping example -------------------------------------------------------------
if(AXOM_ENABLE_MPI AND MFEM_FOUND AND MFEM_USE_MPI
                   AND AXOM_ENABLE_SIDRE AND AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION
                   AND AXOM_ENABLE_KLEE)
    axom_add_executable(
        NAME        quest_shaping_driver_ex
        SOURCES     shaping_driver.cpp
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${quest_example_depends} mfem
        FOLDER      axom/quest/examples
        )

    if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR)
        # 2D shaping tests
        set(_nranks 1)
        if(C2C_FOUND)
            set(_testname quest_shaping_driver_ex_sampling_circles)
            axom_add_test(
                NAME ${_testname}
                COMMAND  quest_shaping_driver_ex 
                            -i ${shaping_data_dir}/circles.yaml
                            --method sampling
                            inline_mesh --min -6 -6 --max 6 6 --res 25 25 -d 2
                NUM_MPI_TASKS ${_nranks})
            # Analytic area for annulus w/ outer/inner radii 5 and 2.5 is ~58.905
            set_tests_properties(${_testname} PROPERTIES 
                PASS_REGULAR_EXPRESSION  "Volume of material 'steel' is 58.")

            set(_testname quest_shaping_driver_ex_sampling_circles_opposite)
            axom_add_test(
                NAME ${_testname}
                COMMAND  quest_shaping_driver_ex 
                        -i ${shaping_data_dir}/circles_opposite.yaml
                        --method sampling
                        inline_mesh --min -6 -6 --max 6 6 --res 25 25 -d 2
                NUM_MPI_TASKS ${_nranks})
            # Analytic area for annulus w/ outer/inner radii 5 and 2.5 is ~65.97
            set_tests_properties(${_testname} PROPERTIES 
                PASS_REGULAR_EXPRESSION  "Volume of material 'steel' is 65.")

            set(_testname quest_shaping_driver_ex_sampling_balls_and_jacks)
            axom_add_test(
                NAME ${_testname}
                COMMAND  quest_shaping_driver_ex 
                        -i ${shaping_data_dir}/balls_and_jacks.yaml
                        --method sampling
                        --background-material air
                        inline_mesh --min -10 0 --max 60 50 --res 25 25 -d 2
                NUM_MPI_TASKS ${_nranks})
            # background: 3500; balls: ~373; jacks: ~230; air ~2897
            set_tests_properties(${_testname} PROPERTIES 
                PASS_REGULAR_EXPRESSION  "Volume of material 'air' is 2,?895.")

            set(_testname quest_shaping_driver_ex_sampling_ball_impact)
            axom_add_test(
                NAME ${_testname}
                COMMAND  quest_shaping_driver_ex 
                        -i ${shaping_data_dir}/ball_impact.yaml
                        --method sampling
                        --background-material air
                        inline_mesh --min -10 0 --max 40 30 --res 25 25 -d 2
                NUM_MPI_TASKS ${_nranks})
            # steel ball: semi-cirle of radius 5.5 has area ~47.516
            set_tests_properties(${_testname} PROPERTIES 
                PASS_REGULAR_EXPRESSION  "Volume of material 'steel' is 47.")

            # Tests 3D mesh against 2D contours
            set(_testname quest_shaping_driver_ex_sampling_ball_impact_3d)
            axom_add_test(
                NAME ${_testname}
                COMMAND  quest_shaping_driver_ex 
                        -i ${shaping_data_dir}/ball_impact.yaml
                        --method sampling
                        --background-material air
                        inline_mesh --min -30 -30 -10 --max 30 30 40 --res 16 16 16 -d3
                NUM_MPI_TASKS ${_nranks})
            # steel ball: cirle of radius 5.5 has analytic volume ~696.91
            set_tests_properties(${_testname} PROPERTIES 
                PASS_REGULAR_EXPRESSION  "Volume of material 'steel' is 6[89][0-9].")
        endif()

        # Tests 2D stl mesh with intersection-based shaper
        set(_testname quest_shaping_driver_ex_star_stl_2d)
        axom_add_test(
            NAME ${_testname}
            COMMAND  quest_shaping_driver_ex
                    -i ${shaping_data_dir}/flat_star.yaml
                    --method intersection
                    inline_mesh --min 0 0 --max 1 1 --res 4 4 -d 2
            NUM_MPI_TASKS ${_nranks})
        # steel triangle: star composed of 16 right triangles with
        # leg length 0.25 has analytic volume 0.5
        set_tests_properties(${_testname} PROPERTIES
            PASS_REGULAR_EXPRESSION  "Volume of material 'steel' is 0.5")

        # Tests 2D c2c contour
        if(C2C_FOUND)
            set(_policies "seq")
            if(RAJA_FOUND)
                blt_list_append(TO _policies ELEMENTS "omp" IF AXOM_ENABLE_OPENMP)
                if(UMPIRE_FOUND)
                  blt_list_append(TO _policies ELEMENTS "cuda" IF AXOM_ENABLE_CUDA)
                  blt_list_append(TO _policies ELEMENTS "hip" IF AXOM_ENABLE_HIP)
                endif()
            endif()

            foreach(_policy ${_policies})
                # sets the number of threads for the omp policy; leave empty for non-omp policies
                set(_num_threads)
                if(_policy STREQUAL "omp")
                    set(_num_threads ${AXOM_TEST_NUM_OMP_THREADS})
                endif()

                set(_testname quest_shaping_driver_ex_semi_circle_2d_${_policy})

                axom_add_test(
                    NAME ${_testname}
                    COMMAND  quest_shaping_driver_ex
                            -i ${shaping_data_dir}/revolved_circle.yaml
                            --policy ${_policy}
                            --segments-per-knot-span 50
                            --method intersection
                            inline_mesh --min 0 0 --max 3 3 --res 3 3 -d 2
                    NUM_MPI_TASKS ${_nranks}
                    NUM_OMP_THREADS ${_num_threads})
                # steel unit semi-circle with an area of pi/2
                set_tests_properties(${_testname} PROPERTIES
                    PASS_REGULAR_EXPRESSION  "Volume of material 'steel' is 1.57")
            endforeach()

            unset(_policies)
        endif()

        # Test 2D MFEM mesh (linearize contours)
        set(_testname quest_shaping_driver_ex_heroic_roses_mfem)
        axom_add_test(
            NAME ${_testname}
            COMMAND  quest_shaping_driver_ex 
                     -i ${shaping_data_dir}/heroic_roses_mfem.yaml
                     --method sampling
                     inline_mesh --min 0 0 --max 300 400 --res 150 200 -d 2
            NUM_MPI_TASKS ${_nranks})
        set_tests_properties(${_testname} PROPERTIES
            PASS_REGULAR_EXPRESSION  "Volume of material 'black' is 36,?595.")

        # Test 2D MFEM mesh (winding number)
        set(_testname quest_shaping_driver_ex_heroic_roses_mfem_cp)
        axom_add_test(
            NAME ${_testname}
            COMMAND  quest_shaping_driver_ex 
                     -i ${shaping_data_dir}/heroic_roses_mfem_cp.yaml
                     --method sampling --sampling windingnumber
                     inline_mesh --min 0 0 --max 300 400 --res 150 200 -d 2
            NUM_MPI_TASKS ${_nranks})
        set_tests_properties(${_testname} PROPERTIES
            PASS_REGULAR_EXPRESSION  "Volume of material 'black' is 36,?654.")

        # 3D shaping tests
        set(_nranks 3)
        set(_testname quest_shaping_driver_ex_sampling_sphere)
        axom_add_test(
            NAME ${_testname}
            COMMAND  quest_shaping_driver_ex 
                    -i ${shaping_data_dir}/sphere.yaml
                    --method sampling
                    inline_mesh --min -6 -6 -6 --max 6 6 6 --res 16 16 16 -d 3
            NUM_MPI_TASKS ${_nranks})
        # sphere of radius 5 has volume: 4/3 pi * r^3 ~523.6
        # input sphere is discretized and mesh is coarse; accuracy is within ~1%
        set_tests_properties(${_testname} PROPERTIES 
            PASS_REGULAR_EXPRESSION  "Volume of material 'steel' is 519.")

        set(_testname quest_shaping_driver_ex_sampling_spheres)
        axom_add_test(
            NAME ${_testname}
            COMMAND  quest_shaping_driver_ex 
                    -i ${shaping_data_dir}/spheres.yaml
                    --method sampling
                    --background-material void
                    inline_mesh --min -6 -6 -6 --max 6 6 6 --res 16 16 16 -d 3
            NUM_MPI_TASKS ${_nranks})
        # bbox volume: 12^3 = 1728; sphere(r=5): ~523.6; sphere(r=2): 33.5
        # expected analytic volume when fully resolved: ~1237.9
        set_tests_properties(${_testname} PROPERTIES 
            PASS_REGULAR_EXPRESSION  "Volume of material 'void' is 1,?239.")

        set(_testname quest_shaping_driver_ex_sampling_plane)
        axom_add_test(
            NAME ${_testname}
            COMMAND  quest_shaping_driver_ex 
                    -i ${shaping_data_dir}/plane.yaml
                    --method sampling
                    --background-material air
                    inline_mesh --min -430 -45 -35 --max 430 855 210 --res 16 16 16 -d 3
            NUM_MPI_TASKS ${_nranks})
        # background: 860*900*245=189,630,000; airplane volume from meshlab: 7,420,578.5
        # expected volume when fully resolved: 182,209,421.5
        # actual value is accurate to about 1% at coarse resolution
        set_tests_properties(${_testname} PROPERTIES 
            PASS_REGULAR_EXPRESSION  "Volume of material 'air' is 182,?227,?963")
    endif()
endif()

# Shaping in-memory example ------------------------------------------------------
if((CONDUIT_FOUND OR
      (AXOM_ENABLE_MPI AND MFEM_FOUND AND MFEM_USE_MPI
        AND AXOM_ENABLE_SIDRE AND AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION))
    AND AXOM_ENABLE_KLEE)
    if(MFEM_FOUND)
        set(optional_dependency, "mfem")
    endif()
    axom_add_executable(
        NAME        quest_shape_in_memory_ex
        SOURCES     quest_shape_in_memory.cpp
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${quest_example_depends} ${optional_dependency}
        FOLDER      axom/quest/examples
        )

    if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR)

        # this test requires MPI when available (but only 1 rank)
        set(_nranks 1)

        # Run the in-memory shaping example on with each enabled policy
        set(_policies "seq")
        if(RAJA_FOUND)
            # Hanging issue on Windows
            if(NOT WIN32)
                blt_list_append(TO _policies ELEMENTS "omp"  IF AXOM_ENABLE_OPENMP)
            endif()
            blt_list_append(TO _policies ELEMENTS "cuda" IF AXOM_ENABLE_CUDA)
            blt_list_append(TO _policies ELEMENTS "hip"  IF AXOM_ENABLE_HIP)
        endif()

        set(_testshapes "tetmesh" "tet" "hex" "sphere" "cyl" "cone" "sor" "all" "plane")
        foreach(_policy ${_policies})
            # sets the number of threads for the omp policy; leave empty for non-omp policies
            set(_num_threads)
            if(_policy STREQUAL "omp")
                set(_num_threads ${AXOM_TEST_NUM_OMP_THREADS})
            endif()

            foreach(_testshape ${_testshapes})
              set(_testname "quest_shape_in_memory_${_policy}_${_testshape}")
              axom_add_test(
                NAME ${_testname}
                COMMAND quest_shape_in_memory_ex
                           --policy ${_policy}
                           --testShape ${_testshape}
                           --refinements 2
                           --scale 0.75 0.75 0.75
                           --dir 0.2 0.4 0.8
                           --meshType bpSidre
                           inline_mesh --min -2 -2 -2 --max 2 2 2 --res 10 10 10
                NUM_MPI_TASKS ${_nranks}
                NUM_OMP_THREADS ${_num_threads})
            endforeach()
        endforeach()

        # Test Conduit and MFEM mesh types with a subset of shapes (so it doesn't take too long)
        set(_testshapes "tetmesh" "sphere" "cyl")
        set(_testMeshTypes "bpConduit")
        blt_list_append(TO _testMeshTypes ELEMENTS "mfem" IF MFEM_FOUND)
        foreach(_policy ${_policies})
            # sets the number of threads for the omp policy; leave empty for non-omp policies
            set(_num_threads)
            if(_policy STREQUAL "omp")
                set(_num_threads ${AXOM_TEST_NUM_OMP_THREADS})
            endif()

            foreach(_testMeshType ${_testMeshTypes})
              foreach(_testshape ${_testshapes})
                set(_testname "quest_shape_in_memory_${_policy}_${_testMeshType}_${_testshape}")
                axom_add_test(
                  NAME ${_testname}
                  COMMAND quest_shape_in_memory_ex
                             --policy ${_policy}
                             --testShape ${_testshape}
                             --refinements 2
                             --scale 0.75 0.75 0.75
                             --dir 0.2 0.4 0.8
                             --meshType ${_testMeshType}
                             inline_mesh --min -2 -2 -2 --max 2 2 2 --res 8 8 8
                  NUM_MPI_TASKS ${_nranks}
                  NUM_OMP_THREADS ${_num_threads})
              endforeach()
            endforeach()
        endforeach()

        # Test 2D Triangle shape
        set(_testshapes "tri")
        set(_testMeshTypes "bpConduit" "bpSidre")
        blt_list_append(TO _testMeshTypes ELEMENTS "mfem" IF MFEM_FOUND)
        foreach(_policy ${_policies})
            foreach(_testMeshType ${_testMeshTypes})
              foreach(_testshape ${_testshapes})
                set(_testname "quest_shape_in_memory_${_policy}_${_testMeshType}_${_testshape}")
                axom_add_test(
                  NAME ${_testname}
                  COMMAND quest_shape_in_memory_ex
                             --policy ${_policy}
                             --testShape ${_testshape}
                             --refinements 2
                             --scale 0.75 0.75 0.75
                             --dir 0.2 0.4 0.8
                             --meshType ${_testMeshType}
                             inline_mesh --min -2 -2 --max 2 2 --res 2 2
                  NUM_MPI_TASKS ${_nranks}
                  NUM_OMP_THREADS ${_num_threads})
              endforeach()
            endforeach()
        endforeach()
    endif()
endif()

# Distributed closest point example -------------------------------------------
if(AXOM_ENABLE_MPI AND AXOM_ENABLE_SIDRE AND HDF5_FOUND)
    axom_add_executable(
            NAME        quest_distributed_distance_query_ex
            SOURCES     quest_distributed_distance_query_example.cpp
            OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
            DEPENDS_ON  ${quest_example_depends}
            FOLDER      axom/quest/examples
            )

    if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR)
        set(_nranks 3)

        # Run the distributed closest point example on N ranks for each enabled policy
        set(_policies "seq")
        if(RAJA_FOUND)
            blt_list_append(TO _policies ELEMENTS "cuda" IF AXOM_ENABLE_CUDA)
            blt_list_append(TO _policies ELEMENTS "hip"  IF AXOM_ENABLE_HIP)
            blt_list_append(TO _policies ELEMENTS "omp"  IF AXOM_ENABLE_OPENMP)
        endif()

        # Non-zero empty-rank probability tests domain underloading case
        set(_meshes "mdmesh.2x1" "mdmesh.2x3" "mdmesh.2x2x1" "mdmeshg.2x2x1")
        # The mdmesh.* files were generated by these commands:
        # src/tools/gen-multidom-structured-mesh.py -ml=0,0 -mu=2,2 -ms=100,100 -dc=2,1 -o mdmesh.2x1
        # src/tools/gen-multidom-structured-mesh.py -ml=0,0 -mu=2,2 -ms=100,100 -dc=2,3 -o mdmesh.2x3
        # src/tools/gen-multidom-structured-mesh.py -ml=0,0,0 -mu=2,2,2 -ms=20,20,15 -dc=2,2,1 -o mdmeshg.2x2x1 --strided
        foreach(_pol ${_policies})
            # sets the number of threads for the omp policy; leave empty for non-omp policies
            set(_num_threads)
            if(_pol STREQUAL "omp")
                set(_num_threads ${AXOM_TEST_NUM_OMP_THREADS})
            endif()

            foreach(_mesh ${_meshes})
                # Determine problem dimension by mesh filename.
                # and set dimension-dependent arguments.
                string(REGEX MATCH "\\.[0-9]+(x[0-9]+)+$" _sizes "${_mesh}")
                string(REGEX MATCHALL "[0-9]+" _sizes ${_sizes})
                list(LENGTH _sizes _ndim)

                if(_ndim EQUAL 2)
                  set(_center 0.7 0.9)
                elseif(_ndim EQUAL 3)
                  set(_center 0.7 0.9 0.5)
                endif()

                set(_test "quest_distributed_closest_point_run_${_ndim}D_${_pol}_${_mesh}")
                axom_add_test(
                    NAME    ${_test}
                    COMMAND quest_distributed_distance_query_ex
                                --mesh-file ${quest_data_dir}/${_mesh}.root
                                --long-point-count 60
                                --center ${_center}
                                --radius 0.9
                                --lat-point-count 30
                                --obj-domain-count-range 0 2
                                --dist-threshold .3
                                --no-random-spacing
                                --check-results
                                --policy ${_pol}
                                --object-file dcp_object_mesh_${_ndim}d_${_pol}_${_mesh}
                                --distance-file dcp_closest_point_2d_${_pol}_${_mesh}
                    NUM_MPI_TASKS ${_nranks}
                    NUM_OMP_THREADS ${_num_threads})
            endforeach()
        endforeach()

        unset(optional_dependency)
        unset(_nranks)
        unset(_policies)
        unset(_test)
    endif()
endif()

# Marching cubes example -------------------------------------------
if(CONDUIT_FOUND AND AXOM_ENABLE_MPI)
    list(APPEND quest_depends_on conduit::conduit
                                 conduit::conduit_mpi)
endif()

if(CONDUIT_FOUND)
    axom_add_executable(
            NAME        quest_marching_cubes_ex
            SOURCES     quest_marching_cubes_example.cpp
            OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
            DEPENDS_ON  ${quest_example_depends}
            FOLDER      axom/quest/examples
            )

    # These examples currently segfault on windows without MPI when loading the test data into conduit
    # Also, the input files are hdf5
    if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR AND NOT WIN32 AND HDF5_FOUND)
        if (ENABLE_MPI)
            set(_nranks 3)
        endif()

        # Run the marching cubes example on N ranks for each enabled policy
        set(_policies "seq")
        if(RAJA_FOUND)
            blt_list_append(TO _policies ELEMENTS "omp" IF AXOM_ENABLE_OPENMP)
            blt_list_append(TO _policies ELEMENTS "cuda" IF AXOM_ENABLE_CUDA)
            blt_list_append(TO _policies ELEMENTS "hip" IF AXOM_ENABLE_HIP)
        endif()

        # Non-zero empty-rank probability tests domain underloading case
        set(_meshes "mdmesh.2x1" "mdmesh.2x3" "mdmesh.2x2x1" "mdmeshg.2x2x1")
        # The amc.* files were generated by these commands:
        # src/tools/gen-multidom-structured-mesh.py -ml=0,0 -mu=2,2 -ms=100,100 -dc=2,1 -o mdmesh.2x1
        # src/tools/gen-multidom-structured-mesh.py -ml=0,0 -mu=2,2 -ms=100,100 -dc=2,3 -o mdmesh.2x3
        # src/tools/gen-multidom-structured-mesh.py -ml=0,0,0 -mu=2,2,2 -ms=20,20,15 -dc=2,2,1 -o mdmesh.2x2x1
        # src/tools/gen-multidom-structured-mesh.py -ml=0,0,0 -mu=2,2,2 -ms=20,20,15 -dc=2,2,1 -o mdmeshg.2x2x1 --strided
        foreach(_pol ${_policies})
            # sets the number of threads for the omp policy; leave empty for non-omp policies
            set(_num_threads)
            if(_pol STREQUAL "omp")
                set(_num_threads ${AXOM_TEST_NUM_OMP_THREADS})
            endif()

            foreach(_mesh ${_meshes})
                # Determine problem dimension by mesh filename.
                # and set dimension-dependent arguments.
                string(REGEX MATCH "\\.[0-9]+(x[0-9]+)+$" _sizes "${_mesh}")
                string(REGEX MATCHALL "[0-9]+" _sizes ${_sizes})
                list(LENGTH _sizes _ndim)

                if(_ndim EQUAL 2)
                  set(_dir 1.0 0.4)
                  set(_center 1.0 0.4)
                  set(_scale 3 3)
                elseif(_ndim EQUAL 3)
                  set(_dir 1.0 0.4 1.2)
                  set(_center 1.0 0.4 1.2)
                  set(_scale 3 3 1.5)
                endif()

                set(_test "quest_marching_cubes_run_${_ndim}D_${_pol}_${_mesh}")
                axom_add_test(
                    NAME    ${_test}
                    COMMAND quest_marching_cubes_ex
                                --policy ${_pol}
                                --mesh-file ${quest_data_dir}/${_mesh}.root
                                --fields-file ${_test}.field
                                --dir ${_dir}
                                --center ${_center}
                                --scale ${_scale}
                                --contourVal 1.25
                                --check-results
                                --verbose
                    NUM_MPI_TASKS ${_nranks}
                    NUM_OMP_THREADS ${_num_threads})
            endforeach()
        endforeach()

        unset(_nranks)
        unset(_policies)
        unset(_test)
    endif()

endif()

# Point in cell example -------------------------------------------------------
if(MFEM_FOUND)
    axom_add_executable(
        NAME        quest_point_in_cell_benchmark_ex
        SOURCES     point_in_cell_benchmark.cpp
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${quest_example_depends} mfem
        FOLDER      axom/quest/examples
        )
endif()

# Delaunay triangulation example ----------------------------------------------
axom_add_executable(
    NAME        quest_delaunay_triangulation_ex
    SOURCES     delaunay_triangulation.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends}
    FOLDER      axom/quest/examples
    )

if(AXOM_ENABLE_TESTS)
    foreach(d 2 3)
        axom_add_test(
            NAME quest_delaunay_${d}_test
            COMMAND quest_delaunay_triangulation_ex -d ${d} -o delaunay_mesh_${d}d -n 100
            )
    endforeach()
endif()

# Scattered interpolation example --------------------------------------------
if(AXOM_ENABLE_SIDRE)
    axom_add_executable(
        NAME        quest_scattered_interpolation_ex
        SOURCES     scattered_interpolation.cpp
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${quest_example_depends}
        FOLDER      axom/quest/examples
        )

    if(AXOM_ENABLE_TESTS)
        foreach(d 2 3)
            axom_add_test(
                NAME quest_scattered_interpolation_${d}_test
                COMMAND quest_scattered_interpolation_ex -d ${d} -n 10000 -q 20000
                )
        endforeach()
    endif()
endif()


# Quest signed distance and inout interface examples (C++ and Fortran) --------

axom_add_executable(
    NAME       quest_signed_distance_interface_ex
    SOURCES    quest_signed_distance_interface.cpp
    OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON ${quest_example_depends}
    FOLDER      axom/quest/examples
    )

axom_add_executable(
    NAME        quest_inout_interface_ex
    SOURCES     quest_inout_interface.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends}
    FOLDER      axom/quest/examples
    )

if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR)
    set(_num_ranks 1)
    if(AXOM_ENABLE_MPI)
        set(_num_ranks 2)
    endif()

    axom_add_test(
        NAME quest_inout_interface_3D_test
        COMMAND quest_inout_interface_ex -i ${quest_data_dir}/sphere_binary.stl
        NUM_MPI_TASKS ${_num_ranks}
        )
    if(C2C_FOUND)
        axom_add_test(
            NAME quest_inout_interface_2D_test
            COMMAND quest_inout_interface_ex -i ${AXOM_DATA_DIR}/contours/unit_circle.contour
            NUM_MPI_TASKS ${_num_ranks}
            )
    endif()
endif()

if (ENABLE_FORTRAN)

    set(quest_fortran_examples quest_signed_distance_interface)
    
    if("${AXOM_CONFIG_NAME}" STREQUAL "GNU-Debug-shared-mpi-openmp")
        # The gcc-11 docker image does not work, turning off for now
        set(_has_inout_fortran_example FALSE)
    elseif("${CMAKE_Fortran_COMPILER}" MATCHES "amdflang$" AND
       CMAKE_BUILD_TYPE MATCHES "(Debug|RelWithDebInfo)")
        # The inout Fortran example fails to compile with hipcc/amdflang in debug configurations
        set(_has_inout_fortran_example FALSE)
    else()
        set(_has_inout_fortran_example TRUE)
        blt_list_append(TO quest_fortran_examples ELEMENTS quest_inout_interface)
    endif()

    foreach(_t ${quest_fortran_examples})
        set(_example_name "${_t}_F_ex")
        axom_add_executable(
            NAME       ${_example_name}
            SOURCES    ${_t}.F
            OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
            DEPENDS_ON ${quest_example_depends}
            FOLDER     axom/quest/examples)

        # When CUDA is enabled, BLT will determine the correct linker, so don't override it here
        if(NOT AXOM_ENABLE_CUDA)
            if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "IntelLLVM")
                set_target_properties(${_example_name}
                                      PROPERTIES LINKER_LANGUAGE Fortran)
            else()
                set_target_properties(${_example_name}
                                      PROPERTIES LINKER_LANGUAGE CXX)
            endif()
        endif()
    endforeach()

    if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR AND ${_has_inout_fortran_example})
        set(_num_ranks 1)
        if(AXOM_ENABLE_MPI)
            set(_num_ranks 2)
        endif()
        axom_add_test(
            NAME quest_inout_interface_F_test
            COMMAND quest_inout_interface_F_ex ${quest_data_dir}/sphere_binary.stl
            NUM_MPI_TASKS ${_num_ranks})
    endif()
endif()

# Quest winding number example ------------------------------------------------

if( MFEM_FOUND)
    axom_add_executable(
        NAME        quest_winding_number_ex
        SOURCES     quest_winding_number.cpp
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  axom mfem cli11 fmt
        FOLDER      axom/quest/examples
        )
endif()

if(MFEM_FOUND AND OPENCASCADE_FOUND)
    axom_add_executable(
        NAME        quest_winding_number_3d_ex
        SOURCES     quest_winding_number_3d.cpp
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  axom mfem cli11 fmt opencascade
        FOLDER      axom/quest/examples
        )

    if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR)
        axom_add_test(
            NAME quest_winding_number_3d_tet_slice
            COMMAND quest_winding_number_3d_ex
                    --input ${quest_data_dir}/step/tet.step
                    --output-prefix wn3d_tet_slice
                    query_mesh --min -0.5 -0.5 --max 0.5 0.5 --res 3 3 --slice-z 0.0 --order 1)
        set_tests_properties(quest_winding_number_3d_tet_slice PROPERTIES
            PASS_REGULAR_EXPRESSION  "4 patches .* 12 trimming curves")
    endif()
endif()


if(OPENCASCADE_FOUND)
    axom_add_executable(
        NAME        quest_step_file_ex
        SOURCES     quest_step_file.cpp
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  axom cli11 fmt opencascade
        FOLDER      axom/quest/examples
    )

    if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR)
        set(_num_ranks 1)
        if(AXOM_ENABLE_MPI)
            set(_num_ranks 2)
        endif()

        axom_add_test(
            NAME quest_step_file_sliced_cylinder
            COMMAND quest_step_file_ex  -f ${quest_data_dir}/step/sliced_cylinder.step
                                        -o step_sliced_cylinder
                                        --output-trimmed vtk
                                        --output-untrimmed stl
                                        --output-svg
                                        --deflection 1e-2
                                        --angular-deflection .1
            NUM_MPI_TASKS ${_num_ranks})
        set_tests_properties(quest_step_file_sliced_cylinder PROPERTIES
            PASS_REGULAR_EXPRESSION  "5 patches .* 18 trimming curves")

        axom_add_test(
            NAME quest_step_file_tet
            COMMAND quest_step_file_ex  -f ${quest_data_dir}/step/tet.step -o step_tet
            NUM_MPI_TASKS ${_num_ranks})
        set_tests_properties(quest_step_file_tet PROPERTIES
            PASS_REGULAR_EXPRESSION  "4 patches .* 12 trimming curves")

        axom_add_test(
            NAME quest_step_file_boxed_sphere
            COMMAND quest_step_file_ex  -f ${quest_data_dir}/step/boxed_sphere.step -o step_boxed_sphere
            NUM_MPI_TASKS ${_num_ranks})
        set_tests_properties(quest_step_file_boxed_sphere PROPERTIES
            PASS_REGULAR_EXPRESSION  "192 patches .* 968 trimming curves")

    endif()
endif()
