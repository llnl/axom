# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Quest benchmarks
#------------------------------------------------------------------------------

set(quest_benchmark_files
    quest_delaunay_circumsphere.cpp
    )

if (ENABLE_BENCHMARKS)
    foreach(test ${quest_benchmark_files})
        get_filename_component(test_name ${test} NAME_WE)
        set(test_name "${test_name}_benchmark")

        axom_add_executable(
            NAME       ${test_name}
            SOURCES    ${test}
            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
            DEPENDS_ON quest gbenchmark
            FOLDER     axom/quest/benchmarks
            )

        blt_add_benchmark(
            NAME    ${test_name}
            COMMAND ${test_name} --benchmark_min_time=0.0001s
            )
    endforeach()

    if (AXOM_ENABLE_KLEE AND AXOM_ENABLE_SIDRE AND CONDUIT_FOUND
        AND RAJA_FOUND AND AXOM_ENABLE_BUMP)
        set(_mesh_clipper_benchmark_depends quest conduit::conduit RAJA)
        blt_list_append(TO _mesh_clipper_benchmark_depends
                        IF AXOM_ENABLE_MPI
                        ELEMENTS conduit::conduit_mpi)

        set(_mesh_clipper_benchmark_num_threads)
        if (AXOM_ENABLE_OPENMP)
            set(_mesh_clipper_benchmark_num_threads ${AXOM_TEST_NUM_OMP_THREADS})
        endif()

        axom_add_executable(
            NAME        quest_mesh_clipper_benchmark
            SOURCES     quest_mesh_clipper_benchmark.cpp
            OUTPUT_DIR  ${TEST_OUTPUT_DIRECTORY}
            DEPENDS_ON  ${_mesh_clipper_benchmark_depends} gbenchmark
            FOLDER      axom/quest/benchmarks
            )

        blt_add_benchmark(
            NAME            quest_mesh_clipper_benchmark
            COMMAND         quest_mesh_clipper_benchmark --benchmark_min_time=0.0001s
            NUM_OMP_THREADS ${_mesh_clipper_benchmark_num_threads}
            )
    endif()
endif()
