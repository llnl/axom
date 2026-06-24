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
endif()
