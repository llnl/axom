# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

axom_add_executable(
    NAME       spin_bvh_benchmark
    SOURCES    spin_bvh_benchmark.cpp
    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
    DEPENDS_ON spin slic gbenchmark
    FOLDER     axom/spin/benchmarks)

set(_spin_benchmark_num_threads)
if(AXOM_ENABLE_OPENMP)
    set(_spin_benchmark_num_threads ${AXOM_TEST_NUM_OMP_THREADS})
endif()

blt_add_benchmark(
    NAME            spin_bvh_benchmark
    COMMAND         spin_bvh_benchmark --benchmark_min_time=0.0001s
    NUM_OMP_THREADS ${_spin_benchmark_num_threads})
