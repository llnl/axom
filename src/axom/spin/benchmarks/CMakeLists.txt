# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Spin benchmarks
#------------------------------------------------------------------------------

axom_add_executable(
    NAME       spin_sfc_codec_benchmark
    SOURCES    spin_sfc_codec.cpp
    HEADERS    ../tests/hilbert_reference.hpp
    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
    DEPENDS_ON spin slic fmt cli11 gbenchmark
    FOLDER     axom/spin/benchmarks
    )

blt_add_benchmark(
    NAME    spin_sfc_codec_benchmark
    COMMAND spin_sfc_codec_benchmark --benchmark_min_time=0.0001s
    )

if (AXOM_ENABLE_TESTS)
    axom_add_test(
        NAME    spin_sfc_codec_candidates
        COMMAND spin_sfc_codec_benchmark --verify-only
        )
endif()
