#!/bin/bash
##############################################################################
# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
##############################################################################

set -x

function or_die () {
    "$@"
    local status=$?
    if [[ $status != 0 ]] ; then
        echo ERROR $status command: $@
        exit $status
    fi
}

function format_duration () {
    local total_seconds="$1"
    local hours=$(( total_seconds / 3600 ))
    local minutes=$(( (total_seconds % 3600) / 60 ))
    local seconds=$(( total_seconds % 60 ))
    printf "%d:%02d:%02d" "$hours" "$minutes" "$seconds"
}

function print_timing_line () {
    local label="$1"
    local seconds="$2"
    if [[ -n "$seconds" ]] ; then
        printf " %-22s %s\n" "$label" "$(format_duration "$seconds")"
    fi
}

function print_timing_summary () {
    local status="$1"
    local total_seconds=$(( $(date +%s) - TOTAL_START_SECONDS ))

    # Disable shell tracing so the summary is readable in CI logs.
    set +x

    echo
    echo "=================================================================="
    echo " Stage Timing Summary"
    echo "=================================================================="
    print_timing_line "configure" "$CONFIGURE_SECONDS"
    print_timing_line "build" "$BUILD_SECONDS"
    print_timing_line "test" "$TEST_SECONDS"
    print_timing_line "benchmarks" "$BENCHMARKS_SECONDS"
    print_timing_line "memcheck" "$MEMCHECK_SECONDS"
    print_timing_line "total" "$total_seconds"
    echo " result                 exit code ${status}"
    echo "=================================================================="
}

function retry_build_verbose () {
    local build_cmd="$1"
    local verbose_build_cmd="$2"

    if eval "$build_cmd" ; then
        return 0
    fi

    echo
    echo "=================================================================="
    echo " Non-verbose build failed. Re-running with verbose output."
    echo "=================================================================="
    echo

    eval "$verbose_build_cmd"
}

echo "~~~~ helpful info ~~~~"
echo "USER="`id -u -n`
echo "PWD="`pwd`
echo "HOST_CONFIG=$HOST_CONFIG"
echo "CMAKE_EXTRA_FLAGS=$CMAKE_EXTRA_FLAGS"
echo "USE_NINJA=$USE_NINJA"
echo "USE_LLD=$USE_LLD"
echo "~~~~~~~~~~~~~~~~~~~~~~"

export BUILD_TYPE=${BUILD_TYPE:-Debug}
export USE_NINJA=${USE_NINJA:-no}
export USE_LLD=${USE_LLD:-no}
TOTAL_START_SECONDS=$(date +%s)
CONFIGURE_SECONDS=""
BUILD_SECONDS=""
TEST_SECONDS=""
BENCHMARKS_SECONDS=""
MEMCHECK_SECONDS=""
trap 'status=$?; set +x; print_timing_summary "$status"' EXIT


if [[ "$DO_BUILD" == "yes" ]] ; then
    echo "~~~~~~ FIND NUMPROCS ~~~~~~~~"
    NUMPROCS=`python3 -c "import os; print(f'{os.cpu_count()}')"`
    NUM_BUILD_PROCS=`python3 -c "import os; print(f'{max(2, os.cpu_count() * 8 // 10)}')"`
    BUILD_GENERATOR_FLAG=""
    BUILD_TOOL="make"
    BUILD_CMD="make -j $NUM_BUILD_PROCS"
    VERBOSE_BUILD_CMD="make -j $NUM_BUILD_PROCS VERBOSE=1"
    if [[ "$USE_NINJA" == "yes" ]] ; then
        BUILD_GENERATOR_FLAG="--ninja"
        BUILD_TOOL="ninja"
        BUILD_CMD="ninja"
        VERBOSE_BUILD_CMD="ninja -v"
    fi
    if [[ "$USE_LLD" == "yes" ]] ; then
        CMAKE_EXTRA_FLAGS="${CMAKE_EXTRA_FLAGS} ${LLD_CMAKE_FLAGS}"
    fi

    echo "~~~~~~ RUNNING CMAKE ~~~~~~~~"
    STAGE_START_SECONDS=$(date +%s)
    python3 ./config-build.py -bp builddir -hc ./host-configs/docker/${HOST_CONFIG} -bt ${BUILD_TYPE} ${BUILD_GENERATOR_FLAG} -DENABLE_GTEST_DEATH_TESTS=ON ${CMAKE_EXTRA_FLAGS}
    STATUS=$?
    CONFIGURE_SECONDS=$(( $(date +%s) - STAGE_START_SECONDS ))
    if [[ $STATUS != 0 ]] ; then
        echo ERROR $STATUS command: python3 ./config-build.py
        exit $STATUS
    fi
    or_die cd builddir

    echo "~~~~~~ BUILDING ~~~~~~~~"
    STAGE_START_SECONDS=$(date +%s)
    retry_build_verbose "${BUILD_CMD}" "${VERBOSE_BUILD_CMD}"
    STATUS=$?
    BUILD_SECONDS=$(( $(date +%s) - STAGE_START_SECONDS ))
    if [[ $STATUS != 0 ]] ; then
        echo ERROR $STATUS command: ${BUILD_CMD}
        exit $STATUS
    fi

    echo "~~~~~~ RUNNING TESTS ~~~~~~~~"
    STAGE_START_SECONDS=$(date +%s)
    ctest --output-on-failure -T Test -VV -j $NUM_BUILD_PROCS
    STATUS=$?
    TEST_SECONDS=$(( $(date +%s) - STAGE_START_SECONDS ))
    if [[ $STATUS != 0 ]] ; then
        echo ERROR $STATUS command: ctest --output-on-failure -T Test -VV -j $NUM_BUILD_PROCS
        exit $STATUS
    fi

    if [[ "${DO_BENCHMARKS}" == "yes" ]] ; then
        echo "~~~~~~ RUNNING BENCHMARKS ~~~~~~~~"
        STAGE_START_SECONDS=$(date +%s)
        ${BUILD_TOOL} -j $NUM_BUILD_PROCS run_benchmarks
        STATUS=$?
        BENCHMARKS_SECONDS=$(( $(date +%s) - STAGE_START_SECONDS ))
        if [[ $STATUS != 0 ]] ; then
            echo ERROR $STATUS command: ${BUILD_TOOL} -j $NUM_BUILD_PROCS run_benchmarks
            exit $STATUS
        fi
    fi

    if [[ "${DO_MEMCHECK}" == "yes" ]] ; then
        echo "~~~~~~ RUNNING MEMCHECK ~~~~~~~~"
        STAGE_START_SECONDS=$(date +%s)
        ctest -T memcheck
        STATUS=$?
        MEMCHECK_SECONDS=$(( $(date +%s) - STAGE_START_SECONDS ))
        if [[ $STATUS != 0 ]] ; then
            echo ERROR $STATUS command: ctest -T memcheck
            exit $STATUS
        fi
    fi
fi
