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

function retry_build_verbose () {
    local build_cmd="$1"
    local verbose_build_cmd="$2"

    if $build_cmd ; then
        return 0
    fi

    echo
    echo "=================================================================="
    echo " Non-verbose build failed. Re-running with verbose output."
    echo "=================================================================="
    echo

    or_die $verbose_build_cmd
}

echo "~~~~ helpful info ~~~~"
echo "USER="`id -u -n`
echo "PWD="`pwd`
echo "HOST_CONFIG=$HOST_CONFIG"
echo "CMAKE_EXTRA_FLAGS=$CMAKE_EXTRA_FLAGS"
echo "USE_NINJA=$USE_NINJA"
echo "~~~~~~~~~~~~~~~~~~~~~~"

export BUILD_TYPE=${BUILD_TYPE:-Debug}
export USE_NINJA=${USE_NINJA:-no}


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
        BUILD_CMD="ninja -j $NUM_BUILD_PROCS"
        VERBOSE_BUILD_CMD="ninja -j $NUM_BUILD_PROCS -v"
    fi

    echo "~~~~~~ RUNNING CMAKE ~~~~~~~~"
    or_die python3 ./config-build.py -bp builddir -hc ./host-configs/docker/${HOST_CONFIG} -bt ${BUILD_TYPE} ${BUILD_GENERATOR_FLAG} -DENABLE_GTEST_DEATH_TESTS=ON ${CMAKE_EXTRA_FLAGS}
    or_die cd builddir

    echo "~~~~~~ BUILDING ~~~~~~~~"
    retry_build_verbose "${BUILD_CMD}" "${VERBOSE_BUILD_CMD}"

    echo "~~~~~~ RUNNING TESTS ~~~~~~~~"
    or_die ctest --output-on-failure -T Test -VV -j $NUM_BUILD_PROCS

    if [[ "${DO_BENCHMARKS}" == "yes" ]] ; then
        echo "~~~~~~ RUNNING BENCHMARKS ~~~~~~~~"
        or_die ${BUILD_TOOL} -j $NUM_BUILD_PROCS run_benchmarks
    fi

    if [[ "${DO_MEMCHECK}" == "yes" ]] ; then
        echo "~~~~~~ RUNNING MEMCHECK ~~~~~~~~"
        or_die ctest -T memcheck
    fi
fi
