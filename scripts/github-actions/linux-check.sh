#!/bin/bash
##############################################################################
# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
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

echo "~~~~ helpful info ~~~~"
echo "USER="`id -u -n`
echo "PWD="`pwd`
echo "HOST_CONFIG=$HOST_CONFIG"
echo "CMAKE_EXTRA_FLAGS=$CMAKE_EXTRA_FLAGS"
echo "~~~~~~~~~~~~~~~~~~~~~~"

echo "~~~~~~ RUNNING CMAKE ~~~~~~~~"
cmake_args="-DCMAKE_BUILD_TYPE=Debug -DBUILD_SHARED_LIBS=ON -DENABLE_CLANGTIDY=OFF"

if [[ "$CHECK_TYPE" == "style" ]] ; then
    CLANGFORMAT_EXECUTABLE=/usr/bin/clang-format
    if [[ ! -f "$CLANGFORMAT_EXECUTABLE" ]]; then
        echo "clang-format not found: $CLANGFORMAT_EXECUTABLE"
        exit 1
    fi    
    cmake_args="$cmake_args -DENABLE_CLANGFORMAT=ON -DCLANGFORMAT_EXECUTABLE=$CLANGFORMAT_EXECUTABLE"
fi

or_die ./config-build.py -hc host-configs/docker/${HOST_CONFIG} -bp build-check-debug -ip install-check-debug $cmake_args
or_die cd build-check-debug

if [[ "$CHECK_TYPE" == "style" ]] ; then
    or_die make VERBOSE=1 clangformat_check
fi

exit 0
