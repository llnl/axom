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

or_die cd axom
git submodule init 
git submodule update 

echo "~~~~ helpful info ~~~~"
echo "USER="`id -u -n`
echo "PWD="`pwd`
echo "HOST_CONFIG=$HOST_CONFIG"
echo "CMAKE_EXTRA_FLAGS=$CMAKE_EXTRA_FLAGS"
echo "~~~~~~~~~~~~~~~~~~~~~~"

echo "~~~~~~~~~ls -al /~~~~~~~~~~"
ls -al /
echo "~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~ls -al /home~~~~~~~~~~"
ls -al /home
echo "~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~ls -al /home/axom~~~~~~~~~~"
ls -al /home/axom
echo "~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~ls -al /home/axom/axom~~~~~~~~~~"
ls -al /home/axom/axom
echo "~~~~~~~~~~~~~~~~~~~~~~"


echo "~~~~~~ RUNNING CMAKE ~~~~~~~~"
or_die ./config-build.py -hc /home/axom/axom/host-configs/docker/${HOST_CONFIG}.cmake ${CMAKE_EXTRA_FLAGS}
or_die cd build-$HOST_CONFIG-debug
echo "~~~~~~ RUNNING make check ~~~~~~~~"
or_die make VERBOSE=1 check

exit 0
