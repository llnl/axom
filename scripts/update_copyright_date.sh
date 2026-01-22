#!/usr/bin/env bash

# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#=============================================================================
# This script is used to update the copyright year in RAJA files that contain
# specific year info.
#
# To use, change the 'sed' commands below as needed to modify the content that
# is changed in each file and what it is changed to.
#
#=============================================================================

echo LICENSE
cp LICENSE LICENSE.sed.bak
sed "s/Copyright (c) 2017-2025/Copyright (c) 2017-2026/" LICENSE.sed.bak > LICENSE
rm LICENSE.sed.bak

echo src/conf.py
cp src/conf.py src/conf.py.sed.bak
sed "s/2017-2025/2017-2026/" src/conf.py.sed.bak > src/conf.py
rm src/conf.py.sed.bak

