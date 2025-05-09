#!/usr/bin/env bash

# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

##======================================================
#
# This script is used to set colors in git
#
##======================================================

# Project configuration instructions: NONE

read -ep 'Enable syntax highlighting of output from git commands? [Y/n] ' color &&
if test "$color" = "y" -o "$color" = "Y"; then
    git config --global color.ui true
else
    git config --global color.ui false
fi

