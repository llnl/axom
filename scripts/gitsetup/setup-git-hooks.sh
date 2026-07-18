#!/usr/bin/env bash

# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#------------------------------------------------------------------------------
# This script installs client-side hooks
#------------------------------------------------------------------------------
basedir=`git rev-parse --show-toplevel`
hooksdir="$basedir/.git/hooks/"
cp -v $basedir/scripts/gitsetup/hooks/commit-msg $hooksdir
cp -v $basedir/scripts/gitsetup/hooks/post-checkout $hooksdir
