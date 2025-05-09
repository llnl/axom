#!/usr/bin/env bash

# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

##======================================================
#
# This script is used to setup git user-settings
#
##======================================================

# Project configuration instructions: NONE

for (( ; ; )); do
    user_name=$(git config user.name || echo '') &&
    user_email=$(git config user.email || echo '') &&
    if test -n "$user_name" -a -n "$user_email"; then
        echo 'Your commits will be recorded as: '"$user_name <$user_email>"'' &&
        read -ep 'Is the author name and e-mail address above correct? [Y/n] ' correct &&
        if test "$correct" != "n" -a "$correct" != "N"; then
            break;
        fi
    fi &&
    read -ep 'Enter full name to associate with your commits, e.g., "John Smith": ' name &&
    read -ep 'Enter corresponding email address, e.g., "smith1@llnl.gov" ' email &&
    git config --global user.name "$name" &&
    git config --global user.email "$email"
done

git config --global push.default current
