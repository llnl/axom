#!/bin/sh
"exec" "python3" "-u" "-B" "$0" "$@"

# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

"""
 file: fix_permissions.py

 description: 
  Fixes the permissions for shared directories

"""

import sys
from llnl_lc_build_tools import *

def main():
    shared_dirs = [get_shared_mirror_dir(),
                   get_shared_devtool_dir(),
                   get_shared_libs_dir()]

    for directory in shared_dirs:
        set_group_and_perms(directory)

    return True


if __name__ == "__main__":
    if main():
        sys.exit(0)
    sys.exit(1)
