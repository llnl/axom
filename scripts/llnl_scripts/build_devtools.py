#!/bin/sh
"exec" "python3" "-u" "-B" "$0" "$@"

# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

"""
 file: build_devtools.py

 description: 
  Builds all Axom Devtools

"""

from llnl_lc_build_tools import *

from argparse import ArgumentParser

import getpass
import os


def parse_args():
    "Parses args from command line"
    parser = ArgumentParser()
    # Location of source directory to build
    parser.add_argument("-d", "--directory",
                        dest="directory",
                        default="",
                        help="Location to build all TPL's, timestamp directory will be created (Defaults to shared location)")

    ###############
    # parse args
    ###############
    opts = parser.parse_args()
    # we want a dict b/c the values could 
    # be passed without using argparse
    opts = vars(opts)
    return opts


def main():
    opts = parse_args()

    # Determine location to do all the building
    if opts["directory"] != "":
        build_dir = opts["directory"]
        if not os.path.exists(build_dir):
            os.makedirs(build_dir)
    else:
        if getpass.getuser() != "atk":
            print("ERROR: Only shared user 'atk' can install into shared directory. Use -d option.")
            return 1
        build_dir = get_shared_devtool_dir()
    build_dir = os.path.abspath(build_dir)

    repo_dir = get_repo_dir()

    try:
        original_wd = os.getcwd()
        os.chdir(repo_dir)

        res = build_devtools(build_dir, get_timestamp())
    finally:
        os.chdir(original_wd)

    return res


if __name__ == "__main__":
    sys.exit(main())
