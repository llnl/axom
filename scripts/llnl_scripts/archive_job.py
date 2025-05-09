#!/bin/sh
"exec" "python3" "-u" "-B" "$0" "$@"

# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

"""
 file: archive_job.py

 description: 
  Archives results from a non-source or TPL job.

  If an exit code is passed to this script, it will be used as the exit code
  for this script.  This is so that archiving the results of an automation
  job does not alter the overall exit status of the job.

  If no "failed.json" or "success.json" does not already exist in the current
  directory then the given exit code of the script is used. Otherwise it defaults
  to failure.

"""

import os
import sys

from llnl_lc_build_tools import *
from argparse import ArgumentParser


def parse_args():
    "Parses args from command line"
    parser = ArgumentParser()
    
    parser.add_argument("-e", "--exitcode",
                        type=int,
                        dest="exitcode",
                        default="0",
                        help="Exit code to be returned from this script (Defaults to 0)")

    parser.add_argument("-n", "--name",
                        dest="name",
                        default="",
                        help="Archive build results under given name")

    parser.add_argument("-f", "--files",
                        dest="files",
                        default="",
                        help="Files to be archived (comma delimited)")

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

    timestamp = get_timestamp()
    archive_dir = pjoin(get_archive_base_dir(), get_system_type())
    archive_dir = pjoin(archive_dir, normalize_job_name(opts["name"]), timestamp)

    cwd = os.getcwd()

    print("[Starting Archiving]")
    print(f"[  Archive Dir: {archive_dir}]")
    print(f"[  Job Dir: {cwd}]")

    # Create a success/failure file if not already created by the job
    if (opts["exitcode"] != 0) and not os.path.exists("failed.json"):
        log_failure(cwd, opts["name"], timestamp)
    elif (opts["exitcode"] == 0) and not os.path.exists("success.json"):
        log_success(cwd, opts["name"], timestamp)
    else:
        log_failure(cwd, opts["name"] + " - No job status given", timestamp)

    # Copy any existing info files
    if not os.path.exists(archive_dir):
        os.makedirs(archive_dir)

    copy_if_exists(pjoin(cwd, "info.json"), archive_dir)
    copy_if_exists(pjoin(cwd, "failed.json"), archive_dir)
    copy_if_exists(pjoin(cwd, "success.json"), archive_dir)

    # Copy explicitly given files
    for file in opts["files"].split(','):
        filepath = file.strip()
        if filepath != "":
            copy_if_exists(file, archive_dir)

    return opts["exitcode"]


if __name__ == "__main__":
    # Always exit with the error code provided by user to
    # preserve for reporting the automated job status
    exitcode = main()
    sys.exit(exitcode)
