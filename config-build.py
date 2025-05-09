#!/bin/sh

# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

"exec" "python3" "-u" "-B" "$0" "$@"

# Python wrapper script for generating the correct cmake line with the
# options specified by the user.
#
# Please keep parser option names as close to possible as the names of
# the cmake options they are wrapping.

from __future__ import print_function

import argparse
import distutils.spawn
import os
import platform
import shutil
import stat
import subprocess
import sys


def extract_cmake_location(file_path):
    if os.path.exists(file_path):
        cmake_line_prefix = "# cmake executable path: "
        file_handle = open(file_path, "r")
        content = file_handle.readlines()
        for line in content:
            if line.lower().startswith(cmake_line_prefix):
                return line.partition(":")[2].strip()
    print("Could not find a cmake entry in host config file.\n"
          "Attempting to find cmake on your path...")
    cmake_path = distutils.spawn.find_executable("cmake")
    print("Found: {0}".format(cmake_path))
    ret_cmake_path = "\"{0}\"".format(cmake_path)
    return ret_cmake_path


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Configure cmake build.",
        epilog="Note: Additional or unrecognized parameters will be passed"
        " directly to cmake."
        " For example, append '-DENABLE_OPENMP=ON' to enable OpenMP."
        " See Axom's QuickStart guide for a list of available CMake options.",
    )

    parser.add_argument(
        "-bp",
        "--buildpath",
        type=str,
        default="",
        help="Specify path for build directory."
        " If not specified, will create in current directory.",
    )

    parser.add_argument(
        "-ip",
        "--installpath",
        type=str,
        default="",
        help="Specify path for installation directory."
        " If not specified, will create in current directory.",
    )

    parser.add_argument(
        "-bt",
        "--buildtype",
        type=str,
        choices=["Release", "Debug", "RelWithDebInfo", "MinSizeRel"],
        default="Debug",
        help="Build type.",
    )

    parser.add_argument(
        "-e",
        "--eclipse",
        action="store_true",
        help="Create an eclipse project file.",
    )

    parser.add_argument(
        "-n", "--ninja", action="store_true", help="Create a Ninja project."
    )

    parser.add_argument(
        "-x", "--xcode", action="store_true", help="Create an xcode project."
    )

    # Add supported versions of MS Visual Studio as needed.
    msvcversions = {'2017': 'Visual Studio 15 2017',
                    '201764': 'Visual Studio 15 2017 Win64',
                    '2019': 'Visual Studio 16 2019',
                    '201964': 'Visual Studio 16 2019',
                    '2022': 'Visual Studio 17 2022',
                    '202264': 'Visual Studio 17 2022'}
    # Newer versions  of MSVC might supply an architecture flag
    generator_archs = {'2019': 'Win32',
                       '201964': 'x64',
                       '2022': 'Win32',
                       '202264': 'x64'}
    parser.add_argument(
        "--msvc",
        type=str,
        choices=msvcversions.keys(),
        help="Create a MS Visual Studio project."
    )

    parser.add_argument(
        "-ecc",
        "--exportcompilercommands",
        action="store_true",
        help="Generate a compilation database."
        " Can be used by the clang tools such as clang-modernize."
        " Will create a 'compile_commands.json' file in build directory.",
    )

    parser.add_argument(
        "-hc",
        "--hostconfig",
        required=True,
        type=str,
        help="Select a specific host-config file to initalize CMake's cache",
    )

    parser.add_argument(
        "--docs-only",
        action="store_true",
        help="Generate a configuration for working on documentation."
        " Disables features not required for building the docs.",
    )

    args, unknown_args = parser.parse_known_args()
    args.msvcversions = msvcversions
    args.generator_archs = generator_archs
    if unknown_args:
        print(
            "[config-build]: Passing the following arguments directly to cmake... %s"
            % unknown_args
        )

    return args, unknown_args


########################
# Find CMake Cache File
########################
def find_host_config(args):
    hostconfigpath = os.path.abspath(args.hostconfig)
    assert os.path.exists(hostconfigpath), (
        "Could not find CMake host config file '%s'." % hostconfigpath
    )
    print("Using host config file: '%s'." % hostconfigpath)
    return hostconfigpath


########################
# Get Platform information from host config name
########################
def get_platform_info(hostconfigpath):
    platform_info = ""
    platform_info = os.path.split(hostconfigpath)[1]
    if platform_info.endswith(".cmake"):
        platform_info = platform_info[:-6]
    return platform_info


#####################
# Setup Build Dir
#####################
def setup_build_dir(args, platform_info):
    if args.buildpath != "":
        # use explicit build path
        buildpath = args.buildpath
    else:
        # use platform info & build type
        pathList = ["build", platform_info]
        if not args.msvc:
            pathList.append(args.buildtype.lower())
        buildpath = "-".join(pathList)

    buildpath = os.path.abspath(buildpath)

    if os.path.exists(buildpath):
        print("Build directory '%s' already exists.  Deleting..." % buildpath)
        shutil.rmtree(buildpath)

    print("Creating build directory '%s'..." % buildpath)
    os.makedirs(buildpath)
    return buildpath


#####################
# Setup Install Dir
#####################
def setup_install_dir(args, platform_info):
    # For install directory, we will clean up old ones, but we don't need to create it, cmake will do that.
    if args.installpath != "":
        installpath = os.path.abspath(args.installpath)
    else:
        # use platform info & build type
        pathList = ["install", platform_info]
        if not args.msvc:
            pathList.append(args.buildtype.lower())
        installpath = "-".join(pathList)

    installpath = os.path.abspath(installpath)

    if os.path.exists(installpath):
        print(
            "Install directory '%s' already exists, deleting..." % installpath
        )
        shutil.rmtree(installpath)

    print("Creating install path '%s'..." % installpath)
    os.makedirs(installpath)
    return installpath


############################
# Check if executable exists
############################
def executable_exists(path):
    if path == "cmake":
        return True
    return os.path.isfile(path) and os.access(path, os.X_OK)


############################
# Build CMake command line
############################
def create_cmake_command_line(
    args, unknown_args, buildpath, hostconfigpath, installpath
):
    cmakeline = extract_cmake_location(hostconfigpath)
    cmakeline = os.path.normpath(cmakeline) # Fixes path for Windows
    assert cmakeline != None, ("No cmake executable found on path")
    assert executable_exists(cmakeline), (
        "['%s'] invalid path to cmake executable or file does not have execute permissions"
        % cmakeline
    )

    # create the ccmake command for convenience
    cmakedir = os.path.dirname(cmakeline)
    ccmake_cmd = os.path.join(cmakedir, "ccmake")
    if executable_exists(ccmake_cmd):
        # write the ccmake command to a file to use for convenience
        ccmake_file = os.path.join(buildpath, "ccmake_cmd")
        with open(ccmake_file, "w") as ccmakefile:
            ccmakefile.write("#!/usr/bin/env bash\n")
            ccmakefile.write(ccmake_cmd)
            ccmakefile.write(" $@")
            ccmakefile.write("\n")

        st = os.stat(ccmake_file)
        os.chmod(ccmake_file, st.st_mode | stat.S_IEXEC)

    # Add cache file option
    cmakeline = '"{0}" -C {1}'.format(cmakeline,hostconfigpath)
    # Add build type (opt or debug); don't add for msvc generator
    if not args.msvc:
        cmakeline += " -DCMAKE_BUILD_TYPE=" + args.buildtype
    # Set install dir
    cmakeline += " -DCMAKE_INSTALL_PREFIX=%s" % installpath

    if args.exportcompilercommands:
        cmakeline += " -DCMAKE_EXPORT_COMPILE_COMMANDS=on"

    if args.eclipse:
        cmakeline += ' -G "Eclipse CDT4 - Unix Makefiles"'

    if args.ninja:
        cmakeline += ' -G "Ninja"'

    if args.xcode:
        cmakeline += ' -G "Xcode"'

    if args.msvc:
        cmakeline += ' -G "%s"' % args.msvcversions[args.msvc]
        if args.generator_archs.get(args.msvc):
            cmakeline += ' -A %s' % args.generator_archs[args.msvc]

    if args.docs_only:
        cmakeline += " -DENABLE_ALL_COMPONENTS=OFF"
        cmakeline += " -DAXOM_ENABLE_TESTS=OFF"
        cmakeline += " -DAXOM_ENABLE_EXAMPLES=OFF"
        cmakeline += " -DAXOM_ENABLE_DOCS=ON"

    if unknown_args:
        cmakeline += " " + " ".join(unknown_args)

    rootdir = os.path.dirname(os.path.abspath(sys.argv[0]))
    cmakeline += " %s " % os.path.join(rootdir, "src")

    # Dump the cmake command to file for convenience
    cmake_file = os.path.join(buildpath, "cmake_cmd")
    with open(cmake_file, "w") as cmdfile:
        cmdfile.write(cmakeline)
        cmdfile.write("\n")

    st = os.stat(cmake_file)
    os.chmod(cmake_file, st.st_mode | stat.S_IEXEC)
    return cmakeline


############################
# Run CMake
############################
def run_cmake(buildpath, cmakeline):
    print("Changing to build directory...")
    os.chdir(buildpath)
    print("Executing CMake line: '%s'" % cmakeline)
    print("")

    returncode = subprocess.call(cmakeline, shell=True)
    if not returncode == 0:
        print(
            "Error: CMake command failed with return code: {0}".format(
                returncode
            )
        )
        return False
    return True


############################
# Main
############################
def main():
    args, unknown_args = parse_arguments()

    hostconfigpath = find_host_config(args)
    platform_info = get_platform_info(hostconfigpath)
    buildpath = setup_build_dir(args, platform_info)
    installpath = setup_install_dir(args, platform_info)

    cmakeline = create_cmake_command_line(
        args, unknown_args, buildpath, hostconfigpath, installpath
    )
    return run_cmake(buildpath, cmakeline)


if __name__ == "__main__":
    exit(0 if main() else 1)
