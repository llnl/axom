#!/bin/bash
# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

# NOTE: This script attempts to automate updating the fmt library in Axom.

function clone_fmt
{
  local ver="$1"

  # Clone repo
  echo "Cloning fmt repo..."
  git clone https://github.com/fmtlib/fmt.git fmt_src
  echo "Checking out version $ver"
  cd fmt_src
  git checkout $ver > /dev/null
  cd ..
}

function copy_headers
{
  # Copy headers into Axom.
  echo "Copying fmt headers..."
  cp fmt_src/include/fmt/*.h fmt
}

function patch_file
{
  # Patch copied file
  cp fmt_src/include/fmt/$1 $1
  if patch -p1 $1 < $2 ; then
    echo "Applied patch $2 to $1."
    echo "Updating patch file $2."
    # Generate diff to update patch.
    diff -u fmt_src/include/fmt/$1 $1 > $2

    mv $1 fmt/$1
  else
    echo "Patch $2 failed for $1. Not generating diff."
    echo "YOU WILL NEED TO PORT THIS FILE MANUALLY AND MAKE A DIFF!"
    rm -f $1
  fi
}

function apply_patches
{
  patch_file base.h       fmt.base.h.patch
  patch_file format-inl.h fmt.runtime_error.patch
  patch_file format.h     fmt.format.h.patch
}

function modify_headers
{
  # Make some AXOM name replacements
  echo "Renaming FMT to Axom in files..."
  cd fmt
  for f in $(ls *.h) ; do
    echo $f
    sed -i "s/FMT_/AXOM_FMT_/g" $f
    sed -i "s/fmt::/axom::fmt::/g" $f
    sed -i "s/namespace fmt/namespace axom::fmt/g" $f
  done
  cd ..
}

function revert
{
  # Revert any changes
  git checkout -- fmt/args.h
  git checkout -- fmt/base.h
  git checkout -- fmt/chrono.h
  git checkout -- fmt/color.h
  git checkout -- fmt/compile.h
  git checkout -- fmt/core.h
  git checkout -- fmt/format.h
  git checkout -- fmt/format-inl.h
  git checkout -- fmt/os.h
  git checkout -- fmt/ostream.h
  git checkout -- fmt/printf.h
  git checkout -- fmt/ranges.h
  git checkout -- fmt/std.h
  git checkout -- fmt/xchar.h

  git checkout -- runtime_error.patch
  git checkout -- format.h.patch
  git checkout -- base.h.patch
}

function cleanup
{
  rm -rf fmt_src
}

function print_help()
{
    cat <<EOF
Usage: $0 [OPTIONS]

Options:
  --version VER   Set the version string (default: ${VERSION})
  --help          Show this help message and exit
EOF
}

function main
{
  # Default version
  VERSION="12.1.0"

  # Parse arguments
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --version)
        if [[ $# -lt 2 ]]; then
          echo "Error: --version requires an argument" >&2
          exit 1
        fi
        VERSION="$2"
        shift 2
        ;;
      --help)
        print_help
        exit 0
        ;;
      *)
        echo "Unknown option: $1" >&2
        echo "Use --help for usage." >&2
        exit 1
        ;;
    esac
  done

  cleanup
  clone_fmt $VERSION
  revert
  copy_headers
  apply_patches
  modify_headers
  cleanup
}

main $@
