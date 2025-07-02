#!/bin/bash
# NOTE: This script attempts to automate updating the fmt library in Axom.

function clone_fmt
{
  # Clone repo
  echo "Cloning fmt repo..."
  git clone https://github.com/fmtlib/fmt.git
}

function copy_headers
{
  # Copy headers into Axom.
  echo "Copying fmt headers..."
  cp fmt/include/fmt/*.h .
}

function patch_file
{
  # Patch copied file
  if patch -p1 $1 < $2 ; then
    echo "Applied patch $2 to $1. Updating patch file $2."
    # Generate diff to update patch.
    AXOM_FMT=$(pwd)
    cd fmt/include/fmt
    git diff format-inl.h $AXOM_FMT/$1 > $AXOM_FMT/$2
    cd $AXOM_FMT
  else
    echo "Patch $2 failed for $1. Not generating diff."
  fi
}

function apply_patches
{
  patch_file format-inl.h runtime_error.patch
  patch_file format.h     hipcc_long_double.patch
}

function modify_headers
{
  # Make some AXOM name replacements
  echo "Renaming FMT to Axom in files..."
  for f in $(ls *.h) ; do
    echo $f
    sed -i "s/FMT_/AXOM_FMT_/g" $f
    sed -i "s/fmt::/axom::fmt::/g" $f
  done
}

function revert
{
  # Revert any changes
  git checkout -- args.h
  git checkout -- base.h
  git checkout -- chrono.h
  git checkout -- color.h
  git checkout -- compile.h
  git checkout -- core.h
  git checkout -- format.h
  git checkout -- format-inl.h
  git checkout -- os.h
  git checkout -- ostream.h
  git checkout -- printf.h
  git checkout -- ranges.h
  git checkout -- std.h
  git checkout -- xchar.h
#  git checkout -- runtime_error.patch
#  git checkout -- hipcc_long_double.patch
}

function cleanup
{
  rm -rf fmt
}

clone_fmt
revert
copy_headers
apply_patches
#modify_headers
#cleanup
