# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#------------------------------------------------------------------------------
# host-config for 'sqa-uno' Windows machine-- using Visual Studio 15
#
# Run the following from a build dir for a 32-bit configuration
#   cmake -G "Visual Studio 15 2017"                       \
#         -C ..\host-configs\other\sqa-uno-MSVC.cmake      \
#         <path-to-axom>
#
# Run the following from a build dir for a 64-bit configuration
#   cmake -G "Visual Studio 15 2017 Win64"               \
#         -C ..\host-configs\other\sqa-uno-MSVC.cmake    \
#         <path-to-axom> 
#
# Build the code from the command line as follows (/m for parallel build in msbuild):
#   cmake --build . --config {Release, Debug, RelWithDebInfo} [-- /m:8]
#
# Test the code as follows (j for parallel testing):
#   ctest -j8 -C {Release,Debug,RelWithDebInfo}
# 
# Install the code from the command line as follows:
#   cmake --build . --config Release --target install
#
#------------------------------------------------------------------------------

# Set the HOME variable (%USERPROFILE% in Windows)
string(REPLACE "\\" "/" HOME "$ENV{USERPROFILE}")

### Setup some devtool/TPL paths

# Sidre will use conduit but not hdf5 until we resolve configuration problems with hdf5
#set(HDF5_DIR    ... CACHE PATH "")
set(CONDUIT_DIR "${HOME}/Projects/install/conduit-no-hdf5" CACHE PATH "")

# Note: Doxygen assumes graphviz 'dot' is in PATH
set(DOXYGEN_EXECUTABLE "${HOME}/Chocolatey/bin/doxygen.exe" CACHE PATH "") 

# Note: uncrustify disabled since this version is different than uberenv version and changes whitespace
# set(UNCRUSTIFY_EXECUTABLE "${HOME}/Code/UniversalIndentGUI/indenters/uncrustify.exe" CACHE PATH "")

# Note: I had to modify the mfem's runtime library flags to link to mfem
#       due to BLT's conversion from \MD to \MT
set(MFEM_DIR "${HOME}/Projects/mfem/install-win" CACHE PATH "")

# Setup MPI
set(ENABLE_MPI ON CACHE BOOL "")
set(MPI_HOME  "C:/Program Files/Microsoft HPC Pack 2008 R2" CACHE PATH "")
set(MPI_GUESS_LIBRARY_NAME "MSMPI" CACHE STRING "")


### Set some additional options
set(ENABLE_FOLDERS ON CACHE BOOL "")
set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

### Disable warning about deprecated TR1 namespace in gtest
set(BLT_DEFINES "/D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING" CACHE STRING "")
