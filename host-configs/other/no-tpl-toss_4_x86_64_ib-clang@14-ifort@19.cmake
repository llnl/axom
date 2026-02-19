# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Minimal host-config (cmake cache file) to build Axom without any pre-built
# third party libraries using the clang compiler for C++ 
# and the intel compiler  for Fortran
#
# Note: Sidre requires Conduit and is therefore disabled. Similarly, the Inlet
# and Klee components both depend on Sidre and must be disabled.
#
# To build the code with this host-config,
# run the following from your build directory:
#
#   cd <axom_root>
#   mkdir build
#   cd build
#   cmake -C ../host-config/other/no-tpl-toss_4_x86_64_ib-clang@14-ifort@19.cmake \
#         -DCMAKE_BUILD_TYPE={Debug,Release}                                      \
#         -DBUILD_SHARED_LIBS={ON,OFF(DEFAULT)}                                   \
#         -DENABLE_EXAMPLES={ON,OFF}                                              \
#         -DENABLE_TESTS={ON,OFF}                                                 \
#         -DCMAKE_INSTALL_PREFIX= /path/to/install/dir                            \
#         ../src
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/bin/cmake
#------------------------------------------------------------------------------

set(ENABLE_FORTRAN ON CACHE BOOL "")

set(_gnu_root   "/usr/tce/packages/gcc/gcc-8.3.1" CACHE PATH "")
set(_clang_root "/usr/tce/packages/clang/clang-14.0.6" CACHE PATH "")
set(_intel_root "/usr/tce/packages/intel-classic-tce/intel-classic-19.0.4" CACHE PATH "")

set(CMAKE_C_COMPILER   "${_clang_root}/bin/clang" CACHE PATH "")
set(CMAKE_CXX_COMPILER "${_clang_root}/bin/clang++" CACHE PATH "")
set(CMAKE_Fortran_COMPILER "${_intel_root}/bin/ifort" CACHE PATH "")

# Tell clang and the intel compiler to use a newer (non-default) version of gcc
# and remove implicit linker paths associated with the default gcc
set(CMAKE_C_FLAGS   "--gcc-toolchain=${_gnu_root} -gcc-name=${_gnu_root}/bin/gcc" CACHE STRING "")
set(CMAKE_CXX_FLAGS "--gcc-toolchain=${_gnu_root} -gcc-name=${_gnu_root}/bin/gcc" CACHE STRING "")
set(BLT_EXE_LINKER_FLAGS " -Wl,-rpath,${COMPILER_HOME}/lib" CACHE STRING "Adds a missing libstdc++ rpath")

set(BLT_CMAKE_IMPLICIT_LINK_DIRECTORIES_EXCLUDE
        ${_clang_root}/release/lib
        /usr/tce/packages/gcc/gcc-4.9.3/lib64/gcc/x86_64-unknown-linux-gnu/4.9.3
        /usr/tce/packages/gcc/gcc-4.9.3/lib64
    CACHE STRING "")

# Sidre requires conduit and hdf5, so disable it in this host-config
# Inlet and Klee both depend on Sidre, so disable them too
set(AXOM_ENABLE_SIDRE OFF CACHE BOOL "")
set(AXOM_ENABLE_INLET OFF CACHE BOOL "")
set(AXOM_ENABLE_KLEE  OFF CACHE BOOL "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------
set(ENABLE_MPI ON CACHE BOOL "")

set(_mpi_root            "/usr/tce/packages/mvapich2/mvapich2-2.3-clang-14.0.6" CACHE PATH "")
set(MPI_C_COMPILER       "${_mpi_root}/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER     "${_mpi_root}/bin/mpicxx" CACHE PATH "")
set(MPI_Fortran_COMPILER "${_mpi_root}/bin/mpif90" CACHE PATH "")

set(MPIEXEC_EXECUTABLE   "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

#------------------------------------------------------------------------------
# Other
#------------------------------------------------------------------------------
set(ENABLE_OPENMP OFF CACHE BOOL "")
set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")
set(ENABLE_DOCS OFF CACHE BOOL "")

set(CLANGFORMAT_EXECUTABLE "/usr/tce/packages/clang/clang-14.0.6/bin/clang-format" CACHE PATH "")
