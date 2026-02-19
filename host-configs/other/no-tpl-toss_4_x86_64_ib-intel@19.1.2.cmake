# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Minimal host-config (cmake cache file) to build Axom without any pre-built
# third party libraries.
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
#   cmake -C ../host-config/other/no-tpl-toss_4_x86_64_ib-intel@19.0.4.cmake \
#         -DCMAKE_BUILD_TYPE={Debug,Release}                                \
#         -DBUILD_SHARED_LIBS={ON,OFF(DEFAULT)}                             \
#         -DENABLE_EXAMPLES={ON,OFF}                                        \
#         -DENABLE_TESTS={ON,OFF}                                           \
#         -DCMAKE_INSTALL_PREFIX= /path/to/install/dir                      \
#         ../src
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/bin/cmake
#------------------------------------------------------------------------------

set(ENABLE_FORTRAN ON CACHE BOOL "")

set(COMPILER_HOME "/usr/tce/packages/intel-classic-tce/intel-classic-19.1.2" )
set(CMAKE_C_COMPILER "${COMPILER_HOME}/bin/icc" CACHE PATH "")
set(CMAKE_CXX_COMPILER "${COMPILER_HOME}/bin/icpc" CACHE PATH "")
set(CMAKE_Fortran_COMPILER "${COMPILER_HOME}/bin/ifort" CACHE PATH "")

# Sidre requires conduit and hdf5, so disable it in this host-config
# Inlet and Klee both depend on Sidre, so disable them too
set(AXOM_ENABLE_SIDRE OFF CACHE BOOL "")
set(AXOM_ENABLE_INLET OFF CACHE BOOL "")
set(AXOM_ENABLE_KLEE OFF CACHE BOOL "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------
set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_HOME             "/usr/tce/packages/mvapich2/mvapich2-2.3.7-intel-classic-19.1.2" CACHE PATH "")
set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpif90" CACHE PATH "")

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

#------------------------------------------------------------------------------
# !---------------------------------------------------------------------------
#------------------------------------------------------------------------------

