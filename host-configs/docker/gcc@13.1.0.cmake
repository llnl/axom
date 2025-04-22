#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/local/bin/cmake
#------------------------------------------------------------------------------

set(CMAKE_PREFIX_PATH "/home/axom/axom_tpls/gcc-13.1.0/blt-0.7.0-sdbpq4uu66agqgyreqkfr5hjebtlscaj;/home/axom/axom_tpls/gcc-13.1.0/caliper-2.12.1-2f52kccafe2q3z5xh72xdwqyoqfjmj5j;/home/axom/axom_tpls/gcc-13.1.0/conduit-git.2eec1ad9548cf6d629e2145181ff7b8cef0a38fa_0.9.3-o3m4pvvff4xazwycwxuqdwxwsevtyxns;/home/axom/axom_tpls/gcc-13.1.0/gmake-4.4.1-sh2qoa2zlsmuuji3nispoi3hx6sqnwkg;/home/axom/axom_tpls/gcc-13.1.0/mfem-4.7.0-tmh67qxbcma4ys6b76njipli6ue3uxex;/home/axom/axom_tpls/gcc-13.1.0/raja-2025.03.0-n4jjwyvdoikji4ywc5qj67lnmaeqtlg6;/home/axom/axom_tpls/gcc-13.1.0/umpire-2025.03.0-bsuzskxalfv75h2k5brt3glxtg5hxmay;/home/axom/axom_tpls/gcc-13.1.0/adiak-0.4.0-x5ti5yx4wcu5jjjny34ugfoo3il4xup6;/home/axom/axom_tpls/gcc-13.1.0/elfutils-0.192-msnl2dc4yclkbpbyfiqflurvhu3j476l;/home/axom/axom_tpls/gcc-13.1.0/libunwind-1.8.1-ohwx2aoeanit3nzd2o2lkg2tap4dp5it;/home/axom/axom_tpls/gcc-13.1.0/hdf5-1.8.23-qw6ua3mffwj2xapopw2sqsymhfdekxkn;/home/axom/axom_tpls/gcc-13.1.0/parmetis-4.0.3-wmb6op3qqxk5l5mf3z4pnzlgaavei6hy;/home/axom/axom_tpls/gcc-13.1.0/hypre-2.24.0-jz7h6o47ufl25jgfcvf453gutnmcfktd;/home/axom/axom_tpls/gcc-13.1.0/camp-2025.03.0-xp6sg6rdf7xp46dw57ovuxzwnewpfomu;/home/axom/axom_tpls/gcc-13.1.0/fmt-11.0.2-umbnbdyf25pufvp3dvsonj4d3pqunh2x;/home/axom/axom_tpls/gcc-13.1.0/libiconv-1.17-7unqwbz22fisucsxznfth72jym5i2uy2;/home/axom/axom_tpls/gcc-13.1.0/xz-5.4.6-hxeb4xtp65kazflulxyapiuuv7bjfakh;/home/axom/axom_tpls/gcc-13.1.0/zstd-1.5.6-6y4avtzsqdwhfxvto2en4t4dip7d2djc;/home/axom/axom_tpls/gcc-13.1.0/pkgconf-2.3.0-xyu57k7pftmjy36ggzxpk76m4ybmatyd;/home/axom/axom_tpls/gcc-13.1.0/zlib-ng-2.2.3-5tvk4dh6opay46nfkqkv457inspaabyj;/home/axom/axom_tpls/gcc-13.1.0/metis-5.1.0-3453ikebr6twodzzjguc2jpsyjsdnx2j;/home/axom/axom_tpls/gcc-13.1.0/gcc-runtime-13.1.0-2hzs4qwswjqszhprvd4mcyuw5lcwkph6" CACHE STRING "")

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH "ON" CACHE STRING "")

set(CMAKE_BUILD_RPATH "/home/axom/axom_tpls/gcc-13.1.0/axom-develop-45yhj5tmag3gqv6ijkdx4ii6fjei4hpv/lib;/home/axom/axom_tpls/gcc-13.1.0/axom-develop-45yhj5tmag3gqv6ijkdx4ii6fjei4hpv/lib64;;" CACHE STRING "")

set(CMAKE_INSTALL_RPATH "/home/axom/axom_tpls/gcc-13.1.0/axom-develop-45yhj5tmag3gqv6ijkdx4ii6fjei4hpv/lib;/home/axom/axom_tpls/gcc-13.1.0/axom-develop-45yhj5tmag3gqv6ijkdx4ii6fjei4hpv/lib64;;" CACHE STRING "")

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: gcc@=13.1.0
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/home/axom/axom_tpls/spack/lib/spack/env/gcc/gcc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/home/axom/axom_tpls/spack/lib/spack/env/gcc/g++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/home/axom/axom_tpls/spack/lib/spack/env/gcc/gfortran" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/bin/gcc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/bin/g++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/bin/gfortran-13" CACHE PATH "")

endif()

set(CMAKE_C_FLAGS "-pthread" CACHE STRING "")

set(CMAKE_CXX_FLAGS "-pthread" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/bin/mpif90" CACHE PATH "")

set(MPIEXEC_EXECUTABLE "/usr/bin/mpirun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING "")

set(ENABLE_MPI ON CACHE BOOL "")

#------------------------------------------------------------------------------
# Hardware
#------------------------------------------------------------------------------

#------------------------------------------------
# Hardware Specifics
#------------------------------------------------

set(ENABLE_OPENMP ON CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

set(TPL_ROOT "/home/axom/axom_tpls/gcc-13.1.0" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-git.2eec1ad9548cf6d629e2145181ff7b8cef0a38fa_0.9.3-o3m4pvvff4xazwycwxuqdwxwsevtyxns" CACHE PATH "")

# C2C not built

set(MFEM_DIR "${TPL_ROOT}/mfem-4.7.0-tmh67qxbcma4ys6b76njipli6ue3uxex" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.23-qw6ua3mffwj2xapopw2sqsymhfdekxkn" CACHE PATH "")

set(LUA_DIR "/usr" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2025.03.0-n4jjwyvdoikji4ywc5qj67lnmaeqtlg6" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2025.03.0-bsuzskxalfv75h2k5brt3glxtg5hxmay" CACHE PATH "")

# OPENCASCADE not built

set(ADIAK_DIR "${TPL_ROOT}/adiak-0.4.0-x5ti5yx4wcu5jjjny34ugfoo3il4xup6" CACHE PATH "")

set(CALIPER_DIR "${TPL_ROOT}/caliper-2.12.1-2f52kccafe2q3z5xh72xdwqyoqfjmj5j" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2025.03.0-xp6sg6rdf7xp46dw57ovuxzwnewpfomu" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

# ClangFormat disabled due to llvm@14 and devtools not in spec

set(ENABLE_CLANGFORMAT OFF CACHE BOOL "")

set(ENABLE_DOCS OFF CACHE BOOL "")



