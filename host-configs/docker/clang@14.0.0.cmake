#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/local/bin/cmake
#------------------------------------------------------------------------------

set(CMAKE_PREFIX_PATH "/home/axom/axom_tpls/clang-14.0.0/blt-0.7.0-sr43mf6inrcrwhao4ckioacpgx5ohupn;/home/axom/axom_tpls/clang-14.0.0/caliper-2.12.1-pxlgplvlumvexjwnuxdlomr5yaciwf7r;/home/axom/axom_tpls/clang-14.0.0/conduit-git.2eec1ad9548cf6d629e2145181ff7b8cef0a38fa_0.9.3-glcvpl5j5fialqd5sdd75zwsezw7wx64;/home/axom/axom_tpls/clang-14.0.0/gmake-4.4.1-qakmu75fp2cyfl36ks3il6e4iej4y7ro;/home/axom/axom_tpls/clang-14.0.0/mfem-4.7.0-hnvseg4gdm6liq7cr5pyyowemgmnjoz4;/home/axom/axom_tpls/clang-14.0.0/raja-2025.03.0-h4mada7t7ncjm3imnmmksncix4kqc7ll;/home/axom/axom_tpls/clang-14.0.0/umpire-2025.03.0-akib7tba4p5zb4id4upcqdiusnmpl6k4;/home/axom/axom_tpls/clang-14.0.0/adiak-0.4.0-lesftqm5cxjgipnpyyyzwhcam4tyjkrz;/home/axom/axom_tpls/clang-14.0.0/elfutils-0.192-5rbp4nbb62gaoswwa3o4sxhsmtpucihz;/home/axom/axom_tpls/clang-14.0.0/libunwind-1.8.1-xbkl6hi2ej7k2iyt2vtl55kkhb7wgz4w;/home/axom/axom_tpls/clang-14.0.0/hdf5-1.8.23-7usjsaiyavaxncd2v7c2bdr5fruy2ept;/home/axom/axom_tpls/clang-14.0.0/parmetis-4.0.3-urgxlhmqfbzdbdfcto6h6ikfdud3h4vj;/home/axom/axom_tpls/clang-14.0.0/hypre-2.24.0-lukvwkcp3qdxs3iathjztwh5gnd6styd;/home/axom/axom_tpls/clang-14.0.0/camp-2025.03.0-ut44gujidwayhspxgne5ni7flacvtr5s;/home/axom/axom_tpls/clang-14.0.0/fmt-11.0.2-yc4tb6krft2dzioamgh7eoxvpkpggjud;/home/axom/axom_tpls/clang-14.0.0/libiconv-1.17-7ejia6pnbidhxea2zksbdyhaojbocoyd;/home/axom/axom_tpls/clang-14.0.0/xz-5.4.6-upnmx7x5jyrstnqy4sscoqhrm5vhdtqv;/home/axom/axom_tpls/clang-14.0.0/zstd-1.5.6-z53wj4a4xffogqutqojer5utkfmetuce;/home/axom/axom_tpls/clang-14.0.0/pkgconf-2.3.0-2op2c4aekepjrqoemhrjbx5t6muhbabo;/home/axom/axom_tpls/clang-14.0.0/zlib-ng-2.2.3-hspizze2o4s3hzjvzllxdkvsqrqitctb;/home/axom/axom_tpls/clang-14.0.0/metis-5.1.0-r6rm2szcrtuf6k3emyax67xewz4v5tec" CACHE STRING "")

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH "ON" CACHE STRING "")

set(CMAKE_BUILD_RPATH "/home/axom/axom_tpls/clang-14.0.0/axom-develop-ewggen3qgytv5h3rf44fl4iahlglc4nv/lib;/home/axom/axom_tpls/clang-14.0.0/axom-develop-ewggen3qgytv5h3rf44fl4iahlglc4nv/lib64;;" CACHE STRING "")

set(CMAKE_INSTALL_RPATH "/home/axom/axom_tpls/clang-14.0.0/axom-develop-ewggen3qgytv5h3rf44fl4iahlglc4nv/lib;/home/axom/axom_tpls/clang-14.0.0/axom-develop-ewggen3qgytv5h3rf44fl4iahlglc4nv/lib64;;" CACHE STRING "")

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: clang@=14.0.0
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/home/axom/axom_tpls/spack/lib/spack/env/clang/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/home/axom/axom_tpls/spack/lib/spack/env/clang/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/home/axom/axom_tpls/spack/lib/spack/env/clang/gfortran" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/bin/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/bin/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/bin/gfortran-11" CACHE PATH "")

endif()

set(CMAKE_C_FLAGS "-fPIC -pthread" CACHE STRING "")

set(CMAKE_CXX_FLAGS "-fPIC -pthread" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

set(BLT_EXE_LINKER_FLAGS " -Wl,-rpath,/usr/lib -Wl,-rpath,/usr/lib64" CACHE STRING "Adds a missing libstdc++ rpath")

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

set(TPL_ROOT "/home/axom/axom_tpls/clang-14.0.0" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-git.2eec1ad9548cf6d629e2145181ff7b8cef0a38fa_0.9.3-glcvpl5j5fialqd5sdd75zwsezw7wx64" CACHE PATH "")

# C2C not built

set(MFEM_DIR "${TPL_ROOT}/mfem-4.7.0-hnvseg4gdm6liq7cr5pyyowemgmnjoz4" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.23-7usjsaiyavaxncd2v7c2bdr5fruy2ept" CACHE PATH "")

set(LUA_DIR "/usr" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2025.03.0-h4mada7t7ncjm3imnmmksncix4kqc7ll" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2025.03.0-akib7tba4p5zb4id4upcqdiusnmpl6k4" CACHE PATH "")

# OPENCASCADE not built

set(ADIAK_DIR "${TPL_ROOT}/adiak-0.4.0-lesftqm5cxjgipnpyyyzwhcam4tyjkrz" CACHE PATH "")

set(CALIPER_DIR "${TPL_ROOT}/caliper-2.12.1-pxlgplvlumvexjwnuxdlomr5yaciwf7r" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2025.03.0-ut44gujidwayhspxgne5ni7flacvtr5s" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

# ClangFormat disabled due to llvm@14 and devtools not in spec

set(ENABLE_CLANGFORMAT OFF CACHE BOOL "")

set(ENABLE_DOCS OFF CACHE BOOL "")



