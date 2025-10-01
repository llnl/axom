#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/local/bin/cmake
#------------------------------------------------------------------------------

set(CMAKE_PREFIX_PATH "/home/axom/axom_tpls/llvm-14.0.0/blt-develop-p6bopj3cxqjef23aka7us72y7x3spzfh;/home/axom/axom_tpls/llvm-14.0.0/caliper-2.12.1-zqwmnodbtvw2wkt22pkufjssuej6biwa;/home/axom/axom_tpls/llvm-14.0.0/conduit-0.9.5-tpyznnm4j23kzrpljqdipsswis4mowob;/home/axom/axom_tpls/llvm-14.0.0/gmake-4.4.1-h2a5czi4uuijaus643wgvwkcwbks3sqj;/home/axom/axom_tpls/llvm-14.0.0/mfem-4.8.0-mlk223ipjs7iijeovplti4lwweyetbhj;/home/axom/axom_tpls/llvm-14.0.0/raja-2025.09.0-6bipn52puk3bpnt2llpthqxtfovgofua;/home/axom/axom_tpls/llvm-14.0.0/umpire-2025.09.0-nrszymzzuznjqvsv32j7mqwifzpcjxsn;/home/axom/axom_tpls/llvm-14.0.0/adiak-0.4.0-jncfkug5buwytnbmhdlg3j2fs3ie6zhf;/home/axom/axom_tpls/llvm-14.0.0/elfutils-0.192-o7hxw5srwdq7prlerzxoyyezideymdb7;/home/axom/axom_tpls/llvm-14.0.0/libunwind-1.8.1-6vft6cdtkxjxhzshq2g4yeuzlq26opeh;/home/axom/axom_tpls/llvm-14.0.0/hdf5-1.8.23-xc7wug34cpu6a7yw6sxon5vm5lkghosx;/home/axom/axom_tpls/llvm-14.0.0/parmetis-4.0.3-rwqzxe2sg3mfnt7mmsxx7263cn5tbpdj;/home/axom/axom_tpls/llvm-14.0.0/hypre-2.24.0-ymqswk5kn7uzbhswzq6oqnoen7dj2hft;/home/axom/axom_tpls/llvm-14.0.0/camp-2025.09.2-xtvodoeqn5pwclptputtnnukoxfpgqaq;/home/axom/axom_tpls/llvm-14.0.0/fmt-11.0.2-hk2r7jbfde722bfehxbaoux3ovcycxs4;/home/axom/axom_tpls/llvm-14.0.0/libiconv-1.18-pth2h6ciltbgluyhrvncolkwxgxidwoa;/home/axom/axom_tpls/llvm-14.0.0/xz-5.6.3-mix5agvh2hbotpbd7qglduhlpbbyzds7;/home/axom/axom_tpls/llvm-14.0.0/zstd-1.5.7-uslpcor5hnodwh6gohoswyfuoidcwirz;/home/axom/axom_tpls/llvm-14.0.0/zlib-ng-2.2.4-av76bkzrf6ef4kx7yiswmg3rfvaygrws;/home/axom/axom_tpls/llvm-14.0.0/metis-5.1.0-6yo5m45psv5bak2rdpaiwshpd2vrwanr;/home/axom/axom_tpls/none-none/gcc-runtime-11.4.0-f6fb6rs3jshnld7slkxw345mcqtghapc;/home/axom/axom_tpls/none-none/compiler-wrapper-1.0-myvi4leujrerocfg4bhvgrlr5ralanoj;/usr/lib/llvm-14" CACHE STRING "")

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH "ON" CACHE STRING "")

set(CMAKE_BUILD_RPATH "/home/axom/axom_tpls/llvm-14.0.0/axom-develop-5vylsrbcbxqsyejrgrh7egcbf77muhxl/lib;/home/axom/axom_tpls/llvm-14.0.0/axom-develop-5vylsrbcbxqsyejrgrh7egcbf77muhxl/lib64;;" CACHE STRING "")

set(CMAKE_INSTALL_RPATH "/home/axom/axom_tpls/llvm-14.0.0/axom-develop-5vylsrbcbxqsyejrgrh7egcbf77muhxl/lib;/home/axom/axom_tpls/llvm-14.0.0/axom-develop-5vylsrbcbxqsyejrgrh7egcbf77muhxl/lib64;;" CACHE STRING "")

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: llvm@14.0.0/ze3omayhuhvvpords3igc43pabd5f3jz
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/home/axom/axom_tpls/none-none/compiler-wrapper-1.0-myvi4leujrerocfg4bhvgrlr5ralanoj/libexec/spack/clang/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/home/axom/axom_tpls/none-none/compiler-wrapper-1.0-myvi4leujrerocfg4bhvgrlr5ralanoj/libexec/spack/clang/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/home/axom/axom_tpls/none-none/compiler-wrapper-1.0-myvi4leujrerocfg4bhvgrlr5ralanoj/libexec/spack/gcc/gfortran" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/lib/llvm-14/bin/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/lib/llvm-14/bin/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/bin/gfortran-11" CACHE PATH "")

endif()

set(CMAKE_C_FLAGS "-fPIC -pthread" CACHE STRING "")

set(CMAKE_CXX_FLAGS "-fPIC -pthread" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

set(BLT_EXE_LINKER_FLAGS " -Wl,-rpath,/usr/lib/llvm-14/lib" CACHE STRING "Adds a missing libstdc++ rpath")

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

set(TPL_ROOT "/home/axom/axom_tpls/llvm-14.0.0" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.9.5-tpyznnm4j23kzrpljqdipsswis4mowob" CACHE PATH "")

# C2C not built

set(MFEM_DIR "${TPL_ROOT}/mfem-4.8.0-mlk223ipjs7iijeovplti4lwweyetbhj" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.23-xc7wug34cpu6a7yw6sxon5vm5lkghosx" CACHE PATH "")

set(LUA_DIR "/usr" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2025.09.0-6bipn52puk3bpnt2llpthqxtfovgofua" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2025.09.0-nrszymzzuznjqvsv32j7mqwifzpcjxsn" CACHE PATH "")

# OPENCASCADE not built

set(ADIAK_DIR "${TPL_ROOT}/adiak-0.4.0-jncfkug5buwytnbmhdlg3j2fs3ie6zhf" CACHE PATH "")

set(CALIPER_DIR "${TPL_ROOT}/caliper-2.12.1-zqwmnodbtvw2wkt22pkufjssuej6biwa" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2025.09.2-xtvodoeqn5pwclptputtnnukoxfpgqaq" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

# ClangFormat disabled since llvm@19 and devtools not in spec

set(ENABLE_CLANGFORMAT OFF CACHE BOOL "")

set(ENABLE_DOCS OFF CACHE BOOL "")


