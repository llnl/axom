# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Defines CMake options for Axom's build system
#------------------------------------------------------------------------------

option(AXOM_ENABLE_ASAN "Enable AddressSanitizer for memory checking (Clang or GCC only)" OFF)
if(AXOM_ENABLE_ASAN)
    if(NOT (C_COMPILER_FAMILY_IS_CLANG OR C_COMPILER_FAMILY_IS_GNU))
        message(FATAL_ERROR "AXOM_ENABLE_ASAN only supports Clang and GCC")
    endif()
endif()

option(AXOM_ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer for undefined behavior detection (Clang or GCC only)" OFF)
if(AXOM_ENABLE_UBSAN)
    if(NOT (C_COMPILER_FAMILY_IS_CLANG OR C_COMPILER_FAMILY_IS_GNU))
        message(FATAL_ERROR "AXOM_ENABLE_UBSAN only supports Clang and GCC")
    endif()
endif()

option(AXOM_ENABLE_SPARSEHASH "Enables Sparsehash." ON)
option(AXOM_ENABLE_ALL_COMPONENTS "Enables all components by default" ON)
option(AXOM_USE_64BIT_INDEXTYPE "Use 64-bit integers for axom::IndexType" ON)
option(AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION "Enable Axom's version of the MFEM SidreDataCollection" ON)

if(NOT CMAKE_CONFIGURATION_TYPES)
    if(CMAKE_BUILD_TYPE MATCHES "(Debug|RelWithDebInfo)")
        option(AXOM_ENABLE_EXPORTS "Add in symbols to demangle axom function names in stacktraces" ON)
    else()
    	option(AXOM_ENABLE_EXPORTS "Add in symbols to demangle axom function names in stacktraces" OFF)
    endif()
endif()

cmake_dependent_option(AXOM_ENABLE_CUDA "Enables Axom with CUDA support" ON "ENABLE_CUDA" OFF)
cmake_dependent_option(AXOM_ENABLE_HIP "Enables Axom with HIP support" ON "ENABLE_HIP" OFF)
cmake_dependent_option(AXOM_ENABLE_MPI "Enables Axom with MPI support" ON "ENABLE_MPI" OFF)
cmake_dependent_option(AXOM_ENABLE_OPENMP "Enables Axom with OPENMP support" ON "ENABLE_OPENMP" OFF)

cmake_dependent_option(AXOM_ENABLE_TESTS "Enables Axom Tests" ON "ENABLE_TESTS" OFF)
cmake_dependent_option(AXOM_ENABLE_PYTHON_TESTS "Enables Axom Python Tests" ON "ENABLE_TESTS" OFF)
cmake_dependent_option(AXOM_ENABLE_DOCS "Enables Axom Docs" ON "ENABLE_DOCS" OFF)
cmake_dependent_option(AXOM_ENABLE_EXAMPLES "Enables Axom Examples" ON "ENABLE_EXAMPLES" OFF)
option(AXOM_ENABLE_TOOLS "Enables Axom Tools" ON)

option(AXOM_ENABLE_TUTORIALS "Builds Axom tutorials as part of the Axom build" ON)
mark_as_advanced(AXOM_ENABLE_TUTORIALS)

#------------------------------------------------------------------------------
# Test execution controls
#------------------------------------------------------------------------------
if(AXOM_ENABLE_OPENMP)
    set(AXOM_TEST_NUM_OMP_THREADS 4 CACHE STRING "Default number of OpenMP threads for tests")
else()
    set(AXOM_TEST_NUM_OMP_THREADS 0 CACHE STRING "Default number of OpenMP threads for tests")
endif()
mark_as_advanced(AXOM_TEST_NUM_OMP_THREADS)

#--------------------------------------------------------------------------
# Option to control whether AXOM_DEFINE compiler define is enabled
#
# Possible values are: "ON", "OFF" and "DEFAULT"
# By default, AXOM_DEBUG is defined in Debug and RelWithDebInfo configurations
#--------------------------------------------------------------------------
set(AXOM_DEBUG_DEFINE "DEFAULT" CACHE STRING "Controls whether AXOM_DEBUG compiler define is enabled")
set_property(CACHE AXOM_DEBUG_DEFINE PROPERTY STRINGS "DEFAULT" "ON" "OFF")

#------------------------------------------------------------------------------
# Option to gradually phase out deprecated types.
# With C++11, some Axom types in src/axom/core/Types.hpp are obsolete.
# They will be removed in steps, as the AXOM_DEPRECATED_TYPES variable
# defaults to WARN, then ERROR, then eventually removed.
#------------------------------------------------------------------------------
set(AXOM_DEPRECATED_TYPES "WARN" CACHE STRING "Controls deprecated types removal phase")
set_property(CACHE AXOM_DEPRECATED_TYPES PROPERTY STRINGS "WARN" "ERROR" "ALLOW")
