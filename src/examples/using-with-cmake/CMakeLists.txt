# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Example that shows how to use an installed instance of Axom 
# in a CMake-based build system.
#
# Configuration variables are stored in a CMake cache file 'host-config.cmake'
# which defines paths to axom and possibly other TPLs.
# It also contains information about the compiler used to build axom.
#------------------------------------------------------------------------------
#
# To build:
#  mkdir build
#  cd build
#  cmake -C ../host-config.cmake ..
#  make
#  ./example
#
#------------------------------------------------------------------------------

if (ENABLE_HIP OR AXOM_ENABLE_HIP OR ENABLE_CUDA OR AXOM_ENABLE_CUDA)
    cmake_minimum_required(VERSION 3.21)
else()
    cmake_minimum_required(VERSION 3.18)
endif()

project(using_with_cmake)

message(STATUS "CMake Version: ${CMAKE_VERSION}")

# _zerotoaxom_docs_start
#------------------------------------------------------------------------------
# Check for AXOM_DIR and use CMake's find_package to import axom's targets
#------------------------------------------------------------------------------
if(NOT DEFINED AXOM_DIR OR NOT EXISTS ${AXOM_DIR}/lib/cmake/axom-config.cmake)
    message(FATAL_ERROR "Missing required 'AXOM_DIR' variable pointing to an installed axom")
endif()

if (ENABLE_CUDA)
    enable_language(CUDA)
endif()

if (ENABLE_HIP)
    enable_language(HIP)
endif()

include(CMakeFindDependencyMacro)

find_dependency(axom REQUIRED
                NO_DEFAULT_PATH 
                PATHS ${AXOM_DIR}/lib/cmake)

# Remove implicitly added link directories added by CMake that are problematic when
# the default system libraries are older than the ones used by the compiler
if(BLT_CMAKE_IMPLICIT_LINK_DIRECTORIES_EXCLUDE)
    list(REMOVE_ITEM CMAKE_C_IMPLICIT_LINK_DIRECTORIES
                     ${BLT_CMAKE_IMPLICIT_LINK_DIRECTORIES_EXCLUDE})
    list(REMOVE_ITEM CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES
                     ${BLT_CMAKE_IMPLICIT_LINK_DIRECTORIES_EXCLUDE})
    list(REMOVE_ITEM CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES
                     ${BLT_CMAKE_IMPLICIT_LINK_DIRECTORIES_EXCLUDE})
endif()

#------------------------------------------------------------------------------
# Set up example target that depends on axom
#------------------------------------------------------------------------------
add_executable(example example.cpp)

# Compile as HIP source if enabled
if (ENABLE_HIP)
    set_source_files_properties(example.cpp PROPERTIES LANGUAGE HIP)
endif()

# setup the axom include path
target_include_directories(example PRIVATE ${AXOM_INCLUDE_DIRS})

# link to axom targets
target_link_libraries(example axom)
target_link_libraries(example axom::fmt)
# _zerotoaxom_docs_end
