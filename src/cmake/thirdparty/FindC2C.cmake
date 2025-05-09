# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Setup C2C
#------------------------------------------------------------------------------
# This file defines:
#  C2C_FOUND - If c2c was found
#  C2C_INCLUDE_DIRS - The c2c include directories
#  C2C_LIBRARY - The c2c library
#------------------------------------------------------------------------------

# First check for C2C_DIR
if (NOT EXISTS "${C2C_DIR}")
    message(FATAL_ERROR "Given C2C_DIR does not exist: ${C2C_DIR}")
endif()

if (NOT IS_DIRECTORY "${C2C_DIR}")
    message(FATAL_ERROR "Given C2C_DIR is not a directory: ${C2C_DIR}")
endif()

# Find includes directory
find_path( C2C_INCLUDE_DIR c2c/C2C.hpp
           PATHS  ${C2C_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

# Find libraries
find_library( C2C_LIBRARY NAMES c2c libc2c
              PATHS ${C2C_DIR}/lib
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set C2C_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(C2C  DEFAULT_MSG
                                  C2C_INCLUDE_DIR
                                  C2C_LIBRARY )

if(NOT C2C_FOUND)
    message(FATAL_ERROR "C2C_DIR is not a path to a valid c2c install: ${C2C_DIR}")
endif()

message(STATUS "c2c includes: ${C2C_INCLUDE_DIR}")
message(STATUS "c2c libraries: ${C2C_LIBRARY}")
