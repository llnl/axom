# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Setup ZFPoly
#------------------------------------------------------------------------------
# This file defines:
#  ZFPOLY_FOUND - If zfpoly was found
#  ZFPOLY_INCLUDE_DIRS - The zfpoly include directories
#  ZFPOLY_LIBRARY - The zfpoly library
#------------------------------------------------------------------------------

# ZFPOLY_DIR is asserted to be an existing directory by the caller
# (see SetupAxomThirdParty.cmake: axom_assert_is_directory).

find_path( ZFPOLY_INCLUDE_DIRS zfpoly.h
           PATHS  ${ZFPOLY_DIR}/include
                  ${ZFPOLY_DIR}
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

# Find libraries
find_library( ZFPOLY_LIBRARY NAMES zfpoly libzfpoly
              PATHS ${ZFPOLY_DIR}/lib
                    ${ZFPOLY_DIR}/lib64
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set ZFPOLY_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(ZFPoly  DEFAULT_MSG
                                  ZFPOLY_INCLUDE_DIRS
                                  ZFPOLY_LIBRARY )

if(NOT ZFPOLY_FOUND)
    message(FATAL_ERROR "ZFPOLY_DIR is not a path to a valid zfpoly install: ${ZFPOLY_DIR}")
endif()

message(STATUS "ZFPoly includes: ${ZFPOLY_INCLUDE_DIRS}")
message(STATUS "ZFPoly library: ${ZFPOLY_LIBRARY}")
