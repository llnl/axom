# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Setup SCR
#------------------------------------------------------------------------------
# This file defines:
#  SCR_FOUND - If SCR was found
#  SCR_INCLUDE_DIRS - The SCR include directories
#  SCR_LIBRARIES - The SCR library
#------------------------------------------------------------------------------

# first Check for SCR_DIR

if(NOT SCR_DIR)
    message(FATAL_ERROR "Could not find SCR. SCR support needs explicit SCR_DIR")
endif()

#find includes
find_path( SCR_INCLUDE_DIR scr.h
           PATHS  ${SCR_DIR}/include/
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

# Add SCR's dependency's include
find_path( KVTREE_INCLUDE_DIR kvtree.h
           PATHS  ${KVTREE_DIR}/include/
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_path( DTCMP_INCLUDE_DIR dtcmp.h
           PATHS  ${DTCMP_DIR}/include/
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_path( SPATH_INCLUDE_DIR spath.h
           PATHS  ${SPATH_DIR}/include/
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_path( AXL_INCLUDE_DIR axl.h
           PATHS  ${AXL_DIR}/include/
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_path( LWGRP_INCLUDE_DIR lwgrp.h
           PATHS  ${LWGRP_DIR}/include/
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_path( ER_INCLUDE_DIR er.h
           PATHS  ${ER_DIR}/include/
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_path( RANKSTR_INCLUDE_DIR rankstr_mpi.h
           PATHS  ${RANKSTR_DIR}/include/
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_path( REDSET_INCLUDE_DIR redset.h
           PATHS  ${REDSET_DIR}/include/
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_path( SHUFFILE_INCLUDE_DIR shuffile.h
           PATHS  ${SHUFFILE_DIR}/include/
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_path( LIBYOGRT_INCLUDE_DIR yogrt.h
           PATHS  ${LIBYOGRT_DIR}/include/
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

set(SCR_INCLUDE_DIRS
    ${SCR_INCLUDE_DIR}
    ${KVTREE_INCLUDE_DIR}
    ${DTCMP_INCLUDE_DIR}
    ${SPATH_INCLUDE_DIR}
    ${AXL_INCLUDE_DIR}
    ${LWGRP_INCLUDE_DIR}
    ${ER_INCLUDE_DIR}
    ${RANKSTR_INCLUDE_DIR}
    ${REDSET_INCLUDE_DIR}
    ${SHUFFILE_INCLUDE_DIR}
    ${LIBYOGRT_INCLUDE_DIR}
    )

set(_library_names
    scr
    kvtree
    dtcmp
    spath
    axl
    lwgrp
    er
    rankstr
    redset
    redset_base
    shuffile
    yogrt
    )

set(_library_paths
   ${SCR_DIR}/lib
   ${KVTREE_DIR}/lib
   ${DTCMP_DIR}/lib
   ${SPATH_DIR}/lib
   ${AXL_DIR}/lib
   ${LWGRP_DIR}/lib
   ${ER_DIR}/lib
   ${RANKSTR_DIR}/lib
   ${REDSET_DIR}/lib
   ${SHUFFILE_DIR}/lib
   ${LIBYOGRT_DIR}/lib
   )

blt_find_libraries(
        FOUND_LIBS SCR_LIBRARIES
        NAMES      ${_library_names}
        REQUIRED   TRUE
        PATHS      ${_library_paths})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set SCR_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(SCR  DEFAULT_MSG
                                  SCR_INCLUDE_DIRS
                                  SCR_LIBRARIES )
