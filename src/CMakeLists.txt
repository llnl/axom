# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Axom project
#------------------------------------------------------------------------------

if (ENABLE_HIP OR AXOM_ENABLE_HIP)
    cmake_minimum_required(VERSION 3.21)
else()
    # Do not bump this over 3.14 due to it changing CMake policies and breaking
    # the CUDA build. We check the required version specifically for this below
    cmake_minimum_required(VERSION 3.14)
endif()

project(axom LANGUAGES C CXX)

if (ENABLE_FORTRAN)
   enable_language(Fortran)
endif()

#------------------------------------------------------------------------------
# Initialize BLT build system
#------------------------------------------------------------------------------
if (DEFINED BLT_SOURCE_DIR)
    # Support having a shared BLT outside of the repository if given a BLT_SOURCE_DIR

    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
        message(FATAL_ERROR "Given BLT_SOURCE_DIR does not contain SetupBLT.cmake")
    endif()
else()
    # Use internal 'blt' submodule path if BLT_SOURCE_DIR not provided
    set(BLT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/cmake/blt" CACHE PATH "")
    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
        message(FATAL_ERROR
            "Cannot locate BLT. "
            "Either run the following command in your git repository: \n"
            "    git submodule update --init --recursive\n"
            "Or add -DBLT_SOURCE_DIR=/path/to/blt to your CMake command." )
    endif()
endif()

if (“${PROJECT_SOURCE_DIR}” STREQUAL “${CMAKE_SOURCE_DIR}”)
    # Set some default BLT options before loading BLT only if not included in
    # another project
    if (NOT BLT_CXX_STD)
        set(BLT_CXX_STD "c++17" CACHE STRING "")
    endif()

    # These are not used in Axom, turn them off
    set(_unused_blt_tools
        CLANGQUERY
        VALGRIND
        ASTYLE
        CMAKEFORMAT
        UNCRUSTIFY
        YAPF)
    foreach(_tool ${_unused_blt_tools})
        set(ENABLE_${_tool} OFF CACHE BOOL "")
    endforeach()

    # These are only used by Axom developers, so turn them off
    # unless an explicit executable path is given
    set(_used_blt_tools
        CLANGFORMAT
        CLANGTIDY
        CLANGAPPLYREPLACEMENTS
        CPPCHECK
        DOXYGEN
        SPHINX)
    foreach(_tool ${_used_blt_tools})
        if(NOT ${_tool}_EXECUTABLE)
            set(ENABLE_${_tool} OFF CACHE BOOL "")
        else()
            set(ENABLE_${_tool} ON CACHE BOOL "")
        endif()
    endforeach()

    set(BLT_REQUIRED_CLANGFORMAT_VERSION  "19" CACHE STRING "")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/.clang-tidy 
                   ${CMAKE_CURRENT_BINARY_DIR}/.clang-tidy
                   COPYONLY)

    # If Axom is the top project and AXOM_ENABLE_TESTS is off, force ENABLE_TESTS to off so
    # gtest doesn't build when it's not needed
    if(DEFINED AXOM_ENABLE_TESTS AND NOT AXOM_ENABLE_TESTS)
        set(ENABLE_TESTS OFF CACHE BOOL "")
    endif()

    # If either AXOM_ENABLE_TESTS or ENABLE_TEST are explicitly turned off by the user,
    # turn off GMock, otherwise turn it on
    if((DEFINED AXOM_ENABLE_TESTS AND NOT AXOM_ENABLE_TESTS) OR
       (DEFINED ENABLE_TESTS AND NOT ENABLE_TESTS))
        set(ENABLE_GMOCK OFF CACHE BOOL "")
    else()
        set(ENABLE_GMOCK ON CACHE BOOL "")
    endif()

    # We use BLT's install targets logic
    set(BLT_EXPORT_THIRDPARTY OFF CACHE BOOL "" FORCE)
endif()

if("${BLT_CXX_STD}" STREQUAL "c++98" OR "${BLT_CXX_STD}" STREQUAL "c++11"
                                     OR "${BLT_CXX_STD}" STREQUAL "c++14" )
    message(FATAL_ERROR "Axom requires BLT_CXX_STD to be 'c++17' or above.")
endif()

include(${BLT_SOURCE_DIR}/SetupBLT.cmake)

#------------------------------------------------------------------------------
# Attempt to set Axom's data directory (used for Axom's tests and examples)
#------------------------------------------------------------------------------
if(NOT AXOM_DATA_DIR)
    # Use internal 'axom_dir' submodule path if AXOM_DATA_DIR not provided
    get_filename_component(_data_dir "${PROJECT_SOURCE_DIR}/../data" ABSOLUTE)
    if(EXISTS ${_data_dir}/README.md)
        set(AXOM_DATA_DIR ${_data_dir} CACHE PATH "")
    endif()
endif()

#------------------------------------------------------------------------------
# Include build system logic and options
#------------------------------------------------------------------------------
include(CMakeDependentOption)
include(cmake/CMakeBasics.cmake)

# Check for minimum CMake version required w/o changing policies like cmake_minimum_required
if(AXOM_ENABLE_CUDA AND ${CMAKE_VERSION} VERSION_LESS 3.18.0)
    message(FATAL_ERROR "Axom requires CMake version 3.18.0+ when CUDA is enabled.")
endif()

axom_add_code_checks()

#------------------------------------------------------------------------------
# Add source directories
#------------------------------------------------------------------------------
add_subdirectory(thirdparty)

include(cmake/AxomVersion.cmake)
add_subdirectory(axom)

if(AXOM_ENABLE_TOOLS)
    add_subdirectory(tools)
endif()

# Installation examples ('Using Axom') and/or tutorials
if(AXOM_ENABLE_EXAMPLES OR AXOM_ENABLE_TUTORIALS)
    add_subdirectory(examples)
endif()

if(AXOM_ENABLE_DOCS)
    if(SPHINX_FOUND)
        blt_add_sphinx_target( axom_docs )
    endif()
    add_subdirectory(docs)
endif()

#------------------------------------------------------------------------------
# Generate header file with configuration options
#------------------------------------------------------------------------------
include(cmake/AxomConfig.cmake)
