# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Tutorial that demonstrates the usage of several Axom components
# in a BLT-based build system.
#
# Configuration variables are stored in a CMake cache file 'host-config.cmake'
# which defines paths to blt, axom and possibly other TPLs.
# It also contains information about the compiler used to build axom.
#------------------------------------------------------------------------------
#
# To build:
#  mkdir build
#  cd build
#  cmake -C ../host-config.cmake ..
#  make
#  ./bin/lesson_*
#
#------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.14)

project(shaping_tutorial)

#------------------------------------------------------------------------------
# Set up BLT with validity checks
#------------------------------------------------------------------------------

# Check that path to BLT is provided and valid
if(NOT DEFINED BLT_SOURCE_DIR OR NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
    message(FATAL_ERROR "Missing required 'BLT_SOURCE_DIR' variable pointing to a valid blt")
endif()
include(${BLT_SOURCE_DIR}/SetupBLT.cmake)


#------------------------------------------------------------------------------
# Check for AXOM_DIR and use CMake's find_package to import axom's targets
#------------------------------------------------------------------------------
if(NOT DEFINED AXOM_DIR OR NOT EXISTS ${AXOM_DIR}/lib/cmake/axom-config.cmake)
    message(FATAL_ERROR "Missing required 'AXOM_DIR' variable pointing to an installed axom")
endif()

include(CMakeFindDependencyMacro)

find_dependency(axom REQUIRED
                NO_DEFAULT_PATH 
                PATHS ${AXOM_DIR}/lib/cmake)

#------------------------------------------------------------------------------
# Set up target executables for tutorial lesson examples
#------------------------------------------------------------------------------
set(tutorial_deps axom axom::fmt axom::cli11)
blt_list_append(TO tutorial_deps ELEMENTS blt::openmp IF ENABLE_OPENMP)
blt_list_append(TO tutorial_deps ELEMENTS blt::cuda IF ENABLE_CUDA)
blt_list_append(TO tutorial_deps ELEMENTS blt::hip IF ENABLE_HIP)

set(_examples lesson_00:check_axom_configuration.cpp
              lesson_01:mesh_metadata_sidre.cpp 
              lesson_02:validated_inlet_metadata.cpp
              lesson_02:improved_inlet_metadata.cpp
              lesson_03:klee_operators_and_validation.cpp)
if(AXOM_USE_MFEM)
    blt_list_append(TO _examples ELEMENTS lesson_04:quest_sampling_shaper.cpp)
endif()

foreach(lesson_and_file ${_examples})
    string(REPLACE ":" ";" _parts "${lesson_and_file}")
    list(GET _parts 0 _lesson)
    list(GET _parts 1 _file_with_ext)
    get_filename_component( _file_no_ext ${_file_with_ext} NAME_WE )

    blt_add_executable(NAME        ${_lesson}_${_file_no_ext}
                       SOURCES     ${_lesson}/${_file_with_ext}
                       DEPENDS_ON  ${tutorial_deps})
endforeach()                   

#------------------------------------------------------------------------------
# Optionally, add tests to run the tutorial lessons
#------------------------------------------------------------------------------
if(ENABLE_TESTS)
    blt_add_test(NAME    lesson_00_check_axom_configuration
                 COMMAND lesson_00_check_axom_configuration)

    blt_add_test(NAME    lesson_01_mesh_metadata_sidre
                 COMMAND lesson_01_mesh_metadata_sidre  --min_x 50 --min_y 350 
                                                        --max_x 50 --max_y 350 
                                                        --res_x  8 --res_y   8 
                                                        --output_blueprint)

    blt_add_test(NAME    lesson_02_validated_inlet_metadata
                 COMMAND lesson_02_validated_inlet_metadata ../lesson_02/input2D.yaml)

    blt_add_test(NAME    lesson_02_improved_inlet_metadata_yaml
                COMMAND lesson_02_improved_inlet_metadata ../lesson_02/input3D.yaml)

    if(AXOM_USE_LUA)
      blt_add_test(NAME    lesson_02_improved_inlet_metadata_lua
                   COMMAND lesson_02_improved_inlet_metadata ../lesson_02/input3D.lua)
    endif()

    blt_add_test(NAME   lesson_03_klee_operators_and_validation
                COMMAND lesson_03_klee_operators_and_validation ../lesson_03/ice_cream.yaml)

    if(AXOM_USE_LUA AND AXOM_USE_MFEM)
      blt_add_test(NAME    lesson_04_quest_sampling_shaper
                   COMMAND lesson_04_quest_sampling_shaper -m ../lesson_04/circle_input.lua 
                                                           -k ../lesson_04/circles.yaml -v)                   
    endif()

endif()

#------------------------------------------------------------------------------
# Optionally, print out information about imported targets
#------------------------------------------------------------------------------
if(EXAMPLE_VERBOSE_OUTPUT)
    blt_print_target_properties(TARGET axom CHILDREN TRUE)  
endif()

