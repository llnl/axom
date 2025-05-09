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

project(radiuss_tutorial)

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
# Set up target executables for the tutorial
#------------------------------------------------------------------------------
set(tutorial_deps axom axom::fmt axom::cli11)
blt_list_append(TO tutorial_deps ELEMENTS blt::openmp IF ENABLE_OPENMP)
blt_list_append(TO tutorial_deps ELEMENTS blt::cuda IF ENABLE_CUDA)
blt_list_append(TO tutorial_deps ELEMENTS blt::hip IF ENABLE_HIP)

blt_add_executable(NAME       lesson_00_check_axom_configuration 
                   SOURCES    lesson_00/check_axom_configuration.cpp
                   DEPENDS_ON ${tutorial_deps})

blt_add_executable(NAME       lesson_01_load_stl_mesh 
                   SOURCES    lesson_01/load_stl_mesh.cpp
                   DEPENDS_ON ${tutorial_deps})

blt_add_executable(NAME       lesson_02_naive_self_intersections
                   SOURCES    lesson_02/naive_self_intersections.cpp
                   DEPENDS_ON ${tutorial_deps})

if(umpire_FOUND AND RAJA_FOUND)
    blt_list_append(TO tutorial_deps ELEMENTS RAJA   IF RAJA_FOUND)
    blt_list_append(TO tutorial_deps ELEMENTS umpire IF umpire_FOUND)

    blt_add_executable(NAME       lesson_03_device_self_intersections
                       SOURCES    lesson_03/device_self_intersections.cpp
                       DEPENDS_ON ${tutorial_deps})

    blt_add_executable(NAME       lesson_04_device_spatial_indexes
                       SOURCES    lesson_04/device_spatial_indexes.cpp
                       DEPENDS_ON ${tutorial_deps})
endif()

#------------------------------------------------------------------------------
# Optionally, add tests to run the tutorial lessons
#------------------------------------------------------------------------------
if(ENABLE_TESTS)
    blt_add_test(NAME    lesson_00_check_axom_configuration
                 COMMAND lesson_00_check_axom_configuration)

    set(_sphere_mesh "${AXOM_DATA_DIR}/quest/sphere_binary.stl")
    set(_plane_problem_mesh "${AXOM_DATA_DIR}/quest/plane_simp_problems.stl")

    if(EXISTS ${_sphere_mesh})
        blt_add_test(
            NAME    lesson_01_load_stl_mesh_sphere
            COMMAND lesson_01_load_stl_mesh -i ${_sphere_mesh})

        blt_add_test(
            NAME    lesson_02_naive_self_intersections_sphere
            COMMAND lesson_02_naive_self_intersections -i ${_sphere_mesh} --use-bounding-boxes)
    endif()

    set (_policies "seq")
    blt_list_append(TO _policies ELEMENTS "omp"  IF ENABLE_OPENMP)
    blt_list_append(TO _policies ELEMENTS "cuda" IF ENABLE_CUDA)
    blt_list_append(TO _policies ELEMENTS "hip"  IF ENABLE_HIP)

    # the plane_mesh can take too long on lesson 3 in Debug configs, so only use it in Release configs
    set(_lesson_03_mesh_name "sphere")
    set(_lesson_03_mesh_path ${_sphere_mesh})
    if(${CMAKE_BUILD_TYPE} MATCHES "^[Rr]elease")
        set(_lesson_03_mesh_name "plane")
        set(_lesson_03_mesh_path ${_plane_problem_mesh})
    endif()
    if(TARGET lesson_03_device_self_intersections AND EXISTS ${_lesson_03_mesh_path})
        foreach(_pol ${_policies})
            blt_add_test(
                NAME    lesson_03_device_self_intersections_${_lesson_03_mesh_name}_${_pol}
                COMMAND lesson_03_device_self_intersections -i ${_lesson_03_mesh_path} --use-bounding-boxes -p ${_pol})
        endforeach()
    endif()

    if(TARGET lesson_04_device_spatial_indexes AND EXISTS ${_plane_problem_mesh})
        foreach(_pol ${_policies})
            blt_add_test(
                NAME    lesson_04_device_spatial_indexes_plane_${_pol}
                COMMAND lesson_04_device_spatial_indexes -i ${_plane_problem_mesh} -p ${_pol})
        endforeach()
    endif()

endif()


#------------------------------------------------------------------------------
# Optionally, print out information about imported targets
#------------------------------------------------------------------------------
if(EXAMPLE_VERBOSE_OUTPUT)
    blt_print_target_properties(TARGET axom CHILDREN TRUE)  
endif()

