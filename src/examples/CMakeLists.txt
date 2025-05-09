# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Examples for using Axom within other build systems: Make, CMake, BLT
#------------------------------------------------------------------------------
# Note: Installed to 'examples/axom' to avoid name conflicts with other libs.
# This will crop up if we use spack to install axom as a python extension.
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# install 'using-with-make' example
#------------------------------------------------------------------------------
# Create a link line string that has all the component libraries in the proper order
# This list is ordered by hand, then we only add components that are enabled
set(_ordered_component_list quest multimat slam mint klee inlet sidre slic lumberjack core)
set(AXOM_ORDERED_LIBS)
foreach(_component ${_ordered_component_list})
  if(_component IN_LIST AXOM_COMPONENTS_ENABLED)
    set(AXOM_ORDERED_LIBS "${AXOM_ORDERED_LIBS} -laxom_${_component}")
  endif()
endforeach()

configure_file(
    using-with-make/Makefile.in
    ${PROJECT_BINARY_DIR}/examples/using-with-make/Makefile)

install(
    FILES
        ${PROJECT_BINARY_DIR}/examples/using-with-make/Makefile
        using-with-make/example.cpp
    DESTINATION 
        examples/axom/using-with-make
)

#------------------------------------------------------------------------------
# configure and install 'using-with-cmake' example
#------------------------------------------------------------------------------
configure_file(
    using-with-cmake/host-config.cmake.in
    ${PROJECT_BINARY_DIR}/examples/using-with-cmake/host-config.cmake)

install(
    FILES
        using-with-cmake/CMakeLists.txt
        using-with-cmake/example.cpp
        ${PROJECT_BINARY_DIR}/examples/using-with-cmake/host-config.cmake
    DESTINATION
        examples/axom/using-with-cmake
)

#------------------------------------------------------------------------------
# configure and install 'using-with-blt' example
#------------------------------------------------------------------------------
configure_file(
    using-with-blt/host-config.cmake.in
    ${PROJECT_BINARY_DIR}/examples/using-with-blt/host-config.cmake)

install(
    FILES
        using-with-blt/CMakeLists.txt
        using-with-blt/example.cpp
        ${PROJECT_BINARY_DIR}/examples/using-with-blt/host-config.cmake
    DESTINATION
        examples/axom/using-with-blt
)


#------------------------------------------------------------------------------
# configure and install 'radiuss_tutorial' example
#------------------------------------------------------------------------------
configure_file(
    radiuss_tutorial/host-config.cmake.in
    ${PROJECT_BINARY_DIR}/examples/radiuss_tutorial/host-config.cmake)

install(
    FILES
        radiuss_tutorial/CMakeLists.txt
        ${PROJECT_BINARY_DIR}/examples/radiuss_tutorial/host-config.cmake
    DESTINATION
        examples/axom/radiuss_tutorial
)

foreach(_lesson lesson_00 lesson_01 lesson_02 lesson_03 lesson_04 patch figs)
    install(DIRECTORY radiuss_tutorial/${_lesson}
            DESTINATION examples/axom/radiuss_tutorial)
endforeach()

