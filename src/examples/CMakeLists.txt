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


#------------------------------------------------------------------------------
# configure and install 'shaping_tutorial' example
#------------------------------------------------------------------------------
configure_file(
    shaping_tutorial/host-config.cmake.in
    ${PROJECT_BINARY_DIR}/examples/shaping_tutorial/host-config.cmake)

install(
    FILES
        shaping_tutorial/CMakeLists.txt
        ${PROJECT_BINARY_DIR}/examples/shaping_tutorial/host-config.cmake
    DESTINATION
        examples/axom/shaping_tutorial
)

foreach(_lesson lesson_00 lesson_01 lesson_02)
    install(DIRECTORY shaping_tutorial/${_lesson}
            DESTINATION examples/axom/shaping_tutorial)
endforeach()

set(AXOM_BUILD_TUTORIALS ON)
if(AXOM_BUILD_TUTORIALS)
    set(shaping_tutorial_deps axom fmt cli11)
    blt_list_append(TO shaping_tutorial_deps ELEMENTS blt::openmp IF ENABLE_OPENMP)
    blt_list_append(TO shaping_tutorial_deps ELEMENTS blt::cuda IF ENABLE_CUDA)
    blt_list_append(TO shaping_tutorial_deps ELEMENTS blt::hip IF ENABLE_HIP)

    foreach(lesson_and_file lesson_01:mesh_metadata_sidre.cpp 
                            lesson_02:validated_inlet_metadata.cpp
                            lesson_02:improved_inlet_metadata.cpp)
        string(REPLACE ":" ";" _parts "${lesson_and_file}")
        list(GET _parts 0 _lesson)
        list(GET _parts 1 _file_with_ext)
        get_filename_component( _file_no_ext ${_file_with_ext} NAME_WE )

        blt_add_executable(NAME       shaping_tutorial_${_lesson}_${_file_no_ext}
                           SOURCES    shaping_tutorial/${_lesson}/${_file_with_ext}
                           DEPENDS_ON ${shaping_tutorial_deps})

    endforeach()
    
endif()
