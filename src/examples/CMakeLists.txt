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
#------------------------------------------------------------------------------
if (AXOM_ENABLE_EXAMPLES)
    set(_ordered_component_list quest mir
                                bump  klee
                                inlet mint  multimat
                                sidre slam
                                slic 
                                sina  lumberjack
                                core)
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
endif()
#------------------------------------------------------------------------------
# configure and install 'using-with-cmake' example
#------------------------------------------------------------------------------
if (AXOM_ENABLE_EXAMPLES)
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
endif()    

#------------------------------------------------------------------------------
# configure and install 'using-with-blt' example
#------------------------------------------------------------------------------
if (AXOM_ENABLE_EXAMPLES)
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
endif()

#------------------------------------------------------------------------------
# configure and install 'radiuss_tutorial' example
# This example requires Quest (which requires Slic, Mint, Primal and Spin)
#------------------------------------------------------------------------------
if(AXOM_ENABLE_TUTORIALS AND AXOM_ENABLE_QUEST)
    configure_file(
        radiuss_tutorial/host-config.cmake.in
        ${PROJECT_BINARY_DIR}/tutorials/radiuss_tutorial/host-config.cmake)

    install(
        FILES
            radiuss_tutorial/CMakeLists.txt
            ${PROJECT_BINARY_DIR}/tutorials/radiuss_tutorial/host-config.cmake
        DESTINATION
            examples/axom/radiuss_tutorial
    )

    foreach(_lesson lesson_00 lesson_01 lesson_02 lesson_03 lesson_04 patch figs)
        install(DIRECTORY radiuss_tutorial/${_lesson}
                DESTINATION examples/axom/radiuss_tutorial)
    endforeach()

    # Also build the tutorial examples; some require RAJA and Umpire
    set(radiuss_tutorial_deps axom fmt cli11)
    blt_list_append(TO radiuss_tutorial_deps ELEMENTS blt::openmp IF ENABLE_OPENMP)
    blt_list_append(TO radiuss_tutorial_deps ELEMENTS blt::cuda IF ENABLE_CUDA)
    blt_list_append(TO radiuss_tutorial_deps ELEMENTS blt::hip IF ENABLE_HIP)
    blt_list_append(TO radiuss_tutorial_deps ELEMENTS RAJA IF RAJA_FOUND)
    blt_list_append(TO radiuss_tutorial_deps ELEMENTS umpire IF umpire_FOUND)

    set(_examples lesson_00:check_axom_configuration.cpp
                  lesson_01:load_stl_mesh.cpp
                  lesson_02:naive_self_intersections.cpp)
    if(RAJA_FOUND AND umpire_FOUND)
        blt_list_append(TO _examples ELEMENTS lesson_03:device_self_intersections.cpp
                                              lesson_04:device_spatial_indexes.cpp)
    endif()

    foreach(lesson_and_file ${_examples})
        string(REPLACE ":" ";" _parts "${lesson_and_file}")
        list(GET _parts 0 _lesson)
        list(GET _parts 1 _file_with_ext)
        get_filename_component( _file_no_ext ${_file_with_ext} NAME_WE )

        blt_add_executable(NAME       radiuss_tutorial_${_lesson}_${_file_no_ext}
                           SOURCES    radiuss_tutorial/${_lesson}/${_file_with_ext}
                           OUTPUT_DIR ${PROJECT_BINARY_DIR}/tutorials
                           DEPENDS_ON ${radiuss_tutorial_deps})
    endforeach()
endif()

#------------------------------------------------------------------------------
# configure and install 'shaping_tutorial' example
# This example requires Klee (which requires Inlet, Sidre and Conduit) as well as Quest
#------------------------------------------------------------------------------
if(AXOM_ENABLE_TUTORIALS AND AXOM_ENABLE_KLEE AND AXOM_ENABLE_QUEST)
    configure_file(
        shaping_tutorial/host-config.cmake.in
        ${PROJECT_BINARY_DIR}/tutorials/shaping_tutorial/host-config.cmake)

    install(
        FILES
            shaping_tutorial/CMakeLists.txt
            ${PROJECT_BINARY_DIR}/tutorials/shaping_tutorial/host-config.cmake
        DESTINATION
            examples/axom/shaping_tutorial
    )

    foreach(_lesson lesson_00 lesson_01 lesson_02 lesson_03 lesson_04)
        install(DIRECTORY shaping_tutorial/${_lesson}
                DESTINATION examples/axom/shaping_tutorial)
    endforeach()

    # Also build the tutorial examples
    set(shaping_tutorial_deps axom fmt cli11)
    blt_list_append(TO shaping_tutorial_deps ELEMENTS blt::openmp IF ENABLE_OPENMP)
    blt_list_append(TO shaping_tutorial_deps ELEMENTS blt::cuda IF ENABLE_CUDA)
    blt_list_append(TO shaping_tutorial_deps ELEMENTS blt::hip IF ENABLE_HIP)
    blt_list_append(TO shaping_tutorial_deps ELEMENTS mfem IF MFEM_FOUND)
    
    set(_examples lesson_00:check_axom_configuration.cpp
                  lesson_01:mesh_metadata_sidre.cpp 
                  lesson_02:validated_inlet_metadata.cpp
                  lesson_02:improved_inlet_metadata.cpp
                  lesson_03:klee_operators_and_validation.cpp)
    if(LUA_FOUND AND MFEM_FOUND)
        blt_list_append(TO _examples ELEMENTS lesson_04:quest_sampling_shaper.cpp)
    endif()

    foreach(lesson_and_file ${_examples})
        string(REPLACE ":" ";" _parts "${lesson_and_file}")
        list(GET _parts 0 _lesson)
        list(GET _parts 1 _file_with_ext)
        get_filename_component( _file_no_ext ${_file_with_ext} NAME_WE )

        blt_add_executable(NAME       shaping_tutorial_${_lesson}_${_file_no_ext}
                           SOURCES    shaping_tutorial/${_lesson}/${_file_with_ext}
                           OUTPUT_DIR ${PROJECT_BINARY_DIR}/tutorials
                           DEPENDS_ON ${shaping_tutorial_deps})
    endforeach()
    
endif()
