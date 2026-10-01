# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
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
# install 'using-with-python' example
#------------------------------------------------------------------------------
if(AXOM_ENABLE_EXAMPLES AND NANOBIND_FOUND AND AXOM_ENABLE_SIDRE)
    install(
        FILES
            using-with-python/README.md
            using-with-python/example.py
        DESTINATION
            examples/axom/using-with-python
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
    blt_list_append(TO radiuss_tutorial_deps ELEMENTS umpire::umpire IF umpire_FOUND)

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

    if(AXOM_ENABLE_TESTS)
        blt_add_test(NAME    shaping_tutorial_lesson_03_klee_operators_and_validation_yaml
                     COMMAND shaping_tutorial_lesson_03_klee_operators_and_validation
                             ${CMAKE_CURRENT_SOURCE_DIR}/shaping_tutorial/lesson_03/ice_cream.yaml)

        if(LUA_FOUND)
            blt_add_test(NAME    shaping_tutorial_lesson_03_klee_operators_and_validation_lua
                         COMMAND shaping_tutorial_lesson_03_klee_operators_and_validation
                                 ${CMAKE_CURRENT_SOURCE_DIR}/shaping_tutorial/lesson_03/ice_cream.lua)

            blt_add_test(NAME    shaping_tutorial_lesson_03_klee_operators_and_validation_lua_initialization
                         COMMAND shaping_tutorial_lesson_03_klee_operators_and_validation
                                 ${CMAKE_CURRENT_SOURCE_DIR}/shaping_tutorial/lesson_03/ice_cream_initialized.lua
                                 --initialization-file
                                 ${CMAKE_CURRENT_SOURCE_DIR}/shaping_tutorial/lesson_03/ice_cream_initialization.lua)

            if(MFEM_FOUND)
                set(_lesson_04_dir "${CMAKE_CURRENT_SOURCE_DIR}/shaping_tutorial/lesson_04")

                set(_lesson_04_yaml_test shaping_tutorial_lesson_04_quest_sampling_shaper_yaml)
                set(_lesson_04_lua_callbacks_test shaping_tutorial_lesson_04_quest_sampling_shaper_lua_callbacks)
                set(_lesson_04_lua_initialization_test shaping_tutorial_lesson_04_quest_sampling_shaper_lua_initialization)
                foreach(_testname
                        ${_lesson_04_yaml_test}
                        ${_lesson_04_lua_callbacks_test}
                        ${_lesson_04_lua_initialization_test})
                    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${_testname}")
                endforeach()

                blt_add_test(NAME    ${_lesson_04_yaml_test}
                             COMMAND shaping_tutorial_lesson_04_quest_sampling_shaper
                                     -m ${_lesson_04_dir}/circle_input.lua
                                     -k ${_lesson_04_dir}/circles.yaml
                             WORKING_DIRECTORY
                                     ${CMAKE_CURRENT_BINARY_DIR}/${_lesson_04_yaml_test})

                blt_add_test(NAME    ${_lesson_04_lua_callbacks_test}
                             COMMAND shaping_tutorial_lesson_04_quest_sampling_shaper
                                     -m ${_lesson_04_dir}/circle_input.lua
                                     -k ${_lesson_04_dir}/circles.lua
                             WORKING_DIRECTORY
                                     ${CMAKE_CURRENT_BINARY_DIR}/${_lesson_04_lua_callbacks_test})

                blt_add_test(NAME    ${_lesson_04_lua_initialization_test}
                             COMMAND shaping_tutorial_lesson_04_quest_sampling_shaper
                                     -m ${_lesson_04_dir}/circle_input.lua
                                     -k ${_lesson_04_dir}/circles_initialized.lua
                                     --lua-init-file ${_lesson_04_dir}/circles_initialization.lua
                             WORKING_DIRECTORY
                                     ${CMAKE_CURRENT_BINARY_DIR}/${_lesson_04_lua_initialization_test})

                set(_testname
                    shaping_tutorial_lesson_04_quest_sampling_shaper_yaml_rejects_lua_initialization)
                blt_add_test(NAME    ${_testname}
                             COMMAND shaping_tutorial_lesson_04_quest_sampling_shaper
                                     -m ${_lesson_04_dir}/circle_input.lua
                                     -k ${_lesson_04_dir}/circles.yaml
                                     --lua-init-file ${_lesson_04_dir}/circles_initialization.lua)
                set_tests_properties(${_testname} PROPERTIES
                                     PASS_REGULAR_EXPRESSION "only supported for Lua input decks")
            endif()
        endif()
    endif()
    
endif()
