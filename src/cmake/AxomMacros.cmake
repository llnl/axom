# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

##------------------------------------------------------------------------------
## axom_add_code_checks()
##
## Adds code checks for all source files recursively in the Axom repository.
## 
## This creates the following parent build targets:
##  check - Runs a non file changing style check and CppCheck
##  style - In-place code formatting
##
## Creates various child build targets that follow this pattern:
##  axom_<check|style>
##  axom_<cppcheck|clangformat>_<check|style>
##------------------------------------------------------------------------------
macro(axom_add_code_checks)

    set(options)
    set(singleValueArgs)
    set(multiValueArgs)

    # Parse the arguments to the macro
    cmake_parse_arguments(arg
         "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

    # Only do code checks if building Axom by itself and not included in
    # another project
    if ("${PROJECT_SOURCE_DIR}" STREQUAL "${CMAKE_SOURCE_DIR}")
        # Create file globbing expressions that only include directories that contain source
        set(_base_dirs "axom" "examples" "thirdparty/tests" "tools")
        set(_ext_expressions "*.cpp" "*.hpp" "*.inl"
                             "*.cxx" "*.hxx" "*.cc" "*.c" "*.h" "*.hh"
                             "*.F" "*.f" "*.f90" "*.F90")

        set(_glob_expressions)
        foreach(_exp ${_ext_expressions})
            foreach(_base_dir ${_base_dirs})
                list(APPEND _glob_expressions "${PROJECT_SOURCE_DIR}/${_base_dir}/${_exp}")
            endforeach()
        endforeach()

        # Glob for list of files to run code checks on
        set(_sources)
        file(GLOB_RECURSE _sources ${_glob_expressions})

        # Filter out exclusions
        set(_exclude_expressions
            "${PROJECT_SOURCE_DIR}/axom/sidre/examples/lulesh2/*"
            "${PROJECT_SOURCE_DIR}/axom/slam/examples/lulesh2.0.3/*"
            "${PROJECT_SOURCE_DIR}/axom/slam/examples/tinyHydro/*")
        foreach(_exp ${_exclude_expressions})
            list(FILTER _sources EXCLUDE REGEX ${_exp})
        endforeach()

        blt_add_code_checks(PREFIX          axom
                            SOURCES         ${_sources}
                            CLANGFORMAT_CFG_FILE ${PROJECT_SOURCE_DIR}/.clang-format
                            CPPCHECK_FLAGS  --enable=all --inconclusive)

        # Set FOLDER property for code check targets
        foreach(_suffix clangformat_check clangformat_style clang_tidy_check clang_tidy_style)
            set(_tgt ${arg_PREFIX}_${_suffix})
            if(TARGET ${_tgt}) 
                set_target_properties(${_tgt} PROPERTIES FOLDER "axom/code_checks")
            endif()
        endforeach()
    endif()

endmacro(axom_add_code_checks)


##------------------------------------------------------------------------------
## axom_add_component( COMPONENT_NAME <name> DEFAULT_STATE [ON/OFF] )
##
## Adds a project component to the build.
##
## Adds a component to the build given the component's name and default state
## (ON/OFF). This macro also adds an "option" so that the user can control,
## which components to build.
##------------------------------------------------------------------------------
set(AXOM_COMPONENTS_FULL    CACHE STRING "List of all components in Axom" FORCE)
set(AXOM_COMPONENTS_ENABLED CACHE STRING "List of all enabled components in Axom" FORCE)
macro(axom_add_component)

    set(options)
    set(singleValueArgs COMPONENT_NAME DEFAULT_STATE )
    set(multiValueArgs)

    # Parse the arguments to the macro
    cmake_parse_arguments(arg
         "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

    # Adds an option so that the user can control whether to build this
    # component.
    # convert the component name to capitals for the AXOM_ENABLE options.
    string(TOUPPER ${arg_COMPONENT_NAME} COMPONENT_NAME_CAPITALIZED)
    string(TOLOWER ${arg_COMPONENT_NAME} COMPONENT_NAME_LOWERED)

    if(NOT "${COMPONENT_NAME_LOWERED}" STREQUAL "core")
        option(AXOM_ENABLE_${COMPONENT_NAME_CAPITALIZED}
              "Enables ${arg_COMPONENT_NAME}"
              ${arg_DEFAULT_STATE})
    endif()

    set(AXOM_COMPONENTS_FULL ${AXOM_COMPONENTS_FULL} ${COMPONENT_NAME_LOWERED}
        CACHE STRING "List of all components in Axom" FORCE)

    if ( AXOM_ENABLE_${COMPONENT_NAME_CAPITALIZED} OR "${COMPONENT_NAME_LOWERED}" STREQUAL "core")
        set(AXOM_COMPONENTS_ENABLED ${AXOM_COMPONENTS_ENABLED} ${COMPONENT_NAME_LOWERED}
            CACHE STRING "List of all enabled components in Axom" FORCE)
        add_subdirectory( ${arg_COMPONENT_NAME} )
    endif()

    unset(COMPONENT_NAME_CAPITALIZED)
    unset(COMPONENT_NAME_LOWERED)
endmacro(axom_add_component)


##------------------------------------------------------------------------------
## axom_add_executable(
##                     NAME        <name>
##                     SOURCES     [source1 [source2 ...]]
##                     HEADERS     [header1 [header2 ...]]
##                     INCLUDES    [dir1 [dir2 ...]]
##                     DEFINES     [define1 [define2 ...]]
##                     DEPENDS_ON  [dep1 [dep2 ...]]
##                     OUTPUT_DIR  [dir]
##                     OUTPUT_NAME [name]
##                     FOLDER      [name])
##
## Wrapper around blt_add_executable
##------------------------------------------------------------------------------
macro(axom_add_executable)

    set(options )
    set(singleValueArgs NAME OUTPUT_DIR OUTPUT_NAME FOLDER)
    set(multiValueArgs HEADERS SOURCES INCLUDES DEFINES DEPENDS_ON)

    # Parse the arguments to the macro
    cmake_parse_arguments(arg
        "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

    # Blanket add openmp as a dependency to get around not propagating
    # openmp to fortran dependencies
    blt_list_append(TO arg_DEPENDS_ON ELEMENTS openmp IF AXOM_ENABLE_OPENMP)

    blt_add_executable(NAME        ${arg_NAME}
                       SOURCES     ${arg_SOURCES}
                       HEADERS     ${arg_HEADERS}
                       INCLUDES    ${arg_INCLUDES}
                       DEFINES     ${arg_DEFINES}
                       DEPENDS_ON  ${arg_DEPENDS_ON} ${axom_device_depends}
                       OUTPUT_DIR  ${arg_OUTPUT_DIR}
                       OUTPUT_NAME ${arg_OUTPUT_NAME}
                       FOLDER      ${arg_FOLDER})

endmacro(axom_add_executable)


##------------------------------------------------------------------------------
## axom_add_library(
##                  NAME         <libname>
##                  SOURCES      [source1 [source2 ...]]
##                  HEADERS      [header1 [header2 ...]]
##                  INCLUDES     [dir1 [dir2 ...]]
##                  DEFINES      [define1 [define2 ...]]
##                  DEPENDS_ON   [dep1 ...] 
##                  OUTPUT_NAME  [name]
##                  OUTPUT_DIR   [dir]
##                  SHARED       [TRUE | FALSE]
##                  OBJECT       [TRUE | FALSE]
##                  CLEAR_PREFIX [TRUE | FALSE]
##                  FOLDER       [name])
##
## Wrapper around blt_add_library
##------------------------------------------------------------------------------
macro(axom_add_library)

    set(options)
    set(singleValueArgs NAME OUTPUT_NAME OUTPUT_DIR SHARED OBJECT CLEAR_PREFIX FOLDER)
    set(multiValueArgs SOURCES HEADERS INCLUDES DEFINES DEPENDS_ON)

    # parse the arguments
    cmake_parse_arguments(arg
        "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN} )

    if("openmp" IN_LIST arg_DEPENDS_ON)
        message(FATAL_ERROR "Do not add 'openmp' to Axom libraries to avoid propegation. It is handled automatically.")
    endif()

    if(DEFINED arg_OUTPUT_NAME)
        set(_output_name "${arg_OUTPUT_NAME}")
    elseif(NOT DEFINED arg_SOURCES)
        # Can't set OUTPUT_NAME on header only libraries
        set(_output_NAME)
    elseif(${arg_NAME} IN_LIST AXOM_COMPONENTS_FULL)
        # Prefix the output name with "axom_" for component libraries
        set(_output_name "axom_${arg_NAME}")
    else()
        set(_output_name "${arg_NAME}")
    endif()

    blt_add_library(NAME        ${arg_NAME}
                    SOURCES     ${arg_SOURCES}
                    HEADERS     ${arg_HEADERS}
                    INCLUDES    ${arg_INCLUDES}
                    DEFINES     ${arg_DEFINES}
                    DEPENDS_ON  ${arg_DEPENDS_ON} ${axom_device_depends}
                    OUTPUT_DIR  ${arg_OUTPUT_DIR}
                    OUTPUT_NAME ${_output_name}
                    FOLDER      ${arg_FOLDER})

    if(AXOM_ENABLE_OPENMP AND (NOT "${arg_SOURCES}" STREQUAL ""))
        # Do not propegate OpenMP due to generator expressions evaluating early
        target_link_libraries(${arg_NAME} PRIVATE openmp)
    endif()

endmacro(axom_add_library)


##------------------------------------------------------------------------------
## axom_add_test(NAME            [name]
##               COMMAND         [command] 
##               NUM_MPI_TASKS   [n]
##               NUM_OMP_THREADS [n]
##               CONFIGURATIONS  [config1 [config2...]])
##
## Wrapper around blt_add_test() that handles functionality 
## that Axom applies to all tests.
##
## Note that NUM_OMP_THREADS delegates to the corresponding argument 
## in blt_add_test() and sets the OpenMP environment variable OMP_NUM_THREADS.
## When AXOM_ENABLE_OPENMP is set and NUM_OMP_THREADS is not provided, 
## this macros also sets the environment variable OMP_NUM_THREADS=1.
##------------------------------------------------------------------------------
macro(axom_add_test)

    set(options )
    set(singleValueArgs NAME NUM_MPI_TASKS NUM_OMP_THREADS)
    set(multiValueArgs COMMAND CONFIGURATIONS)

    # Parse the arguments to the macro
    cmake_parse_arguments(arg
        "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN} )


    blt_add_test(NAME            ${arg_NAME}
                 COMMAND         ${arg_COMMAND}
                 NUM_MPI_TASKS   ${arg_NUM_MPI_TASKS}
                 NUM_OMP_THREADS ${arg_NUM_OMP_THREADS}
                 CONFIGURATIONS  ${arg_CONFIGURATIONS} )

    #--------------------------------------------------------------------------
    # Newer versions of OpenMPI require OMPI_MCA_rmaps_base_oversubscribe=1
    # to run with more tasks than actual cores
    # Since this is an OpenMPI specific env var, it shouldn't interfere
    # with other mpi implementations.
    #--------------------------------------------------------------------------
    set_property(TEST ${arg_NAME}
                 APPEND
                 PROPERTY ENVIRONMENT  "OMPI_MCA_rmaps_base_oversubscribe=1")

    #--------------------------------------------------------------------------
    # Cap OpenMP parallelism for tests that do not explicitly
    # specify NUM_OMP_THREADS to avoid accidental oversubscription
    #--------------------------------------------------------------------------
    if(AXOM_ENABLE_OPENMP AND (NOT arg_NUM_OMP_THREADS))
        set_property(TEST ${arg_NAME}
                     APPEND 
                     PROPERTY ENVIRONMENT OMP_NUM_THREADS=1)
    endif()

endmacro(axom_add_test)


#------------------------------------------------------------------------------
# axom_assert_is_directory(DIR_VARIABLE <variable that holds the prefix>)
#
# Asserts that the given DIR_VARIABLE's value is a directory and exists.
# Fails with a helpful message when it doesn't.
#------------------------------------------------------------------------------
macro(axom_assert_is_directory)

    set(options)
    set(singleValueArgs DIR_VARIABLE)
    set(multiValueArgs)

    # Parse the arguments to the macro
    cmake_parse_arguments(arg
         "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

    if (NOT EXISTS "${${arg_DIR_VARIABLE}}")
        message(FATAL_ERROR "Given ${arg_DIR_VARIABLE} does not exist: ${${arg_DIR_VARIABLE}}")
    endif()

    if (NOT IS_DIRECTORY "${${arg_DIR_VARIABLE}}")
        message(FATAL_ERROR "Given ${arg_DIR_VARIABLE} is not a directory: ${${arg_DIR_VARIABLE}}")
    endif()

endmacro(axom_assert_is_directory)

#------------------------------------------------------------------------------
# axom_assert_find_succeeded(PROJECT_NAME <project name>
#                             TARGET       <found target>
#                             DIR_VARIABLE <variable that holds the prefix>)
#
# Asserts that the given PROJECT_NAME's TARGET exists.
# Fails with a helpful message when it doesn't.
#------------------------------------------------------------------------------
macro(axom_assert_find_succeeded)

    set(options)
    set(singleValueArgs DIR_VARIABLE PROJECT_NAME TARGET)
    set(multiValueArgs)

    # Parse the arguments to the macro
    cmake_parse_arguments(arg
         "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

    message(STATUS "Checking for expected ${arg_PROJECT_NAME} target '${arg_TARGET}'")
    if (NOT TARGET ${arg_TARGET})
        message(FATAL_ERROR "${arg_PROJECT_NAME} failed to load: ${${arg_DIR_VARIABLE}}")
    else()
        message(STATUS "${arg_PROJECT_NAME} loaded: ${${arg_DIR_VARIABLE}}")
    endif()

endmacro(axom_assert_find_succeeded)

##------------------------------------------------------------------------------
## convert_to_native_escaped_file_path( path output )
##
## This macro converts a cmake path to a platform specific string literal
## usable in C++.  (For example, on windows C:/Path will become C:\\Path)
##------------------------------------------------------------------------------

macro(convert_to_native_escaped_file_path path output)
    file(TO_NATIVE_PATH ${path} ${output})
    string(REPLACE "\\" "\\\\"  ${output} "${${output}}")
endmacro()


##------------------------------------------------------------------------------
## axom_check_code_compiles
##
## This macro checks if a snippet of C++ code compiles.
##
## SOURCE_STRING The source snippet to compile.
## Must be a valid C++ program with a main() function.
## Note: This parameter should be passed in as a quoted string variable. Otherwise,
## cmake will convert the string into a list and lose the semicolons.
## E.g. axom_check_code_compiles(SOURCE_STRING "${str_var}" ...)
##
## CODE_COMPILES A boolean variable the contains the compilation result.
##
## VERBOSE_OUTPUT Optional parameter to output debug information (Default: off)
##------------------------------------------------------------------------------
macro(axom_check_code_compiles)

    set(options)
    set(singleValueArgs CODE_COMPILES VERBOSE_OUTPUT)
    set(multiValueArgs SOURCE_STRING)

    # Parse the arguments to the macro
    cmake_parse_arguments(arg
         "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

    # Check the arguments
    if(NOT DEFINED arg_SOURCE_STRING)
        message(FATAL_ERROR "[axom_check_code_compiles] 'SOURCE_STRING' is a required parameter")
    endif()
    if(NOT DEFINED arg_CODE_COMPILES)
        message(FATAL_ERROR "[axom_check_code_compiles] 'CODE_COMPILES' is a required parameter")
    endif()

    if(NOT DEFINED arg_VERBOSE_OUTPUT)
        set(arg_VERBOSE_OUTPUT FALSE)
    endif()

    if(${arg_VERBOSE_OUTPUT})
        message(STATUS "[axom_check_code_compiles] Attempting to compile source string: \n${arg_SOURCE_STRING}")
    endif()

    # Write string as temp file, try to compile it and then remove file
    string(RANDOM LENGTH 5 _rand)
    set(_fname ${CMAKE_CURRENT_BINARY_DIR}/_axomCheckCompiles${_rand}.cpp)
    file(WRITE ${_fname} "${arg_SOURCE_STRING}")
    try_compile(${arg_CODE_COMPILES}
                ${CMAKE_CURRENT_BINARY_DIR}/CMakeTmp
                SOURCES ${_fname}
                CXX_STANDARD ${CMAKE_CXX_STANDARD}
                OUTPUT_VARIABLE _res)
    file(REMOVE ${_fname})

    if(${arg_VERBOSE_OUTPUT})
        message(STATUS "[axom_check_code_compiles] Compiler output: \n${_res}\n")

        if(${arg_CODE_COMPILES})
            message(STATUS "[axom_check_code_compiles] The code snippet successfully compiled")
        else()
            message(STATUS "[axom_check_code_compiles] The code snippet failed to compile")
        endif()
    endif()

    # clear the variables set within the macro
    unset(_fname)
    unset(_res)

endmacro(axom_check_code_compiles)


##------------------------------------------------------------------------------
## axom_component_requires
##
## This macro checks for the required dependencies of the given component
##
## NAME - The name of the component that we are checking the dependencies of
##
## COMPONENTS - Internal required components
##
## TPLS - Third party required libraries
##------------------------------------------------------------------------------
macro(axom_component_requires)

    set(options)
    set(singleValueArgs NAME)
    set(multiValueArgs COMPONENTS TPLS)

    # Parse the arguments to the macro
    cmake_parse_arguments(arg
         "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

    foreach(_dep ${arg_COMPONENTS})
        string(TOUPPER ${_dep} _ucdep)
        if(NOT AXOM_ENABLE_${_ucdep})
            message(FATAL_ERROR "${arg_NAME} requires ${_dep}. Set AXOM_ENABLE_${_ucdep} to ON.")
        endif()
    endforeach()

    foreach(_dep ${arg_TPLS})
        string(TOUPPER ${_dep} _ucdep)
        if(NOT ${_ucdep}_FOUND)
            message(FATAL_ERROR "${arg_NAME} requires ${_dep}. Set ${_ucdep}_DIR to location of a ${_dep} install.")
        endif()
    endforeach()

endmacro(axom_component_requires)

##------------------------------------------------------------------------------
## axom_install_component
##
## This macro installs libraries, fortran modules, headers, and exports the CMake
## target while preserving the directory stucture.  This macro assumes the following:
##
##    * CMake Target to install is the same as NAME
##    * If there is a Fortran module built, it is named axom_<NAME>.mod
##
## NAME - The name of the component that we are installing.
##
## HEADERS - Headers to be installed
##
##------------------------------------------------------------------------------
macro(axom_install_component)

    set(options )
    set(singleValueArgs NAME)
    set(multiValueArgs HEADERS)

    # Parse the arguments to the macro
    cmake_parse_arguments(arg
         "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

    set(_header_base_dir include/axom/${arg_NAME})

    foreach( _file ${arg_HEADERS} )
        get_filename_component( _dir ${_file} DIRECTORY )
        install(FILES ${_file} DESTINATION ${_header_base_dir}/${_dir} )
    endforeach()

    if(ENABLE_FORTRAN)
        set(_mod ${CMAKE_Fortran_MODULE_DIRECTORY}/axom_${arg_NAME}.mod)
        # TODO: Remove optional once all components have fortran wrappers
        install(FILES ${_mod} DESTINATION lib/fortran OPTIONAL)
    endif()

endmacro(axom_install_component)


##------------------------------------------------------------------------------
## axom_write_unified_header
##
## This macro writes the unified header (axom/<lowered NAME>.hpp) to the build directory for the
## given component NAME with the given HEADERS included inside of it.
##
## NAME - The name of the component for the unified header.
##
## HEADERS - Headers to be included in the header.
##
##------------------------------------------------------------------------------
macro(axom_write_unified_header)

    set(options )
    set(singleValueArgs NAME)
    set(multiValueArgs HEADERS EXCLUDE)

    # Parse the arguments to the macro
    cmake_parse_arguments(arg
         "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

    string(TOUPPER ${arg_NAME} _ucname)
    string(TOLOWER ${arg_NAME} _lcname)
    set(_header ${PROJECT_BINARY_DIR}/include/axom/${_lcname}.hpp)
    set(_tmp_header ${_header}.tmp)

    file(WRITE ${_tmp_header} "\/\/ Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
\/\/ other Axom Project Developers. See the top-level LICENSE file for details.
\/\/
\/\/ SPDX-License-Identifier: (BSD-3-Clause)
\n
")

    file(APPEND ${_tmp_header} "#ifndef AXOM_UNIFIED_${_ucname}_HPP\n")
    file(APPEND ${_tmp_header} "#define AXOM_UNIFIED_${_ucname}_HPP\n\n")

    file(APPEND ${_tmp_header} "#include \"axom\/config.hpp\"\n\n")

    foreach(_file ${arg_HEADERS})
        set(_headerPath "axom\/${_lcname}\/${_file}")

        if(${_file} IN_LIST arg_EXCLUDE)
            continue()
        elseif(${_headerPath} MATCHES "(\/detail\/)|(\/internal\/)")
            continue()
        else()
            file(APPEND ${_tmp_header} "#include \"${_headerPath}\"\n")
        endif()
    endforeach()

    file(APPEND ${_tmp_header} "\n#endif // AXOM_UNIFIED_${_ucname}_HPP\n")

    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_tmp_header} ${_header})
    execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${_tmp_header})

    install(FILES       ${_header}
            DESTINATION include/axom)

endmacro(axom_write_unified_header)

##------------------------------------------------------------------------------
## axom_configure_file
##
## This macro is a thin wrapper over the builtin configure_file command.
## It has the same arguments/options as configure_file but introduces an
## intermediate file that is only copied to the target file if the target differs
## from the intermediate.
##------------------------------------------------------------------------------
macro(axom_configure_file _source _target)
    set(_tmp_target ${_target}.tmp)
    configure_file(${_source} ${_tmp_target} ${ARGN})
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_tmp_target} ${_target})
    execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${_tmp_target})
endmacro(axom_configure_file)
