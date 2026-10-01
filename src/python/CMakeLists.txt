# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Compile Axom's Python bindings against an installed Axom and its third-party libraries.
#------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.21)
project(axom_python LANGUAGES C CXX)

# Enable with both build-time flags; see src/python/README.md:
#   -C cmake.define.AXOM_PYTHON_STABLE_ABI=ON  -C wheel.py-api=cp312
option(AXOM_PYTHON_STABLE_ABI
       "Build the extension against Python's stable ABI (abi3); needs Python >= 3.12" OFF)

if(AXOM_PYTHON_STABLE_ABI AND CMAKE_VERSION VERSION_LESS 3.26)
    message(FATAL_ERROR
            "AXOM_PYTHON_STABLE_ABI=ON requires CMake >= 3.26 because "
            "FindPython's Development.SABIModule component is needed "
            "to locate Python stable ABI development artifacts. "
            "Current CMake is ${CMAKE_VERSION}. "
            "Build with CMake >= 3.26 or leave AXOM_PYTHON_STABLE_ABI off.")
endif()

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.26)
    find_package(Python 3.10 REQUIRED
                 COMPONENTS Interpreter Development.Module
                 OPTIONAL_COMPONENTS Development.SABIModule)
else()
    find_package(Python 3.10 REQUIRED
                 COMPONENTS Interpreter Development.Module)
endif()

if(AXOM_PYTHON_STABLE_ABI AND NOT Python_Development.SABIModule_FOUND)
    message(FATAL_ERROR
            "AXOM_PYTHON_STABLE_ABI=ON requires CMake >= 3.26 and Python's "
            "Development.SABIModule component (CPython >= 3.12), "
            "which was not found for ${Python_EXECUTABLE} (version ${Python_VERSION}). "
            "Build with a 3.12+ interpreter that provides stable ABI development artifacts, "
            "or leave AXOM_PYTHON_STABLE_ABI off.")
endif()

# pip/uv builds may run outside the shell's working directory.
# Require absolute paths so find_package resolves the intended Axom install.
foreach(_axom_py_dir_var AXOM_DIR axom_DIR)
    if(DEFINED ${_axom_py_dir_var} AND NOT IS_ABSOLUTE "${${_axom_py_dir_var}}")
        message(FATAL_ERROR
                "${_axom_py_dir_var}='${${_axom_py_dir_var}}' must be an absolute path. "
                "The build ran from the current CMake directory '${CMAKE_CURRENT_SOURCE_DIR}', "
                "rather than the shell's working directory.")
    endif()
endforeach()

# AXOM_DIR accepts the install prefix or the directory containing axom-config.cmake.
# CMake's axom_DIR takes precedence when set.
if(DEFINED AXOM_DIR AND NOT DEFINED axom_DIR)
    if(EXISTS "${AXOM_DIR}/lib/cmake/axom-config.cmake")
        set(_axom_py_package_dir "${AXOM_DIR}/lib/cmake")
    elseif(EXISTS "${AXOM_DIR}/axom-config.cmake")
        set(_axom_py_package_dir "${AXOM_DIR}")
    else()
        message(FATAL_ERROR
                "AXOM_DIR='${AXOM_DIR}' is not an Axom install: "
                "neither '${AXOM_DIR}/lib/cmake/axom-config.cmake' "
                "nor '${AXOM_DIR}/axom-config.cmake' exists. "
                "Set AXOM_DIR to the Axom install prefix.")
    endif()
    set(axom_DIR "${_axom_py_package_dir}" CACHE PATH "Axom CMake package directory")
    unset(_axom_py_package_dir)
endif()

find_package(axom CONFIG REQUIRED)
find_package(nanobind CONFIG REQUIRED)   # supplied via build-system.requires

# Compile the bindings with the C++ standard of the Axom install.
if(NOT AXOM_CXX_STANDARD)
    set(AXOM_CXX_STANDARD 20)
endif()

if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD ${AXOM_CXX_STANDARD})
elseif(CMAKE_CXX_STANDARD LESS AXOM_CXX_STANDARD)
    message(FATAL_ERROR
            "CMAKE_CXX_STANDARD is C++${CMAKE_CXX_STANDARD}, "
            "but the Axom install at '${axom_DIR}' was built with C++${AXOM_CXX_STANDARD}.")
endif()

# Match how BLT compiles Axom
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT DEFINED CMAKE_CXX_EXTENSIONS)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

set(_axom_py_conduit_cmake_dir "${Conduit_DIR}")
if(NOT _axom_py_conduit_cmake_dir AND AXOM_CONDUIT_DIR)
    set(_axom_py_conduit_cmake_dir "${AXOM_CONDUIT_DIR}/lib/cmake/conduit")
endif()

# Use a separate override because ConduitConfig.cmake sets CONDUIT_PYTHON_MODULE_DIR
# as a normal variable, shadowing caller-supplied cache values.
set(AXOM_PYTHON_CONDUIT_MODULE_DIR "" CACHE PATH
    "Directory holding the same-build conduit Python package; overrides the value recorded by the Axom install")

if(AXOM_PYTHON_CONDUIT_MODULE_DIR)
    set(_axom_py_conduit_python_module_dir "${AXOM_PYTHON_CONDUIT_MODULE_DIR}")
elseif(AXOM_CONDUIT_PYTHON_MODULE_DIR)
    set(_axom_py_conduit_python_module_dir "${AXOM_CONDUIT_PYTHON_MODULE_DIR}")
else()
    set(_axom_py_conduit_python_module_dir "${CONDUIT_PYTHON_MODULE_DIR}")
endif()

if(_axom_py_conduit_python_module_dir
   AND NOT IS_ABSOLUTE "${_axom_py_conduit_python_module_dir}")
    get_filename_component(_axom_py_conduit_python_module_dir
                           "${AXOM_CONDUIT_DIR}/${_axom_py_conduit_python_module_dir}" ABSOLUTE)
endif()

if(NOT _axom_py_conduit_python_module_dir
   OR NOT EXISTS "${_axom_py_conduit_python_module_dir}/conduit")
    message(FATAL_ERROR
            "Could not find Conduit's Python package for this Axom install.\n"
            "  Looked for a directory containing a 'conduit' package at: "
            "'${_axom_py_conduit_python_module_dir}'\n"
            "  AXOM_CONDUIT_DIR                = ${AXOM_CONDUIT_DIR}\n"
            "  AXOM_CONDUIT_PYTHON_MODULE_DIR  = ${AXOM_CONDUIT_PYTHON_MODULE_DIR}\n"
            "Either rebuild Axom against a Conduit configured with Python support, or "
            "pass the path explicitly:\n"
            "  -C cmake.define.AXOM_PYTHON_CONDUIT_MODULE_DIR=<path>")
endif()

# Keep the wheel metadata version aligned with the Axom install it links.
if(DEFINED SKBUILD_PROJECT_VERSION AND DEFINED AXOM_VERSION_MAJOR)
    set(_axom_installed_version
        "${AXOM_VERSION_MAJOR}.${AXOM_VERSION_MINOR}.${AXOM_VERSION_PATCH}")
    if(NOT SKBUILD_PROJECT_VERSION VERSION_EQUAL _axom_installed_version)
        message(FATAL_ERROR
                "Axom version mismatch: the wheel's metadata version is ${SKBUILD_PROJECT_VERSION} "
                "but the Axom install it would link against is ${_axom_installed_version} (axom_DIR=${axom_DIR}). "
                "Build the wheel from the source tree that produced the install, "
                "or point AXOM_DIR at an install built from this source tree.")
    endif()
endif()

#------------------------------------------------------------------------------
# Compile the binding source from the Sidre component.
#------------------------------------------------------------------------------
set(_sidre_binding_sources
    "${CMAKE_CURRENT_SOURCE_DIR}/../axom/sidre/nanobind_sidre.cpp")

# Match the in-tree Sidre binding domain so future Axom modules can share bound
# C++ types when built together.
set(_axom_nb_module_args NB_DOMAIN axom)
if(AXOM_PYTHON_STABLE_ABI)
    list(APPEND _axom_nb_module_args STABLE_ABI)
endif()
nanobind_add_module(_sidre ${_axom_nb_module_args} ${_sidre_binding_sources})

target_link_libraries(_sidre PRIVATE axom::sidre conduit::conduit_python)

set_target_properties(_sidre PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

#------------------------------------------------------------------------------
# Generate configuration helpers for this Axom and Conduit install.
#------------------------------------------------------------------------------
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/axom-python-host-config.cmake.in"
               "${CMAKE_CURRENT_BINARY_DIR}/axom-python-host-config.cmake"
               @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/axom-python-env.sh.in"
               "${CMAKE_CURRENT_BINARY_DIR}/axom-python-env.sh"
               @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/axom-python-host-config.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/axom-python-env.sh"
        DESTINATION axom/share)

# Prefix the filename to avoid replacing another package's conduit.pth.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/axom-conduit.pth.in"
               "${CMAKE_CURRENT_BINARY_DIR}/axom-conduit.pth"
               @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/axom-conduit.pth"
        DESTINATION ".")

#------------------------------------------------------------------------------
# TODO: Test wheel builds with HIP and CUDA.
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# nanobind_add_stub imports _sidre to generate its type stub.
# Disable it if the build host cannot import the module:
#   -C cmake.define.AXOM_PYTHON_GENERATE_STUB=OFF
#------------------------------------------------------------------------------
option(AXOM_PYTHON_GENERATE_STUB
       "Generate the _sidre.pyi type stub at build time (imports the module)" ON)

if(AXOM_PYTHON_GENERATE_STUB)
    nanobind_add_stub(
        _sidre_stub
        MODULE _sidre
        OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/_sidre.pyi"
        PYTHON_PATH $<TARGET_FILE_DIR:_sidre>
        DEPENDS _sidre)

    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/_sidre.pyi"
            DESTINATION axom/sidre)
endif()

#------------------------------------------------------------------------------
# Pure-Python package files ship via [tool.scikit-build] wheel.packages.
#------------------------------------------------------------------------------
install(TARGETS _sidre LIBRARY DESTINATION axom/sidre)
