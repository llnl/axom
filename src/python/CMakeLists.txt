# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Thin, binding-only build of Axom's Python extension module(s).
#
# This project compiles Axom's Python bindings against an installed Axom.
# It does not build Axom or its third-party libraries.
#------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.21)
project(axom_python LANGUAGES C CXX)

# Enable with both build-time flags; see src/python/README.md:
#   -C cmake.define.AXOM_PYTHON_STABLE_ABI=ON  -C wheel.py-api=cp312
option(AXOM_PYTHON_STABLE_ABI
       "Build the extension against Python's stable ABI (abi3); needs Python >= 3.12" OFF)

find_package(Python 3.10 REQUIRED
             COMPONENTS Interpreter Development.Module
             OPTIONAL_COMPONENTS Development.SABIModule)

if(AXOM_PYTHON_STABLE_ABI AND NOT Python_Development.SABIModule_FOUND)
    message(FATAL_ERROR
            "AXOM_PYTHON_STABLE_ABI=ON requires Python's Development.SABIModule "
            "component (CPython >= 3.12), which was not found for "
            "${Python_EXECUTABLE} (version ${Python_VERSION}). Build with a 3.12+ "
            "interpreter or leave AXOM_PYTHON_STABLE_ABI off.")
endif()

if(DEFINED AXOM_DIR AND NOT DEFINED axom_DIR)
    set(axom_DIR "${AXOM_DIR}" CACHE PATH "Axom CMake package directory")
endif()

find_package(axom CONFIG REQUIRED)
find_package(nanobind CONFIG REQUIRED)   # supplied via build-system.requires

set(_axom_py_conduit_cmake_dir "")
if(DEFINED Conduit_DIR)
    set(_axom_py_conduit_cmake_dir "${Conduit_DIR}")
elseif(DEFINED CONDUIT_INSTALL_PREFIX)
    set(_axom_py_conduit_cmake_dir "${CONDUIT_INSTALL_PREFIX}/lib/cmake/conduit")
elseif(DEFINED AXOM_CONDUIT_DIR)
    set(_axom_py_conduit_cmake_dir "${AXOM_CONDUIT_DIR}/lib/cmake/conduit")
endif()

set(_axom_py_conduit_prefix "")
if(DEFINED CONDUIT_INSTALL_PREFIX)
    set(_axom_py_conduit_prefix "${CONDUIT_INSTALL_PREFIX}")
elseif(DEFINED AXOM_CONDUIT_DIR)
    set(_axom_py_conduit_prefix "${AXOM_CONDUIT_DIR}")
elseif(DEFINED CONDUIT_DIR)
    set(_axom_py_conduit_prefix "${CONDUIT_DIR}")
endif()

set(_axom_py_conduit_python_module_dir "")
if(DEFINED CONDUIT_PYTHON_MODULE_DIR)
    if(IS_ABSOLUTE "${CONDUIT_PYTHON_MODULE_DIR}")
        set(_axom_py_conduit_python_module_dir "${CONDUIT_PYTHON_MODULE_DIR}")
    elseif(_axom_py_conduit_prefix)
        set(_axom_py_conduit_python_module_dir
            "${_axom_py_conduit_prefix}/${CONDUIT_PYTHON_MODULE_DIR}")
    endif()
elseif(_axom_py_conduit_prefix)
    foreach(_axom_py_conduit_python_module_dir_candidate
            "${_axom_py_conduit_prefix}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages"
            "${_axom_py_conduit_prefix}/python-modules")
        if(EXISTS "${_axom_py_conduit_python_module_dir_candidate}/conduit")
            set(_axom_py_conduit_python_module_dir
                "${_axom_py_conduit_python_module_dir_candidate}")
            break()
        endif()
    endforeach()
endif()

if(NOT _axom_py_conduit_python_module_dir
   OR NOT EXISTS "${_axom_py_conduit_python_module_dir}/conduit")
    message(FATAL_ERROR
            "Could not determine the Conduit Python module directory for this "
            "Axom install. Pass -C cmake.define.CONDUIT_PYTHON_MODULE_DIR=<path> "
            "where <path> contains the same-build conduit Python package.")
endif()

# Keep the wheel metadata version aligned with the Axom install it links.
if(DEFINED SKBUILD_PROJECT_VERSION AND DEFINED AXOM_VERSION_MAJOR)
    set(_axom_installed_version
        "${AXOM_VERSION_MAJOR}.${AXOM_VERSION_MINOR}.${AXOM_VERSION_PATCH}")
    if(NOT SKBUILD_PROJECT_VERSION VERSION_EQUAL _axom_installed_version)
        message(FATAL_ERROR
                "Axom version mismatch: the wheel's metadata version is "
                "${SKBUILD_PROJECT_VERSION} (from src/cmake/AxomVersion.cmake in this "
                "source tree) but the Axom install it would link against is "
                "${_axom_installed_version} (AXOM_DIR=${axom_DIR}). Build the wheel from "
                "the source tree that produced the install, or point AXOM_DIR at an "
                "install built from this source tree.")
    endif()
endif()

#------------------------------------------------------------------------------
# Binding sources. The translation unit lives with its component, not in this project directory.
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
# Generated runtime/development helpers tied to the Axom/Conduit install.
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

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/conduit.pth.in"
               "${CMAKE_CURRENT_BINARY_DIR}/conduit.pth"
               @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/conduit.pth"
        DESTINATION ".")

#------------------------------------------------------------------------------
# HIP: replicate the in-tree special-case from src/axom/sidre/CMakeLists.txt.
# CMake treats MODULE libraries separately from executables, so HIP executable
# link flags are not applied to the module automatically.
#------------------------------------------------------------------------------
if(AXOM_ENABLE_HIP)
    set_source_files_properties(${_sidre_binding_sources} PROPERTIES LANGUAGE HIP)
    string(REPLACE " " ";" _axom_py_module_link_flags "${CMAKE_EXE_LINKER_FLAGS}")
    target_link_options(_sidre PRIVATE ${_axom_py_module_link_flags})
    unset(_axom_py_module_link_flags)
endif()

#------------------------------------------------------------------------------
# nanobind_add_stub imports the module ('import _sidre') to introspect it.
# Turn this off on hosts where the build-time import is problematic:
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
