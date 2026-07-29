# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Thin, binding-only build of Axom's Python extension module(s).
#
# This project compiles Axom's nanobind translation unit(s) against an already-installed Axom and Conduit.
# It does NOT build Axom or its third-party libraries. 
#
# It is driven by the scikit-build-core project in the sibling pyproject.toml (see src/python/README.md),
# but is a valid standalone CMake project as well.
#
# When driven by scikit-build-core (uv build / uv pip install),
# point find_package at the installs with -DAXOM_DIR / -DConduit_DIR:
# scikit-build-core force-sets CMAKE_PREFIX_PATH to the isolated build env (to locate its own nanobind),
# so a user-supplied CMAKE_PREFIX_PATH is ignored.
# A standalone cmake invocation has no such layer and uses CMAKE_PREFIX_PATH directly:
#
#   cmake -S src/python -B build/py -DCMAKE_PREFIX_PATH="$AXOM_INSTALL;$CONDUIT_INSTALL;$(python -m nanobind --cmake_dir)"
#   cmake --build   build/py
#   cmake --install build/py --prefix <site-packages-shaped-dir>
#------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.21)
project(axom_python LANGUAGES C CXX)

# Optional: build against Python's stable ABI (abi3) so one wheel serves every CPython >= 3.12 on the machine.
# Opt-in (default OFF) to provide per-Python wheel tags by default.
# Enable with BOTH build-time flags (they must agree; see pyproject.toml and src/python/README.md):
#   -C cmake.define.AXOM_PYTHON_STABLE_ABI=ON  -C wheel.py-api=cp312
# nanobind silently builds a non-stable module below Python 3.12,
# so an abi3 wheel needs a 3.12+ interpreter (which also provides Development.SABIModule).
option(AXOM_PYTHON_STABLE_ABI
       "Build the extension against Python's stable ABI (abi3); needs Python >= 3.12" OFF)

# Development.Module (not Development) so we link libpython-free extensions,
# matching the in-tree build's discovery.
# Development.SABIModule is requested as an OPTIONAL component (nanobind's recommended pattern) 
# so a stable-ABI build (AXOM_PYTHON_STABLE_ABI=ON, which needs Python >= 3.12) can find it,
# without a conditional component list.
find_package(Python 3.10 REQUIRED
             COMPONENTS Interpreter Development.Module
             OPTIONAL_COMPONENTS Development.SABIModule)

# An optional component that is silently absent would turn a requested abi3 build
# into a confusing failure inside nanobind (or a non-stable module), so fail here
# with the actionable message instead.
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

# The wheel's version metadata is read from this checkout's src/cmake/AxomVersion.cmake 
# (see [[tool.dynamic-metadata]] in pyproject.toml),
# but the extension is compiled and linked against the Axom install found above.
#
# Building from a checkout that does not match the install would produce a wheel
# whose axom.__version__ misreports its own binary, so they must agree.
# (AXOM_VERSION_* come from axom-config.cmake; SKBUILD_PROJECT_VERSION is injected by scikit-build-core,
# so this check is active for wheel builds and skipped for a standalone cmake invocation.)
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

# Build the extension under the shared 'axom' nanobind domain so that C++ types bound in one Axom module
# (e.g. a sidre::Group*) are recognized by another module from the same build.
# nanobind only shares type bindings across modules that agree on domain *and* nanobind ABI, compiler, and build mode,
# hence the one-build/one-wheel rule. This mirrors src/axom/sidre/CMakeLists.txt.
#
# STABLE_ABI is appended only when AXOM_PYTHON_STABLE_ABI is set
# nanobind degrades to a non-stable build below 3.12.
set(_axom_nb_module_args NB_DOMAIN axom)
if(AXOM_PYTHON_STABLE_ABI)
    list(APPEND _axom_nb_module_args STABLE_ABI)
endif()
nanobind_add_module(_sidre ${_axom_nb_module_args} ${_sidre_binding_sources})

# conduit::conduit_python provides conduit_python.hpp and is needed only by the binding TU, not by libsidre.
target_link_libraries(_sidre PRIVATE axom::sidre conduit::conduit_python)

# Bake the linked libraries' locations into the module's rpath so it resolves libsidre/libconduit/HDF5
# from the Axom/Conduit install with no LD_LIBRARY_PATH.
# pyproject also sets CMAKE_INSTALL_RPATH_USE_LINK_PATH for the wheel;
# this makes the standalone `cmake --install` path behave the same.
set_target_properties(_sidre PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

#------------------------------------------------------------------------------
# Generated runtime/development helpers.
#
# A thin wheel is tied to the Axom/Conduit install it was built against. Ship the
# corresponding Conduit Python path and CMake configuration hints with the wheel
# so users do not have to rediscover them by hand.
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

if(_axom_py_conduit_python_module_dir)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/conduit.pth.in"
                   "${CMAKE_CURRENT_BINARY_DIR}/conduit.pth"
                   @ONLY)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/conduit.pth"
            DESTINATION ".")
endif()

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
# Type stubs (PEP 561).
#
# nanobind_add_stub imports the module ('import _sidre') to introspect it.
#
# Turn this off on hosts where the build-time import is problematic:
#   -C cmake.define.AXOM_PYTHON_GENERATE_STUB=OFF
# The wheel still ships the checked-in package stub (axom/sidre/__init__.pyi)
# and py.typed via wheel.packages, so it stays a typed package either way;
# only the detailed _sidre.pyi (which __init__.pyi re-exports) is then absent.
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
# Install only build products here. 
# The pure-Python package files 
# (axom/__init__.py, axom/py.typed, axom/sidre/__init__.py, axom/sidre/__init__.pyi)
# ship as pure Python via [tool.scikit-build] wheel.packages in pyproject.toml,
# so they are not installed by CMake (avoiding double-packaging).
#------------------------------------------------------------------------------
install(TARGETS _sidre LIBRARY DESTINATION axom/sidre)
