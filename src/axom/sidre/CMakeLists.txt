# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# SIDRE -- SI(mulation) D(ata) RE(pository)
#
# Also includes functionality for SPIO -- S(idre) P(arallel) I/O
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Check necessary dependencies
#------------------------------------------------------------------------------
axom_component_requires(NAME       Sidre
                        COMPONENTS SLIC
                        TPLS       Conduit )

#------------------------------------------------------------------------------
# Add the sidre sources
#------------------------------------------------------------------------------
set(sidre_headers
    core/Array.hpp
    core/Buffer.hpp
    core/ConduitMemory.hpp
    core/Group.hpp
    core/DataStore.hpp
    core/View.hpp
    core/Attribute.hpp
    core/AttrValues.hpp
    core/SidreTypes.hpp
    core/SidreDataTypeIds.h )

set(sidre_sources
    core/Buffer.cpp
    core/ConduitMemory.cpp
    core/Group.cpp
    core/DataStore.cpp
    core/View.cpp
    core/Attribute.cpp
    core/AttrValues.cpp )

# Add spio headers and sources when MPI is available
if(AXOM_ENABLE_MPI)
   list(APPEND sidre_headers
        spio/IOManager.hpp
        spio/IOBaton.hpp )
   list(APPEND sidre_sources
        spio/IOManager.cpp
        spio/IOBaton.cpp )
endif()

if(SHROUD_FOUND)
    add_subdirectory(interface)
    add_subdirectory(spio/interface)
endif()

if (ENABLE_FORTRAN)
    # Sidre's C API is only built to provide a Fortran interface.
    list(APPEND sidre_headers
        interface/sidre.h
        interface/SidreTypes.h
        # generated headers
        interface/c_fortran/wrapSidre.h
        interface/c_fortran/wrapDataStore.h
        interface/c_fortran/wrapGroup.h
        interface/c_fortran/wrapBuffer.h
        interface/c_fortran/wrapView.h
        interface/c_fortran/typesSidre.h
    )
    list(APPEND sidre_sources
        # generated source
        interface/c_fortran/wrapSidre.cpp
        interface/c_fortran/wrapDataStore.cpp
        interface/c_fortran/wrapGroup.cpp
        interface/c_fortran/wrapBuffer.cpp
        interface/c_fortran/wrapView.cpp
        interface/c_fortran/utilSidre.cpp
        interface/c_fortran/wrapfsidre.F
    )

    # Add spio's C API when MPI is available
    if(AXOM_ENABLE_MPI)
        list(APPEND sidre_headers
             spio/interface/c_fortran/wrapIOManager.h
             spio/interface/c_fortran/typesSPIO.h
        )
        list(APPEND sidre_sources
             spio/interface/c_fortran/wrapfspio.f
             spio/interface/c_fortran/wrapIOManager.cpp
             spio/interface/c_fortran/utilSPIO.cpp
        )
    endif()
endif()

if (MFEM_FOUND)
    list(APPEND sidre_headers core/MFEMSidreDataCollection.hpp)
    list(APPEND sidre_sources core/MFEMSidreDataCollection.cpp)
endif()

#------------------------------------------------------------------------------
# Build and install the library
#------------------------------------------------------------------------------
set(sidre_depends
    slic
    conduit::conduit
    )

blt_list_append(TO sidre_depends ELEMENTS hdf5 IF HDF5_FOUND)
blt_list_append(TO sidre_depends ELEMENTS sparsehash IF SPARSEHASH_FOUND)
if(MFEM_FOUND)
    blt_list_append(TO sidre_depends ELEMENTS mfem)
    blt_list_append(TO sidre_depends ELEMENTS fmt)
endif()

# Include additional dependencies for spio when MPI is available
if(AXOM_ENABLE_MPI)
    list(APPEND sidre_depends conduit::conduit_mpi fmt)
    blt_list_append(TO sidre_depends ELEMENTS scr IF SCR_FOUND)
endif()

axom_add_library(NAME       sidre
                 SOURCES    ${sidre_sources}
                 HEADERS    ${sidre_headers}
                 DEPENDS_ON ${sidre_depends}
                 FOLDER     axom/sidre)

axom_write_unified_header(NAME    sidre
                          HEADERS ${sidre_headers})

axom_install_component(NAME    sidre
                       HEADERS ${sidre_headers}
                       )

if(ENABLE_FORTRAN AND AXOM_ENABLE_MPI)
    install(FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/axom_spio.mod DESTINATION lib/fortran)
endif()


if(NANOBIND_FOUND)
    # Python bindings for Sidre.
    # The pure-Python package scaffolding lives once under src/python/src/

    # site-packages-shaped install directory for Axom's Python package(s).
    # Keep this relative to the install prefix so `cmake --install --prefix`
    # relocates the Python package along with Axom's other install artifacts.
    set(_axom_python_install_default
        "lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages")
    set(_axom_python_install_description
        "Install destination for Axom's Python package(s), relative to the install prefix")
    set(AXOM_PYTHON_MODULE_INSTALL_PREFIX
        "${_axom_python_install_default}"
        CACHE PATH
        "${_axom_python_install_description}")

    if(IS_ABSOLUTE "${AXOM_PYTHON_MODULE_INSTALL_PREFIX}")
        file(RELATIVE_PATH _axom_python_install_relpath
             "${CMAKE_INSTALL_PREFIX}"
             "${AXOM_PYTHON_MODULE_INSTALL_PREFIX}")
        if(NOT _axom_python_install_relpath MATCHES "^\\.\\.")
            set(AXOM_PYTHON_MODULE_INSTALL_PREFIX
                "${_axom_python_install_relpath}"
                CACHE PATH
                "${_axom_python_install_description}"
                FORCE)
        endif()
        unset(_axom_python_install_relpath)
    endif()
    unset(_axom_python_install_default)
    unset(_axom_python_install_description)

    # Root of the staged package tree in the build directory.
    # The build-tree interpreter (and run_python_with_axom.sh, via _PYEXT_DIR)
    # puts this single directory on PYTHONPATH and then 'import axom.sidre' works.
    set(_axom_py_build_root "${PROJECT_BINARY_DIR}/python")
    set(_pysidre_pkg_src    "${CMAKE_CURRENT_SOURCE_DIR}/../../python/src")

    nanobind_add_module(_sidre nanobind_sidre.cpp)
    # conduit::conduit_python provides conduit_python.hpp
    # and is needed only by the binding translation unit, not by libsidre
    target_link_libraries(_sidre PRIVATE sidre conduit::conduit_python)

    # Place the built extension directly into the staged package tree so the
    # build tree is import-ready without an extra copy step.
    set_target_properties(_sidre PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${_axom_py_build_root}/axom/sidre")

    # Use HIP executable linker flags for the python module
    # (CMake treats modules separately from executables,
    # executable flags not automatically applied)
    if(AXOM_ENABLE_HIP)
        # Make CMake compile the binding file with the HIP compiler.
        set_source_files_properties(nanobind_sidre.cpp PROPERTIES LANGUAGE HIP)

        string (REPLACE " " ";" MODULE_LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
        target_link_options(_sidre PRIVATE ${MODULE_LINK_FLAGS})
    endif()

    # Stage the pure-Python package scaffolding into the build tree at configure
    # time (axom/ namespace root + py.typed, axom/sidre/ re-export, pysidre shim).
    axom_configure_file("${_pysidre_pkg_src}/axom/__init__.py"
                        "${_axom_py_build_root}/axom/__init__.py" COPYONLY)
    axom_configure_file("${_pysidre_pkg_src}/axom/py.typed"
                        "${_axom_py_build_root}/axom/py.typed" COPYONLY)
    axom_configure_file("${_pysidre_pkg_src}/axom/sidre/__init__.py"
                        "${_axom_py_build_root}/axom/sidre/__init__.py" COPYONLY)
    axom_configure_file("${_pysidre_pkg_src}/pysidre/__init__.py"
                        "${_axom_py_build_root}/pysidre/__init__.py" COPYONLY)

    # Type stubs (PEP 561). nanobind_add_stub imports the module by its bare name
    # ('import _sidre'), so the directory holding the built extension must be on
    # PYTHON_PATH, along with the module's runtime deps (conduit for Node interop,
    # numpy for ndarray returns). On an interpreter that already has conduit/numpy
    # these extra entries are harmless.
    set(_sidre_stub_pythonpath "${_axom_py_build_root}/axom/sidre")
    if(CONDUIT_PYTHON_MODULE_DIR)
        list(APPEND _sidre_stub_pythonpath ${CONDUIT_PYTHON_MODULE_DIR})
    endif()
    if(PY_NUMPY_DIR)
        list(APPEND _sidre_stub_pythonpath ${PY_NUMPY_DIR})
    endif()

    nanobind_add_stub(
        _sidre_stub
        MODULE _sidre
        OUTPUT "${_axom_py_build_root}/axom/sidre/_sidre.pyi"
        PYTHON_PATH ${_sidre_stub_pythonpath}
        DEPENDS _sidre)

    #--------------------------------------------------------------------------
    # Install the package tree into the site-packages-shaped prefix.
    #--------------------------------------------------------------------------
    install(TARGETS _sidre
            LIBRARY DESTINATION "${AXOM_PYTHON_MODULE_INSTALL_PREFIX}/axom/sidre")

    install(FILES "${_axom_py_build_root}/axom/sidre/_sidre.pyi"
                  "${_pysidre_pkg_src}/axom/sidre/__init__.py"
            DESTINATION "${AXOM_PYTHON_MODULE_INSTALL_PREFIX}/axom/sidre")

    # Namespace-root package files install once (not per component).
    install(FILES "${_pysidre_pkg_src}/axom/__init__.py"
                  "${_pysidre_pkg_src}/axom/py.typed"
            DESTINATION "${AXOM_PYTHON_MODULE_INSTALL_PREFIX}/axom")

    # Deprecation shim for the historical top-level 'pysidre' module.
    install(FILES "${_pysidre_pkg_src}/pysidre/__init__.py"
            DESTINATION "${AXOM_PYTHON_MODULE_INSTALL_PREFIX}/pysidre")
endif()


#------------------------------------------------------------------------------
# Add tests and examples
#------------------------------------------------------------------------------
if(AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

if(AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()
