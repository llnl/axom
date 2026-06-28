vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO llnl/conduit
    REF v0.9.5
    SHA512 c0b8e92e6eb3c42a2efc8656c488b62c2fbacbe015df7820d8eae8752deea8271db65f2d129df251a7eb97402352ff2260e3ac8e38a3c37c8606ab1096511757
    HEAD_REF develop
    PATCHES 
        "./setup_deps_vcpkg_triplet.patch"
        "./python_install_prefix_slashes.patch"
)

set(_python_options -DENABLE_PYTHON=OFF)
if("python" IN_LIST FEATURES)
    include("${CURRENT_INSTALLED_DIR}/share/python3/vcpkg-port-config.cmake")

    # Conduit needs numpy while building its python module.
    # Use vcpkg's helper python environment for those build-time packages.
    vcpkg_get_vcpkg_installed_python(_conduit_python)
    if(CMAKE_HOST_WIN32)
        vcpkg_execute_required_process(
            COMMAND "${_conduit_python}" -I -m ensurepip --upgrade
            WORKING_DIRECTORY "${CURRENT_BUILDTREES_DIR}"
            LOGNAME "ensurepip-conduit-python-${TARGET_TRIPLET}")
        vcpkg_execute_required_process(
            COMMAND "${_conduit_python}" -I -m pip install
                    --disable-pip-version-check
                    --no-warn-script-location
                    virtualenv
            WORKING_DIRECTORY "${CURRENT_BUILDTREES_DIR}"
            LOGNAME "pip-install-conduit-virtualenv-${TARGET_TRIPLET}")
    endif()
    x_vcpkg_get_python_packages(
        PYTHON_VERSION "3"
        PYTHON_EXECUTABLE "${_conduit_python}"
        PACKAGES
            "numpy==2.4.2"
            "setuptools==80.9.0"
        OUT_PYTHON_VAR _conduit_python_with_packages)

    set(_conduit_config_python "${_conduit_python_with_packages}")
    if(VCPKG_TARGET_IS_WINDOWS)
        # On Windows, configure Conduit with vcpkg's installed python.exe
        # so the generated config matches the interpreter Axom will use later.
        # Provide the helper environment's numpy path through PYTHONPATH for configure.
        set(_conduit_config_python "${CURRENT_INSTALLED_DIR}/tools/python3/python.exe")
        execute_process(
            COMMAND "${_conduit_python_with_packages}" -I -c
                    "import pathlib, numpy; print(pathlib.Path(numpy.__file__).resolve().parents[1])"
            OUTPUT_VARIABLE _conduit_pythonpath
            RESULT_VARIABLE _conduit_pythonpath_result
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(NOT _conduit_pythonpath_result EQUAL 0)
            message(FATAL_ERROR "Failed to determine Conduit build-time Python package path.")
        endif()
        set(ENV{PYTHONPATH} "${_conduit_pythonpath}")
    endif()

    set(_conduit_python_library "")
    if(VCPKG_TARGET_IS_WINDOWS)
        set(_conduit_python_library
            "${CURRENT_INSTALLED_DIR}/lib/python${PYTHON3_VERSION_MAJOR}${PYTHON3_VERSION_MINOR}.lib")
    endif()

    set(_python_module_install_prefix "${CURRENT_PACKAGES_DIR}/${PYTHON3_SITE}")

    set(_python_options
        -DENABLE_PYTHON=ON
        -DPYTHON_EXECUTABLE=${_conduit_config_python}
        -DPYTHON_MODULE_INSTALL_PREFIX=${_python_module_install_prefix})

    if(_conduit_python_library)
        list(APPEND _python_options -DPYTHON_LIBRARY=${_conduit_python_library})
    endif()
endif()

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}/src"
    OPTIONS 
        -DBLT_SOURCE_DIR=${CURRENT_HOST_INSTALLED_DIR}/share/blt
        -DCONDUIT_ENABLE_TESTS=OFF
        -DENABLE_COVERAGE=OFF
        -DENABLE_DOCS=OFF
        -DENABLE_EXAMPLES=OFF
        ${_python_options}
        -DENABLE_TESTS=OFF
        -DENABLE_UTILS=OFF
        -DHDF5_DIR=${CURRENT_INSTALLED_DIR}
        -DCONDUIT_INSTALL_CONFIG_DIR="share/conduit"
        -DCONDUIT_INSTALL_CMAKE_MODULE_DIR="share"
)

vcpkg_cmake_install()
vcpkg_cmake_config_fixup(
        CONFIG_PATH  lib/cmake/conduit
        TOOLS_PATH   tools/conduit)
vcpkg_copy_pdbs()

if("python" IN_LIST FEATURES)
    # The helper environment above is only for building Conduit.
    # Also stage numpy into vcpkg's installed python site so downstream imports work.
    vcpkg_execute_required_process(
        COMMAND "${_conduit_python_with_packages}" -I -m pip install
                --disable-pip-version-check
                --no-warn-script-location
                --target "${CURRENT_PACKAGES_DIR}/${PYTHON3_SITE}"
                "numpy==2.4.2"
        WORKING_DIRECTORY "${CURRENT_BUILDTREES_DIR}"
        LOGNAME "pip-install-conduit-python-packages-${TARGET_TRIPLET}")

    # NumPy's nested pkgconfig file is not useful to vcpkg consumers and trips
    # vcpkg's pkgconfig layout audit.
    file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/${PYTHON3_SITE}/numpy/_core/lib/pkgconfig")

    if(VCPKG_TARGET_IS_WINDOWS)
        # Imported extension modules need vcpkg's bin directory in the Windows
        # DLL search path; keep the handles alive on os to preserve the entries.
        file(WRITE "${CURRENT_PACKAGES_DIR}/${PYTHON3_SITE}/conduit_vcpkg_dll_path.pth"
             "import os,sys,pathlib; os._vcpkg_dll_dirs=getattr(os,'_vcpkg_dll_dirs',[])+[os.add_dll_directory(str(p/'bin')) for r in {pathlib.Path(sys.prefix).resolve(), pathlib.Path(sys.base_prefix).resolve()} for p in (r,*r.parents) if hasattr(os,'add_dll_directory') and (p/'bin').is_dir()]\n")
    endif()

    # Conduit exports absolute package-staging paths for custom python installs.
    # Rewrite them so consumers resolve modules from the installed vcpkg prefix.
    set(_conduit_config "${CURRENT_PACKAGES_DIR}/share/conduit/ConduitConfig.cmake")
    file(READ "${_conduit_config}" _conduit_config_contents)
    string(REGEX REPLACE
        "set\\(CONDUIT_INSTALL_PREFIX \"[^\"]*\"\\)"
        "set(CONDUIT_INSTALL_PREFIX \"\${VCPKG_IMPORT_PREFIX}\")"
        _conduit_config_contents
        "${_conduit_config_contents}")
    file(WRITE "${_conduit_config}" "${_conduit_config_contents}")
    vcpkg_replace_string("${_conduit_config}"
        "set(CONDUIT_PYTHON_MODULE_DIR \"${CURRENT_PACKAGES_DIR}/${PYTHON3_SITE}\")"
        "set(CONDUIT_PYTHON_MODULE_DIR \"${PYTHON3_SITE}\")"
        IGNORE_UNCHANGED)
    vcpkg_replace_string("${_conduit_config}"
        "set(CONDUIT_PYTHON_MODULE_DIR \"\${VCPKG_IMPORT_PREFIX}/${PYTHON3_SITE}\")"
        "set(CONDUIT_PYTHON_MODULE_DIR \"${PYTHON3_SITE}\")"
        IGNORE_UNCHANGED)
    vcpkg_replace_string("${_conduit_config}"
        "set(CONDUIT_PYTHON_MODULE_CUSTOM_PREFIX \"ON\")"
        "set(CONDUIT_PYTHON_MODULE_CUSTOM_PREFIX \"OFF\")"
        IGNORE_UNCHANGED)
endif()


# Conduit's generated runner and makefile config embed package/build paths.
# Axom uses its own Python wrapper and CMake package config instead.
file(REMOVE
    "${CURRENT_PACKAGES_DIR}/bin/run_python_with_conduit.sh"
    "${CURRENT_PACKAGES_DIR}/debug/bin/run_python_with_conduit.sh"
    "${CURRENT_PACKAGES_DIR}/share/conduit/conduit_config.mk")

file(REMOVE_RECURSE
    "${CURRENT_PACKAGES_DIR}/debug/include"
    "${CURRENT_PACKAGES_DIR}/debug/share")

# Remove exe files -- vcpkg doesn't like them
# (Future): It might be possible to move them to the vcpkg 'tools' directory
foreach(_dir "bin" "debug/bin")
    file(GLOB _exe "${CURRENT_PACKAGES_DIR}/${_dir}/*.exe")
    if(_exe)
        file(REMOVE ${_exe})
    endif()
endforeach()

if(VCPKG_LIBRARY_LINKAGE STREQUAL static)
    file(REMOVE_RECURSE
        "${CURRENT_PACKAGES_DIR}/bin"
        "${CURRENT_PACKAGES_DIR}/debug/bin")
endif()

vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE")
