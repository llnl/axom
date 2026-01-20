vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO llnl/raja
    REF v2025.12.0
    SHA512 3668df960d4f1f5fc7dc72e8e6b522a8e7b5da9c2fcefa5a975c2b311050bdf46a58662b9344522ac93765c83f9f63ee9f9d028d7a6ea13c61f72877a11b6929
    #PATCHES
    #    "./vcpkg_raja_openmp_forall.patch"
    #    "./vcpkg_raja_threadutils.patch"
)

set(_is_shared TRUE)
if(VCPKG_LIBRARY_LINKAGE STREQUAL "static")
    set(_is_shared FALSE)
else()
    set(_extra_cxx_flags "/DRAJASHAREDDLL_EXPORTS")
endif()

vcpkg_check_features(OUT_FEATURE_OPTIONS FEATURE_OPTIONS
    FEATURES
        openmp       ENABLE_OPENMP
)

vcpkg_configure_cmake(
    SOURCE_PATH ${SOURCE_PATH}
    PREFER_NINJA
    OPTIONS 
        -DBLT_SOURCE_DIR:PATH=${CURRENT_INSTALLED_DIR}/share/blt
        -Dcamp_DIR:PATH=${CURRENT_INSTALLED_DIR}
        -DENABLE_ALL_WARNINGS:BOOL=OFF
        -DENABLE_COVERAGE:BOOL=OFF
        -DENABLE_EXAMPLES:BOOL=OFF
        -DENABLE_TESTS:BOOL=OFF
        -DRAJA_ENABLE_RUNTIME_PLUGINS:BOOL=OFF
        -DRAJA_ENABLE_TESTS:BOOL=OFF
        -DRAJA_ENABLE_EXAMPLES:BOOL=OFF
        -DRAJA_ENABLE_EXERCISES:BOOL=OFF
        -DRAJA_ENABLE_REPRODUCERS:BOOL=OFF
        -DRAJA_ENABLE_DOCUMENTATION:BOOL=OFF
        -DRAJA_ENABLE_BENCHMARKS:BOOL=OFF
        -DBUILD_SHARED_LIBS:BOOL=${_is_shared}
        -DBLT_CXX_FLAGS:STRING=${_extra_cxx_flags}
        ${FEATURE_OPTIONS}
)

vcpkg_install_cmake()
vcpkg_fixup_cmake_targets(CONFIG_PATH lib/cmake/raja
                          TARGET_PATH share/raja)
vcpkg_copy_pdbs()


## shuffle the output directories to make vcpkg happy
# Remove extraneous debug header files
file(REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/debug/include)
file(REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/debug/share)

message(STATUS "CURRENT_PACKAGES_DIR -- ${CURRENT_PACKAGES_DIR}")

if(VCPKG_LIBRARY_LINKAGE STREQUAL static)
    # Note: Not tested
    file(REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/bin ${CURRENT_PACKAGES_DIR}/debug/bin)
else()
    set(_config_dir "${CURRENT_PACKAGES_DIR}/share/raja")

    # Move dll files from lib to bin directory
    file(MAKE_DIRECTORY ${CURRENT_PACKAGES_DIR}/bin )
    file(MAKE_DIRECTORY ${CURRENT_PACKAGES_DIR}/debug/bin )

    file(RENAME ${CURRENT_PACKAGES_DIR}/lib/RAJA.dll
                ${CURRENT_PACKAGES_DIR}/bin/RAJA.dll)

    file(RENAME ${CURRENT_PACKAGES_DIR}/debug/lib/RAJA.dll
                ${CURRENT_PACKAGES_DIR}/debug/bin/RAJA.dll)

    # Update paths to dlls in CMake config files
    foreach(_c  debug release)
        set(_f ${_config_dir}/RAJATargets-${_c}.cmake)
        file(READ ${_f} _fdata)
        string(REPLACE "lib/RAJA.dll" "bin/RAJA.dll" _fdata "${_fdata}")
        file(WRITE  ${_f} "${_fdata}")
    endforeach()
endif()

# Put the license file where vcpkg expects it
file(INSTALL ${SOURCE_PATH}/LICENSE DESTINATION ${CURRENT_PACKAGES_DIR}/share/raja RENAME copyright)
