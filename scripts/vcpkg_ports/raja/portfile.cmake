vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO llnl/raja
    REF b92047352b30ebd87012eb95da5fe70a2927ff9f
    SHA512 c3fd435bab2999d06e44b6b2a3a2dfdbe9bb583b880c186de40a14c4c07f7ef7a376230e809b427b8d3b8c158438b11a0a9f860f0fbf5745bf46d2fab8384ff5
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
