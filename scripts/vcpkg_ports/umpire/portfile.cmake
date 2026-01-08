vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO llnl/umpire
    REF v2025.12.0
    SHA512 3dbb321b3b79ae60ecad4ef81ee40e98e41b370da77dfbd9bd037b118d05da9ea67a8aabf82f60d4302ab40897a0f58b75ad910cce3ec831984b6bd47530760b
    HEAD_REF develop
)


vcpkg_check_features(OUT_FEATURE_OPTIONS FEATURE_OPTIONS
    FEATURES
        openmp       ENABLE_OPENMP
)

set(_is_shared TRUE)
if(VCPKG_LIBRARY_LINKAGE STREQUAL "static")
    set(_is_shared FALSE)
else()
    list(APPEND FEATURE_OPTIONS -DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=ON)
endif()


vcpkg_configure_cmake(
    SOURCE_PATH ${SOURCE_PATH}
    PREFER_NINJA
    OPTIONS 
        -DBLT_SOURCE_DIR:PATH=${CURRENT_INSTALLED_DIR}/share/blt
        -Dcamp_DIR:PATH=${CURRENT_INSTALLED_DIR}
        -Dfmt_DIR:PATH=${CURRENT_INSTALLED_DIR}
        -DENABLE_ALL_WARNINGS:BOOL=OFF
        -DENABLE_WARNINGS_AS_ERRORS:BOOL=OFF
        -DENABLE_COVERAGE:BOOL=OFF
        -DENABLE_EXAMPLES:BOOL=OFF
        -DENABLE_TESTS:BOOL=OFF
        -DENABLE_BENCHMARKS:BOOL=OFF
        -DUMPIRE_ENABLE_FILESYSTEM:BOOL=ON
        -DBLT_CXX_STD:STRING=c++17
        -DBLT_OPENMP_LINK_FLAGS:STRING=" "
        -DUMPIRE_ENABLE_TOOLS:BOOL=OFF
        -DUMPIRE_ENABLE_TESTS:BOOL=OFF
        -DUMPIRE_ENABLE_BENCHMARKS:BOOL=OFF
        -DUMPIRE_ENABLE_DEVELOPER_BENCHMARKS:BOOL=OFF
        -DBUILD_SHARED_LIBS:BOOL=${_is_shared}
        -DOpenMP_RUNTIME_MSVC:STRING="experimental"
        ${FEATURE_OPTIONS}
)

vcpkg_install_cmake()
vcpkg_fixup_cmake_targets(CONFIG_PATH lib/cmake/umpire
                          TARGET_PATH share/umpire)
vcpkg_copy_pdbs()


## shuffle the output directories to make vcpkg happy
# Remove extraneous debug header files
file(REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/debug/include)
file(REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/debug/share)

message(STATUS "CURRENT_PACKAGES_DIR -- ${CURRENT_PACKAGES_DIR}")

if(VCPKG_LIBRARY_LINKAGE STREQUAL static)
    # Note: Not tested
    file(REMOVE_RECURSE ${CURRENT_PACKAGES_DIR}/bin ${CURRENT_PACKAGES_DIR}/debug/bin)
endif()


# Put the license file where vcpkg expects it
file(INSTALL ${SOURCE_PATH}/LICENSE DESTINATION ${CURRENT_PACKAGES_DIR}/share/umpire RENAME copyright)
