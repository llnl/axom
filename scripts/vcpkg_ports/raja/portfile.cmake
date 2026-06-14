vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO llnl/raja
    REF v2025.12.1
    SHA512 fca1d5d336cb552bbee1f33da69734c3c0dd53632e04db3ae76f08e1160801676e085c41b52fa1d1e7bd91b8b880338c30396a03c13dd32f6941c327c52cf5e8
    PATCHES
        cuda-13-mem-advise-location.patch
        windows-cuda-msvc.patch
)

set(_extra_cxx_flags "")
set(_is_shared TRUE)
set(FEATURE_OPTIONS_DEBUG "")
set(FEATURE_OPTIONS_RELEASE "")

vcpkg_check_features(OUT_FEATURE_OPTIONS FEATURE_OPTIONS
    FEATURES
        cuda         ENABLE_CUDA
        openmp       ENABLE_OPENMP
)

if(VCPKG_TARGET_IS_WINDOWS AND "cuda" IN_LIST FEATURES)
    vcpkg_check_linkage(ONLY_STATIC_LIBRARY)
endif()

if(VCPKG_LIBRARY_LINKAGE STREQUAL "static" OR (VCPKG_TARGET_IS_WINDOWS AND "cuda" IN_LIST FEATURES))
    set(_is_shared FALSE)
endif()

if("cuda" IN_LIST FEATURES)
    include("${CURRENT_INSTALLED_DIR}/share/cuda/vcpkg_find_cuda.cmake")
    vcpkg_find_cuda(OUT_CUDA_TOOLKIT_ROOT _cuda_toolkit_root)

    if((NOT DEFINED CUDA_ARCHITECTURES OR "${CUDA_ARCHITECTURES}" STREQUAL "") AND DEFINED ENV{CUDA_ARCHITECTURES})
        set(CUDA_ARCHITECTURES "$ENV{CUDA_ARCHITECTURES}")
    endif()

    if(VCPKG_TARGET_IS_WINDOWS)
        string(APPEND _extra_cxx_flags " /Zc:preprocessor /utf-8")
        list(APPEND FEATURE_OPTIONS_DEBUG
            "-DBLT_CUDA_FLAGS:STRING=-Xcompiler=/Zc:preprocessor -Xcompiler=/utf-8 -Xcompiler=/MDd"
        )
        list(APPEND FEATURE_OPTIONS_RELEASE
            "-DBLT_CUDA_FLAGS:STRING=-Xcompiler=/Zc:preprocessor -Xcompiler=/utf-8 -Xcompiler=/MD"
        )
    endif()

    if(DEFINED CUDA_ARCHITECTURES AND NOT "${CUDA_ARCHITECTURES}" STREQUAL "")
        list(APPEND FEATURE_OPTIONS
            "-DCMAKE_CUDA_ARCHITECTURES=${CUDA_ARCHITECTURES}"
        )
    endif()

    if(EXISTS "${_cuda_toolkit_root}/include/cccl/cub/cub.cuh")
        list(APPEND FEATURE_OPTIONS
            "-DCUB_DIR=${_cuda_toolkit_root}/include/cccl"
        )
    elseif(EXISTS "${_cuda_toolkit_root}/include/cub/cub.cuh")
        list(APPEND FEATURE_OPTIONS
            "-DCUB_DIR=${_cuda_toolkit_root}/include"
        )
    endif()
endif()

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
    OPTIONS_DEBUG
        ${FEATURE_OPTIONS_DEBUG}
    OPTIONS_RELEASE
        ${FEATURE_OPTIONS_RELEASE}
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

if(NOT _is_shared)
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
