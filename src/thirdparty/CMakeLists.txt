# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Axom third-party builtin libs are all currently header only.
# We have tests here that make sure they link and can be used.
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Register the sol library for parsing lua
#------------------------------------------------------------------------------

if (LUA_FOUND)
  blt_add_library(NAME    sol
                  HEADERS axom/sol.hpp axom/sol_forward.hpp)

  target_include_directories(sol SYSTEM INTERFACE
              $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/thirdparty>)

  install(TARGETS              sol
          EXPORT               axom-targets
          INCLUDES DESTINATION include)
  install(FILES        ${PROJECT_SOURCE_DIR}/thirdparty/axom/sol.hpp
                       ${PROJECT_SOURCE_DIR}/thirdparty/axom/sol_forward.hpp
          DESTINATION include/axom )

  set(SOL_FOUND TRUE CACHE INTERNAL "")
else()
  set(SOL_FOUND FALSE CACHE INTERNAL "")
endif()
mark_as_advanced(SOL_FOUND)

message(STATUS "SOL Enabled: ${SOL_FOUND}")

#------------------------------------------------------------------------------
# Register the built-in CLI11 library for command line option parsing
#------------------------------------------------------------------------------

blt_add_library(NAME    cli11
                HEADERS axom/CLI11.hpp)

target_include_directories(cli11 SYSTEM INTERFACE
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/thirdparty>)

install(TARGETS              cli11
        EXPORT               axom-targets
        INCLUDES DESTINATION include)
install(FILES       ${PROJECT_SOURCE_DIR}/thirdparty/axom/CLI11.hpp
        DESTINATION include/axom )

set(CLI11_FOUND TRUE CACHE INTERNAL "")
mark_as_advanced(CLI11_FOUND)

#------------------------------------------------------------------------------
# Register the built-in fmt library for string formatting
#------------------------------------------------------------------------------

# NOTE: When updating fmt you must do the following changes to properly guard
#  our fmt from other fmt's.
# 
# 1) Change the namespace from 'fmt' to 'axom::fmt'.  This is in fmt/core.h:203-204,
#    as well as, search and replace `fmt::` to `axom::fmt`.
# 2) Search and replace compiler define's prefix from `FMT_` to `AXOM_FMT_`.
# 3) Keep fmt's headers installed under axom's include directory.

# Headers that Axom adds to fmt
set(axom_fmt_headers
    axom/fmt.hpp
    )

# Headers that come from fmt
set(fmt_headers
    axom/fmt/args.h
    axom/fmt/chrono.h
    axom/fmt/color.h
    axom/fmt/compile.h
    axom/fmt/core.h
    axom/fmt/format.h
    axom/fmt/format-inl.h
    axom/fmt/os.h
    axom/fmt/ostream.h
    axom/fmt/printf.h
    axom/fmt/ranges.h
    axom/fmt/std.h
    axom/fmt/xchar.h
    )

blt_add_library(NAME    fmt
                HEADERS ${fmt_headers} ${axom_fmt_headers})

target_include_directories(fmt SYSTEM INTERFACE
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/thirdparty>)

target_include_directories(fmt INTERFACE 
            $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include> )

# Setup some variables for fmt in Axom's config.hpp
set(AXOM_FMT_EXCEPTIONS FALSE PARENT_SCOPE)
set(AXOM_FMT_HEADER_ONLY TRUE PARENT_SCOPE)

install(TARGETS              fmt
        EXPORT               axom-targets
        INCLUDES DESTINATION include)
# Install our helper header to axom
install(FILES   ${axom_fmt_headers}
        DESTINATION include/axom
        )
# Install fmt's headers to axom/fmt
install(DIRECTORY   ${PROJECT_SOURCE_DIR}/thirdparty/axom/fmt
        DESTINATION include/axom
        )

set(FMT_FOUND TRUE CACHE INTERNAL "")
mark_as_advanced(FMT_FOUND)

#------------------------------------------------------------------------------
# Register the built-in sparse hash library
#------------------------------------------------------------------------------

if (AXOM_ENABLE_SPARSEHASH)

    set(SPARSEHASH_HASHFUN_HEADER <functional> CACHE INTERNAL "")
    set(SPARSEHASH_HASHFUN_NAMESPACE std CACHE INTERNAL "")
    set(SPARSEHASH_FOUND TRUE CACHE INTERNAL "")

    mark_as_advanced(
      SPARSEHASH_FOUND 
      SPARSEHASH_HASHFUN_HEADER 
      SPARSEHASH_HASHFUN_NAMESPACE
    )
    
    # Explicitly list sparsehash files for copy-headers,
    # Separate into two lists to maintain sparsehash's hierarchy
    set(sparsehash_headers
      axom/sparsehash/README
      axom/sparsehash/dense_hash_map
      axom/sparsehash/dense_hash_set
      axom/sparsehash/sparse_hash_map
      axom/sparsehash/sparse_hash_set
      axom/sparsehash/sparsetable
      axom/sparsehash/template_util.h
      axom/sparsehash/type_traits.h     
     )

    set(sparsehash_headers_internal
      axom/sparsehash/internal/densehashtable.h
      axom/sparsehash/internal/hashtable-common.h
      axom/sparsehash/internal/libc_allocator_with_realloc.h
      axom/sparsehash/internal/sparseconfig.h
      axom/sparsehash/internal/sparsehashtable.h
     )
     
    blt_add_library(NAME sparsehash
                    HEADERS ${sparsehash_headers})
     
    target_include_directories(sparsehash SYSTEM INTERFACE
                $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/thirdparty>)

    target_include_directories(sparsehash INTERFACE 
                $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include> )

    # Disable warning introduced in gcc@8.1+ related to how sparsehash casts memory
    # The double guarding for compiler family helps when compiling libraries with GNU
    # then compiling the main code with another compiler (this is happens for example when
    # running clang-query)
    if(C_COMPILER_FAMILY_IS_GNU AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "8.1")
        blt_add_target_compile_flags(TO sparsehash FLAGS
            $<$<AND:$<CXX_COMPILER_ID:GNU>,$<COMPILE_LANGUAGE:CXX>>:-Wno-class-memaccess>)
    endif()

    install(TARGETS              sparsehash
            EXPORT               axom-targets
            INCLUDES DESTINATION include)

    install(DIRECTORY   ${PROJECT_SOURCE_DIR}/thirdparty/axom/sparsehash
            DESTINATION include/axom)
endif()

#------------------------------------------------------------------------------
# Add tests
#------------------------------------------------------------------------------
if (AXOM_ENABLE_TESTS)
    add_subdirectory(tests)
endif()
