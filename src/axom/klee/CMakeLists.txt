# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
axom_component_requires(NAME       Klee
                        COMPONENTS PRIMAL SLIC INLET SIDRE)

#------------------------------------------------------------------------------
# Specify all headers/sources
#------------------------------------------------------------------------------
set(klee_headers
    AffineMatrixVisitor.hpp
    Dimensions.hpp
    Geometry.hpp
    GeometryOperators.hpp
    Shape.hpp
    ShapeSet.hpp
    Units.hpp
    KleeError.hpp
    io/IO.hpp
    )

set(klee_internal_headers
    io/GeometryOperatorsIO.hpp
    io/IOUtil.hpp
    )

set(klee_sources
    AffineMatrixVisitor.cpp
    Geometry.cpp
    GeometryOperators.cpp
    KleeError.cpp
    Shape.cpp
    ShapeSet.cpp
    Units.cpp
    io/GeometryOperatorsIO.cpp
    io/IO.cpp
    io/IOUtil.cpp
    )

#------------------------------------------------------------------------------
# Make/Install the library
#------------------------------------------------------------------------------
axom_add_library(NAME        klee
                 SOURCES     ${klee_sources}
                 HEADERS     ${klee_headers} ${klee_internal_headers}
                 DEPENDS_ON  inlet
                 FOLDER      axom/klee
                 )

axom_write_unified_header( NAME klee
                           HEADERS ${klee_headers}
                           )

axom_install_component( NAME    klee
                        HEADERS ${klee_headers}
                        )

if(NANOBIND_FOUND)
    nanobind_add_module(pyklee nanobind_klee.cpp)
    target_link_libraries(pyklee PRIVATE klee inlet sidre primal fmt)

    if(AXOM_ENABLE_HIP)
        set_source_files_properties(nanobind_klee.cpp PROPERTIES LANGUAGE HIP)

        string(REPLACE " " ";" MODULE_LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
        target_link_options(pyklee PRIVATE ${MODULE_LINK_FLAGS})
    endif()

    install(TARGETS pyklee LIBRARY DESTINATION lib)
endif()


#------------------------------------------------------------------------------
# Add tests
#------------------------------------------------------------------------------
if (AXOM_ENABLE_TESTS AND ENABLE_GMOCK)
  add_subdirectory(tests)
endif()
