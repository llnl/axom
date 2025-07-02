# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
axom_component_requires(NAME       Klee
                        COMPONENTS PRIMAL SLIC INLET SIDRE)

#------------------------------------------------------------------------------
# Specify all headers/sources
#------------------------------------------------------------------------------
set(klee_headers
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


#------------------------------------------------------------------------------
# Add tests
#------------------------------------------------------------------------------
if (AXOM_ENABLE_TESTS AND ENABLE_GMOCK)
  add_subdirectory(tests)
endif()
