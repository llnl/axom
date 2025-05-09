# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# MultiMat component
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Check necessary dependencies
#------------------------------------------------------------------------------
axom_component_requires(NAME       MultiMat
                        COMPONENTS SLAM SLIC)

#------------------------------------------------------------------------------
# Specify all headers/sources
#------------------------------------------------------------------------------
set(multimat_headers
    multimat.hpp
    mmfield.hpp
    mmsubfield.hpp
    )

set(multimat_sources
    multimat.cpp
    )

#------------------------------------------------------------------------------
# Specify multimat dependencies
#------------------------------------------------------------------------------
set(multimat_depends_on slam)

blt_list_append( TO multimat_depends_on ELEMENTS RAJA IF RAJA_FOUND )
blt_list_append( TO multimat_depends_on ELEMENTS umpire IF UMPIRE_FOUND )

#------------------------------------------------------------------------------
# Build and install the library
#------------------------------------------------------------------------------
axom_add_library(
    NAME        multimat
    SOURCES     ${multimat_sources}
    HEADERS     ${multimat_headers}
    DEPENDS_ON  ${multimat_depends_on}
    FOLDER      axom/multimat
    )

axom_write_unified_header(NAME    multimat
                          HEADERS ${multimat_headers}
                          )

axom_install_component(NAME    multimat
                       HEADERS ${multimat_headers}
                       )

#------------------------------------------------------------------------------
# Add tests, benchmarks and examples
#------------------------------------------------------------------------------
if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

if (AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()
