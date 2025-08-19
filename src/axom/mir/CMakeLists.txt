# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# MIR -- Material Interface Reconstruction
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Check necessary dependencies
#------------------------------------------------------------------------------
axom_component_requires(NAME       MIR
                        COMPONENTS SLIC SLAM PRIMAL BUMP
                        TPLS       Conduit)

#------------------------------------------------------------------------------
# Specify all headers/sources
#------------------------------------------------------------------------------

# A serial reference implementation (remove eventually)
set(mir_reference_headers
    reference/CellClipper.hpp
    reference/CellData.hpp
    reference/CellGenerator.hpp
    reference/InterfaceReconstructor.hpp
    reference/MIRMesh.hpp
    reference/MIRMeshTypes.hpp
    reference/MIRUtilities.hpp
    reference/ZooClippingTables.hpp
    )
set(mir_reference_sources
    reference/MIRMesh.cpp
    reference/InterfaceReconstructor.cpp
    reference/ZooClippingTables.cpp
    reference/CellData.cpp
    reference/CellClipper.cpp
    reference/CellGenerator.cpp
   )

# Future classes (evaluate whether to keep them)
set(mir_future_headers
    future/ClipFieldFilter.hpp
    future/ClipFieldFilterDevice.hpp
    )
set(mir_future_sources
    future/ClipFieldFilter.cpp
    )

set(mir_headers
    ${mir_reference_headers}
    detail/elvira_detail.hpp
    detail/elvira_impl.hpp
    detail/equiz_detail.hpp
    ElviraAlgorithm.hpp
    EquiZAlgorithm.hpp
    MIRAlgorithm.hpp
    )

set(mir_sources
    ${mir_reference_sources}
    MIRAlgorithm.cpp
    )

#------------------------------------------------------------------------------
# Build and install the library
#------------------------------------------------------------------------------
set(mir_depends_on core bump slic slam primal conduit::conduit)
blt_list_append( TO mir_depends_on ELEMENTS RAJA IF RAJA_FOUND )

axom_add_library(
    NAME        mir
    SOURCES     ${mir_sources}
    HEADERS     ${mir_headers}
    DEPENDS_ON  ${mir_depends_on}
    FOLDER      axom/mir)

axom_write_unified_header(NAME    mir
                          HEADERS ${mir_headers} )

axom_install_component(NAME    mir
                       HEADERS ${mir_headers} )

#------------------------------------------------------------------------------
# Add tests, benchmarks and examples
#------------------------------------------------------------------------------
if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

if (AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()


