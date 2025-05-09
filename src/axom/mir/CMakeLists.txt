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
                        COMPONENTS SLIC SLAM PRIMAL
                        TPLS       Conduit RAJA)

#------------------------------------------------------------------------------
# Specify all headers/sources
#------------------------------------------------------------------------------

# A serial reference implementation (remove eventually)
set(mir_reference_headers
    reference/MIRMesh.hpp
    reference/MIRMeshTypes.hpp
    reference/ZooClippingTables.hpp
    reference/InterfaceReconstructor.hpp
    reference/CellData.hpp
    reference/CellClipper.hpp
    reference/MIRUtilities.hpp
    reference/CellGenerator.hpp
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
    future/ClipFieldFilterDevice.hpp
    future/ClipFieldFilter.hpp
    )
set(mir_future_sources
    future/ClipFieldFilter.cpp
    )

set(mir_headers
    ${mir_reference_headers}
    MeshTester.hpp

    clipping/BlendGroupBuilder.hpp
    clipping/ClipCases.h
    clipping/ClipField.hpp
    clipping/ClipOptions.hpp
    clipping/ClipTableManager.hpp
    detail/elvira_impl.hpp
    detail/elvira_detail.hpp
    ElviraAlgorithm.hpp
    EquiZAlgorithm.hpp
    MIRAlgorithm.hpp
    MIROptions.hpp
    Options.hpp
    utilities/blueprint_utilities.hpp
    utilities/BlendData.hpp
    utilities/CoordsetBlender.hpp
    utilities/CoordsetSlicer.hpp
    utilities/ExtractZones.hpp
    utilities/ExtrudeMesh.hpp
    utilities/FieldBlender.hpp
    utilities/FieldSlicer.hpp
    utilities/MakeZoneCenters.hpp
    utilities/MakeZoneVolumes.hpp
    utilities/MakeUnstructured.hpp
    utilities/MatsetSlicer.hpp
    utilities/MergeMeshes.hpp
    utilities/NodeToZoneRelationBuilder.hpp
    utilities/PrimalAdaptor.hpp
    utilities/RecenterField.hpp
    utilities/SelectedZones.hpp
    utilities/TopologyMapper.hpp
    utilities/utilities.hpp
    utilities/VariableShape.hpp
    utilities/ZoneListBuilder.hpp
    views/dispatch_coordset.hpp
    views/dispatch_material.hpp
    views/dispatch_rectilinear_topology.hpp
    views/dispatch_structured_topology.hpp
    views/dispatch_topology.hpp
    views/dispatch_uniform_topology.hpp
    views/dispatch_unstructured_topology.hpp
    views/dispatch_utilities.hpp
    views/ExplicitCoordsetView.hpp
    views/MaterialView.hpp
    views/NodeArrayView.hpp
    views/RectilinearCoordsetView.hpp
    views/Shapes.hpp
    views/StridedStructuredIndexing.hpp
    views/StructuredIndexing.hpp
    views/StructuredTopologyView.hpp
    views/UniformCoordsetView.hpp
    views/UnstructuredTopologyMixedShapeView.hpp
    views/UnstructuredTopologyPolyhedralView.hpp
    views/UnstructuredTopologySingleShapeView.hpp
    views/view_traits.hpp
    )

set(mir_sources
    ${mir_reference_sources}

    clipping/ClipCasesHex.cpp
    clipping/ClipCasesPyr.cpp
    clipping/ClipCasesQua.cpp
    clipping/ClipCasesTet.cpp
    clipping/ClipCasesTri.cpp
    clipping/ClipCasesWdg.cpp
    MeshTester.cpp
    MIRAlgorithm.cpp
    utilities/blueprint_utilities.cpp
    views/MaterialView.cpp
    views/UnstructuredTopologyMixedShapeView.cpp
    )

#------------------------------------------------------------------------------
# Build and install the library
#------------------------------------------------------------------------------
set(mir_depends_on core slic slam primal conduit::conduit RAJA)

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


