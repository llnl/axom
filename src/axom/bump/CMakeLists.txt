# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# BUMP -- Blueprint Utilities for Mesh Processing
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Check necessary dependencies
#------------------------------------------------------------------------------
axom_component_requires(NAME       BUMP
                        COMPONENTS SLIC SPIN PRIMAL
                        TPLS       Conduit)

#------------------------------------------------------------------------------
# Specify all headers/sources
#------------------------------------------------------------------------------

set(bump_headers
    clipping/BlendGroupBuilder.hpp
    clipping/ClipCases.h
    clipping/ClipField.hpp
    clipping/ClipOptions.hpp
    clipping/ClipTableManager.hpp
    data/MeshTester.hpp
    io/save.hpp
    utilities/blueprint_utilities.hpp
    utilities/conduit_memory.hpp
    utilities/conduit_traits.hpp
    utilities/utilities.hpp
    views/BasicIndexing.hpp
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
    BlendData.hpp
    CoordsetBlender.hpp
    CoordsetExtents.hpp
    CoordsetSlicer.hpp
    ExtractZonesAndMatsetPolyhedral.hpp
    ExtractZones.hpp
    ExtrudeMesh.hpp
    FieldBlender.hpp
    FieldSlicer.hpp
    HashNaming.hpp
    IndexingPolicies.hpp
    MakePointMesh.hpp
    MakePolyhedralTopology.hpp
    MakeUnstructured.hpp
    MakeZoneCenters.hpp
    MakeZoneVolumes.hpp
    MatsetSlicer.hpp
    MergeCoordsetPoints.hpp
    MergeMeshes.hpp
    MergePolyhedralFaces.hpp
    MinMax.hpp
    NodeToZoneRelationBuilder.hpp
    Options.hpp
    PrimalAdaptor.hpp
    RecenterField.hpp
    SelectedZones.hpp
    TopologyMapper.hpp
    Unique.hpp
    VariableShape.hpp
    ZoneListBuilder.hpp
    )

set(bump_sources
    clipping/ClipCasesHex.cpp
    clipping/ClipCasesPyr.cpp
    clipping/ClipCasesQua.cpp
    clipping/ClipCasesTet.cpp
    clipping/ClipCasesTri.cpp
    clipping/ClipCasesWdg.cpp
    data/MeshTester.cpp
    io/save.cpp
    utilities/conduit_traits.cpp
    utilities/blueprint_utilities.cpp
    views/MaterialView.cpp
    views/UnstructuredTopologyMixedShapeView.cpp
    )

#------------------------------------------------------------------------------
# Build and install the library
#------------------------------------------------------------------------------
set(bump_depends_on core slic spin primal conduit::conduit)
blt_list_append( TO bump_depends_on ELEMENTS RAJA IF RAJA_FOUND )

axom_add_library(
    NAME        bump
    SOURCES     ${bump_sources}
    HEADERS     ${bump_headers}
    DEPENDS_ON  ${bump_depends_on}
    FOLDER      axom/bump)

axom_write_unified_header(NAME    bump
                          HEADERS ${bump_headers} )

axom_install_component(NAME    bump
                       HEADERS ${bump_headers} )

#------------------------------------------------------------------------------
# Add tests, benchmarks and examples
#------------------------------------------------------------------------------
if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

