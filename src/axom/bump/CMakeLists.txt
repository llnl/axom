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
    extraction/BlendGroupBuilder.hpp
    extraction/tables/clipping/ClipCases.h
    extraction/tables/cutting/CutCases.hpp
    extraction/ClipField.hpp
    extraction/CutField.hpp
    extraction/ExtractorOptions.hpp
    extraction/FieldOptions.hpp
    extraction/ClipTableManager.hpp
    extraction/CutTableManager.hpp
    extraction/ExtractionConstants.hpp
    extraction/FieldIntersector.hpp
    extraction/PlaneIntersector.hpp
    extraction/PlaneSlice.hpp
    extraction/Table.hpp
    extraction/TableBasedExtractor.hpp
    extraction/TableManager.hpp
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
    ExtractZones.hpp
    ExtractZonesAndMatsetPolyhedral.hpp
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
    MapBasedNaming.hpp
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
    extraction/ClipTableManager.cpp
    extraction/CutTableManager.cpp
    extraction/Table.cpp
    extraction/TableManager.cpp
    extraction/tables/clipping/ClipCasesHex.cpp
    extraction/tables/clipping/ClipCasesPoly5.cpp
    extraction/tables/clipping/ClipCasesPoly6.cpp
    extraction/tables/clipping/ClipCasesPoly7.cpp
    extraction/tables/clipping/ClipCasesPoly8.cpp
    extraction/tables/clipping/ClipCasesPyr.cpp
    extraction/tables/clipping/ClipCasesQua.cpp
    extraction/tables/clipping/ClipCasesTet.cpp
    extraction/tables/clipping/ClipCasesTri.cpp
    extraction/tables/clipping/ClipCasesWdg.cpp
    extraction/tables/cutting/CutCasesHex.cpp
    extraction/tables/cutting/CutCasesPoly5.cpp
    extraction/tables/cutting/CutCasesPoly6.cpp
    extraction/tables/cutting/CutCasesPoly7.cpp
    extraction/tables/cutting/CutCasesPoly8.cpp
    extraction/tables/cutting/CutCasesPyramid.cpp
    extraction/tables/cutting/CutCasesQua.cpp
    extraction/tables/cutting/CutCasesTet.cpp
    extraction/tables/cutting/CutCasesTri.cpp
    extraction/tables/cutting/CutCasesWedge.cpp
    data/MeshTester.cpp
    io/save.cpp
    utilities/conduit_traits.cpp
    utilities/blueprint_utilities.cpp
    views/MaterialView.cpp
    views/UnstructuredTopologyMixedShapeView.cpp
    views/dispatch_utilities.cpp
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

