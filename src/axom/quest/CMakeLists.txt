# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Quest component
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Check necessary dependencies
#
# Note: Quest also optionally depends on Klee, when Klee is enabled
#------------------------------------------------------------------------------
axom_component_requires(NAME       QUEST
                        COMPONENTS MINT PRIMAL SLAM SLIC SPIN)

#------------------------------------------------------------------------------
# Specify headers, sources and dependencies
#------------------------------------------------------------------------------
set( quest_headers

    Delaunay.hpp
    detail/DelaunayElementFinder.hpp
    detail/DelaunayImpl.hpp
    detail/DelaunayInsertionHelper.hpp
    detail/DelaunayPointLocation.hpp
    detail/DelaunayValidation.hpp
    LinearizeCurves.hpp
    SignedDistance.hpp

    ## All-nearest-neighbors query
    AllNearestNeighbors.hpp
    detail/AllNearestNeighbors_detail.hpp

    ## Shaping
    Discretize.hpp
    detail/Discretize_detail.hpp

    ## In/out query
    InOutOctree.hpp
    detail/inout/BlockData.hpp
    detail/inout/MeshWrapper.hpp
    detail/inout/InOutOctreeMeshDumper.hpp
    detail/inout/InOutOctreeStats.hpp
    detail/inout/InOutOctreeValidator.hpp

    # Mesh tester
    MeshTester.hpp
    detail/MeshTester_detail.hpp

    # PointInCell
    PointInCell.hpp
    detail/PointFinder.hpp
    detail/PointInCellMeshWrapper_mfem.hpp

    ## File I/O
    io/STLReader.hpp
    io/STLWriter.hpp
    io/ProEReader.hpp

    ## Generalized Winding Numbers
    FastApproximateGWN.hpp

    ## quest interface
    interface/internal/mpicomm_wrapper.hpp
    interface/internal/QuestHelpers.hpp
    interface/inout.hpp
    interface/signed_distance.hpp

    util/mesh_helpers.hpp
    MeshViewUtil.hpp
    )

set( quest_sources

    ## Signed-distance query
    SignedDistance.cpp

    ## All-nearest-neighbors query
    AllNearestNeighbors.cpp

    ## Discretize shapes
    Discretize.cpp

    LinearizeCurves.cpp

    ## Mesh tester
    MeshTester.cpp

    ## File I/O
    io/STLReader.cpp
    io/STLWriter.cpp
    io/ProEReader.cpp

    ## quest interface
    interface/internal/QuestHelpers.cpp
    interface/inout.cpp
    interface/signed_distance.cpp

    util/mesh_helpers.cpp
   )

if(MFEM_FOUND)
    list(APPEND quest_headers io/MFEMReader.hpp)
    list(APPEND quest_headers GWNMethods.hpp)
    list(APPEND quest_sources io/MFEMReader.cpp)
endif()


set( quest_depends_on
    spin
    mint
    fmt
    )

blt_list_append(TO quest_depends_on IF AXOM_ENABLE_SIDRE ELEMENTS sidre)
blt_list_append(TO quest_depends_on IF AXOM_ENABLE_BUMP ELEMENTS bump)
blt_list_append(TO quest_depends_on IF SPARSEHASH_FOUND ELEMENTS sparsehash)
blt_list_append(TO quest_depends_on IF MFEM_FOUND ELEMENTS mfem)

if(AXOM_ENABLE_SIDRE)
    list(APPEND quest_headers    ScatteredInterpolation.hpp )
endif()

if(CONDUIT_FOUND AND AXOM_ENABLE_MPI)
    blt_list_append(
        TO       quest_headers
        ELEMENTS DistributedClosestPoint.hpp
                 detail/DistributedClosestPointImpl.hpp
                   )
    blt_list_append(
        TO       quest_sources
        ELEMENTS DistributedClosestPoint.cpp
                   )
    list(APPEND quest_depends_on conduit::conduit
                                 conduit::conduit_mpi)
endif()

blt_list_append(
    TO       quest_headers
    ELEMENTS MarchingCubes.hpp detail/MarchingCubesSingleDomain.hpp detail/MarchingCubesImpl.hpp
    IF       CONDUIT_FOUND
               )

blt_list_append(
    TO       quest_sources
    ELEMENTS MarchingCubes.cpp detail/MarchingCubesSingleDomain.cpp
    IF       CONDUIT_FOUND
               )

blt_list_append( TO quest_depends_on ELEMENTS conduit::conduit IF CONDUIT_FOUND )

if(AXOM_ENABLE_KLEE AND AXOM_ENABLE_SIDRE)
    if(MFEM_FOUND OR CONDUIT_FOUND)
        list(APPEND quest_headers Shaper.hpp
                                  DiscreteShape.hpp
                                  IntersectionShaper.hpp
                                  detail/shaping/shaping_helpers.hpp
                                  )
        list(APPEND quest_sources Shaper.cpp
                                  DiscreteShape.cpp
                                  detail/shaping/shaping_helpers.cpp
                                  )
        list(APPEND quest_depends_on klee)
    endif()

    # Geometry clipping requires Conduit, Sidre, and RAJA.
    # (TetMeshClipper additionally requires the Bump component.)
    if(CONDUIT_FOUND AND RAJA_FOUND)
        list(APPEND quest_headers ShapeMesh.hpp
                                  MeshClipper.hpp
                                  MeshClipperStrategy.hpp
                                  detail/clipping/Plane3DClipper.hpp
                                  detail/clipping/TetClipper.hpp
                                  detail/clipping/HexClipper.hpp
                                  detail/clipping/SphereClipper.hpp
                                  detail/clipping/MonotonicZSORClipper.hpp
                                  detail/clipping/SORClipper.hpp
                                  detail/clipping/MeshClipperImpl.hpp
                                  util/make_clipper_strategy.hpp)

        list(APPEND quest_sources ShapeMesh.cpp
                                  MeshClipper.cpp
                                  MeshClipperStrategy.cpp
                                  detail/clipping/Plane3DClipper.cpp
                                  detail/clipping/TetClipper.cpp
                                  detail/clipping/HexClipper.cpp
                                  detail/clipping/SphereClipper.cpp
                                  detail/clipping/MonotonicZSORClipper.cpp
                                  detail/clipping/SORClipper.cpp
                                  detail/clipping/MeshClipperImpl.hpp
                                  util/make_clipper_strategy.cpp)

        # MeshClipperImpl uses RAJA algorithms directly.
        list(APPEND quest_depends_on RAJA)
        # MeshClipperStrategy and the clippers depend on Klee geometry.
        list(APPEND quest_depends_on klee)

        if(AXOM_ENABLE_BUMP)
            list(APPEND quest_headers detail/clipping/TetMeshClipper.hpp)
            list(APPEND quest_sources detail/clipping/TetMeshClipper.cpp)
        endif()
    endif()

    if(MFEM_FOUND)
        list(APPEND quest_headers SamplingShaper.hpp 
                                  detail/shaping/InOutSampler.hpp
                                  detail/shaping/PrimitiveSampler.hpp
                                  detail/shaping/WindingNumberSampler.hpp
                                  )
    endif()
endif()

if(C2C_FOUND)
    list(APPEND quest_headers       io/C2CReader.hpp)
    list(APPEND quest_sources       io/C2CReader.cpp)
    list(APPEND quest_depends_on    c2c)
endif()

if(OPENCASCADE_FOUND)
    list(APPEND quest_headers       io/STEPReader.hpp)
    list(APPEND quest_sources       io/STEPReader.cpp)
    list(APPEND quest_depends_on    opencascade)
endif()

if (AXOM_ENABLE_MPI)
    list(APPEND quest_headers io/PSTLReader.hpp
                              io/PProEReader.hpp)
    blt_list_append(TO quest_headers IF C2C_FOUND ELEMENTS io/PC2CReader.hpp)
    blt_list_append(TO quest_headers IF OPENCASCADE_FOUND ELEMENTS io/PSTEPReader.hpp)

    list(APPEND quest_sources io/PSTLReader.cpp
                              io/PProEReader.cpp)
    blt_list_append(TO quest_sources IF C2C_FOUND ELEMENTS io/PC2CReader.cpp)
    blt_list_append(TO quest_sources IF OPENCASCADE_FOUND ELEMENTS io/PSTEPReader.cpp)
endif()

if (SHROUD_FOUND)
    add_subdirectory(interface)
endif()

if (ENABLE_FORTRAN)
    ## append shroud generated headers/sources to corresponding lists
    list(APPEND quest_headers interface/c_fortran/wrapQUEST.h)
    list(APPEND quest_headers interface/c_fortran/typesQUEST.h)
    list(APPEND quest_sources interface/c_fortran/wrapQUEST.cpp)
    list(APPEND quest_sources interface/c_fortran/utilQUEST.cpp)
    list(APPEND quest_sources interface/c_fortran/wrapfquest.F)
endif()

if (BUILD_SHARED_LIBS)
    # Install setup.py file used to build quest Python module.
    configure_file(interface/python/setup.py.in
                   ${CMAKE_CURRENT_BINARY_DIR}/interface/python/setup.py)
    configure_file(interface/python/quest_test.py.in
                   ${CMAKE_CURRENT_BINARY_DIR}/interface/python/quest_test.py)
endif()

#------------------------------------------------------------------------------
# Make and install the library
#------------------------------------------------------------------------------
axom_add_library(
    NAME        quest
    SOURCES     ${quest_sources}
    HEADERS     ${quest_headers}
    DEPENDS_ON  ${quest_depends_on}
    FOLDER      axom/quest
    )

axom_write_unified_header(NAME    quest
                          HEADERS ${quest_headers}
                          )

axom_install_component(NAME    quest
                       HEADERS ${quest_headers}
                       )

#------------------------------------------------------------------------------
# Add tests and examples
#------------------------------------------------------------------------------
if (AXOM_ENABLE_EXAMPLES)
    add_subdirectory(examples)
endif()

if (ENABLE_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

if (AXOM_ENABLE_TESTS)
    add_subdirectory(tests)
endif()
