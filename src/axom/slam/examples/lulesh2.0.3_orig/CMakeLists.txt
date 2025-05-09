# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Original version of lulesh 2.0, with a few light edits
#------------------------------------------------------------------------------

set(lulesh_headers
        lulesh.h
        )

set(lulesh_sources
        lulesh.cc
        lulesh-comm.cc
        lulesh-init.cc
        lulesh-util.cc
        lulesh-viz.cc
        )

# Add some extra configuration flags for MPI, OpenMP and allowing unused
# parameters (in 3rd party source code)
set(luleshExtraCompileFlags ${AXOM_DISABLE_UNUSED_PARAMETER_WARNINGS})
if (NOT AXOM_ENABLE_OPENMP )
    set(luleshExtraCompileFlags "${luleshExtraCompileFlags} ${AXOM_DISABLE_OMP_PRAGMA_WARNINGS}")
endif()

set(lulesh_depends_on)
blt_list_append(TO lulesh_depends_on ELEMENTS mpi IF AXOM_ENABLE_MPI)

set_source_files_properties(${lulesh_sources} ${lulesh_headers}
    PROPERTIES COMPILE_FLAGS " ${luleshExtraCompileFlags} " )

if ( lulesh_depends_on )
    axom_add_executable(
        NAME        slam_lulesh_orig_ex
        SOURCES     ${lulesh_sources}
        HEADERS     ${lulesh_headers}
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${lulesh_depends_on}
        FOLDER      axom/slam/examples )
else()
    axom_add_executable(
        NAME        slam_lulesh_orig_ex
        SOURCES     ${lulesh_sources}
        HEADERS     ${lulesh_headers}
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        FOLDER      axom/slam/examples )
endif()

if ( AXOM_ENABLE_MPI )
    target_compile_definitions(slam_lulesh_orig_ex PRIVATE USE_MPI=1)
else()
    target_compile_definitions(slam_lulesh_orig_ex PRIVATE USE_MPI=0)
endif()
