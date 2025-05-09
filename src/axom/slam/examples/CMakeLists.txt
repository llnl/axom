# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Slam examples
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Single source standalone examples
#------------------------------------------------------------------------------
# The following flags might be helpful for analyzing performance:
#   for TemplateEx:        -S -fverbose-asm
#   for UnstructMeshField: -fwhole-program -funroll-loops
#------------------------------------------------------------------------------

set(example_sources
    HandleMesh.cpp
    PolicyPrototype.cpp
    ShockTube.cpp
    UserDocs.cpp
    )

foreach(example_source ${example_sources})
    get_filename_component(example_name ${example_source} NAME_WE)
    axom_add_executable(
        NAME        slam_${example_name}_ex
        SOURCES     ${example_source}
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  slam
        FOLDER      axom/slam/examples )

    if(AXOM_ENABLE_TESTS)
        axom_add_test(
            NAME        slam_${example_name}
            COMMAND     slam_${example_name}_ex )
    endif()
endforeach()


#------------------------------------------------------------------------------
# Slam version of tinyHydro
#------------------------------------------------------------------------------
add_subdirectory(tinyHydro)


#------------------------------------------------------------------------------
# Unstructured hexahedral mesh example
#------------------------------------------------------------------------------
axom_add_executable(
    NAME        slam_unstructMesh_ex
    SOURCES     UnstructMeshField.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  slam fmt
    FOLDER      axom/slam/examples )

if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR)
    axom_add_test(
        NAME        slam_unstructMesh
        COMMAND     slam_unstructMesh_ex ${AXOM_DATA_DIR}/slam )
endif()


#------------------------------------------------------------------------------
# Slam version of lulesh 2.0
# Note: Disabled on windows builds due to Lulesh's usage of unix APIs.
#------------------------------------------------------------------------------
if(NOT WIN32)
    add_subdirectory(lulesh2.0.3)

    # Original lulesh example
    set(_buildOrigLulesh TRUE)          # ---- Should we build original Lulesh?
    if ( ${_buildOrigLulesh} )
       add_subdirectory(lulesh2.0.3_orig)
    endif()
    unset(_buildOrigLulesh)
endif()
