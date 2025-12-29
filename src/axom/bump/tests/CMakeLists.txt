# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Bump unit tests
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Specify list of tests
#------------------------------------------------------------------------------

set(gtest_bump_tests
    bump_clipfield.cpp
    bump_cutfield.cpp
    bump_coordset_extents.cpp
    bump_make_polyhedral_topology.cpp
    bump_mergemeshes.cpp
    bump_mesh_operations.cpp
    bump_node_to_zone_relation.cpp
    bump_planeslice.cpp
    bump_slicers.cpp
    bump_topology_mapper.cpp
    bump_utilities.cpp
    bump_views.cpp
    bump_views_indexing.cpp
    )

set(bump_tests_depends_on 
    slic
    bump 
    gtest)

#------------------------------------------------------------------------------
# Add gtest based tests
#------------------------------------------------------------------------------
foreach(test ${gtest_bump_tests})
    get_filename_component( test_name ${test} NAME_WE )
    axom_add_executable( NAME ${test_name}_test
                         SOURCES ${test}
                         OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                         DEPENDS_ON ${bump_tests_depends_on}
                         FOLDER axom/bump/tests )

    axom_add_test( NAME ${test_name}
                   COMMAND ${test_name}_test )
endforeach()

